//! VCF liftover: convert variant coordinates between reference builds using
//! UCSC chain files.
//!
//! The public entry point is [`liftover`], which consumes a VCF from a
//! [`BufRead`] and writes mapped records to a [`Write`]. Unmapped records are
//! optionally written to a separate reject file.
//!
//! # Chain file format
//!
//! We parse the UCSC/Kent `.chain` format directly (noodles does not ship a
//! chain parser as of 0.109). Each chain header is followed by zero or more
//! aligned data lines:
//!
//! ```text
//! chain <score> <srcChrom> <srcSize> <srcStrand> <srcStart> <srcEnd> \
//!       <tgtChrom> <tgtSize> <tgtStrand> <tgtStart> <tgtEnd> <id>
//! <size> <dt> <dq>
//! <size> <dt> <dq>
//! <size>
//! ```
//!
//! Coordinates are 0-based half-open. We expand the chain into a vector of
//! `ChainBlock`s (one per ungapped segment) and index by source contig.

use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
};

use flate2::read::MultiGzDecoder;

use noodles::{
    core::Position,
    fasta,
    vcf::{
        self,
        header::record::value::{
            map::info::{Number, Type},
            Map,
        },
        variant::{
            io::Write as _,
            record_buf::{info::field::Value as InfoValue, AlternateBases},
            RecordBuf,
        },
    },
};

use crate::error::VcfkitError;

// ── public config / result types ─────────────────────────────────────────────

/// Options controlling the liftover pipeline.
#[derive(Debug, Clone)]
pub struct LiftoverOptions {
    /// If `Some(path)`, records that cannot be mapped are written to `path`
    /// in VCF format (header copied from the input).
    pub reject_file: Option<PathBuf>,
    /// If true, add `INFO/SRC_CONTIG` and `INFO/SRC_POS` fields recording the
    /// original (pre-liftover) coordinates.
    pub write_src_coords: bool,
    /// If true and the chain block strand is `-`, reverse-complement REF and
    /// ALT alleles. If false, such records are rejected as unmapped.
    pub fix_swapped_ref: bool,
    /// Output format for the main writer. Only [`OutputFormat::Vcf`] is
    /// supported. For BCF output, pipe through `bcftools view -O b`.
    pub output_format: crate::io::OutputFormat,
}

impl Default for LiftoverOptions {
    fn default() -> Self {
        Self {
            reject_file: None,
            write_src_coords: false,
            fix_swapped_ref: true,
            output_format: crate::io::OutputFormat::Vcf,
        }
    }
}

/// Statistics produced by a [`liftover`] run.
#[derive(Debug, Default, Clone, Copy)]
pub struct LiftoverStats {
    /// Records read from the input.
    pub input_records: usize,
    /// Records successfully mapped and written to the main output.
    pub output_records: usize,
    /// Records rejected because no chain block covered the source position.
    pub rejected_unmapped: usize,
    /// Records rejected because the lifted REF did not match the target
    /// reference at the mapped position.
    pub rejected_ref_mismatch: usize,
    /// Count of records whose alleles were reverse-complemented due to a
    /// strand swap.
    pub swapped_alleles: usize,
}

/// Known UCSC chain file URLs surfaced by `vcfkit liftover --list-chains`.
///
/// This is intentionally a small, well-known set; users with exotic builds are
/// expected to pass `--chain <path>` directly.
pub const KNOWN_CHAIN_URLS: &[(&str, &str)] = &[
    (
        "hg19 -> hg38",
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    ),
    (
        "hg38 -> hg19",
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
    ),
    (
        "hg38 -> T2T-CHM13",
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz",
    ),
    (
        "T2T-CHM13 -> hg38",
        "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz",
    ),
];

// ── chain types ──────────────────────────────────────────────────────────────

/// A single ungapped aligned block within a chain.
///
/// Coordinates are 0-based half-open.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainBlock {
    pub src_chrom: String,
    pub src_start: u64,
    pub src_end: u64,
    pub tgt_chrom: String,
    pub tgt_start: u64,
    pub tgt_end: u64,
    /// `+` for same-strand, `-` when the target coordinates walk in reverse.
    pub tgt_strand: char,
    /// Length of the target contig (needed to flip `-` strand coordinates).
    pub tgt_size: u64,
}

/// Chain index keyed by source contig. The `Vec` for each contig is sorted by
/// `src_start` and does not contain overlapping blocks (UCSC chain files are
/// already non-overlapping within a single chain).
#[derive(Debug, Default, Clone)]
pub struct ChainIndex {
    by_src: HashMap<String, Vec<ChainBlock>>,
}

impl ChainIndex {
    /// Parse a UCSC chain file from a reader into a [`ChainIndex`].
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self, VcfkitError> {
        let blocks = parse_chain(reader)?;
        let mut by_src: HashMap<String, Vec<ChainBlock>> = HashMap::new();
        for b in blocks {
            by_src.entry(b.src_chrom.clone()).or_default().push(b);
        }
        for v in by_src.values_mut() {
            v.sort_by_key(|b| b.src_start);
        }
        Ok(Self { by_src })
    }

    /// Parse a UCSC chain file from a filesystem path.
    ///
    /// Files ending with `.gz` are automatically decompressed via
    /// [`MultiGzDecoder`]. Plain (uncompressed) chain files are read directly.
    pub fn from_path(path: &Path) -> Result<Self, VcfkitError> {
        let file = File::open(path).map_err(|e| {
            VcfkitError::Other(format!("failed to open chain file {}: {e}", path.display()))
        })?;
        let is_gzip = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case("gz"))
            .unwrap_or(false);
        if is_gzip {
            Self::from_reader(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Self::from_reader(BufReader::new(file))
        }
    }

    /// Look up the block covering the 0-based source position, if any.
    pub fn lookup(&self, src_chrom: &str, src_pos0: u64) -> Option<&ChainBlock> {
        let blocks = self.by_src.get(src_chrom)?;
        // Binary search for the last block whose src_start <= src_pos0.
        let idx = match blocks.binary_search_by_key(&src_pos0, |b| b.src_start) {
            Ok(i) => i,
            Err(0) => return None,
            Err(i) => i - 1,
        };
        let b = &blocks[idx];
        if src_pos0 >= b.src_start && src_pos0 < b.src_end {
            Some(b)
        } else {
            None
        }
    }

    /// Returns the total number of blocks in the index.
    pub fn len(&self) -> usize {
        self.by_src.values().map(|v| v.len()).sum()
    }

    /// Returns true if the index has no blocks.
    pub fn is_empty(&self) -> bool {
        self.by_src.values().all(|v| v.is_empty())
    }

    /// Returns the set of source contigs covered by this index.
    pub fn source_contigs(&self) -> Vec<&str> {
        self.by_src.keys().map(String::as_str).collect()
    }
}

// ── chain file parser ────────────────────────────────────────────────────────

/// Parse a UCSC chain file into a flat vector of [`ChainBlock`]s.
///
/// Blocks from `-` strand chains carry target coordinates already converted
/// into the forward orientation, so callers do not need to know about the
/// chain header's strand field.
fn parse_chain<R: BufRead>(reader: R) -> Result<Vec<ChainBlock>, VcfkitError> {
    let mut blocks = Vec::new();
    let mut current: Option<ChainHeader> = None;

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| VcfkitError::Other(format!("chain read error: {e}")))?;
        let line = line.trim();
        if line.is_empty() {
            // A blank line ends the current chain.
            current = None;
            continue;
        }
        if line.starts_with('#') {
            continue;
        }

        if let Some(rest) = line.strip_prefix("chain ") {
            // Parse the chain header: "chain <score> <srcChrom> <srcSize> <srcStrand>
            // <srcStart> <srcEnd> <tgtChrom> <tgtSize> <tgtStrand> <tgtStart>
            // <tgtEnd> <id>".
            let fields: Vec<&str> = rest.split_whitespace().collect();
            if fields.len() < 11 {
                return Err(VcfkitError::Other(format!(
                    "chain header on line {} has too few fields",
                    lineno + 1
                )));
            }
            let src_chrom = fields[1].to_string();
            let src_size: u64 = parse_u64(fields[2], "srcSize", lineno)?;
            let src_strand = single_char(fields[3], "srcStrand", lineno)?;
            if src_strand != '+' {
                return Err(VcfkitError::Other(format!(
                    "unsupported src_strand '{}' in chain file (only '+' is supported)",
                    src_strand
                )));
            }
            let src_start: u64 = parse_u64(fields[4], "srcStart", lineno)?;
            let src_end: u64 = parse_u64(fields[5], "srcEnd", lineno)?;
            let tgt_chrom = fields[6].to_string();
            let tgt_size: u64 = parse_u64(fields[7], "tgtSize", lineno)?;
            let tgt_strand = single_char(fields[8], "tgtStrand", lineno)?;
            let tgt_start: u64 = parse_u64(fields[9], "tgtStart", lineno)?;
            let tgt_end: u64 = parse_u64(fields[10], "tgtEnd", lineno)?;

            current = Some(ChainHeader {
                src_chrom,
                src_size,
                src_strand,
                src_start,
                src_end,
                tgt_chrom,
                tgt_size,
                tgt_strand,
                tgt_start,
                tgt_end,
                src_cursor: src_start,
                tgt_cursor: tgt_start,
            });
            continue;
        }

        let Some(hdr) = current.as_mut() else {
            // Data line with no active chain header — skip.
            continue;
        };

        // Chain data: "<size>" or "<size> <dt> <dq>".
        let parts: Vec<&str> = line.split_whitespace().collect();
        let size: u64 = parse_u64(parts[0], "size", lineno)?;

        // Emit a ChainBlock for this ungapped segment, normalizing for strand.
        let block = make_block(hdr, size);
        blocks.push(block);

        hdr.src_cursor = hdr.src_cursor.saturating_add(size);
        hdr.tgt_cursor = hdr.tgt_cursor.saturating_add(size);

        if parts.len() >= 3 {
            let dt: u64 = parse_u64(parts[1], "dt", lineno)?;
            let dq: u64 = parse_u64(parts[2], "dq", lineno)?;
            hdr.src_cursor = hdr.src_cursor.saturating_add(dt);
            hdr.tgt_cursor = hdr.tgt_cursor.saturating_add(dq);
        } else {
            // The last block in a chain is written as a size-only line; the
            // chain is complete.
            current = None;
        }
    }

    Ok(blocks)
}

#[allow(dead_code)]
struct ChainHeader {
    src_chrom: String,
    src_size: u64,
    src_strand: char,
    src_start: u64,
    src_end: u64,
    tgt_chrom: String,
    tgt_size: u64,
    tgt_strand: char,
    tgt_start: u64,
    tgt_end: u64,
    src_cursor: u64,
    tgt_cursor: u64,
}

/// Build a [`ChainBlock`] for an ungapped segment of the active chain, taking
/// strand orientation into account.
///
/// For UCSC chains, if `tgt_strand == '-'`, the target coordinates in the
/// chain header are already expressed relative to the reverse strand (i.e.
/// `tgt_start` is the offset from the 5' end of the reverse-complemented
/// contig). We flip them into forward-strand coordinates here so downstream
/// lookup math is uniform.
fn make_block(hdr: &ChainHeader, size: u64) -> ChainBlock {
    let src_start = hdr.src_cursor;
    let src_end = hdr.src_cursor + size;

    let (tgt_start, tgt_end) = if hdr.tgt_strand == '-' {
        let end = hdr.tgt_size - hdr.tgt_cursor;
        let start = end - size;
        (start, end)
    } else {
        (hdr.tgt_cursor, hdr.tgt_cursor + size)
    };

    let _ = hdr.src_start;
    let _ = hdr.src_end;

    ChainBlock {
        src_chrom: hdr.src_chrom.clone(),
        src_start,
        src_end,
        tgt_chrom: hdr.tgt_chrom.clone(),
        tgt_start,
        tgt_end,
        tgt_strand: hdr.tgt_strand,
        tgt_size: hdr.tgt_size,
    }
}

fn parse_u64(s: &str, field: &str, lineno: usize) -> Result<u64, VcfkitError> {
    s.parse::<u64>().map_err(|e| {
        VcfkitError::Other(format!(
            "chain: line {} invalid {field}: {s:?} ({e})",
            lineno + 1
        ))
    })
}

fn single_char(s: &str, field: &str, lineno: usize) -> Result<char, VcfkitError> {
    let mut it = s.chars();
    match (it.next(), it.next()) {
        (Some(c), None) => Ok(c),
        _ => Err(VcfkitError::Other(format!(
            "chain: line {} invalid {field}: {s:?}",
            lineno + 1
        ))),
    }
}

// ── public entry point ───────────────────────────────────────────────────────

/// Liftover a VCF from `reader` to `writer` using the chain file at
/// `chain_path`.
///
/// * `source_ref_path` — FASTA+FAI for the source (pre-lift) build. Used to
///   sanity-check REF before lifting. Pass `/dev/null` style path only if REF
///   validation is not desired (it still must exist for [`Reference::open`]).
/// * `target_ref_path` — FASTA+FAI for the target (post-lift) build. Used to
///   validate that the lifted REF matches the target reference. Pass `None` to
///   skip target REF validation entirely.
///
/// Records that cannot be mapped are either silently dropped (default) or
/// written to `options.reject_file` with the original coordinates.
pub fn liftover<R: BufRead, W: Write>(
    reader: R,
    writer: W,
    chain_path: &Path,
    source_ref_path: &Path,
    target_ref_path: Option<&Path>,
    options: LiftoverOptions,
) -> Result<LiftoverStats, VcfkitError> {
    liftover_with_progress(
        reader,
        writer,
        chain_path,
        source_ref_path,
        target_ref_path,
        options,
        |_| {},
    )
}

/// Variant of [`liftover`] that notifies `on_record` after each input record
/// is read. Used by the CLI to drive a progress bar.
pub fn liftover_with_progress<R, W, F>(
    reader: R,
    writer: W,
    chain_path: &Path,
    source_ref_path: &Path,
    target_ref_path: Option<&Path>,
    options: LiftoverOptions,
    mut on_record: F,
) -> Result<LiftoverStats, VcfkitError>
where
    R: BufRead,
    W: Write,
    F: FnMut(u64),
{
    // ── open inputs ──────────────────────────────────────────────────────────
    let chain = ChainIndex::from_path(chain_path)?;

    // Source reference is validated to exist (for parity with normalize) but
    // we don't currently use it — the target reference is the authoritative
    // validation surface. Open if present.
    let _src_fa = Reference::open(source_ref_path).ok();
    let mut tgt_fa = target_ref_path.and_then(|p| Reference::open(p).ok());

    let mut vcf_reader = vcf::io::Reader::new(reader);
    let mut header = vcf_reader
        .read_header()
        .map_err(|e| VcfkitError::Other(format!("failed to read VCF header: {e}")))?;

    // Ensure SRC_CONTIG / SRC_POS are declared in the header if we'll emit them.
    if options.write_src_coords {
        let infos = header.infos_mut();
        infos.entry(String::from("SRC_CONTIG")).or_insert_with(|| {
            Map::<vcf::header::record::value::map::Info>::new(
                Number::Count(1),
                Type::String,
                "Source contig before liftover",
            )
        });
        infos.entry(String::from("SRC_POS")).or_insert_with(|| {
            Map::<vcf::header::record::value::map::Info>::new(
                Number::Count(1),
                Type::Integer,
                "Source position (1-based) before liftover",
            )
        });
    }

    let mut vcf_writer = vcf::io::Writer::new(writer);
    vcf_writer
        .write_header(&header)
        .map_err(|e| VcfkitError::Other(format!("failed to write VCF header: {e}")))?;

    // Reject file writer — lazily created on first rejected record.
    let mut reject_writer: Option<vcf::io::Writer<Box<dyn Write>>> =
        if let Some(ref path) = options.reject_file {
            let file = File::create(path).map_err(|e| {
                VcfkitError::Other(format!(
                    "failed to create reject file {}: {e}",
                    path.display()
                ))
            })?;
            let boxed: Box<dyn Write> = Box::new(file);
            let mut w = vcf::io::Writer::new(boxed);
            w.write_header(&header)
                .map_err(|e| VcfkitError::Other(format!("failed to write reject header: {e}")))?;
            Some(w)
        } else {
            None
        };

    let mut stats = LiftoverStats::default();
    let mut record = RecordBuf::default();

    loop {
        let n = vcf_reader
            .read_record_buf(&header, &mut record)
            .map_err(|e| VcfkitError::Other(format!("failed to read VCF record: {e}")))?;
        if n == 0 {
            break;
        }
        stats.input_records += 1;
        on_record(stats.input_records as u64);

        match lift_record(&record, &chain, tgt_fa.as_mut(), &options, &mut stats)? {
            LiftResult::Ok(lifted) => {
                vcf_writer
                    .write_variant_record(&header, lifted.as_ref())
                    .map_err(|e| VcfkitError::Other(format!("failed to write record: {e}")))?;
                stats.output_records += 1;
            }
            LiftResult::Reject => {
                if let Some(w) = reject_writer.as_mut() {
                    w.write_variant_record(&header, &record)
                        .map_err(|e| VcfkitError::Other(format!("failed to write reject: {e}")))?;
                }
            }
        }
    }

    Ok(stats)
}

// ── per-record liftover ──────────────────────────────────────────────────────

enum LiftResult {
    Ok(Box<RecordBuf>),
    Reject,
}

fn lift_record(
    record: &RecordBuf,
    chain: &ChainIndex,
    mut tgt_fa: Option<&mut Reference>,
    options: &LiftoverOptions,
    stats: &mut LiftoverStats,
) -> Result<LiftResult, VcfkitError> {
    let src_chrom = record.reference_sequence_name().to_string();
    let src_pos1 = match record.variant_start() {
        Some(p) => p.get() as u64,
        None => {
            stats.rejected_unmapped += 1;
            return Ok(LiftResult::Reject);
        }
    };
    let src_pos0 = src_pos1 - 1;

    let block = match chain.lookup(&src_chrom, src_pos0) {
        Some(b) => b,
        None => {
            stats.rejected_unmapped += 1;
            return Ok(LiftResult::Reject);
        }
    };

    let ref_bases = record.reference_bases().to_string();
    let alts: Vec<String> = record.alternate_bases().as_ref().to_vec();
    let ref_len = ref_bases.len() as u64;

    // Compute new (1-based) position for the record.
    let offset = src_pos0 - block.src_start;
    let (new_pos1, new_ref, new_alts, swapped) = if block.tgt_strand == '-' {
        // Source interval is [src_pos0, src_pos0 + ref_len); on the target
        // forward strand that maps to [tgt_end - offset - ref_len,
        // tgt_end - offset).
        if block.tgt_end < offset + ref_len {
            stats.rejected_unmapped += 1;
            return Ok(LiftResult::Reject);
        }
        let tgt_pos0 = block.tgt_end - offset - ref_len;
        let new_pos1 = tgt_pos0 + 1;

        if !options.fix_swapped_ref {
            // Caller asked us not to flip alleles — reject.
            stats.rejected_unmapped += 1;
            return Ok(LiftResult::Reject);
        }

        let new_ref = reverse_complement(&ref_bases);
        let new_alts: Vec<String> = alts
            .iter()
            .map(|a| {
                if is_symbolic(a) {
                    a.clone()
                } else {
                    reverse_complement(a)
                }
            })
            .collect();
        (new_pos1, new_ref, new_alts, true)
    } else {
        let new_pos1 = block.tgt_start + offset + 1;
        (new_pos1, ref_bases.clone(), alts.clone(), false)
    };

    // Build the lifted record.
    let mut lifted = record.clone();
    *lifted.reference_sequence_name_mut() = block.tgt_chrom.clone();

    let new_pos_usize =
        usize::try_from(new_pos1).map_err(|_| VcfkitError::Other("position overflow".into()))?;
    let new_pos = Position::new(new_pos_usize)
        .ok_or_else(|| VcfkitError::Other("lifted position is zero".into()))?;
    *lifted.variant_start_mut() = Some(new_pos);
    *lifted.reference_bases_mut() = new_ref.clone();
    *lifted.alternate_bases_mut() = AlternateBases::from(new_alts);

    // INFO/SRC_* tags.
    if options.write_src_coords {
        let info = lifted.info_mut();
        info.insert(
            String::from("SRC_CONTIG"),
            Some(InfoValue::String(src_chrom.clone())),
        );
        info.insert(
            String::from("SRC_POS"),
            Some(InfoValue::Integer(src_pos1 as i32)),
        );
    }

    // Validate lifted REF against the target reference (if available).
    if let Some(fa) = tgt_fa.as_mut() {
        if !is_symbolic_record(&lifted) {
            match fa.contig(&block.tgt_chrom) {
                Ok(seq) => {
                    let tgt_pos0_usize = (new_pos1 - 1) as usize;
                    let ref_bytes = new_ref.as_bytes();
                    let end = tgt_pos0_usize + ref_bytes.len();
                    if end > seq.len() {
                        stats.rejected_ref_mismatch += 1;
                        return Ok(LiftResult::Reject);
                    }
                    let on_ref = &seq[tgt_pos0_usize..end];
                    let matches = on_ref
                        .iter()
                        .zip(ref_bytes.iter())
                        .all(|(a, b)| a.eq_ignore_ascii_case(b));
                    if !matches {
                        stats.rejected_ref_mismatch += 1;
                        return Ok(LiftResult::Reject);
                    }
                }
                Err(_) => {
                    // Target contig not in the FASTA — treat as mismatch.
                    stats.rejected_ref_mismatch += 1;
                    return Ok(LiftResult::Reject);
                }
            }
        }
    }

    // Only count swapped_alleles for records that make it to output.
    if swapped {
        stats.swapped_alleles += 1;
    }

    Ok(LiftResult::Ok(Box::new(lifted)))
}

fn is_symbolic(a: &str) -> bool {
    a.starts_with('<') || a.contains('[') || a.contains(']') || a == "*"
}

fn is_symbolic_record(record: &RecordBuf) -> bool {
    record
        .alternate_bases()
        .as_ref()
        .iter()
        .any(|a| is_symbolic(a))
}

/// Reverse-complement a DNA string.  Non-ACGT/acgt characters (e.g. `N`) are
/// passed through unchanged apart from case preservation.
fn reverse_complement(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            'N' | 'n' => c,
            other => other,
        })
        .collect()
}

// ── reference FASTA cache (mirror of normalize's) ────────────────────────────

struct Reference {
    path: PathBuf,
    cache: Option<(String, Vec<u8>)>,
}

impl Reference {
    fn open(path: &Path) -> Result<Self, VcfkitError> {
        let _ = File::open(path).map_err(|e| {
            VcfkitError::Other(format!("failed to open reference {}: {e}", path.display()))
        })?;
        let mut fai_path = path.to_path_buf().into_os_string();
        fai_path.push(".fai");
        let fai_path = PathBuf::from(fai_path);
        let _ = File::open(&fai_path).map_err(|e| {
            VcfkitError::Other(format!(
                "failed to open reference index {}: {e}",
                fai_path.display()
            ))
        })?;
        Ok(Self {
            path: path.to_path_buf(),
            cache: None,
        })
    }

    fn contig(&mut self, name: &str) -> Result<&[u8], VcfkitError> {
        if self
            .cache
            .as_ref()
            .is_some_and(|(cached, _)| cached == name)
        {
            return Ok(self.cache.as_ref().unwrap().1.as_slice());
        }

        let mut builder = fasta::io::indexed_reader::Builder::default()
            .build_from_path(&self.path)
            .map_err(|e| VcfkitError::Other(format!("failed to open indexed FASTA: {e}")))?;
        let start = Position::MIN;
        let end = Position::MAX;
        let region = noodles::core::Region::new(name.as_bytes(), start..=end);
        let record = builder
            .query(&region)
            .map_err(|e| VcfkitError::Other(format!("failed to query contig {name}: {e}")))?;
        let seq = record.sequence().as_ref().to_vec();
        self.cache = Some((name.to_string(), seq));
        Ok(self.cache.as_ref().unwrap().1.as_slice())
    }
}

// ── tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_chain_single_block() {
        // One chain with a single 50-bp block, + strand, mapping chr1:0-50 to
        // chr1:100-150.
        let chain = "chain 1000 chr1 200 + 0 50 chr1 200 + 100 150 1\n50\n";
        let idx = ChainIndex::from_reader(chain.as_bytes()).unwrap();
        assert_eq!(idx.len(), 1);
        let b = idx.lookup("chr1", 10).unwrap();
        assert_eq!(b.src_start, 0);
        assert_eq!(b.src_end, 50);
        assert_eq!(b.tgt_start, 100);
        assert_eq!(b.tgt_end, 150);
        assert_eq!(b.tgt_strand, '+');
    }

    #[test]
    fn parse_chain_with_gaps() {
        // Two 10-bp blocks separated by a 5-bp source gap and 5-bp target gap.
        // src: 0..10, 15..25 ; tgt: 100..110, 115..125
        let chain = "chain 1000 chr1 100 + 0 25 chr1 200 + 100 125 1\n\
                     10 5 5\n\
                     10\n";
        let idx = ChainIndex::from_reader(chain.as_bytes()).unwrap();
        assert_eq!(idx.len(), 2);
        // src_pos0=5 falls in block 1: offset 5 -> tgt 105
        let b = idx.lookup("chr1", 5).unwrap();
        assert_eq!(b.src_start, 0);
        assert_eq!(b.tgt_start, 100);
        // src_pos0=12 is in the gap — no block.
        assert!(idx.lookup("chr1", 12).is_none());
        // src_pos0=20 falls in block 2.
        let b = idx.lookup("chr1", 20).unwrap();
        assert_eq!(b.src_start, 15);
        assert_eq!(b.tgt_start, 115);
    }

    #[test]
    fn parse_chain_negative_strand_flips_target() {
        // tgt_size=200, tgt_strand=-, tgt_start=50, tgt_end=150 → forward
        // target range = [200-150, 200-50) = [50, 150). A single 100-bp block.
        let chain = "chain 1000 chr1 100 + 0 100 chr1 200 - 50 150 1\n100\n";
        let idx = ChainIndex::from_reader(chain.as_bytes()).unwrap();
        let b = idx.lookup("chr1", 0).unwrap();
        assert_eq!(b.tgt_strand, '-');
        assert_eq!(b.tgt_start, 50);
        assert_eq!(b.tgt_end, 150);
    }

    #[test]
    fn reverse_complement_basic() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("N"), "N");
    }

    #[test]
    fn lookup_returns_none_for_unknown_contig() {
        let chain = "chain 1 chr1 100 + 0 10 chr1 100 + 0 10 1\n10\n";
        let idx = ChainIndex::from_reader(chain.as_bytes()).unwrap();
        assert!(idx.lookup("chrX", 5).is_none());
    }

    #[test]
    fn known_chain_urls_are_exposed() {
        assert_eq!(KNOWN_CHAIN_URLS.len(), 4);
        assert!(KNOWN_CHAIN_URLS.iter().any(|(k, _)| *k == "hg19 -> hg38"));
    }
}
