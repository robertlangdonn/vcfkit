//! VCF normalization: left-alignment of indels, multi-allelic splitting, and
//! REF-against-reference checks.
//!
//! The public entry point is [`normalize`], which consumes a VCF from a
//! [`BufRead`] and writes normalized records to a [`Write`]. The algorithm
//! follows Tan et al. 2015 for left-alignment and the bcftools convention for
//! splitting multi-allelic sites (INFO `Number=A/R/G` fields are re-sliced per
//! allele; `Number=1/./0` are copied verbatim).

use std::{
    fs::File,
    io::{BufRead, Write},
    path::{Path, PathBuf},
};

use noodles::{
    core::Position,
    fasta,
    vcf::{
        self,
        header::record::value::map::{format, info},
        variant::{
            RecordBuf,
            io::Write as _,
            record_buf::{
                AlternateBases, Samples,
                info::field::{Value as InfoValue, value::Array as InfoArray},
                samples::sample::{Value as SampleValue, value::Array as SampleArray},
            },
        },
    },
};

use crate::error::VcfkitError;

// ── public config / result types ─────────────────────────────────────────────

/// How to handle a mismatch between the VCF REF column and the reference FASTA.
///
/// Note: only the first base of REF is compared against the reference FASTA,
/// consistent with bcftools `norm -c` behavior.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RefCheck {
    /// Don't check REF against the reference.
    Ignore,
    /// Warn on stderr; count the mismatch but continue.
    #[default]
    Warn,
    /// Return an error on the first mismatch.
    Error,
}

/// Options controlling the normalization pipeline.
#[derive(Debug, Clone)]
pub struct NormalizeOptions {
    /// If true, split multi-allelic records into biallelic records.
    pub split_multiallelics: bool,
    /// If true, left-align indels using Tan et al. 2015.
    pub left_align: bool,
    /// How to respond when the VCF REF does not match the reference FASTA.
    /// Only the first base of REF is validated (consistent with bcftools behavior).
    pub check_ref: RefCheck,
    /// Output format to write to `writer`.
    pub output_format: crate::io::OutputFormat,
}

impl Default for NormalizeOptions {
    fn default() -> Self {
        Self {
            split_multiallelics: true,
            left_align: true,
            check_ref: RefCheck::Warn,
            output_format: crate::io::OutputFormat::Vcf,
        }
    }
}

/// Statistics produced by a `normalize` run.
#[derive(Debug, Default, Clone, Copy)]
pub struct NormalizeStats {
    /// Records read from the input.
    pub input_records: usize,
    /// Records written to the output (may exceed input after splitting).
    pub output_records: usize,
    /// Count of records whose position/REF/ALT changed during left-alignment.
    pub left_aligned: usize,
    /// Count of input sites that were split into biallelic records.
    pub split_sites: usize,
    /// Count of records whose REF mismatched the reference FASTA.
    pub ref_mismatches: usize,
    /// Count of records skipped because their POS exceeded the declared contig length.
    pub out_of_bounds: usize,
}

// ── public entry point ───────────────────────────────────────────────────────

/// Normalize a VCF read from `reader` and write the result to `writer`.
///
/// `reference_path` must point at a bgzf-indexed or plain FASTA alongside a
/// `.fai` index (SAMtools style).
///
/// Currently only emits plain VCF; `options.output_format` is reserved for a
/// future BCF writer.
pub fn normalize<R: BufRead, W: Write>(
    reader: R,
    writer: W,
    reference_path: &Path,
    options: NormalizeOptions,
) -> Result<NormalizeStats, VcfkitError> {
    normalize_with_progress(reader, writer, reference_path, options, |_| {})
}

/// Variant of [`normalize`] that invokes `on_record` after each *input*
/// record is consumed, passing the running count. Used by the CLI to drive a
/// progress bar; bench/test code should prefer [`normalize`].
pub fn normalize_with_progress<R, W, F>(
    reader: R,
    writer: W,
    reference_path: &Path,
    options: NormalizeOptions,
    mut on_record: F,
) -> Result<NormalizeStats, VcfkitError>
where
    R: BufRead,
    W: Write,
    F: FnMut(u64),
{
    let mut vcf_reader = vcf::io::Reader::new(reader);
    let header = vcf_reader
        .read_header()
        .map_err(|e| VcfkitError::Other(format!("failed to read VCF header: {e}")))?;

    let mut vcf_writer = vcf::io::Writer::new(writer);
    vcf_writer
        .write_header(&header)
        .map_err(|e| VcfkitError::Other(format!("failed to write VCF header: {e}")))?;

    let mut fasta = if options.check_ref != RefCheck::Ignore || options.left_align {
        Some(Reference::open(reference_path)?)
    } else {
        None
    };

    let mut stats = NormalizeStats::default();
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

        // Out-of-bounds check: if the VCF header declares a contig length and
        // the record's POS exceeds it, warn and skip rather than aborting.
        if let Some(pos) = record.variant_start() {
            let pos = pos.get();
            let chrom = record.reference_sequence_name();
            if let Some(contig_map) = header.contigs().get(chrom) {
                if let Some(contig_length) = contig_map.length() {
                    if pos > contig_length {
                        tracing::warn!(
                            "position {} out of bounds for contig {} (length {}); skipping",
                            pos, chrom, contig_length
                        );
                        stats.out_of_bounds += 1;
                        continue;
                    }
                }
            }
        }

        let processed = process_record(&record, &header, fasta.as_mut(), &options, &mut stats)?;

        for rec in processed {
            vcf_writer
                .write_variant_record(&header, &rec)
                .map_err(|e| VcfkitError::Other(format!("failed to write record: {e}")))?;
            stats.output_records += 1;
        }
    }

    Ok(stats)
}

// ── core algorithm ───────────────────────────────────────────────────────────

/// Per-record pipeline: REF check (once per site) → split → left-align.
fn process_record(
    record: &RecordBuf,
    header: &vcf::Header,
    mut fasta: Option<&mut Reference>,
    options: &NormalizeOptions,
    stats: &mut NormalizeStats,
) -> Result<Vec<RecordBuf>, VcfkitError> {
    // Symbolic or breakend ALTs — pass through unchanged.
    if has_symbolic_alt(record) {
        return Ok(vec![record.clone()]);
    }

    // 1. REF check — run ONCE on the original record before any splitting so
    //    that a triallelic site counts as a single mismatch (not one per split).
    let ref_mismatch = if options.check_ref != RefCheck::Ignore {
        if let Some(fa) = fasta.as_deref_mut() {
            check_ref(record, fa, options.check_ref)?
        } else {
            None
        }
    } else {
        None
    };
    if let Some(ref msg) = ref_mismatch {
        stats.ref_mismatches += 1;
        eprintln!("{msg}");
    }

    // 2. Split multi-allelics (if enabled and needed).
    let split_records: Vec<RecordBuf> = if options.split_multiallelics
        && record.alternate_bases().as_ref().len() > 1
    {
        stats.split_sites += 1;
        split_multiallelic(record, header)
    } else {
        vec![record.clone()]
    };

    // 3. For each (possibly-split) record, left-align.
    let mut out = Vec::with_capacity(split_records.len());
    for mut rec in split_records {
        if options.left_align {
            if let Some(fa) = fasta.as_deref_mut() {
                if left_align_record(&mut rec, fa)? {
                    stats.left_aligned += 1;
                }
            }
        }

        out.push(rec);
    }

    Ok(out)
}

/// Returns true if any ALT allele is a symbolic allele (`<…>`) or a breakend.
fn has_symbolic_alt(record: &RecordBuf) -> bool {
    record
        .alternate_bases()
        .as_ref()
        .iter()
        .any(|a| a.starts_with('<') || a.contains('[') || a.contains(']'))
}

// ── left-alignment (Tan et al. 2015) ─────────────────────────────────────────

/// Tan et al. 2015 left-alignment. Returns true when the record was modified.
fn left_align_record(record: &mut RecordBuf, fasta: &mut Reference) -> Result<bool, VcfkitError> {
    let mut pos = match record.variant_start() {
        Some(p) => p.get(),
        None => return Ok(false),
    };
    let mut r: Vec<u8> = record.reference_bases().as_bytes().to_vec();
    let mut a: Vec<u8> = match record.alternate_bases().as_ref().first() {
        Some(alt) => alt.as_bytes().to_vec(),
        None => return Ok(false),
    };

    // Only left-align indels (len(REF) != len(ALT)), and skip symbolic/breakend.
    if r.len() == a.len() {
        return Ok(false);
    }

    let chrom = record.reference_sequence_name().to_string();
    let seq = fasta.contig(&chrom)?;

    let start_pos = pos;
    let start_r = r.clone();
    let start_a = a.clone();

    while pos > 1 && !r.is_empty() && !a.is_empty() && r.last() == a.last() {
        let prev_base = match seq.get(pos - 2).copied() {
            Some(b) => b.to_ascii_uppercase(),
            None => break,
        };
        r.pop();
        a.pop();
        r.insert(0, prev_base);
        a.insert(0, prev_base);
        pos -= 1;
    }

    // Trim common prefix beyond the anchor base (bcftools: keep at least 1 base
    // on each side). This is not strictly Tan 2015, but matches the tests that
    // expect a minimal anchored representation.
    while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        pos += 1;
    }

    let changed = pos != start_pos || r != start_r || a != start_a;
    if changed {
        let new_pos = Position::new(pos).ok_or_else(|| {
            VcfkitError::Other(format!("invalid position {pos} after left-alignment"))
        })?;
        *record.variant_start_mut() = Some(new_pos);
        *record.reference_bases_mut() =
            String::from_utf8(r).map_err(|e| VcfkitError::Other(e.to_string()))?;
        let alt_str = String::from_utf8(a).map_err(|e| VcfkitError::Other(e.to_string()))?;
        *record.alternate_bases_mut() = AlternateBases::from(vec![alt_str]);
    }
    Ok(changed)
}

// ── REF check ────────────────────────────────────────────────────────────────

fn check_ref(
    record: &RecordBuf,
    fasta: &mut Reference,
    mode: RefCheck,
) -> Result<Option<String>, VcfkitError> {
    let chrom = record.reference_sequence_name().to_string();
    let seq = fasta.contig(&chrom)?;

    let pos = match record.variant_start() {
        Some(p) => p.get(),
        None => return Ok(None),
    };
    let vcf_ref = record.reference_bases();
    if vcf_ref.is_empty() {
        return Ok(None);
    }
    let want = vcf_ref.as_bytes()[0].to_ascii_uppercase();
    // If the position is beyond the FASTA sequence, treat as out-of-bounds.
    // For Error mode we still return an error; for Warn/Ignore we return None
    // (the record will have already been skipped at the loop level if the VCF
    // header declared a contig length, but we handle the FASTA-only case here).
    let have = match seq.get(pos - 1).copied() {
        Some(b) => b.to_ascii_uppercase(),
        None => match mode {
            RefCheck::Error => {
                return Err(VcfkitError::Other(format!(
                    "position {pos} out of bounds for contig {chrom}"
                )));
            }
            _ => return Ok(None),
        },
    };
    if want == have {
        return Ok(None);
    }

    let msg = format!(
        "REF mismatch at {chrom}:{pos}: VCF has {}, reference has {}",
        want as char, have as char
    );

    match mode {
        RefCheck::Ignore => Ok(None),
        RefCheck::Warn => Ok(Some(msg)),
        RefCheck::Error => Err(VcfkitError::Other(msg)),
    }
}

// ── multi-allelic splitting ──────────────────────────────────────────────────

/// Splits a record with N ALT alleles into N biallelic records, handling
/// `Number=A`, `Number=R`, and `Number=G` INFO/FORMAT fields per the VCF spec.
fn split_multiallelic(record: &RecordBuf, header: &vcf::Header) -> Vec<RecordBuf> {
    let alts = record.alternate_bases().as_ref().to_vec();
    let n_alts = alts.len();
    let mut out = Vec::with_capacity(n_alts);

    for (alt_idx, alt) in alts.iter().enumerate() {
        // allele_index = alt_idx + 1 (0 = REF)
        let allele_index = alt_idx + 1;

        let mut new = record.clone();
        *new.alternate_bases_mut() = AlternateBases::from(vec![alt.clone()]);

        // Rewrite INFO fields.
        let mut new_info = noodles::vcf::variant::record_buf::Info::default();
        for (key, value) in record.info().as_ref().iter() {
            let number = header
                .infos()
                .get(key)
                .map(info_number)
                .unwrap_or(FieldNumber::Other);
            let new_value = split_info_value(value.as_ref(), number, n_alts, allele_index);
            new_info.insert(key.clone(), new_value);
        }
        *new.info_mut() = new_info;

        // Rewrite FORMAT fields per sample.
        let keys = record.samples().keys().clone();
        let old_values = record.samples().clone();
        let mut new_values: Vec<Vec<Option<SampleValue>>> =
            Vec::with_capacity(old_values.values().count());

        for sample in old_values.values() {
            let mut row: Vec<Option<SampleValue>> = Vec::with_capacity(keys.as_ref().len());
            for (i, key) in keys.as_ref().iter().enumerate() {
                let number = header
                    .formats()
                    .get(key)
                    .map(format_number)
                    .unwrap_or(FieldNumber::Other);
                let v = sample.values().get(i).cloned().flatten();
                let new_v = split_sample_value(v.as_ref(), number, n_alts, allele_index, key);
                row.push(new_v);
            }
            new_values.push(row);
        }
        *new.samples_mut() = Samples::new(keys, new_values);

        out.push(new);
    }

    out
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum FieldNumber {
    /// `Number=A` (one value per ALT).
    Alt,
    /// `Number=R` (one value per REF + ALT).
    Ref,
    /// `Number=G` (one value per genotype).
    Genotype,
    /// Anything else (fixed count, unknown, flag, etc.).
    Other,
}

fn info_number(map: &noodles::vcf::header::record::value::Map<info::Info>) -> FieldNumber {
    use noodles::vcf::header::record::value::map::info::Number;
    match map.number() {
        Number::AlternateBases => FieldNumber::Alt,
        Number::ReferenceAlternateBases => FieldNumber::Ref,
        Number::Samples => FieldNumber::Genotype,
        _ => FieldNumber::Other,
    }
}

fn format_number(map: &noodles::vcf::header::record::value::Map<format::Format>) -> FieldNumber {
    use noodles::vcf::header::record::value::map::format::Number;
    match map.number() {
        Number::AlternateBases => FieldNumber::Alt,
        Number::ReferenceAlternateBases => FieldNumber::Ref,
        Number::Samples => FieldNumber::Genotype,
        _ => FieldNumber::Other,
    }
}

/// Returns the g-indices (within the parent record's G-ordered value list) to
/// keep for a biallelic split with `allele_index` (1-based into the original
/// allele vector `[REF, ALT1, ALT2, …]`).
///
/// The VCF G ordering for ploidy P with K alleles enumerates genotypes in
/// "colex" order; for diploid the sequence is (0/0, 0/1, 1/1, 0/2, 1/2, 2/2,
/// 0/3, 1/3, 2/3, 3/3, …). After splitting, the biallelic record has alleles
/// [REF, ALT_k], so we keep the g-indices that correspond to genotypes whose
/// allele set is {0, k}: namely (0/0), (0/k), (k/k).
fn diploid_g_indices_to_keep(n_alts: usize, allele_index: usize) -> Vec<usize> {
    // gt_index(a, b) with a <= b is b*(b+1)/2 + a.
    let gi = |a: usize, b: usize| {
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        b * (b + 1) / 2 + a
    };
    let n_alleles = n_alts + 1;
    // Sanity: callers should only use this when n_alleles >= 2.
    let k = allele_index;
    if k >= n_alleles {
        return vec![0];
    }
    vec![gi(0, 0), gi(0, k), gi(k, k)]
}

fn split_info_value(
    value: Option<&InfoValue>,
    number: FieldNumber,
    n_alts: usize,
    allele_index: usize,
) -> Option<InfoValue> {
    let v = value?;
    match number {
        FieldNumber::Alt => pick_info_array(v, &[allele_index - 1]),
        FieldNumber::Ref => pick_info_array(v, &[0, allele_index]),
        FieldNumber::Genotype => {
            let idxs = diploid_g_indices_to_keep(n_alts, allele_index);
            pick_info_array(v, &idxs)
        }
        FieldNumber::Other => Some(v.clone()),
    }
}

fn split_sample_value(
    value: Option<&SampleValue>,
    number: FieldNumber,
    n_alts: usize,
    allele_index: usize,
    _key: &str,
) -> Option<SampleValue> {
    let v = value?;
    match number {
        FieldNumber::Alt => pick_sample_array(v, &[allele_index - 1]),
        FieldNumber::Ref => pick_sample_array(v, &[0, allele_index]),
        FieldNumber::Genotype => {
            let idxs = diploid_g_indices_to_keep(n_alts, allele_index);
            pick_sample_array(v, &idxs)
        }
        FieldNumber::Other => Some(v.clone()),
    }
}

fn pick_info_array(value: &InfoValue, indices: &[usize]) -> Option<InfoValue> {
    // Slice an array-valued INFO field, preserving its scalar element type.
    let arr = match value {
        InfoValue::Array(a) => a,
        // Non-array scalar with Number=A/R/G is unusual but treat it as pass-through.
        _ => return Some(value.clone()),
    };
    let picked = match arr {
        InfoArray::Integer(xs) => InfoArray::Integer(pick(xs, indices)),
        InfoArray::Float(xs) => InfoArray::Float(pick(xs, indices)),
        InfoArray::Character(xs) => InfoArray::Character(pick(xs, indices)),
        InfoArray::String(xs) => InfoArray::String(pick(xs, indices)),
    };
    Some(InfoValue::Array(picked))
}

fn pick_sample_array(value: &SampleValue, indices: &[usize]) -> Option<SampleValue> {
    let arr = match value {
        SampleValue::Array(a) => a,
        _ => return Some(value.clone()),
    };
    let picked = match arr {
        SampleArray::Integer(xs) => SampleArray::Integer(pick(xs, indices)),
        SampleArray::Float(xs) => SampleArray::Float(pick(xs, indices)),
        SampleArray::Character(xs) => SampleArray::Character(pick(xs, indices)),
        SampleArray::String(xs) => SampleArray::String(pick(xs, indices)),
    };
    Some(SampleValue::Array(picked))
}

fn pick<T: Clone>(xs: &[Option<T>], indices: &[usize]) -> Vec<Option<T>> {
    indices
        .iter()
        .map(|&i| xs.get(i).cloned().unwrap_or(None))
        .collect()
}

// ── reference FASTA cache ────────────────────────────────────────────────────

struct Reference {
    path: PathBuf,
    cache: Option<(String, Vec<u8>)>,
}

impl Reference {
    fn open(path: &Path) -> Result<Self, VcfkitError> {
        // Probe once that the file and index exist.
        let _ = File::open(path).map_err(|e| {
            VcfkitError::Other(format!("failed to open reference {}: {e}", path.display()))
        })?;
        let mut fai_path = path.to_path_buf();
        let mut with_fai = fai_path.as_os_str().to_os_string();
        with_fai.push(".fai");
        fai_path = PathBuf::from(with_fai);
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

    /// Loads the full sequence for a contig (cached — most VCFs are sorted by
    /// contig so only one sequence is resident at a time).
    fn contig(&mut self, name: &str) -> Result<&[u8], VcfkitError> {
        if self
            .cache
            .as_ref()
            .is_some_and(|(cached_name, _)| cached_name == name)
        {
            return Ok(self.cache.as_ref().unwrap().1.as_slice());
        }

        // Build an indexed reader, query the whole contig.
        let mut builder = fasta::io::indexed_reader::Builder::default()
            .build_from_path(&self.path)
            .map_err(|e| VcfkitError::Other(format!("failed to open indexed FASTA: {e}")))?;

        // Build a "whole contig" region. Lacking an explicit length, we use
        // Position::MAX; IndexedReader clamps to the contig length internally.
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

    fn corpus_dir() -> PathBuf {
        // $CARGO_MANIFEST_DIR points at crates/vcfkit-core — walk up twice.
        let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        p.pop(); // crates/
        p.pop(); // workspace root
        p.push("tests/corpus/synthetic");
        p
    }

    fn normalize_to_string(input: &[u8], options: NormalizeOptions) -> (String, NormalizeStats) {
        let mut out = Vec::new();
        let fa = corpus_dir().join("mini_ref.fa");
        let stats = normalize(input, &mut out, &fa, options).expect("normalize ok");
        (String::from_utf8(out).expect("utf8"), stats)
    }

    // ── left-alignment ────────────────────────────────────────────────────────

    #[test]
    fn left_align_matches_corpus_norm_tags() {
        let input = std::fs::read(corpus_dir().join("indels_unnormalized.vcf")).unwrap();
        let opts = NormalizeOptions {
            split_multiallelics: false,
            left_align: true,
            check_ref: RefCheck::Ignore,
            output_format: crate::io::OutputFormat::Vcf,
        };
        let (out, stats) = normalize_to_string(&input, opts);

        // Parse output records and compare against the NORM_* tags.
        let mut reader = vcf::io::Reader::new(out.as_bytes());
        let header = reader.read_header().unwrap();
        let mut rec = RecordBuf::default();
        let mut checked = 0usize;
        while reader.read_record_buf(&header, &mut rec).unwrap() > 0 {
            let pos = rec.variant_start().unwrap().get();
            let r = rec.reference_bases().to_string();
            let a = rec.alternate_bases().as_ref()[0].clone();

            let info = rec.info();
            let norm_pos = match info.get("NORM_POS").flatten() {
                Some(InfoValue::Integer(n)) => *n as usize,
                other => panic!("missing/bad NORM_POS: {other:?}"),
            };
            let norm_ref = match info.get("NORM_REF").flatten() {
                Some(InfoValue::String(s)) => s.clone(),
                other => panic!("missing/bad NORM_REF: {other:?}"),
            };
            let norm_alt = match info.get("NORM_ALT").flatten() {
                Some(InfoValue::String(s)) => s.clone(),
                other => panic!("missing/bad NORM_ALT: {other:?}"),
            };
            assert_eq!(
                (pos, r.clone(), a.clone()),
                (norm_pos, norm_ref.clone(), norm_alt.clone()),
                "left-alignment mismatch for record (expected NORM_POS={norm_pos}, \
                 NORM_REF={norm_ref}, NORM_ALT={norm_alt})",
            );
            checked += 1;
        }
        assert_eq!(checked, 5, "expected 5 indel records");
        assert_eq!(stats.input_records, 5);
        assert_eq!(stats.output_records, 5);
        assert!(stats.left_aligned > 0);
    }

    // ── multi-allelic splitting ───────────────────────────────────────────────

    #[test]
    fn split_multi_allelic_number_a_and_r() {
        let input = std::fs::read(corpus_dir().join("multi_allelic.vcf")).unwrap();
        let opts = NormalizeOptions {
            split_multiallelics: true,
            left_align: false,
            check_ref: RefCheck::Ignore,
            output_format: crate::io::OutputFormat::Vcf,
        };
        let (out, stats) = normalize_to_string(&input, opts);
        assert_eq!(stats.input_records, 5);
        // 4 biallelic-after-split + 1 triallelic = 2+2+3+2+2 = 11 output records
        assert_eq!(stats.output_records, 11);
        assert_eq!(stats.split_sites, 5);

        // Spot-check the first record: A -> T,G with AF=0.3,0.2 AD=50,30,20.
        let mut reader = vcf::io::Reader::new(out.as_bytes());
        let header = reader.read_header().unwrap();
        let mut rec = RecordBuf::default();
        let _ = reader.read_record_buf(&header, &mut rec).unwrap();
        // Record 1: ALT=T, AF=[0.3], AD=[50,30]
        assert_eq!(rec.alternate_bases().as_ref(), &[String::from("T")]);
        match rec.info().get("AF").flatten() {
            Some(InfoValue::Array(InfoArray::Float(v))) => {
                assert_eq!(v, &vec![Some(0.3_f32)]);
            }
            other => panic!("expected AF float array, got {other:?}"),
        }
        match rec.info().get("AD").flatten() {
            Some(InfoValue::Array(InfoArray::Integer(v))) => {
                assert_eq!(v, &vec![Some(50), Some(30)]);
            }
            other => panic!("expected AD int array, got {other:?}"),
        }

        let _ = reader.read_record_buf(&header, &mut rec).unwrap();
        // Record 2: ALT=G, AF=[0.2], AD=[50,20]
        assert_eq!(rec.alternate_bases().as_ref(), &[String::from("G")]);
        match rec.info().get("AF").flatten() {
            Some(InfoValue::Array(InfoArray::Float(v))) => {
                assert_eq!(v, &vec![Some(0.2_f32)]);
            }
            other => panic!("expected AF float array, got {other:?}"),
        }
        match rec.info().get("AD").flatten() {
            Some(InfoValue::Array(InfoArray::Integer(v))) => {
                assert_eq!(v, &vec![Some(50), Some(20)]);
            }
            other => panic!("expected AD int array, got {other:?}"),
        }
    }

    // ── SNP pass-through ──────────────────────────────────────────────────────

    #[test]
    fn snps_pass_through_unchanged() {
        let input = std::fs::read(corpus_dir().join("basic.vcf")).unwrap();
        let opts = NormalizeOptions {
            split_multiallelics: true,
            left_align: true,
            check_ref: RefCheck::Ignore,
            output_format: crate::io::OutputFormat::Vcf,
        };
        let (_out, stats) = normalize_to_string(&input, opts);
        assert_eq!(stats.input_records, 5);
        assert_eq!(stats.output_records, 5);
        assert_eq!(stats.split_sites, 0);
        assert_eq!(stats.left_aligned, 0);
    }

    // ── symbolic ALT pass-through ─────────────────────────────────────────────

    #[test]
    fn symbolic_alts_pass_through() {
        let input = std::fs::read(corpus_dir().join("empty_alt.vcf")).unwrap();
        let opts = NormalizeOptions {
            split_multiallelics: true,
            left_align: true,
            check_ref: RefCheck::Ignore,
            output_format: crate::io::OutputFormat::Vcf,
        };
        let (_out, stats) = normalize_to_string(&input, opts);
        assert_eq!(stats.input_records, stats.output_records);
        assert_eq!(stats.split_sites, 0);
        assert_eq!(stats.left_aligned, 0);
    }

    // ── g-index math ──────────────────────────────────────────────────────────

    #[test]
    fn diploid_g_indices_triallelic() {
        // Alleles [REF, A, B] -> split on A (k=1) keeps {0/0, 0/1, 1/1} = indices 0,1,2
        assert_eq!(diploid_g_indices_to_keep(2, 1), vec![0, 1, 2]);
        // Split on B (k=2) keeps {0/0, 0/2, 2/2} = indices 0,3,5
        assert_eq!(diploid_g_indices_to_keep(2, 2), vec![0, 3, 5]);
    }
}
