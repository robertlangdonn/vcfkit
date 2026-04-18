//! Integration tests for `vcfkit liftover`.
//!
//! Every test works against inline chain and FASTA fixtures written to the
//! process's temp directory, so no real genome files are required. Tests that
//! exercise the installed `vcfkit` binary resolve it via `CARGO_BIN_EXE_vcfkit`.

mod common;

use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::atomic::{AtomicU64, Ordering},
};

use vcfkit_core::{
    io::OutputFormat,
    liftover::{liftover, ChainIndex, LiftoverOptions, LiftoverStats},
};

use crate::common::diff::parse_vcf_records;

// ── test infrastructure ──────────────────────────────────────────────────────

static UNIQ: AtomicU64 = AtomicU64::new(0);

/// Create a unique scratch directory under the OS temp dir and return it.
fn scratch_dir(tag: &str) -> PathBuf {
    let seq = UNIQ.fetch_add(1, Ordering::SeqCst);
    let pid = std::process::id();
    let dir = std::env::temp_dir().join(format!("vcfkit-liftover-{tag}-{pid}-{seq}"));
    fs::create_dir_all(&dir).expect("create scratch dir");
    dir
}

/// Write a file with the given `name` and `contents` inside `dir`, returning
/// the absolute path.
fn write_file(dir: &Path, name: &str, contents: &[u8]) -> PathBuf {
    let p = dir.join(name);
    fs::write(&p, contents).unwrap_or_else(|e| panic!("write {name}: {e}"));
    p
}

/// Return the path to the built `vcfkit` binary.
fn vcfkit_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_vcfkit"))
}

/// Write a tiny FASTA + FAI pair representing one or more contigs.  Sequences
/// are padded (if necessary) up to the next 60-bp boundary with 'N' so each
/// FAI entry describes a clean `linebases=60`/`linewidth=61` record.  The
/// returned path points at the `.fa`; the `.fa.fai` sidecar sits next to it.
fn write_mini_fasta(dir: &Path, name: &str, contigs: &[(&str, &[u8])]) -> PathBuf {
    let mut fasta = Vec::new();
    let mut fai_entries: Vec<String> = Vec::new();

    for (contig_name, seq) in contigs {
        // Pad the sequence to a multiple of 60 bp with 'N'.
        let line_blen = 60usize;
        let line_width = line_blen + 1;
        let padded_len = seq.len().div_ceil(line_blen) * line_blen;
        let mut padded = Vec::with_capacity(padded_len);
        padded.extend_from_slice(seq);
        padded.resize(padded_len, b'N');

        let header_line = format!(">{contig_name}\n");
        let offset = fasta.len() + header_line.len();
        fasta.extend_from_slice(header_line.as_bytes());
        for chunk in padded.chunks(line_blen) {
            fasta.extend_from_slice(chunk);
            fasta.push(b'\n');
        }

        fai_entries.push(format!(
            "{contig_name}\t{}\t{}\t{}\t{}",
            padded_len, offset, line_blen, line_width
        ));
    }

    let fasta_path = write_file(dir, name, &fasta);
    let fai_path_str = format!("{}.fai", fasta_path.display());
    let fai_path = PathBuf::from(&fai_path_str);
    fs::write(&fai_path, fai_entries.join("\n") + "\n").expect("write fai");
    fasta_path
}

/// Build a 60-bp sequence cycling A, C, G, T, ...  Deterministic helper.
fn cycle_acgt(len: usize) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ACGT";
    (0..len).map(|i| ALPHABET[i % 4]).collect()
}

/// Reverse-complement a byte sequence.
fn rc(bytes: &[u8]) -> Vec<u8> {
    bytes
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => *other,
        })
        .collect()
}

/// A minimal VCF body given a list of (chrom, pos, ref, alt) tuples.
fn make_vcf(records: &[(&str, u64, &str, &str)], contig_headers: &[(&str, usize)]) -> Vec<u8> {
    let mut v = String::new();
    v.push_str("##fileformat=VCFv4.2\n");
    v.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    for (c, len) in contig_headers {
        v.push_str(&format!("##contig=<ID={c},length={len}>\n"));
    }
    v.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    v.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n");
    for (chrom, pos, r, a) in records {
        v.push_str(&format!(
            "{chrom}\t{pos}\t.\t{r}\t{a}\t50\tPASS\t.\tGT\t0/1\n"
        ));
    }
    v.into_bytes()
}

/// Helper that runs the core `liftover` fn and returns (stdout, stats).
fn run_core(
    input: &[u8],
    chain_path: &Path,
    source_ref: &Path,
    target_ref: &Path,
    options: LiftoverOptions,
) -> (String, LiftoverStats) {
    let mut out = Vec::new();
    let stats = liftover(input, &mut out, chain_path, source_ref, target_ref, options)
        .expect("liftover ok");
    (String::from_utf8(out).expect("utf-8"), stats)
}

// ── 1. same-strand contig-preserving liftover ────────────────────────────────

#[test]
fn same_strand_same_contig_shifts_pos() {
    // chr1:1-50 → chr1:101-150 (same strand). FASTA contigs are 200 bp each,
    // enough to cover both ranges. Source/target both have the same cycled
    // sequence so that lifted REF matches.
    let dir = scratch_dir("same-strand");
    let seq = cycle_acgt(200);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &seq)]);

    // For a chain-only mapping chr1:0-50 → chr1:100-150 the target positions
    // don't contain the source sequence unless the two FASTAs agree at those
    // offsets. Since our cycled alphabet is period-4, pos_0 and pos_0+100
    // carry the same base. Good.
    let chain = "chain 1000 chr1 200 + 0 50 chr1 200 + 100 150 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    // Record at chr1:10 (1-based, pos0=9). REF base at pos0=9 is seq[9]
    // where seq = A,C,G,T,A,C,G,T,A,C,... → seq[9] = 'C'.
    assert_eq!(seq[9], b'C');
    let input = make_vcf(&[("chr1", 10, "C", "A")], &[("chr1", 200)]);
    let (out, stats) = run_core(
        &input,
        &chain_path,
        &src_fa,
        &tgt_fa,
        LiftoverOptions::default(),
    );

    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].chrom, "chr1");
    assert_eq!(records[0].pos, 110, "offset 9 -> target pos 100+9+1 = 110");
    assert_eq!(records[0].ref_allele, "C");
    assert_eq!(records[0].alt_alleles, vec!["A".to_string()]);
    assert_eq!(stats.input_records, 1);
    assert_eq!(stats.output_records, 1);
    assert_eq!(stats.rejected_unmapped, 0);
}

// ── 2. contig rename ─────────────────────────────────────────────────────────

#[test]
fn contig_rename_chr1_to_chr2() {
    let dir = scratch_dir("rename");
    let seq = cycle_acgt(260);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq[..200])]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr2", &seq)]);

    // chr1:0-50 → chr2:200-250.  Because the cycled alphabet is period-4, REF
    // at offset 0 in chr1 == REF at offset 200 in chr2 ('A').
    let chain = "chain 1000 chr1 200 + 0 50 chr2 260 + 200 250 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    // chr1:5 (pos0=4) → chr2: 200+4+1 = 205. REF at pos0=4 is seq[4] = 'A'.
    let input = make_vcf(&[("chr1", 5, "A", "G")], &[("chr1", 200), ("chr2", 260)]);
    let (out, stats) = run_core(
        &input,
        &chain_path,
        &src_fa,
        &tgt_fa,
        LiftoverOptions::default(),
    );

    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].chrom, "chr2");
    assert_eq!(records[0].pos, 205);
    assert_eq!(records[0].ref_allele, "A");
    assert_eq!(stats.output_records, 1);
}

// ── 3. negative-strand block reverse-complements REF/ALT ─────────────────────

#[test]
fn negative_strand_reverse_complements_alleles() {
    // Pick a target sequence whose reverse-complement at the mapped position
    // equals the source REF.
    //
    // Let source seq (chr1, 200 bp) = cycle ACGT.
    //   source pos0=10 → base 'G'.
    // Chain: chr1:0-50 → chr1:50-100 on '-' strand, tgt_size=200.
    //
    // Per the -strand normalization, forward target range is:
    //   [tgt_size - tgt_end, tgt_size - tgt_start) = [200-100, 200-50) = [100, 150)
    // And for src offset 10 with ref_len=1, lifted forward pos0 =
    //   tgt_end - offset - ref_len = 150 - 10 - 1 = 139 → 1-based 140.
    //
    // The lifted REF must equal reverse_complement(source_ref_base), which is
    // 'C'. So tgt seq at pos0=139 must be 'C'. We construct tgt accordingly.
    let dir = scratch_dir("neg");
    let src_seq = cycle_acgt(200);
    assert_eq!(src_seq[10], b'G');
    let mut tgt_seq = cycle_acgt(200);
    tgt_seq[139] = b'C';

    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &src_seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &tgt_seq)]);

    let chain = "chain 1000 chr1 200 + 0 50 chr1 200 - 50 150 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    // Source: chr1:11 REF=G ALT=A.  After flip: REF=C, ALT=T.
    let input = make_vcf(&[("chr1", 11, "G", "A")], &[("chr1", 200)]);
    let (out, stats) = run_core(
        &input,
        &chain_path,
        &src_fa,
        &tgt_fa,
        LiftoverOptions::default(),
    );

    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].pos, 140);
    assert_eq!(records[0].ref_allele, "C");
    assert_eq!(records[0].alt_alleles, vec!["T".to_string()]);
    assert_eq!(stats.swapped_alleles, 1);
    assert_eq!(stats.output_records, 1);
    // rc round-trip sanity check
    assert_eq!(rc(b"G"), b"C");
}

// ── 4. unmapped variant routes to reject file ────────────────────────────────

#[test]
fn unmapped_variant_goes_to_reject_file() {
    let dir = scratch_dir("reject");
    let seq = cycle_acgt(200);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &seq)]);

    // Chain covers only chr1:0-50 (source).  A variant at chr1:100 (pos0=99)
    // falls outside the chain and must be rejected.
    let chain = "chain 1000 chr1 200 + 0 50 chr1 200 + 0 50 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());
    let reject_path = dir.join("rejects.vcf");

    let input = make_vcf(&[("chr1", 100, "A", "G")], &[("chr1", 200)]);
    let opts = LiftoverOptions {
        reject_file: Some(reject_path.clone()),
        ..Default::default()
    };
    let (out, stats) = run_core(&input, &chain_path, &src_fa, &tgt_fa, opts);

    assert_eq!(stats.input_records, 1);
    assert_eq!(stats.output_records, 0);
    assert_eq!(stats.rejected_unmapped, 1);

    let main_records = parse_vcf_records(&out);
    assert!(main_records.is_empty());

    let reject_text = fs::read_to_string(&reject_path).expect("read reject file");
    let reject_records = parse_vcf_records(&reject_text);
    assert_eq!(reject_records.len(), 1);
    assert_eq!(reject_records[0].chrom, "chr1");
    assert_eq!(reject_records[0].pos, 100);
}

// ── 5. --write-src-coords adds INFO fields ───────────────────────────────────

#[test]
fn write_src_coords_populates_info_fields() {
    let dir = scratch_dir("src-coords");
    let seq = cycle_acgt(200);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &seq)]);

    let chain = "chain 1000 chr1 200 + 0 50 chr1 200 + 100 150 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    // chr1:5 -> pos 105. REF at pos0=4 is 'A'.
    let input = make_vcf(&[("chr1", 5, "A", "T")], &[("chr1", 200)]);
    let opts = LiftoverOptions {
        write_src_coords: true,
        ..Default::default()
    };
    let (out, stats) = run_core(&input, &chain_path, &src_fa, &tgt_fa, opts);
    assert_eq!(stats.output_records, 1);

    // Header must have SRC_CONTIG and SRC_POS lines.
    assert!(out.contains("##INFO=<ID=SRC_CONTIG"), "header: {out}");
    assert!(out.contains("##INFO=<ID=SRC_POS"), "header: {out}");

    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 1);
    assert_eq!(
        records[0].info.get("SRC_CONTIG").map(|s| s.as_str()),
        Some("chr1")
    );
    assert_eq!(
        records[0].info.get("SRC_POS").map(|s| s.as_str()),
        Some("5")
    );
    assert_eq!(records[0].pos, 105);
}

// ── 6. --list-chains prints URLs and exits 0 ─────────────────────────────────

#[test]
fn list_chains_prints_known_urls_and_exits_0() {
    let out = Command::new(vcfkit_bin())
        .args(["liftover", "--list-chains"])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("spawn vcfkit");

    assert!(out.status.success(), "--list-chains should succeed");
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    assert!(stdout.contains("hg19ToHg38.over.chain.gz"));
    assert!(stdout.contains("hg38ToHg19.over.chain.gz"));
    assert!(stdout.contains("hg38ToHs1.over.chain.gz"));
    assert!(stdout.contains("hs1ToHg38.over.chain.gz"));

    let line_count = stdout.lines().filter(|l| !l.is_empty()).count();
    assert_eq!(line_count, 4, "expected 4 non-empty lines in --list-chains");
}

// ── 7. stats counts are accurate (mapped + unmapped mix) ─────────────────────

#[test]
fn stats_accurately_count_mapped_and_unmapped() {
    let dir = scratch_dir("stats");
    let seq = cycle_acgt(200);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &seq)]);

    let chain = "chain 1000 chr1 200 + 0 50 chr1 200 + 100 150 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    // 3 mappable (pos 1, 10, 25) + 2 unmapped (pos 80, 180).
    // Reference bases (cycled ACGT, 0-indexed):
    //   pos 1 → pos0 0 → seq[0] = 'A'
    //   pos 10 → pos0 9 → seq[9] = 'C'
    //   pos 25 → pos0 24 → seq[24] = 'A'
    assert_eq!(seq[0], b'A');
    assert_eq!(seq[9], b'C');
    assert_eq!(seq[24], b'A');
    let input = make_vcf(
        &[
            ("chr1", 1, "A", "C"),
            ("chr1", 10, "C", "A"),
            ("chr1", 25, "A", "G"),
            ("chr1", 80, "A", "G"),
            ("chr1", 180, "A", "G"),
        ],
        &[("chr1", 200)],
    );
    let (_out, stats) = run_core(
        &input,
        &chain_path,
        &src_fa,
        &tgt_fa,
        LiftoverOptions::default(),
    );
    assert_eq!(stats.input_records, 5);
    assert_eq!(stats.output_records, 3);
    assert_eq!(stats.rejected_unmapped, 2);
}

// ── 8. REF mismatch at target is rejected ─────────────────────────────────────

#[test]
fn ref_mismatch_at_target_is_rejected() {
    let dir = scratch_dir("ref-mismatch");
    let src_seq = cycle_acgt(200);
    let mut tgt_seq = cycle_acgt(200);
    // Corrupt a single base at pos0=110 in target so the lifted REF won't match.
    tgt_seq[110] = b'N';
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &src_seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &tgt_seq)]);

    let chain = "chain 1000 chr1 200 + 0 50 chr1 200 + 100 150 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    // pos 11 -> lifted pos 111 (pos0 110). source REF at pos0=10 is 'G'. Target
    // has 'N' there, so mismatch.
    assert_eq!(src_seq[10], b'G');
    let input = make_vcf(&[("chr1", 11, "G", "A")], &[("chr1", 200)]);
    let (_out, stats) = run_core(
        &input,
        &chain_path,
        &src_fa,
        &tgt_fa,
        LiftoverOptions::default(),
    );
    assert_eq!(stats.output_records, 0);
    assert_eq!(stats.rejected_ref_mismatch, 1);
}

// ── 9. no_fix_swapped_ref rejects strand-flipping records ────────────────────

#[test]
fn no_fix_swapped_ref_rejects_negative_strand() {
    let dir = scratch_dir("no-fix-swap");
    let seq = cycle_acgt(200);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &seq)]);

    let chain = "chain 1000 chr1 200 + 0 50 chr1 200 - 50 150 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    let input = make_vcf(&[("chr1", 11, "G", "A")], &[("chr1", 200)]);
    let opts = LiftoverOptions {
        fix_swapped_ref: false,
        ..Default::default()
    };
    let (out, stats) = run_core(&input, &chain_path, &src_fa, &tgt_fa, opts);
    assert_eq!(stats.output_records, 0);
    assert_eq!(stats.rejected_unmapped, 1);
    let recs = parse_vcf_records(&out);
    assert!(recs.is_empty());
}

// ── 10. unknown source contig is unmapped ────────────────────────────────────

#[test]
fn unknown_source_contig_is_unmapped() {
    let dir = scratch_dir("unknown-contig");
    let seq = cycle_acgt(200);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq), ("chrX", &seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &seq)]);

    let chain = "chain 1 chr1 200 + 0 50 chr1 200 + 100 150 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    let input = make_vcf(&[("chrX", 5, "A", "G")], &[("chrX", 200)]);
    let (_out, stats) = run_core(
        &input,
        &chain_path,
        &src_fa,
        &tgt_fa,
        LiftoverOptions::default(),
    );
    assert_eq!(stats.rejected_unmapped, 1);
    assert_eq!(stats.output_records, 0);
}

// ── 11. chain with a gap rejects variants in the gap ─────────────────────────

#[test]
fn variants_in_chain_gaps_are_rejected() {
    let dir = scratch_dir("gap");
    let seq = cycle_acgt(200);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &seq)]);

    // Two 10-bp blocks with a 5-bp source gap:
    //   src 0..10 → tgt 100..110
    //   src 15..25 → tgt 115..125
    let chain = "chain 1000 chr1 200 + 0 25 chr1 200 + 100 125 1\n10 5 5\n10\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    // pos 5 → pos0=4 → seq[4] = 'A'
    // pos 12 → pos0=11 → seq[11] = 'T' (but in the gap — won't be checked)
    // pos 20 → pos0=19 → seq[19] = 'T'
    assert_eq!(seq[4], b'A');
    assert_eq!(seq[19], b'T');
    let input = make_vcf(
        &[
            ("chr1", 5, "A", "G"),  // in block 1 -> pos 104
            ("chr1", 12, "T", "G"), // in the gap -> rejected
            ("chr1", 20, "T", "A"), // in block 2 -> pos 120
        ],
        &[("chr1", 200)],
    );
    let (out, stats) = run_core(
        &input,
        &chain_path,
        &src_fa,
        &tgt_fa,
        LiftoverOptions::default(),
    );
    // Two of those should map; one (pos=12) falls in the gap.
    assert_eq!(stats.output_records, 2);
    assert_eq!(stats.rejected_unmapped, 1);

    let recs = parse_vcf_records(&out);
    let positions: Vec<u64> = recs.iter().map(|r| r.pos).collect();
    // pos 5 (offset 4 in block 0..10) → tgt 100 + 4 + 1 = 105
    // pos 20 (offset 4 in block 15..25) → tgt 115 + 4 + 1 = 120
    assert_eq!(positions, vec![105, 120]);
}

// ── 12. ChainIndex::lookup rejects positions at exact block boundary ─────────

#[test]
fn chain_index_lookup_is_half_open() {
    let chain = "chain 1 chr1 100 + 0 10 chr1 100 + 0 10 1\n10\n";
    let idx = ChainIndex::from_reader(chain.as_bytes()).unwrap();
    assert!(idx.lookup("chr1", 0).is_some());
    assert!(idx.lookup("chr1", 9).is_some());
    assert!(idx.lookup("chr1", 10).is_none(), "end is exclusive");
}

// ── 13. CLI invocation: read stdin, write stdout ─────────────────────────────

#[test]
fn cli_reads_stdin_and_writes_stdout() {
    let dir = scratch_dir("cli-pipe");
    let seq = cycle_acgt(200);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &seq)]);
    let chain = "chain 1000 chr1 200 + 0 50 chr1 200 + 100 150 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    let input = make_vcf(&[("chr1", 5, "A", "T")], &[("chr1", 200)]);

    let mut child = Command::new(vcfkit_bin())
        .args([
            "liftover",
            "-s",
            src_fa.to_str().unwrap(),
            "-t",
            tgt_fa.to_str().unwrap(),
            "-c",
            chain_path.to_str().unwrap(),
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn vcfkit");
    child
        .stdin
        .as_mut()
        .expect("stdin")
        .write_all(&input)
        .expect("write stdin");
    let out = child.wait_with_output().expect("wait vcfkit");
    assert!(
        out.status.success(),
        "vcfkit liftover exited {:?}: {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    let records = parse_vcf_records(&stdout);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].pos, 105);
}

// ── 14. CLI: --write-src-coords adds the INFO ────────────────────────────────

#[test]
fn cli_write_src_coords_end_to_end() {
    let dir = scratch_dir("cli-src-coords");
    let seq = cycle_acgt(200);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &seq)]);
    let chain = "chain 1000 chr1 200 + 0 50 chr1 200 + 100 150 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    // pos 10 -> pos0=9 -> seq[9]='C'
    assert_eq!(seq[9], b'C');
    let input = make_vcf(&[("chr1", 10, "C", "A")], &[("chr1", 200)]);
    let out_path = dir.join("out.vcf");
    let input_path = write_file(&dir, "in.vcf", &input);

    let status = Command::new(vcfkit_bin())
        .args([
            "liftover",
            "-s",
            src_fa.to_str().unwrap(),
            "-t",
            tgt_fa.to_str().unwrap(),
            "-c",
            chain_path.to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
            "--write-src-coords",
            input_path.to_str().unwrap(),
        ])
        .status()
        .expect("spawn vcfkit");
    assert!(status.success());

    let text = fs::read_to_string(&out_path).expect("read output");
    assert!(text.contains("##INFO=<ID=SRC_CONTIG"));
    assert!(text.contains("##INFO=<ID=SRC_POS"));
    let recs = parse_vcf_records(&text);
    assert_eq!(recs.len(), 1);
    assert_eq!(recs[0].info.get("SRC_POS").map(|s| s.as_str()), Some("10"));
    assert_eq!(
        recs[0].info.get("SRC_CONTIG").map(|s| s.as_str()),
        Some("chr1")
    );
}

// ── 15. symbolic ALTs are passed through untouched ───────────────────────────

#[test]
fn symbolic_alts_pass_through_on_positive_strand() {
    let dir = scratch_dir("symbolic");
    let seq = cycle_acgt(200);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr1", &seq)]);
    let chain = "chain 1000 chr1 200 + 0 50 chr1 200 + 100 150 1\n50\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    // Need to add ALT declaration in the mini VCF. Just build it by hand.
    let vcf = b"##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=200>\n\
##ALT=<ID=DEL,Description=\"Deletion\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t10\t.\tT\t<DEL>\t50\tPASS\t.\tGT\t0/1\n";
    let (out, stats) = run_core(
        vcf,
        &chain_path,
        &src_fa,
        &tgt_fa,
        LiftoverOptions {
            output_format: OutputFormat::Vcf,
            ..Default::default()
        },
    );
    assert_eq!(stats.output_records, 1);
    let recs = parse_vcf_records(&out);
    assert_eq!(recs[0].alt_alleles, vec!["<DEL>".to_string()]);
    assert_eq!(recs[0].pos, 110);
}

// ── 16. multi-block `-` strand chain: tgt_cursor advances correctly ───────────
//
// Chain layout:
//   chain 1000 chr1 50 + 0 25 chr2 100 - 0 23 1
//   10 5 3
//   10
//
// Block 1 (size=10): src=[0,10)
//   tgt_cursor=0 → forward tgt_end = 100-0 = 100, tgt_start = 100-10 = 90
//   → block maps src[0,10) → fwd-tgt[90,100)
//
// After gap (dt=5, dq=3): src_cursor=15, tgt_cursor=13
//
// Block 2 (size=10): src=[15,25)
//   tgt_cursor=13 → forward tgt_end = 100-13 = 87, tgt_start = 87-10 = 77
//   → block maps src[15,25) → fwd-tgt[77,87)
//
// Variant at src_pos0=17 (1-based pos=18), offset=2 from block2.src_start=15:
//   lifted fwd pos0 = tgt_end - offset - ref_len = 87 - 2 - 1 = 84  → 1-based 85
//
// We omit real FASTA files so REF validation is skipped (Reference::open fails
// on non-existent paths, returning None via .ok()).

#[test]
fn multi_block_neg_strand_tgt_cursor_advances_correctly() {
    let dir = scratch_dir("neg-multiblock");

    // chain: two 10-bp blocks, gap (dt=5, dq=3) between them, '-' strand target
    let chain = "chain 1000 chr1 50 + 0 25 chr2 100 - 0 23 1\n10 5 3\n10\n";
    let chain_path = write_file(&dir, "x.chain", chain.as_bytes());

    // Build a minimal target FASTA with chr2 (100 bp) so REF validation can run.
    // Block 2 maps src[15,25) → fwd-tgt[77,87). Variant at src_pos0=17 (offset=2)
    // lifts to tgt_pos0=84. REF on source = rev_comp(tgt[84]).
    //
    // We set tgt_seq[84] = 'C' and pick source REF = rev_comp("C") = "G".
    let mut tgt_seq = cycle_acgt(120);
    tgt_seq[84] = b'C';
    let src_seq = cycle_acgt(50);
    let src_fa = write_mini_fasta(&dir, "src.fa", &[("chr1", &src_seq)]);
    let tgt_fa = write_mini_fasta(&dir, "tgt.fa", &[("chr2", &tgt_seq)]);

    // src pos 18 (1-based) = pos0 17, in block2 (src_start=15), offset=2.
    // Source REF = "G" (will be rev_comp'd to "C" to match tgt[84]).
    let input = make_vcf(&[("chr1", 18, "G", "A")], &[("chr1", 50)]);
    let (out, stats) = run_core(
        &input,
        &chain_path,
        &src_fa,
        &tgt_fa,
        LiftoverOptions::default(),
    );

    assert_eq!(stats.input_records, 1, "one record in");
    assert_eq!(stats.output_records, 1, "should be lifted (not rejected)");
    assert_eq!(
        stats.swapped_alleles, 1,
        "alleles should be swapped for '-' strand"
    );

    let recs = parse_vcf_records(&out);
    assert_eq!(recs.len(), 1);
    assert_eq!(recs[0].chrom, "chr2");
    // tgt_end - offset - ref_len + 1 (1-based) = 87 - 2 - 1 + 1 = 85
    assert_eq!(
        recs[0].pos, 85,
        "second block maps to lower forward-strand position"
    );
    assert_eq!(recs[0].ref_allele, "C", "rev_comp(G) = C");
    assert_eq!(
        recs[0].alt_alleles,
        vec!["T".to_string()],
        "rev_comp(A) = T"
    );
}
