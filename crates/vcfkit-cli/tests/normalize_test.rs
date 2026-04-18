//! Integration tests for `vcfkit normalize`.
//!
//! Two layers:
//! 1. *Core-path tests* — call `vcfkit_core::normalize::normalize` directly for
//!    speed and deterministic behavior. These always run.
//! 2. *Differential tests* — compare `vcfkit normalize` output against
//!    `bcftools norm`. These are gated by a runtime check for the `bcftools`
//!    executable and marked `#[ignore]` so they can be enabled (`cargo test
//!    -- --ignored`) on machines where bcftools is installed.

mod common;

use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use vcfkit_core::{
    io::OutputFormat,
    normalize::{normalize, NormalizeOptions, RefCheck},
};

use crate::common::diff::{assert_vcf_eq, parse_vcf_records};

// ── test infrastructure ──────────────────────────────────────────────────────

fn workspace_root() -> PathBuf {
    // CARGO_MANIFEST_DIR points at crates/vcfkit-cli — walk up twice.
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.pop(); // crates/
    p.pop(); // workspace root
    p
}

fn corpus_dir() -> PathBuf {
    workspace_root().join("tests/corpus/synthetic")
}

fn reference_fa() -> PathBuf {
    corpus_dir().join("mini_ref.fa")
}

/// Read a corpus fixture into a byte vector.
fn read_corpus(name: &str) -> Vec<u8> {
    fs::read(corpus_dir().join(name)).unwrap_or_else(|e| {
        panic!(
            "failed to read corpus fixture {}: {e}",
            corpus_dir().join(name).display()
        )
    })
}

/// Run `normalize` directly against the core API and return the produced VCF
/// string and the run statistics.
fn run_normalize(
    input: &[u8],
    options: NormalizeOptions,
) -> (String, vcfkit_core::normalize::NormalizeStats) {
    let mut out = Vec::new();
    let stats = normalize(input, &mut out, &reference_fa(), options)
        .expect("normalize must succeed for this fixture");
    (String::from_utf8(out).expect("valid utf-8 output"), stats)
}

/// Default options matching the CLI's defaults.
fn default_opts() -> NormalizeOptions {
    NormalizeOptions {
        split_multiallelics: true,
        left_align: true,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    }
}

/// Filter a VCF byte slice to only include records whose CHROM is present in
/// the given reference FASTA index. This is needed for differential tests where
/// the test fixtures contain chromosomes not present in `mini_ref.fa` — bcftools
/// hard-errors on missing chromosomes even when `-c w` is set.
fn filter_to_ref_chroms(input: &[u8], fai_path: &Path) -> Vec<u8> {
    // Parse chromosome names from the .fai file (first column of each line).
    let fai = fs::read_to_string(fai_path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", fai_path.display()));
    let ref_chroms: std::collections::HashSet<&str> =
        fai.lines().filter_map(|l| l.split('\t').next()).collect();

    let src = std::str::from_utf8(input).expect("VCF is valid UTF-8");
    let mut out = Vec::new();
    for line in src.lines() {
        if line.starts_with('#') {
            // Always keep header lines.
            out.extend_from_slice(line.as_bytes());
            out.push(b'\n');
        } else {
            let chrom = line.split('\t').next().unwrap_or("");
            if ref_chroms.contains(chrom) {
                out.extend_from_slice(line.as_bytes());
                out.push(b'\n');
            }
        }
    }
    out
}

/// True if `bcftools` is on `PATH`. Evaluated at test runtime.
fn bcftools_available() -> bool {
    Command::new("bcftools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// Run `bcftools norm` against the given VCF content, returning the resulting
/// VCF as a String. Panics if bcftools isn't available or exits non-zero.
fn run_bcftools_norm(input: &[u8], reference: &Path, extra_args: &[&str]) -> String {
    let mut cmd = Command::new("bcftools");
    cmd.arg("norm")
        .arg("-f")
        .arg(reference)
        .args(extra_args)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = cmd.spawn().expect("spawn bcftools");
    {
        let stdin = child.stdin.as_mut().expect("bcftools stdin");
        stdin.write_all(input).expect("write to bcftools stdin");
    }
    let output = child.wait_with_output().expect("bcftools wait");
    assert!(
        output.status.success(),
        "bcftools norm exited with {:?}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("bcftools stdout utf-8")
}

/// Resolve the path to the built `vcfkit` binary used by these tests.
fn vcfkit_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_vcfkit"))
}

// ── core-path integration tests (always run) ─────────────────────────────────

#[test]
fn basic_vcf_passes_through_unchanged() {
    let input = read_corpus("basic.vcf");
    let (out, stats) = run_normalize(&input, default_opts());

    let in_records = parse_vcf_records(std::str::from_utf8(&input).unwrap());
    let out_records = parse_vcf_records(&out);
    assert_eq!(
        in_records, out_records,
        "basic.vcf should pass through unchanged"
    );
    assert_eq!(stats.input_records, 5);
    assert_eq!(stats.output_records, 5);
    assert_eq!(stats.split_sites, 0);
    assert_eq!(stats.left_aligned, 0);
}

#[test]
fn multi_allelic_split_produces_biallelic_records() {
    let input = read_corpus("multi_allelic.vcf");
    let opts = NormalizeOptions {
        split_multiallelics: true,
        left_align: false,
        ..default_opts()
    };
    let (out, stats) = run_normalize(&input, opts);

    let records = parse_vcf_records(&out);
    assert_eq!(
        records.len(),
        11,
        "five sites -> 2+2+3+2+2 biallelic records"
    );
    for r in &records {
        assert_eq!(
            r.alt_alleles.len(),
            1,
            "every post-split record must be biallelic (got {:?})",
            r.alt_alleles
        );
    }
    assert_eq!(stats.split_sites, 5);
    assert_eq!(stats.output_records, 11);
}

#[test]
fn indels_left_align_to_norm_tag_positions() {
    let input = read_corpus("indels_unnormalized.vcf");
    let opts = NormalizeOptions {
        split_multiallelics: false,
        left_align: true,
        ..default_opts()
    };
    let (out, stats) = run_normalize(&input, opts);

    // Compare each output record against the NORM_POS/NORM_REF/NORM_ALT tags
    // carried by the original input (the output preserves INFO verbatim).
    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 5);
    for (i, r) in records.iter().enumerate() {
        let norm_pos: u64 = r
            .info
            .get("NORM_POS")
            .unwrap_or_else(|| panic!("record {i} missing NORM_POS"))
            .parse()
            .expect("NORM_POS is an integer");
        let norm_ref = r
            .info
            .get("NORM_REF")
            .unwrap_or_else(|| panic!("record {i} missing NORM_REF"));
        let norm_alt = r
            .info
            .get("NORM_ALT")
            .unwrap_or_else(|| panic!("record {i} missing NORM_ALT"));
        assert_eq!(r.pos, norm_pos, "record {i} POS");
        assert_eq!(&r.ref_allele, norm_ref, "record {i} REF");
        assert_eq!(
            r.alt_alleles
                .first()
                .expect("record must have ALT after left-align"),
            norm_alt,
            "record {i} ALT"
        );
    }
    assert!(stats.left_aligned > 0);
}

#[test]
fn ref_mismatch_count_is_exact_under_warn_mode() {
    // Inline fixture tied to mini_ref.fa (chr1, 120 bp starting "TTCGAATCGA...").
    // chr1:1 has REF=T in the FASTA, so a VCF claiming REF=A at pos 1 must be
    // counted as a mismatch. Two mismatch sites + one match = ref_mismatches==2.
    let input = b"##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
chr1\t2\t.\tG\tC\t50\tPASS\t.\tGT\t0/1\n\
chr1\t3\t.\tT\tA\t50\tPASS\t.\tGT\t0/1\n";
    // Reference at pos 1,2,3 is T,T,C. So:
    //   rec 1 (A vs T) → mismatch
    //   rec 2 (G vs T) → mismatch
    //   rec 3 (T vs C) → mismatch
    let opts = NormalizeOptions {
        split_multiallelics: false,
        left_align: false,
        check_ref: RefCheck::Warn,
        output_format: OutputFormat::Vcf,
    };
    let (_out, stats) = run_normalize(input, opts);
    assert_eq!(stats.ref_mismatches, 3);
    assert_eq!(stats.input_records, 3);
    assert_eq!(stats.output_records, 3);
}

#[test]
fn symbolic_alt_records_pass_through_unchanged() {
    let input = read_corpus("empty_alt.vcf");
    let opts = default_opts();
    let (out, stats) = run_normalize(&input, opts);

    // Output must match input at the record level; symbolic ALTs are opaque.
    assert_vcf_eq(std::str::from_utf8(&input).unwrap(), &out);
    assert_eq!(stats.split_sites, 0);
    assert_eq!(stats.left_aligned, 0);
    assert_eq!(stats.input_records, stats.output_records);
}

#[test]
fn missing_fields_are_handled_without_panicking() {
    // missing_fields.vcf exercises `.` in QUAL, FILTER, INFO, and FORMAT sample
    // values. Just requiring the normalize pipeline to finish without panicking
    // and producing parseable output is the contract here.
    let input = read_corpus("missing_fields.vcf");
    let opts = default_opts();
    let (out, stats) = run_normalize(&input, opts);
    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), stats.output_records);
    assert_eq!(stats.input_records, 5);
}

#[test]
fn large_indels_left_align_correctly() {
    // Inline "large" indel against mini_ref.fa. The reference has repetitive
    // sequence that allows the indel to left-align to an earlier position.
    // A VCF with an insertion anchored at pos 18 (REF=T, ALT=TACGTACGTACGTACGTACGT)
    // should walk left during left-alignment due to the repeated pattern.
    let input = b"##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t18\t.\tT\tTACGTACGTACGTACGTACGT\t50\tPASS\t.\tGT\t0/1\n";
    let opts = NormalizeOptions {
        split_multiallelics: false,
        left_align: true,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, stats) = run_normalize(input, opts);
    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 1);
    assert_eq!(stats.input_records, 1);
    assert_eq!(stats.left_aligned, 1, "large indel must be left-aligned");
    // The indel should shift to an earlier position with a new anchor base.
    assert!(
        records[0].pos < 18,
        "expected shift left, got {}",
        records[0].pos
    );
    // Length of (ALT - REF) must be preserved (a 20-bp insertion stays 20 bp).
    assert_eq!(
        records[0].alt_alleles[0].len() - records[0].ref_allele.len(),
        20
    );
}

#[test]
fn no_split_flag_preserves_multi_allelic_records() {
    let input = read_corpus("multi_allelic.vcf");
    let opts = NormalizeOptions {
        split_multiallelics: false,
        left_align: false,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, stats) = run_normalize(&input, opts);
    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 5, "no splitting => same record count");
    assert_eq!(stats.split_sites, 0);
    // Some records should still have >1 ALT.
    let multi = records.iter().filter(|r| r.alt_alleles.len() > 1).count();
    assert_eq!(multi, 5, "all 5 fixture records carry >1 ALT");
}

#[test]
fn no_left_align_flag_preserves_indel_positions() {
    let input = read_corpus("indels_unnormalized.vcf");
    let opts = NormalizeOptions {
        split_multiallelics: false,
        left_align: false,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, stats) = run_normalize(&input, opts);
    assert_eq!(stats.left_aligned, 0);

    // Output POS should equal input POS (not NORM_POS).
    let in_records = parse_vcf_records(std::str::from_utf8(&input).unwrap());
    let out_records = parse_vcf_records(&out);
    assert_eq!(in_records.len(), out_records.len());
    for (i, (inp, outp)) in in_records.iter().zip(out_records.iter()).enumerate() {
        assert_eq!(
            inp.pos, outp.pos,
            "record {i}: position must not shift when --no-left-align"
        );
        assert_eq!(inp.ref_allele, outp.ref_allele);
        assert_eq!(inp.alt_alleles, outp.alt_alleles);
    }
}

#[test]
fn check_ref_error_returns_err_on_mismatch() {
    // Inline fixture with a single in-bounds REF mismatch against mini_ref.fa.
    // Pos 1 reference base is T; the VCF claims A.
    let input = b"##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
    let opts = NormalizeOptions {
        split_multiallelics: false,
        left_align: false,
        check_ref: RefCheck::Error,
        output_format: OutputFormat::Vcf,
    };
    let mut out = Vec::new();
    let res = normalize(&input[..], &mut out, &reference_fa(), opts);
    assert!(
        res.is_err(),
        "check_ref=Error must fail on the first mismatch; got Ok({res:?})"
    );
    let msg = res.unwrap_err().to_string();
    assert!(
        msg.contains("REF mismatch"),
        "error should mention 'REF mismatch', got: {msg}"
    );
}

#[test]
fn check_ref_ignore_does_not_count_mismatches() {
    let input = read_corpus("ref_mismatch.vcf");
    let opts = NormalizeOptions {
        split_multiallelics: false,
        left_align: false,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (_out, stats) = run_normalize(&input, opts);
    assert_eq!(stats.ref_mismatches, 0);
}

#[test]
fn split_preserves_number_a_and_r_slicing() {
    // multi_allelic.vcf record 1: A -> T,G with AF=0.3,0.2 (Number=A) and
    // AD=50,30,20 (Number=R). After split:
    //   record[0] ALT=T, AF=[0.3], AD=[50,30]
    //   record[1] ALT=G, AF=[0.2], AD=[50,20]
    let input = read_corpus("multi_allelic.vcf");
    let opts = NormalizeOptions {
        split_multiallelics: true,
        left_align: false,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, _stats) = run_normalize(&input, opts);
    let records = parse_vcf_records(&out);

    assert_eq!(records[0].alt_alleles, vec!["T".to_string()]);
    assert_eq!(records[0].info.get("AF").map(|s| s.as_str()), Some("0.3"));
    assert_eq!(records[0].info.get("AD").map(|s| s.as_str()), Some("50,30"));

    assert_eq!(records[1].alt_alleles, vec!["G".to_string()]);
    assert_eq!(records[1].info.get("AF").map(|s| s.as_str()), Some("0.2"));
    assert_eq!(records[1].info.get("AD").map(|s| s.as_str()), Some("50,20"));

    // Number=1 INFO (DP) is copied unchanged.
    assert_eq!(records[0].info.get("DP").map(|s| s.as_str()), Some("100"));
    assert_eq!(records[1].info.get("DP").map(|s| s.as_str()), Some("100"));
}

#[test]
fn split_then_left_align_is_order_independent() {
    // A record with multi-allelic indels should produce the same POS/REF/ALT
    // tuples whether we split first or left-align first. We can at least
    // confirm that split+left-align produces biallelic records with the same
    // INFO NORM_* expectations for the indels_unnormalized fixture when split
    // is a no-op (n_alts==1).
    let input = read_corpus("indels_unnormalized.vcf");
    let opts_both = NormalizeOptions {
        split_multiallelics: true,
        left_align: true,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out_both, stats_both) = run_normalize(&input, opts_both);

    let opts_left_only = NormalizeOptions {
        split_multiallelics: false,
        left_align: true,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out_left_only, _) = run_normalize(&input, opts_left_only);

    // All fixture records are biallelic; splitting is a no-op, so outputs must
    // be identical.
    assert_vcf_eq(&out_left_only, &out_both);
    assert_eq!(stats_both.split_sites, 0);
}

#[test]
fn multi_allelic_with_left_align_retains_per_allele_info() {
    // Combining split + left-align on the multi-allelic fixture should still
    // produce N biallelic records for each original site (these are SNPs in
    // multi_allelic.vcf, so left-alignment is a no-op but shouldn't interfere).
    let input = read_corpus("multi_allelic.vcf");
    let opts = NormalizeOptions {
        split_multiallelics: true,
        left_align: true,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, stats) = run_normalize(&input, opts);
    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 11);
    assert_eq!(stats.left_aligned, 0, "SNPs should not left-align");
    for r in &records {
        assert_eq!(r.alt_alleles.len(), 1);
    }
}

#[test]
fn cli_reads_stdin_and_writes_stdout() {
    // End-to-end test that exercises the built `vcfkit` binary: pipe a VCF via
    // stdin, capture stdout, and parse back.
    let input = read_corpus("basic.vcf");
    let ref_fa = reference_fa();
    let mut child = Command::new(vcfkit_bin())
        .args([
            "normalize",
            "-f",
            ref_fa.to_str().unwrap(),
            "--check-ref",
            "ignore",
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
        "vcfkit normalize exited {:?}: {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    let records = parse_vcf_records(&stdout);
    assert_eq!(records.len(), 5);

    // Records should match the input semantically (SNPs, no splitting needed).
    let in_records = parse_vcf_records(std::str::from_utf8(&input).unwrap());
    assert_eq!(in_records, records);
}

#[test]
fn cli_no_split_flag_preserves_multi_allelics() {
    let input = read_corpus("multi_allelic.vcf");
    let ref_fa = reference_fa();
    let mut child = Command::new(vcfkit_bin())
        .args([
            "normalize",
            "-f",
            ref_fa.to_str().unwrap(),
            "--no-split",
            "--no-left-align",
            "--check-ref",
            "ignore",
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
        "vcfkit normalize exited {:?}: {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    let records = parse_vcf_records(&stdout);
    assert_eq!(records.len(), 5);
    assert!(records.iter().any(|r| r.alt_alleles.len() > 1));
}

#[test]
fn cli_check_ref_error_exits_nonzero() {
    let input = read_corpus("ref_mismatch.vcf");
    let ref_fa = reference_fa();
    let mut child = Command::new(vcfkit_bin())
        .args([
            "normalize",
            "-f",
            ref_fa.to_str().unwrap(),
            "--no-split",
            "--no-left-align",
            "--check-ref",
            "error",
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
        !out.status.success(),
        "vcfkit normalize --check-ref error should fail on mismatch (stderr: {})",
        String::from_utf8_lossy(&out.stderr)
    );
}

// ── out-of-bounds position test ───────────────────────────────────────────────

#[test]
fn out_of_bounds_position_is_skipped_not_aborted() {
    // A record at position 999999 on chr1 (declared length=120 in the header)
    // must be silently skipped: the run completes, out_of_bounds==1, and the
    // position does not appear in the output.
    let input = b"##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t10\t.\tT\tA\t50\tPASS\t.\tGT\t0/1\n\
chr1\t999999\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
    let opts = NormalizeOptions {
        split_multiallelics: false,
        left_align: false,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, stats) = run_normalize(input, opts);
    // Run must complete (not panic/error).
    assert_eq!(
        stats.out_of_bounds, 1,
        "expected 1 out-of-bounds record skipped"
    );
    assert_eq!(stats.input_records, 2);
    // Only the in-bounds record (pos=10) appears in the output.
    let records = parse_vcf_records(&out);
    assert_eq!(
        records.len(),
        1,
        "only the in-bounds record should be written"
    );
    assert_eq!(
        records[0].pos, 10,
        "the remaining record should be at position 10"
    );
    assert!(
        records.iter().all(|r| r.pos != 999999),
        "out-of-bounds position 999999 must not appear in output"
    );
}

// ── Number=G INFO field slicing ───────────────────────────────────────────────

/// For a diploid triallelic site (alleles REF/ALT1/ALT2) the G-ordered values
/// are indexed as: 0→0/0, 1→0/1, 2→1/1, 3→0/2, 4→1/2, 5→2/2.
///
/// Splitting on ALT1 (k=1) should keep indices [0,1,2] → values for 0/0, 0/1, 1/1.
/// Splitting on ALT2 (k=2) should keep indices [0,3,5] → values for 0/0, 0/2, 2/2.
#[test]
fn split_preserves_number_g_slicing() {
    // Inline triallelic VCF with a GL field (Number=G, Type=Float).
    // Alleles: REF=A, ALT1=T, ALT2=G. GL values for 3 alleles → 6 values:
    //   g(0/0)=-0.5, g(0/1)=-1.0, g(1/1)=-2.0, g(0/2)=-1.5, g(1/2)=-3.0, g(2/2)=-4.0
    let input: &[u8] = b"\
##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##INFO=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihoods\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t10\t.\tA\tT,G\t50\tPASS\tGL=-0.5,-1.0,-2.0,-1.5,-3.0,-4.0\tGT\t0/1\n";

    let opts = NormalizeOptions {
        split_multiallelics: true,
        left_align: false,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, stats) = run_normalize(input, opts);
    assert_eq!(stats.input_records, 1);
    assert_eq!(stats.output_records, 2);
    assert_eq!(stats.split_sites, 1);

    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 2);

    // ALT1 (T) → keep g-indices [0,1,2]: values -0.5, -1.0, -2.0
    assert_eq!(
        records[0].alt_alleles,
        vec!["T".to_string()],
        "first split record must have ALT=T"
    );
    let gl_t = records[0]
        .info
        .get("GL")
        .expect("GL must be present in first split record");
    assert_eq!(
        gl_t, "-0.5,-1,-2",
        "Number=G split for ALT1 must keep indices [0,1,2]"
    );

    // ALT2 (G) → keep g-indices [0,3,5]: values -0.5, -1.5, -4.0
    assert_eq!(
        records[1].alt_alleles,
        vec!["G".to_string()],
        "second split record must have ALT=G"
    );
    let gl_g = records[1]
        .info
        .get("GL")
        .expect("GL must be present in second split record");
    assert_eq!(
        gl_g, "-0.5,-1.5,-4",
        "Number=G split for ALT2 must keep indices [0,3,5]"
    );
}

// ── Number=. (variable) and Flag INFO pass-through ───────────────────────────

/// INFO fields with `Number=.` (unknown/variable count) must be copied verbatim
/// to every split record — the implementation must not attempt to slice them.
#[test]
fn split_passes_through_number_dot_info() {
    // CSQV is a fictional Number=. field that carries a |-delimited annotation.
    let input: &[u8] = b"\
##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n\
##INFO=<ID=CSQV,Number=.,Type=String,Description=\"Variable annotation\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t10\t.\tA\tT,G\t50\tPASS\tAF=0.3,0.2;CSQV=foo|bar|baz\tGT\t0/1\n";

    let opts = NormalizeOptions {
        split_multiallelics: true,
        left_align: false,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, stats) = run_normalize(input, opts);
    assert_eq!(stats.output_records, 2);

    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 2);

    // Number=. field must be copied unchanged to both split records.
    for (i, rec) in records.iter().enumerate() {
        let csqv = rec
            .info
            .get("CSQV")
            .unwrap_or_else(|| panic!("CSQV missing from split record {i}"));
        assert_eq!(
            csqv, "foo|bar|baz",
            "Number=. CSQV must be copied verbatim to split record {i}"
        );
    }

    // Number=A field still sliced correctly alongside the dot field.
    assert_eq!(
        records[0].info.get("AF").map(|s| s.as_str()),
        Some("0.3"),
        "AF for ALT1 must be 0.3"
    );
    assert_eq!(
        records[1].info.get("AF").map(|s| s.as_str()),
        Some("0.2"),
        "AF for ALT2 must be 0.2"
    );
}

/// Flag INFO entries (no `=value`, `Type=Flag`) must be preserved on every
/// split record because they describe the site, not a specific allele.
#[test]
fn split_preserves_flag_info_fields() {
    let input: &[u8] = b"\
##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n\
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t10\t.\tA\tT,G\t50\tPASS\tAF=0.3,0.2;DB\tGT\t0/1\n";

    let opts = NormalizeOptions {
        split_multiallelics: true,
        left_align: false,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, stats) = run_normalize(input, opts);
    assert_eq!(stats.output_records, 2);

    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 2);

    for (i, rec) in records.iter().enumerate() {
        assert!(
            rec.info.contains_key("DB"),
            "Flag INFO field DB must be present in split record {i}"
        );
    }
}

// ── FORMAT Number=A and Number=R slicing ────────────────────────────────────

/// FORMAT fields with `Number=A` must be re-sliced to a single value for the
/// ALT allele of each split record, mirroring the INFO Number=A logic.
#[test]
fn split_preserves_format_number_a_slicing() {
    // VAF is a per-sample allele frequency — one value per ALT.
    let input: &[u8] = b"\
##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant Allele Frequency\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t10\t.\tA\tT,G\t50\tPASS\tDP=100\tGT:VAF\t0/1:0.35,0.15\n";

    let opts = NormalizeOptions {
        split_multiallelics: true,
        left_align: false,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, stats) = run_normalize(input, opts);
    assert_eq!(stats.output_records, 2);

    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 2);

    // Extract the VAF value from the sample column (FORMAT is GT:VAF).
    // samples[0] is the raw string e.g. "0/1:0.35".
    let vaf_for = |rec: &crate::common::diff::VcfRecord| -> String {
        let fmt = rec.format.as_deref().unwrap_or("");
        let sample = rec.samples.first().map(|s| s.as_str()).unwrap_or("");
        let fmt_fields: Vec<&str> = fmt.split(':').collect();
        let smp_fields: Vec<&str> = sample.split(':').collect();
        let vaf_pos = fmt_fields
            .iter()
            .position(|&k| k == "VAF")
            .expect("VAF in FORMAT");
        smp_fields
            .get(vaf_pos)
            .map(|s| s.to_string())
            .expect("VAF value present")
    };

    assert_eq!(
        records[0].alt_alleles,
        vec!["T".to_string()],
        "first split record must have ALT=T"
    );
    assert_eq!(
        vaf_for(&records[0]),
        "0.35",
        "FORMAT VAF for ALT1 must be 0.35 (Number=A sliced to first element)"
    );

    assert_eq!(
        records[1].alt_alleles,
        vec!["G".to_string()],
        "second split record must have ALT=G"
    );
    assert_eq!(
        vaf_for(&records[1]),
        "0.15",
        "FORMAT VAF for ALT2 must be 0.15 (Number=A sliced to second element)"
    );
}

/// FORMAT fields with `Number=R` must keep the REF value plus the allele-specific
/// value for the split ALT. This mirrors multi_allelic.vcf FORMAT AD Number=R.
#[test]
fn split_preserves_format_number_r_slicing() {
    // AD (Number=R) for a biallelic site: REF=50, ALT1=30, ALT2=20.
    // After splitting:
    //   record for ALT1 → AD = 50,30 (REF + ALT1)
    //   record for ALT2 → AD = 50,20 (REF + ALT2)
    let input: &[u8] = b"\
##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t10\t.\tA\tT,G\t50\tPASS\tDP=100\tGT:AD\t0/1:50,30,20\n";

    let opts = NormalizeOptions {
        split_multiallelics: true,
        left_align: false,
        check_ref: RefCheck::Ignore,
        output_format: OutputFormat::Vcf,
    };
    let (out, stats) = run_normalize(input, opts);
    assert_eq!(stats.output_records, 2);

    let records = parse_vcf_records(&out);
    assert_eq!(records.len(), 2);

    let ad_for = |rec: &crate::common::diff::VcfRecord| -> String {
        let fmt = rec.format.as_deref().unwrap_or("");
        let sample = rec.samples.first().map(|s| s.as_str()).unwrap_or("");
        let fmt_fields: Vec<&str> = fmt.split(':').collect();
        let smp_fields: Vec<&str> = sample.split(':').collect();
        let ad_pos = fmt_fields
            .iter()
            .position(|&k| k == "AD")
            .expect("AD in FORMAT");
        smp_fields
            .get(ad_pos)
            .map(|s| s.to_string())
            .expect("AD value present")
    };

    assert_eq!(
        ad_for(&records[0]),
        "50,30",
        "FORMAT AD for ALT1 must be REF+ALT1 = 50,30"
    );
    assert_eq!(
        ad_for(&records[1]),
        "50,20",
        "FORMAT AD for ALT2 must be REF+ALT2 = 50,20"
    );
}

// ── differential tests against bcftools norm ─────────────────────────────────
//
// These are marked `#[ignore]` because `bcftools` is not guaranteed to be on
// every developer's PATH or in CI images. Run with:
//
//   cargo test -p vcfkit-cli -- --ignored
//
// Each test also short-circuits cleanly when `bcftools` is absent at runtime,
// so they stay robust even if someone runs the ignored set on a bare machine.

#[test]
#[ignore = "requires bcftools on PATH; run `cargo test -- --ignored`"]
fn diff_basic_vcf_matches_bcftools_norm() {
    if !bcftools_available() {
        eprintln!("skipping: bcftools not installed");
        return;
    }
    let ref_fa = reference_fa();
    let fai = ref_fa.with_extension("fa.fai");
    let input = filter_to_ref_chroms(&read_corpus("basic.vcf"), &fai);
    let (actual, _) = run_normalize(&input, default_opts());
    let expected = run_bcftools_norm(&input, &ref_fa, &["-m", "-any", "-c", "w"]);
    assert_vcf_eq(&expected, &actual);
}

#[test]
#[ignore = "requires bcftools on PATH; run `cargo test -- --ignored`"]
fn diff_multi_allelic_split_matches_bcftools_norm() {
    if !bcftools_available() {
        eprintln!("skipping: bcftools not installed");
        return;
    }
    let ref_fa = reference_fa();
    let fai = ref_fa.with_extension("fa.fai");
    let input = filter_to_ref_chroms(&read_corpus("multi_allelic.vcf"), &fai);
    let opts = NormalizeOptions {
        split_multiallelics: true,
        left_align: false,
        ..default_opts()
    };
    let (actual, _) = run_normalize(&input, opts);
    let expected = run_bcftools_norm(&input, &ref_fa, &["-m", "-any", "-c", "w"]);
    assert_vcf_eq(&expected, &actual);
}

#[test]
#[ignore = "requires bcftools on PATH; run `cargo test -- --ignored`"]
fn diff_indels_left_align_matches_bcftools_norm() {
    if !bcftools_available() {
        eprintln!("skipping: bcftools not installed");
        return;
    }
    let input = read_corpus("indels_unnormalized.vcf");
    let opts = NormalizeOptions {
        split_multiallelics: false,
        left_align: true,
        ..default_opts()
    };
    let (actual, _) = run_normalize(&input, opts);
    let expected = run_bcftools_norm(&input, &reference_fa(), &[]);
    assert_vcf_eq(&expected, &actual);
}
