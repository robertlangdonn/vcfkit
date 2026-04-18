//! Integration tests for `vcfkit filter` and the underlying
//! `vcfkit_core::filter` API.

mod common;

use std::{
    fs,
    io::Write,
    path::PathBuf,
    process::{Command, Stdio},
};

use vcfkit_core::filter::{FilterExpression, FilterOptions, filter};

use crate::common::diff::parse_vcf_records;

// ── test infrastructure ──────────────────────────────────────────────────────

fn workspace_root() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.pop();
    p.pop();
    p
}

fn corpus_dir() -> PathBuf {
    workspace_root().join("tests/corpus/synthetic")
}

fn read_corpus(name: &str) -> Vec<u8> {
    fs::read(corpus_dir().join(name)).unwrap_or_else(|e| {
        panic!(
            "failed to read corpus fixture {}: {e}",
            corpus_dir().join(name).display()
        )
    })
}

fn run_filter(input: &[u8], expr: &str, invert: bool) -> (String, vcfkit_core::filter::FilterStats) {
    let ast = FilterExpression::parse(expr)
        .unwrap_or_else(|e| panic!("expression {expr:?} should parse: {e}"));
    let opts = FilterOptions {
        invert,
        ..Default::default()
    };
    let mut out = Vec::new();
    let stats = filter(input, &mut out, ast, opts).expect("filter should succeed");
    (String::from_utf8(out).expect("utf-8"), stats)
}

fn vcfkit_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_vcfkit"))
}

// A small inline fixture with INFO fields for numeric tests.
const TYPED_VCF: &[u8] = b"##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FILTER=<ID=LowQual,Description=\"Low quality\">\n\
##contig=<ID=chr1,length=248956422>\n\
##contig=<ID=chr17,length=83257441>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n\
chr1\t100\t.\tA\tT\t55\tPASS\tAF=0.02;DP=40;CSQ=missense_variant|GENE1\tGT:DP\t0/1:40\n\
chr1\t200\t.\tC\tG\t18\tLowQual\tAF=0.30;DP=30;CSQ=synonymous_variant|GENE1\tGT:DP\t0/1:30\n\
chr17\t150\t.\tG\tA\t72\tPASS\tAF=0.003;DP=80;CSQ=missense_variant|BRCA1\tGT:DP\t1/1:80\n\
chr17\t250\t.\tT\tC\t25\tLowQual\tAF=0.45;DP=9;CSQ=intron_variant|BRCA1\tGT:DP\t0/1:9\n\
chr1\t300\t.\tA\tG\t90\tPASS\tAF=0.001;DP=100;CSQ=missense_variant|GENE2\tGT:DP\t1/1:100\n";

// ── core tests ───────────────────────────────────────────────────────────────

/// 1. `INFO/AF < 0.05` — keeps records below the threshold.
#[test]
fn af_below_threshold() {
    let (out, stats) = run_filter(TYPED_VCF, "INFO/AF < 0.05", false);
    let recs = parse_vcf_records(&out);
    assert_eq!(stats.input_records, 5);
    // AFs: 0.02 ✓, 0.30 ✗, 0.003 ✓, 0.45 ✗, 0.001 ✓
    assert_eq!(stats.output_records, 3);
    assert_eq!(recs.len(), 3);
    for r in &recs {
        let af: f64 = r.info["AF"].parse().expect("AF is a float");
        assert!(af < 0.05, "expected AF<0.05, got {af}");
    }
}

/// 2. `FILTER == 'PASS'` — exact string match.
#[test]
fn filter_equals_pass() {
    let (out, stats) = run_filter(TYPED_VCF, "FILTER == 'PASS'", false);
    let recs = parse_vcf_records(&out);
    assert_eq!(stats.output_records, 3);
    for r in &recs {
        assert_eq!(r.filter, vec!["PASS".to_string()]);
    }
}

/// 3. `INFO/CSQ ~ 'missense'` — substring match.
#[test]
fn csq_contains_missense() {
    let (out, stats) = run_filter(TYPED_VCF, "INFO/CSQ ~ 'missense'", false);
    let recs = parse_vcf_records(&out);
    assert_eq!(stats.output_records, 3);
    for r in &recs {
        assert!(
            r.info["CSQ"].contains("missense"),
            "CSQ did not contain 'missense': {}",
            r.info["CSQ"]
        );
    }
}

/// 4. `CHROM == 'chr1' && INFO/DP > 10` — AND combination.
#[test]
fn chrom_and_dp() {
    let (out, stats) = run_filter(TYPED_VCF, "CHROM == 'chr1' && INFO/DP > 10", false);
    let recs = parse_vcf_records(&out);
    // chr1 records: pos 100 DP=40 ✓, pos 200 DP=30 ✓, pos 300 DP=100 ✓
    assert_eq!(stats.output_records, 3);
    for r in &recs {
        assert_eq!(r.chrom, "chr1");
    }
}

/// 5. `INFO/AF < 0.01 || CHROM == 'chr17'` — OR combination.
#[test]
fn af_or_chrom() {
    let (out, stats) = run_filter(
        TYPED_VCF,
        "INFO/AF < 0.01 || CHROM == 'chr17'",
        false,
    );
    let recs = parse_vcf_records(&out);
    // AF<0.01: pos 300 (0.001) ✓, pos 17:150 (0.003) ✓. Plus all chr17 records.
    // chr17 records: pos 150, pos 250 → 2. pos 300 adds 1. Total unique: 3.
    assert_eq!(stats.output_records, 3);
    // Every record must be either chr17 or have AF<0.01.
    for r in &recs {
        let af: f64 = r.info["AF"].parse().unwrap();
        assert!(r.chrom == "chr17" || af < 0.01);
    }
}

/// 6. `!` (NOT) operator.
#[test]
fn not_operator() {
    let (out, stats) = run_filter(TYPED_VCF, "!(FILTER == 'PASS')", false);
    let recs = parse_vcf_records(&out);
    // 2 LowQual records
    assert_eq!(stats.output_records, 2);
    for r in &recs {
        assert_ne!(r.filter, vec!["PASS".to_string()]);
    }
}

/// 7. `--invert` flag inverts the filter.
#[test]
fn invert_flag_inverts_filter() {
    let (out_plain, _) = run_filter(TYPED_VCF, "FILTER == 'PASS'", false);
    let (out_inv, stats_inv) = run_filter(TYPED_VCF, "FILTER == 'PASS'", true);
    let plain = parse_vcf_records(&out_plain);
    let inverted = parse_vcf_records(&out_inv);

    assert_eq!(plain.len() + inverted.len(), 5);
    assert_eq!(stats_inv.filtered_out, 3);
    for r in &inverted {
        assert_ne!(r.filter, vec!["PASS".to_string()]);
    }
}

/// 8. Missing INFO field evaluates to false.
#[test]
fn missing_info_field_does_not_pass() {
    // Use missing_fields.vcf; records lack AF in most rows, so
    // `INFO/AF < 1.0` should only match rows where AF is present and < 1.0.
    let input = read_corpus("missing_fields.vcf");
    let (out, _stats) = run_filter(&input, "INFO/AF < 1.0", false);
    let recs = parse_vcf_records(&out);
    // In missing_fields.vcf, only one record has an AF value (AF=.) and it's
    // missing — so the filter should select zero records.
    assert_eq!(recs.len(), 0, "expected no records to pass when INFO field is missing");
    for r in &recs {
        assert!(
            r.info.contains_key("AF") && !r.info["AF"].is_empty() && r.info["AF"] != ".",
            "unexpected row without AF: {:?}",
            r
        );
    }
}

/// 9. `QUAL > 30` numeric comparison.
#[test]
fn qual_numeric_comparison() {
    let (out, stats) = run_filter(TYPED_VCF, "QUAL > 30", false);
    let recs = parse_vcf_records(&out);
    // QUALs: 55, 18, 72, 25, 90 → keep 55, 72, 90.
    assert_eq!(stats.output_records, 3);
    for r in &recs {
        assert!(r.qual.unwrap() > 30.0);
    }
}

/// 10. `FILTER != 'PASS'`.
#[test]
fn filter_not_equals_pass() {
    let (out, _stats) = run_filter(TYPED_VCF, "FILTER != 'PASS'", false);
    let recs = parse_vcf_records(&out);
    assert_eq!(recs.len(), 2);
    for r in &recs {
        assert_ne!(r.filter, vec!["PASS".to_string()]);
    }
}

/// 11. Parentheses grouping.
#[test]
fn parentheses_grouping() {
    // Without parens, `&&` would bind tighter than `||`:
    //   CHROM=='chr17' && INFO/AF<0.05  -> chr17 only (AF=0.003) = 1 record
    //   FILTER=='PASS' || above         -> all PASS (3) union chr17 hit (already in PASS)
    // With parens below: (FILTER=='PASS' || CHROM=='chr17') && INFO/AF<0.05
    //   keep records that are PASS or chr17 AND AF<0.05.
    //   pos 100 (PASS, AF=0.02) ✓, pos 17:150 (PASS+chr17, 0.003) ✓,
    //   pos 17:250 (LowQual+chr17, 0.45) ✗, pos 300 (PASS, 0.001) ✓ → 3.
    let (out, _stats) = run_filter(
        TYPED_VCF,
        "(FILTER == 'PASS' || CHROM == 'chr17') && INFO/AF < 0.05",
        false,
    );
    let recs = parse_vcf_records(&out);
    assert_eq!(recs.len(), 3);
}

/// 12. All records from basic.vcf pass through with a permissive expression.
#[test]
fn basic_permissive_expression_passes_all() {
    let input = read_corpus("basic.vcf");
    let (out, stats) = run_filter(&input, "POS > 0", false);
    let recs = parse_vcf_records(&out);
    assert_eq!(stats.input_records, 5);
    assert_eq!(stats.output_records, 5);
    assert_eq!(recs.len(), 5);
}

/// 13. Parse error returns a helpful message (not a panic).
#[test]
fn parse_error_returns_helpful_message() {
    let err = FilterExpression::parse("INFO/AF << 3").unwrap_err();
    let msg = err.to_string();
    assert!(
        msg.contains("invalid filter expression"),
        "error should include 'invalid filter expression': {msg}"
    );
    assert!(msg.contains('^'), "error should include a caret: {msg}");
    assert!(
        msg.contains("INFO/AF << 3"),
        "error should echo the expression: {msg}"
    );
}

/// 14. FORMAT/GT field access on the first sample.
#[test]
fn format_gt_field_access() {
    // TYPED_VCF has GTs: 0/1, 0/1, 1/1, 0/1, 1/1 → `FORMAT/GT == '1/1'` = 2 records.
    let (out, stats) = run_filter(TYPED_VCF, "FORMAT/GT == '1/1'", false);
    let recs = parse_vcf_records(&out);
    assert_eq!(stats.output_records, 2);
    for r in &recs {
        assert!(
            r.samples.first().map(|s| s.starts_with("1/1")).unwrap_or(false),
            "expected 1/1 GT, got {:?}",
            r.samples
        );
    }
}

/// 15. multiline_info.vcf corpus: filter `INFO/AF < 0.05 && INFO/CSQ ~ 'missense'`.
#[test]
fn multiline_info_compound_filter() {
    let input = read_corpus("multiline_info.vcf");
    let (out, stats) = run_filter(
        &input,
        "INFO/AF < 0.05 && INFO/CSQ ~ 'missense'",
        false,
    );
    let recs = parse_vcf_records(&out);
    assert!(stats.input_records > 0);
    assert_eq!(stats.output_records, recs.len());
    // Every kept record must satisfy both predicates.
    for r in &recs {
        let af: f64 = r.info["AF"].parse().expect("AF is a float");
        assert!(af < 0.05, "AF={af} should be < 0.05");
        assert!(
            r.info["CSQ"].contains("missense"),
            "CSQ {:?} should contain 'missense'",
            r.info["CSQ"]
        );
    }
    // From the fixture: missense records with AF<0.05 are:
    //   chr1:925952 AF=0.03 missense ✓
    //   chr1:2000000 AF=0.02 missense ✓
    //   chr17:41244936 AF=0.01 missense ✓
    //   chr17:41276077 AF=0.04 splice_region ✗
    //   chr17:7674220 AF=0.03 missense ✓
    //   chr17:41194312 AF=0.06 missense ✗ (AF too high)
    // → 4 records expected.
    assert_eq!(recs.len(), 4, "expected 4 AF<0.05 AND missense records");
}

/// 16. `~` on a multi-valued INFO array: any element containing the substring
/// should match.
///
/// The inline VCF has `CSQ=intron_variant,missense_variant` for one record and
/// `CSQ=synonymous_variant` for another. Filtering with `INFO/CSQ ~ 'missense'`
/// should return only the first record.
#[test]
fn csq_contains_matches_any_array_element() {
    // Build a small inline VCF with a comma-separated multi-value CSQ field.
    let vcf: &[u8] = b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=248956422>\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tT\t.\t.\tCSQ=intron_variant,missense_variant\n\
chr1\t200\t.\tC\tG\t.\t.\tCSQ=synonymous_variant\n";

    let (out, stats) = run_filter(vcf, "INFO/CSQ ~ 'missense'", false);
    let recs = parse_vcf_records(&out);
    assert_eq!(stats.input_records, 2);
    assert_eq!(
        recs.len(),
        1,
        "expected exactly 1 record: the one with 'missense_variant' in CSQ"
    );
    assert!(
        recs[0].info["CSQ"].contains("missense"),
        "CSQ field did not contain 'missense': {}",
        recs[0].info["CSQ"]
    );
}

/// 17. Bare `!` negation without parentheses.
///
/// `! FILTER == 'PASS'` is equivalent to `!(FILTER == 'PASS')` because `!`
/// applies to the next atom. With 3 PASS and 2 LowQual records in `TYPED_VCF`,
/// we expect the 2 non-PASS records.
#[test]
fn bare_not_negation_without_parens() {
    // `!` binds to the next unary/atom, so `! FILTER == 'PASS'` should parse
    // as `!(FILTER == 'PASS')` and keep the two LowQual records.
    let (out, stats) = run_filter(TYPED_VCF, "! FILTER == 'PASS'", false);
    let recs = parse_vcf_records(&out);
    assert_eq!(
        recs.len(),
        2,
        "expected 2 non-PASS records, got {} (stats: {stats:?})",
        recs.len()
    );
    for r in &recs {
        assert_ne!(
            r.filter,
            vec!["PASS".to_string()],
            "record should not have PASS filter"
        );
    }
}

// ── CLI-level smoke test ─────────────────────────────────────────────────────

#[test]
fn cli_reads_stdin_and_writes_stdout() {
    let input = TYPED_VCF;
    let mut child = Command::new(vcfkit_bin())
        .args(["filter", "-e", "FILTER == 'PASS'"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn vcfkit");
    child
        .stdin
        .as_mut()
        .expect("stdin")
        .write_all(input)
        .expect("write stdin");
    let out = child.wait_with_output().expect("wait vcfkit");
    assert!(
        out.status.success(),
        "vcfkit filter exited {:?}: {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    let recs = parse_vcf_records(&stdout);
    assert_eq!(recs.len(), 3);
}

#[test]
fn cli_invert_flag_passes_non_matching() {
    let input = TYPED_VCF;
    let mut child = Command::new(vcfkit_bin())
        .args(["filter", "-e", "FILTER == 'PASS'", "--invert"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn vcfkit");
    child
        .stdin
        .as_mut()
        .expect("stdin")
        .write_all(input)
        .expect("write stdin");
    let out = child.wait_with_output().expect("wait vcfkit");
    assert!(out.status.success());
    let stdout = String::from_utf8(out.stdout).expect("utf-8");
    let recs = parse_vcf_records(&stdout);
    assert_eq!(recs.len(), 2);
    for r in &recs {
        assert_ne!(r.filter, vec!["PASS".to_string()]);
    }
}

#[test]
fn cli_bad_expression_exits_nonzero() {
    let input = TYPED_VCF;
    let mut child = Command::new(vcfkit_bin())
        .args(["filter", "-e", "INFO/AF << 3"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn vcfkit");
    child
        .stdin
        .as_mut()
        .expect("stdin")
        .write_all(input)
        .expect("write stdin");
    let out = child.wait_with_output().expect("wait vcfkit");
    assert!(
        !out.status.success(),
        "vcfkit filter should fail on a bad expression"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("invalid filter expression"),
        "stderr should mention 'invalid filter expression': {stderr}"
    );
}
