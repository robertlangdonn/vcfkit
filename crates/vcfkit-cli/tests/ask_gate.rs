//! Integration tests for --ask confidence gate and mock translation path.
//!
//! All tests use VCFKIT_MOCK_TRANSLATION to bypass the Anthropic API entirely.

use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_vcfkit");

// Minimal VCF with one variant.
const MINIMAL_VCF: &str = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t50\tPASS\tAF=0.005
";

fn write_temp_vcf(contents: &str) -> tempfile::NamedTempFile {
    use std::io::Write;
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(contents.as_bytes()).unwrap();
    f
}

fn high_conf_mock() -> String {
    r#"{"expression":"INFO/AF < 0.01","reasoning":"Rare variants.","confidence":0.9,"caveats":[]}"#
        .to_string()
}

fn low_conf_mock() -> String {
    r#"{"expression":"INFO/AF < 0.01","reasoning":"Best guess.","confidence":0.3,"caveats":["Field not found"]}"#
        .to_string()
}

/// --ask --yes with high confidence should run the filter and exit 0.
#[test]
fn ask_yes_high_confidence_runs_filter() {
    let vcf = write_temp_vcf(MINIMAL_VCF);
    let out = Command::new(BIN)
        .env("VCFKIT_MOCK_TRANSLATION", high_conf_mock())
        .args([
            "filter",
            "--ask",
            "rare variants",
            "--yes",
            vcf.path().to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("#CHROM"), "expected VCF header in output");
}

/// --ask --yes with low confidence (< 50%) should exit non-zero and print an
/// error telling the user to add --accept-low-confidence.
#[test]
fn ask_yes_low_confidence_blocked() {
    let vcf = write_temp_vcf(MINIMAL_VCF);
    let out = Command::new(BIN)
        .env("VCFKIT_MOCK_TRANSLATION", low_conf_mock())
        .args([
            "filter",
            "--ask",
            "some ambiguous query",
            "--yes",
            vcf.path().to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "expected non-zero exit for low-confidence + --yes"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("accept-low-confidence") || stderr.contains("confidence"),
        "expected confidence error in stderr, got: {stderr}"
    );
}

/// --ask --yes --accept-low-confidence with low confidence should run.
#[test]
fn ask_yes_accept_low_confidence_runs() {
    let vcf = write_temp_vcf(MINIMAL_VCF);
    let out = Command::new(BIN)
        .env("VCFKIT_MOCK_TRANSLATION", low_conf_mock())
        .args([
            "filter",
            "--ask",
            "some ambiguous query",
            "--yes",
            "--accept-low-confidence",
            vcf.path().to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("#CHROM"), "expected VCF header in output");
}

/// --ask --yes at exactly 50% confidence should be blocked (boundary case).
/// Previously slipped through because the gate used strict < 0.5.
#[test]
fn ask_yes_exactly_50_percent_is_blocked() {
    let vcf = write_temp_vcf(MINIMAL_VCF);
    let boundary_mock = r#"{"expression":"FORMAT/AD > 19","reasoning":"Workaround.","confidence":0.5,"caveats":["Matches wrong records"]}"#;
    let out = Command::new(BIN)
        .env("VCFKIT_MOCK_TRANSLATION", boundary_mock)
        .args([
            "filter",
            "--ask",
            "variants with exactly 20 alt reads",
            "--yes",
            vcf.path().to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "exactly 50% confidence should be blocked, not allowed through"
    );
}

/// -e and --ask are mutually exclusive.
#[test]
fn expression_and_ask_are_mutually_exclusive() {
    let vcf = write_temp_vcf(MINIMAL_VCF);
    let out = Command::new(BIN)
        .env("VCFKIT_MOCK_TRANSLATION", high_conf_mock())
        .args([
            "filter",
            "-e",
            "QUAL > 30",
            "--ask",
            "rare variants",
            vcf.path().to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "expected non-zero exit when both -e and --ask are given"
    );
}
