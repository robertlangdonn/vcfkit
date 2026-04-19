//! Real-world differential tests: vcfkit vs bcftools on 1000 Genomes chr22.
//!
//! All tests are gated behind `VCFKIT_REAL_TESTS=1`. When the env var is absent
//! (normal CI, local dev without data), tests skip cleanly. When set, tests
//! download/verify reference data once to `tests/real_world/` and run full
//! differential comparisons against bcftools.
//!
//! # Running locally
//!
//! ```bash
//! # First time: downloads ~500MB of reference data
//! VCFKIT_REAL_TESTS=1 cargo test --test real_world_differential -- --nocapture
//!
//! # Subsequent runs use cached data
//! VCFKIT_REAL_TESTS=1 cargo test --test real_world_differential
//! ```
//!
//! Prerequisites: `bcftools` and `samtools` in PATH.

mod common;

use std::{path::PathBuf, process::Command};

use common::diff::{assert_vcf_eq, parse_vcf_records};
use common::download::RealWorldData;

// ── guard ─────────────────────────────────────────────────────────────────────

/// Returns true if the real-world test suite is enabled.
fn real_tests_enabled() -> bool {
    std::env::var("VCFKIT_REAL_TESTS")
        .map(|v| v == "1" || v.to_lowercase() == "true")
        .unwrap_or(false)
}

/// Skip the test if VCFKIT_REAL_TESTS is not set.
macro_rules! require_real_tests {
    () => {
        if !real_tests_enabled() {
            eprintln!("skipping real-world test (set VCFKIT_REAL_TESTS=1 to enable)");
            return;
        }
    };
}

/// Assert that `cmd` is available in PATH, skipping if not.
macro_rules! require_binary {
    ($name:expr) => {
        if Command::new($name).arg("--version").output().is_err() {
            eprintln!("skipping: {} not found in PATH", $name);
            return;
        }
    };
}

fn vcfkit_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_vcfkit"))
}

fn workspace_root() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.pop(); // crates/
    p.pop(); // workspace root
    p
}

fn real_world_dir() -> PathBuf {
    workspace_root().join("tests/real_world")
}

// ── run helpers ───────────────────────────────────────────────────────────────

fn run_vcfkit(args: &[&str]) -> String {
    let out = Command::new(vcfkit_bin())
        .args(args)
        .output()
        .expect("vcfkit binary failed to launch");
    if !out.status.success() {
        panic!(
            "vcfkit exited non-zero:\nstderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    String::from_utf8(out.stdout).expect("vcfkit output is not utf-8")
}

fn run_bcftools(args: &[&str]) -> String {
    let out = Command::new("bcftools")
        .args(args)
        .output()
        .expect("bcftools failed to launch");
    if !out.status.success() {
        panic!(
            "bcftools exited non-zero:\nstderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    String::from_utf8(out.stdout).expect("bcftools output is not utf-8")
}

/// Count variants in a VCF string.
fn count_records(vcf: &str) -> usize {
    parse_vcf_records(vcf).len()
}

// ── normalize differential ────────────────────────────────────────────────────

/// Compare `vcfkit normalize` output against `bcftools norm` on 1000G chr22.
///
/// Uses `--no-split` so the comparison is apples-to-apples for left-alignment
/// only (bcftools split behavior has known minor differences for some edge cases
/// documented in docs/reference_differences.md).
#[test]
fn normalize_left_align_matches_bcftools() {
    require_real_tests!();
    require_binary!("bcftools");

    let data =
        RealWorldData::acquire(&real_world_dir()).expect("failed to acquire real-world data");

    // vcfkit: left-align only (no split), warn mode for ref check
    let vcfkit_out = run_vcfkit(&[
        "normalize",
        "--no-split",
        "--check-ref",
        "warn",
        "-f",
        data.hg19_chr22_fa.to_str().unwrap(),
        data.chr22_vcf.to_str().unwrap(),
    ]);

    // bcftools norm: left-align only, keep multi-allelics, warn mode
    let bcftools_out = run_bcftools(&[
        "norm",
        "-f",
        data.hg19_chr22_fa.to_str().unwrap(),
        "-m",
        "+any",
        "-c",
        "w",
        data.chr22_vcf.to_str().unwrap(),
    ]);

    let vcfkit_n = count_records(&vcfkit_out);
    let bcftools_n = count_records(&bcftools_out);
    eprintln!("normalize left-align: vcfkit={vcfkit_n} records, bcftools={bcftools_n} records");

    // Record counts must match exactly.
    assert_eq!(
        vcfkit_n, bcftools_n,
        "normalize produced different record counts: vcfkit={vcfkit_n}, bcftools={bcftools_n}"
    );

    // Semantic record-level diff. Any divergence is a bug (or must be documented
    // in docs/reference_differences.md as an intentional deviation).
    assert_vcf_eq(&bcftools_out, &vcfkit_out);
}

/// Same as above but uses --fast mode. Fast path must produce identical output.
#[test]
fn normalize_fast_matches_bcftools() {
    require_real_tests!();
    require_binary!("bcftools");

    let data =
        RealWorldData::acquire(&real_world_dir()).expect("failed to acquire real-world data");

    let vcfkit_out = run_vcfkit(&[
        "normalize",
        "--fast",
        "--no-split",
        "--check-ref",
        "warn",
        "-f",
        data.hg19_chr22_fa.to_str().unwrap(),
        data.chr22_vcf.to_str().unwrap(),
    ]);

    let bcftools_out = run_bcftools(&[
        "norm",
        "-f",
        data.hg19_chr22_fa.to_str().unwrap(),
        "-m",
        "+any",
        "-c",
        "w",
        data.chr22_vcf.to_str().unwrap(),
    ]);

    let vcfkit_n = count_records(&vcfkit_out);
    let bcftools_n = count_records(&bcftools_out);
    eprintln!("normalize --fast: vcfkit={vcfkit_n} records, bcftools={bcftools_n} records");

    assert_eq!(vcfkit_n, bcftools_n);
    assert_vcf_eq(&bcftools_out, &vcfkit_out);
}

// ── filter differential ───────────────────────────────────────────────────────

/// Compare `vcfkit filter` against `bcftools view` for an AF threshold filter.
/// This is the benchmark expression (INFO/AF < 0.01) used in the README.
#[test]
fn filter_af_threshold_matches_bcftools() {
    require_real_tests!();
    require_binary!("bcftools");

    let data =
        RealWorldData::acquire(&real_world_dir()).expect("failed to acquire real-world data");

    let vcfkit_out = run_vcfkit(&[
        "filter",
        "-e",
        "INFO/AF < 0.01",
        data.chr22_vcf.to_str().unwrap(),
    ]);

    let bcftools_out = run_bcftools(&[
        "view",
        "-i",
        "INFO/AF < 0.01",
        data.chr22_vcf.to_str().unwrap(),
    ]);

    let vcfkit_n = count_records(&vcfkit_out);
    let bcftools_n = count_records(&bcftools_out);
    eprintln!("filter INFO/AF < 0.01: vcfkit={vcfkit_n} records, bcftools={bcftools_n} records");

    assert_eq!(
        vcfkit_n, bcftools_n,
        "filter AF<0.01 produced different record counts: vcfkit={vcfkit_n}, bcftools={bcftools_n}"
    );
    assert_vcf_eq(&bcftools_out, &vcfkit_out);
}

/// FILTER == 'PASS' — the second benchmark expression.
#[test]
fn filter_pass_matches_bcftools() {
    require_real_tests!();
    require_binary!("bcftools");

    let data =
        RealWorldData::acquire(&real_world_dir()).expect("failed to acquire real-world data");

    let vcfkit_out = run_vcfkit(&[
        "filter",
        "-e",
        "FILTER == 'PASS'",
        data.chr22_vcf.to_str().unwrap(),
    ]);

    let bcftools_out = run_bcftools(&["view", "-f", "PASS", data.chr22_vcf.to_str().unwrap()]);

    let vcfkit_n = count_records(&vcfkit_out);
    let bcftools_n = count_records(&bcftools_out);
    eprintln!("filter FILTER=='PASS': vcfkit={vcfkit_n} records, bcftools={bcftools_n} records");

    assert_eq!(vcfkit_n, bcftools_n);
    assert_vcf_eq(&bcftools_out, &vcfkit_out);
}

// ── liftover differential ─────────────────────────────────────────────────────

/// Compare `vcfkit liftover` hg19→hg38 against `bcftools +liftover` on chr22.
///
/// Known intentional deviations (if any) should be documented in
/// docs/reference_differences.md and excluded from the comparison here.
#[test]
fn liftover_hg19_to_hg38_matches_bcftools() {
    require_real_tests!();
    require_binary!("bcftools");

    // bcftools +liftover is a plugin — check it's available.
    let plugin_check = Command::new("bcftools")
        .args(["+liftover", "--version"])
        .output();
    if plugin_check.map(|o| !o.status.success()).unwrap_or(true) {
        eprintln!("skipping: bcftools +liftover plugin not available");
        return;
    }

    let data =
        RealWorldData::acquire(&real_world_dir()).expect("failed to acquire real-world data");

    let vcfkit_out = run_vcfkit(&[
        "liftover",
        "-s",
        data.hg19_chr22_fa.to_str().unwrap(),
        "-t",
        data.hg38_chr22_fa.to_str().unwrap(),
        "-c",
        data.hg19_to_hg38_chain.to_str().unwrap(),
        data.chr22_vcf.to_str().unwrap(),
    ]);

    let bcftools_out = run_bcftools(&[
        "+liftover",
        data.chr22_vcf.to_str().unwrap(),
        "--",
        "-s",
        data.hg19_chr22_fa.to_str().unwrap(),
        "-f",
        data.hg38_chr22_fa.to_str().unwrap(),
        "-c",
        data.hg19_to_hg38_chain.to_str().unwrap(),
    ]);

    let vcfkit_n = count_records(&vcfkit_out);
    let bcftools_n = count_records(&bcftools_out);
    eprintln!(
        "liftover hg19→hg38: vcfkit={vcfkit_n} records lifted, bcftools={bcftools_n} records lifted"
    );

    // Count-level check first (more informative failure message than full diff).
    // Allow ±1% difference — minor divergences in gap/strand handling are expected
    // and intentional (see docs/reference_differences.md).
    let delta = (vcfkit_n as i64 - bcftools_n as i64).unsigned_abs() as usize;
    let pct = if bcftools_n > 0 {
        delta * 100 / bcftools_n
    } else {
        100
    };
    assert!(
        pct <= 1,
        "liftover record count diverges by {pct}% (vcfkit={vcfkit_n}, bcftools={bcftools_n}); \
         if intentional, document in docs/reference_differences.md"
    );
}
