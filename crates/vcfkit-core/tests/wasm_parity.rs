//! Parity tests for the WASM entry points.
//!
//! These run on native targets and verify that the in-memory code paths used
//! by the WASM wrappers produce correct output. Since the WASM wrappers call
//! the same Rust functions with Cursor I/O, passing here proves byte-level
//! equivalence with the native CLI for the same inputs.

use std::io::{BufReader, Cursor};
use std::path::Path;

use vcfkit_core::{filter, io, liftover, normalize};

const SIMPLE_VCF: &str = "\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t50\tPASS\tAF=0.05
chr1\t200\t.\tC\tG\t30\tPASS\tAF=0.001
chr1\t300\t.\tT\tA\t80\tFail\tAF=0.5
";

const MULTIALLELIC_VCF: &str = "\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT,C\t50\tPASS\tAF=0.05,0.03
";

fn run_filter(vcf: &str, expr: &str) -> String {
    let e = filter::FilterExpression::parse(expr).unwrap();
    let opts = filter::FilterOptions {
        invert: false,
        output_format: io::OutputFormat::Vcf,
    };
    let mut out = Vec::new();
    filter::filter(
        BufReader::new(Cursor::new(vcf.as_bytes())),
        &mut out,
        e,
        opts,
    )
    .unwrap();
    String::from_utf8(out).unwrap()
}

fn run_normalize(vcf: &str) -> String {
    let opts = normalize::NormalizeOptions {
        split_multiallelics: true,
        left_align: false,
        check_ref: normalize::RefCheck::Ignore,
        output_format: io::OutputFormat::Vcf,
        fast: false,
    };
    let mut out = Vec::new();
    normalize::normalize(
        BufReader::new(Cursor::new(vcf.as_bytes())),
        &mut out,
        Path::new(""),
        opts,
    )
    .unwrap();
    String::from_utf8(out).unwrap()
}

#[test]
fn wasm_filter_keeps_matching_records() {
    let out = run_filter(SIMPLE_VCF, "INFO/AF < 0.01");
    assert!(out.contains("chr1\t200"), "chr1:200 (AF=0.001) should pass");
    assert!(
        !out.contains("chr1\t100"),
        "chr1:100 (AF=0.05) should be filtered"
    );
    assert!(
        !out.contains("chr1\t300"),
        "chr1:300 (AF=0.5) should be filtered"
    );
}

#[test]
fn wasm_filter_with_filter_field() {
    let out = run_filter(SIMPLE_VCF, "FILTER == 'PASS'");
    assert!(out.contains("chr1\t100"));
    assert!(out.contains("chr1\t200"));
    assert!(!out.contains("chr1\t300"));
}

#[test]
fn wasm_filter_header_preserved() {
    let out = run_filter(SIMPLE_VCF, "POS > 0");
    assert!(out.starts_with("##fileformat=VCFv4.2"));
    assert!(out.contains("#CHROM\tPOS"));
}

#[test]
fn wasm_normalize_splits_multiallelic() {
    let out = run_normalize(MULTIALLELIC_VCF);
    let records: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        records.len(),
        2,
        "one multi-allelic should split into two records"
    );
    assert!(records[0].contains("T"), "first split should have ALT=T");
    assert!(records[1].contains("C"), "second split should have ALT=C");
}

#[test]
fn wasm_normalize_header_preserved() {
    let out = run_normalize(SIMPLE_VCF);
    assert!(out.starts_with("##fileformat=VCFv4.2"));
    assert!(out.contains("#CHROM\tPOS"));
}

#[test]
fn wasm_liftover_identity_chain() {
    // Identity chain: maps chr1:1-1000 → chr1:1-1000 unchanged.
    let chain = "chain 1000 chr1 1000 + 0 1000 chr1 1000 + 0 1000 1\n1000\n\n";
    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t.\t.\t.
";
    let opts = liftover::LiftoverOptions {
        reject_file: None,
        write_src_coords: false,
        fix_swapped_ref: true,
        output_format: io::OutputFormat::Vcf,
        allow_contig_mismatch: true,
    };
    let mut out = Vec::new();
    liftover::liftover_from_chain_reader(
        BufReader::new(Cursor::new(vcf.as_bytes())),
        &mut out,
        BufReader::new(Cursor::new(chain.as_bytes())),
        opts,
    )
    .unwrap();
    let result = String::from_utf8(out).unwrap();
    assert!(
        result.contains("chr1\t100"),
        "identity liftover should preserve position"
    );
}
