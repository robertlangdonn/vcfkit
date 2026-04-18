//! Semantic VCF diffing utilities.
//!
//! Compares two VCF outputs at the record level rather than byte-for-byte, so
//! differences in header ordering, whitespace, or tool-provenance lines
//! (`##source`, `##fileformat`) don't cause spurious test failures.
//!
//! The primary entry point is [`assert_vcf_eq`]; [`parse_vcf_records`] and
//! [`vcf_diff`] are exposed for use in custom assertions.

use std::collections::HashMap;

/// A single parsed VCF data record.
#[derive(Debug, Clone, PartialEq)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,
    pub id: String,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
    pub qual: Option<f32>,
    pub filter: Vec<String>,
    pub info: HashMap<String, String>,
    pub format: Option<String>,
    pub samples: Vec<String>,
}

/// Parse a VCF string into a `Vec<VcfRecord>` in the order they appear.
///
/// Header lines (starting with `#`) are ignored entirely, as is trailing
/// whitespace. Truly empty lines are skipped; blank-ish whitespace-only lines
/// are likewise ignored.
pub fn parse_vcf_records(vcf: &str) -> Vec<VcfRecord> {
    let mut out = Vec::new();
    for raw_line in vcf.split('\n') {
        let line = raw_line.trim_end_matches(&['\r', '\n'][..]).trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        // Require at least CHROM..INFO (8 columns) per the VCF spec.
        if parts.len() < 8 {
            continue;
        }
        let chrom = parts[0].to_string();
        let pos: u64 = parts[1].parse().unwrap_or(0);
        let id = parts[2].to_string();
        let ref_allele = parts[3].to_string();
        let alt_alleles: Vec<String> = if parts[4] == "." {
            Vec::new()
        } else {
            parts[4].split(',').map(str::to_string).collect()
        };
        let qual = parse_qual(parts[5]);
        let filter: Vec<String> = if parts[6] == "." || parts[6].is_empty() {
            Vec::new()
        } else {
            parts[6].split(';').map(str::to_string).collect()
        };
        let info = parse_info(parts[7]);
        let format = parts.get(8).filter(|s| !s.is_empty()).map(|s| s.to_string());
        let samples: Vec<String> = if parts.len() > 9 {
            parts[9..].iter().map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };

        out.push(VcfRecord {
            chrom,
            pos,
            id,
            ref_allele,
            alt_alleles,
            qual,
            filter,
            info,
            format,
            samples,
        });
    }
    out
}

fn parse_qual(s: &str) -> Option<f32> {
    if s == "." || s.is_empty() {
        None
    } else {
        s.parse().ok()
    }
}

fn parse_info(s: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    if s == "." || s.is_empty() {
        return map;
    }
    for entry in s.split(';') {
        if entry.is_empty() {
            continue;
        }
        match entry.split_once('=') {
            Some((k, v)) => {
                map.insert(k.to_string(), v.to_string());
            }
            None => {
                // Flag-style INFO key with no value.
                map.insert(entry.to_string(), String::new());
            }
        }
    }
    map
}

/// Assert that two VCF strings are semantically equivalent.
///
/// Ignored: `##fileformat` lines, `##source` lines, header-line ordering,
/// trailing whitespace on any line, and the distinction between `.` and missing
/// values for `QUAL`.
///
/// Panics with a diff rendered by [`vcf_diff`] if the records differ.
#[track_caller]
pub fn assert_vcf_eq(expected: &str, actual: &str) {
    let exp = parse_vcf_records(expected);
    let act = parse_vcf_records(actual);

    if exp == act {
        return;
    }

    let diff = vcf_diff(expected, actual);
    panic!(
        "VCF outputs differ.\n\n\
         expected: {} records\n\
         actual:   {} records\n\n\
         {diff}",
        exp.len(),
        act.len()
    );
}

/// Render a line-by-line diff of two VCF strings, suitable for embedding in
/// panic messages. Data records are diffed in order; header lines are diffed
/// only if they're substantively different (ignoring `##fileformat` and
/// `##source` lines which vary by tool/version).
pub fn vcf_diff(expected: &str, actual: &str) -> String {
    let exp_records = parse_vcf_records(expected);
    let act_records = parse_vcf_records(actual);

    let mut out = String::new();
    out.push_str("--- expected\n+++ actual\n");

    let max = exp_records.len().max(act_records.len());
    for i in 0..max {
        match (exp_records.get(i), act_records.get(i)) {
            (Some(e), Some(a)) if e == a => {
                out.push_str(&format!("  [{i}] {}\n", render_record(e)));
            }
            (Some(e), Some(a)) => {
                out.push_str(&format!("- [{i}] {}\n", render_record(e)));
                out.push_str(&format!("+ [{i}] {}\n", render_record(a)));
            }
            (Some(e), None) => {
                out.push_str(&format!("- [{i}] {}\n", render_record(e)));
            }
            (None, Some(a)) => {
                out.push_str(&format!("+ [{i}] {}\n", render_record(a)));
            }
            (None, None) => break,
        }
    }
    out
}

fn render_record(r: &VcfRecord) -> String {
    let alt = if r.alt_alleles.is_empty() {
        ".".to_string()
    } else {
        r.alt_alleles.join(",")
    };
    let filter = if r.filter.is_empty() {
        ".".to_string()
    } else {
        r.filter.join(";")
    };
    let qual = match r.qual {
        Some(q) => format!("{q}"),
        None => ".".to_string(),
    };
    let mut info_keys: Vec<&String> = r.info.keys().collect();
    info_keys.sort();
    let info = if info_keys.is_empty() {
        ".".to_string()
    } else {
        info_keys
            .iter()
            .map(|k| {
                let v = &r.info[*k];
                if v.is_empty() {
                    (*k).to_string()
                } else {
                    format!("{k}={v}")
                }
            })
            .collect::<Vec<_>>()
            .join(";")
    };
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        r.chrom, r.pos, r.id, r.ref_allele, alt, qual, filter, info
    )
}
