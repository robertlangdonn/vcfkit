//! WASM bindings for vcfkit-core operations.
//!
//! Compiled only when targeting `wasm32`. Build with:
//!
//! ```text
//! wasm-pack build crates/vcfkit-core --target web \
//!   --out-dir ../../web/public/wasm --release
//! ```
//!
//! Limitations vs the native CLI:
//! - `normalize_vcf` splits multi-allelic records but does **not** left-align
//!   indels (no FASTA available in the browser).
//! - `liftover_vcf` skips REF validation against the target reference.

use std::io::{BufReader, Cursor};
use std::path::Path;

use wasm_bindgen::prelude::*;

use crate::{filter, io, liftover, normalize};

#[wasm_bindgen(start)]
pub fn init() {
    console_error_panic_hook::set_once();
}

/// Filter a VCF (as a string) by an expression. Returns filtered VCF including header.
///
/// # Example (JavaScript)
/// ```js
/// import init, { filter_vcf } from './vcfkit_core.js';
/// await init();
/// const out = filter_vcf(vcfText, "INFO/AF < 0.01 && FILTER == 'PASS'");
/// ```
#[wasm_bindgen]
pub fn filter_vcf(vcf_input: &str, expression: &str) -> Result<String, JsError> {
    let expr = filter::FilterExpression::parse(expression)
        .map_err(|e| JsError::new(&format!("invalid expression: {e}")))?;
    let options = filter::FilterOptions {
        invert: false,
        output_format: io::OutputFormat::Vcf,
    };
    let reader = BufReader::new(Cursor::new(vcf_input.as_bytes()));
    let mut output = Vec::new();
    filter::filter(reader, &mut output, expr, options).map_err(|e| JsError::new(&e.to_string()))?;
    String::from_utf8(output).map_err(|e| JsError::new(&e.to_string()))
}

/// Split multi-allelic VCF records. Left-alignment is skipped (requires a
/// reference FASTA, which is not available in WASM). Returns normalized VCF.
#[wasm_bindgen]
pub fn normalize_vcf(vcf_input: &str) -> Result<String, JsError> {
    let options = normalize::NormalizeOptions {
        split_multiallelics: true,
        left_align: false,
        check_ref: normalize::RefCheck::Ignore,
        output_format: io::OutputFormat::Vcf,
        fast: false,
    };
    let reader = BufReader::new(Cursor::new(vcf_input.as_bytes()));
    let mut output = Vec::new();
    // Path is never opened when left_align=false and check_ref=Ignore.
    normalize::normalize(reader, &mut output, Path::new(""), options)
        .map_err(|e| JsError::new(&e.to_string()))?;
    String::from_utf8(output).map_err(|e| JsError::new(&e.to_string()))
}

/// Liftover a VCF using a chain file provided as bytes. REF validation is
/// skipped (no FASTA). Unmapped records are silently dropped. Returns lifted VCF.
#[wasm_bindgen]
pub fn liftover_vcf(vcf_input: &str, chain_bytes: &[u8]) -> Result<String, JsError> {
    let options = liftover::LiftoverOptions {
        reject_file: None,
        write_src_coords: false,
        fix_swapped_ref: true,
        output_format: io::OutputFormat::Vcf,
        allow_contig_mismatch: true,
    };
    let reader = BufReader::new(Cursor::new(vcf_input.as_bytes()));
    let mut output = Vec::new();
    liftover::liftover_from_chain_reader(
        reader,
        &mut output,
        BufReader::new(Cursor::new(chain_bytes)),
        options,
    )
    .map_err(|e| JsError::new(&e.to_string()))?;
    String::from_utf8(output).map_err(|e| JsError::new(&e.to_string()))
}
