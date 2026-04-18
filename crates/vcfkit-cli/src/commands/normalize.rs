//! `vcfkit normalize` — left-align indels, split multi-allelic records, and
//! optionally check REF against the reference FASTA.

use std::{
    fs::File,
    io::{self, BufReader, BufWriter},
};

use anyhow::Context;

use vcfkit_core::{
    io::OutputFormat,
    normalize::{NormalizeOptions, RefCheck, normalize_with_progress},
};

use crate::output::ProgressReporter;
use crate::{CheckRefMode, NormalizeArgs};

/// Run the normalize subcommand.
pub fn run(args: &NormalizeArgs, quiet: bool) -> anyhow::Result<()> {
    // Build options from CLI args.
    let options = NormalizeOptions {
        split_multiallelics: !args.no_split,
        left_align: !args.no_left_align,
        check_ref: match args.check_ref {
            CheckRefMode::Ignore => RefCheck::Ignore,
            CheckRefMode::Warn => RefCheck::Warn,
            CheckRefMode::Error => RefCheck::Error,
        },
        output_format: args
            .output
            .as_deref()
            .map(vcfkit_core::io::format_from_path)
            .unwrap_or(OutputFormat::Vcf),
    };

    // Validate reference up front so the user sees the clearest possible
    // error before we touch the input stream.
    if !args.reference.exists() {
        return Err(anyhow::anyhow!(
            "failed to load reference FASTA '{}': file not found",
            args.reference.display()
        ));
    }

    let reporter = ProgressReporter::new(None, quiet);
    let on_record = |_n: u64| reporter.inc();

    // Open input (path or stdin).
    let stats = match (args.input.as_deref(), args.output.as_deref()) {
        (Some(in_path), Some(out_path)) => {
            let reader = BufReader::new(File::open(in_path).with_context(|| {
                format!("failed to open input file '{}'", in_path.display())
            })?);
            let writer = BufWriter::new(File::create(out_path).with_context(|| {
                format!("failed to create output file '{}'", out_path.display())
            })?);
            normalize_with_progress(reader, writer, &args.reference, options, on_record)
                .with_context(|| "normalize failed")?
        }
        (Some(in_path), None) => {
            let reader = BufReader::new(File::open(in_path).with_context(|| {
                format!("failed to open input file '{}'", in_path.display())
            })?);
            let stdout = io::stdout();
            let writer = BufWriter::new(stdout.lock());
            normalize_with_progress(reader, writer, &args.reference, options, on_record)
                .with_context(|| "normalize failed")?
        }
        (None, Some(out_path)) => {
            let stdin = io::stdin();
            let reader = BufReader::new(stdin.lock());
            let writer = BufWriter::new(File::create(out_path).with_context(|| {
                format!("failed to create output file '{}'", out_path.display())
            })?);
            normalize_with_progress(reader, writer, &args.reference, options, on_record)
                .with_context(|| "normalize failed")?
        }
        (None, None) => {
            let stdin = io::stdin();
            let stdout = io::stdout();
            let reader = BufReader::new(stdin.lock());
            let writer = BufWriter::new(stdout.lock());
            normalize_with_progress(reader, writer, &args.reference, options, on_record)
                .with_context(|| "normalize failed")?
        }
    };

    reporter.finish("normalize complete");

    if !quiet {
        eprintln!(
            "normalize: {} in, {} out ({} left-aligned, {} multi-allelic sites split, {} REF mismatches)",
            stats.input_records,
            stats.output_records,
            stats.left_aligned,
            stats.split_sites,
            stats.ref_mismatches,
        );
    }

    Ok(())
}
