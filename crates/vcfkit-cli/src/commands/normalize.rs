//! `vcfkit normalize` — left-align indels, split multi-allelic records, and
//! optionally check REF against the reference FASTA.

use std::{
    fs::File,
    io::{self, BufReader, BufWriter},
};

use anyhow::Context;

use vcfkit_core::{
    io::OutputFormat,
    normalize::{NormalizeOptions, RefCheck, normalize},
};

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

    // Open input (path or stdin).
    let stats = match (args.input.as_deref(), args.output.as_deref()) {
        (Some(in_path), Some(out_path)) => {
            let reader = BufReader::new(
                File::open(in_path)
                    .with_context(|| format!("opening input {}", in_path.display()))?,
            );
            let writer = BufWriter::new(
                File::create(out_path)
                    .with_context(|| format!("creating output {}", out_path.display()))?,
            );
            normalize(reader, writer, &args.reference, options)?
        }
        (Some(in_path), None) => {
            let reader = BufReader::new(
                File::open(in_path)
                    .with_context(|| format!("opening input {}", in_path.display()))?,
            );
            let stdout = io::stdout();
            let writer = BufWriter::new(stdout.lock());
            normalize(reader, writer, &args.reference, options)?
        }
        (None, Some(out_path)) => {
            let stdin = io::stdin();
            let reader = BufReader::new(stdin.lock());
            let writer = BufWriter::new(
                File::create(out_path)
                    .with_context(|| format!("creating output {}", out_path.display()))?,
            );
            normalize(reader, writer, &args.reference, options)?
        }
        (None, None) => {
            let stdin = io::stdin();
            let stdout = io::stdout();
            let reader = BufReader::new(stdin.lock());
            let writer = BufWriter::new(stdout.lock());
            normalize(reader, writer, &args.reference, options)?
        }
    };

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
