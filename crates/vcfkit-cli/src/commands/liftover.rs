//! `vcfkit liftover` — convert variants between reference builds using a
//! UCSC chain file.

use std::{
    fs::File,
    io::{self, BufReader, BufWriter},
};

use anyhow::Context;
use tracing::warn;

use vcfkit_core::{
    io::OutputFormat,
    liftover::{liftover_with_progress, LiftoverOptions, KNOWN_CHAIN_URLS},
};

use crate::output::ProgressReporter;
use crate::LiftoverArgs;

/// Print the known chain-file URLs to stdout and return.
pub fn print_known_chains() {
    for (name, url) in KNOWN_CHAIN_URLS {
        println!("{name}\t{url}");
    }
}

/// Run the liftover subcommand.
pub fn run(args: &LiftoverArgs, quiet: bool) -> anyhow::Result<()> {
    if args.list_chains {
        print_known_chains();
        return Ok(());
    }

    // After the early-return above, source_ref and chain are guaranteed
    // present by clap's required_unless_present attributes. target_ref is
    // fully optional.
    let source_ref = args
        .source_ref
        .as_deref()
        .context("internal: --source-ref missing after list_chains check")?;
    let target_ref = args.target_ref.as_deref();
    let chain_path = args
        .chain
        .as_deref()
        .context("internal: --chain missing after list_chains check")?;

    if !source_ref.exists() {
        return Err(anyhow::anyhow!(
            "failed to load reference FASTA '{}': file not found",
            source_ref.display()
        ));
    }
    if let Some(tref) = target_ref {
        if !tref.exists() {
            return Err(anyhow::anyhow!(
                "failed to load reference FASTA '{}': file not found",
                tref.display()
            ));
        }
    } else {
        warn!("no target reference provided — REF alleles will not be validated after liftover");
    }
    if !chain_path.exists() {
        return Err(anyhow::anyhow!(
            "failed to open chain file '{}': file not found",
            chain_path.display()
        ));
    }

    let options = LiftoverOptions {
        reject_file: args.reject.clone(),
        write_src_coords: args.write_src_coords,
        fix_swapped_ref: !args.no_fix_swapped_ref,
        output_format: args
            .output
            .as_deref()
            .map(vcfkit_core::io::format_from_path)
            .unwrap_or(OutputFormat::Vcf),
    };

    let reporter = ProgressReporter::new_with_flags(None, quiet, args.no_progress);
    let on_record = |_n: u64| reporter.inc();

    let stats =
        match (args.input.as_deref(), args.output.as_deref()) {
            (Some(in_path), Some(out_path)) => {
                let reader = BufReader::new(File::open(in_path).with_context(|| {
                    format!("failed to open input file '{}'", in_path.display())
                })?);
                let writer = BufWriter::new(File::create(out_path).with_context(|| {
                    format!("failed to create output file '{}'", out_path.display())
                })?);
                liftover_with_progress(
                    reader, writer, chain_path, source_ref, target_ref, options, on_record,
                )
                .with_context(|| "liftover failed")?
            }
            (Some(in_path), None) => {
                let reader = BufReader::new(File::open(in_path).with_context(|| {
                    format!("failed to open input file '{}'", in_path.display())
                })?);
                let stdout = io::stdout();
                let writer = BufWriter::new(stdout.lock());
                liftover_with_progress(
                    reader, writer, chain_path, source_ref, target_ref, options, on_record,
                )
                .with_context(|| "liftover failed")?
            }
            (None, Some(out_path)) => {
                let stdin = io::stdin();
                let reader = BufReader::new(stdin.lock());
                let writer = BufWriter::new(File::create(out_path).with_context(|| {
                    format!("failed to create output file '{}'", out_path.display())
                })?);
                liftover_with_progress(
                    reader, writer, chain_path, source_ref, target_ref, options, on_record,
                )
                .with_context(|| "liftover failed")?
            }
            (None, None) => {
                let stdin = io::stdin();
                let stdout = io::stdout();
                let reader = BufReader::new(stdin.lock());
                let writer = BufWriter::new(stdout.lock());
                liftover_with_progress(
                    reader, writer, chain_path, source_ref, target_ref, options, on_record,
                )
                .with_context(|| "liftover failed")?
            }
        };

    reporter.finish("liftover complete");

    if !quiet {
        eprintln!(
            "liftover: {} in, {} out ({} unmapped, {} REF mismatches, {} strand-flipped)",
            stats.input_records,
            stats.output_records,
            stats.rejected_unmapped,
            stats.rejected_ref_mismatch,
            stats.swapped_alleles,
        );
    }

    Ok(())
}
