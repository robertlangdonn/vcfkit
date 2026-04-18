//! `vcfkit filter` — keep VCF records that match an expression.

use std::{
    fs::File,
    io::{self, BufReader, BufWriter},
};

use anyhow::{Context, anyhow};

use vcfkit_core::{
    filter::{FilterExpression, FilterOptions, filter_with_progress},
    io::OutputFormat,
};

use crate::FilterArgs;
use crate::output::ProgressReporter;

/// Run the filter subcommand.
pub fn run(args: &FilterArgs, quiet: bool) -> anyhow::Result<()> {
    // Parse the filter expression up front so CLI errors surface before we
    // open input/output files.
    let expression = FilterExpression::parse(&args.expression)
        .map_err(|e| anyhow!("{e}"))
        .with_context(|| format!("failed to parse expression: {}", args.expression))?;

    let options = FilterOptions {
        invert: args.invert,
        output_format: args
            .output
            .as_deref()
            .map(vcfkit_core::io::format_from_path)
            .unwrap_or(OutputFormat::Vcf),
    };

    let reporter = ProgressReporter::new_with_flags(None, quiet, args.no_progress);
    let on_record = |_n: u64| reporter.inc();

    let stats = match (args.input.as_deref(), args.output.as_deref()) {
        (Some(in_path), Some(out_path)) => {
            let reader = BufReader::new(File::open(in_path).with_context(|| {
                format!("failed to open input file '{}'", in_path.display())
            })?);
            let writer = BufWriter::new(File::create(out_path).with_context(|| {
                format!("failed to create output file '{}'", out_path.display())
            })?);
            filter_with_progress(reader, writer, expression, options, on_record)
                .with_context(|| "filter failed")?
        }
        (Some(in_path), None) => {
            let reader = BufReader::new(File::open(in_path).with_context(|| {
                format!("failed to open input file '{}'", in_path.display())
            })?);
            let stdout = io::stdout();
            let writer = BufWriter::new(stdout.lock());
            filter_with_progress(reader, writer, expression, options, on_record)
                .with_context(|| "filter failed")?
        }
        (None, Some(out_path)) => {
            let stdin = io::stdin();
            let reader = BufReader::new(stdin.lock());
            let writer = BufWriter::new(File::create(out_path).with_context(|| {
                format!("failed to create output file '{}'", out_path.display())
            })?);
            filter_with_progress(reader, writer, expression, options, on_record)
                .with_context(|| "filter failed")?
        }
        (None, None) => {
            let stdin = io::stdin();
            let stdout = io::stdout();
            let reader = BufReader::new(stdin.lock());
            let writer = BufWriter::new(stdout.lock());
            filter_with_progress(reader, writer, expression, options, on_record)
                .with_context(|| "filter failed")?
        }
    };

    reporter.finish("filter complete");

    if !quiet {
        eprintln!(
            "filter: {} in, {} out ({} filtered out)",
            stats.input_records, stats.output_records, stats.filtered_out,
        );
    }

    Ok(())
}
