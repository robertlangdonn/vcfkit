//! `vcfkit filter` — keep VCF records that match an expression.

use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Write},
};

use anyhow::{anyhow, Context};
use owo_colors::OwoColorize;

use vcfkit_core::{
    filter::{extract_header_schema, filter_with_progress, FilterExpression, FilterOptions},
    io::OutputFormat,
};

use crate::english::{self, HeaderSchema};
use crate::output::ProgressReporter;
use crate::FilterArgs;

/// Run the filter subcommand.
pub fn run(args: &FilterArgs, quiet: bool) -> anyhow::Result<()> {
    super::reject_bcf_output(args.output.as_deref())?;

    // Determine the expression — either direct (-e) or translated (--english).
    let expression_str: String = if let Some(query) = &args.english {
        resolve_english(query, args, quiet)?
    } else {
        args.expression
            .clone()
            .expect("clap group ensures one of -e/--english is set")
    };

    let expression = FilterExpression::parse(&expression_str)
        .map_err(|e| anyhow!("{e}"))
        .with_context(|| format!("failed to parse expression: {expression_str}"))?;

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

    let stats =
        match (args.input.as_deref(), args.output.as_deref()) {
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

/// Translate `query` to a filter expression via the Anthropic API, show it to
/// the user for confirmation (unless `--yes`), and return the confirmed expression.
fn resolve_english(query: &str, args: &FilterArgs, quiet: bool) -> anyhow::Result<String> {
    // --english requires a file path — stdin can't be read twice.
    let in_path = args.input.as_deref().ok_or_else(|| {
        anyhow!("--english requires an input file path (stdin is not supported with --english)")
    })?;

    // Read the VCF header to build the schema for the LLM prompt.
    if !quiet {
        eprint!("Reading VCF header...");
        let _ = io::stderr().flush();
    }
    let header_reader = BufReader::new(
        File::open(in_path).with_context(|| format!("failed to open '{}'", in_path.display()))?,
    );
    let vcf_schema =
        extract_header_schema(header_reader).with_context(|| "failed to read VCF header")?;

    if !quiet {
        eprintln!(
            " done ({} INFO, {} FORMAT fields)",
            vcf_schema.info_fields.len(),
            vcf_schema.format_fields.len()
        );
        eprint!("Translating query via Anthropic API...");
        let _ = io::stderr().flush();
    }

    // Build HeaderSchema from vcfkit-core's VcfHeaderSchema.
    let schema = HeaderSchema {
        info_fields: vcf_schema
            .info_fields
            .iter()
            .map(|f| english::FieldDef {
                id: f.id.clone(),
                number: f.number.clone(),
                ty: f.ty.clone(),
                description: f.description.clone(),
            })
            .collect(),
        format_fields: vcf_schema
            .format_fields
            .iter()
            .map(|f| english::FieldDef {
                id: f.id.clone(),
                number: f.number.clone(),
                ty: f.ty.clone(),
                description: f.description.clone(),
            })
            .collect(),
        contigs: vcf_schema.contigs.clone(),
    };

    // Run the async translation in a single-threaded tokio runtime.
    let rt = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()
        .map_err(|e| anyhow!("failed to create async runtime: {e}"))?;
    let translation = rt
        .block_on(english::translate(query, &schema))
        .map_err(|e| anyhow!("{e}"))?;

    if !quiet {
        eprintln!(
            " done (model: {}, confidence: {:.0}%)",
            translation.model,
            translation.confidence * 100.0
        );
    }

    // Display the translation and ask for confirmation (unless --yes).
    eprintln!();
    eprintln!("  Query:      {}", query.bold());
    eprintln!("  Expression: {}", translation.expression.green().bold());
    eprintln!("  Reasoning:  {}", translation.reasoning);
    if !translation.caveats.is_empty() {
        for caveat in &translation.caveats {
            eprintln!("  {}: {}", "Caveat".yellow(), caveat);
        }
    }
    eprintln!();

    if args.yes {
        return Ok(translation.expression);
    }

    // Interactive confirmation.
    eprint!("Run this filter? [Y/n/edit] ");
    let _ = io::stderr().flush();

    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .with_context(|| "failed to read confirmation")?;

    match input.trim().to_lowercase().as_str() {
        "" | "y" | "yes" => Ok(translation.expression),
        "n" | "no" => Err(anyhow!("filter cancelled")),
        "edit" | "e" => open_in_editor(translation.expression),
        other => Err(anyhow!(
            "unrecognised response '{other}'; expected Y, n, or edit"
        )),
    }
}

/// Open `expression` in `$EDITOR`, let the user modify it, and return the
/// edited value. Falls back to `vi` when `$EDITOR` is not set.
fn open_in_editor(expression: String) -> anyhow::Result<String> {
    let tmp = std::env::temp_dir().join("vcfkit_english_edit.tmp");
    std::fs::write(&tmp, &expression)
        .with_context(|| format!("failed to write temp file '{}'", tmp.display()))?;

    let editor = std::env::var("EDITOR").unwrap_or_else(|_| "vi".to_string());
    let status = std::process::Command::new(&editor)
        .arg(&tmp)
        .status()
        .with_context(|| format!("failed to launch editor '{editor}'"))?;
    if !status.success() {
        return Err(anyhow!("editor exited with non-zero status"));
    }

    let edited =
        std::fs::read_to_string(&tmp).with_context(|| "failed to read edited expression")?;
    let edited = edited.trim().to_string();
    if edited.is_empty() {
        return Err(anyhow!("edited expression is empty"));
    }

    // Validate the edited expression before running.
    FilterExpression::parse(&edited)
        .map_err(|e| anyhow!("edited expression failed to parse: {e}"))?;

    Ok(edited)
}
