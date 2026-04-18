//! `vcfkit liftover` — convert variants between reference builds using a
//! UCSC chain file.

use std::{
    fs::File,
    io::{self, BufReader, BufWriter},
};

use anyhow::Context;

use vcfkit_core::{
    io::OutputFormat,
    liftover::{KNOWN_CHAIN_URLS, LiftoverOptions, liftover},
};

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

    // After the early-return above, these three are guaranteed present by
    // clap's required_unless_present attributes. Failing here would indicate
    // a CLI schema bug.
    let source_ref = args
        .source_ref
        .as_deref()
        .context("internal: --source-ref missing after list_chains check")?;
    let target_ref = args
        .target_ref
        .as_deref()
        .context("internal: --target-ref missing after list_chains check")?;
    let chain_path = args
        .chain
        .as_deref()
        .context("internal: --chain missing after list_chains check")?;

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
            liftover(reader, writer, chain_path, source_ref, target_ref, options)?
        }
        (Some(in_path), None) => {
            let reader = BufReader::new(
                File::open(in_path)
                    .with_context(|| format!("opening input {}", in_path.display()))?,
            );
            let stdout = io::stdout();
            let writer = BufWriter::new(stdout.lock());
            liftover(reader, writer, chain_path, source_ref, target_ref, options)?
        }
        (None, Some(out_path)) => {
            let stdin = io::stdin();
            let reader = BufReader::new(stdin.lock());
            let writer = BufWriter::new(
                File::create(out_path)
                    .with_context(|| format!("creating output {}", out_path.display()))?,
            );
            liftover(reader, writer, chain_path, source_ref, target_ref, options)?
        }
        (None, None) => {
            let stdin = io::stdin();
            let stdout = io::stdout();
            let reader = BufReader::new(stdin.lock());
            let writer = BufWriter::new(stdout.lock());
            liftover(reader, writer, chain_path, source_ref, target_ref, options)?
        }
    };

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
