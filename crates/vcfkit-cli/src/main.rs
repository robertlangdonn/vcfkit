use std::io;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use anyhow::Result;
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use clap_complete::{generate, Shell};
use tracing_subscriber::{fmt, EnvFilter};

mod commands;
mod output;
mod telemetry;

// ---------------------------------------------------------------------------
// Top-level CLI
// ---------------------------------------------------------------------------

/// vcfkit — bioinformatics VCF toolkit
#[derive(Debug, Parser)]
#[command(name = "vcfkit", version, about, long_about = None)]
pub struct Cli {
    /// Increase verbosity (-v = info, -vv = debug, -vvv = trace)
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,

    /// Suppress non-error output
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    pub quiet: bool,

    /// Disable telemetry for this run
    #[arg(long, global = true)]
    pub no_telemetry: bool,

    /// Control colour output
    #[arg(long, value_name = "WHEN", default_value = "auto", global = true)]
    pub color: ColorWhen,

    #[command(subcommand)]
    pub command: Commands,
}

// ---------------------------------------------------------------------------
// Colour flag
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, ValueEnum)]
pub enum ColorWhen {
    Auto,
    Always,
    Never,
}

// ---------------------------------------------------------------------------
// Subcommands enum
// ---------------------------------------------------------------------------

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Normalize a VCF (left-align, split multi-allelics)
    Normalize(NormalizeArgs),

    /// Convert variants between reference builds
    Liftover(LiftoverArgs),

    /// Filter variants by expression
    Filter(FilterArgs),

    /// Generate shell completion scripts
    Completions(CompletionsArgs),
}

// ---------------------------------------------------------------------------
// normalize
// ---------------------------------------------------------------------------

#[derive(Debug, Parser)]
pub struct NormalizeArgs {
    /// Reference genome FASTA (required)
    #[arg(short = 'f', long, value_name = "FASTA", required = true)]
    pub reference: PathBuf,

    /// Output file (default: stdout)
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Do not split multi-allelic sites
    #[arg(long)]
    pub no_split: bool,

    /// Do not left-align indels
    #[arg(long)]
    pub no_left_align: bool,

    /// How to handle reference mismatches
    #[arg(long, value_name = "MODE", default_value = "warn")]
    pub check_ref: CheckRefMode,

    /// Suppress the progress bar even when stderr is a TTY
    #[arg(long)]
    pub no_progress: bool,

    /// Input VCF/BCF (default: stdin)
    #[arg(value_name = "INPUT")]
    pub input: Option<PathBuf>,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum CheckRefMode {
    Ignore,
    Warn,
    Error,
}

// ---------------------------------------------------------------------------
// liftover
// ---------------------------------------------------------------------------

#[derive(Debug, Parser)]
pub struct LiftoverArgs {
    /// Source reference genome (required unless --list-chains)
    #[arg(
        short = 's',
        long,
        value_name = "FASTA",
        required_unless_present = "list_chains"
    )]
    pub source_ref: Option<PathBuf>,

    /// Target reference genome (optional; when omitted REF alleles are not
    /// validated against the target build after liftover)
    #[arg(short = 't', long, value_name = "FASTA")]
    pub target_ref: Option<PathBuf>,

    /// Chain file (required unless --list-chains)
    #[arg(
        short = 'c',
        long,
        value_name = "FILE",
        required_unless_present = "list_chains"
    )]
    pub chain: Option<PathBuf>,

    /// Write rejected variants here
    #[arg(short = 'r', long, value_name = "FILE")]
    pub reject: Option<PathBuf>,

    /// Output file (default: stdout)
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Print known chain file URLs and exit
    #[arg(long)]
    pub list_chains: bool,

    /// Add INFO/SRC_CONTIG and INFO/SRC_POS to each mapped record
    #[arg(long)]
    pub write_src_coords: bool,

    /// Do not flip alleles when a chain block is on the opposite strand;
    /// reject such records instead.
    #[arg(long)]
    pub no_fix_swapped_ref: bool,

    /// Suppress the progress bar even when stderr is a TTY
    #[arg(long)]
    pub no_progress: bool,

    /// Input VCF/BCF (default: stdin)
    #[arg(value_name = "INPUT")]
    pub input: Option<PathBuf>,
}

// ---------------------------------------------------------------------------
// filter
// ---------------------------------------------------------------------------

#[derive(Debug, Parser)]
pub struct FilterArgs {
    /// Filter expression (required)
    #[arg(short = 'e', long, value_name = "EXPR", required = true)]
    pub expression: String,

    /// Keep variants NOT matching the expression
    #[arg(long)]
    pub invert: bool,

    /// Output file (default: stdout)
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Suppress the progress bar even when stderr is a TTY
    #[arg(long)]
    pub no_progress: bool,

    /// Input VCF/BCF (default: stdin)
    #[arg(value_name = "INPUT")]
    pub input: Option<PathBuf>,
}

// ---------------------------------------------------------------------------
// completions
// ---------------------------------------------------------------------------

#[derive(Debug, Parser)]
pub struct CompletionsArgs {
    /// Shell to generate completions for
    #[arg(value_name = "SHELL")]
    pub shell: Shell,
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() -> ExitCode {
    let cli = Cli::parse();

    // Set up tracing based on --verbose count (--quiet forces warn level)
    let level = if cli.quiet {
        "warn"
    } else {
        match cli.verbose {
            0 => "warn",
            1 => "info",
            2 => "debug",
            _ => "trace",
        }
    };

    fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(level)),
        )
        .with_writer(io::stderr)
        .init();

    // Resolve telemetry setting once, up front. We only prompt for actual
    // pipeline commands (not `completions` — that's usually automated).
    let telemetry_enabled = if cli.no_telemetry {
        false
    } else {
        match &cli.command {
            Commands::Completions(_) => false,
            _ => {
                let mut cfg = telemetry::TelemetryConfig::load();
                cfg.ensure_prompted()
            }
        }
    };

    let command_name = match &cli.command {
        Commands::Normalize(_) => "normalize",
        Commands::Liftover(_) => "liftover",
        Commands::Filter(_) => "filter",
        Commands::Completions(_) => "completions",
    };

    let started = Instant::now();
    let result: Result<()> = match cli.command {
        Commands::Normalize(args) => commands::normalize::run(&args, cli.quiet),
        Commands::Liftover(args) => commands::liftover::run(&args, cli.quiet),
        Commands::Filter(args) => commands::filter::run(&args, cli.quiet),
        Commands::Completions(args) => {
            let mut cmd = Cli::command();
            let name = cmd.get_name().to_string();
            generate(args.shell, &mut cmd, name, &mut io::stdout());
            Ok(())
        }
    };

    let success = result.is_ok();
    let duration = started.elapsed();

    if telemetry_enabled {
        telemetry::send_event(command_name, duration, success);
    }

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("Error: {:#}", e);
            ExitCode::FAILURE
        }
    }
}
