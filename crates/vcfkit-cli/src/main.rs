mod commands;
mod output;
mod telemetry;

use std::io;
use std::path::PathBuf;

use anyhow::Result;
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use clap_complete::{generate, Shell};
use tracing_subscriber::{EnvFilter, fmt};

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
    #[arg(short, long, global = true)]
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
    /// Source reference genome (required)
    #[arg(short = 's', long, value_name = "FASTA", required = true)]
    pub source_ref: PathBuf,

    /// Target reference genome (required)
    #[arg(short = 't', long, value_name = "FASTA", required = true)]
    pub target_ref: PathBuf,

    /// Chain file (required)
    #[arg(short = 'c', long, value_name = "FILE", required = true)]
    pub chain: PathBuf,

    /// Write rejected variants here
    #[arg(short = 'r', long, value_name = "FILE")]
    pub reject: Option<PathBuf>,

    /// Output file (default: stdout)
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Print known chain file URLs and exit
    #[arg(long)]
    pub list_chains: bool,

    /// Input VCF/BCF (default: stdin)
    #[arg(value_name = "INPUT")]
    pub input: Option<PathBuf>,
}

// ---------------------------------------------------------------------------
// filter
// ---------------------------------------------------------------------------

#[derive(Debug, Parser)]
pub struct FilterArgs {
    /// Filter expression (required unless --list-fields)
    #[arg(short = 'e', long, value_name = "EXPR")]
    pub expression: Option<String>,

    /// Keep variants NOT matching the expression
    #[arg(short = 'v', long)]
    pub invert: bool,

    /// Output file (default: stdout)
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

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

fn main() -> Result<()> {
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

    match cli.command {
        Commands::Normalize(_args) => {
            eprintln!("normalize: not yet implemented");
        }
        Commands::Liftover(args) => {
            if args.list_chains {
                eprintln!("liftover --list-chains: not yet implemented");
            } else {
                eprintln!("liftover: not yet implemented");
            }
        }
        Commands::Filter(_args) => {
            eprintln!("filter: not yet implemented");
        }
        Commands::Completions(args) => {
            let mut cmd = Cli::command();
            let name = cmd.get_name().to_string();
            generate(args.shell, &mut cmd, name, &mut io::stdout());
        }
    }

    Ok(())
}
