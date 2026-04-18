# vcfkit — Alpha Testing

Thanks for taking a look. This document is everything you need to try vcfkit and tell me what broke.

## What it is

vcfkit is a fast command-line toolkit for three VCF operations bioinformaticians run constantly:

- **normalize** — left-align indels, split multi-allelic sites
- **liftover** — convert variants between genome builds (hg19 ↔ hg38 ↔ T2T-CHM13)
- **filter** — keep variants matching an expression (`INFO/AF < 0.01 && FILTER == 'PASS'`)

It's a single static binary, no dependencies, no htslib, no Docker. Installs in one command.

## Install

```bash
# Requires Rust 1.75+
cargo install vcfkit-cli   # installs the `vcfkit` binary

# Verify
vcfkit --version
```

If you don't have Rust: `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`

## Try it

```bash
# normalize — left-align indels and split multi-allelics
vcfkit normalize -f reference.fasta input.vcf > normalized.vcf

# liftover — convert hg19 → hg38
vcfkit liftover --list-chains   # shows known chain file URLs
vcfkit liftover -s hg19.fa -t hg38.fa -c hg19ToHg38.over.chain.gz input.vcf > lifted.vcf

# filter — expression-based variant selection
vcfkit filter -e "INFO/AF < 0.01" input.vcf
vcfkit filter -e "CHROM == 'chr17' && QUAL > 30 && FILTER == 'PASS'" input.vcf
vcfkit filter -e "INFO/CSQ ~ 'missense'" input.vcf

# all three pipe cleanly
cat input.vcf | vcfkit normalize -f ref.fa | vcfkit filter -e "QUAL > 30" > output.vcf
```

## What's working

- All three operations on real VCFs (tested on synthetic data; needs more real-world validation)
- stdin/stdout piping
- Progress bars with rate + ETA (hidden automatically when piping)
- Helpful error messages with file/line context
- Shell completions: `vcfkit completions bash > vcfkit.bash`
- Opt-in anonymous telemetry (off by default, asks on first run)

## Known issues / not yet done

- BCF output not supported yet — pipe through `bcftools view -O b` to convert
- No bgzip output (`.vcf.gz`) yet
- LLM natural-language filter coming in a later phase (`vcfkit filter -e "rare coding variants on chr17"`)
- `vcfkit filter` benchmarks at 4× faster than bcftools on 1000G chr22 (1.1M variants)

## How to give feedback

Open a GitHub issue: https://github.com/robertlangdonn/vcfkit/issues

Or email directly. Most useful things to report:
- It crashed (include the command, input file type/size, and the error message)
- Output differs from bcftools (include the bcftools command you'd use for the same operation)
- A flag or behavior you expected to exist that doesn't
- Performance — faster or slower than your current tool?
