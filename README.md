# vcfkit

Fast VCF toolkit for bioinformaticians. Three operations every pipeline needs — normalize, liftover, filter — as a single static binary with zero dependencies.

> **Status:** v0.3.0-alpha.2 — early preview. Suitable for research pipelines and evaluation. Not validated for production clinical use. Known behavioral differences from bcftools are documented in [docs/known_differences.md](docs/known_differences.md).

On 1000 Genomes chr22 (1.1M variants), vcfkit's fast paths beat bcftools by ~4×:

| Operation | vcfkit | bcftools | speedup |
|-----------|--------|----------|---------|
| `filter -e 'INFO/AF < 0.01'` | **422 ms** | 1,695 ms | **4.0×** |
| `normalize --fast --no-split` | **682 ms** | 2,820 ms | **4.1×** |
| `normalize` (standard, with noodles) | 6,481 ms | 2,820 ms | 0.43× |
| `liftover` | 6,713 ms | — | ~164K rec/s |

Standard normalize (without `--fast`) is currently **2.3× slower** than bcftools — the fast path is opt-in because it only handles biallelic SNPs/MNPs. Liftover has no bcftools equivalent on this platform (`bcftools +liftover` plugin unavailable). Measured 2026-04-19 on macOS aarch64 with bcftools 1.23.1. See [BENCHMARKS.md](BENCHMARKS.md) for full methodology.

```
vcfkit normalize -f ref.fa input.vcf > normalized.vcf
vcfkit liftover -s hg19.fa -t hg38.fa -c hg19ToHg38.over.chain.gz input.vcf > lifted.vcf
vcfkit filter -e "INFO/AF < 0.01 && FILTER == 'PASS'" input.vcf > rare_variants.vcf
```

## Install

```bash
cargo install vcfkit-cli   # installs the `vcfkit` binary
```

Or download a pre-built binary from [Releases](https://github.com/robertlangdonn/vcfkit/releases).  
No Rust required for the pre-built binary.

## Credits

vcfkit exists because of decades of work by others:

**[htslib](https://github.com/samtools/htslib)** and **[bcftools](https://github.com/samtools/bcftools)** — the reference implementations for VCF/BCF processing. Created and maintained by the Wellcome Sanger Institute; primary authorship by Heng Li (original author), Petr Danecek (bcftools lead), and hundreds of contributors over 15+ years. vcfkit's normalization behavior was developed by reading bcftools source. Differential tests validate against bcftools output — if vcfkit and bcftools diverge, vcfkit is wrong by default.

**[noodles](https://github.com/zaeleus/noodles)** by Michael Macias — the pure-Rust VCF, BCF, FASTA, and chain file I/O primitives that vcfkit builds on. Without noodles this project would not exist in its current form.

**[Tan, Abecasis, Kang 2015](https://doi.org/10.1093/bioinformatics/btv112)** — "Unified representation of genetic variants," *Bioinformatics* 31(13):2202–2204. The normalization algorithm implemented in `vcfkit normalize`.

**[UCSC Genome Browser](https://genome.ucsc.edu/)** — chain files and reference FASTAs used by `vcfkit liftover`.

### AI assistance

Portions of this codebase were written with assistance from Claude (Anthropic). The AI wrote code; a human verified correctness and owns the result. See [CREDITS.md](CREDITS.md) for full attribution details.

vcfkit's contribution: a modern CLI UX, single-binary distribution, measured performance improvements on specific hot paths via raw-line parsing (the same approach htslib uses in C, applied in Rust), and — in future releases — WASM and natural-language filter queries. It does not replace bcftools. See [BENCHMARKS.md](BENCHMARKS.md) for methodology.

## Usage

### normalize

Left-align indels and split multi-allelic sites:

```bash
vcfkit normalize -f reference.fasta input.vcf -o output.vcf

# Options
-f, --reference <FASTA>     Reference genome (required)
-o, --output <FILE>         Output file (default: stdout)
    --no-split              Keep multi-allelic sites
    --no-left-align         Skip left-alignment
    --check-ref <MODE>      ignore | warn | error  (default: warn)
```

With `--no-split` (multi-allelic records preserved), equivalent to `bcftools norm -f ref.fa -c w`. Without `--no-split`, equivalent to `bcftools norm -f ref.fa -m-any -c w` for biallelic records; multi-allelic indels are currently passed through without left-alignment (see [known differences](docs/known_differences.md)).

### liftover

Convert variants between genome builds using UCSC chain files:

```bash
vcfkit liftover \
  -s hg19.fa -t hg38.fa \
  -c hg19ToHg38.over.chain.gz \
  input.vcf -o output.vcf

# Known chain file URLs
vcfkit liftover --list-chains
```

Unmapped variants are written to a reject file (`-r rejects.vcf`) instead of silently dropped.

### filter

Keep variants matching an expression:

```bash
vcfkit filter -e "INFO/AF < 0.01" input.vcf
vcfkit filter -e "CHROM == 'chr17' && QUAL > 30" input.vcf
vcfkit filter -e "INFO/CSQ ~ 'missense'" input.vcf
vcfkit filter -e "FILTER == 'PASS'" --invert input.vcf   # keep non-PASS
```

Supported fields: `INFO/*`, `FORMAT/*`, `CHROM`, `POS`, `QUAL`, `FILTER`  
Operators: `<` `<=` `>` `>=` `==` `!=` `&&` `||` `!` `~` `!~`

Multi-allelic INFO fields (e.g. `AF=0.12,0.003`) use any-element semantics — the filter matches if any value satisfies the condition.

#### Natural-language filter (`--ask`)

Translate plain English to a filter expression via Anthropic's Claude API:

```bash
export ANTHROPIC_API_KEY=sk-ant-...

vcfkit filter --ask "rare PASS variants" input.vcf
vcfkit filter --ask "rare PASS variants" --yes input.vcf   # skip confirmation
vcfkit filter -a "rare PASS variants" --yes input.vcf      # short flag
```

The LLM sees only the VCF header schema (field names, types, descriptions) and your query. Variant data never leaves your machine. The translated expression is shown for review before filtering runs, unless `--yes` is passed. When confidence is below 50%, `--yes` is blocked; add `--accept-low-confidence` to override.

## Piping

All three operations read from stdin and write to stdout by default:

```bash
cat input.vcf \
  | vcfkit normalize -f ref.fa \
  | vcfkit filter -e "QUAL > 30" \
  > output.vcf
```

## Performance

The filter fast path reads raw VCF lines and only parses fields the expression references — matching records are written as raw bytes without re-serialization. See [BENCHMARKS.md](BENCHMARKS.md) for full results and methodology.

## Building from source

```bash
git clone https://github.com/robertlangdonn/vcfkit
cd vcfkit
cargo build --release
./target/release/vcfkit --version
```

Requires Rust 1.75+. No C dependencies, no htslib.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for release history.

## License

MIT
