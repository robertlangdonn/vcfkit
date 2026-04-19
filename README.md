# vcfkit

Fast VCF toolkit for bioinformaticians. Three operations every pipeline needs — normalize, liftover, filter — as a single static binary with zero dependencies.

On 1000 Genomes chr22 (1.1M variants), vcfkit's fast paths beat bcftools by ~4×:

| Operation | vcfkit | bcftools | speedup |
|-----------|--------|----------|---------|
| `filter -e 'INFO/AF < 0.01'` | **422 ms** | 1,695 ms | **4.0×** |
| `normalize --fast --no-split` | **682 ms** | 2,820 ms | **4.1×** |

Measured 2026-04-19 on macOS aarch64 with bcftools 1.23.1. See [BENCHMARKS.md](BENCHMARKS.md) for full methodology.

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

## Credits and prior art

vcfkit exists because of decades of work by others:

**[htslib](https://github.com/samtools/htslib)** and **[bcftools](https://github.com/samtools/bcftools)** — the reference implementations for VCF/BCF processing. Created and maintained by the Wellcome Sanger Institute; primary authorship by Heng Li (original author), Petr Danecek (bcftools lead), and hundreds of contributors over 15+ years. vcfkit's differential tests validate against bcftools output — if vcfkit and bcftools diverge, vcfkit is wrong by default.

**[noodles](https://github.com/zaeleus/noodles)** by Michael Macias — the pure-Rust VCF, BCF, FASTA, and chain file I/O primitives that vcfkit builds on. Without noodles this project would not exist in its current form.

**[Tan, Abecasis, Kang 2015](https://doi.org/10.1093/bioinformatics/btv112)** — "Unified representation of genetic variants," *Bioinformatics* 31(13):2202–2204. The normalization algorithm implemented in `vcfkit normalize`.

**[UCSC Genome Browser](https://genome.ucsc.edu/)** — chain files and reference FASTAs used by `vcfkit liftover`.

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

Equivalent to `bcftools norm -f ref.fa -m-any -c w`.

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

## License

MIT
