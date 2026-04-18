# vcfkit

Fast VCF toolkit for bioinformaticians. Three operations every pipeline needs — normalize, liftover, filter — as a single static binary with zero dependencies.

```
vcfkit normalize -f ref.fa input.vcf > normalized.vcf
vcfkit liftover -s hg19.fa -t hg38.fa -c hg19ToHg38.over.chain.gz input.vcf > lifted.vcf
vcfkit filter -e "INFO/AF < 0.01 && FILTER == 'PASS'" input.vcf > rare_variants.vcf
```

## Install

```bash
cargo install vcfkit
```

Or download a pre-built binary from [Releases](https://github.com/robertlangdonn/vcfkit/releases).

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

Equivalent to `bcftools norm -f ref.fa -m-any -c w`, but faster.

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

## Piping

All three operations read from stdin and write to stdout by default:

```bash
cat input.vcf \
  | vcfkit normalize -f ref.fa \
  | vcfkit filter -e "QUAL > 30" \
  > output.vcf
```

## Performance

Target: ≥1.5× faster than bcftools on 1M-variant inputs. See [BENCHMARKS.md](BENCHMARKS.md).

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
