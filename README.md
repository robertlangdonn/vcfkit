# vcfkit

Fast VCF toolkit for bioinformaticians — normalize, liftover, filter — as a single static binary with zero dependencies.

**[vcfkit.dev](https://vcfkit.dev)** · [Docs](https://vcfkit.dev/introduction) · [Install](https://vcfkit.dev/install) · [Changelog](CHANGELOG.md)

> **v0.3.0-alpha.4** — suitable for research pipelines and evaluation. Not validated for clinical use. [Known differences from bcftools →](docs/known_differences.md)

---

## What it does

Three operations every VCF pipeline needs, rewritten in Rust:

| Command | What it does |
|---------|-------------|
| `normalize` | Left-align indels, split multi-allelic sites, validate against reference FASTA |
| `liftover` | Convert between genome builds (hg19, hg38, T2T-CHM13) using UCSC chain files |
| `filter` | Keep variants matching expressions over INFO, FORMAT, CHROM, POS, QUAL, FILTER |

```bash
vcfkit normalize -f ref.fa input.vcf > normalized.vcf
vcfkit liftover -s hg19.fa -t hg38.fa -c hg19ToHg38.over.chain.gz input.vcf > lifted.vcf
vcfkit filter -e "INFO/AF < 0.01 && FILTER == 'PASS'" input.vcf > rare_variants.vcf
```

All three read from stdin and write to stdout by default — pipe them freely.

## Performance

Measured on 1000 Genomes chr22 (1,103,547 variants), macOS aarch64, bcftools 1.23.1:

| Operation | vcfkit | bcftools | Speedup |
|-----------|--------|----------|---------|
| `filter -e 'INFO/AF < 0.01'` | **422 ms** | 1,695 ms | **4.0×** |
| `normalize --fast --no-split` | **682 ms** | 2,820 ms | **4.1×** |
| `normalize` (standard path) | 6,481 ms | 2,820 ms | 0.43× |
| `liftover` | 6,713 ms | — | ~164K rec/s |

The fast path applies to biallelic SNPs/MNPs (~80% of typical VCFs). Standard normalize uses a pure-Rust parser — correct on all inputs but slower than bcftools' C implementation. See [BENCHMARKS.md](BENCHMARKS.md) for full methodology.

## Install

```bash
# Cargo
cargo install vcfkit-cli

# Pre-built binaries (macOS, Linux, Windows) — no Rust required
# See https://vcfkit.dev/install or GitHub Releases
```

Homebrew tap is planned. Follow the repo or watch [Releases](https://github.com/robertlangdonn/vcfkit/releases) for when it lands.

Full install instructions (including Windows, ARM, shell completions): **[vcfkit.dev/install](https://vcfkit.dev/install)**

## Usage

### filter

```bash
# Expression filter
vcfkit filter -e "INFO/AF < 0.01" input.vcf
vcfkit filter -e "QUAL > 30 && FILTER == 'PASS'" input.vcf
vcfkit filter -e "CHROM == 'chr17' && POS >= 43044295" input.vcf
vcfkit filter -e "INFO/CSQ ~ 'missense'" input.vcf
vcfkit filter -e "FILTER == 'PASS'" --invert input.vcf   # keep non-PASS

# Natural-language filter via Claude (requires ANTHROPIC_API_KEY)
vcfkit filter --ask "rare PASS variants on chromosome 17" input.vcf
vcfkit filter --ask "rare PASS variants" --yes input.vcf   # skip confirmation
```

Fields: `INFO/*`, `FORMAT/*`, `CHROM`, `POS`, `QUAL`, `FILTER`  
Operators: `<` `<=` `>` `>=` `==` `!=` `&&` `||` `!` `~` (contains) `!~`

With `--ask`: the LLM sees only your query and the VCF header schema. Variant data never leaves your machine. The translated expression is shown for review before filtering runs.

### normalize

```bash
vcfkit normalize -f hg38.fa input.vcf > normalized.vcf
vcfkit normalize --fast -f hg38.fa input.vcf > normalized.vcf   # 4× faster for SNPs/MNPs
vcfkit normalize -f hg38.fa --no-split input.vcf                 # keep multi-allelic
vcfkit normalize -f hg38.fa --check-ref error input.vcf          # strict REF validation
```

### liftover

```bash
vcfkit liftover -s hg19.fa -t hg38.fa -c hg19ToHg38.over.chain.gz input.vcf
vcfkit liftover ... -r rejects.vcf input.vcf   # keep unmapped records
```

Download chain files: `vcfkit liftover --list-chains`

## Try it in the browser

**[vcfkit.dev](https://vcfkit.dev)** — all three operations run in WebAssembly in your browser. Nothing uploads anywhere.

## Build from source

```bash
git clone https://github.com/robertlangdonn/vcfkit
cd vcfkit
cargo build --release
./target/release/vcfkit --version
```

Requires Rust 1.75+. No C dependencies, no htslib.

## Correctness

Validated against bcftools in differential tests on 1000 Genomes chr22 (1.1M real variants). Tests run nightly in CI. Known divergences are documented in [docs/known_differences.md](docs/known_differences.md).

## Credits

vcfkit builds on the work of:

- **[htslib](https://github.com/samtools/htslib) / [bcftools](https://github.com/samtools/bcftools)** — Wellcome Sanger Institute, Heng Li, Petr Danecek, and contributors. The reference implementation vcfkit validates against.
- **[noodles](https://github.com/zaeleus/noodles)** by Michael Macias — the pure-Rust VCF/BCF/FASTA I/O library vcfkit builds on.
- **[Tan, Abecasis, Kang 2015](https://doi.org/10.1093/bioinformatics/btv112)** — the normalization algorithm implemented in `vcfkit normalize`.
- **[UCSC Genome Browser](https://genome.ucsc.edu/)** — chain files for liftover.

Portions of this codebase were written with assistance from Claude (Anthropic). See [CREDITS.md](CREDITS.md) for full attribution.

## License

[MIT](LICENSE)
