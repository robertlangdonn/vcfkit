# vcfkit Benchmarks

## Running Microbenchmarks (Criterion)

```bash
# Run all Criterion microbenchmarks
cargo bench

# Compile benchmarks without running them (used in CI)
cargo bench --no-run

# Run a specific benchmark by name
cargo bench normalize_50_snp_records
cargo bench liftover_50_records
cargo bench filter_50_records_af_lt_0.05
```

HTML reports are written to `target/criterion/` and can be opened in a browser.

## Running the E2E Benchmark Script

```bash
bash benches/e2e/run.sh
```

Prerequisites: `bcftools` and `hyperfine` must be in `PATH`. If either is missing
the script exits cleanly with an install hint so it does not break CI.

Results are written to `benches/e2e/report.md`.

## E2E Results (1000 Genomes chr22, 1.1M variants, macOS aarch64)

Measured 2026-04-19 with bcftools 1.23.1 and hyperfine 1.20.0 (5 runs, 1 warmup).  
Input: `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
extracted to plain sites-only VCF (149 MB, 1,103,547 records).

### filter (`INFO/AF < 0.01`)

| Command | Mean time | vs bcftools |
|---------|-----------|-------------|
| `vcfkit filter -e 'INFO/AF < 0.01'` | **422 ms** | **4.0× faster** |
| `bcftools view -i 'INFO/AF < 0.01'` | 1,695 ms | — |

### filter (`FILTER == 'PASS'`)

| Command | Mean time | vs bcftools |
|---------|-----------|-------------|
| `vcfkit filter -e "FILTER == 'PASS'"` | **462 ms** | **3.5× faster** |
| `bcftools view -f PASS` | 1,610 ms | — |

### normalize

| Command | Mean time | vs bcftools |
|---------|-----------|-------------|
| `vcfkit normalize --fast --no-split` | **682 ms** | **4.1× faster** |
| `bcftools norm -m +any` | 2,820 ms | — |
| `vcfkit normalize --no-split` (noodles path) | 6,481 ms | 2.3× slower |

The standard noodles path is slower because `noodles` allocates a `BTreeMap + Vec + String` per record. The `--fast` path reads raw bytes for biallelic SNPs/MNPs (≈80% of chr22) and falls back to noodles only for multi-allelics and indels.

### liftover

| Command | Mean time | throughput |
|---------|-----------|------------|
| `vcfkit liftover` (hg19 → hg38) | 6,713 ms | ~164K rec/s |

`bcftools +liftover` is a plugin that requires manual compilation and was not available for comparison. The chain lookup is O(log n) per record (binary search into a sorted `Vec<ChainBlock>` per source contig).

## How the fast paths work

### filter fast path (4.0× measured speedup)

The filter hot loop bypasses noodles record parsing entirely:

1. **Read raw lines** with `BufRead::read_line` into a reusable `String` buffer — no per-record allocation.
2. **Lazy field access** via `FastRecord`, which holds `&str` slices into the line buffer. Only the fields referenced by the expression are parsed.
3. **Pass-through writes** — matching records are written as raw bytes without re-serialization through noodles.
4. **INFO metadata** from the header (field types, Number=A/R/G) is extracted once before the loop and reused.

### normalize fast path (`--fast`, 4.1× measured speedup)

Same approach as filter for biallelic SNPs/MNPs (≈80% of typical VCFs):

1. Detect SNP/MNP from raw tab-split columns — no noodles parse.
2. Optional REF check: single byte lookup into cached contig sequence.
3. Write raw bytes unchanged.
4. Multi-allelic records and indels (when left-align is enabled) automatically fall back to the full noodles pipeline.

Correctness is validated by parity tests: fast mode and noodles mode produce identical output on all corpus VCFs.

bcftools (htslib) parses VCF lazily in C with field offsets. vcfkit's fast paths do the same in Rust.

---

*Regenerate: `bash benches/e2e/run.sh`*
