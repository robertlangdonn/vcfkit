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

Measured 2026-04-18 with bcftools 1.23.1 and hyperfine 1.20.0.  
Input: `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
extracted to plain VCF (sites-only, 149 MB, 1,103,547 records).

### filter (`INFO/AF < 0.01`)

| Command | Mean time | vs bcftools |
|---------|-----------|-------------|
| `vcfkit filter -e 'INFO/AF < 0.01'` | **390 ms** | **4.2× faster** |
| `bcftools view -i 'INFO/AF < 0.01'` | 1,635 ms | — |

### filter (`FILTER == 'PASS'`)

| Command | Mean time | vs bcftools |
|---------|-----------|-------------|
| `vcfkit filter -e "FILTER == 'PASS'"` | **462 ms** | **3.5× faster** |
| `bcftools view -f PASS` | 1,610 ms | — |

### normalize (pending — v0.1.3 target)

`vcfkit normalize --fast` was added in the post-v0.1.2 work. Benchmark numbers for
normalize vs `bcftools norm` on chr22 have not yet been collected.

Expected: ~4× speedup for SNP-heavy VCFs (chr22 is ~80% SNPs), similar to filter,
since both use the same raw-line + lazy-parse approach.

### liftover (pending)

Liftover benchmarks require a real hg19→hg38 chain file and a chr22 hg19 VCF.
The chain lookup is O(log n) by construction (binary search into a sorted `Vec<ChainBlock>`
per source contig). End-to-end timing vs `bcftools +liftover` will be added in v0.1.3.

## How the fast paths work

### filter fast path (active, 4× measured speedup)

The filter hot loop bypasses noodles record parsing entirely:

1. **Read raw lines** with `BufRead::read_line` into a reusable `String` buffer — no per-record allocation.
2. **Lazy field access** via `FastRecord`, which holds `&str` slices into the line buffer. Only the fields referenced by the expression are parsed.
3. **Pass-through writes** — matching records are written as raw bytes without re-serialization through noodles.
4. **INFO metadata** from the header (field types, Number=A/R/G) is extracted once before the loop and reused.

### normalize fast path (active via `--fast`, benchmarks pending)

Same approach as filter for biallelic SNPs/MNPs (≈80% of typical VCFs):

1. Detect SNP/MNP from raw tab-split columns — no noodles parse.
2. Optional REF check: single byte lookup into cached contig sequence.
3. Write raw bytes unchanged.
4. Multi-allelic records and indels (when left-align is enabled) automatically fall back to the full noodles pipeline.

Correctness is validated by parity tests: fast mode and noodles mode produce identical output on all corpus VCFs.

bcftools (htslib) parses VCF lazily in C with field offsets. vcfkit's fast paths do the same in Rust.

---

*Regenerate: `bash benches/e2e/run.sh`*
