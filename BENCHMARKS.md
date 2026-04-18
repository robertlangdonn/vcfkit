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

## How the fast path works

The filter hot loop bypasses noodles record parsing entirely:

1. **Read raw lines** with `BufRead::read_line` into a reusable `String` buffer — no per-record allocation.
2. **Lazy field access** via `FastRecord`, which holds `&str` slices into the line buffer. Only the fields referenced by the expression are parsed.
3. **Pass-through writes** — matching records are written as raw bytes without re-serialization through noodles.
4. **INFO metadata** from the header (field types, Number=A/R/G) is extracted once before the loop and reused.

bcftools (htslib) parses VCF lazily in C with field offsets. vcfkit's fast path does the same in Rust.

Correctness is validated: both tools produce identical variant counts on the same input (910,425 of 1,103,547 variants with AF < 0.01 on chr22).

---

*Regenerate: `bash benches/e2e/run.sh`*
