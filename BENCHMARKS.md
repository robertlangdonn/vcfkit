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

## Performance Target

vcfkit aims to stay within **1.5×** of `bcftools` throughput for equivalent operations.

## E2E Results (placeholder — run the script to fill)

### normalize

| Command | Mean time |
|---------|-----------|
| `vcfkit normalize` | _not yet measured_ |
| `bcftools norm` | _not yet measured_ |

### filter

| Command | Mean time |
|---------|-----------|
| `vcfkit filter` | _not yet measured_ |
| `bcftools view -i` | _not yet measured_ |

---

*Regenerate: `bash benches/e2e/run.sh`*
