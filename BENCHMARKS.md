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
| `vcfkit filter -e 'INFO/AF < 0.01'` | 4.38 s | 2.6× slower |
| `bcftools view -i 'INFO/AF < 0.01'` | 1.69 s | — |

### filter (`FILTER == 'PASS'`)

| Command | Mean time | vs bcftools |
|---------|-----------|-------------|
| `vcfkit filter -e "FILTER == 'PASS'"` | 4.87 s | 3.0× slower |
| `bcftools view -f PASS` | 1.62 s | — |

## Performance Notes

Current throughput on plain VCF is ~2.6–3× slower than bcftools. The bottleneck is
`noodles::vcf::io::Reader::read_record_buf`, which parses every field into owned Rust
types on each record. bcftools/htslib parses lazily and reuses scratch buffers.

Planned fix (Phase 2): switch to noodles lazy record reader so only the fields touched
by the expression are parsed. This should close the gap substantially.

Correctness is verified: both tools produce identical variant counts on the same input.

---

*Regenerate: `bash benches/e2e/run.sh`*
