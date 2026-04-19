# Changelog

## [0.1.5] — 2026-04-19

### Fixed

- **BCF output: error instead of silent VCF fallback** — `vcfkit normalize/filter/liftover -o out.bcf` previously wrote a VCF file named `.bcf` without any warning. Now returns a clear error with a bcftools workaround. BCF write support is planned for v0.2.

## [0.1.4] — 2026-04-19

### Fixed

- **normalize: multi-allelic indel corruption** — `left_align_record` was applying the biallelic left-alignment algorithm using only the first ALT allele, then overwriting all ALTs with that single result. For a record like `GA → GAA,G`, vcfkit was producing `G → GA` (losing the deletion allele entirely). Fix: skip left-alignment for multi-allelic records; bcftools norm without `-m` also passes these through unchanged. Found and confirmed by the real-world differential test harness against 1000 Genomes chr22 (875 affected records out of 1.1M).

### Changed

- Differential test harness: corrected bcftools norm command from `-m +any` (merge) to no `-m` flag (left-align only, matching vcfkit `--no-split` semantics). The `+any` flag was merging co-located records into multi-allelics, making the comparison apples-to-oranges.

## [0.1.3] — 2026-04-19

### Added

- **normalize `--fast` flag** — Raw-line fast path for biallelic SNPs/MNPs (≈80% of typical VCFs): reads raw bytes, skips noodles parse, falls back to full pipeline for multi-allelics and indels. Measured **4.1× faster** than `bcftools norm` on 1000G chr22 (682ms vs 2,820ms). Use with `vcfkit normalize --fast`.
- **liftover: b37/UCSC contig name mismatch detection** — Detects when the VCF uses b37 names ("1", "2", ...) but the chain uses UCSC names ("chr1", "chr2", ...) and errors with a clear `bcftools annotate --rename-chrs` hint. Use `--allow-contig-mismatch` to suppress and process anyway.
- **Real-world differential test harness** — `VCFKIT_REAL_TESTS=1 cargo test --test real_world_differential` runs 5 bcftools-vs-vcfkit comparisons on 1000G chr22. Gated behind an env var so it skips in normal CI. Nightly GitHub Actions workflow at `.github/workflows/nightly-real-tests.yml`.

### Performance

| Operation | vcfkit | bcftools | speedup |
|-----------|--------|----------|---------|
| `filter -e 'INFO/AF < 0.01'` | **422 ms** | 1,695 ms | **4.0×** |
| `normalize --fast --no-split` | **682 ms** | 2,820 ms | **4.1×** |

Measured on 1000G chr22, 1.1M variants, macOS aarch64, bcftools 1.23.1. See BENCHMARKS.md.

## [0.1.2] — 2026-04-18

### Fixed

- **liftover: gzip chain file support** — All UCSC chain files are `.gz`. v0.1.1 read them as raw bytes, silently producing zero lifted variants. Now auto-detects `.gz` extension and decompresses with `flate2`. This fixes liftover on every real-world UCSC chain file.
- **normalize: FASTA cache miss bug** — When the VCF uses b37 contig names ("22") but the FASTA uses UCSC names ("chr22"), v0.1.1 called `build_from_path()` on every record for the missing contig, causing 87× excess file I/O. Fixed with a `HashSet<absent_contigs>` — one syscall per absent contig, then fast path forever.
- **normalize: clearer .fai missing error** — Error message now includes `hint: create it with 'samtools faidx ref.fa'`.
- **normalize: contig-skip warning** — First time a contig is skipped (not in reference FASTA), a `WARN` is emitted: `contig "22" not found in reference FASTA — skipping normalization for this contig`. Subsequent records on the same contig are silent.
- **liftover: `--target-ref` now optional** — v0.1.1 required `--target-ref` even though it's only used for REF validation of lifted variants. It now warns when absent rather than failing.
- **filter: multi-allelic INFO fields** — `INFO/AF < 0.01` on `AF=0.09,0.003` now matches (any-element semantics), matching bcftools behavior. v0.1.1 only checked the first value.

### Added

- 5 new normalize tests covering `Number=G`, `Number=.`, `Flag` fields, `FORMAT Number=A`, `FORMAT Number=R` (all pass through the split pipeline correctly).

## [0.1.1] — 2026-04-12

- Cross-platform binary releases via cargo-dist (Linux x86_64, macOS arm64/x86_64, Windows x86_64)
- Published to crates.io as `vcfkit-cli` + `vcfkit-core`
- Filter fast path: 4× faster than bcftools on 1000 Genomes chr22 (390ms vs 1,635ms)
- Workspace metadata: homepage, keywords, categories

## [0.1.0] — 2026-04-10

- Initial release: `normalize`, `liftover`, `filter` operations
- Single static binary, no htslib, no C dependencies
- stdin/stdout piping, progress bars, shell completions, opt-in telemetry
