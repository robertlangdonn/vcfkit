# Changelog

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
