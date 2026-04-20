# Changelog

## [0.3.0-alpha.4] — 2026-04-20

### Fixed

- **Confidence gate boundary**: gate now fires at confidence ≤ 50%, not strictly
  < 50%. Previously the LLM could return exactly 50% for a known-compromised
  expression (e.g. the FORMAT/AD any-element trap) and pass `--yes` silently.
  `LOW_CONFIDENCE_THRESHOLD` is now 0.51; the user-facing error message still
  reads "below 50% threshold". Regression test added (`ask_yes_exactly_50_percent_is_blocked`).

- **Prompt: array-indexing workaround example**: added a second compromised-translation
  worked example to the system prompt showing the FORMAT/AD any-element case at
  confidence 0.25. Includes the explicit rule: "setting confidence to exactly 0.5
  is not a safe compromise — any expression that the caveats describe as 'matches
  records the user did not intend' must be explicitly below the gate."

  Found during v0.3.0-alpha.3 dogfood on GiAB HG001: "variants with exactly 20 alt
  reads" returned confidence 50% and matched 421,921 records instead of ~thousands.

---

## [0.3.0-alpha.3] — 2026-04-20

### Changed

- **`--ask` confidence calibration** — added explicit calibration rules to the
  system prompt so the LLM uses confidence levels consistently:
  - ≥ 0.8: expression fully answers the query
  - 0.5–0.8: answers with reasonable caveats
  - < 0.5: compromised expression (proxy field, over-matches, language limitation)

  Added a trigger rule: any caveat of the form "this may match records the user did
  not intend" or "a more complete expression would be X" forces confidence below 0.5,
  ensuring the gate fires rather than silently running a wrong answer.

- Added a worked example of a compromised translation ("biallelic SNPs") to the
  system prompt, teaching calibration by demonstration.

### Fixed

- Queries like "biallelic SNPs" and "variants with exactly 20 alt reads" previously
  returned ~55–60% confidence and ran silently under `--yes`. After this fix they
  should gate at < 50%, prompting the user to review before running.

### Known limitation tracked

- GitHub issue #2: `--ask` cannot express per-element array access (`FORMAT/AD[1]`,
  `INFO/AF[0]`). Queries needing alt-allele depth require a bcftools workaround.
  Planned for v0.4 filter expression language extension.

---

## [0.3.0-alpha.2] — 2026-04-19

### Changed

- **`--english` renamed to `--ask`** (short flag `-a`). The old name is removed.
- **Confidence gate**: `--ask --yes` on a translation with confidence below 50% now
  exits non-zero with an error. Add `--accept-low-confidence` to override.
- **Mock path for testing**: set `VCFKIT_MOCK_TRANSLATION=<json>` to bypass the API
  entirely (useful for CI and offline smoke-testing).
- Eval binary renamed from `english_eval` to `ask_eval`; corpus file renamed from
  `tests/english_filter_corpus.yaml` to `tests/ask_corpus.yaml`.

### Added

- `--accept-low-confidence` flag: proceed with a low-confidence translation when
  `--yes` is also set.
- 4 integration tests in `crates/vcfkit-cli/tests/ask_gate.rs` covering the mock
  path, confidence gate, and flag exclusivity.

---

## [0.3.0-alpha.1] — 2026-04-19

### Added

- **Natural-language filter queries** via `vcfkit filter --english "..."`.
  Translates plain English like `"rare missense variants"` into a deterministic
  filter expression via Anthropic's Claude API, shows the expression for
  confirmation, then runs it. `--yes` skips the prompt for scripting.

  ```
  $ vcfkit filter --english "rare PASS variants" input.vcf

    Query:      rare PASS variants
    Expression: INFO/AF < 0.01 && FILTER == 'PASS'
    Reasoning:  Rare is conventionally AF < 1%. FILTER == PASS ensures all
                caller filters passed.

  Run this filter? [Y/n/edit]
  ```

- **`english_eval` binary** — measures translation accuracy against
  `tests/english_filter_corpus.yaml` (25 cases). Target: ≥85% pass rate.
  Run with `ANTHROPIC_API_KEY=... VCFKIT_EVAL_CONFIRM=1 cargo run --bin english_eval`.

### Privacy

The LLM receives only the VCF header schema (INFO/FORMAT field names, types,
descriptions) and your query text. Variant data (CHROM, POS, REF, ALT,
genotypes) never leaves your machine. The translated expression runs through
the existing deterministic parser — the LLM cannot cause arbitrary behaviour.

### Requirements

- `--english` requires `ANTHROPIC_API_KEY` to be set.
- Default model: `claude-haiku-4-5`. Override with `VCFKIT_LLM_MODEL`.
- Default timeout: 30 s. Override with `VCFKIT_LLM_TIMEOUT_SECS`.
- `--english` requires an input file path (stdin not supported — the header
  must be readable without consuming the variant data stream).

### CLI-only

`--english` is not available in the WASM browser demo (API key cannot be
safely exposed client-side). The demo continues to support the `-e` expression
syntax for all three operations.

### No breaking changes to `-e` / existing operations.

---

## [0.1.6] — 2026-04-19

### Testing infrastructure

- WASM parity test suite expanded from 4 to 15 cases — filter now covers QUAL, CHROM,
  INFO/AF any-element semantics, compound `&&`, FILTER field, OR (`||`), negation (`!`),
  regex (`~`), and POS range; normalize covers multi-allelic split, SNV passthrough,
  and Number=G PL field; liftover covers identity, partial coverage, and position offset.
- CI verifies committed WASM artifacts in `web/public/wasm/` match a fresh `wasm-pack`
  build. Stale WASM fails the build with a rebuild hint.
- Nightly CI workflow runs differential tests against `bcftools` on 1000 Genomes chr22
  (1.1M variants) for all three operations. Any divergence fails the build within 24 hours.
- Adversarial poly-A tract multi-allelic test added with committed fixtures
  (`tests/corpus/synthetic/multiallelic_polyA.vcf`). Confirms the `--no-split` shortcut
  matches `bcftools norm` on a deliberately stressful input.

### Documentation

- `TESTING.md` — five-level validation strategy with run commands for each level.
- `docs/known_differences.md` — expanded with adversarial test reference, root cause
  analysis, and v0.2 fix plan. GitHub issue #1 tracks the v0.2 implementation.

### No runtime behaviour changes.

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
