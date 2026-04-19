# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

**vcfkit** is a fast, pure-Rust CLI toolkit for bioinformaticians performing three VCF operations:
- `normalize` — Left-align indels, split multi-allelic sites, validate against reference genome
- `liftover` — Convert variants between genome builds (hg19 ↔ hg38 ↔ T2T-CHM13)
- `filter` — Select variants using deterministic expressions or natural-language via `--ask`

The positioning: what `ripgrep` is to `grep` — faster, better UX, single static binary, zero dependencies.

## Commands

```bash
cargo build --release
cargo install --path crates/vcfkit-cli   # install locally
cargo test --workspace                   # all tests (skips #[ignore] by default)
cargo test --workspace -- --ignored --include-ignored  # include differential tests (needs bcftools)
cargo clippy --workspace -- -D warnings
cargo fmt --check
cargo bench                              # Criterion microbenchmarks
bash benches/e2e/run.sh                  # E2E vs bcftools (needs bcftools + hyperfine)

# Run a single test by name
cargo test -p vcfkit-cli normalize_splits_multiallelic

# Run the ask_gate integration tests (no API key needed — uses mock)
cargo test -p vcfkit-cli --test ask_gate

# Run the eval harness against the real Anthropic API (costs money)
ANTHROPIC_API_KEY=sk-ant-... VCFKIT_EVAL_CONFIRM=1 cargo run --bin ask_eval

# WASM build (outputs to web/public/wasm/)
wasm-pack build crates/vcfkit-core --target web --out-dir ../../web/public/wasm --release

# WASM parity tests (run from repo root, after WASM build)
node tests/wasm_runtime/run.mjs

# Web dev server (site lives at vcfkit.dev on Cloudflare Pages)
cd web && npm install && npm run dev
```

**No Makefile or justfile** — pure `cargo` commands only.

## Architecture

Two-crate Rust workspace. The split keeps `vcfkit-core` usable from WASM and future bindings without pulling in CLI concerns.

```
crates/
  vcfkit-core/        # Library — pure logic, no CLI, compiles to both native and wasm32
    src/
      io.rs           # VCF/BCF read/write, format auto-detection, OutputFormat enum
      normalize.rs    # Left-align (Tan 2015 algorithm), multi-allelic split (Number=A/R/G/.), ref check
      liftover.rs     # UCSC chain file parser, strand handling, b37/UCSC contig mismatch detection
      filter.rs       # nom expression parser + evaluator (INFO/*, FORMAT/*, CHROM, POS, QUAL, FILTER)
                      # also contains extract_header_schema() — builds LLM prompt schema from VCF header
      wasm.rs         # #[wasm_bindgen] wrappers: filter_vcf, normalize_vcf, liftover_vcf
      error.rs        # VcfkitError (thiserror)
  vcfkit-cli/         # Binary — clap subcommands, progress bars, telemetry
    src/
      commands/       # normalize.rs, liftover.rs, filter.rs — wire clap args to core functions
      ask.rs          # Anthropic API client: translate() async fn, AskError, HeaderSchema,
                      # build_system_prompt(), extract_json(), VCFKIT_MOCK_TRANSLATION mock path
      output.rs       # ProgressReporter (indicatif, auto-hides when output is piped)
      telemetry.rs    # Opt-in usage telemetry, config at ~/.config/vcfkit/config.toml
    bin/
      ask_eval.rs     # Eval harness: runs ask_corpus.yaml through the API, reports pass rate
    tests/
      ask_gate.rs     # Integration tests for --ask: mock path, confidence gate, flag exclusivity

web/                  # Astro 6 + Starlight 0.38 + React 18 docs site at vcfkit.dev
  src/
    components/       # Demo.tsx (main interactive demo), VcfEditor.tsx (CodeMirror 6),
                      # ResultPanel.tsx, CliEquivalent.tsx, ExamplePicker.tsx
    lib/
      wasm-loader.ts  # Singleton WASM init; import via ensureWasm()
      vcf-examples.ts # Built-in example VCFs per operation
      format-utils.ts # countRecords, truncateForDemo (10K record cap)
    content/docs/     # MDX pages: index, install, commands/*, benchmarks, known-differences, privacy, credits
  public/wasm/        # Committed WASM artifacts (vcfkit_core_bg.wasm + vcfkit_core.js)
                      # CI fails if these are stale relative to vcfkit-core source

tests/
  corpus/synthetic/   # Hand-crafted VCF edge cases + mini_ref.fa
  ask_corpus.yaml     # 25 eval cases for --ask translation accuracy (target: ≥85% pass rate)
  wasm_runtime/       # Node.js WASM parity tests (run.mjs, fixtures/, expected/)
  real_world/         # 1000G chr22 data (downloaded on first run, then cached)
```

## Key Constraints

**Read `_docs/01_tech_stack.md` before adding any dependency.** The stack is settled.

**Correctness is validated by differential tests against bcftools.** `#[ignore]` tests in `normalize_test.rs` and `liftover_test.rs` call bcftools directly when it's on PATH. Any divergence is a bug.

**`normalize --fast`** is a raw-byte fast path for biallelic SNPs/MNPs (~80% of typical VCFs). It bypasses noodles parsing. The non-fast path uses full noodles parsing. Both must produce identical output.

**BCF write is not supported.** Writing to a `.bcf` path returns a clear error with a bcftools workaround. BCF write is planned for v0.2.

**`liftover` needs `--source-ref` from the CLI but not from WASM.** `liftover_from_chain_reader` (used by WASM) accepts a `BufRead` chain and does no ref validation.

**normalize in WASM** always uses `left_align=false, check_ref=Ignore` — the ref path is never opened so `Path::new("")` is safe to pass.

**`--ask` lives entirely in `vcfkit-cli`, never in `vcfkit-core`.** The WASM build is unaffected. `--ask` requires an input file path — stdin is not supported because the header must be read separately from the variant stream. Not available in the browser demo.

**`--ask` confidence gate:** `--yes` is blocked when translation confidence < 50%. Add `--accept-low-confidence` to override. Integration tests use `VCFKIT_MOCK_TRANSLATION=<json>` to bypass the API entirely — set this env var to a `TranslationPayload` JSON string.

**`ask_eval.rs` uses `#[path = "../ask.rs"] mod ask;`** because `vcfkit-cli` has no `lib.rs`. This pattern must be preserved if `ask.rs` is ever renamed or moved.

## WASM artifacts

When you change `vcfkit-core`, rebuild and commit the WASM bundle:

```bash
wasm-pack build crates/vcfkit-core --target web --out-dir ../../web/public/wasm --release
git add web/public/wasm/
git commit -m "chore(wasm): rebuild bundle for <change>"
```

CI (`wasm` job) diffs the committed artifacts against a fresh build and fails if they're stale.

**Exception:** pure documentation changes, `web/` changes, or `vcfkit-cli`-only changes do not require a WASM rebuild — only changes to `crates/vcfkit-core/src/` do.

## Phase Status

- **Phase 1 (Core CLI):** ✅ Complete — normalize/liftover/filter, 100+ tests, Criterion benchmarks, real-world differential tests vs bcftools on 1000G chr22
- **Phase 2 (WASM + Web):** ✅ Complete — wasm-bindgen wrappers, Astro/Starlight site at vcfkit.dev on Cloudflare Pages, 11-test WASM parity suite
- **Phase 3 (LLM filter):** ✅ v0.3.0-alpha.2 shipped — `vcfkit filter --ask "<query>"` / `-a`, Anthropic Claude API (`ANTHROPIC_API_KEY`), confidence gate, `ask_eval` harness vs `tests/ask_corpus.yaml`. **Not published to crates.io yet** — pending dogfood validation.

## Testing levels

See `TESTING.md` for full details. Summary:

| Level | Command | Notes |
|-------|---------|-------|
| Unit + integration | `cargo test --workspace` | Always runs in CI |
| `--ask` gate (mock) | `cargo test -p vcfkit-cli --test ask_gate` | No API key needed |
| Synthetic differential vs bcftools | `cargo test -- --ignored --include-ignored` | Needs `bcftools` on PATH |
| Real-world differential (1000G chr22) | `VCFKIT_REAL_TESTS=1 cargo test --release --test real_world_differential` | Runs nightly in CI |
| WASM parity | `node tests/wasm_runtime/run.mjs` | Runs in CI after WASM build |
| LLM eval (real API) | `ANTHROPIC_API_KEY=... VCFKIT_EVAL_CONFIRM=1 cargo run --bin ask_eval` | Costs ~$0.025 per run |

## Dependencies (Settled)

| Crate | Purpose |
|-------|---------|
| `noodles` 0.109 | VCF/BCF/FASTA/chain I/O |
| `clap` v4 + `clap_complete` | CLI with derive macros + shell completions |
| `nom` v7 | Expression parser for `filter` |
| `anyhow` / `thiserror` | Error handling (CLI layer vs library) |
| `tracing` + `tracing-subscriber` | Structured logging |
| `owo-colors` | Coloured terminal output |
| `indicatif` | Progress bars with ETA |
| `serde` + `serde_json` | Telemetry serialisation |
| `serde_yaml` | YAML parsing for `ask_eval` corpus |
| `flate2` | gzip chain file decompression in liftover |
| `tokio` | Async runtime for Anthropic API calls in `--ask` |
| `reqwest` 0.12 (rustls-tls) | HTTP client for Anthropic API in `--ask` |
| `wasm-bindgen` 0.2 | WASM bindings (wasm32 target only) |
| `console_error_panic_hook` 0.1 | Better panic messages in browser (wasm32 only) |

No Python, no htslib, no C dependencies. Single static binary is a hard requirement.
