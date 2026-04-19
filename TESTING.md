# Testing strategy

vcfkit's correctness is validated at five levels.

## 1. Unit and integration tests

`cargo test --workspace` — 100+ tests covering individual functions, edge cases,
INFO field Number= types, header preservation, and CLI end-to-end behaviour.
Runs in CI on every push.

## 2. Synthetic differential tests vs bcftools

Five `#[ignore]` tests that invoke both vcfkit and bcftools on synthetic VCFs
and compare output. Covers Number=1, Number=A, Number=R, Number=G, and Number=.
INFO field splitting. Requires `bcftools` on PATH.

```bash
cargo test --workspace -- --ignored --include-ignored
```

## 3. Real-world differential tests vs bcftools

Five tests that run both tools on 1000 Genomes chr22 (~1.1 M variants) and
semantically diff the output. Downloads the data on first run; subsequent runs
use the cached copy. Gated behind `VCFKIT_REAL_TESTS=1`.

```bash
VCFKIT_REAL_TESTS=1 cargo test --release --test real_world_differential -- --nocapture --test-threads=1
```

Runs nightly in CI via `.github/workflows/nightly-real-tests.yml`.

## 4. WASM parity tests

Runs the WASM build under Node.js and compares output against expected files
generated from the native CLI. Currently 11 tests:

| Category | Tests |
|----------|-------|
| filter | QUAL threshold, CHROM equality, INFO/AF any-element, compound &&, FILTER field |
| normalize | multi-allelic split, SNV passthrough, Number=G PL field |
| liftover | identity chain, partial coverage, position offset |

```bash
# Build WASM first
wasm-pack build crates/vcfkit-core --target web --out-dir ../../web/public/wasm --release

# Run from repo root
node tests/wasm_runtime/run.mjs
```

Runs in CI on every push. CI also verifies that the committed WASM artifacts
match a fresh build — stale WASM fails the `wasm` job.

## 5. WASM in-process tests

`crates/vcfkit-core/tests/wasm_parity.rs` — six Rust integration tests that
exercise the same code paths as the WASM bindings using in-process I/O
(no WASM toolchain needed). Runs in `cargo test --workspace`.

## Running everything

```bash
# Levels 1 + 2 + 5
cargo test --workspace -- --ignored --include-ignored

# Level 3 (requires bcftools, downloads ~150 MB on first run)
VCFKIT_REAL_TESTS=1 cargo test --release --test real_world_differential -- --nocapture --test-threads=1

# Level 4 (requires a WASM build in web/public/wasm/)
node tests/wasm_runtime/run.mjs
```
