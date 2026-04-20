# Contributing to vcfkit

## Development setup

Requires Rust 1.75+. No C dependencies, no htslib.

```bash
git clone https://github.com/robertlangdonn/vcfkit
cd vcfkit
cargo build --release
cargo test --workspace
```

For differential tests against bcftools (optional):

```bash
# Install bcftools, then:
cargo test --workspace -- --ignored --include-ignored
```

## Before submitting a PR

```bash
cargo test --workspace          # all tests must pass
cargo clippy --workspace -- -D warnings
cargo fmt --check
```

If you changed anything in `crates/vcfkit-core/src/`, rebuild the WASM bundle:

```bash
wasm-pack build crates/vcfkit-core --target web --out-dir ../../web/public/wasm --release
git add web/public/wasm/
```

CI fails if the committed WASM artifacts are stale relative to the Rust source.

## What to work on

- Bug reports with a reproducible test case are the highest-value contributions.
- Check [known differences from bcftools](docs/known_differences.md) — those are documented gaps, not bugs.
- Open an issue before starting large changes (new commands, API changes) to avoid wasted effort.

## Correctness standard

vcfkit validates against bcftools. Any divergence on real data is a bug. Differential tests in `normalize_test.rs` and `liftover_test.rs` call bcftools directly when it's on PATH. New normalization or liftover changes must pass these.

## Code style

- No `unsafe` without a comment explaining why.
- Errors use `thiserror` in `vcfkit-core`, `anyhow` in `vcfkit-cli`.
- `--ask` stays entirely in `vcfkit-cli`, never in `vcfkit-core` (keeps WASM clean).
- Single static binary is a hard requirement — no new C dependencies, no Python.

## License

By contributing, you agree your contributions are licensed under the [MIT License](LICENSE).
