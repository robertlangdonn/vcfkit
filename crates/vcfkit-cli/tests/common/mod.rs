//! Shared test helpers for vcfkit-cli integration tests.

// Not every integration test uses every helper here; silence the dead-code
// warnings that would otherwise fire from the per-test-binary analysis.
#![allow(dead_code)]

pub mod diff;
pub mod download;
