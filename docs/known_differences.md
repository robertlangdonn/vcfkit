# Known Differences Between vcfkit and bcftools

This document tracks intentional behavioral differences between `vcfkit normalize`
and `bcftools norm`. Each entry notes the version where the limitation was introduced,
a concrete example, and the planned fix.

---

## 1. Multi-allelic indel left-alignment (v0.1.x — planned fix in v0.2)

### Behavior

`vcfkit normalize` passes multi-allelic indel records through **unchanged** when
`--no-split` is used (or when splitting is disabled). It does **not** attempt to
left-align multi-allelic indels.

`bcftools norm` (without `-m`) will left-align multi-allelic indels jointly — i.e.
it applies the Tan 2015 algorithm treating all ALTs together as a group, choosing
the leftmost position where all alleles remain valid.

### Concrete Example

Reference (chr22, hg19 b37) around position 16404838 is a poly-A run:

```
...GGGAAAAAAAAAAAAAAAA...
         ^16404838
```

Input VCF record that is **not** fully left-aligned (shifted right by 2):

```
22  16404840  .  AA  AAA,A  100  PASS  DP=100
```

**bcftools norm output** — shifts the record left by 2 into the poly-A run:

```
22  16404838  .  GA  GAA,G  100  PASS  DP=100
```

**vcfkit normalize --no-split output** — record passes through unchanged:

```
22  16404840  .  AA  AAA,A  100  PASS  DP=100
```

Note: the 1000 Genomes Project chr22 data contains this exact record at position
16404838 (`GA -> GAA,G`), which is already fully left-aligned. The divergence
only manifests when the input record is not yet at the leftmost valid position.

### Root Cause

The biallelic Tan 2015 algorithm in `vcfkit-core/src/normalize.rs` operates on a
single (REF, ALT) pair. Applying it to only the first ALT of a multi-allelic record
would silently drop the remaining ALTs — a data-corrupting bug fixed in v0.1.4 by
skipping left-alignment for multi-allelic records entirely.

Proper joint multi-allelic left-alignment requires selecting a single anchor
position where all ALT alleles can be simultaneously represented, which is a
non-trivial extension of the algorithm.

### Planned Fix

v0.2 will implement joint Tan 2015 left-alignment across all ALTs: find the
leftmost position P such that `trim_and_extend(REF, ALT_k, P)` is valid for every
k simultaneously, then rewrite the record at P.

---

## 2. No other known differences

The differential tests in `tests/normalize_test.rs` compare vcfkit output against
`bcftools norm` on the synthetic corpus (basic SNPs, multi-allelic SNPs, unnormalized
indels) and all pass. No additional behavioral differences have been identified.
