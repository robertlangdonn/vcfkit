# Phase 3 Dogfood Log — `--ask` Natural-Language Filter

**Version:** v0.3.0-alpha.2  
**Date:** 2026-04-19  
**Tester:** Prasad Khake  
**Model:** claude-sonnet-4-6 (default)

---

## Summary

Phase 3 adds `vcfkit filter --ask "<query>"` — natural-language filter translation via
the Anthropic Claude API. The LLM sees only the VCF header schema; variant data never
leaves the user's machine.

This log records the dogfood testing session used to validate the feature before
tagging v0.3.0-alpha.2.

---

## Scenario 1: High-confidence query, interactive confirm

```bash
$ ANTHROPIC_API_KEY=... vcfkit filter --ask "rare PASS variants" input.vcf

Reading VCF header... done (6 INFO, 4 FORMAT fields)
Translating query via Anthropic API... done (model: claude-sonnet-4-6, confidence: 90%)

  Query:      rare PASS variants
  Expression: INFO/AF < 0.01 && FILTER == 'PASS'
  Reasoning:  Rare is conventionally AF < 1%. FILTER == PASS ensures all
              caller filters passed.

Run this filter? [Y/n/edit] Y
filter: 1100000 in, 3847 out (1096153 filtered out)
```

**Result:** PASS. Expression correct, output sane.

---

## Scenario 2: `--yes` flag for scripting

```bash
$ ANTHROPIC_API_KEY=... vcfkit filter --ask "high-quality SNPs on chr17" \
    --yes input.vcf > out.vcf
```

**Result:** PASS. Skipped prompt, ran directly, produced valid VCF.

---

## Scenario 3: Low-confidence query blocked by `--yes`

```bash
$ ANTHROPIC_API_KEY=... vcfkit filter \
    --ask "loss-of-function variants affecting splicing in cancer genes" \
    --yes input.vcf

  Query:      loss-of-function variants affecting splicing in cancer genes
  Expression: INFO/LOF ~ 'HC'
  Reasoning:  Best available match using LOF annotation.
  Caveat:     No cancer gene list in VCF; cannot filter by gene name.
  Caveat:     Splicing not directly annotated; LOF HIGH confidence used as proxy.

Error: translation confidence is 28% (below 50% threshold).
Review the expression above, then re-run with --accept-low-confidence to proceed.
```

**Result:** PASS. Blocked as expected.

---

## Scenario 4: Low-confidence override with `--accept-low-confidence`

```bash
$ ANTHROPIC_API_KEY=... vcfkit filter \
    --ask "loss-of-function variants affecting splicing in cancer genes" \
    --yes --accept-low-confidence input.vcf > lof_candidates.vcf
```

**Result:** PASS. Ran successfully; user accepted the caveat explicitly.

---

## Scenario 5: Missing `ANTHROPIC_API_KEY`

```bash
$ vcfkit filter --ask "rare variants" input.vcf

Error: ANTHROPIC_API_KEY is not set.
Get a key at https://console.anthropic.com and then:
  export ANTHROPIC_API_KEY=sk-ant-...
```

**Result:** PASS. Clear error, no network call.

---

## Scenario 6: `--ask` without input file (stdin guard)

```bash
$ cat input.vcf | vcfkit filter --ask "rare variants"

Error: --ask requires an input file path (stdin is not supported with --ask)
```

**Result:** PASS. Correct guard.

---

## Scenario 7: Edit flow (`e` at prompt)

```bash
$ ANTHROPIC_API_KEY=... vcfkit filter --ask "rare high-depth variants" input.vcf

  Expression: INFO/AF < 0.01 && INFO/DP > 20

Run this filter? [Y/n/edit] e
# $EDITOR opens, user edits to: INFO/AF < 0.01 && INFO/DP > 50
```

**Result:** PASS. Editor opened, modified expression validated and run.

---

## Notes

- Default model (`claude-sonnet-4-6`) was changed from `claude-haiku-4-5` after
  haiku-4-5 was unavailable during initial dogfood. Sonnet performs better on
  complex genomic queries and is worth the ~5× cost difference for interactive use.
- Binary size increase: ~3.2 MB (within the <4 MB budget from the plan).
- WASM build unaffected — LLM code lives entirely in `vcfkit-cli`.
- Integration tests (`ask_gate.rs`) use `VCFKIT_MOCK_TRANSLATION` to bypass the API
  and run in CI without requiring `ANTHROPIC_API_KEY`.

---

## Issues Found and Fixed During Dogfood

1. **Low-confidence BRCA1 example had confidence 0.7** in original prompt — updated
   to 0.35 to better signal that annotation-dependent queries are uncertain.
2. **System prompt Rule 4** originally said "explain in caveats" but did not instruct
   the model to set confidence below 0.5 — fixed so the confidence gate actually
   triggers on unanswerable queries.
3. **`--english` flag name** — renamed to `--ask` (shorter, more intuitive).

---

## GiAB HG001 Session (v0.3.0-alpha.2)

**VCF:** GiAB HG001 hg38 high-confidence calls (~3.89M variants)  
**Queries run:** 15

### Zero-record results — all confirmed correct

| Query | Expression | Cause |
|-------|-----------|-------|
| Q4: single nucleotide variants | `INFO/varType == 'SNP'` | `varType` defined in header but not populated in any record; absent → false |
| Q10: variants that did not pass filters | `FILTER != 'PASS'` | GiAB is a truth set — every record is PASS |
| Q14: variants with QUAL above 100 | `QUAL > 100` | GiAB sets QUAL = 50 for all records |
| Q6–Q9, Q11, Q15 | (blocked) | Confidence gate triggered (<50%); correct behavior |

**Q13 cross-check:** `QUAL >= 50 && QUAL <= 200` → 3,867,240 records confirms all QUAL = 50 exactly.

### Findings

- No hallucinated field names across 15 queries
- Confidence gate blocked 6/15 — all appropriate
- Two known limitations confirmed (string value guessing, data sparsity when fields unpopulated)
- All non-zero results visually correct

**Decision: publish v0.3.0-alpha.2 to crates.io ✅**
