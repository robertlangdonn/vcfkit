# vcfkit Performance Audit

Audited 2026-04-18. Focus: per-record heap allocations in hot loops.

## Already fixed

- **normalize: FASTA cache miss** — `absent_contigs: HashSet<String>` prevents repeated `build_from_path()` calls. Was 87× excess syscalls on b37/UCSC mismatch. Fixed in v0.1.2.
- **filter: fast path** — `FastRecord` bypasses noodles entirely. Uses `BufRead::read_line`, zero-copy `&str` slices, raw byte writes. 4× faster than bcftools on 1000G chr22.

## Hot loop allocations — live issues

### normalize.rs

| Line | Pattern | Every record? | Notes |
|------|---------|--------------|-------|
| 265-269 | `to_vec()` for REF + ALT | Yes | Two Vec allocs per record before left-align check |
| 276, 328 | `chrom.to_string()` | Yes | String per record for contig lookup |
| 283-284 | `r.clone()` + `a.clone()` | Yes | Clone of already-allocated Vecs |
| 388 | `record.alternate_bases().to_vec()` | Per multi-allelic | Vec copy before split |
| 404-405 | `keys.clone()`, `old_values.clone()` | Per ALT during split | Clones entire samples struct |
| 406-410 | `Vec::with_capacity(...)` (2×) | Per ALT, per sample | Nested Vec rebuilds |
| 625, 636 | `name.to_string()` | Per new contig | String for `absent_contigs` + cache key |

**Biggest win available:** Hoist the REF/ALT Vec allocation outside the left-align check; only allocate if the record is actually an indel candidate. For SNPs (majority of typical VCFs) this would eliminate 2 Vec allocs per record.

### liftover.rs

| Line | Pattern | Every record? | Notes |
|------|---------|--------------|-------|
| 538 | `chrom.to_string()` | Yes | String for chain lookup |
| 556 | `ref_bases.to_string()` | Yes | String copy of REF |
| 557 | `alts.to_vec()` + `alt.clone()` | Yes | Vec + String per ALT |
| 579-590 | `reverse_complement()` returning `String` | Per minus-strand record | Unavoidable without different data model |
| 593 | `ref_bases.clone()`, `alts.clone()` | Per forward-strand record | Clones right after alloc |
| 597 | `record.clone()` | Yes | Deep clone of full noodles Record |
| 612-618 | `String::from("SRC_CONTIG")`, `String::from("SRC_POS")` | Per record (if option set) | Two static string keys allocated per record |

**Biggest win available:** `record.clone()` at line 597 is the most expensive single allocation — copies all INFO/FORMAT/sample data. Ideally we'd construct a new minimal Record instead. Second biggest: `alts.to_vec()` creates a new Vec + clones each alt String even for single-allelic records.

### filter.rs — fast path

The fast path is already good. Remaining allocations:

| Line | Pattern | Every record? | Notes |
|------|---------|--------------|-------|
| 403 | `rec.chrom.to_string()` | Per CHROM access | Could compare `&str` directly in many cases |
| 508, 513 | `raw.split(',').collect::<Vec<_>>()` | Per multi-value INFO | Vec per multi-allelic INFO lookup |
| 559 | `joined.join(",")` | Per string INFO array | String for regex matching |
| 896-901 | `.map(|v| v.to_string()).collect::<Vec<_>>()` | Per FloatArray ~ compare | Vec + N Strings for regex on float arrays |

**Note:** The fast path only parses fields the expression actually references. If your expression is `INFO/AF < 0.01`, CHROM/FILTER/FORMAT are never touched.

### filter.rs — noodles fallback

The fallback path is used when noodles was already invoked (for headers, BCF, etc.). Its allocations are numerous but it's not the hot path:

- `info_array_all`: 3 separate `Vec<f64>` collection paths
- `format_genotype`: `String::new()` rebuilt per genotype
- `compare` string fallback: `scalar_to_string()` on both sides

**Not prioritizing** — the fast path is used for all plain VCF; the fallback only runs when noodles has already paid its overhead.

## Prioritized next actions

1. **normalize fast path** (Step 3 in roadmap) — Same approach as filter. Read raw line, only parse fields the operation needs. Expected: match/beat bcftools on chr22. The REF/ALT Vec alloc issue above is moot once we do this; we'd avoid noodles entirely for the hot path.

2. **liftover `record.clone()`** — Consider building the lifted record from scratch using only the fields that change (CHROM, POS, REF, ALT) rather than deep-cloning and mutating. Depends on noodles API surface.

3. **liftover static INFO keys** — `String::from("SRC_CONTIG")` / `String::from("SRC_POS")` should use `static` or `Lazy<String>` if the `write_src_coords` option is common.

4. **normalize split path** — `old_values.clone()` cloning the entire samples struct per ALT is expensive for multi-sample VCFs with many FORMAT fields. Worth profiling with a real multi-sample VCF.

## Methodology

- Grep for `String::new`, `to_string()`, `to_vec()`, `.clone()`, `Vec::new`, `Vec::with_capacity`, `File::open`, `build_from_path` in hot loop context
- Cross-referenced with `cargo bench` baselines and hyperfine E2E numbers
- Confirmed fast path is already in use for `filter`; normalize and liftover still use full noodles per-record path
