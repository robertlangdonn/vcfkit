#!/usr/bin/env bash
# benches/e2e/run.sh — compare vcfkit vs bcftools using hyperfine
#
# Usage: bash benches/e2e/run.sh
#
# Requires: bcftools, hyperfine (both in PATH)
# If either is missing the script exits 0 with an install hint so CI is not broken.
#
# Output: benches/e2e/report.md  (markdown table of timing results)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
E2E_DIR="${SCRIPT_DIR}"
REPORT="${E2E_DIR}/report.md"
CORPUS_DIR="${REPO_ROOT}/tests/corpus/synthetic"

# ── Dependency checks ─────────────────────────────────────────────────────────

check_dep() {
    if ! command -v "$1" &>/dev/null; then
        echo "INFO: '$1' not found — skipping E2E benchmarks."
        echo "      Install hint: $2"
        exit 0
    fi
}

check_dep bcftools "brew install bcftools  OR  conda install -c bioconda bcftools"
check_dep hyperfine "brew install hyperfine  OR  cargo install hyperfine"

# ── Build vcfkit release binary ───────────────────────────────────────────────

VCFKIT_BIN="${REPO_ROOT}/target/release/vcfkit"

echo "==> Building vcfkit (release)..."
(cd "${REPO_ROOT}" && cargo build --release --quiet)

if [[ ! -x "${VCFKIT_BIN}" ]]; then
    echo "ERROR: vcfkit binary not found at ${VCFKIT_BIN}"
    exit 1
fi

# ── Ensure a reference FASTA exists with a .fai index ────────────────────────

REF_FA="${CORPUS_DIR}/mini_ref.fa"
if [[ ! -f "${REF_FA}.fai" ]]; then
    echo "==> Indexing mini_ref.fa with samtools faidx..."
    if command -v samtools &>/dev/null; then
        samtools faidx "${REF_FA}"
    else
        echo "INFO: samtools not found; skipping faidx (normalize bench may fail)"
    fi
fi

# ── Prepare a slightly larger synthetic VCF for benchmarking ─────────────────
# We replicate the basic.vcf 100× to get ~500 records for a stable measurement.

BENCH_VCF="${E2E_DIR}/bench_input.vcf"

if [[ ! -f "${BENCH_VCF}" ]]; then
    echo "==> Generating bench_input.vcf (~500 records)..."
    BASIC="${CORPUS_DIR}/basic.vcf"
    # Copy header from basic.vcf
    grep '^#' "${BASIC}" > "${BENCH_VCF}"
    # Append data lines 100 times, adjusting POS to avoid collisions
    local_i=0
    while IFS= read -r line; do
        if [[ "${line}" == "#"* ]]; then continue; fi
        for offset in $(seq 0 100 9900); do
            # Extract POS field (column 2), add offset
            IFS=$'\t' read -ra fields <<< "${line}"
            fields[1]=$(( fields[1] + offset + local_i ))
            printf '%s\n' "$(IFS=$'\t'; echo "${fields[*]}")"
        done
        (( local_i++ )) || true
    done < "${BASIC}" >> "${BENCH_VCF}"
fi

# ── Benchmark: normalize ──────────────────────────────────────────────────────

echo ""
echo "==> Benchmarking: vcfkit normalize vs bcftools norm"

NORM_RESULTS="${E2E_DIR}/normalize_bench.json"

hyperfine \
    --warmup 2 \
    --min-runs 5 \
    --export-json "${NORM_RESULTS}" \
    --command-name "vcfkit normalize" \
        "${VCFKIT_BIN} normalize --ref ${REF_FA} ${BENCH_VCF} -o /dev/null" \
    --command-name "bcftools norm" \
        "bcftools norm -f ${REF_FA} ${BENCH_VCF} -o /dev/null" \
    || true   # don't fail if bcftools errors on missing index

# ── Benchmark: filter ─────────────────────────────────────────────────────────

echo ""
echo "==> Benchmarking: vcfkit filter vs bcftools view"

# Use a VCF with INFO/AF fields for the filter bench
FILTER_VCF="${CORPUS_DIR}/multi_allelic.vcf"

FILTER_RESULTS="${E2E_DIR}/filter_bench.json"

hyperfine \
    --warmup 2 \
    --min-runs 5 \
    --export-json "${FILTER_RESULTS}" \
    --command-name "vcfkit filter" \
        "${VCFKIT_BIN} filter 'INFO/AF < 0.3' ${FILTER_VCF} -o /dev/null" \
    --command-name "bcftools view -i" \
        "bcftools view -i 'AF<0.3' ${FILTER_VCF} -o /dev/null" \
    || true

# ── Generate report.md ────────────────────────────────────────────────────────

echo ""
echo "==> Generating ${REPORT}..."

parse_mean() {
    local json_file="$1"
    local cmd_name="$2"
    if [[ ! -f "${json_file}" ]]; then
        echo "N/A"
        return
    fi
    # Extract mean from hyperfine JSON: results[].mean (seconds)
    python3 -c "
import json, sys
data = json.load(open('${json_file}'))
for r in data.get('results', []):
    if r.get('command', '') == '${cmd_name}':
        ms = r['mean'] * 1000
        stddev = r.get('stddev', 0) * 1000
        print(f'{ms:.1f} ms ± {stddev:.1f} ms')
        sys.exit(0)
print('N/A')
" 2>/dev/null || echo "N/A"
}

VCFKIT_NORM=$(parse_mean "${NORM_RESULTS}" "vcfkit normalize")
BCFTOOLS_NORM=$(parse_mean "${NORM_RESULTS}" "bcftools norm")
VCFKIT_FILT=$(parse_mean "${FILTER_RESULTS}" "vcfkit filter")
BCFTOOLS_FILT=$(parse_mean "${FILTER_RESULTS}" "bcftools view -i")

cat > "${REPORT}" <<EOF
# E2E Benchmark Results

Generated: $(date -u +"%Y-%m-%d %Human:%M:%S UTC" 2>/dev/null || date -u)

Corpus: \`tests/corpus/synthetic/\`
Target: vcfkit should be within **1.5×** of bcftools throughput.

## normalize

| Command | Mean time |
|---------|-----------|
| \`vcfkit normalize\` | ${VCFKIT_NORM} |
| \`bcftools norm\` | ${BCFTOOLS_NORM} |

## filter

| Command | Mean time |
|---------|-----------|
| \`vcfkit filter\` | ${VCFKIT_FILT} |
| \`bcftools view -i\` | ${BCFTOOLS_FILT} |

---
*Run \`bash benches/e2e/run.sh\` to regenerate.*
EOF

echo ""
echo "Report written to: ${REPORT}"
cat "${REPORT}"
