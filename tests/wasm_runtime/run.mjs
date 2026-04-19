/**
 * WASM parity tests: verify that the wasm-bindgen functions produce output
 * byte-identical to the native CLI for the same inputs.
 *
 * Run: node tests/wasm_runtime/run.mjs
 * (from repo root, after building the WASM package)
 */

import { readFile } from 'fs/promises';
import { fileURLToPath } from 'url';
import { dirname, join, resolve } from 'path';

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(__dirname, '../..');

// The wasm-pack output is in web/public/wasm/
const wasmPkgDir = join(repoRoot, 'web/public/wasm');

async function loadFile(relPath) {
  return readFile(join(__dirname, relPath), 'utf-8');
}

function normalizeVcf(text) {
  // Strip trailing whitespace per line for comparison stability
  return text
    .split('\n')
    .map((l) => l.trimEnd())
    .filter((l, i, arr) => !(i === arr.length - 1 && l === ''))
    .join('\n');
}

async function main() {
  // Import the WASM module from its built location
  const { default: init, filter_vcf, normalize_vcf, liftover_vcf } = await import(
    join(wasmPkgDir, 'vcfkit_core.js')
  );

  // Node requires a BufferSource or WebAssembly.Module — not a bare file path.
  const wasmBytes = await readFile(join(wasmPkgDir, 'vcfkit_core_bg.wasm'));
  await init({ module_or_path: wasmBytes.buffer });

  const encoder = new TextEncoder();

  const tests = [
    // ── filter ────────────────────────────────────────────────────────────────
    {
      name: 'filter: QUAL > 50',
      run: async () => {
        const input = await loadFile('fixtures/simple.vcf');
        const expected = await loadFile('expected/filter_qual_gt_50.vcf');
        const result = filter_vcf(input, 'QUAL > 50');
        return { result, expected };
      },
    },
    {
      name: 'filter: CHROM == chr1',
      run: async () => {
        const input = await loadFile('fixtures/simple.vcf');
        const expected = await loadFile('expected/filter_chr1.vcf');
        const result = filter_vcf(input, "CHROM == 'chr1'");
        return { result, expected };
      },
    },
    {
      name: 'filter: INFO/AF < 0.25 (any-element semantics on multi-allelic)',
      run: async () => {
        const input = await loadFile('fixtures/multi_allelic.vcf');
        const expected = await loadFile('expected/filter_af_lt_0.25.vcf');
        const result = filter_vcf(input, 'INFO/AF < 0.25');
        return { result, expected };
      },
    },
    {
      name: 'filter: compound QUAL > 50 && FILTER == PASS',
      run: async () => {
        const input = await loadFile('fixtures/simple.vcf');
        const expected = await loadFile('expected/filter_qual_and_pass.vcf');
        const result = filter_vcf(input, "QUAL > 50 && FILTER == 'PASS'");
        return { result, expected };
      },
    },
    {
      name: 'filter: FILTER == PASS (mixed filter values)',
      run: async () => {
        const input = await loadFile('fixtures/filter_pass_mixed.vcf');
        const expected = await loadFile('expected/filter_pass_only.vcf');
        const result = filter_vcf(input, "FILTER == 'PASS'");
        return { result, expected };
      },
    },

    {
      name: 'filter: OR expression (INFO/AF > 0.35 || QUAL > 60)',
      run: async () => {
        const input = await loadFile('fixtures/multi_allelic.vcf');
        const expected = await loadFile('expected/filter_or.vcf');
        const result = filter_vcf(input, 'INFO/AF > 0.35 || QUAL > 60');
        return { result, expected };
      },
    },
    {
      name: 'filter: negation !(FILTER == PASS)',
      run: async () => {
        const input = await loadFile('fixtures/filter_pass_mixed.vcf');
        const expected = await loadFile('expected/filter_negation.vcf');
        const result = filter_vcf(input, "!(FILTER == 'PASS')");
        return { result, expected };
      },
    },
    {
      name: 'filter: regex INFO/GENE ~ BRCA',
      run: async () => {
        const input = await loadFile('fixtures/filter_anno.vcf');
        const expected = await loadFile('expected/filter_gene_regex.vcf');
        const result = filter_vcf(input, "INFO/GENE ~ 'BRCA'");
        return { result, expected };
      },
    },
    {
      name: 'filter: POS range (POS >= 200 && POS <= 300)',
      run: async () => {
        const input = await loadFile('fixtures/simple.vcf');
        const expected = await loadFile('expected/filter_pos_range.vcf');
        const result = filter_vcf(input, 'POS >= 200 && POS <= 300');
        return { result, expected };
      },
    },

    // ── normalize ─────────────────────────────────────────────────────────────
    {
      name: 'normalize: split multi-allelic',
      run: async () => {
        const input = await loadFile('fixtures/multi_allelic.vcf');
        const expected = await loadFile('expected/normalize_split.vcf');
        const result = normalize_vcf(input);
        return { result, expected };
      },
    },
    {
      name: 'normalize: SNV-only passthrough (no multi-allelic)',
      run: async () => {
        const input = await loadFile('fixtures/normalize_snv_only.vcf');
        const expected = await loadFile('expected/normalize_snv_passthrough.vcf');
        const result = normalize_vcf(input);
        return { result, expected };
      },
    },
    {
      name: 'normalize: Number=G INFO field split across alleles',
      run: async () => {
        const input = await loadFile('fixtures/normalize_number_g.vcf');
        const expected = await loadFile('expected/normalize_number_g_split.vcf');
        const result = normalize_vcf(input);
        return { result, expected };
      },
    },

    // ── liftover ──────────────────────────────────────────────────────────────
    {
      name: 'liftover: identity chain',
      run: async () => {
        // A synthetic identity chain: maps chr1:0-1000000 → chr1:0-1000000 unchanged.
        const chain = 'chain 1000000 chr1 248956422 + 0 1000000 chr1 248956422 + 0 1000000 1\n1000000\n\n';
        const input = await loadFile('fixtures/simple.vcf');
        const result = liftover_vcf(input, encoder.encode(chain));
        return { result, expected: null, check: () => {
          if (!result.includes('chr1\t100')) throw new Error('chr1:100 missing from liftover output');
          if (!result.includes('chr1\t200')) throw new Error('chr1:200 missing from liftover output');
          if (result.includes('chr2\t')) throw new Error('chr2 records should be unmapped and dropped');
        }};
      },
    },
    {
      name: 'liftover: partial chain (only chr1:0-250 covered, rest unmapped)',
      run: async () => {
        // Chain covers only chr1 positions 0–250; chr1:300+ and all chr2 are unmapped.
        const chain = 'chain 250 chr1 248956422 + 0 250 chr1 248956422 + 0 250 1\n250\n\n';
        const input = await loadFile('fixtures/simple.vcf');
        const result = liftover_vcf(input, encoder.encode(chain));
        return { result, expected: null, check: () => {
          if (!result.includes('chr1\t100')) throw new Error('chr1:100 should be mapped');
          if (!result.includes('chr1\t200')) throw new Error('chr1:200 should be mapped');
          if (result.includes('chr1\t300')) throw new Error('chr1:300 is outside chain and should be unmapped');
          if (result.includes('chr2\t')) throw new Error('chr2 records should be unmapped');
        }};
      },
    },
    {
      name: 'liftover: offset chain (chr1 positions shifted +1000)',
      run: async () => {
        // Chain maps chr1:0-1000000 → chr1:1000-1001000 (shift all positions by +1000).
        const chain = 'chain 1000000 chr1 248956422 + 0 1000000 chr1 248956422 + 1000 1001000 1\n1000000\n\n';
        const input = await loadFile('fixtures/simple.vcf');
        const result = liftover_vcf(input, encoder.encode(chain));
        return { result, expected: null, check: () => {
          if (!result.includes('chr1\t1100')) throw new Error('chr1:100 should lift to chr1:1100');
          if (!result.includes('chr1\t1200')) throw new Error('chr1:200 should lift to chr1:1200');
          if (!result.includes('chr1\t1300')) throw new Error('chr1:300 should lift to chr1:1300');
          if (result.includes('chr2\t')) throw new Error('chr2 records not in chain should be unmapped');
        }};
      },
    },
  ];

  let passed = 0;
  let failed = 0;

  for (const { name, run } of tests) {
    try {
      const { result, expected, check } = await run();
      if (check) {
        check();
        console.log(`PASS: ${name}`);
        passed++;
      } else if (normalizeVcf(result) === normalizeVcf(expected)) {
        console.log(`PASS: ${name}`);
        passed++;
      } else {
        console.error(`FAIL: ${name}`);
        const resultLines = normalizeVcf(result).split('\n');
        const expectedLines = normalizeVcf(expected).split('\n');
        for (let i = 0; i < Math.max(resultLines.length, expectedLines.length); i++) {
          if (resultLines[i] !== expectedLines[i]) {
            console.error(`  Line ${i + 1}:`);
            console.error(`    expected: ${JSON.stringify(expectedLines[i])}`);
            console.error(`    got:      ${JSON.stringify(resultLines[i])}`);
            if (i > 5) { console.error('  ... (truncated)'); break; }
          }
        }
        failed++;
      }
    } catch (e) {
      console.error(`ERROR: ${name}: ${e.message}`);
      failed++;
    }
  }

  console.log(`\n${passed}/${passed + failed} tests passed`);
  process.exit(failed > 0 ? 1 : 0);
}

main().catch((e) => {
  console.error('Fatal:', e);
  process.exit(1);
});
