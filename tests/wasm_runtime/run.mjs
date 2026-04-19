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

  const tests = [
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
      name: 'normalize: split multi-allelic',
      run: async () => {
        const input = await loadFile('fixtures/multi_allelic.vcf');
        const expected = await loadFile('expected/normalize_split.vcf');
        const result = normalize_vcf(input);
        return { result, expected };
      },
    },
    {
      name: 'liftover: identity chain',
      run: async () => {
        // A synthetic identity chain: maps chr1:0-1000000 → chr1:0-1000000 unchanged.
        const chain = 'chain 1000000 chr1 248956422 + 0 1000000 chr1 248956422 + 0 1000000 1\n1000000\n\n';
        const encoder = new TextEncoder();
        const input = await loadFile('fixtures/simple.vcf');
        const result = liftover_vcf(input, encoder.encode(chain));
        // All chr1 records should appear in output; chr2 records unmapped (dropped).
        return { result, expected: null, check: () => {
          if (!result.includes('chr1\t100')) throw new Error('chr1:100 missing from liftover output');
          if (!result.includes('chr1\t200')) throw new Error('chr1:200 missing from liftover output');
          if (result.includes('chr2\t')) throw new Error('chr2 records should be unmapped and dropped');
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
