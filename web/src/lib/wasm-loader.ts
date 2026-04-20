export interface VcfkitWasm {
  filter_vcf: (vcf: string, expression: string) => string;
  normalize_vcf: (vcf: string) => string;
  liftover_vcf: (vcf: string, chainBytes: Uint8Array) => string;
}

let mod: VcfkitWasm | null = null;
let initPromise: Promise<VcfkitWasm> | null = null;

export async function ensureWasm(): Promise<VcfkitWasm> {
  if (mod) return mod;
  if (initPromise) return initPromise;
  initPromise = (async () => {
    // Use Function constructor so Vite's static analyzer never sees the
    // /public path — avoids "Cannot import non-asset file" in dev mode.
    // The build is unaffected (rollupOptions.external handles it there).
    const dynamicImport = new Function('u', 'return import(u)');
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const m: any = await dynamicImport('/wasm/vcfkit_core.js');
    // Pass explicit WASM URL so the browser knows where to fetch it.
    await m.default('/wasm/vcfkit_core_bg.wasm');
    mod = m as VcfkitWasm;
    return mod;
  })();
  return initPromise;
}
