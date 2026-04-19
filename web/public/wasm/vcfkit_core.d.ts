/* tslint:disable */
/* eslint-disable */

/**
 * Filter a VCF (as a string) by an expression. Returns filtered VCF including header.
 *
 * # Example (JavaScript)
 * ```js
 * import init, { filter_vcf } from './vcfkit_core.js';
 * await init();
 * const out = filter_vcf(vcfText, "INFO/AF < 0.01 && FILTER == 'PASS'");
 * ```
 */
export function filter_vcf(vcf_input: string, expression: string): string;

export function init(): void;

/**
 * Liftover a VCF using a chain file provided as bytes. REF validation is
 * skipped (no FASTA). Unmapped records are silently dropped. Returns lifted VCF.
 */
export function liftover_vcf(vcf_input: string, chain_bytes: Uint8Array): string;

/**
 * Split multi-allelic VCF records. Left-alignment is skipped (requires a
 * reference FASTA, which is not available in WASM). Returns normalized VCF.
 */
export function normalize_vcf(vcf_input: string): string;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly filter_vcf: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly liftover_vcf: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly normalize_vcf: (a: number, b: number, c: number) => void;
    readonly init: () => void;
    readonly __wbindgen_export: (a: number, b: number, c: number) => void;
    readonly __wbindgen_export2: (a: number, b: number) => number;
    readonly __wbindgen_export3: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_add_to_stack_pointer: (a: number) => number;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
