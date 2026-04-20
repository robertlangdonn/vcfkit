import { useState, useCallback, useRef } from 'react';
import { VcfEditor } from './VcfEditor';
import { ResultPanel } from './ResultPanel';
import { CliEquivalent } from './CliEquivalent';
import { ExamplePicker } from './ExamplePicker';
import { EXAMPLES, type VcfExample } from '../lib/vcf-examples';
import { ensureWasm } from '../lib/wasm-loader';
import { countRecords, truncateForDemo, MAX_DEMO_RECORDS } from '../lib/format-utils';

type OperationId = 'filter' | 'normalize' | 'liftover';

interface RunStats {
  inputRecords: number;
  outputRecords: number;
  durationMs: number;
}

interface DemoState {
  input: string;
  output: string;
  expression: string;
  chainFile: string;
  status: 'idle' | 'running' | 'success' | 'error';
  error: string | null;
  stats: RunStats | null;
  truncated: boolean;
  selectedExampleId: string | null;
}

const DEFAULT_STATE: DemoState = {
  input: EXAMPLES.filter[0].vcf,
  output: '',
  expression: EXAMPLES.filter[0].expression ?? '',
  chainFile: '',
  status: 'idle',
  error: null,
  stats: null,
  truncated: false,
  selectedExampleId: EXAMPLES.filter[0].id,
};

const TAB_LABELS: Record<OperationId, string> = {
  filter: 'filter',
  normalize: 'normalize',
  liftover: 'liftover',
};

export function Demo() {
  const [activeTab, setActiveTab] = useState<OperationId>('filter');
  const [states, setStates] = useState<Record<OperationId, DemoState>>({
    filter: { ...DEFAULT_STATE },
    normalize: { ...DEFAULT_STATE, input: EXAMPLES.normalize[0].vcf, expression: '', selectedExampleId: EXAMPLES.normalize[0].id },
    liftover: { ...DEFAULT_STATE, input: EXAMPLES.liftover[0].vcf, expression: '', chainFile: '', selectedExampleId: EXAMPLES.liftover[0].id },
  });

  const fileInputRef = useRef<HTMLInputElement>(null);

  const state = states[activeTab];

  const update = useCallback(
    (tab: OperationId, patch: Partial<DemoState>) => {
      setStates((prev) => ({ ...prev, [tab]: { ...prev[tab], ...patch } }));
    },
    [],
  );

  const handleExampleSelect = useCallback(
    (ex: VcfExample) => {
      update(activeTab, {
        input: ex.vcf,
        expression: ex.expression ?? state.expression,
        chainFile: ex.chain ?? state.chainFile,
        output: '',
        status: 'idle',
        error: null,
        stats: null,
        selectedExampleId: ex.id,
      });
    },
    [activeTab, state.expression, state.chainFile, update],
  );

  const handleFileUpload = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const file = e.target.files?.[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = () => {
        const text = reader.result as string;
        const { text: truncated, truncated: wasTruncated } = truncateForDemo(text);
        update(activeTab, {
          input: truncated,
          truncated: wasTruncated,
          output: '',
          status: 'idle',
          error: null,
          stats: null,
        });
      };
      reader.readAsText(file);
      // Reset so re-uploading the same file triggers onChange
      e.target.value = '';
    },
    [activeTab, update],
  );

  const handleRun = useCallback(async () => {
    const s = states[activeTab];
    update(activeTab, { status: 'running', error: null, output: '', stats: null });

    try {
      const wasm = await ensureWasm();
      const t0 = performance.now();
      let output: string;
      const inputRecords = countRecords(s.input);

      switch (activeTab) {
        case 'filter': {
          if (!s.expression.trim()) {
            update(activeTab, { status: 'error', error: 'Expression is required.' });
            return;
          }
          output = wasm.filter_vcf(s.input, s.expression);
          break;
        }
        case 'normalize': {
          output = wasm.normalize_vcf(s.input);
          break;
        }
        case 'liftover': {
          if (!s.chainFile.trim()) {
            update(activeTab, {
              status: 'error',
              error: 'Paste a UCSC chain file before running.',
            });
            return;
          }
          const encoder = new TextEncoder();
          output = wasm.liftover_vcf(s.input, encoder.encode(s.chainFile));
          break;
        }
      }

      const durationMs = performance.now() - t0;
      const outputRecords = countRecords(output);

      update(activeTab, {
        status: 'success',
        output,
        stats: { inputRecords, outputRecords, durationMs },
      });
    } catch (err) {
      update(activeTab, {
        status: 'error',
        error: err instanceof Error ? err.message : String(err),
      });
    }
  }, [activeTab, states, update]);

  return (
    <div className="demo-root">
      {/* Privacy note */}
      <div className="privacy-note">
        <svg width="11" height="11" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" aria-hidden="true" style={{ flexShrink: 0 }}>
          <rect x="3" y="11" width="18" height="11" rx="2" ry="2"/>
          <path d="M7 11V7a5 5 0 0 1 10 0v4"/>
        </svg>
        Your VCF never leaves this tab — all processing runs in WebAssembly locally.
      </div>

      {/* Tabs */}
      <div className="demo-tabs" role="tablist">
        {(Object.keys(TAB_LABELS) as OperationId[]).map((tab) => (
          <button
            key={tab}
            role="tab"
            aria-selected={activeTab === tab}
            className={`demo-tab ${activeTab === tab ? 'active' : ''}`}
            onClick={() => setActiveTab(tab)}
          >
            {TAB_LABELS[tab]}
          </button>
        ))}
      </div>

      {/* Toolbar: example picker + file upload */}
      <div className="demo-toolbar">
        <ExamplePicker operation={activeTab} onSelect={handleExampleSelect} selectedId={state.selectedExampleId} />
        <div style={{ display: 'flex', gap: '6px', alignItems: 'center' }}>
          <button
            className="toolbar-btn"
            onClick={() => fileInputRef.current?.click()}
            title="Upload a .vcf file"
          >
            Upload .vcf
          </button>
          <input
            ref={fileInputRef}
            type="file"
            accept=".vcf,.vcf.gz"
            style={{ display: 'none' }}
            onChange={handleFileUpload}
          />
        </div>
      </div>

      {state.truncated && (
        <div className="truncation-notice">
          File truncated to {MAX_DEMO_RECORDS.toLocaleString()} records for browser performance.
          Use the CLI for large files:{' '}
          <code>cargo install vcfkit-cli</code>
        </div>
      )}

      {/* Main editing area */}
      <div className="demo-editors">
        <div className="demo-editor-col">
          <span className="editor-label">Input VCF</span>
          <VcfEditor
            value={state.input}
            onChange={(v) => update(activeTab, { input: v })}
          />
        </div>
        <div className="demo-editor-col">
          <ResultPanel
            value={state.output}
            stats={state.stats}
            status={state.status}
            error={state.error}
            operation={activeTab}
          />
        </div>
      </div>

      {/* Operation-specific controls */}
      <div className="demo-controls">
        {activeTab === 'filter' && (
          <div style={{ flex: 1 }}>
            <label className="control-label" htmlFor="filter-expr">
              Expression
            </label>
            <input
              id="filter-expr"
              type="text"
              className="control-input"
              value={state.expression}
              onChange={(e) => update('filter', { expression: e.target.value })}
              placeholder="INFO/AF < 0.01 && FILTER == 'PASS'"
              onKeyDown={(e) => e.key === 'Enter' && handleRun()}
            />
          </div>
        )}

        {activeTab === 'normalize' && (
          <div style={{ flex: 1 }}>
            <p className="control-note">
              Multi-allelic sites will be split into biallelic records. Left-alignment requires a
              reference FASTA and is not available in the browser — use{' '}
              <code>vcfkit normalize -f ref.fa</code> in the CLI.
            </p>
          </div>
        )}

        {activeTab === 'liftover' && (
          <div style={{ flex: 1 }}>
            <label className="control-label" htmlFor="chain-input">
              Chain file content{' '}
              <a
                href="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/"
                target="_blank"
                rel="noopener noreferrer"
                style={{ fontSize: '11px', opacity: 0.7 }}
              >
                Download from UCSC →
              </a>
            </label>
            <textarea
              id="chain-input"
              className="control-textarea"
              value={state.chainFile}
              onChange={(e) => update('liftover', { chainFile: e.target.value })}
              placeholder="Paste the contents of a UCSC .chain file here…"
              rows={4}
            />
          </div>
        )}

        <button
          className="run-btn"
          onClick={handleRun}
          disabled={state.status === 'running'}
        >
          {state.status === 'running' ? 'Running…' : `Run ${activeTab}`}
        </button>
      </div>

      {/* CLI equivalent bar */}
      <CliEquivalent operation={activeTab} expression={state.expression} />
    </div>
  );
}

export default Demo;
