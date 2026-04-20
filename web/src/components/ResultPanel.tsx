import { VcfEditor } from './VcfEditor';
import { downloadVcf, formatDuration } from '../lib/format-utils';

interface Stats {
  inputRecords: number;
  outputRecords: number;
  durationMs: number;
}

interface Props {
  value: string;
  stats: Stats | null;
  status: 'idle' | 'running' | 'success' | 'error';
  error: string | null;
  operation: 'filter' | 'normalize' | 'liftover';
}

export function ResultPanel({ value, stats, status, error, operation }: Props) {
  const handleDownload = () => downloadVcf(value, `vcfkit-${operation}-output.vcf`);

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '8px' }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
        <span style={{ fontSize: '12px', color: 'var(--sl-color-gray-3)' }}>Output</span>
        {value && (
          <button
            onClick={handleDownload}
            style={{
              fontSize: '11px',
              padding: '3px 10px',
              background: 'var(--sl-color-gray-6)',
              border: '1px solid var(--sl-color-gray-5)',
              borderRadius: '4px',
              cursor: 'pointer',
              color: 'var(--sl-color-text)',
            }}
          >
            ↓ Download .vcf
          </button>
        )}
      </div>

      <div style={{ position: 'relative' }}>
        <VcfEditor value={value} readOnly placeholder="Output appears here after running…" height="260px" />
        {!value && status === 'idle' && (
          <div style={{
            position: 'absolute', inset: 0,
            display: 'flex', alignItems: 'center', justifyContent: 'center',
            pointerEvents: 'none',
          }}>
            <span style={{ fontSize: '12px', color: 'var(--sl-color-gray-4)', fontFamily: 'ui-monospace, monospace' }}>
              Press Run to see output
            </span>
          </div>
        )}
      </div>

      {status === 'running' && (
        <div className="result-status">
          <span style={{ opacity: 0.6 }}>Running…</span>
        </div>
      )}

      {status === 'success' && stats && (
        <div className="result-status result-ok">
          ✓{' '}
          {operation === 'filter' && (
            <>
              {stats.outputRecords.toLocaleString()} of {stats.inputRecords.toLocaleString()} records
              passed — {formatDuration(stats.durationMs)}
            </>
          )}
          {operation === 'normalize' && (
            <>
              {stats.inputRecords.toLocaleString()} in → {stats.outputRecords.toLocaleString()} out
              (multi-allelic split) — {formatDuration(stats.durationMs)}
            </>
          )}
          {operation === 'liftover' && (
            <>
              {stats.outputRecords.toLocaleString()} of {stats.inputRecords.toLocaleString()} records
              lifted — {formatDuration(stats.durationMs)}
            </>
          )}
        </div>
      )}

      {status === 'error' && error && (
        <div className="result-status result-err">
          <div>❌ {error}</div>
          {error.includes('&&') === false && error.includes('&') && (
            <div className="error-hint">Hint: use &amp;&amp; for AND, not &amp;</div>
          )}
        </div>
      )}
    </div>
  );
}
