import { EXAMPLES, type VcfExample } from '../lib/vcf-examples';

interface Props {
  operation: 'filter' | 'normalize' | 'liftover';
  onSelect: (example: VcfExample) => void;
  selectedId?: string | null;
}

export function ExamplePicker({ operation, onSelect, selectedId }: Props) {
  const examples = EXAMPLES[operation];

  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: '8px', flexWrap: 'wrap' }}>
      <span style={{ display: 'flex', alignItems: 'center', gap: '4px', fontSize: '11px', color: 'var(--sl-color-gray-3)', flexShrink: 0, fontWeight: 600, textTransform: 'uppercase', letterSpacing: '0.06em' }}>
        <svg width="11" height="11" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" aria-hidden="true"><polygon points="12 2 15.09 8.26 22 9.27 17 14.14 18.18 21.02 12 17.77 5.82 21.02 7 14.14 2 9.27 8.91 8.26 12 2"/></svg>
        Examples
      </span>
      <div style={{ display: 'flex', gap: '5px', flexWrap: 'wrap' }}>
        {examples.map((ex) => {
          const isActive = ex.id === selectedId;
          return (
            <button
              key={ex.id}
              onClick={() => onSelect(ex)}
              style={{
                fontSize: '11px',
                padding: '4px 12px',
                background: isActive ? 'color-mix(in srgb, var(--sl-color-accent) 18%, transparent)' : 'transparent',
                border: `1px solid ${isActive ? 'var(--sl-color-accent)' : 'var(--sl-color-gray-5)'}`,
                borderRadius: '20px',
                cursor: 'pointer',
                color: isActive ? 'var(--sl-color-accent-high, var(--sl-color-accent))' : 'var(--sl-color-gray-3)',
                whiteSpace: 'nowrap',
                fontWeight: isActive ? '600' : '400',
                transition: 'color 0.15s, border-color 0.15s, background 0.15s',
              }}
            >
              {ex.label}
            </button>
          );
        })}
      </div>
    </div>
  );
}
