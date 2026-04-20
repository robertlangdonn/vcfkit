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
      <span style={{ fontSize: '12px', color: 'var(--sl-color-gray-3)', flexShrink: 0 }}>
        Example:
      </span>
      <div style={{ display: 'flex', gap: '6px', flexWrap: 'wrap' }}>
        {examples.map((ex) => {
          const isActive = ex.id === selectedId;
          return (
            <button
              key={ex.id}
              onClick={() => onSelect(ex)}
              style={{
                fontSize: '11px',
                padding: '3px 10px',
                background: isActive ? 'color-mix(in srgb, var(--sl-color-accent) 18%, transparent)' : 'var(--sl-color-gray-7)',
                border: `1px solid ${isActive ? 'var(--sl-color-accent)' : 'var(--sl-color-gray-5)'}`,
                borderRadius: '4px',
                cursor: 'pointer',
                color: isActive ? 'var(--sl-color-accent-high, var(--sl-color-accent))' : 'var(--sl-color-text)',
                whiteSpace: 'nowrap',
                fontWeight: isActive ? '600' : '400',
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
