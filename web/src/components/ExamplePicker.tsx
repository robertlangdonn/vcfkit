import { EXAMPLES, type VcfExample } from '../lib/vcf-examples';

interface Props {
  operation: 'filter' | 'normalize' | 'liftover';
  onSelect: (example: VcfExample) => void;
}

export function ExamplePicker({ operation, onSelect }: Props) {
  const examples = EXAMPLES[operation];

  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: '8px', flexWrap: 'wrap' }}>
      <span style={{ fontSize: '12px', color: 'var(--sl-color-gray-3)', flexShrink: 0 }}>
        Load example:
      </span>
      <div style={{ display: 'flex', gap: '6px', flexWrap: 'wrap' }}>
        {examples.map((ex) => (
          <button
            key={ex.id}
            onClick={() => onSelect(ex)}
            style={{
              fontSize: '11px',
              padding: '3px 10px',
              background: 'var(--sl-color-gray-7)',
              border: '1px solid var(--sl-color-gray-5)',
              borderRadius: '4px',
              cursor: 'pointer',
              color: 'var(--sl-color-text)',
              whiteSpace: 'nowrap',
            }}
          >
            {ex.label}
          </button>
        ))}
      </div>
    </div>
  );
}
