import { useEffect, useRef } from 'react';
import { EditorState } from '@codemirror/state';
import {
  EditorView,
  keymap,
  lineNumbers,
  highlightActiveLine,
  placeholder as placeholderExt,
} from '@codemirror/view';
import { defaultKeymap, history, historyKeymap } from '@codemirror/commands';
import { oneDark } from '@codemirror/theme-one-dark';

interface Props {
  value: string;
  onChange?: (v: string) => void;
  placeholder?: string;
  readOnly?: boolean;
  minHeight?: string;
}

export function VcfEditor({
  value,
  onChange,
  placeholder = 'Paste VCF content here, or load an example above…',
  readOnly = false,
  minHeight = '280px',
}: Props) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewRef = useRef<EditorView | null>(null);
  const onChangeRef = useRef(onChange);
  onChangeRef.current = onChange;

  useEffect(() => {
    if (!containerRef.current) return;

    const state = EditorState.create({
      doc: value,
      extensions: [
        lineNumbers(),
        highlightActiveLine(),
        history(),
        keymap.of([...defaultKeymap, ...historyKeymap]),
        EditorView.editable.of(!readOnly),
        EditorView.lineWrapping,
        placeholderExt(placeholder),
        EditorView.theme({
          '&': { fontSize: '12px', minHeight },
          '.cm-scroller': {
            fontFamily: 'ui-monospace, "Cascadia Code", "Fira Mono", monospace',
            overflow: 'auto',
          },
          '.cm-content': { padding: '8px 0' },
          '.cm-line': { padding: '0 8px' },
        }),
        oneDark,
        EditorView.updateListener.of((update) => {
          if (update.docChanged) {
            onChangeRef.current?.(update.state.doc.toString());
          }
        }),
      ],
    });

    viewRef.current = new EditorView({ state, parent: containerRef.current });
    return () => {
      viewRef.current?.destroy();
      viewRef.current = null;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Sync external value changes (e.g. "Load example" button)
  useEffect(() => {
    const view = viewRef.current;
    if (!view) return;
    const current = view.state.doc.toString();
    if (current !== value) {
      view.dispatch({
        changes: { from: 0, to: current.length, insert: value },
      });
    }
  }, [value]);

  return (
    <div
      ref={containerRef}
      className="vcf-editor-wrap"
      style={{ border: '1px solid var(--sl-color-gray-5)', borderRadius: '6px', overflow: 'hidden' }}
    />
  );
}
