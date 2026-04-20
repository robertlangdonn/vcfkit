import { useState } from 'react';
import { cliCommand } from '../lib/format-utils';

interface Props {
  operation: 'filter' | 'normalize' | 'liftover';
  expression?: string;
}

export function CliEquivalent({ operation, expression }: Props) {
  const [copied, setCopied] = useState(false);
  const cmd = cliCommand(operation, expression);

  const handleCopy = async () => {
    await navigator.clipboard.writeText(cmd);
    setCopied(true);
    setTimeout(() => setCopied(false), 1500);
  };

  return (
    <div className="cli-bar">
      <span className="cli-label">CLI equivalent</span>
      <code className="cli-code">{cmd}</code>
      <button className="cli-copy" onClick={handleCopy} title="Copy to clipboard" aria-label="Copy command">
        {copied ? (
          <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" aria-hidden="true"><polyline points="20 6 9 17 4 12"/></svg>
        ) : (
          <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" aria-hidden="true"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg>
        )}
      </button>
    </div>
  );
}
