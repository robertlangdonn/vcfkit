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
      <button className="cli-copy" onClick={handleCopy} title="Copy to clipboard">
        {copied ? '✓' : 'Copy'}
      </button>
    </div>
  );
}
