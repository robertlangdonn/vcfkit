export function countRecords(vcf: string): number {
  return vcf.split('\n').filter((l) => l && !l.startsWith('#')).length;
}

export function formatDuration(ms: number): string {
  if (ms < 1000) return `${Math.round(ms)}ms`;
  return `${(ms / 1000).toFixed(1)}s`;
}

export function downloadVcf(content: string, filename = 'vcfkit-output.vcf') {
  const blob = new Blob([content], { type: 'text/plain' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

export function cliCommand(op: 'filter' | 'normalize' | 'liftover', expression = ''): string {
  switch (op) {
    case 'filter':
      return `vcfkit filter -e "${expression || 'INFO/AF < 0.01'}" input.vcf`;
    case 'normalize':
      return `vcfkit normalize -f reference.fasta input.vcf`;
    case 'liftover':
      return `vcfkit liftover -s hg19.fa -t hg38.fa -c hg19ToHg38.over.chain.gz input.vcf`;
  }
}

export const MAX_DEMO_RECORDS = 10_000;

export function truncateForDemo(vcf: string): { text: string; truncated: boolean } {
  const lines = vcf.split('\n');
  const headers: string[] = [];
  const records: string[] = [];
  for (const line of lines) {
    if (line.startsWith('#')) headers.push(line);
    else if (line) records.push(line);
  }
  if (records.length <= MAX_DEMO_RECORDS) {
    return { text: vcf, truncated: false };
  }
  return {
    text: [...headers, ...records.slice(0, MAX_DEMO_RECORDS)].join('\n') + '\n',
    truncated: true,
  };
}
