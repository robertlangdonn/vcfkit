export interface VcfExample {
  id: string;
  label: string;
  vcf: string;
  expression?: string;
  chain?: string;
}

export const EXAMPLES: Record<'filter' | 'normalize' | 'liftover', VcfExample[]> = {
  filter: [
    {
      id: 'rare-variants',
      label: 'Rare variants (AF < 0.01)',
      expression: "INFO/AF < 0.01",
      vcf: `##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t925952\trs3107975\tC\tT\t50\tPASS\tAF=0.45;DP=120
chr1\t931271\trs9988001\tG\tA\t30\tPASS\tAF=0.006;DP=98
chr1\t941119\trs2298214\tT\tC\t80\tPASS\tAF=0.72;DP=215
chr1\t944296\trs4422948\tC\tA\t90\tPASS\tAF=0.003;DP=188
chr7\t117548628\trs113993960\tA\tT\t60\tPASS\tAF=0.009;DP=142
`,
    },
    {
      id: 'pass-qual',
      label: 'High-quality PASS (QUAL > 30 && FILTER == PASS)',
      expression: "QUAL > 30 && FILTER == 'PASS'",
      vcf: `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t55\tPASS\tDP=200
chr1\t200\t.\tC\tG\t20\tPASS\tDP=80
chr1\t300\t.\tT\tA\t80\tLowQual\tDP=220
chr1\t400\t.\tG\tC\t40\tPASS\tDP=195
chr1\t500\t.\tA\tG\t10\tLowQual\tDP=40
`,
    },
    {
      id: 'brca1',
      label: 'BRCA1 region (CHROM + POS range)',
      expression: "CHROM == 'chr17' && POS >= 43044295 && POS <= 43125483",
      vcf: `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000000\t.\tA\tT\t60\tPASS\tGENE=OTHER
chr17\t43000000\t.\tC\tG\t70\tPASS\tGENE=BRCA1_5prime
chr17\t43050000\trs1234\tT\tA\t90\tPASS\tGENE=BRCA1
chr17\t43080000\t.\tG\tC\t85\tPASS\tGENE=BRCA1
chr17\t43200000\t.\tA\tG\t75\tPASS\tGENE=BRCA1_3prime
chr22\t100000\t.\tC\tT\t65\tPASS\tGENE=OTHER2
`,
    },
    {
      id: 'multiallelic-any',
      label: 'Multi-allelic AF (any-element semantics)',
      expression: "INFO/AF < 0.01",
      vcf: `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency per ALT">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT,C\t50\tPASS\tAF=0.05,0.003
chr1\t200\t.\tG\tA,T\t60\tPASS\tAF=0.12,0.15
chr1\t300\t.\tC\tG,A,T\t70\tPASS\tAF=0.3,0.008,0.2
`,
    },
  ],

  normalize: [
    {
      id: 'multiallelic-split',
      label: 'Split multi-allelic sites',
      vcf: `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT,C\t50\tPASS\tAF=0.05,0.03;DP=200
chr1\t200\t.\tG\tA\t60\tPASS\tAF=0.12;DP=180
chr1\t300\t.\tC\tG,A,T\t70\tPASS\tAF=0.3,0.1,0.2;DP=220
`,
    },
    {
      id: 'number-ar',
      label: 'Number=R INFO fields (allele-indexed)',
      vcf: `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF per ALT">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allele depth (REF + ALTs)">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT,C\t50\tPASS\tAF=0.05,0.03;AD=150,10,6;DP=166
chr1\t200\t.\tG\tA,T,C\t60\tPASS\tAF=0.3,0.1,0.2;AD=120,30,10,20;DP=180
`,
    },
  ],

  liftover: [
    {
      id: 'simple-snv',
      label: 'Simple SNV (paste your chain file)',
      vcf: `##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10000\t.\tA\tT\t.\t.\t.
chr1\t50000\t.\tC\tG\t.\t.\t.
chr2\t30000\t.\tG\tA\t.\t.\t.
`,
      chain: '',
    },
  ],
};
