# Synthetic VCF Test Corpus

Hand-crafted VCF 4.2 fixtures for differential testing of vcfkit operations (normalize, liftover, filter).

| File | Description |
|------|-------------|
| `basic.vcf` | 5 clean biallelic SNPs on chr1/chr2, FILTER=PASS, QUAL>30 — the happy path |
| `multi_allelic.vcf` | 5 records with 2–3 ALT alleles; includes `Number=A` (AF) and `Number=R` (AD) INFO fields for split testing |
| `indels_unnormalized.vcf` | 5 indels written in right-aligned form against `mini_ref.fa`; NORM_* INFO tags record expected left-aligned form |
| `mini_ref.fa` | 120 bp FASTA reference (chr1, repeating ACGT) used by normalization tests |
| `mini_ref.fa.fai` | SAMtools-style index for `mini_ref.fa` |
| `ref_mismatch.vcf` | 3 variants where REF does not match the reference genome; ACTUAL_REF tag records the correct base |
| `missing_fields.vcf` | Variants with `.` QUAL, `.` FILTER, `.` INFO subfields, and missing FORMAT sample values |
| `empty_alt.vcf` | 6 structural variants with symbolic ALT alleles (`<DEL>`, `<DUP>`, `<INS>`, `<CNV>`); must pass through normalization unchanged |
| `mixed_filters.vcf` | 10 variants covering FILTER values PASS, LowQual, DP5, and LowQual;DP5 (multi-value) |
| `multiline_info.vcf` | 12 variants with rich INFO (AF, DP, MQ, CSQ annotation) for filter expression testing |
| `chr_prefix.vcf` | Same variants written with both `chr1`/`chr2`/`chrX` and `1`/`2`/`X` contig styles |
| `large_indel.vcf` | 50 bp insertion and deletion on chr1 and chr2 for normalization edge cases and benchmarking |
| `hg19_coords.vcf` | 5 variants at known hg19 positions on chr1/chr17/chrX for liftover testing |
