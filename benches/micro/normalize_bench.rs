use criterion::{black_box, criterion_group, criterion_main, Criterion};
use vcfkit_core::normalize::{normalize, NormalizeOptions, RefCheck};

// ── inline FASTA fixture ─────────────────────────────────────────────────────
//
// chr1 is 120 bp (matches the corpus mini_ref.fa).
const MINI_REF_FA: &str = ">chr1\n\
TTCGAATCGAACGTACGTGCCATAATCGACGTACGTATCGTTCGAATCGATCGAATCGAT\n\
CGTTACGTATCGAATCGATCGTTACGTATCGAATCGATCGTTACGTATCGAATCGATCGG\n";

// fai index for the above: <name> <len> <offset> <bases/line> <bytes/line>
// offset 6 (len of ">chr1\n"), bases 60, bytes 61
const MINI_REF_FAI: &str = "chr1\t120\t6\t60\t61\n";

// ── inline VCF fixture (50 records) ──────────────────────────────────────────
//
// Mix of SNPs, biallelic indels (multi-allelic), and records with INFO fields.
// All positions are within the 120-bp chr1 above.
const VCF_50: &str = "\
##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t2\t.\tT\tA\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t3\t.\tC\tG\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t5\t.\tA\tT\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t7\t.\tT\tC\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t9\t.\tC\tA\t50\tPASS\tAF=0.3;DP=25\tGT\t0/1\n\
chr1\t11\t.\tG\tT\t55\tPASS\tAF=0.4;DP=40\tGT\t1/1\n\
chr1\t13\t.\tA\tC\t60\tPASS\tAF=0.6;DP=50\tGT\t0/1\n\
chr1\t15\t.\tC\tT\t45\tPASS\tAF=0.2;DP=20\tGT\t0/1\n\
chr1\t17\t.\tG\tA\t70\tPASS\tAF=0.7;DP=60\tGT\t1/1\n\
chr1\t19\t.\tT\tG\t35\tPASS\tAF=0.1;DP=15\tGT\t0/1\n\
chr1\t21\t.\tA\tT\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t23\t.\tT\tC\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t25\t.\tA\tG\t50\tPASS\tAF=0.4;DP=35\tGT\t0/1\n\
chr1\t27\t.\tG\tC\t65\tPASS\tAF=0.6;DP=45\tGT\t1/1\n\
chr1\t29\t.\tT\tA\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t31\t.\tC\tT\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t33\t.\tA\tG\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t35\t.\tT\tC\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t37\t.\tG\tA\t50\tPASS\tAF=0.3;DP=28\tGT\t0/1\n\
chr1\t39\t.\tT\tG\t55\tPASS\tAF=0.4;DP=32\tGT\t0/1\n\
chr1\t41\t.\tA\tC\t60\tPASS\tAF=0.2;DP=18\tGT\t0/1\n\
chr1\t43\t.\tG\tT\t45\tPASS\tAF=0.6;DP=48\tGT\t1/1\n\
chr1\t45\t.\tT\tA\t70\tPASS\tAF=0.7;DP=55\tGT\t0/1\n\
chr1\t47\t.\tC\tG\t35\tPASS\tAF=0.1;DP=10\tGT\t0/1\n\
chr1\t49\t.\tA\tT\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t51\t.\tT\tC\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t53\t.\tG\tA\t50\tPASS\tAF=0.4;DP=35\tGT\t0/1\n\
chr1\t55\t.\tC\tT\t65\tPASS\tAF=0.6;DP=45\tGT\t1/1\n\
chr1\t57\t.\tA\tG\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t59\t.\tT\tA\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t61\t.\tC\tG\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t63\t.\tG\tC\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t65\t.\tT\tG\t50\tPASS\tAF=0.3;DP=28\tGT\t0/1\n\
chr1\t67\t.\tA\tC\t55\tPASS\tAF=0.4;DP=32\tGT\t0/1\n\
chr1\t69\t.\tG\tT\t60\tPASS\tAF=0.2;DP=18\tGT\t0/1\n\
chr1\t71\t.\tC\tA\t45\tPASS\tAF=0.6;DP=48\tGT\t1/1\n\
chr1\t73\t.\tT\tG\t70\tPASS\tAF=0.7;DP=55\tGT\t0/1\n\
chr1\t75\t.\tA\tC\t35\tPASS\tAF=0.1;DP=10\tGT\t0/1\n\
chr1\t77\t.\tG\tA\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t79\t.\tT\tC\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t81\t.\tA\tT\t50\tPASS\tAF=0.4;DP=35\tGT\t0/1\n\
chr1\t83\t.\tC\tG\t65\tPASS\tAF=0.6;DP=45\tGT\t1/1\n\
chr1\t85\t.\tG\tA\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t87\t.\tT\tG\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t89\t.\tA\tC\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t91\t.\tT\tA\t50\tPASS\tAF=0.5;DP=30\tGT\t0/1\n\
chr1\t93\t.\tC\tT\t50\tPASS\tAF=0.3;DP=28\tGT\t0/1\n\
chr1\t95\t.\tG\tC\t55\tPASS\tAF=0.4;DP=32\tGT\t0/1\n\
chr1\t97\t.\tA\tG\t60\tPASS\tAF=0.2;DP=18\tGT\t0/1\n\
chr1\t99\t.\tT\tA\t45\tPASS\tAF=0.6;DP=48\tGT\t1/1\n\
";

/// Write the inline FASTA and FAI to a temporary directory, return the fa path.
fn write_mini_ref(dir: &std::path::Path) -> std::path::PathBuf {
    let fa_path = dir.join("mini_ref.fa");
    let fai_path = dir.join("mini_ref.fa.fai");
    std::fs::write(&fa_path, MINI_REF_FA).expect("write mini_ref.fa");
    std::fs::write(&fai_path, MINI_REF_FAI).expect("write mini_ref.fa.fai");
    fa_path
}

fn bench_normalize_snps(c: &mut Criterion) {
    let dir = std::env::temp_dir().join("vcfkit_bench_normalize");
    std::fs::create_dir_all(&dir).expect("create temp dir");
    let ref_path = write_mini_ref(&dir);

    c.bench_function("normalize_50_snp_records", |b| {
        b.iter(|| {
            let cursor = std::io::Cursor::new(black_box(VCF_50.as_bytes()));
            let mut output = Vec::new();
            normalize(
                cursor,
                &mut output,
                black_box(&ref_path),
                NormalizeOptions {
                    split_multiallelics: true,
                    left_align: true,
                    check_ref: RefCheck::Ignore,
                    output_format: vcfkit_core::io::OutputFormat::Vcf,
                },
            )
            .unwrap()
        })
    });
}

criterion_group!(benches, bench_normalize_snps);
criterion_main!(benches);
