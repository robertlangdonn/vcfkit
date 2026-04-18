use criterion::{black_box, criterion_group, criterion_main, Criterion};
use vcfkit_core::liftover::{liftover, LiftoverOptions};

// ── inline reference FASTAs ───────────────────────────────────────────────────
//
// Source reference: chr1 120 bp (same as normalize bench / corpus mini_ref).
const SRC_FA: &str = ">chr1\n\
TTCGAATCGAACGTACGTGCCATAATCGACGTACGTATCGTTCGAATCGATCGAATCGAT\n\
CGTTACGTATCGAATCGATCGTTACGTATCGAATCGATCGTTACGTATCGAATCGATCGG\n";
const SRC_FAI: &str = "chr1\t120\t6\t60\t61\n";

// Target reference: chr2 200 bp — all 'A' so any SNP REF lifted here matches.
// We set the first 120 bases to match chr1 exactly (identity mapping test).
const TGT_FA: &str = ">chr2\n\
TTCGAATCGAACGTACGTGCCATAATCGACGTACGTATCGTTCGAATCGATCGAATCGAT\n\
CGTTACGTATCGAATCGATCGTTACGTATCGAATCGATCGTTACGTATCGAATCGATCGG\n\
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\
AAAAAAAAAAAAAAAAAAA\n";
const TGT_FAI: &str = "chr2\t200\t6\t60\t61\n";

// ── inline chain fixture ──────────────────────────────────────────────────────
//
// Identity mapping: chr1:0-100 → chr2:0-100 (+ strand, single 100-bp block).
const CHAIN: &str = "chain 1000 chr1 120 + 0 100 chr2 200 + 0 100 1\n100\n\n";

// ── inline VCF fixture (50 SNP records on chr1 within 0-100) ─────────────────
const VCF_50: &str = "\
##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=120>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t2\t.\tT\tA\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t3\t.\tC\tG\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t5\t.\tA\tT\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t7\t.\tT\tC\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t9\t.\tC\tA\t50\tPASS\tAF=0.3\tGT\t0/1\n\
chr1\t11\t.\tG\tT\t55\tPASS\tAF=0.4\tGT\t1/1\n\
chr1\t13\t.\tA\tC\t60\tPASS\tAF=0.6\tGT\t0/1\n\
chr1\t15\t.\tC\tT\t45\tPASS\tAF=0.2\tGT\t0/1\n\
chr1\t17\t.\tG\tA\t70\tPASS\tAF=0.7\tGT\t1/1\n\
chr1\t19\t.\tT\tG\t35\tPASS\tAF=0.1\tGT\t0/1\n\
chr1\t21\t.\tA\tT\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t23\t.\tT\tC\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t25\t.\tA\tG\t50\tPASS\tAF=0.4\tGT\t0/1\n\
chr1\t27\t.\tG\tC\t65\tPASS\tAF=0.6\tGT\t1/1\n\
chr1\t29\t.\tT\tA\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t31\t.\tC\tT\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t33\t.\tA\tG\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t35\t.\tT\tC\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t37\t.\tG\tA\t50\tPASS\tAF=0.3\tGT\t0/1\n\
chr1\t39\t.\tT\tG\t55\tPASS\tAF=0.4\tGT\t0/1\n\
chr1\t41\t.\tA\tC\t60\tPASS\tAF=0.2\tGT\t0/1\n\
chr1\t43\t.\tG\tT\t45\tPASS\tAF=0.6\tGT\t1/1\n\
chr1\t45\t.\tT\tA\t70\tPASS\tAF=0.7\tGT\t0/1\n\
chr1\t47\t.\tC\tG\t35\tPASS\tAF=0.1\tGT\t0/1\n\
chr1\t49\t.\tA\tT\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t51\t.\tT\tC\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t53\t.\tG\tA\t50\tPASS\tAF=0.4\tGT\t0/1\n\
chr1\t55\t.\tC\tT\t65\tPASS\tAF=0.6\tGT\t1/1\n\
chr1\t57\t.\tA\tG\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t59\t.\tT\tA\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t61\t.\tC\tG\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t63\t.\tG\tC\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t65\t.\tT\tG\t50\tPASS\tAF=0.3\tGT\t0/1\n\
chr1\t67\t.\tA\tC\t55\tPASS\tAF=0.4\tGT\t0/1\n\
chr1\t69\t.\tG\tT\t60\tPASS\tAF=0.2\tGT\t0/1\n\
chr1\t71\t.\tC\tA\t45\tPASS\tAF=0.6\tGT\t1/1\n\
chr1\t73\t.\tT\tG\t70\tPASS\tAF=0.7\tGT\t0/1\n\
chr1\t75\t.\tA\tC\t35\tPASS\tAF=0.1\tGT\t0/1\n\
chr1\t77\t.\tG\tA\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t79\t.\tT\tC\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t81\t.\tA\tT\t50\tPASS\tAF=0.4\tGT\t0/1\n\
chr1\t83\t.\tC\tG\t65\tPASS\tAF=0.6\tGT\t1/1\n\
chr1\t85\t.\tG\tA\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t87\t.\tT\tG\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t89\t.\tA\tC\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t91\t.\tT\tA\t50\tPASS\tAF=0.5\tGT\t0/1\n\
chr1\t93\t.\tC\tT\t50\tPASS\tAF=0.3\tGT\t0/1\n\
chr1\t95\t.\tG\tC\t55\tPASS\tAF=0.4\tGT\t0/1\n\
chr1\t97\t.\tA\tG\t60\tPASS\tAF=0.2\tGT\t0/1\n\
chr1\t99\t.\tT\tA\t45\tPASS\tAF=0.6\tGT\t1/1\n\
";

fn write_files(
    dir: &std::path::Path,
) -> (std::path::PathBuf, std::path::PathBuf, std::path::PathBuf) {
    std::fs::create_dir_all(dir).expect("create temp dir");
    let chain_path = dir.join("bench.chain");
    let src_fa_path = dir.join("src.fa");
    let tgt_fa_path = dir.join("tgt.fa");
    std::fs::write(&chain_path, CHAIN).expect("write chain");
    std::fs::write(&src_fa_path, SRC_FA).expect("write src.fa");
    std::fs::write(dir.join("src.fa.fai"), SRC_FAI).expect("write src.fa.fai");
    std::fs::write(&tgt_fa_path, TGT_FA).expect("write tgt.fa");
    std::fs::write(dir.join("tgt.fa.fai"), TGT_FAI).expect("write tgt.fa.fai");
    (chain_path, src_fa_path, tgt_fa_path)
}

fn bench_liftover(c: &mut Criterion) {
    let dir = std::env::temp_dir().join("vcfkit_bench_liftover");
    let (chain_path, src_fa_path, tgt_fa_path) = write_files(&dir);

    c.bench_function("liftover_50_records", |b| {
        b.iter(|| {
            let cursor = std::io::Cursor::new(black_box(VCF_50.as_bytes()));
            let mut output = Vec::new();
            liftover(
                cursor,
                &mut output,
                black_box(&chain_path),
                black_box(&src_fa_path),
                Some(black_box(tgt_fa_path.as_path())),
                LiftoverOptions::default(),
            )
            .unwrap()
        })
    });
}

criterion_group!(benches, bench_liftover);
criterion_main!(benches);
