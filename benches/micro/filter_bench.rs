use criterion::{black_box, criterion_group, criterion_main, Criterion};
use vcfkit_core::filter::{filter, FilterExpression, FilterOptions};

// ── inline VCF fixture (50 records with INFO/AF fields) ──────────────────────
//
// Half the records have AF < 0.05 (rare variants) and half have AF >= 0.05.
// This exercises both the matching and non-matching paths through the filter.
const VCF_50: &str = "\
##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=248956422>\n\
##contig=<ID=chr2,length=242193529>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">\n\
##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=\"VQSLOD score\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t100\trs001\tA\tT\t50\tPASS\tAF=0.01;DP=100;VQSLOD=5.2\tGT:GQ\t0/1:45\n\
chr1\t200\trs002\tC\tG\t60\tPASS\tAF=0.50;DP=120;VQSLOD=8.1\tGT:GQ\t1/1:55\n\
chr1\t300\trs003\tG\tA\t75\tPASS\tAF=0.02;DP=80;VQSLOD=3.4\tGT:GQ\t0/1:70\n\
chr1\t400\trs004\tT\tC\t45\tPASS\tAF=0.35;DP=90;VQSLOD=7.2\tGT:GQ\t0/1:40\n\
chr1\t500\trs005\tA\tG\t90\tPASS\tAF=0.04;DP=200;VQSLOD=4.8\tGT:GQ\t1/1:85\n\
chr1\t600\trs006\tC\tT\t55\tPASS\tAF=0.60;DP=150;VQSLOD=9.0\tGT:GQ\t0/1:50\n\
chr1\t700\trs007\tG\tC\t65\tPASS\tAF=0.008;DP=60;VQSLOD=2.1\tGT:GQ\t0/1:60\n\
chr1\t800\trs008\tT\tA\t40\tPASS\tAF=0.45;DP=110;VQSLOD=6.5\tGT:GQ\t1/1:35\n\
chr1\t900\trs009\tA\tC\t70\tPASS\tAF=0.03;DP=130;VQSLOD=4.2\tGT:GQ\t0/1:65\n\
chr1\t1000\trs010\tC\tG\t80\tPASS\tAF=0.55;DP=170;VQSLOD=8.7\tGT:GQ\t1/1:75\n\
chr1\t1100\trs011\tG\tT\t50\tPASS\tAF=0.015;DP=75;VQSLOD=3.8\tGT:GQ\t0/1:45\n\
chr1\t1200\trs012\tT\tC\t60\tPASS\tAF=0.40;DP=95;VQSLOD=7.5\tGT:GQ\t0/1:55\n\
chr1\t1300\trs013\tA\tG\t75\tPASS\tAF=0.025;DP=85;VQSLOD=3.1\tGT:GQ\t0/1:70\n\
chr1\t1400\trs014\tC\tA\t45\tPASS\tAF=0.30;DP=140;VQSLOD=6.8\tGT:GQ\t1/1:40\n\
chr1\t1500\trs015\tG\tC\t90\tPASS\tAF=0.03;DP=190;VQSLOD=4.5\tGT:GQ\t0/1:85\n\
chr2\t100\trs016\tA\tT\t55\tPASS\tAF=0.70;DP=145;VQSLOD=9.2\tGT:GQ\t1/1:50\n\
chr2\t200\trs017\tC\tG\t65\tPASS\tAF=0.01;DP=65;VQSLOD=2.5\tGT:GQ\t0/1:60\n\
chr2\t300\trs018\tG\tA\t40\tPASS\tAF=0.45;DP=105;VQSLOD=7.0\tGT:GQ\t1/1:35\n\
chr2\t400\trs019\tT\tC\t70\tPASS\tAF=0.04;DP=125;VQSLOD=4.0\tGT:GQ\t0/1:65\n\
chr2\t500\trs020\tA\tG\t80\tPASS\tAF=0.55;DP=165;VQSLOD=8.5\tGT:GQ\t1/1:75\n\
chr2\t600\trs021\tC\tT\t50\tPASS\tAF=0.009;DP=70;VQSLOD=2.8\tGT:GQ\t0/1:45\n\
chr2\t700\trs022\tG\tC\t60\tPASS\tAF=0.38;DP=90;VQSLOD=7.3\tGT:GQ\t0/1:55\n\
chr2\t800\trs023\tT\tA\t75\tPASS\tAF=0.02;DP=80;VQSLOD=3.6\tGT:GQ\t0/1:70\n\
chr2\t900\trs024\tA\tC\t45\tPASS\tAF=0.32;DP=135;VQSLOD=6.6\tGT:GQ\t1/1:40\n\
chr2\t1000\trs025\tC\tG\t90\tPASS\tAF=0.035;DP=185;VQSLOD=4.3\tGT:GQ\t0/1:85\n\
chr2\t1100\trs026\tG\tT\t55\tPASS\tAF=0.65;DP=140;VQSLOD=9.1\tGT:GQ\t1/1:50\n\
chr2\t1200\trs027\tT\tC\t65\tPASS\tAF=0.012;DP=60;VQSLOD=2.3\tGT:GQ\t0/1:60\n\
chr2\t1300\trs028\tA\tG\t40\tPASS\tAF=0.42;DP=100;VQSLOD=7.1\tGT:GQ\t1/1:35\n\
chr2\t1400\trs029\tC\tA\t70\tPASS\tAF=0.045;DP=120;VQSLOD=4.1\tGT:GQ\t0/1:65\n\
chr2\t1500\trs030\tG\tC\t80\tPASS\tAF=0.52;DP=160;VQSLOD=8.4\tGT:GQ\t1/1:75\n\
chr1\t2000\trs031\tA\tT\t50\tPASS\tAF=0.007;DP=68;VQSLOD=2.2\tGT:GQ\t0/1:45\n\
chr1\t2100\trs032\tC\tG\t60\tPASS\tAF=0.36;DP=88;VQSLOD=7.4\tGT:GQ\t0/1:55\n\
chr1\t2200\trs033\tG\tA\t75\tPASS\tAF=0.022;DP=82;VQSLOD=3.3\tGT:GQ\t0/1:70\n\
chr1\t2300\trs034\tT\tC\t45\tPASS\tAF=0.28;DP=132;VQSLOD=6.4\tGT:GQ\t1/1:40\n\
chr1\t2400\trs035\tA\tG\t90\tPASS\tAF=0.038;DP=182;VQSLOD=4.4\tGT:GQ\t0/1:85\n\
chr1\t2500\trs036\tC\tT\t55\tPASS\tAF=0.62;DP=138;VQSLOD=9.0\tGT:GQ\t1/1:50\n\
chr1\t2600\trs037\tG\tC\t65\tPASS\tAF=0.014;DP=58;VQSLOD=2.6\tGT:GQ\t0/1:60\n\
chr1\t2700\trs038\tT\tA\t40\tPASS\tAF=0.44;DP=98;VQSLOD=7.2\tGT:GQ\t1/1:35\n\
chr1\t2800\trs039\tA\tC\t70\tPASS\tAF=0.042;DP=118;VQSLOD=4.2\tGT:GQ\t0/1:65\n\
chr1\t2900\trs040\tC\tG\t80\tPASS\tAF=0.50;DP=158;VQSLOD=8.3\tGT:GQ\t1/1:75\n\
chr2\t2000\trs041\tG\tT\t50\tPASS\tAF=0.006;DP=66;VQSLOD=2.0\tGT:GQ\t0/1:45\n\
chr2\t2100\trs042\tT\tC\t60\tPASS\tAF=0.34;DP=86;VQSLOD=7.5\tGT:GQ\t0/1:55\n\
chr2\t2200\trs043\tA\tG\t75\tPASS\tAF=0.018;DP=78;VQSLOD=3.2\tGT:GQ\t0/1:70\n\
chr2\t2300\trs044\tC\tA\t45\tPASS\tAF=0.26;DP=128;VQSLOD=6.2\tGT:GQ\t1/1:40\n\
chr2\t2400\trs045\tG\tC\t90\tPASS\tAF=0.033;DP=178;VQSLOD=4.6\tGT:GQ\t0/1:85\n\
chr2\t2500\trs046\tT\tA\t55\tPASS\tAF=0.60;DP=136;VQSLOD=8.9\tGT:GQ\t1/1:50\n\
chr2\t2600\trs047\tA\tC\t65\tPASS\tAF=0.016;DP=56;VQSLOD=2.7\tGT:GQ\t0/1:60\n\
chr2\t2700\trs048\tC\tG\t40\tPASS\tAF=0.46;DP=96;VQSLOD=7.3\tGT:GQ\t1/1:35\n\
chr2\t2800\trs049\tG\tT\t70\tPASS\tAF=0.048;DP=116;VQSLOD=4.7\tGT:GQ\t0/1:65\n\
chr2\t2900\trs050\tT\tC\t80\tPASS\tAF=0.48;DP=156;VQSLOD=8.2\tGT:GQ\t1/1:75\n\
";

fn bench_filter_af(c: &mut Criterion) {
    let expr = FilterExpression::parse("INFO/AF < 0.05").expect("parse filter expr");

    c.bench_function("filter_50_records_af_lt_0.05", |b| {
        b.iter(|| {
            let cursor = std::io::Cursor::new(black_box(VCF_50.as_bytes()));
            let mut output = Vec::new();
            filter(
                cursor,
                &mut output,
                black_box(expr.clone()),
                FilterOptions::default(),
            )
            .unwrap()
        })
    });
}

fn bench_filter_compound(c: &mut Criterion) {
    let expr =
        FilterExpression::parse("INFO/AF < 0.05 && INFO/DP > 50").expect("parse filter expr");

    c.bench_function("filter_50_records_compound_expr", |b| {
        b.iter(|| {
            let cursor = std::io::Cursor::new(black_box(VCF_50.as_bytes()));
            let mut output = Vec::new();
            filter(
                cursor,
                &mut output,
                black_box(expr.clone()),
                FilterOptions::default(),
            )
            .unwrap()
        })
    });
}

criterion_group!(benches, bench_filter_af, bench_filter_compound);
criterion_main!(benches);
