//! Evaluation harness for the --english filter translation feature.
//!
//! Loads tests/english_filter_corpus.yaml, runs each query through the
//! Anthropic API, and reports how many accepted expressions were matched.
//!
//! Usage:
//!   ANTHROPIC_API_KEY=sk-ant-... VCFKIT_EVAL_CONFIRM=1 \
//!     cargo run --bin english_eval
//!
//! Gated behind VCFKIT_EVAL_CONFIRM=1 to prevent accidental API spend.

use std::time::Duration;

use serde::Deserialize;

#[path = "../english.rs"]
mod english;

#[derive(Debug, Deserialize)]
struct CorpusEntry {
    query: String,
    accepted_expressions: Vec<String>,
    #[serde(default)]
    notes: String,
}

#[tokio::main]
async fn main() {
    let corpus_path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/english_filter_corpus.yaml"
    );
    let corpus_text = std::fs::read_to_string(corpus_path).expect("failed to read corpus");
    let entries: Vec<CorpusEntry> =
        serde_yaml::from_str(&corpus_text).expect("failed to parse corpus YAML");

    let estimated_cost = entries.len() as f64 * 0.001;
    eprintln!(
        "Corpus: {} cases. Estimated cost: ${:.2}.",
        entries.len(),
        estimated_cost
    );
    eprintln!(
        "Model:  {} (override with VCFKIT_LLM_MODEL)",
        english::DEFAULT_MODEL
    );

    let confirmed = std::env::var("VCFKIT_EVAL_CONFIRM")
        .map(|v| v == "1")
        .unwrap_or(false);
    if !confirmed {
        eprintln!("\nSet VCFKIT_EVAL_CONFIRM=1 to proceed.");
        eprintln!("Waiting 5s for Ctrl-C...");
        std::thread::sleep(Duration::from_secs(5));
        eprintln!("Proceeding.");
    }

    // Provide a representative schema so the model uses real field names.
    let schema = english::HeaderSchema {
        info_fields: vec![
            field("AF", "A", "Float", "Allele frequency"),
            field("DP", "1", "Integer", "Total read depth"),
            field(
                "CSQ",
                ".",
                "String",
                "VEP consequence annotations (pipe-separated)",
            ),
            field("GENE", "1", "String", "Gene symbol"),
            field("INDEL", "0", "Flag", "Variant is an indel"),
            field("SVTYPE", "1", "String", "Structural variant type"),
            field("VT", "1", "String", "Variant type (SNP, INDEL, SV, etc.)"),
            field(
                "LOF",
                ".",
                "String",
                "Loss-of-function annotation (HC=high-confidence)",
            ),
            field("MULTI_ALLELIC", "0", "Flag", "Site is multi-allelic"),
        ],
        format_fields: vec![
            field("GT", "1", "String", "Genotype"),
            field("DP", "1", "Integer", "Read depth"),
            field("GQ", "1", "Integer", "Genotype quality"),
            field("AD", "R", "Integer", "Allelic depths"),
        ],
        contigs: vec!["chr1".into(), "chr2".into(), "chr17".into(), "chr22".into()],
    };

    let mut passed = 0usize;
    let mut failed = 0usize;
    let mut errored = 0usize;

    for (i, entry) in entries.iter().enumerate() {
        eprint!("[{}/{}] \"{}\" ... ", i + 1, entries.len(), entry.query);

        match english::translate(&entry.query, &schema).await {
            Ok(t) => {
                let got = normalise(&t.expression);
                let matched = entry
                    .accepted_expressions
                    .iter()
                    .any(|a| normalise(a) == got);
                if matched {
                    eprintln!("PASS");
                    eprintln!(
                        "    expr: {} (conf={:.0}%, model={})",
                        t.expression,
                        t.confidence * 100.0,
                        t.model
                    );
                    passed += 1;
                } else {
                    eprintln!("FAIL");
                    eprintln!("    got:       {}", t.expression);
                    eprintln!("    accepted:  {}", entry.accepted_expressions.join(" | "));
                    eprintln!("    reasoning: {}", t.reasoning);
                    if !t.caveats.is_empty() {
                        eprintln!("    caveats:   {}", t.caveats.join("; "));
                    }
                    if !entry.notes.is_empty() {
                        eprintln!("    notes:     {}", entry.notes);
                    }
                    failed += 1;
                }
            }
            Err(e) => {
                eprintln!("ERROR: {e}");
                errored += 1;
            }
        }
    }

    let total = entries.len();
    let pct = (passed * 100).checked_div(total).unwrap_or(0);
    let target_cases = (total * 85).div_ceil(100);
    eprintln!();
    eprintln!("Results: {passed}/{total} passed ({pct}%), {failed} failed, {errored} errored");
    eprintln!("Target:  >=85% ({target_cases} cases needed)");

    if pct < 85 {
        eprintln!("BELOW TARGET — review failed cases and iterate on the system prompt.");
        std::process::exit(1);
    } else {
        eprintln!("AT OR ABOVE TARGET.");
    }
}

fn field(id: &str, number: &str, ty: &str, description: &str) -> english::FieldDef {
    english::FieldDef {
        id: id.into(),
        number: number.into(),
        ty: ty.into(),
        description: description.into(),
    }
}

fn normalise(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}
