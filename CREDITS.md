# Credits

## htslib and bcftools

**[htslib](https://github.com/samtools/htslib)** and **[bcftools](https://github.com/samtools/bcftools)**  
Created and maintained by the Wellcome Sanger Institute. Primary authors: Heng Li (original author), Petr Danecek (bcftools lead), and hundreds of contributors over 15+ years.

vcfkit owes htslib/bcftools a great deal:

- The VCF/BCF format semantics and edge-case handling that vcfkit implements are documented almost entirely through bcftools source code and behavior.
- vcfkit's correctness is validated via differential tests against `bcftools norm`, `bcftools view`, and `bcftools filter`. Any divergence is treated as a bug in vcfkit.
- The fast-path raw-line parsing strategy in `filter` is the same technique htslib applies in C.

## noodles

**[noodles](https://github.com/zaeleus/noodles)** by Michael Macias  
Pure-Rust I/O for VCF, BCF, FASTA, BAM, and chain files. vcfkit uses noodles for all record parsing, writing, and format auto-detection. Without noodles this project would not exist in its current form.

## Normalization algorithm

**Tan G, Abecasis GR, Kang HM (2015).** "Unified representation of genetic variants."  
*Bioinformatics* 31(13):2202–2204. DOI: [10.1093/bioinformatics/btv112](https://doi.org/10.1093/bioinformatics/btv112)

The left-alignment and parsimony algorithm implemented in `vcfkit normalize`.

## UCSC Genome Browser

**[UCSC Genome Browser](https://genome.ucsc.edu/)**  
Chain files and reference FASTAs used by `vcfkit liftover` (hg19 ↔ hg38 ↔ T2T-CHM13).

## VCF specification

**[VCFv4.3 specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)**  
Maintained at [samtools.github.io/hts-specs](https://samtools.github.io/hts-specs/). The authoritative reference for format semantics.

## Rust dependencies

| Crate | Version | Purpose |
|---|---|---|
| [clap](https://github.com/clap-rs/clap) | v4 | CLI argument parsing with derive macros and shell completions |
| [nom](https://github.com/rust-bakery/nom) | v7 | Parser combinator for the `filter` expression language |
| [anyhow](https://github.com/dtolnay/anyhow) / [thiserror](https://github.com/dtolnay/thiserror) | — | Error handling: anyhow in the CLI, thiserror in the library |
| [tracing](https://github.com/tokio-rs/tracing) | — | Structured, leveled diagnostic logging |
| [owo-colors](https://github.com/jam1garner/owo-colors) | — | Colored terminal output with no_color support |
| [indicatif](https://github.com/console-rs/indicatif) | — | Progress bars with ETA; hidden automatically when output is piped |
| [serde](https://serde.rs/) / [serde_json](https://github.com/serde-rs/json) | — | Serialization for opt-in telemetry payloads |

## AI assistance

Significant portions of this code were written with AI coding assistants (Anthropic Claude). Every change was reviewed by a human and validated against bcftools via differential tests on 1000 Genomes data.
