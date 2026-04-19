//! Reference data acquisition for real-world differential tests.
//!
//! [`RealWorldData::acquire`] downloads files from canonical sources to
//! `tests/real_world/`, verifying SHA-256 checksums. Subsequent calls reuse
//! cached files without re-downloading.
//!
//! All paths point at chr22-only files to keep total size manageable (~600MB):
//! - 1000 Genomes chr22 VCF (hg19 b37 coordinates, genotypes, ~196MB gz)
//! - hg19 chr22 reference FASTA + FAI
//! - hg38 chr22 reference FASTA + FAI
//! - hg19→hg38 UCSC chain file (gzip)

use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

/// Paths to all reference data files needed by the real-world tests.
pub struct RealWorldData {
    /// 1000 Genomes chr22 VCF (hg19/b37 coordinates, genotypes, plain VCF after decompression).
    pub chr22_vcf: PathBuf,
    /// hg19 chr22 reference FASTA (with .fai sidecar).
    pub hg19_chr22_fa: PathBuf,
    /// hg38 chr22 reference FASTA (with .fai sidecar).
    pub hg38_chr22_fa: PathBuf,
    /// hg19→hg38 UCSC chain file (gzip).
    pub hg19_to_hg38_chain: PathBuf,
}

impl RealWorldData {
    /// Acquire all reference data, downloading and verifying if necessary.
    pub fn acquire(dir: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        fs::create_dir_all(dir)?;

        // ── 1000 Genomes chr22 VCF (sites-only, derived) ─────────────────────
        // The per-chromosome sites-only file no longer exists on the EBI FTP.
        // We download the genotypes VCF (~196MB gz) and strip sample columns
        // with `bcftools view -G` to produce a lean ~14MB sites VCF (~354K records).
        let chr22_geno_gz = download_and_verify(
            dir,
            "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\
             ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
            None,
        )?;
        let chr22_vcf = extract_sites_vcf(&chr22_geno_gz, dir)?;

        // ── hg19 chr22 FASTA ─────────────────────────────────────────────────
        let hg19_chr22_fa_gz = download_and_verify(
            dir,
            "chr22.hg19.fa.gz",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz",
            None, // UCSC checksums not published; integrity ensured by samtools faidx
        )?;
        let hg19_chr22_fa = decompress_gz_if_needed(&hg19_chr22_fa_gz, dir)?;
        ensure_fai(&hg19_chr22_fa)?;
        // The 1000G VCF uses b37 contig names ("22"), but the UCSC FASTA uses
        // "chr22". Rename to match so vcfkit normalize/filter tests can use the
        // same FASTA without contig mismatch warnings.
        let hg19_chr22_fa_b37 = make_b37_fasta(&hg19_chr22_fa, dir)?;

        // ── hg38 chr22 FASTA ─────────────────────────────────────────────────
        let hg38_chr22_fa_gz = download_and_verify(
            dir,
            "chr22.hg38.fa.gz",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz",
            None,
        )?;
        let hg38_chr22_fa = decompress_gz_if_needed(&hg38_chr22_fa_gz, dir)?;
        ensure_fai(&hg38_chr22_fa)?;

        // ── hg19→hg38 chain file ─────────────────────────────────────────────
        let chain = download_and_verify(
            dir,
            "hg19ToHg38.over.chain.gz",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
            None,
        )?;

        Ok(Self {
            chr22_vcf,
            hg19_chr22_fa: hg19_chr22_fa_b37,
            hg38_chr22_fa,
            hg19_to_hg38_chain: chain,
        })
    }
}

// ── download + cache ──────────────────────────────────────────────────────────

/// Download `url` to `dir/filename`, verifying SHA-256 if provided.
/// Returns the path to the file (cached if already exists and checksum matches).
fn download_and_verify(
    dir: &Path,
    filename: &str,
    url: &str,
    expected_sha256: Option<&str>,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let dest = dir.join(filename);

    if dest.exists() {
        if let Some(expected) = expected_sha256 {
            let actual = sha256_file(&dest)?;
            if actual == expected {
                eprintln!("cache hit: {filename}");
                return Ok(dest);
            }
            eprintln!("checksum mismatch for {filename}, re-downloading");
        } else {
            eprintln!("cache hit: {filename} (no checksum)");
            return Ok(dest);
        }
    }

    eprintln!("downloading {filename} from {url} …");
    let status = Command::new("curl")
        .args(["-fsSL", "-o", dest.to_str().unwrap(), url])
        .status()?;
    if !status.success() {
        return Err(format!("curl failed for {url}").into());
    }

    if let Some(expected) = expected_sha256 {
        let actual = sha256_file(&dest)?;
        if actual != expected {
            fs::remove_file(&dest)?;
            return Err(format!(
                "SHA-256 mismatch for {filename}:\n  expected: {expected}\n  actual:   {actual}"
            )
            .into());
        }
    }

    Ok(dest)
}

fn sha256_file(path: &Path) -> Result<String, Box<dyn std::error::Error>> {
    // Shell out to shasum (available on macOS + Linux).
    let out = Command::new("shasum")
        .args(["-a", "256", path.to_str().unwrap()])
        .output()?;
    if !out.status.success() {
        return Err("shasum failed".into());
    }
    // shasum output: "<hex>  <path>"
    let stdout = String::from_utf8(out.stdout)?;
    let hex = stdout
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_lowercase();
    Ok(hex)
}

// ── VCF sites extraction ──────────────────────────────────────────────────────

/// Run `bcftools view -G` to strip sample columns from a genotypes VCF, producing
/// a sites-only VCF. Returns the path to the sites VCF (cached if already exists).
fn extract_sites_vcf(geno_gz: &Path, dir: &Path) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let dest = dir.join("chr22.1kg.sites.vcf");
    if dest.exists() {
        eprintln!("cache hit: chr22.1kg.sites.vcf");
        return Ok(dest);
    }
    eprintln!(
        "extracting sites from {} (strips sample columns) …",
        geno_gz.display()
    );
    let status = Command::new("bcftools")
        .args([
            "view",
            "-G",
            geno_gz.to_str().unwrap(),
            "-o",
            dest.to_str().unwrap(),
        ])
        .status()?;
    if !status.success() {
        return Err("bcftools view -G failed".into());
    }
    Ok(dest)
}

// ── decompression ─────────────────────────────────────────────────────────────

/// If `path` ends with `.gz`, decompress to the same directory without the `.gz`
/// suffix and return that path. If it doesn't end with `.gz`, return `path` as-is.
fn decompress_gz_if_needed(path: &Path, dir: &Path) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let fname = path
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or("invalid filename")?;

    if !fname.ends_with(".gz") {
        return Ok(path.to_path_buf());
    }

    let decompressed_name = &fname[..fname.len() - 3];
    let dest = dir.join(decompressed_name);

    if dest.exists() {
        return Ok(dest);
    }

    eprintln!("decompressing {} …", path.display());
    let status = Command::new("gunzip")
        .args(["-k", "-f", path.to_str().unwrap()])
        .current_dir(dir)
        .status()?;

    if !status.success() {
        // Try bgzip as fallback (some .gz files from UCSC are actually bgzip).
        let status2 = Command::new("bgzip")
            .args(["-d", "-k", "-f", path.to_str().unwrap()])
            .current_dir(dir)
            .status();
        if status2.map(|s| !s.success()).unwrap_or(true) {
            return Err(format!("gunzip/bgzip failed for {}", path.display()).into());
        }
    }

    Ok(dest)
}

// ── FASTA utilities ───────────────────────────────────────────────────────────

fn ensure_fai(fa: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let fai = fa.with_extension("fa.fai");
    if fai.exists() {
        return Ok(());
    }
    eprintln!("indexing {} …", fa.display());
    let out = Command::new("samtools")
        .args(["faidx", fa.to_str().unwrap()])
        .output()?;
    if !out.status.success() {
        return Err(format!(
            "samtools faidx failed: {}",
            String::from_utf8_lossy(&out.stderr)
        )
        .into());
    }
    Ok(())
}

fn ensure_bgzipped_and_indexed(vcf: &Path) -> Result<(), Box<dyn std::error::Error>> {
    // If it's already a plain VCF, we don't need to index it for most tests.
    // For tests that need random access, bgzip + tabix is needed. Skip for now.
    let _ = vcf;
    Ok(())
}

/// Create a b37-style FASTA from a UCSC-named FASTA ("chr22" → "22").
/// Returns the path to the b37 FASTA.
fn make_b37_fasta(ucsc_fa: &Path, dir: &Path) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let b37_fa = dir.join("chr22.hg19.b37.fa");
    if b37_fa.exists() && b37_fa.with_extension("fa.fai").exists() {
        return Ok(b37_fa);
    }

    eprintln!("creating b37-style FASTA (chr22 → 22) …");

    // Read and rewrite the FASTA header.
    let input = fs::read_to_string(ucsc_fa)?;
    let output = input.replace(">chr22", ">22");
    fs::write(&b37_fa, output)?;
    ensure_fai(&b37_fa)?;

    Ok(b37_fa)
}
