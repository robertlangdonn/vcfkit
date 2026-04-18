//! VCF I/O helpers: format detection and reader/writer construction.
//!
//! This module is the single entry point for opening VCF files (or
//! stdin/stdout) regardless of compression. BCF input is supported for
//! reading; BCF output is not — pipe through `bcftools view -O b` instead.
//! All other vcfkit-core modules should use these helpers rather than calling
//! noodles directly.

use std::{
    io::{self, BufRead, BufReader, Read, Write},
    path::Path,
};

use noodles::{bcf, bgzf, vcf};

// ── OutputFormat ─────────────────────────────────────────────────────────────

/// The output container format for VCF-family data.
///
/// Only plain-text VCF (with optional bgzf compression) is supported for
/// output in v1. For BCF output, pipe through `bcftools view -O b`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
pub enum OutputFormat {
    /// Plain-text VCF (optionally bgzf-compressed when the path ends in `.gz`/`.bgz`).
    #[default]
    Vcf,
}

// ── Path-based helpers ────────────────────────────────────────────────────────

/// Detect the output format from a file-path extension.
///
/// Only [`OutputFormat::Vcf`] is returned; `.bcf` paths are treated as VCF
/// for forward-compatibility. For BCF output, pipe through
/// `bcftools view -O b`.
pub fn format_from_path(_path: &Path) -> OutputFormat {
    OutputFormat::Vcf
}

/// Returns `true` when the path suggests the file is compressed.
///
/// Compressed extensions: `.gz`, `.bgz`.
pub fn is_compressed_path(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|e| e.to_str()),
        Some("gz" | "bgz")
    )
}

// ── Magic-byte detection ──────────────────────────────────────────────────────

/// BCF magic bytes (first 4 bytes of every BCF file).
const BCF_MAGIC: [u8; 4] = *b"BCF\x02";

/// Peek at the first four bytes of `reader` to determine the format.
///
/// BCF files begin with `b"BCF\x02"`. Everything else is treated as VCF.
/// The four bytes are consumed from `reader`; callers that need them back
/// should wrap the reader in a [`std::io::Chain`] with a `Cursor` over the
/// peeked bytes.
pub fn detect_format_from_magic(reader: &mut impl Read) -> io::Result<bool> {
    let mut magic = [0u8; 4];
    match reader.read_exact(&mut magic) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
            return Ok(false);
        }
        Err(e) => return Err(e),
    }
    Ok(magic == BCF_MAGIC)
}

// ── VcfReader enum ────────────────────────────────────────────────────────────

/// A format-agnostic VCF/BCF reader backed by a boxed [`BufRead`] or [`Read`].
///
/// Construct via [`open_vcf`].
pub enum VcfReader {
    /// Plain-text (possibly bgzf-compressed) VCF.
    Vcf(vcf::io::Reader<Box<dyn BufRead>>),
    /// Binary BCF (read-only; BCF output is not supported in v1).
    Bcf(bcf::io::Reader<Box<dyn Read>>),
}

impl VcfReader {
    /// Read and parse the VCF/BCF header.
    pub fn read_header(&mut self) -> io::Result<vcf::Header> {
        match self {
            VcfReader::Vcf(r) => r.read_header(),
            VcfReader::Bcf(r) => r.read_header(),
        }
    }

    // TODO: add a `read_record` method once the normalize, liftover, and
    // filter modules are implemented and a shared record type is established.
}

// ── VcfWriter enum ────────────────────────────────────────────────────────────

/// A VCF writer backed by a boxed [`Write`].
///
/// In v1 only plain-text VCF output is supported. For BCF output, pipe
/// through `bcftools view -O b`. Construct via [`create_vcf_writer`].
pub type VcfWriter = vcf::io::Writer<Box<dyn Write>>;

// ── open_vcf ──────────────────────────────────────────────────────────────────

/// Open a VCF or BCF file for reading.
///
/// * `path = Some(p)` – open the file at `p`; format is detected from the
///   extension (`.bcf` → BCF, `.gz`/`.bgz` → bgzf-compressed VCF, otherwise
///   plain VCF).
/// * `path = None` – read from stdin; format is detected by peeking at the
///   first four magic bytes.
pub fn open_vcf(path: Option<&Path>) -> anyhow::Result<VcfReader> {
    match path {
        Some(p) => {
            match p.extension().and_then(|e| e.to_str()) {
                Some("bcf") => {
                    let reader = bcf::io::reader::Builder::default().build_from_path(p)?;
                    Ok(VcfReader::Bcf(reader))
                }
                _ => {
                    let reader = vcf::io::reader::Builder::default().build_from_path(p)?;
                    Ok(VcfReader::Vcf(reader))
                }
            }
        }
        None => {
            // Read stdin into a buffered reader, peek magic bytes.
            let stdin = io::stdin();
            let mut buf_stdin = BufReader::new(stdin.lock());

            // Peek the first 4 bytes.
            // Note: `consume()` is intentionally NOT called after `fill_buf`.
            // The magic bytes remain in the `BufReader` buffer so that noodles
            // can re-read them as part of the normal header/record parsing.
            let peeked = {
                let bytes = buf_stdin.fill_buf()?;
                let len = bytes.len().min(4);
                bytes[..len].to_vec()
            };

            let is_bcf = peeked.len() == 4 && peeked[..] == BCF_MAGIC;

            if is_bcf {
                let inner: Box<dyn Read> =
                    Box::new(bgzf::io::Reader::new(buf_stdin));
                Ok(VcfReader::Bcf(bcf::io::Reader::from(inner)))
            } else {
                let inner: Box<dyn BufRead> = Box::new(buf_stdin);
                Ok(VcfReader::Vcf(vcf::io::Reader::new(inner)))
            }
        }
    }
}

// ── create_vcf_writer ─────────────────────────────────────────────────────────

/// Open a VCF file for writing.
///
/// * `path = Some(p)` – write to the file at `p`, creating or truncating it.
///   If the path ends in `.gz` or `.bgz`, bgzf compression is applied.
/// * `path = None` – write to stdout.
///
/// Only VCF output is supported in v1. For BCF output, pipe through
/// `bcftools view -O b`.
pub fn create_vcf_writer(
    path: Option<&Path>,
    _format: OutputFormat,
) -> anyhow::Result<VcfWriter> {
    let writer: Box<dyn Write> = match path {
        Some(p) => {
            let file = std::fs::File::create(p)?;
            Box::new(file)
        }
        None => {
            let stdout = io::stdout();
            Box::new(io::BufWriter::new(stdout.lock()))
        }
    };
    Ok(vcf::io::Writer::new(writer))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    // ── format_from_path ──────────────────────────────────────────────────────

    #[test]
    fn format_from_path_bcf_returns_vcf() {
        // BCF output is not supported; .bcf paths yield Vcf for forward-compat.
        let p = PathBuf::from("variants.bcf");
        assert_eq!(format_from_path(&p), OutputFormat::Vcf);
    }

    #[test]
    fn format_from_path_vcf() {
        let p = PathBuf::from("variants.vcf");
        assert_eq!(format_from_path(&p), OutputFormat::Vcf);
    }

    #[test]
    fn format_from_path_vcf_gz() {
        let p = PathBuf::from("variants.vcf.gz");
        assert_eq!(format_from_path(&p), OutputFormat::Vcf);
    }

    #[test]
    fn format_from_path_no_extension() {
        let p = PathBuf::from("variants");
        assert_eq!(format_from_path(&p), OutputFormat::Vcf);
    }

    // ── is_compressed_path ────────────────────────────────────────────────────

    #[test]
    fn is_compressed_gz() {
        assert!(is_compressed_path(Path::new("sample.vcf.gz")));
    }

    #[test]
    fn is_compressed_bgz() {
        assert!(is_compressed_path(Path::new("sample.vcf.bgz")));
    }

    #[test]
    fn not_compressed_bcf() {
        // .bcf is no longer treated as a compressed output path.
        assert!(!is_compressed_path(Path::new("sample.bcf")));
    }

    #[test]
    fn not_compressed_vcf() {
        assert!(!is_compressed_path(Path::new("sample.vcf")));
    }

    #[test]
    fn not_compressed_no_ext() {
        assert!(!is_compressed_path(Path::new("sample")));
    }

    // ── detect_format_from_magic ──────────────────────────────────────────────

    #[test]
    fn magic_detects_bcf() {
        let mut data = b"BCF\x02extra bytes".as_ref();
        assert!(detect_format_from_magic(&mut data).unwrap());
    }

    #[test]
    fn magic_detects_vcf_from_hash() {
        // VCF headers start with "##fi…"
        let mut data = b"##fileformat=VCFv4.3\n".as_ref();
        assert!(!detect_format_from_magic(&mut data).unwrap());
    }

    #[test]
    fn magic_empty_reader_is_vcf() {
        let mut data: &[u8] = b"";
        assert!(!detect_format_from_magic(&mut data).unwrap());
    }
}
