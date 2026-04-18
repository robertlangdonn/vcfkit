//! VCF/BCF I/O helpers: format detection and reader/writer construction.
//!
//! This module is the single entry point for opening VCF and BCF files
//! (or stdin/stdout) regardless of compression. All other vcfkit-core
//! modules should use these helpers rather than calling noodles directly.

use std::{
    io::{self, BufRead, BufReader, Read, Write},
    path::Path,
};

use noodles::{bcf, bgzf, vcf};

// ── OutputFormat ─────────────────────────────────────────────────────────────

/// The output container format for VCF-family data.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
pub enum OutputFormat {
    /// Plain-text VCF (optionally bgzf-compressed when the path ends in `.gz`/`.bgz`).
    #[default]
    Vcf,
    /// Binary BCF (always bgzf-compressed per the BCF spec).
    Bcf,
}

// ── Path-based helpers ────────────────────────────────────────────────────────

/// Detect the output format from a file-path extension.
///
/// * `.bcf` → [`OutputFormat::Bcf`]
/// * anything else → [`OutputFormat::Vcf`]
pub fn format_from_path(path: &Path) -> OutputFormat {
    match path.extension().and_then(|e| e.to_str()) {
        Some("bcf") => OutputFormat::Bcf,
        _ => OutputFormat::Vcf,
    }
}

/// Returns `true` when the path suggests the file is compressed.
///
/// Compressed extensions: `.gz`, `.bgz`, `.bcf` (BCF is always bgzf-compressed).
pub fn is_compressed_path(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|e| e.to_str()),
        Some("gz" | "bgz" | "bcf")
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
pub fn detect_format_from_magic(reader: &mut impl Read) -> io::Result<OutputFormat> {
    let mut magic = [0u8; 4];
    match reader.read_exact(&mut magic) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
            // File is shorter than 4 bytes — treat as VCF.
            return Ok(OutputFormat::Vcf);
        }
        Err(e) => return Err(e),
    }
    if magic == BCF_MAGIC {
        Ok(OutputFormat::Bcf)
    } else {
        Ok(OutputFormat::Vcf)
    }
}

// ── VcfReader enum ────────────────────────────────────────────────────────────

/// A format-agnostic VCF/BCF reader backed by a boxed [`BufRead`] or [`Read`].
///
/// Construct via [`open_vcf`].
pub enum VcfReader {
    /// Plain-text (possibly bgzf-compressed) VCF.
    Vcf(vcf::io::Reader<Box<dyn BufRead>>),
    /// Binary BCF.
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

/// A format-agnostic VCF/BCF writer backed by a boxed [`Write`].
///
/// Construct via [`create_vcf_writer`].
pub enum VcfWriter {
    /// Plain-text (possibly bgzf-compressed) VCF.
    Vcf(vcf::io::Writer<Box<dyn Write>>),
    /// Binary BCF.
    Bcf(bcf::io::Writer<Box<dyn Write>>),
}

impl VcfWriter {
    /// Write the VCF/BCF header.
    pub fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        match self {
            VcfWriter::Vcf(w) => w.write_header(header),
            VcfWriter::Bcf(w) => w.write_header(header),
        }
    }

    // TODO: add a `write_record` method once the normalize, liftover, and
    // filter modules are implemented and a shared record type is established.
}

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
            let fmt = format_from_path(p);
            match fmt {
                OutputFormat::Bcf => {
                    let reader = bcf::io::reader::Builder::default().build_from_path(p)?;
                    Ok(VcfReader::Bcf(reader))
                }
                OutputFormat::Vcf => {
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

            let fmt = if peeked.len() == 4 && peeked[..] == BCF_MAGIC {
                OutputFormat::Bcf
            } else {
                OutputFormat::Vcf
            };

            // Box the stdin as a `dyn Read` / `dyn BufRead`.
            match fmt {
                OutputFormat::Bcf => {
                    let inner: Box<dyn Read> =
                        Box::new(bgzf::io::Reader::new(buf_stdin));
                    Ok(VcfReader::Bcf(bcf::io::Reader::from(inner)))
                }
                OutputFormat::Vcf => {
                    let inner: Box<dyn BufRead> = Box::new(buf_stdin);
                    Ok(VcfReader::Vcf(vcf::io::Reader::new(inner)))
                }
            }
        }
    }
}

// ── create_vcf_writer ─────────────────────────────────────────────────────────

/// Open a VCF or BCF file for writing.
///
/// * `path = Some(p)` – write to the file at `p`, creating or truncating it.
/// * `path = None` – write to stdout.
///
/// The `format` argument controls whether VCF or BCF is written. When writing
/// VCF to a path ending in `.gz`/`.bgz`, bgzf compression is applied
/// automatically.
pub fn create_vcf_writer(
    path: Option<&Path>,
    format: OutputFormat,
) -> anyhow::Result<VcfWriter> {
    match path {
        Some(p) => match format {
            OutputFormat::Bcf => {
                let writer = bcf::io::writer::Builder::default().build_from_path(p)?;
                Ok(VcfWriter::Bcf(writer))
            }
            OutputFormat::Vcf => {
                let writer = vcf::io::writer::Builder::default().build_from_path(p)?;
                Ok(VcfWriter::Vcf(writer))
            }
        },
        None => {
            let stdout = io::stdout();
            match format {
                OutputFormat::Bcf => {
                    let inner: Box<dyn Write> =
                        Box::new(bgzf::io::Writer::new(stdout.lock()));
                    Ok(VcfWriter::Bcf(bcf::io::Writer::from(inner)))
                }
                OutputFormat::Vcf => {
                    let inner: Box<dyn Write> = Box::new(io::BufWriter::new(stdout.lock()));
                    Ok(VcfWriter::Vcf(vcf::io::Writer::new(inner)))
                }
            }
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    // ── format_from_path ──────────────────────────────────────────────────────

    #[test]
    fn format_from_path_bcf() {
        let p = PathBuf::from("variants.bcf");
        assert_eq!(format_from_path(&p), OutputFormat::Bcf);
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
    fn is_compressed_bcf() {
        assert!(is_compressed_path(Path::new("sample.bcf")));
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
        assert_eq!(
            detect_format_from_magic(&mut data).unwrap(),
            OutputFormat::Bcf
        );
    }

    #[test]
    fn magic_detects_vcf_from_hash() {
        // VCF headers start with "##fi…"
        let mut data = b"##fileformat=VCFv4.3\n".as_ref();
        assert_eq!(
            detect_format_from_magic(&mut data).unwrap(),
            OutputFormat::Vcf
        );
    }

    #[test]
    fn magic_empty_reader_is_vcf() {
        let mut data: &[u8] = b"";
        assert_eq!(
            detect_format_from_magic(&mut data).unwrap(),
            OutputFormat::Vcf
        );
    }
}
