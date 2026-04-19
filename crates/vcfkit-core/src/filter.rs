//! VCF record filtering by expression.
//!
//! The public entry points are:
//!
//! * [`FilterExpression::parse`] — parse a filter expression string into an AST.
//! * [`FilterExpression::evaluate`] — evaluate an expression against a record.
//! * [`filter`] — stream records through a parsed expression, writing those
//!   that match to the output.
//!
//! # Expression language (Phase 1)
//!
//! Fields:
//!
//! * `INFO/<key>` — access an INFO field by key.
//! * `FORMAT/<key>` — access the first sample's FORMAT field by key.
//! * `CHROM` — chromosome name (string).
//! * `POS` — position (integer).
//! * `QUAL` — quality score (float or `.`).
//! * `FILTER` — filter string (e.g. `PASS`).
//!
//! Operators:
//!
//! * Comparison: `<`, `<=`, `>`, `>=`, `==`, `!=`
//! * Logical: `&&`, `||`, `!`, parentheses
//! * Substring: `~` (contains), `!~` (not contains)
//! * Existence: bare field reference evaluates to true if the field is present
//!   and non-missing.
//!
//! Literals:
//!
//! * Numbers: `3.14`, `42`, `-1`
//! * Strings: `'PASS'`, `'chr17'`
//!
//! Type coercion follows the header schema — INFO/FORMAT fields declared with
//! `Type=Float` are parsed as f64 for numeric comparison; `Type=Integer` as
//! i64; anything else (including `FILTER`) as string. Comparisons against a
//! missing value (`.`) evaluate to false.

use std::{
    collections::HashMap,
    io::{BufRead, Write},
};

use noodles::vcf::{
    self,
    variant::{
        record_buf::{
            info::field::{value::Array as InfoArray, Value as InfoValue},
            samples::sample::{
                value::{Array as SampleArray, Genotype},
                Value as SampleValue,
            },
        },
        RecordBuf,
    },
};

use crate::error::VcfkitError;

// ── public API ───────────────────────────────────────────────────────────────

/// Options controlling the [`filter`] pipeline.
#[derive(Debug, Clone)]
pub struct FilterOptions {
    /// If true, invert the filter: keep records that *don't* match.
    pub invert: bool,
    /// Output format to write.
    pub output_format: crate::io::OutputFormat,
}

// Kept as a manual impl (rather than `#[derive(Default)]`) to match the
// spec exactly and to document the default values in one place.
#[allow(clippy::derivable_impls)]
impl Default for FilterOptions {
    fn default() -> Self {
        Self {
            invert: false,
            output_format: crate::io::OutputFormat::default(),
        }
    }
}

/// Statistics produced by a [`filter`] run.
#[derive(Debug, Default, Clone, Copy)]
pub struct FilterStats {
    /// Records read from the input.
    pub input_records: usize,
    /// Records written to the output.
    pub output_records: usize,
    /// Records that did not match the filter (i.e. were discarded).
    pub filtered_out: usize,
}

/// A parsed filter expression.
#[derive(Debug, Clone)]
pub struct FilterExpression {
    ast: Expr,
}

impl FilterExpression {
    /// Parse an expression string into an AST. Returns an error with a
    /// human-readable message (including a caret pointing at the failure) on
    /// syntax errors.
    pub fn parse(expr: &str) -> Result<Self, VcfkitError> {
        parser::parse(expr).map(|ast| FilterExpression { ast })
    }

    /// Evaluate the expression against a VCF record.
    ///
    /// Returns `Ok(true)` when the record passes the filter, `Ok(false)` when
    /// it does not, and `Err` only for evaluator bugs (well-formed expressions
    /// over well-formed records never return `Err`).
    pub fn evaluate(&self, record: &RecordBuf, header: &vcf::Header) -> Result<bool, VcfkitError> {
        eval_expr(&self.ast, record, header)
    }
}

/// Stream VCF records from `reader`, evaluate `expression` against each, and
/// write matching records (or non-matching, when `options.invert`) to
/// `writer`.
pub fn filter<R: BufRead, W: Write>(
    reader: R,
    writer: W,
    expression: FilterExpression,
    options: FilterOptions,
) -> Result<FilterStats, VcfkitError> {
    filter_with_progress(reader, writer, expression, options, |_| {})
}

/// Variant of [`filter`] that notifies `on_record` after each input record is
/// read. Used by the CLI to drive a progress bar.
pub fn filter_with_progress<R, W, F>(
    reader: R,
    mut writer: W,
    expression: FilterExpression,
    options: FilterOptions,
    mut on_record: F,
) -> Result<FilterStats, VcfkitError>
where
    R: BufRead,
    W: Write,
    F: FnMut(u64),
{
    // Phase 1: read the header with noodles so we have INFO type metadata.
    let mut vcf_reader = vcf::io::Reader::new(reader);
    let header = vcf_reader
        .read_header()
        .map_err(|e| VcfkitError::Other(format!("failed to read VCF header: {e}")))?;

    // Build a lightweight INFO-metadata map from the noodles header.
    let info_meta = build_info_meta(&header);

    // Write the header to output using noodles, then release the borrow.
    {
        let mut vcf_writer = vcf::io::Writer::new(&mut writer);
        vcf_writer
            .write_header(&header)
            .map_err(|e| VcfkitError::Other(format!("failed to write VCF header: {e}")))?;
    }

    // Phase 2: fast raw-line loop — no record allocation per line.
    // `vcf::io::Reader` wraps a `BufRead`; we recover it via `into_inner`.
    let mut raw_reader = vcf_reader.into_inner();
    let mut line = String::with_capacity(4096);
    let mut stats = FilterStats::default();

    loop {
        line.clear();
        let n = raw_reader
            .read_line(&mut line)
            .map_err(|e| VcfkitError::Other(format!("read error: {e}")))?;
        if n == 0 {
            break;
        }

        // Skip any residual header lines (should not appear after read_header,
        // but be safe).
        if line.starts_with('#') {
            continue;
        }

        stats.input_records += 1;
        on_record(stats.input_records as u64);

        let Some(rec) = FastRecord::parse(&line) else {
            continue;
        };

        let matches = eval_fast(&expression.ast, &rec, &info_meta)?;
        let keep = if options.invert { !matches } else { matches };

        if keep {
            writer
                .write_all(rec.raw.as_bytes())
                .map_err(|e| VcfkitError::Other(format!("write error: {e}")))?;
            stats.output_records += 1;
        } else {
            stats.filtered_out += 1;
        }
    }

    Ok(stats)
}

// ── fast record type ─────────────────────────────────────────────────────────

/// A zero-copy view of the first eight (mandatory) VCF columns, plus the
/// FORMAT and first-sample columns when present.  All slices borrow from the
/// raw line buffer — no heap allocation is needed to parse a record.
struct FastRecord<'a> {
    chrom: &'a str,
    pos_str: &'a str,
    _id: &'a str,
    _ref: &'a str,
    _alt: &'a str,
    qual_str: &'a str,
    filter_str: &'a str,
    info_str: &'a str,
    /// FORMAT column (column 8), or `""` when absent.
    format_str: &'a str,
    /// First sample column (column 9), or `""` when absent.
    sample_str: &'a str,
    /// The entire raw line including the trailing newline (used for
    /// pass-through writes so we never re-serialise).
    raw: &'a str,
}

impl<'a> FastRecord<'a> {
    /// Split a tab-delimited VCF data line into field slices.
    /// Returns `None` for lines that have fewer than 8 columns.
    fn parse(line: &'a str) -> Option<Self> {
        // We only need to split up to column 9 (FORMAT + first sample).
        // Using a manual byte scan avoids allocating a Vec<&str>.
        let mut iter = line.splitn(11, '\t');
        let chrom = iter.next()?;
        let pos_str = iter.next()?;
        let _id = iter.next()?;
        let _ref = iter.next()?;
        let _alt = iter.next()?;
        let qual_str = iter.next()?;
        let filter_str = iter.next()?;
        let info_str_raw = iter.next()?;
        // Trim any trailing newline from the INFO field when there are no
        // further columns (8-column VCF with no FORMAT/sample).
        let info_str = info_str_raw.trim_end_matches(&['\n', '\r'][..]);

        let format_str = iter
            .next()
            .unwrap_or("")
            .trim_end_matches(&['\n', '\r'][..]);
        let sample_str = iter
            .next()
            .unwrap_or("")
            .trim_end_matches(&['\n', '\r'][..]);

        Some(FastRecord {
            chrom,
            pos_str,
            _id,
            _ref,
            _alt,
            qual_str,
            filter_str,
            info_str,
            format_str,
            sample_str,
            raw: line,
        })
    }

    /// Look up a key in the INFO column using a linear scan.
    ///
    /// INFO format: `KEY1=VAL1;KEY2=VAL2;FLAG1;…`
    ///
    /// Returns the raw value string for `KEY=VAL` entries, or `"1"` for bare
    /// flag keys, or `None` when the key is absent.
    fn info_raw(&self, key: &str) -> Option<&'a str> {
        if self.info_str == "." || self.info_str.is_empty() {
            return None;
        }
        for entry in self.info_str.split(';') {
            match entry.split_once('=') {
                Some((k, v)) if k == key => {
                    // Value of "." means explicitly missing.
                    if v == "." {
                        return None;
                    }
                    return Some(v);
                }
                Some(_) => {}
                None => {
                    // Bare flag.
                    if entry == key {
                        return Some("1");
                    }
                }
            }
        }
        None
    }

    /// Look up a key in the FORMAT/first-sample columns.
    ///
    /// FORMAT: `GT:DP:GQ`   SAMPLE1: `0/1:40:55`
    ///
    /// Returns `None` when FORMAT or sample columns are absent or the key is
    /// not found.
    fn format_raw(&self, key: &str) -> Option<&'a str> {
        if self.format_str.is_empty() || self.sample_str.is_empty() {
            return None;
        }
        // Find the index of `key` in the FORMAT column, then extract the
        // corresponding value from the first sample column.
        let idx = self.format_str.split(':').position(|k| k == key)?;
        let val = self.sample_str.split(':').nth(idx)?;
        if val == "." {
            None
        } else {
            Some(val)
        }
    }
}

// ── INFO field metadata extracted from the noodles header ────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InfoNumberKind {
    /// `Number=1` — always a single value.
    One,
    /// `Number=A`, `R`, `G`, or `.` — potentially multiple comma-separated values.
    Many,
}

#[derive(Debug, Clone, Copy)]
struct InfoFieldMeta {
    ty: noodles::vcf::header::record::value::map::info::Type,
    number: InfoNumberKind,
}

fn build_info_meta(header: &vcf::Header) -> HashMap<String, InfoFieldMeta> {
    use noodles::vcf::header::record::value::map::info::Number;

    header
        .infos()
        .iter()
        .map(|(key, map)| {
            let number = match map.number() {
                Number::Count(1) => InfoNumberKind::One,
                _ => InfoNumberKind::Many,
            };
            (
                key.to_string(),
                InfoFieldMeta {
                    ty: map.ty(),
                    number,
                },
            )
        })
        .collect()
}

// ── VCF header schema extraction (used by --english CLI feature) ─────────────

/// Plain-data description of a single INFO or FORMAT field.
#[derive(Debug, Clone)]
pub struct VcfFieldDef {
    pub id: String,
    pub number: String,
    pub ty: String,
    pub description: String,
}

/// Schema extracted from a VCF header for use in LLM prompts.
#[derive(Debug, Clone, Default)]
pub struct VcfHeaderSchema {
    pub info_fields: Vec<VcfFieldDef>,
    pub format_fields: Vec<VcfFieldDef>,
    pub contigs: Vec<String>,
}

/// Read only the VCF header from `reader` and return a plain-data schema
/// suitable for building LLM prompts. Variant rows are not read.
pub fn extract_header_schema(
    reader: impl std::io::BufRead,
) -> Result<VcfHeaderSchema, VcfkitError> {
    use noodles::vcf::header::record::value::map::info::Number as InfoNumber;

    let mut vcf_reader = vcf::io::Reader::new(reader);
    let header = vcf_reader
        .read_header()
        .map_err(|e| VcfkitError::Other(format!("failed to read VCF header: {e}")))?;

    let info_fields = header
        .infos()
        .iter()
        .map(|(key, map)| VcfFieldDef {
            id: key.to_string(),
            number: match map.number() {
                InfoNumber::Count(n) => n.to_string(),
                InfoNumber::AlternateBases => "A".to_string(),
                InfoNumber::ReferenceAlternateBases => "R".to_string(),
                InfoNumber::Samples => "G".to_string(),
                InfoNumber::Unknown => ".".to_string(),
            },
            ty: format!("{:?}", map.ty()),
            description: map.description().to_string(),
        })
        .collect();

    let format_fields = header
        .formats()
        .iter()
        .map(|(key, map)| {
            use noodles::vcf::header::record::value::map::format::Number as FmtNumber;
            VcfFieldDef {
                id: key.to_string(),
                number: match map.number() {
                    FmtNumber::Count(n) => n.to_string(),
                    FmtNumber::AlternateBases => "A".to_string(),
                    FmtNumber::ReferenceAlternateBases => "R".to_string(),
                    FmtNumber::Samples => "G".to_string(),
                    FmtNumber::LocalAlternateBases => "LA".to_string(),
                    FmtNumber::LocalReferenceAlternateBases => "LR".to_string(),
                    FmtNumber::LocalSamples => "LG".to_string(),
                    FmtNumber::Ploidy => "P".to_string(),
                    FmtNumber::BaseModifications => "M".to_string(),
                    FmtNumber::Unknown => ".".to_string(),
                },
                ty: format!("{:?}", map.ty()),
                description: map.description().to_string(),
            }
        })
        .collect();

    let contigs = header.contigs().keys().map(|k| k.to_string()).collect();

    Ok(VcfHeaderSchema {
        info_fields,
        format_fields,
        contigs,
    })
}

// ── fast evaluator ───────────────────────────────────────────────────────────

fn eval_fast(
    e: &Expr,
    rec: &FastRecord<'_>,
    meta: &HashMap<String, InfoFieldMeta>,
) -> Result<bool, VcfkitError> {
    match e {
        Expr::And(a, b) => Ok(eval_fast(a, rec, meta)? && eval_fast(b, rec, meta)?),
        Expr::Or(a, b) => Ok(eval_fast(a, rec, meta)? || eval_fast(b, rec, meta)?),
        Expr::Not(x) => Ok(!eval_fast(x, rec, meta)?),
        Expr::Compare(lhs, op, rhs) => {
            let lv = fast_load_operand(lhs, rec, meta);
            let rv = fast_load_operand(rhs, rec, meta);
            compare(&lv, *op, &rv)
        }
        Expr::Exists(field) => Ok(!matches!(
            fast_load_field(field, rec, meta),
            Scalar::Missing
        )),
    }
}

fn fast_load_operand(
    o: &Operand,
    rec: &FastRecord<'_>,
    meta: &HashMap<String, InfoFieldMeta>,
) -> Scalar {
    match o {
        Operand::Field(f) => fast_load_field(f, rec, meta),
        Operand::Number(n) => Scalar::Float(*n),
        Operand::String(s) => Scalar::String(s.clone()),
    }
}

fn fast_load_field(
    f: &Field,
    rec: &FastRecord<'_>,
    meta: &HashMap<String, InfoFieldMeta>,
) -> Scalar {
    match f {
        Field::Chrom => Scalar::String(rec.chrom.to_string()),
        Field::Pos => match rec.pos_str.parse::<i64>() {
            Ok(p) => Scalar::Integer(p),
            Err(_) => Scalar::Missing,
        },
        Field::Qual => {
            let q = rec.qual_str.trim_end_matches(&['\n', '\r'][..]);
            if q == "." || q.is_empty() {
                Scalar::Missing
            } else {
                match q.parse::<f64>() {
                    Ok(v) => Scalar::Float(v),
                    Err(_) => Scalar::Missing,
                }
            }
        }
        Field::Filter => {
            let f = rec.filter_str;
            if f == "." || f.is_empty() {
                Scalar::Missing
            } else {
                Scalar::String(f.to_string())
            }
        }
        Field::Info(key) => fast_load_info(rec, meta, key),
        Field::Format(key) => fast_load_format(rec, key),
    }
}

fn fast_load_info(
    rec: &FastRecord<'_>,
    meta: &HashMap<String, InfoFieldMeta>,
    key: &str,
) -> Scalar {
    use noodles::vcf::header::record::value::map::info::Type as InfoType;

    let raw = match rec.info_raw(key) {
        Some(v) => v,
        None => return Scalar::Missing,
    };

    // "1" is returned for bare FLAG keys — treat as integer 1.
    if raw == "1" {
        // Confirm it's truly a flag by checking for bare key existence or
        // declared Flag type.  Since info_raw returns "1" only for bare keys,
        // we can just check: if the header says Flag, or the value is literally
        // "1" with no meta entry, return Integer(1).
        if let Some(m) = meta.get(key) {
            if m.ty == InfoType::Flag {
                return Scalar::Integer(1);
            }
        } else {
            // No header declaration but bare key — treat as flag.
            return Scalar::Integer(1);
        }
    }

    let field_meta = meta.get(key);

    // If the value contains commas and the field is declared as multi-valued,
    // parse into a FloatArray for "any element" semantics.
    let is_many = field_meta
        .map(|m| m.number == InfoNumberKind::Many)
        .unwrap_or(false);

    if is_many && raw.contains(',') {
        return fast_parse_array(raw, field_meta);
    }

    // Single value — parse according to declared type.
    if let Some(m) = field_meta {
        match m.ty {
            InfoType::Integer => {
                if let Ok(n) = raw.parse::<i64>() {
                    return Scalar::Integer(n);
                }
            }
            InfoType::Float => {
                if let Ok(x) = raw.parse::<f64>() {
                    return Scalar::Float(x);
                }
            }
            InfoType::Flag => return Scalar::Integer(1),
            _ => {}
        }
    } else {
        // No header declaration: try numeric coercion.
        if let Ok(n) = raw.parse::<i64>() {
            return Scalar::Integer(n);
        }
        if let Ok(x) = raw.parse::<f64>() {
            return Scalar::Float(x);
        }
    }

    Scalar::String(raw.to_string())
}

/// Parse a comma-separated multi-value INFO field into a `FloatArray` (for
/// numeric types) or a joined `String` (for string types).
///
/// Mirrors the semantics of `info_array_all` exactly.
fn fast_parse_array(raw: &str, field_meta: Option<&InfoFieldMeta>) -> Scalar {
    use noodles::vcf::header::record::value::map::info::Type as InfoType;

    let elements: Vec<&str> = raw.split(',').collect();

    if let Some(m) = field_meta {
        match m.ty {
            InfoType::Integer => {
                let vals: Option<Vec<f64>> = elements
                    .iter()
                    .filter(|s| **s != ".")
                    .map(|s| s.parse::<i64>().ok().map(|n| n as f64))
                    .collect();
                if let Some(v) = vals {
                    if v.is_empty() {
                        return Scalar::Missing;
                    }
                    return Scalar::FloatArray(v);
                }
            }
            InfoType::Float => {
                let vals: Option<Vec<f64>> = elements
                    .iter()
                    .filter(|s| **s != ".")
                    .map(|s| s.parse::<f64>().ok())
                    .collect();
                if let Some(v) = vals {
                    if v.is_empty() {
                        return Scalar::Missing;
                    }
                    return Scalar::FloatArray(v);
                }
            }
            _ => {}
        }
    } else {
        // No declared type — try to coerce all elements to float.
        let vals: Option<Vec<f64>> = elements
            .iter()
            .filter(|s| **s != ".")
            .map(|s| s.parse::<f64>().ok())
            .collect();
        if let Some(v) = vals {
            if !v.is_empty() {
                return Scalar::FloatArray(v);
            }
        }
    }

    // Fall back to joined string (for `~` substring matching on string arrays).
    let joined: Vec<&str> = elements.iter().filter(|s| **s != ".").copied().collect();
    if joined.is_empty() {
        Scalar::Missing
    } else {
        Scalar::String(joined.join(","))
    }
}

fn fast_load_format(rec: &FastRecord<'_>, key: &str) -> Scalar {
    match rec.format_raw(key) {
        None => Scalar::Missing,
        Some(val) => {
            // Attempt numeric coercion without header metadata.
            if let Ok(n) = val.parse::<i64>() {
                return Scalar::Integer(n);
            }
            if let Ok(x) = val.parse::<f64>() {
                return Scalar::Float(x);
            }
            Scalar::String(val.to_string())
        }
    }
}

// ── AST ──────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq)]
enum Expr {
    /// Logical AND.
    And(Box<Expr>, Box<Expr>),
    /// Logical OR.
    Or(Box<Expr>, Box<Expr>),
    /// Logical NOT.
    Not(Box<Expr>),
    /// Comparison `lhs <op> rhs`.
    Compare(Operand, CmpOp, Operand),
    /// Bare field reference: true iff the field is present and non-missing.
    Exists(Field),
}

#[derive(Debug, Clone, PartialEq)]
enum Operand {
    Field(Field),
    Number(f64),
    String(String),
}

#[derive(Debug, Clone, PartialEq)]
enum Field {
    /// `INFO/<key>`.
    Info(String),
    /// `FORMAT/<key>` (first sample).
    Format(String),
    Chrom,
    Pos,
    Qual,
    Filter,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CmpOp {
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
    /// String contains.
    Contains,
    /// String not-contains.
    NotContains,
}

// ── evaluator ────────────────────────────────────────────────────────────────

/// A scalar value extracted from a record, typed loosely so we can compare
/// across INFO/FORMAT type boundaries.
#[derive(Debug, Clone, PartialEq)]
enum Scalar {
    /// A missing value (`.`) or a field that is not present in the record.
    Missing,
    Integer(i64),
    Float(f64),
    String(String),
    /// Multiple numeric values from an INFO array field (Number=A/R/G/etc).
    /// Numeric comparisons use "any element satisfies" semantics (bcftools default).
    FloatArray(Vec<f64>),
}

fn eval_expr(e: &Expr, rec: &RecordBuf, hdr: &vcf::Header) -> Result<bool, VcfkitError> {
    match e {
        Expr::And(a, b) => Ok(eval_expr(a, rec, hdr)? && eval_expr(b, rec, hdr)?),
        Expr::Or(a, b) => Ok(eval_expr(a, rec, hdr)? || eval_expr(b, rec, hdr)?),
        Expr::Not(x) => Ok(!eval_expr(x, rec, hdr)?),
        Expr::Compare(lhs, op, rhs) => {
            let lv = load_operand(lhs, rec, hdr);
            let rv = load_operand(rhs, rec, hdr);
            compare(&lv, *op, &rv)
        }
        Expr::Exists(field) => Ok(!matches!(load_field(field, rec, hdr), Scalar::Missing)),
    }
}

fn load_operand(o: &Operand, rec: &RecordBuf, hdr: &vcf::Header) -> Scalar {
    match o {
        Operand::Field(f) => load_field(f, rec, hdr),
        Operand::Number(n) => Scalar::Float(*n),
        Operand::String(s) => Scalar::String(s.clone()),
    }
}

fn load_field(f: &Field, rec: &RecordBuf, hdr: &vcf::Header) -> Scalar {
    match f {
        Field::Chrom => Scalar::String(rec.reference_sequence_name().to_string()),
        Field::Pos => rec
            .variant_start()
            .map(|p| Scalar::Integer(p.get() as i64))
            .unwrap_or(Scalar::Missing),
        Field::Qual => rec
            .quality_score()
            .map(|q| Scalar::Float(q as f64))
            .unwrap_or(Scalar::Missing),
        Field::Filter => {
            let filters = rec.filters().as_ref();
            if filters.is_empty() {
                Scalar::Missing
            } else {
                // Join multiple filter strings with ';' as per the VCF spec.
                let joined: Vec<String> = filters.iter().cloned().collect();
                Scalar::String(joined.join(";"))
            }
        }
        Field::Info(key) => load_info(rec, hdr, key),
        Field::Format(key) => load_format(rec, hdr, key),
    }
}

fn load_info(rec: &RecordBuf, hdr: &vcf::Header, key: &str) -> Scalar {
    let value = match rec.info().get(key) {
        Some(Some(v)) => v,
        _ => return Scalar::Missing,
    };

    // Infer a scalar type from the header, but fall back to whatever the value
    // carries if the header doesn't declare the field.
    let declared_type = hdr.infos().get(key).map(|m| m.ty());

    info_value_to_scalar(value, declared_type)
}

fn info_value_to_scalar(
    value: &InfoValue,
    declared_type: Option<noodles::vcf::header::record::value::map::info::Type>,
) -> Scalar {
    use noodles::vcf::header::record::value::map::info::Type as InfoType;

    match value {
        InfoValue::Integer(n) => Scalar::Integer(*n as i64),
        InfoValue::Float(x) => Scalar::Float(*x as f64),
        InfoValue::Flag => Scalar::Integer(1),
        InfoValue::Character(c) => Scalar::String(c.to_string()),
        InfoValue::String(s) => {
            // Try to coerce to the declared numeric type first.
            if let Some(t) = declared_type {
                match t {
                    InfoType::Integer => {
                        if let Ok(n) = s.parse::<i64>() {
                            return Scalar::Integer(n);
                        }
                    }
                    InfoType::Float => {
                        if let Ok(x) = s.parse::<f64>() {
                            return Scalar::Float(x);
                        }
                    }
                    _ => {}
                }
            }
            Scalar::String(s.clone())
        }
        InfoValue::Array(arr) => info_array_all(arr, declared_type),
    }
}

/// For multi-value INFO arrays, returns a `FloatArray` so that numeric
/// comparisons use "any element satisfies" semantics (matching bcftools).
/// String arrays are joined for `~` matching.
fn info_array_all(
    arr: &InfoArray,
    declared_type: Option<noodles::vcf::header::record::value::map::info::Type>,
) -> Scalar {
    match arr {
        InfoArray::Integer(v) => {
            let vals: Vec<f64> = v
                .iter()
                .filter_map(|o| o.as_ref().map(|n| *n as f64))
                .collect();
            if vals.is_empty() {
                Scalar::Missing
            } else {
                Scalar::FloatArray(vals)
            }
        }
        InfoArray::Float(v) => {
            let vals: Vec<f64> = v
                .iter()
                .filter_map(|o| o.as_ref().map(|n| *n as f64))
                .collect();
            if vals.is_empty() {
                Scalar::Missing
            } else {
                Scalar::FloatArray(vals)
            }
        }
        InfoArray::Character(v) => v
            .first()
            .and_then(|o| o.as_ref())
            .map(|c| Scalar::String(c.to_string()))
            .unwrap_or(Scalar::Missing),
        InfoArray::String(v) => {
            let elements: Vec<&str> = v
                .iter()
                .filter_map(|o| o.as_ref().map(|s| s.as_str()))
                .collect();
            if elements.is_empty() {
                return Scalar::Missing;
            }
            // Try numeric coercion — if all elements parse as floats, use FloatArray.
            if let Some(t) = declared_type {
                use noodles::vcf::header::record::value::map::info::Type as InfoType;
                match t {
                    InfoType::Integer => {
                        let nums: Option<Vec<f64>> = elements
                            .iter()
                            .map(|s| s.parse::<i64>().ok().map(|n| n as f64))
                            .collect();
                        if let Some(vals) = nums {
                            return Scalar::FloatArray(vals);
                        }
                    }
                    InfoType::Float => {
                        let nums: Option<Vec<f64>> =
                            elements.iter().map(|s| s.parse::<f64>().ok()).collect();
                        if let Some(vals) = nums {
                            return Scalar::FloatArray(vals);
                        }
                    }
                    _ => {}
                }
            }
            // Return all elements joined so `~` can match any of them.
            Scalar::String(elements.join(","))
        }
    }
}

fn load_format(rec: &RecordBuf, _hdr: &vcf::Header, key: &str) -> Scalar {
    // Look up the key's column in the FORMAT keys, then pick the first
    // sample's value in that column.
    let samples = rec.samples();
    let series = match samples.select(key) {
        Some(s) => s,
        None => return Scalar::Missing,
    };
    let first = match series.get(0) {
        Some(Some(v)) => v,
        _ => return Scalar::Missing,
    };
    sample_value_to_scalar(first)
}

fn sample_value_to_scalar(v: &SampleValue) -> Scalar {
    match v {
        SampleValue::Integer(n) => Scalar::Integer(*n as i64),
        SampleValue::Float(x) => Scalar::Float(*x as f64),
        SampleValue::Character(c) => Scalar::String(c.to_string()),
        SampleValue::String(s) => Scalar::String(s.clone()),
        SampleValue::Genotype(g) => Scalar::String(format_genotype(g)),
        SampleValue::Array(arr) => match arr {
            SampleArray::Integer(v) => v
                .first()
                .and_then(|o| o.as_ref())
                .map(|n| Scalar::Integer(*n as i64))
                .unwrap_or(Scalar::Missing),
            SampleArray::Float(v) => v
                .first()
                .and_then(|o| o.as_ref())
                .map(|n| Scalar::Float(*n as f64))
                .unwrap_or(Scalar::Missing),
            SampleArray::Character(v) => v
                .first()
                .and_then(|o| o.as_ref())
                .map(|c| Scalar::String(c.to_string()))
                .unwrap_or(Scalar::Missing),
            SampleArray::String(v) => v
                .first()
                .and_then(|o| o.as_ref())
                .map(|s| Scalar::String(s.clone()))
                .unwrap_or(Scalar::Missing),
        },
    }
}

/// Render a genotype back into its VCF textual form (e.g. `0/1`, `1|0`,
/// `./.`). The first allele's phasing is emitted as a leading separator only
/// when phased, per the VCF spec.
fn format_genotype(g: &Genotype) -> String {
    use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;

    let alleles = g.as_ref();
    let mut out = String::new();
    for (i, a) in alleles.iter().enumerate() {
        if i > 0 {
            match a.phasing() {
                Phasing::Phased => out.push('|'),
                Phasing::Unphased => out.push('/'),
            }
        } else if matches!(a.phasing(), Phasing::Phased) {
            out.push('|');
        }
        match a.position() {
            Some(p) => out.push_str(&p.to_string()),
            None => out.push('.'),
        }
    }
    out
}

/// Compare two scalars. A missing operand always evaluates to `false`.
/// For `FloatArray` (multi-value INFO fields), returns true if ANY element
/// satisfies the condition — matching bcftools default semantics.
fn compare(lhs: &Scalar, op: CmpOp, rhs: &Scalar) -> Result<bool, VcfkitError> {
    // Any comparison with a missing operand is false.
    if matches!(lhs, Scalar::Missing) || matches!(rhs, Scalar::Missing) {
        return Ok(false);
    }

    // FloatArray: any-element semantics for numeric ops, joined string for ~.
    if let Scalar::FloatArray(vals) = lhs {
        return match op {
            CmpOp::Contains | CmpOp::NotContains => {
                let joined = vals
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                compare(&Scalar::String(joined), op, rhs)
            }
            _ => {
                if let Some(b) = to_f64(rhs) {
                    Ok(vals.iter().any(|a| cmp_f64(*a, op, b)))
                } else {
                    Ok(false)
                }
            }
        };
    }

    match op {
        CmpOp::Contains => Ok(match (lhs, rhs) {
            (Scalar::String(a), Scalar::String(b)) => a.contains(b.as_str()),
            _ => false,
        }),
        CmpOp::NotContains => Ok(match (lhs, rhs) {
            (Scalar::String(a), Scalar::String(b)) => !a.contains(b.as_str()),
            _ => false,
        }),
        _ => {
            // Numeric vs numeric (or either side numeric) -> numeric compare.
            // String vs string -> string compare.
            if let (Some(a), Some(b)) = (to_f64(lhs), to_f64(rhs)) {
                Ok(cmp_f64(a, op, b))
            } else {
                let sa = scalar_to_string(lhs);
                let sb = scalar_to_string(rhs);
                Ok(match op {
                    CmpOp::Lt => sa < sb,
                    CmpOp::Le => sa <= sb,
                    CmpOp::Gt => sa > sb,
                    CmpOp::Ge => sa >= sb,
                    CmpOp::Eq => sa == sb,
                    CmpOp::Ne => sa != sb,
                    CmpOp::Contains | CmpOp::NotContains => {
                        return Err(VcfkitError::Other(
                            "internal: contains op in numeric comparison".to_string(),
                        ));
                    }
                })
            }
        }
    }
}

fn cmp_f64(a: f64, op: CmpOp, b: f64) -> bool {
    match op {
        CmpOp::Lt => a < b,
        CmpOp::Le => a <= b,
        CmpOp::Gt => a > b,
        CmpOp::Ge => a >= b,
        CmpOp::Eq => a == b,
        CmpOp::Ne => a != b,
        CmpOp::Contains | CmpOp::NotContains => false,
    }
}

fn to_f64(s: &Scalar) -> Option<f64> {
    match s {
        Scalar::Integer(n) => Some(*n as f64),
        Scalar::Float(x) => Some(*x),
        Scalar::String(s) => s.parse().ok(),
        Scalar::Missing | Scalar::FloatArray(_) => None,
    }
}

fn scalar_to_string(s: &Scalar) -> String {
    match s {
        Scalar::Integer(n) => n.to_string(),
        Scalar::Float(x) => x.to_string(),
        Scalar::String(s) => s.clone(),
        Scalar::FloatArray(v) => v
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(","),
        Scalar::Missing => String::new(),
    }
}

// ── parser (nom) ─────────────────────────────────────────────────────────────

mod parser {
    use super::{CmpOp, Expr, Field, Operand};
    use crate::error::VcfkitError;
    use nom::{
        branch::alt,
        bytes::complete::{tag, take_while1},
        character::complete::{char, digit1, multispace0},
        combinator::{map, opt, recognize},
        multi::many0,
        sequence::{delimited, pair, preceded, terminated, tuple},
        IResult,
    };

    /// Public parse entry point. On failure, builds an error message with a
    /// caret pointing at the offending position.
    pub(super) fn parse(input: &str) -> Result<Expr, VcfkitError> {
        let trimmed = input.trim();
        match terminated(or_expr, multispace0)(trimmed) {
            Ok(("", ast)) => Ok(ast),
            Ok((rest, _)) => Err(build_error(
                input,
                trimmed,
                rest,
                "unexpected trailing input",
            )),
            Err(e) => Err(nom_error(input, trimmed, e)),
        }
    }

    fn build_error(original: &str, trimmed: &str, rest: &str, msg: &str) -> VcfkitError {
        // Find offset of `rest` inside the trimmed input, then map to the
        // original (accounting for leading whitespace that was trimmed).
        let offset_in_trimmed = trimmed.len().saturating_sub(rest.len());
        let leading_ws = original.len() - original.trim_start().len();
        let col = leading_ws + offset_in_trimmed;
        let caret: String = (0..col).map(|_| ' ').collect::<String>() + "^";
        VcfkitError::Other(format!(
            "invalid filter expression: {msg}\n  {original}\n  {caret}"
        ))
    }

    fn nom_error(
        original: &str,
        trimmed: &str,
        err: nom::Err<nom::error::Error<&str>>,
    ) -> VcfkitError {
        match err {
            nom::Err::Incomplete(_) => VcfkitError::Other(format!(
                "invalid filter expression: incomplete input\n  {original}"
            )),
            nom::Err::Error(e) | nom::Err::Failure(e) => {
                build_error(original, trimmed, e.input, "parse error")
            }
        }
    }

    // ── grammar ──────────────────────────────────────────────────────────────
    //
    //   or_expr    := and_expr ( "||" and_expr )*
    //   and_expr   := unary ( "&&" unary )*
    //   unary      := "!" unary | atom
    //   atom       := "(" or_expr ")" | comparison | existence
    //   comparison := operand cmp_op operand
    //   existence  := field  (bare field reference)
    //   operand    := field | number | string
    //   field      := "INFO/" ident | "FORMAT/" ident | "CHROM" | "POS" | "QUAL" | "FILTER"
    //   number     := -?[0-9]+ ( "." [0-9]+ )?
    //   string     := "'" ... "'"
    //   cmp_op     := "<=" | ">=" | "==" | "!=" | "<" | ">" | "~" | "!~"

    fn or_expr(i: &str) -> IResult<&str, Expr> {
        let (i, first) = and_expr(i)?;
        let (i, rest) = many0(preceded(ws(tag("||")), and_expr))(i)?;
        Ok((
            i,
            rest.into_iter()
                .fold(first, |acc, x| Expr::Or(Box::new(acc), Box::new(x))),
        ))
    }

    fn and_expr(i: &str) -> IResult<&str, Expr> {
        let (i, first) = unary(i)?;
        let (i, rest) = many0(preceded(ws(tag("&&")), unary))(i)?;
        Ok((
            i,
            rest.into_iter()
                .fold(first, |acc, x| Expr::And(Box::new(acc), Box::new(x))),
        ))
    }

    fn unary(i: &str) -> IResult<&str, Expr> {
        alt((
            map(preceded(ws(char('!')), unary), |e| Expr::Not(Box::new(e))),
            atom,
        ))(i)
    }

    fn atom(i: &str) -> IResult<&str, Expr> {
        alt((
            delimited(ws(char('(')), or_expr, ws(char(')'))),
            comparison_or_existence,
        ))(i)
    }

    /// Try a comparison first; if that fails, fall back to a bare field
    /// reference as an existence check.
    fn comparison_or_existence(i: &str) -> IResult<&str, Expr> {
        alt((comparison, map(ws(field), Expr::Exists)))(i)
    }

    fn comparison(i: &str) -> IResult<&str, Expr> {
        map(
            tuple((ws(operand), ws(cmp_op), ws(operand))),
            |(a, op, b)| Expr::Compare(a, op, b),
        )(i)
    }

    fn cmp_op(i: &str) -> IResult<&str, CmpOp> {
        // Order matters: longer operators first.
        alt((
            map(tag("<="), |_| CmpOp::Le),
            map(tag(">="), |_| CmpOp::Ge),
            map(tag("=="), |_| CmpOp::Eq),
            map(tag("!="), |_| CmpOp::Ne),
            map(tag("!~"), |_| CmpOp::NotContains),
            map(tag("~"), |_| CmpOp::Contains),
            map(tag("<"), |_| CmpOp::Lt),
            map(tag(">"), |_| CmpOp::Gt),
        ))(i)
    }

    fn operand(i: &str) -> IResult<&str, Operand> {
        alt((
            map(string_lit, Operand::String),
            map(number, Operand::Number),
            map(field, Operand::Field),
        ))(i)
    }

    fn field(i: &str) -> IResult<&str, Field> {
        alt((
            map(preceded(tag("INFO/"), ident), |k: &str| {
                Field::Info(k.to_string())
            }),
            map(preceded(tag("FORMAT/"), ident), |k: &str| {
                Field::Format(k.to_string())
            }),
            map(tag("CHROM"), |_| Field::Chrom),
            map(tag("POS"), |_| Field::Pos),
            map(tag("QUAL"), |_| Field::Qual),
            map(tag("FILTER"), |_| Field::Filter),
        ))(i)
    }

    fn ident(i: &str) -> IResult<&str, &str> {
        recognize(pair(
            take_while1(|c: char| c.is_ascii_alphabetic() || c == '_'),
            take_while1_or_empty(|c: char| c.is_ascii_alphanumeric() || c == '_'),
        ))(i)
    }

    /// `take_while` that is allowed to match zero characters.
    fn take_while1_or_empty<F>(f: F) -> impl Fn(&str) -> IResult<&str, &str>
    where
        F: Fn(char) -> bool,
    {
        move |i: &str| {
            let end = i
                .char_indices()
                .find(|(_, c)| !f(*c))
                .map(|(idx, _)| idx)
                .unwrap_or(i.len());
            Ok((&i[end..], &i[..end]))
        }
    }

    fn number(i: &str) -> IResult<&str, f64> {
        let (i, raw) = recognize(tuple((
            opt(char('-')),
            digit1,
            opt(preceded(char('.'), digit1)),
        )))(i)?;
        let n = raw.parse::<f64>().map_err(|_| {
            nom::Err::Error(nom::error::Error::new(i, nom::error::ErrorKind::Float))
        })?;
        Ok((i, n))
    }

    fn string_lit(i: &str) -> IResult<&str, String> {
        let (i, _) = char('\'')(i)?;
        let (i, body) = take_while1_or_empty(|c: char| c != '\'')(i)?;
        let (i, _) = char('\'')(i)?;
        Ok((i, body.to_string()))
    }

    fn ws<'a, F, O>(inner: F) -> impl FnMut(&'a str) -> IResult<&'a str, O>
    where
        F: FnMut(&'a str) -> IResult<&'a str, O>,
    {
        delimited(multispace0, inner, multispace0)
    }
}

// ── tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_comparison() {
        let e = FilterExpression::parse("INFO/AF < 0.05").unwrap();
        match e.ast {
            Expr::Compare(Operand::Field(Field::Info(k)), CmpOp::Lt, Operand::Number(n)) => {
                assert_eq!(k, "AF");
                assert!((n - 0.05).abs() < 1e-9);
            }
            other => panic!("unexpected AST: {other:?}"),
        }
    }

    #[test]
    fn parse_and_or_precedence() {
        // && binds tighter than ||
        let e = FilterExpression::parse("CHROM == 'chr1' || CHROM == 'chr2' && POS > 10").unwrap();
        match e.ast {
            Expr::Or(_, rhs) => match *rhs {
                Expr::And(_, _) => {}
                other => panic!("rhs should be And, got {other:?}"),
            },
            other => panic!("top should be Or, got {other:?}"),
        }
    }

    #[test]
    fn parse_not_operator() {
        let e = FilterExpression::parse("!(FILTER == 'PASS')").unwrap();
        match e.ast {
            Expr::Not(_) => {}
            other => panic!("expected Not, got {other:?}"),
        }
    }

    #[test]
    fn parse_substring_operator() {
        let e = FilterExpression::parse("INFO/CSQ ~ 'missense'").unwrap();
        match e.ast {
            Expr::Compare(Operand::Field(Field::Info(k)), CmpOp::Contains, Operand::String(s)) => {
                assert_eq!(k, "CSQ");
                assert_eq!(s, "missense");
            }
            other => panic!("unexpected AST: {other:?}"),
        }
    }

    #[test]
    fn parse_existence_bare_field() {
        let e = FilterExpression::parse("INFO/AF").unwrap();
        match e.ast {
            Expr::Exists(Field::Info(k)) => assert_eq!(k, "AF"),
            other => panic!("unexpected AST: {other:?}"),
        }
    }

    #[test]
    fn parse_error_has_caret() {
        let err = FilterExpression::parse("INFO/AF <> 3").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("invalid filter expression"), "got: {msg}");
        assert!(msg.contains('^'), "expected a caret, got: {msg}");
    }

    #[test]
    fn parse_parentheses_group() {
        let e = FilterExpression::parse("(POS > 100 && POS < 200) || CHROM == 'chrX'").unwrap();
        match e.ast {
            Expr::Or(lhs, _) => match *lhs {
                Expr::And(_, _) => {}
                other => panic!("lhs should be And, got {other:?}"),
            },
            other => panic!("top should be Or, got {other:?}"),
        }
    }

    #[test]
    fn parse_negative_number() {
        let e = FilterExpression::parse("POS > -1").unwrap();
        match e.ast {
            Expr::Compare(_, CmpOp::Gt, Operand::Number(n)) => {
                assert!((n - -1.0).abs() < 1e-9);
            }
            other => panic!("unexpected AST: {other:?}"),
        }
    }
}
