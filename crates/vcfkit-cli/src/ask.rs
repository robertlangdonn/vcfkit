//! Natural-language filter translation via Anthropic's Claude API.
//!
//! When a user runs `vcfkit filter --ask "<query>"`, this module translates
//! the query into a deterministic filter expression. The expression goes
//! through the existing parser before running, so the LLM cannot cause
//! arbitrary behavior.
//!
//! The LLM sees only the VCF header schema (INFO/FORMAT definitions,
//! contig list) and the query text. Variant data never leaves the user's
//! machine.

use std::time::Duration;

use serde::{Deserialize, Serialize};

pub const DEFAULT_MODEL: &str = "claude-sonnet-4-6";
const DEFAULT_TIMEOUT_SECS: u64 = 30;
const MAX_TOKENS: u32 = 1024;
const API_URL: &str = "https://api.anthropic.com/v1/messages";
const ANTHROPIC_VERSION: &str = "2023-06-01";

// ── public types ─────────────────────────────────────────────────────────────

pub struct Translation {
    pub expression: String,
    pub reasoning: String,
    pub confidence: f64,
    pub caveats: Vec<String>,
    pub model: String,
}

#[derive(Debug, thiserror::Error)]
pub enum AskError {
    #[error(
        "ANTHROPIC_API_KEY is not set.\n\
         Get a key at https://console.anthropic.com and then:\n\
           export ANTHROPIC_API_KEY=sk-ant-..."
    )]
    MissingApiKey,

    #[error("Anthropic API returned HTTP {status}:\n{body}")]
    ApiError { status: u16, body: String },

    #[error("Network error contacting Anthropic API: {0}")]
    NetworkError(#[from] reqwest::Error),

    #[error(
        "Request timed out after {0}s.\n\
         Override with: export VCFKIT_LLM_TIMEOUT_SECS=60"
    )]
    Timeout(u64),

    #[error(
        "Translation produced an expression that failed to parse:\n\
         \n  {expression}\n\
         \nParser error: {error}\n\
         \nTry rephrasing your query, or use -e to write the expression directly."
    )]
    InvalidTranslation { expression: String, error: String },

    #[error("Unexpected API response format: {0}")]
    ResponseParse(String),
}

pub struct HeaderSchema {
    pub info_fields: Vec<FieldDef>,
    pub format_fields: Vec<FieldDef>,
    pub contigs: Vec<String>,
}

pub struct FieldDef {
    pub id: String,
    pub number: String,
    pub ty: String,
    pub description: String,
}

// ── Anthropic API types ───────────────────────────────────────────────────────

#[derive(Serialize)]
struct AnthropicRequest<'a> {
    model: &'a str,
    max_tokens: u32,
    temperature: f32,
    system: &'a str,
    messages: Vec<Message<'a>>,
}

#[derive(Serialize)]
struct Message<'a> {
    role: &'a str,
    content: &'a str,
}

#[derive(Deserialize)]
struct AnthropicResponse {
    content: Vec<ContentBlock>,
    model: String,
}

#[derive(Deserialize)]
struct ContentBlock {
    #[serde(rename = "type")]
    block_type: String,
    text: Option<String>,
}

#[derive(Deserialize)]
struct TranslationPayload {
    expression: String,
    reasoning: String,
    confidence: f64,
    #[serde(default)]
    caveats: Vec<String>,
}

struct ApiCallResult {
    content: String,
    model: String,
}

// ── public entry point ────────────────────────────────────────────────────────

/// Translate a natural-language query into a vcfkit filter expression.
///
/// Checks `VCFKIT_MOCK_TRANSLATION` first — if set, that JSON is used directly
/// (bypasses the API entirely; useful for testing and offline use).
///
/// Otherwise calls the Anthropic API. If the returned expression fails to
/// parse, retries once with the parser error appended to the prompt.
pub async fn translate(query: &str, schema: &HeaderSchema) -> Result<Translation, AskError> {
    // Mock path — set VCFKIT_MOCK_TRANSLATION to a JSON string to bypass the
    // API. Useful for integration tests and offline smoke-testing.
    if let Ok(mock_json) = std::env::var("VCFKIT_MOCK_TRANSLATION") {
        let payload: TranslationPayload = serde_json::from_str(&mock_json)
            .map_err(|e| AskError::ResponseParse(format!("VCFKIT_MOCK_TRANSLATION parse: {e}")))?;
        return Ok(Translation {
            expression: payload.expression,
            reasoning: payload.reasoning,
            confidence: payload.confidence,
            caveats: payload.caveats,
            model: "mock".to_string(),
        });
    }

    let api_key = std::env::var("ANTHROPIC_API_KEY").map_err(|_| AskError::MissingApiKey)?;

    let model = std::env::var("VCFKIT_LLM_MODEL").unwrap_or_else(|_| DEFAULT_MODEL.to_string());

    let timeout_secs = std::env::var("VCFKIT_LLM_TIMEOUT_SECS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(DEFAULT_TIMEOUT_SECS);

    let system_prompt = build_system_prompt(schema);

    let first = call_api(&api_key, &model, &system_prompt, query, timeout_secs).await?;
    let payload = parse_response(&first.content)?;

    match validate_expression(&payload.expression) {
        Ok(()) => Ok(Translation {
            expression: payload.expression,
            reasoning: payload.reasoning,
            confidence: payload.confidence,
            caveats: payload.caveats,
            model: first.model,
        }),
        Err(parser_error) => {
            let retry_prompt = format!(
                "{}\n\nYour previous expression `{}` was rejected by the parser with: {}\n\
                 Please produce a corrected expression.",
                query, payload.expression, parser_error
            );
            let retry = call_api(
                &api_key,
                &model,
                &system_prompt,
                &retry_prompt,
                timeout_secs,
            )
            .await?;
            let retry_payload = parse_response(&retry.content)?;
            validate_expression(&retry_payload.expression).map_err(|e| {
                AskError::InvalidTranslation {
                    expression: retry_payload.expression.clone(),
                    error: e,
                }
            })?;
            Ok(Translation {
                expression: retry_payload.expression,
                reasoning: retry_payload.reasoning,
                confidence: retry_payload.confidence,
                caveats: retry_payload.caveats,
                model: retry.model,
            })
        }
    }
}

// ── API call ──────────────────────────────────────────────────────────────────

async fn call_api(
    api_key: &str,
    model: &str,
    system_prompt: &str,
    user_prompt: &str,
    timeout_secs: u64,
) -> Result<ApiCallResult, AskError> {
    let client = reqwest::Client::builder()
        .timeout(Duration::from_secs(timeout_secs))
        .build()?;

    let body = AnthropicRequest {
        model,
        max_tokens: MAX_TOKENS,
        temperature: 0.0,
        system: system_prompt,
        messages: vec![Message {
            role: "user",
            content: user_prompt,
        }],
    };

    let response = client
        .post(API_URL)
        .header("x-api-key", api_key)
        .header("anthropic-version", ANTHROPIC_VERSION)
        .header("content-type", "application/json")
        .json(&body)
        .send()
        .await
        .map_err(|e| {
            if e.is_timeout() {
                AskError::Timeout(timeout_secs)
            } else {
                AskError::NetworkError(e)
            }
        })?;

    if !response.status().is_success() {
        let status = response.status().as_u16();
        let body = response
            .text()
            .await
            .unwrap_or_else(|_| "<no body>".to_string());
        return Err(AskError::ApiError { status, body });
    }

    let parsed: AnthropicResponse = response
        .json()
        .await
        .map_err(|e| AskError::ResponseParse(e.to_string()))?;

    let content = parsed
        .content
        .into_iter()
        .find(|b| b.block_type == "text")
        .and_then(|b| b.text)
        .ok_or_else(|| AskError::ResponseParse("no text block in response".to_string()))?;

    Ok(ApiCallResult {
        content,
        model: parsed.model,
    })
}

// ── response parsing ──────────────────────────────────────────────────────────

fn parse_response(content: &str) -> Result<TranslationPayload, AskError> {
    let json_str = extract_json(content);
    serde_json::from_str(&json_str)
        .map_err(|e| AskError::ResponseParse(format!("{e}; raw response: {content}")))
}

/// Pull a JSON object out of a response that may be wrapped in prose or
/// markdown code fences.
pub fn extract_json(content: &str) -> String {
    let trimmed = content.trim();
    if let Some(stripped) = trimmed.strip_prefix("```json") {
        if let Some(end) = stripped.rfind("```") {
            return stripped[..end].trim().to_string();
        }
    }
    if let Some(stripped) = trimmed.strip_prefix("```") {
        if let Some(end) = stripped.rfind("```") {
            return stripped[..end].trim().to_string();
        }
    }
    if let (Some(start), Some(end)) = (trimmed.find('{'), trimmed.rfind('}')) {
        if start <= end {
            return trimmed[start..=end].to_string();
        }
    }
    trimmed.to_string()
}

fn validate_expression(expr: &str) -> Result<(), String> {
    vcfkit_core::filter::FilterExpression::parse(expr)
        .map(|_| ())
        .map_err(|e| e.to_string())
}

// ── system prompt ─────────────────────────────────────────────────────────────

pub fn build_system_prompt(schema: &HeaderSchema) -> String {
    let mut p = String::from(
        "You translate natural-language queries about genomic variants into vcfkit filter expressions.\n\n\
         ## Filter expression grammar\n\
         - Fields: INFO/<name>, FORMAT/<name>, CHROM, POS, QUAL, FILTER\n\
         - Comparison operators: < <= > >= == !=\n\
         - Logical operators: && || ! (parentheses allowed)\n\
         - Substring operators: ~ (contains) !~ (does not contain)\n\
         - Values: unquoted numbers (1.5, 42, -1), single-quoted strings ('PASS', 'chr17')\n\
         - Examples: INFO/AF < 0.01  |  QUAL > 30 && FILTER == 'PASS'  |  CHROM == 'chr1'\n\n\
         ## Available fields in this VCF\n",
    );

    if !schema.info_fields.is_empty() {
        p.push_str("\n### INFO fields\n");
        for f in &schema.info_fields {
            p.push_str(&format!(
                "  INFO/{} (Number={}, Type={}): {}\n",
                f.id, f.number, f.ty, f.description
            ));
        }
    }
    if !schema.format_fields.is_empty() {
        p.push_str("\n### FORMAT fields\n");
        for f in &schema.format_fields {
            p.push_str(&format!(
                "  FORMAT/{} (Number={}, Type={}): {}\n",
                f.id, f.number, f.ty, f.description
            ));
        }
    }
    if !schema.contigs.is_empty() {
        let display = if schema.contigs.len() > 12 {
            format!(
                "{} ... ({} total)",
                schema.contigs[..12].join(", "),
                schema.contigs.len()
            )
        } else {
            schema.contigs.join(", ")
        };
        p.push_str(&format!("\n### Contigs\n  {display}\n"));
    }

    p.push_str(
        "\n## Output format\n\
         Output ONLY a JSON object — no prose before or after:\n\
         {\n  \
           \"expression\": \"<valid filter expression using only the fields listed above>\",\n  \
           \"reasoning\": \"<1-2 sentences explaining the translation>\",\n  \
           \"confidence\": <0.0 to 1.0>,\n  \
           \"caveats\": [\"<caveat>\", ...]\n\
         }\n\n\
         ## Rules\n\
         1. Only reference fields listed above. NEVER invent field names.\n\
         2. For multi-allelic INFO fields, any-element semantics apply automatically — INFO/AF < 0.01 matches if any allele has AF < 0.01.\n\
         3. For ambiguous thresholds (e.g. 'rare', 'high quality'), use standard genomic conventions and note them in caveats.\n\
         4. If the query cannot be answered with available fields, explain in caveats, set confidence below 0.5, and produce the nearest valid expression.\n\
         5. expression must never be empty.\n\n\
         ## Examples\n\n\
         Query: \"rare variants\"\n\
         {\"expression\": \"INFO/AF < 0.01\", \"reasoning\": \"Rare variants conventionally means allele frequency below 1%.\", \"confidence\": 0.9, \"caveats\": [\"Using AF < 0.01 as the rarity threshold; adjust as needed.\"]}\n\n\
         Query: \"high-quality PASS variants on chromosome 17\"\n\
         {\"expression\": \"QUAL > 30 && FILTER == 'PASS' && CHROM == 'chr17'\", \"reasoning\": \"QUAL > 30 is a standard quality cutoff (Phred scale). FILTER == PASS ensures all caller filters passed.\", \"confidence\": 0.95, \"caveats\": []}\n\n\
         Query: \"missense variants in BRCA1\"\n\
         {\"expression\": \"INFO/CSQ ~ 'missense_variant' && INFO/CSQ ~ 'BRCA1'\", \"reasoning\": \"Assumes VEP CSQ annotation. Substring match finds records where annotation contains both terms.\", \"confidence\": 0.35, \"caveats\": [\"Requires VEP INFO/CSQ annotation field.\", \"Match is substring-based and may capture related terms.\"]}\n"
    );

    p
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_plain_json() {
        let s =
            r#"{"expression": "QUAL > 30", "reasoning": "ok", "confidence": 0.9, "caveats": []}"#;
        assert_eq!(extract_json(s), s);
    }

    #[test]
    fn extract_fenced_json() {
        let s = "```json\n{\"a\": 1}\n```";
        assert_eq!(extract_json(s), "{\"a\": 1}");
    }

    #[test]
    fn extract_plain_fence() {
        let s = "```\n{\"a\": 1}\n```";
        assert_eq!(extract_json(s), "{\"a\": 1}");
    }

    #[test]
    fn extract_json_from_prose() {
        let s = "Here is the result:\n\n{\"expression\": \"QUAL > 30\"}\n\nHope that helps.";
        let extracted = extract_json(s);
        assert!(extracted.starts_with('{'));
        assert!(extracted.ends_with('}'));
        assert!(extracted.contains("QUAL > 30"));
    }

    #[test]
    fn validate_good_expressions() {
        assert!(validate_expression("QUAL > 30").is_ok());
        assert!(validate_expression("INFO/AF < 0.01 && FILTER == 'PASS'").is_ok());
        assert!(validate_expression("CHROM == 'chr17' || CHROM == 'chr22'").is_ok());
    }

    #[test]
    fn validate_bad_expressions() {
        assert!(validate_expression("QUAL >>> 30").is_err());
        assert!(validate_expression("INVENTED_OPERATOR INFO/AF").is_err());
    }

    #[test]
    fn build_prompt_includes_info_fields() {
        let schema = HeaderSchema {
            info_fields: vec![FieldDef {
                id: "AF".to_string(),
                number: "A".to_string(),
                ty: "Float".to_string(),
                description: "Allele frequency".to_string(),
            }],
            format_fields: vec![],
            contigs: vec![],
        };
        let prompt = build_system_prompt(&schema);
        assert!(prompt.contains("INFO/AF"));
        assert!(prompt.contains("Allele frequency"));
    }

    #[test]
    fn build_prompt_truncates_long_contig_list() {
        let schema = HeaderSchema {
            info_fields: vec![],
            format_fields: vec![],
            contigs: (1..=25).map(|i| format!("chr{i}")).collect(),
        };
        let prompt = build_system_prompt(&schema);
        assert!(prompt.contains("25 total"));
    }
}
