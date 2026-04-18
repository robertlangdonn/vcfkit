//! Opt-in telemetry scaffolding.
//!
//! We persist an `enabled: Option<bool>` tri-state (None = not yet asked)
//! in `~/.config/vcfkit/config.toml`. First-run, we prompt on stderr and
//! save the user's choice. Payloads are sent fire-and-forget via `curl` in
//! a background thread — zero added HTTP dependencies, and we never block
//! the CLI exit path.
//!
//! The payload deliberately omits anything identifying: no file paths,
//! no hostnames, no content, just command name, OS/arch, duration, and a
//! success bit.

use std::io::{self, BufRead, IsTerminal, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::thread;
use std::time::Duration;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

/// Compile-time default endpoint for telemetry events, overridable at build
/// time via the `VCFKIT_TELEMETRY_URL` environment variable, and again at
/// runtime by setting the same env var.
const DEFAULT_TELEMETRY_URL: &str = match option_env!("VCFKIT_TELEMETRY_URL") {
    Some(u) => u,
    None => "https://telemetry.vcfkit.dev/v1/event",
};

/// Persisted config (lives at `~/.config/vcfkit/config.toml`).
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct TelemetryConfig {
    /// `None` means the user has not yet been asked.
    pub enabled: Option<bool>,
}

impl TelemetryConfig {
    /// Resolve the on-disk config path honouring `$XDG_CONFIG_HOME`.
    pub fn config_path() -> PathBuf {
        let base = std::env::var_os("XDG_CONFIG_HOME")
            .map(PathBuf::from)
            .filter(|p| !p.as_os_str().is_empty())
            .or_else(|| std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".config")))
            .unwrap_or_else(|| PathBuf::from(".config"));
        base.join("vcfkit").join("config.toml")
    }

    /// Read the config from disk; return default (enabled = None) if missing
    /// or unparseable.
    pub fn load() -> Self {
        let path = Self::config_path();
        let Ok(text) = std::fs::read_to_string(&path) else {
            return Self::default();
        };
        parse_config(&text)
    }

    /// Write the config atomically-ish (create parent dir, overwrite file).
    pub fn save(&self) -> Result<()> {
        let path = Self::config_path();
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).with_context(|| {
                format!("failed to create config directory {}", parent.display())
            })?;
        }
        let body = render_config(self);
        std::fs::write(&path, body)
            .with_context(|| format!("failed to write config file {}", path.display()))?;
        Ok(())
    }

    /// Prompt the user (first run only) and return the final effective
    /// enabled flag. Writes the user's answer to disk so we never ask twice.
    ///
    /// When stdin isn't a TTY (e.g. the user is piping a VCF in) we skip the
    /// prompt and return `false` for this run without persisting — we'll ask
    /// again next time they run interactively.
    pub fn ensure_prompted(&mut self) -> bool {
        if let Some(enabled) = self.enabled {
            return enabled;
        }
        if !io::stdin().is_terminal() || !io::stderr().is_terminal() {
            // No way to ask; stay off for this run, ask again next time.
            return false;
        }
        let answer = prompt_user();
        self.enabled = Some(answer);
        // Best-effort: a failure to persist shouldn't break the main command.
        let _ = self.save();
        answer
    }
}

/// Show the opt-in prompt and read a single line from stdin.
fn prompt_user() -> bool {
    let _ = writeln!(
        io::stderr(),
        "vcfkit can send anonymous usage statistics to help improve the tool.\n\
         This includes: command name, input size bucket, duration, success/error.\n\
         It does NOT include: file contents, file paths, hostnames, or any identifying data.\n"
    );
    let _ = write!(io::stderr(), "Enable telemetry? [y/N]: ");
    let _ = io::stderr().flush();

    let mut line = String::new();
    let stdin = io::stdin();
    if stdin.lock().read_line(&mut line).is_err() {
        return false;
    }
    matches!(line.trim(), "y" | "Y" | "yes" | "YES" | "Yes")
}

// ── payload + send ───────────────────────────────────────────────────────────

/// Send an event in a detached background thread.
///
/// Fire-and-forget: if `curl` is missing, if the body can't be serialized, if
/// the thread panics — we swallow every error. Telemetry must never affect
/// the user-visible exit path.
pub fn send_event(command: &str, duration: Duration, success: bool) {
    let event = OwnedEvent {
        version: env!("CARGO_PKG_VERSION").to_string(),
        os: std::env::consts::OS.to_string(),
        arch: std::env::consts::ARCH.to_string(),
        command: command.to_string(),
        duration_ms: duration.as_millis() as u64,
        success,
    };

    thread::spawn(move || {
        let _ = post_event(&event);
    });
}

/// Owned companion to `Event<'a>` so we can move it into a background thread.
#[derive(Debug, Serialize)]
struct OwnedEvent {
    version: String,
    os: String,
    arch: String,
    command: String,
    duration_ms: u64,
    success: bool,
}

fn post_event(event: &OwnedEvent) -> Result<()> {
    let url = std::env::var("VCFKIT_TELEMETRY_URL")
        .unwrap_or_else(|_| DEFAULT_TELEMETRY_URL.to_string());
    let body = serde_json::to_string(event)?;

    // Only send if `curl` is on PATH; otherwise skip silently.
    if Command::new("curl")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_err()
    {
        return Ok(());
    }

    let _ = Command::new("curl")
        .arg("--silent")
        .arg("--max-time")
        .arg("5")
        .arg("-X")
        .arg("POST")
        .arg("-H")
        .arg("Content-Type: application/json")
        .arg("--data")
        .arg(body)
        .arg(url)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn(); // fire-and-forget — don't wait for curl to finish

    Ok(())
}

// ── tiny hand-rolled TOML read/write ─────────────────────────────────────────
//
// We have no runtime TOML dependency, and the config is a single table with
// one boolean key. The parser below handles that minimal surface (and ignores
// comments/blank lines) so users can hand-edit without surprise.

fn render_config(cfg: &TelemetryConfig) -> String {
    let mut out = String::from("[telemetry]\n");
    match cfg.enabled {
        Some(true) => out.push_str("enabled = true\n"),
        Some(false) => out.push_str("enabled = false\n"),
        None => out.push_str("# enabled = false\n"),
    }
    out
}

fn parse_config(text: &str) -> TelemetryConfig {
    let mut cfg = TelemetryConfig::default();
    let mut in_telemetry = false;

    for raw in text.lines() {
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('[') && line.ends_with(']') {
            in_telemetry = line == "[telemetry]";
            continue;
        }
        if !in_telemetry {
            continue;
        }
        let Some((key, value)) = line.split_once('=') else {
            continue;
        };
        let key = key.trim();
        let value = value.trim();
        if key == "enabled" {
            cfg.enabled = match value {
                "true" => Some(true),
                "false" => Some(false),
                _ => continue,
            };
        }
    }
    cfg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_enabled_true() {
        let cfg = parse_config("[telemetry]\nenabled = true\n");
        assert_eq!(cfg.enabled, Some(true));
    }

    #[test]
    fn parse_enabled_false() {
        let cfg = parse_config("[telemetry]\nenabled = false\n");
        assert_eq!(cfg.enabled, Some(false));
    }

    #[test]
    fn parse_missing_section() {
        let cfg = parse_config("# nothing here\n");
        assert_eq!(cfg.enabled, None);
    }

    #[test]
    fn parse_ignores_unrelated_section() {
        let cfg = parse_config("[other]\nenabled = true\n");
        assert_eq!(cfg.enabled, None);
    }

    #[test]
    fn render_roundtrip_true() {
        let before = TelemetryConfig {
            enabled: Some(true),
        };
        let text = render_config(&before);
        let after = parse_config(&text);
        assert_eq!(after.enabled, Some(true));
    }

    #[test]
    fn render_roundtrip_false() {
        let before = TelemetryConfig {
            enabled: Some(false),
        };
        let text = render_config(&before);
        let after = parse_config(&text);
        assert_eq!(after.enabled, Some(false));
    }

    #[test]
    fn config_path_uses_xdg() {
        // Sanity: path ends with the expected suffix regardless of
        // environment, and is non-empty.
        let p = TelemetryConfig::config_path();
        let s = p.to_string_lossy();
        assert!(s.ends_with("vcfkit/config.toml"), "got {}", s);
    }
}
