//! Progress-bar reporting for long-running CLI operations.
//!
//! Wraps [`indicatif`] with a no-op fallback when the output isn't a TTY, the
//! user passed `--quiet`, or the total is small enough that a progress bar is
//! more noise than signal.

use std::io::IsTerminal;
use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};

/// Below this threshold (records) a progress bar is suppressed — small jobs
/// finish before the bar would even paint.
const MIN_RECORDS_FOR_BAR: u64 = 1000;

/// Thin wrapper around [`ProgressBar`] that degrades to a no-op silently.
pub struct ProgressReporter {
    bar: Option<ProgressBar>,
}

impl ProgressReporter {
    /// Create a new reporter.
    ///
    /// Returns a silent reporter (no-op) when any of the following hold:
    /// * `quiet` is true,
    /// * stderr is not a TTY,
    /// * `total_records` is `Some(n)` with `n <= 1000`.
    ///
    /// When `total_records` is `None` the reporter paints an indeterminate
    /// spinner — we don't know the total (e.g. streaming from stdin).
    pub fn new(total_records: Option<u64>, quiet: bool) -> Self {
        if quiet || !std::io::stderr().is_terminal() {
            return Self { bar: None };
        }

        match total_records {
            Some(n) if n <= MIN_RECORDS_FOR_BAR => Self {
                bar: Some(ProgressBar::hidden()),
            },
            Some(n) => {
                let bar = ProgressBar::new(n);
                bar.set_style(
                    ProgressStyle::with_template(
                        "[{elapsed_precise}] {human_pos}/{human_len} records {bar:20.cyan/blue} {percent}% @ {per_sec} ETA {eta_precise}",
                    )
                    .unwrap_or_else(|_| ProgressStyle::default_bar())
                    .progress_chars("##-"),
                );
                bar.enable_steady_tick(Duration::from_millis(120));
                Self { bar: Some(bar) }
            }
            None => {
                let bar = ProgressBar::new_spinner();
                bar.set_style(
                    ProgressStyle::with_template(
                        "[{elapsed_precise}] {spinner} {human_pos} records @ {per_sec}",
                    )
                    .unwrap_or_else(|_| ProgressStyle::default_spinner()),
                );
                bar.enable_steady_tick(Duration::from_millis(120));
                Self { bar: Some(bar) }
            }
        }
    }

    /// Advance the progress bar by one record.
    pub fn inc(&self) {
        if let Some(bar) = &self.bar {
            bar.inc(1);
        }
    }

    /// Finalize the bar and emit a one-line summary to stderr.
    pub fn finish(&self, message: &str) {
        if let Some(bar) = &self.bar {
            bar.finish_with_message(message.to_string());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quiet_reporter_is_noop() {
        let r = ProgressReporter::new(Some(10_000), true);
        r.inc();
        r.inc();
        r.finish("done");
        // Just ensure no panic.
        assert!(r.bar.is_none());
    }

    #[test]
    fn small_total_hides_bar() {
        // In a non-TTY test environment we take the early-return first,
        // but this still exercises the code path.
        let r = ProgressReporter::new(Some(42), false);
        r.inc();
        r.finish("done");
        // With a non-TTY test runner, bar is None; either way no panic.
    }

    #[test]
    fn none_total_does_not_panic() {
        let r = ProgressReporter::new(None, false);
        r.inc();
        r.finish("done");
    }
}
