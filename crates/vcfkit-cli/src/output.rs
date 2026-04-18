//! Progress-bar reporting for long-running CLI operations.
//!
//! Wraps [`indicatif`] with a no-op fallback when the output isn't a TTY, the
//! user passed `--quiet`, `--no-progress`, or the total is small enough that a
//! progress bar is more noise than signal.
//!
//! The draw target is **always** stderr so progress output never contaminates
//! a piped VCF stream on stdout.

use std::io::IsTerminal;
use std::time::Duration;

use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};

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
    /// * `quiet` or `no_progress` is true,
    /// * stderr is not a TTY,
    /// * `total_records` is `Some(n)` with `n <= 1000`.
    ///
    /// When `total_records` is `None` the reporter paints an indeterminate
    /// spinner — we don't know the total (e.g. streaming from stdin).
    ///
    /// The draw target is explicitly set to stderr so that the bar never
    /// writes to stdout and cannot corrupt piped VCF output.
    pub fn new_with_flags(total_records: Option<u64>, quiet: bool, no_progress: bool) -> Self {
        if quiet || no_progress || !std::io::stderr().is_terminal() {
            return Self { bar: None };
        }

        match total_records {
            Some(n) if n <= MIN_RECORDS_FOR_BAR => Self { bar: None },
            Some(n) => {
                let bar = ProgressBar::new(n);
                bar.set_draw_target(ProgressDrawTarget::stderr());
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
                bar.set_draw_target(ProgressDrawTarget::stderr());
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
        let r = ProgressReporter::new_with_flags(Some(10_000), true, false);
        r.inc();
        r.inc();
        r.finish("done");
        // Just ensure no panic.
        assert!(r.bar.is_none());
    }

    #[test]
    fn no_progress_flag_forces_noop() {
        let r = ProgressReporter::new_with_flags(Some(10_000), false, true);
        r.inc();
        r.finish("done");
        assert!(r.bar.is_none());
    }

    #[test]
    fn small_total_hides_bar() {
        // Both the non-TTY early-return and the small-record branch now
        // return bar: None — assert that in either case bar is None.
        let r = ProgressReporter::new_with_flags(Some(42), false, false);
        r.inc();
        r.finish("done");
        assert!(r.bar.is_none());
    }

    #[test]
    fn none_total_does_not_panic() {
        let r = ProgressReporter::new_with_flags(None, false, false);
        r.inc();
        r.finish("done");
    }
}
