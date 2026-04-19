pub mod filter;
pub mod liftover;
pub mod normalize;

use std::path::Path;

/// Reject BCF output paths before doing any I/O — BCF writing is not yet
/// implemented and silently falling back to VCF would surprise users.
pub(crate) fn reject_bcf_output(path: Option<&Path>) -> anyhow::Result<()> {
    if let Some(p) = path {
        if p.extension().and_then(|e| e.to_str()) == Some("bcf") {
            anyhow::bail!(
                "BCF output is not yet supported (planned for v0.2).\n\
                 Workaround: write VCF and convert with bcftools:\n\
                 \n  vcfkit {} ... | bcftools view -O b -o {}\n",
                std::env::args().nth(1).unwrap_or_default(),
                p.display()
            );
        }
    }
    Ok(())
}
