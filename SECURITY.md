# Security Policy

## Supported Versions

vcfkit is in alpha. Security fixes are applied to the latest release only.

| Version | Supported |
|---------|-----------|
| 0.3.x   | ✅ |
| < 0.3   | ❌ |

## Reporting a Vulnerability

**Please do not open a public GitHub issue for security vulnerabilities.**

Email: prasadkhake@gmail.com

Include:
- Description of the vulnerability
- Steps to reproduce
- Potential impact
- Any suggested fix (optional)

You'll receive a response within 72 hours. If the issue is confirmed, a fix will be released as quickly as possible and you'll be credited in the release notes (unless you prefer to remain anonymous).

## Scope

vcfkit processes local files and writes to stdout. It has no server component and no persistent state beyond `~/.config/vcfkit/config.toml`.

The one network operation is `vcfkit filter --ask`, which sends your query text and VCF header schema to the Anthropic API. Variant data (CHROM, POS, REF, ALT, genotypes) is never sent. See [Privacy](https://vcfkit.dev/privacy) for full details.
