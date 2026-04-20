# vcfkit.dev — Design Notes

Round 1 design direction proposals. Review and pick one option per section
before any code is touched.

---

## 1. Color Palette

Three options. Each has a signature color and one accent, with dark-mode
values. The current Starlight default is blue — all three options move away
from that.

---

### Option A — Muted Teal

```
Signature  #2E9688   (dark teal)
Accent     #56D4C8   (bright teal, for highlights / active states)
Dark bg    #0E1F1E
Dark text  #D6EFED
```

Light mode counterparts (if Starlight light mode is ever enabled):
```
Light bg   #F2FAFA
Light text #0E2523
```

**Rationale:** Teal reads as "genomics-adjacent" — it's the color family
UCSC Genome Browser uses for annotation tracks. Rare among dev tools (most
pick blue or green). Calm and readable. The two shades give enough range
for buttons, active states, and code highlights without looking garish.

---

### Option B — Warm Amber

```
Signature  #C97C3A   (warm amber, accessible on dark)
Accent     #F0A55C   (lighter for hover states and highlights)
Dark bg    #131009
Dark text  #F0E8DC
```

Light mode:
```
Light bg   #FDF8F2
Light text #2A1A08
```

**Rationale:** Nobody in bioinformatics uses amber. It's warm, distinctive,
and creates strong contrast on dark backgrounds. Feels handcrafted rather
than default-framework. The warm tone also makes the site feel approachable
and human — not cold/corporate.

---

### Option C — Muted Violet

```
Signature  #7455BA   (muted violet)
Accent     #A27FE8   (lighter for highlights)
Dark bg    #0F0D1A
Dark text  #E8E3F5
```

Light mode:
```
Light bg   #F8F5FF
Light text #1A1228
```

**Rationale:** Academic/research-adjacent. The purple family is used by
Ensembl Genome Browser and feels "scientific institution" without being
stuffy. Creates a sense of depth. Among the three options, most distinctive
at a glance.

---

**→ Pick one: A (teal), B (amber), or C (violet)**

---

## 2. Typography

Two pairings. Both use free/open-source fonts available from Google Fonts
or with permissive licenses. No licensing risk.

---

### Pairing 1 — Geist Sans + JetBrains Mono

| Role     | Font            | Weight |
|----------|-----------------|--------|
| Headline | Geist Sans      | 700    |
| Body     | Geist Sans      | 400    |
| Code     | JetBrains Mono  | 400    |

**Rationale:** Geist Sans (Vercel's open-source typeface) has tight spacing
and a technical feel without looking like a startup clone. It reads
professionally in dense technical contexts. JetBrains Mono has optical
ligatures for `>=`, `!=`, `&&`, `||` — these appear constantly in VCF
filter expressions, so they render cleanly.

---

### Pairing 2 — IBM Plex Sans + IBM Plex Mono

| Role     | Font          | Weight |
|----------|---------------|--------|
| Headline | IBM Plex Sans | 600    |
| Body     | IBM Plex Sans | 400    |
| Code     | IBM Plex Mono | 400    |

**Rationale:** IBM Plex has a slightly more "scientific paper" quality.
The Sans and Mono are designed as a matched pair — they share optical
metrics so inline code snippets feel native to the prose. Slightly more
formal, which fits the "validated against bcftools" tone.

---

**→ Pick one: 1 (Geist + JetBrains) or 2 (IBM Plex)**

---

## 3. Wordmark

Two treatments. Both are CSS-only — no custom fonts required beyond the
typography pick above.

---

### Treatment 1 — Monospace brackets

```
[vcfkit]
```

The brackets are rendered in the accent color at regular weight; `vcfkit`
is semibold. Uses the code font so the name visually encodes "this is a
CLI tool."

```css
/* mockup */
.wordmark {
  font-family: var(--font-mono);
  font-weight: 600;
  letter-spacing: 0.02em;
  color: var(--color-text);
}
.wordmark-bracket {
  font-weight: 400;
  color: var(--color-accent);
}
/* renders as: [vcfkit] where brackets are accent-colored */
```

**Rationale:** The bracket instantly signals CLI. It's subtle — just two
characters. No custom glyph needed.

---

### Treatment 2 — Two-weight split

```
vcf kit
```

`vcf` is rendered at normal weight in a muted text color; `kit` is bold in
the primary text color. A single-pixel vertical rule or zero-width spacer
separates them. Uses the headline font.

```css
/* mockup */
.wordmark-vcf {
  font-family: var(--font-sans);
  font-weight: 400;
  color: var(--color-text-muted);
  letter-spacing: -0.01em;
}
.wordmark-kit {
  font-family: var(--font-sans);
  font-weight: 700;
  color: var(--color-text);
  letter-spacing: -0.01em;
}
/* renders: vcf (muted) + kit (bold), no space */
```

**Rationale:** `vcf` is the file format everyone in the audience knows;
`kit` is the differentiator. The weight split draws the eye to `-kit` as
the operative word without being loud.

---

**→ Pick one: 1 (brackets) or 2 (weight split)**

---

## 4. Information Architecture

### Landing page flow

One scrollable page. Sections in order:

```
1. Hero
   — headline, benchmark stat strip, three CTAs

2. Demo (full-width)
   — immediately below the fold, the main attraction

3. What is this
   — two-paragraph plain-English positioning

4. Three operations
   — cards: normalize · liftover · filter

5. Natural-language queries (--ask)
   — AI feature callout, realistic example, confirmation flow

6. Benchmarks
   — table with methodology link and nightly CI link

7. Credits
   — bcftools, noodles, Tan et al., UCSC — with descriptions

8. Trust signals
   — four-bullet block: correctness, status, clinical disclaimer, known differences

9. Get in touch
   — GitHub Issues, Discussions, email placeholder

10. Footer
    — single line: Built by Prasad Khake · MIT · GitHub · crates.io · Docs
```

### Nav (always visible, top)

```
vcfkit [wordmark]          Docs ▾   GitHub   crates.io
```

Docs dropdown:
- Install
- normalize
- liftover
- filter
- Benchmarks
- Known differences
- Privacy
- Credits

No sidebar on the landing page (it's a marketing page, not docs).
Sidebar appears on all /commands/*, /install, /benchmarks, etc. pages.

### Footer

```
Built by Prasad Khake · MIT License · GitHub · crates.io · Docs
© 2026
```

"Prasad Khake" links to personal site/profile (placeholder — provide URL).

---

## 5. Mobile Breakpoints

| Breakpoint | Layout changes |
|------------|----------------|
| 1280px     | Full layout. Two-column demo (input / output side by side). |
| 1024px     | Demo collapses to stacked (input above, output below). Three operations cards stay in row. |
| 768px      | Nav collapses to hamburger. Operations cards stack 1-column. Hero text scales down. |
| 480px      | Phone view. Demo shows a banner: "For the best experience, open on a larger screen" (demo still functional). Font sizes reduce. Footer goes to 2 lines. |

---

## 6. Subtle Bioinformatics Touches

Five candidates. Pick 3–5 to include.

---

### Touch 1 — Monospace genomic coordinates

Any time a genomic coordinate appears in text (`chr17:43,050,000`,
`REF→ALT`), render it in the mono font with a slight `background: var(--code-bg)`
pill. Not a code block — just enough to signal "this is a position."

Example: "BRCA1 spans `chr17:43,044,295–43,125,483` on the reference genome."

---

### Touch 2 — ACGT base accent stripe

A single four-segment horizontal stripe used once: in the hero section as a
decorative underline under the wordmark, or as a section separator.
Segments: A (accent green), C (blue), G (amber), T (red), each 25% width,
about 3px tall, rounded caps.

Used exactly once. Not animated. Could also appear as the active tab
underline in the demo.

```css
.acgt-stripe {
  height: 3px;
  border-radius: 2px;
  background: linear-gradient(
    to right,
    #4ADE80 0% 25%,    /* A — green */
    #60A5FA 25% 50%,   /* C — blue  */
    #FBBF24 50% 75%,   /* G — amber */
    #F87171 75% 100%   /* T — red   */
  );
}
```

---

### Touch 3 — Terminal-style section separators

Between major landing-page sections, a thin `1px` rule styled to look like
a terminal divider — `────────` in a dim color, rendered with `::before`
pseudo-element or as a `<hr>` styled with repeating dashes. Subtle.

---

### Touch 4 — File extension badges

When showing filenames inline (`.vcf`, `.bcf`, `.vcf.gz`), render the
extension as a small pill badge in the code font. Color-coded:
- `.vcf` — neutral (gray)
- `.vcf.gz` — muted blue (compressed)
- `.bcf` — muted amber (binary)

---

### Touch 5 — Chromosome-style run progress

While the demo is running (WASM processing), the run button's inner span
shows a subtle left-to-right fill animation — like a chromosome ideogram
coverage bar. Contained to the button, ~2px tall at the bottom of the
button, filling over ~200ms.

Not a spinner (spinners feel anxious). A quiet fill feels like progress.

---

**→ Pick which touches to include (recommend 3–4 max for subtlety)**

---

## Summary — what you pick

| Decision | Options |
|----------|---------|
| Color palette | A (teal) · B (amber) · C (violet) |
| Typography | 1 (Geist + JetBrains) · 2 (IBM Plex) |
| Wordmark | 1 (brackets) · 2 (weight split) |
| Bio touches | pick 3–4 from the 5 above |
| Personal site URL | for the footer "Prasad Khake" link |
| Contact email | for the "Get in touch" section |

Once you reply with picks, Round 2 begins: the landing page rebuild.
No code has been changed yet.
