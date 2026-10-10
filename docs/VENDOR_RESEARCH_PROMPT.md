# Prompt: vendor sequence-rule research for CLIPPR

*Hand this to a web-research agent verbatim. Self-contained — it assumes no knowledge of the
project.*

> **Revision 2.** Revision 1 was given to an agent that returned a well-formatted report in
> which one cited URL was a 404, another contained none of the four quotes attributed to it,
> and the "verbatim" quotes were paraphrases with the numbers preserved. The numbers were
> mostly correct, which made it worse: it looked checked and was not.
>
> This revision does not ask more politely. It requires an anchor string per claim that is
> **mechanically checked against the live page** by
> `validation/check_vendor_claims.py`. Rows whose anchor is not found on the page are
> rejected automatically. Read "How your work will be checked" before you start.

---

## Your task

Find every published rule, threshold and specification from **IDT** and **Twist Bioscience**
that bears on whether a synthetic DNA sequence will be accepted, rejected, flagged, or priced
differently. Report them in the machine-checkable format below.

## How your work will be checked

Every claim you submit carries a `verbatim_anchor`: a contiguous run of **8 to 20 words copied
character-for-character** from the page at `url`. After you deliver, a script fetches each URL
and searches for each anchor. Outcomes:

| Result | Meaning |
|---|---|
| anchor found on the page | claim accepted |
| page fetched, **anchor absent** | claim **rejected as unsupported**, and every other claim sharing that URL is re-examined |
| page unreachable | claim marked **unverified** — not a failure, flagged for a human |

So:

- **Copy, do not retype.** Normalised whitespace is fine; changed words are not. "homopolymeric
  runs of 10 or more As and Ts" is an anchor. "homopolymeric runs of 10 or more consecutive As
  or Ts" is a fabrication, even though it means the same thing and the number is right.
- **Pick a distinctive anchor.** "GC content" appears on every page and proves nothing. Choose
  a phrase that occurs once.
- **A claim you cannot anchor does not go in the table.** Put it in "Leads" instead, clearly
  marked. A short table of solid rows beats a long table that fails the check.
- **If a page is behind a sign-in wall, blocked, or will not load, say so and move on.** That
  is a reportable finding. Do not substitute a search-engine summary of the page and cite the
  page. That is the specific failure that produced revision 1.
- **If the content is in a linked PDF rather than on the page, cite the PDF's URL**, not the
  landing page that links to it.

Inventing a plausible URL, or reconstructing a quote from memory of what a page probably says,
wastes more of our time than returning nothing. We will find it.

## Context: what we order, and why the details matter

We design synthetic genes and split them into fragments for Golden Gate assembly.

- **Vendor: IDT. Product: oPools Oligo Pools, 50 pmol scale.**
- **Our fragments are 257–264 nt**, ordered as single-stranded DNA oligos in one pooled tube.
- Typical order: ~42 oligos, ~11,000 bases total.
- Each fragment carries a fixed wrapper at both ends: a pad, a BsaI site (`GGTCTC`), a spacer
  and a 4-nt overhang. Those bases cannot change.
- We can change **synonymous codons** in the coding region, and nothing else.

**Twist is a comparator, not our vendor.** It publishes more numbers than IDT, which is why we
care. A Twist rule is not binding on an IDT order and must never be presented as though it
were.

## The substantive failure mode

Our existing GC rule — local GC 35–65% over any 50 bp window, enforced as a hard rejection —
matches exactly one document: **Twist's codon-optimiser guidance**. Wrong vendor, wrong
product, and explicitly a design *target* rather than an acceptance criterion. Meanwhile:

- Twist's **stated acceptance** threshold for a 50 bp window is 10–90%, five times wider.
- Twist's **only hard rules** are "avoid homopolymers ≥ 14 bp" and "no CcdB".
- Twist's **Clonal Genes TechNote** uses a different mechanism entirely: a *spread* rule
  (highest minus lowest 50 bp window ≤ 50 points), not a band.
- **IDT publishes no GC window rule at all** for any product we have found.

Four Twist documents, four different 50 bp GC rules. We adopted the strictest by accident and
enforced it as a gate. Hence the obsession with which product and which enforcement level.

## Output format

Deliver **two files**.

### File 1 — `vendor_claims.jsonl`, one JSON object per line

```json
{"vendor":"IDT","product":"oPools Oligo Pools","parameter":"oligo length","value":"40-350 bases","enforcement":"hard","url":"https://...","date_read":"2026-10-11","verbatim_anchor":"copied word for word from the page","note":""}
```

Fields, all required except `note`:

- `vendor` — IDT, Twist, GenScript, Eurofins, Azenta, Thermo Fisher
- `product` — the exact product the rule governs. **oPools, Ultramer, gBlocks, eBlocks,
  Clonal Genes, Express Genes, Complex Genes and the codon optimiser have different rules.**
  A gBlocks rule says nothing about oPools. If a page does not say which product, write
  `"unspecified"`.
- `parameter` — e.g. "local GC window", "homopolymer run", "secondary structure"
- `value` — the numbers, with units and window size
- `enforcement` — exactly one of:
  - `hard` — the vendor states the order is refused
  - `complexity_flag` — affects eligibility, price tier or turnaround, not acceptance
  - `design_guidance` — recommended, not enforced
  - `optimiser_target` — what their design tool aims for
  - `scored` — fed to a model with no published threshold
  - `level_not_stated` — the page gives a number without saying what it does
  Use the vendor's own words to decide. Do not promote guidance to a rule because it has a
  number in it. `level_not_stated` is a perfectly good answer and is preferred to a guess.
- `url` — the page carrying the anchor. If it is a PDF, the PDF's URL.
- `date_read` — ISO date
- `verbatim_anchor` — 8–20 words, copied exactly
- `note` — optional; contradictions, caveats, scope

### File 2 — `vendor_report.md`

1. **Contradictions** — every case where two documents disagree. Quote both, with both URLs.
   Do **not** resolve them: do not pick the stricter, the looser or the newer. We decide.
2. **Absences** — what you looked for and did not find, and the specific pages you checked.
   *"No GC rule for oPools after checking the product page, the spec table, the FAQ list and
   the ordering guidelines"* is a valuable result.
3. **Blocked** — pages you could not reach, and why (sign-in, 403, redirect loop, JS-only).
4. **Corrections** — anything in "What we currently believe" below that is wrong or outdated.
5. **Leads** — claims you could not anchor. Clearly marked as unverified, with where you think
   the source is.

Do **not** recommend which rule we should adopt. Report what exists.

## What to find

### A. IDT — our actual vendor, highest priority

1. **oPools Oligo Pools** — every published specification:
   - length range, oligos per pool, scales
   - **any sequence-content rule at all**: GC, homopolymer, repeat, secondary structure,
     complexity. We believe none exists. Confirm or refute, and name the pages you checked.
   - QC, quantification, dropout rate, and any statement linking sequence complexity to yield
   - pricing tiers, modification fees (especially 5′ phosphorylation), turnaround, storage
2. **The IDT complexity checker** — what it screens for; any thresholds published outside the
   sign-in wall; whether it applies to oPools or only to gBlocks and genes.
3. **IDT SciTools API** — public documentation? Does it expose complexity scoring, or
   OligoAnalyzer hairpin/dimer ΔG, programmatically? Endpoints, auth, rate limits.
4. **OligoAnalyzer** — the documented maximum input length **per function** (Analyze, Hairpin,
   Self-Dimer, Hetero-Dimer), and any stated validity range. We have observed it accept a
   257-base sequence for Analyze while stating its melt-temperature model is accurate only for
   8–60 bases. A third-party claim of a 255-base ceiling is unverified; settle it.
5. **gBlocks / eBlocks / Ultramer** — their rules, labelled as theirs, for comparison only.
6. **Biosecurity and sequence screening** — what IDT refuses on policy grounds, and how
   screening works.

### B. Twist Bioscience — comparator

7. Every distinct GC rule it publishes, across **all** products and documents, each with its
   product and enforcement level. We know of four and expect more.
8. Its hard rules, complexity flags, "Impossible" scoring factors, and how it describes the
   accept/reject decision.
9. Homopolymer, repeat and secondary-structure thresholds, per product.
10. The Clonal Genes design TechNote (DOC-001081 or current revision) — every rule in it.
    **Cite the PDF URL, not the landing page**; the landing page carries no rules.

### C. Secondary structure — the open question that matters most

11. **Does any vendor publish a secondary-structure threshold for oligo pools?** We have found
    only Twist's "strong hairpins (ΔG < −8) in the first 48 base pairs", which is
    gene-synthesis optimiser guidance, not an oPools rule. We are borrowing it for want of
    anything better and would rather not.
12. Any vendor statement linking GC content to secondary structure, or calling GC a proxy
    for it.
13. Any published data relating sequence complexity to oligo-pool dropout rate.

### D. Alternatives — brief

14. **GenScript, Azenta/GENEWIZ, Eurofins, Thermo Fisher**: does any publish *numeric*
    acceptance criteria for oligo pools or long oligos, rather than "contact us"? Include the
    maximum oligo length each supports — if a vendor caps below 264 nt it cannot make our
    fragments at all, which settles it. One paragraph each.
15. Vendors offering iGEM team discounts or sponsorship, with terms and the source. Prefer
    igem.org's own partner pages over the vendor's marketing.

## What we currently believe

Treat all of this as **unverified**. Correcting any row is a useful result. Rows marked
*first-hand* were read directly from the live page by the project team on 2026-10-11; if you
disagree with one, say so explicitly and show the anchor.

| Claim | Status |
|---|---|
| IDT oPools: 40–350 bases; 2–384 oligos at 50 pmol, 10–2,000 at 10 pmol, 100–20,000 at 1 pmol | first-hand |
| IDT oPools: no published GC, complexity or structure rule | first-hand, want independent confirmation |
| IDT oPools: QC "None"; average dropout 0.8% over ~900,000 oligos; dropout "can vary by length and sequence complexity" | first-hand |
| IDT gBlocks: GC 25–75%; homopolymer 10+ A/T or 6+ G/C; acceptance "complex and multifactorial" | first-hand, from the **sequence-motifs FAQ**, not the product page |
| IDT Ultramer caps at 200 bases | first-hand |
| IDT pricing: €109 first 3,300 bases, then €0.038 / €0.025 / €0.013 per base; €1.63/oligo for 5′ phosphorylation (EU site) | first-hand |
| Twist: only hard rules are homopolymer ≥14 bp and no CcdB | first-hand |
| Twist Express complexity flag: 50 bp window below 10% or above 90% | first-hand |
| Twist codon optimiser: global 25–65%; local 35–65% over 50 bp; hairpins ΔG < −8 in the first 48 bp | first-hand |
| Twist Clonal Genes TechNote: global 25–65%; 50 bp window *spread* ≤ 50 points; homopolymer ≤13 nt; avoid repeats ≥20 bp or Tm ≥60 °C | first-hand, from the PDF |

**Access note:** IDT's site redirects automated fetchers in a loop; a real browser session
reaches it. `www.idtdna.com` redirects to `eu.idtdna.com` from Europe. If you cannot reach a
page, name it in "Blocked" rather than substituting a summary of it.

## One last thing

If you are about to write "vendors generally recommend" or "industry standard is" — stop. That
sentence has no source behind it and is how we got here.
