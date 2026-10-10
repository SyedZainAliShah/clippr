# Which vendor should the sequence policy be modelled on?

*Researched 11 October 2026 from primary vendor pages. Supersedes the single unsourced band.*

The GC band refusing all 200 de novo designs matches exactly one published document, and it is
the least applicable one we could have chosen. This records what each vendor actually publishes,
how our designs score against each, and what to model on.

**Verification levels.** **[F]** fetched and read directly. **[S]** quoted by a search engine
from the vendor's page; IDT's site returns a redirect loop to automated fetches, so its numbers
could not be page-verified here and should be confirmed in IDT's own order-entry checker.

---

## 1. What we actually order

**IDT oPools DNA, 50 pmol.** Not gene synthesis. `ordering.CURRENT_OPOOL_50PMOL` models it:
40–350 bases per oligo, 2–384 oligos per pool, sourced and dated.

Our fragments are **257–264 nt**, comfortably inside that range.

**IDT publishes no GC or complexity limit for oPools.** So there is no vendor sequence
requirement for our actual product. Everything below is guidance borrowed from adjacent
products.

## 2. What each document actually says

| Document | Global GC | Local 50 bp rule | Homopolymer | Level |
|---|---|---|---|---|
| **Twist Clonal Genes TechNote** DOC-001081 REV4 | 25–65% | **spread** between highest and lowest window ≤ 50 pts | ≤ 13 nt | design guide **[F]** |
| **Twist codon-optimiser FAQ** ← *our band* | 25–65% | **35–65% absolute** | not stated | optimiser guidance **[F]** |
| **Twist Express** complexity flag | — | 10–90% absolute | run > 30 bp | complexity flag, not rejection **[F]** |
| **Twist "the only hard rules"** | — | — | **≥ 14 bp forbidden** | hard **[F]** |
| **Twist "Impossible" score factors** | 25–65% | — | < 10 | ML-scored **[F]** |
| **IDT gBlocks** sequence motifs FAQ | 25–75% | none published | 10+ A/T, 6+ G/C | guidance **[S]** |
| **IDT oPools** (our product) | **none published** | **none published** | **none published** | — **[F]** |

Two things to notice. Twist's **only** stated hard rules are a 14 bp homopolymer and no CcdB.
And Twist's *acceptance* threshold for a 50 bp window — 10–90% — is roughly five times wider
than the band we enforce as a rejection gate.

Twist also states that rejection is scored by a machine-learning model and that *"most of the
time, there is no single reason for the rejection of a sequence."* No threshold set predicts
acceptance, which is why `STRICT_LEGACY.unresolved` has always said so.

## 3. How our 1050 ordered fragments score

Measured on the pinned corpus:

| Rule | Result |
|---|---|
| Twist TechNote: global GC 25–65% | **0 of 1050 violate** (range 0.565–0.611) |
| Twist TechNote: window spread ≤ 50 pts | **0 of 1050 violate** (range 0.080–0.200) |
| Twist TechNote: homopolymer ≤ 13 nt | **0 of 1050 violate** (longest 4) |
| Twist Express flag: window 10–90% | **0 of 1050 violate** |
| Twist hard rule: homopolymer ≥ 14 bp | **0 of 1050 violate** |
| IDT gBlocks: global 25–75% | **0 of 1050 violate** |
| IDT gBlocks: G/C run ≥ 6 | **0 of 1050 violate** (longest 4) |
| **our rule: window 35–65% absolute** | **277 of 1050 violate** |

**Every published vendor rule we can find accepts every fragment we refuse.** The only
document that rejects them is the one our band was modelled on, and that document is Twist's
codon-optimiser guidance — a different vendor, a different product, and explicitly a design
target rather than an acceptance criterion.

## 4. Recommendation

**Model on IDT, because IDT is who we order from.** Concretely:

1. **Keep `ordering.CURRENT_OPOOL_50PMOL` as the acceptance model.** It is the right product,
   properly sourced, and our fragments satisfy it.
2. **Adopt the Twist Clonal Genes TechNote as the declared *design* guide** — global 25–65%
   plus the 50 bp spread rule. It is a genuine published design document, we pass it
   completely, and the spread rule is a better fit for oligo pools than an absolute band:
   pool synthesis uniformity depends on fragments resembling each other, which is what a
   spread limit controls and an absolute band does not.
3. **Demote the 35–65% absolute window from hard to a declared target.** It is not a vendor
   requirement. Keeping it visible is reasonable engineering conservatism; enforcing it as a
   rejection gate is not supportable by any document found.
4. **Never describe any of this as vendor acceptance.** Both vendors score acceptance with
   proprietary models and say so.

### The counter-argument, which is real

Our fragments are 257–264 nt oligos in a pool, not cloned genes. Oligo-pool synthesis may be
more GC-sensitive than gene synthesis, and neither vendor publishes oPools sequence limits for
us to check. Third-party oligo-pool guides suggest 40–60% as optimal and 30–70% as acceptable;
our fragments sit at 0.565–0.611 global — inside the acceptable range, near the top of the
optimal one.

So a conservative band is defensible **as our own engineering choice**. What is not defensible
is calling it a vendor rule, or keeping it hard while no vendor document supports it.

### What this does not decide

Whether to switch the default. That interacts with the wrapper-aware result
(`PANEL_WRAPPER_AWARE_2026-10-11.md`): with the wrapper-aware constraint on, all 200 designs
satisfy even the strict band, so the band could be kept hard **and** the route made usable.
Loosening the band and fixing the optimiser are alternative solutions to the same problem, and
doing both is not required.

## 5. Sources

| | |
|---|---|
| Twist Clonal Genes TechNote | DOC-001081_TechNote_TwistTipVectorDesign-REV4 **[F]**, local copy |
| Twist gene synthesis resources and FAQs | https://www.twistbioscience.com/products/genes/resources **[F]** |
| Twist design guidelines FAQ | https://www.twistbioscience.com/faq/gene-synthesis/are-there-any-sequence-limitationsdesign-guidelines-genes-which-i-should-follow **[F]** |
| Twist gene synthesis FAQ | https://www.twistbioscience.com/faq/gene-synthesis **[F]** |
| Twist codon optimisation tool FAQ | https://www.twistbioscience.com/faq/using-your-twist-account/what-does-twist-codon-optimization-tool-do **[F]** |
| IDT gBlocks sequence motifs FAQ | https://www.idtdna.com/pages/support/faqs/what-types-of-sequence-motifs-should-be-avoided-when-ordering-gblocks-gene-fragments- **[S]** |
| IDT oPools product page | https://www.idtdna.com/pages/products/custom-dna-rna/dna-oligos/custom-dna-oligos/opools-oligo-pools **[F]**, read 2026-09-15 for the ordering profile |

All pages read 11 October 2026 unless stated. Vendor pages change; re-read before publication,
and confirm IDT's numbers in its order-entry complexity checker, which is the authoritative
accept/reject tool.
