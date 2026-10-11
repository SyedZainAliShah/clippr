# Which vendor should the sequence policy be modelled on?

*Researched 11 October 2026 from primary vendor pages. Supersedes the single unsourced band.*

The GC band refusing all 200 de novo designs matches exactly one published document, and it is
the least applicable one we could have chosen. This records what each vendor actually publishes,
how our designs score against each, and what to model on.

**Verification levels.** **[F]** fetched and read directly. **[S]** quoted by a search engine only.
IDT's site defeats simple fetches with a redirect loop; its pages here were read in a real
browser instead, so every IDT figure below is **[F]**. Confirm current numbers in IDT's own
order-entry complexity checker, which it names as the authoritative accept/reject tool.

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
| **IDT gBlocks** sequence motifs FAQ | 25–75% | **none published** | 10+ A/T, 6+ G/C | guidance **[F]** |
| **IDT Codon Optimization Tool** | **none published** | **none published** | **none published** | algorithm + checker **[F]** |
| **IDT oPools** (our product) | **none published** | **none published** | **none published** | — **[F]** |

Three things to notice.

**Twist's only stated hard rules** are a 14 bp homopolymer and no CcdB, and its *acceptance*
threshold for a 50 bp window — 10–90% — is roughly five times wider than the band we enforce
as a rejection gate.

**IDT publishes no GC window rule anywhere.** Not for oPools, not for gBlocks, not for its
Codon Optimization Tool, not in its DECODED explainer. Its only numeric GC guidance is a
global 25–75% for gBlocks, a different product. There is no IDT equivalent of our 35–65% over
50 nt — the rule is ours alone. Verbatim, from the gBlocks sequence-motifs FAQ:

> "Sequences that cause problems for gBlocks Gene Fragment synthesis are those having
> extremely low or high GC content (less than 25% and greater than 75%), homopolymeric runs
> of 10 or more As and Ts or 6 or more Gs and Cs. Other structural motifs such as repeats or
> hairpins may also influence our ability to construct gBlocks Gene Fragments. The decision
> to accept or reject an order for a particular gBlocks Gene Fragment sequence is complex and
> multifactorial."

**IDT's own explanation says GC is a proxy, not the thing itself.** Its DECODED article on
complexity screening states that *"Regions with extreme GC content may have higher propensity
to form stable secondary structures,"* and describes screening as secondary-structure
prediction, GC analysis and melting-temperature calculation together. The quantity that
actually obstructs synthesis is **structure**; GC is a cheap stand-in for it.

We check the stand-in and not the thing. PPROPT has no hairpin, self-dimer or cross-dimer
check of any kind, which `ordering.CURRENT_OPOOL_50PMOL.unresolved_rules` already admits:
*"secondary-structure and synthesis-difficulty judgements are the vendor's."*

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

## 3a. Verified against IDT's live pages, 11 October 2026

Read while signed in, so these are the pages as a customer sees them.

**Our oPools profile is exactly right.** IDT's specification table:

| Scale (pmol/oligo) | Oligos per pool | Oligo length |
|---|---|---|
| 1 | 100–20,000 | 40–350 |
| 10 | 10–2,000 | 40–350 |
| **50** | **2–384** | **40–350** |

`ordering.CURRENT_OPOOL_50PMOL` carries 40–350 nt and 2–384 oligos at 50 pmol. Transcribed
2026-09-15, still accurate today, and the scale-dependent count ranges its `unresolved_rules`
warned about are now confirmed rather than assumed.

**Our fragments need oPools specifically.** At 257–264 nt they exceed Ultramer's 200-base
ceiling, so the product choice is forced, not preferred.

**The list prices are confirmed and now dated.** €109.00 for the first 3,300 bases, then
€0.038, €0.025 and €0.013 per base in tiers, and €1.63 per oligo for 5′ phosphorylation —
matching the reference implementation's transcription exactly.
`docs/reference_implementation_findings.md` §3.8 observes that PPROPT refuses to price while
the reference's figures *"are undated too — but they are there, and ours are not."* That
objection is now answerable: these figures carry a date.

**oPools has no published GC, complexity or structure rule**, and no QC: the page lists
*"QC and quantification: None"*, an average dropout of 0.8% across 900,000 oligos, and notes
dropout *"can vary by length and sequence complexity."* Complexity affects yield, not
acceptance.

## 3b. Secondary structure, measured for the first time

Folded all 1,050 ordered fragments in both arms with ViennaRNA 2.7.2 under the Mathews 2004
DNA parameters — `validation/experiments/secondary_structure_scan.py`.

| | baseline | wrapper-aware |
|---|---|---|
| first 48 nt ΔG, median | −5.0 | −5.1 |
| first 48 nt ΔG, worst | −10.7 | −10.7 |
| whole fragment ΔG, median | −44.7 | −44.9 |
| **breaching Twist's −8 in the first 48 nt** | **14 of 1050** | **16 of 1050** |

**The GC fix is structure-neutral.** Medians moved 0.1–0.2 kcal/mol, which is noise. The
wrapper-aware work did not trade secondary structure for GC — a risk nobody had checked.

**But roughly 1.4% of fragments carry a hairpin stronger than −8 kcal/mol in the first 48 nt,
in the shipped corpus, today.** That is the only published structure threshold found anywhere,
and we have been shipping fragments past it without looking.

Four reasons this is a flag and not a failure: ViennaRNA is not IDT's implementation; Twist's
−8 rule belongs to its codon optimiser for gene synthesis, not to oPools; a folding prediction
is not a measurement; and IDT publishes no structure rule for oPools at all, treating
complexity as a yield question. It is nonetheless a quantified gap where there was previously
no number.

Note that neither IDT's OligoAnalyzer nor primer3 can evaluate these fragments — both cap
thermodynamic calculations at 60 bases, against our 257–264 nt. A folding engine is the only
instrument that applies.

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
   requirement — no IDT document contains it, and the one Twist document that does frames it
   as optimisation guidance. Keeping it visible is reasonable engineering conservatism;
   enforcing it as a rejection gate is not supportable by any document found.
4. **Add a secondary-structure check, which is the real gap.** Both vendors screen structure;
   neither publishes a GC window rule we match. Tuning the proxy while ignoring the quantity
   it proxies for is the wrong priority. IDT's OligoAnalyzer computes hairpin and dimer ΔG and
   accepts our full 257 nt fragments, and IDT exposes an API, so this can be scripted rather
   than left as a manual step.
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
| IDT gBlocks sequence motifs FAQ | https://eu.idtdna.com/pages/support/faqs/what-types-of-sequence-motifs-should-be-avoided-when-ordering-gblocks-gene-fragments- **[F]** |
| IDT Codon Optimization Tool | https://eu.idtdna.com/pages/tools/codon-optimization-tool **[F]** — no numeric thresholds published; the tool is behind sign-in |
| IDT DECODED, codon optimisation basics | https://eu.idtdna.com/page/support-and-education/decoded-plus/codon-optimization-the-basics-explained/ **[F]** |
| IDT oPools product page | https://www.idtdna.com/pages/products/custom-dna-rna/dna-oligos/custom-dna-oligos/opools-oligo-pools **[F]**, read 2026-09-15 for the ordering profile |

All pages read 11 October 2026 unless stated. Vendor pages change; re-read before publication,
and confirm IDT's numbers in its order-entry complexity checker, which is the authoritative
accept/reject tool.
