# The GC band, and the constraint model behind it

*Written 2026-09-16, after the matched oracle benchmark showed the reference implementation
ahead on mean codon adaptation. The case for a change, with the measurements that support it
and the ones that do not.*

**Implemented 2026-09-16:** the profile mechanism this document argues for now exists — see
`docs/synthesis_policy.md`. The **default is unchanged**, so every number here still holds.
What the implementation settles is that the policy is one object reaching both the solver and
the validator; what it deliberately does **not** settle is which profile should be the default.

---

## 0. The short version

CLIPPR enforces `GC_BAND = (0.35, 0.65)` over every 50-nt window of a synthesised sequence.
That constant is the **single largest determinant of our headline codon-adaptation figure** —
worth about 0.16 CAI — and it has **no recorded source in this repository**.

**Correction, 2026-09-16.** An earlier version of this document said the band was "tighter than
any published vendor guideline we can find." That is **false**, and the error was mine. Twist's
own gene-synthesis guidance states verbatim: *"We avoid fitted sequences that create global GC%
of less than 25% or more than 65% and local GC windows (50 bp) of less than 35% or more than
65%."* That is our band and our window, exactly. I had compared against the reference
implementation's *transcription* of Twist's rules rather than reading Twist's page, having
written two paragraphs earlier in this same document warning that the transcription should not
be trusted for exactly this purpose.

What survives the correction, and it is still the point:

- **The band has no recorded source here.** Matching a published guideline after the fact does
  not tell us where the constant came from, and Twist must not be retro-fitted as its origin.
- **It is a Twist *codon-optimisation* guideline, and we order from IDT.** A vendor name is not
  a policy; the product is. Twist's own page separates this guidance from its sequence
  acceptance and complexity rules.
- **We enforce the local half and not the global half.** Twist pairs the 50-bp window with a
  global 25–65% band. CLIPPR has no global band at all — so we are not simply "stricter", we
  are *differently* scoped, which is harder to defend than either.
- **We enforce as hard what Twist frames as optimisation guidance.** That reading is worth the
  0.16 CAI on its own.

This is a decision to take deliberately, not a defect to patch. It is also **not** a one-line
change — see §7.1.

---

## 1. What prompted it

The matched oracle benchmark scores both systems' delivered sequences with one scorer, over the
translated span they share:

| | CLIPPR | reference |
|---|---:|---:|
| mean CAI, aligned spans | 0.621004 | **0.658672** |
| modules ahead | **22** | 20 |

The reference is ahead on the mean while CLIPPR wins the per-module split, so their advantage
is concentrated rather than uniform. The question this document answers is *why*, and whether
it is something to copy.

## 2. The measurement

**My first version of this experiment changed two variables at once** — the GC band *and* the
homopolymer cap — and reported the result as "only the band changes". It was not a controlled
experiment. The orchestrator separated them in a 2×2 over all 42 modules, same solver, same
four seeds, same locked interfaces:

| local GC band | homopolymer cap | first 12, best coding CAI | all 42, best coding CAI | best coding choices failing the full-substrate checks |
|---|---:|---:|---:|---:|
| 35–65% | 4 | 0.692434 | 0.657195 | 3 / 42 |
| 35–65% | 3 | 0.688232 | 0.653145 | 11 / 42 |
| 15–85% | 4 | **0.854766** | **0.819921** | **0 / 42** |
| 15–85% | 3 | 0.854766 | 0.818415 | 8 / 42 |

**The GC effect survives the correction.** Holding the cap at 4, the first-12 figures reproduce
exactly: 0.692434 → 0.854766. The homopolymer cap contributes almost nothing by comparison.

Three things the wider experiment shows that mine could not:

- **Coding-span CAI is not a deliverable.** Scoring the best coding candidate is not the same as
  producing an inventory that satisfies the whole synthesis contract. Under full-substrate
  filtering the strict cap-4 run returns **0.618830** — reproducing today's shipped recoding
  result — while the broad cap-4 run returns **0.819921** with all 42 selected substrates
  passing that experiment's own band, homopolymer and site checks.
- **The 12-module cohort hid a scope problem.** In the cap-3 columns, eight wide-band coding
  winners fail the full-substrate homopolymer target. A clean first-12 result would have missed
  every one.
- **This varies hard bounds against hard bounds.** It says nothing about soft penalties, and it
  does not isolate every difference between the two systems.

Corroboration from the other direction: in the matched benchmark, **36 of 42** of the
reference's ordered sequences breach our GC window. (Each is reported at 0.660 because the
scanner stops at the first offending window — that is a violation, not that sequence's worst
window.)

**The optimiser is not the difference** — on this axis. Both maximise a mean log relative
adaptiveness. But that is a shared *component*, not a shared objective: theirs also carries
weighted GC, homopolymer, repeat and library-similarity penalties in the same scalar. Calling
the two optimisers equivalent would go beyond what was measured.

## 3. Where our band came from

    src/clippr/recoding.py:33
    GC_BAND = (0.35, 0.65)
    GC_WINDOW = 50
    MAX_HOMOPOLYMER = 4

The comment above it is careful and correct about **scope** — it explains that the band applies
to the whole module rather than the optimiser's sub-span, which was a real defect once. It says
nothing about **provenance**. There is no citation anywhere in the package.

Nor does the vendor profile it presumably serves. `CURRENT_OPOOL_50PMOL` — IDT oPools, 50 pmol,
transcribed from the published specification — constrains:

- oligo length 40–350 nt
- 2–384 oligos per pool
- one scale

and states **no GC rule at all**. Its `unresolved_rules` explicitly defer "secondary-structure
and synthesis-difficulty judgements" to the vendor.

So the constant enforces a rule that our own vendor model does not contain.

## 4. What the reference does instead, and what is worth taking

The reference carries **five named vendor profiles**, each with a product URL, and selects one:

| profile | global GC | window | homopolymer target | machine-hard |
|---|---|---|---|---|
| Twist · Express / Low complexity | 0.25–0.65 | 0.20–0.80 / 50 nt | 3 | ≤ 13 |
| Twist · Standard guidelines | 0.25–0.65 | 0.15–0.85 / 50 nt | 3 | ≤ 13 |
| Twist · Complex Genes tolerant | 0.25–0.65 | 0.10–0.90 / 50 nt | 3 | ≤ 30 |
| **IDT · gBlocks / eBlocks conservative** | 0.25–0.70 | 0.20–0.80 / 40 nt | 3 | — |
| Generic · conservative (default) | 0.30–… | 0.25–0.75 / 50 nt | 3 | — |
| **CLIPPR** | **none** | **0.35–0.65 / 50 nt** | 4 | — |

**Read that table with care — it is the reference's transcription, not the vendors' pages.**
Trusting it is how this document originally reached a false conclusion. Twist's own page gives
**35–65% over 50 bp** for codon optimisation, which no row above reproduces, and which is
exactly CLIPPR's band.

What the primary pages support, product by product:

| primary source | supports | does **not** support |
|---|---|---|
| Twist gene resources (read 2026-09-16) | global 25–65%, local 35–65% over 50 bp, as **codon-optimisation** guidance | treating it as an acceptance rule, or as CLIPPR's historical source |
| Twist Complex Genes | 50-bp GC outside 10–90% among high-complexity criteria | universal acceptance bounds |
| IDT gBlocks FAQ *(orchestrator's reading; the page redirect-loops for me)* | GC below 25% or above 75% can cause synthesis problems; acceptance is multifactorial | applying gBlocks rules to oPools |
| IDT oPools product page | length and scale-dependent pool counts; **no numeric GC threshold** | that any sequence will be accepted |

The actionable point is narrower than the one I made, and better: **we enforce the local half of
Twist's pair and not the global half**, we enforce it as hard where Twist frames it as
optimisation guidance, and we order from IDT, whose oPools page states no GC rule at all.

**They separate three things we collapse into one:**

| their field | meaning | ours |
|---|---|---|
| `synthesis` | design targets, optimised against as soft penalties | enforced as hard |
| `machine_hard_constraints` | what is actually rejected — e.g. Twist's real homopolymer rule is < 14, not ≤ 3 | — |
| `manual_review_rules` | human or vendor judgement, e.g. "no CcdB" | — |

We treat a design target as a rejection criterion. That single conflation is the measured gap.

## 5. The deeper difference: soft penalties versus hard constraints

Their objective is one scalar, minimised by simulated annealing:

    score = w_codon·codon_score
          − w_global_gc·gc_penalty²      − w_local_gc·local_gc_penalty
          − w_homopolymer·excess²        − w_internal_repeat·repeat_penalty
          − w_library_similarity·shared_kmer_penalty

Only a forbidden restriction site is an absolute veto. Everything else is a **cost**, and the
weights are configuration.

Ours hands GC, homopolymer and enzyme sites to DNA Chisel as **constraints**, with codon
adaptation as the objective. A candidate one window-percentage-point over the band is discarded
outright.

The consequence is not cosmetic:

- **Theirs trades.** It will accept a hot window for a large codon gain, then mark the result
  `qc_status: WARNING` and deliver it. The user sees the compromise.
- **Ours refuses.** If no candidate satisfies the band, the module is returned *unchanged* — its
  original, poorly adapted sequence. Our "unchanged" count is not an absence of opportunity; it
  is the band vetoing every improvement.

That is a second, independent contribution to the CAI difference, on top of the band width
itself, and it is invisible in any single number.

## 6. What this does **not** establish

- **That a wider band is safe.** No sequence from either system has been synthesised. Twist and
  IDT both score submissions with proprietary models; a published guideline is not an
  acceptance guarantee, and neither is a looser one a licence.
- **That the reference's output would manufacture.** Their own QC flags many of these as
  `WARNING`, and their profile notes say final acceptance is ML-scored by the vendor.
- **That CAI predicts expression.** It does not, and no part of this document should be read as
  saying a higher number is a better protein. `docs/objectives.md` §8 states the limit.
- **That 0.35–0.65 is wrong.** It matches Twist's published codon-optimisation guidance
  exactly. What is established is that it is *undocumented here*, *enforced as hard where its
  matching guideline is advisory*, *applied without the global band that guideline pairs it
  with*, and *expensive*.
- **That the two optimisers are equivalent.** They share a codon component, not a complete
  objective.
- **That a soft-penalty model would do better.** This experiment varies hard bounds against hard
  bounds. Hard versus soft at the *same* thresholds is a separate experiment and has not been
  run.

## 7. The proposal

In order, smallest defensible change first.

**7.1 — It is not a one-line change, and I said it was.** `recoding.GC_BAND` is the *validator's*
band; `codons.optimize_cds` carries its own `gc_bounds=(0.35, 0.65)` default and
`recode_inventory` never passes one. Patching the constant alone moves the four-seed result to
**0.657195**, not the fully broadened **0.819921** — it changes which candidates the final check
accepts while the solver still searches the narrow space.

That is the real finding under the one I claimed: **the solver and the validator hold the same
threshold in two places and neither knows about the other.** A band is not a constant to widen;
it is a policy that has to reach every stage that judges a sequence.

(One further care: in that probe all three previously-unchanged modules improved, which supports
a constraint-veto reading *for those modules and those attempts*. "Every unchanged module means
the band vetoed every possible improvement" is a wider claim and was not tested.)

**7.2 — Source the band, or declare it.** Either attach a citation to `GC_BAND`, or relabel it
as a CLIPPR engineering choice and publish the codon adaptation it costs. Matching Twist's
published guidance is not the same as having been derived from it, and the repository still
records no source.

**7.3 — Add a global band.** Every vendor document leads with global GC. We check only windows,
which is the weaker half of the published guidance.

**7.4 — Separate targets from rejections.** Adopt the three-way split: what we optimise toward,
what we refuse to ship, and what only a vendor can judge. A module breaching a *target* should
be delivered with a recorded warning, not silently returned unchanged.

**7.5 — Make the profile selectable.** IDT and Twist have materially different rules, and we
order from IDT. A user ordering elsewhere is currently held to a constant that describes
neither vendor.

**7.6 — Report the trade, every time.** Whatever band is chosen, the recoding report should
state the CAI forgone to satisfy it. That converts an invisible ceiling into a declared
decision, which is the whole point of this document.

---

## Appendix — reproducing the measurement

    # the controlled band experiment in §2
    python - <<'PY'
    import sys, json, statistics; sys.path.insert(0, "src")
    from Bio.Seq import Seq
    from clippr import inventories as inv, constants as C
    from clippr.codons import complete_table, optimize_cds
    from clippr.objectives import codon_adaptation
    from clippr.recoding import locked_interface_sites
    table = complete_table(json.load(open("data/codon_tables/kazusa_3055.json")), 1)
    dep = inv.load_deposited("data/grasp_supp/Table S1.xlsx")
    for name, kw in {"ours": dict(gc_bounds=(0.35, 0.65), max_homopolymer=4),
                     "theirs": dict(gc_bounds=(0.15, 0.85), max_homopolymer=3)}.items():
        cais = []
        for m in sorted(dep.modules)[:12]:
            r = dep.modules[m]; s, e = r.coding_interval
            protein = str(Seq(r.dna[s:e]).translate())
            best = codon_adaptation(r.dna[s:e], table)["cai"]
            for seed in (42, 43, 44, 45):
                got = optimize_cds(protein, locked_sites=locked_interface_sites(r.dna, r.frame),
                                   codon_table=table, genetic_code=1, seed=seed,
                                   enzymes=C.DEFAULT_ENZYME_PROFILE, unique_kmer_size=None, **kw)
                if got.get("constraints_ok"):
                    best = max(best, codon_adaptation(got["cds"], table)["cai"])
            cais.append(best)
        print(name, round(statistics.mean(cais), 6))
    PY

Expected: `ours 0.692434`, `theirs 0.854766`.

The vendor-profile table in §4 is read from the reference checkout's
`grasp_library/synthesis_vendors.py`, which carries the product URLs it cites. Those rules are
**their transcription** of Twist and IDT documentation; before any of them is adopted here, the
primary vendor pages should be read directly and cited from source.
