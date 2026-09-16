# The GC band, and the constraint model behind it

*Written 2026-09-16, after the matched oracle benchmark showed the reference implementation
ahead on mean codon adaptation. Nothing in this document has been implemented. It is the case
for a change, with the measurements that support it and the ones that do not.*

---

## 0. The short version

CLIPPR enforces `GC_BAND = (0.35, 0.65)` over every 50-nt window of a synthesised sequence.
That constant is the **single largest determinant of our headline codon-adaptation figure**, it
is **tighter than any published vendor guideline we can find**, and it has **no recorded
source**.

It is also enforced as a *hard* constraint, where every vendor document that mentions such a
band describes a *design guideline*. The difference between those two readings is worth about
0.16 CAI — more than four times the gap that prompted this document.

This is a decision to be taken deliberately, not a defect to be quietly patched. Widening the
band is one line and would put CLIPPR ahead of the reference on this metric immediately. That
is precisely why it should be sourced and published rather than edited.

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

One controlled experiment. Same solver (DNA Chisel, `CodonOptimize(method="use_best_codon")`),
same four seeds, same 12 deposited modules, same locked interfaces. **Only the band changes.**

| band | mean CAI |
|---|---:|
| CLIPPR — 0.35–0.65 per 50 nt, homopolymer ≤ 4 | 0.692434 |
| reference — 0.15–0.85 per 50 nt, homopolymer ≤ 3 | **0.854766** |

The band is worth **≈ 0.16 CAI**. The gap it was invoked to explain is 0.037.

This is corroborated from the other direction: in the matched benchmark, **36 of 42** of the
reference's ordered sequences breach our GC window, every one of them at 0.660. They are buying
codon adaptation with GC we forbid, and the effect size is not marginal.

**The optimiser is not the difference.** Both systems maximise the same quantity — a mean log
relative adaptiveness. Neither is cleverer than the other at it.

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

Two observations, both actionable:

**Our window is tighter than every profile in that table**, including the two labelled
*conservative*, and tighter than the tolerant profile by a factor of four in width. Meanwhile we
enforce **no global band at all** — the opposite emphasis to every vendor document cited, all
of which lead with a global 25–65% or 25–70% range and treat the window as a local-complexity
guideline.

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
- **That 0.35–0.65 is wrong.** It may be a perfectly good engineering choice. What is
  established is that it is *undocumented*, *tighter than every guideline we can cite*, and
  *expensive* — and that nobody currently knows which of those three facts was intended.

## 7. The proposal

In order, smallest defensible change first.

**7.1 — Source the band, or declare it.** Either attach a citation to `GC_BAND`, or relabel it
as a CLIPPR engineering choice and publish the codon adaptation it costs. The cost is now
measurable to six decimal places; there is no excuse for it being implicit.

**7.2 — Add a global band.** Every vendor document leads with global GC. We check only windows,
which is the weaker half of the published guidance.

**7.3 — Separate targets from rejections.** Adopt the three-way split: what we optimise toward,
what we refuse to ship, and what only a vendor can judge. A module breaching a *target* should
be delivered with a recorded warning, not silently returned unchanged.

**7.4 — Make the profile selectable.** IDT and Twist have materially different rules, and we
order from IDT. A user ordering elsewhere is currently held to a constant that describes
neither vendor.

**7.5 — Report the trade, every time.** Whatever band is chosen, the recoding report should
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
