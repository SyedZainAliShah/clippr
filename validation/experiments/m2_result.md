# M2 — joint full-CDS library selection: result

**Verdict: DO NOT KEEP.** The preregistered gate required improvement on at least two of
three panels. It improved one.

Configuration frozen before any result existed: `m2_joint_selection.json`,
sha256 `ca5582455298044a008ce0d0ce939c94b4d795059a4e944a1572c788cd8939d9`.

## What was measured

Primary objective, declared in advance: minimise the worst exact DNA tract shared between
two distinct library members. Flag threshold 50 nt, unchanged from `HR_THRESHOLD_NT`.

| panel | independent | scaffold diversified | joint greedy (registered) | joint exhaustive (diagnostic) |
|---|---|---|---|---|
| 9S | 107 nt, 15/15 flagged | 77 nt, **2/15** | **73 nt**, 15/15 | 66 nt, 15/15 |
| 14S | 120 nt, 15/15 | **48 nt, 0/15** | 65 nt, 15/15 | 65 nt, 15/15 |
| 19S | 92 nt, 15/15 | **56 nt, 2/15** | 65 nt, 15/15 | 65 nt, 15/15 |

The registered method is deterministic greedy. Exhaustive enumeration is reported as a
diagnostic at this panel size and is **not** the production method; at 9S it would also have
improved, at 14S and 19S it would not.

Banks: 18 targets across three panels, four seeds each, **4 of 4 distinct valid candidates
for every target, zero failures, zero validation problems**. Every candidate was checked for
translation, locked junction overhangs, constraint satisfaction and synthesis QC before
being banked.

## Why it lost, measured rather than assumed

An earlier explanation here claimed every candidate for a target shared one scaffold
encoding, so selection could only shuffle the repeat body. That was false. The bank holds
2–4 distinct 69 nt prefixes per target (8, 9 and 11 across the panels), and an assignment
giving all six members *distinct* prefixes exists in every panel.

The real constraint is that distinct is not dissimilar:

    15 distinct prefixes in the bank
    longest shared run between two DISTINCT prefixes   min 59, median 59, max 68 nt
    distinct prefix pairs still sharing >= 50 nt        105 of 105
    the same measure for scaffold_encoding(0..5)       min 6, median 10, max 17 nt

Every pair of members selectable from these banks must therefore share at least 59 nt in its
prefix — distinct prefixes share ≥59, identical ones share 69 — so **no selector restricted
to these banks can beat the 48 nt and 56 nt baselines**. The bank is the binding constraint,
not the selector. The worst residual tract starts at 0/0 in all three panels, consistent
with that.

## Scope of the conclusion

Seed variation did not provide sufficient separation **in this configuration**. An objective
that requires library members to differ should measure or constrain the separation rather
than assume sampling produces it.

This is not a claim that stochastic design never diversifies. A different optimiser,
constraint set or candidate distribution was not tested, and nothing here measures
recombination frequency, expression or any biological outcome — the objective is a sequence
measurement.

## Cost

Bank construction dominated: 19S banking took roughly 25 minutes, against seconds for the
selection itself. Reporting a single per-method runtime would make joint selection look cheap
by hiding where its cost was paid, so the stages are recorded separately in
`work/m2/method_comparison.json`.

The exhaustive diagnostic remains expensive at 19S — 1530.7 s, down from 2000.1 s after the
two-pass repair. That is a 23% reduction rather than the order of magnitude intended: ranking
now uses the precomputed distance table, but a large number of assignments tie on the primary
metric and the tie-break pass profiles all of them. The redundant work was removed; the tied
work was not.

## Independent reproduction

The orchestrator reproduced the greedy and exhaustive values and the flagged-pair counts on
all three panels with a separate exact-substring implementation that never called
`homology.longest_shared` or this selector, and independently reproduced the prefix
distribution. The 9S diversified baseline was regenerated independently; the 14S and 19S
diversified baselines were read from these results rather than regenerated.

## Commands

    python validation/experiments/build_banks.py
    python validation/experiments/compare_methods.py

Artefacts: `work/m2/candidate_bank.jsonl`, `bank_summary.json`, `method_comparison.json`.
