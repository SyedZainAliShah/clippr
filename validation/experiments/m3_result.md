# M3 — fixed-interface reusable inventory: result

**Verdict: the pilot passes all five gate criteria. The limited user-facing route is
authorised.**

The fifth criterion — release preparation complete enough to protect the internal freeze —
was outstanding when the first four were settled, and this document said so. It was closed
later in the same pass: `validation/release_check.py` now reports 10 of 10, with the M1
product exports and the manifest integration that criterion was waiting on. The verdict here
was revised on that evidence, not on a reread of the pilot.

The objective was declared before the pilot ran: recode the kit modules for one host,
maximising codon adaptation, without changing any protein or any four-base interface.

## 1. What was built

The 21 modules the six reference targets need, recoded for *Chlamydomonas reinhardtii*
nuclear expression under genetic code 1, with both interfaces frozen and the protein fixed.

Source inventory: user-supplied Table S1, sha256
`6324ede3e3863ced60f12cc822ebe68f530e9ebc6a5ebc45c490a8816930ef35`. (An earlier draft
quoted `1e0a3c74…` here, which is the package *source fingerprint*, not this file's hash.) The other 21 modules in the
kit were not touched, no host search was performed, and no interface was redesigned.

The pilot rests on one measured fact, not an assumption: across all 200 corpus targets **no
module ever appears at more than one reading frame** (4 modules at frame 0, 28 at frame 1, 8
at frame 2; an earlier draft transposed the last two counts, though the computation behind
them was always the one stated here). That is a property of the kit's architecture, not of codon arithmetic — only 18
of 42 module contributions are codon multiples. Because each module has exactly one frame,
and because the partial codons at each end always fall inside the four-base interface,
freezing the interfaces also freezes every boundary-spanning base. A neighbour's codons
cannot be disturbed.

## 2. Exact before/after metrics

| | min | median | max |
|---|---:|---:|---:|
| CAI before | 0.1784 | 0.2206 | 0.3033 |
| CAI after | 0.5076 | 0.6450 | 0.7199 |
| gain | +0.3009 | +0.4203 | +0.4956 |

**21 of 21 modules improved.** 17 of the 21 are reused by more than one of the six targets,
and all 17 improved; the smallest gain among the reused modules is +0.3009. This is not one
trivial change reported as broad host adaptation — the distribution is tight and every module
moved.

Module-level synthesis QC: **21/21 PASS, unchanged before and after.** No mandatory
constraint regressed.

## 3. Every module invariant, recomputed rather than trusted

The pilot calls the production optimiser and records whatever `constraints_ok` it returned —
the optimiser grading its own work. `m3_independent_qc.py` recomputes each promise from the
original and recoded DNA, importing neither `clippr.qc` nor `clippr.codons`, and deriving the
reading frames from the kit layout rather than reading them from the pilot.

**21 of 21 intended modules pass every invariant.** All 21 required by the six targets are in
the denominator; none was skipped.

| checked independently | result |
|---|---|
| length unchanged | 21/21 |
| both four-base interfaces byte-identical to the deposited original | 21/21 |
| in-frame protein identical to the original's | 21/21, no stop codons |
| no forbidden site (BsaI, BbsI, SapI, either strand) | 21/21 |
| GC inside 0.35–0.65 in every 50 nt window | 21/21 |
| no homopolymer run over 4 | 21/21 |
| CAI reproduced within 5×10⁻⁴ of the pilot's figure | 21/21 |
| reading frame derived here agrees with the pilot's | 21/21 |

Improvement-direction disagreements with the pilot: **0**.

## 4. Assembly and reuse, checked independently

`m3_validate_inventory.py` writes its own join rather than calling `products.reconstruct`, so
a fault in the production reconstructor cannot hide a fault in the inventory.

- **6/6 products assemble**, every four-base overlap verified, zero failures.
- Assembled proteins identical to those from the original inventory; no internal stop codons.
- **No junction-spanning forbidden site** — BsaI, BbsI or SapI, either strand — appears in
  any product. This is the check per-module optimisation cannot do, because neither
  optimiser ever sees its neighbour's sequence.
- One versioned sequence per module identity; every product referencing a module uses the
  same record.
- **Interfaces identical to the deposited originals on 21/21.** Compatible overlaps are not
  the promise: two modules could agree perfectly with each other and both have drifted from
  the kit, and an overlap check alone would not notice.

Three corruption probes, each judged by the same function the real inventory goes through —
a probe scored by its own private test would say nothing about the gate:

| probe | result |
|---|---|
| external entry base altered on the first module | rejected |
| junction overhang mutated on the second module | rejected |
| BsaI site planted inside a module | rejected |

## 5. What this cannot do, measured rather than asserted

The longest exact DNA tract shared between two kit products went **326 nt → 317 nt**. That is
not a diversification result and was never expected to be one. Products that reuse the same
module records contain the same DNA; recoding a shared module changes every product
containing it identically. The 9 nt that did move came from modules that are not shared
between the two worst-case products.

M3 is a reusable-inventory capability. It is not a route to library diversification, and the
two must never be promised together.

A product-level repeat count is likewise not a defect in a module: reuse necessarily repeats
module DNA. Module synthesis constraints and product repeat measurements are reported apart
for that reason.

## 6. Reference comparison (E6)

Run against the pinned local checkout, **v0.1.15 at commit `8882759a`**, in its own
environment — not the v0.1.5 wheel that every earlier oracle comparison in this project used.
The two are different references and the older one does not contain this workflow.

The reference was given our Chlamydomonas codon table in its own `codon,aa,frequency` format,
so both systems optimised toward the same host. Module identities join directly once naming
is normalised (`pPR-1_1A_5N_AATG` ↔ `1A_5N_AATG_v1`).

**The spans are not defined the same way.** The reference's CDS runs 0–2 codons longer,
because it includes the codon that spans a junction and is completed by bases the neighbouring
module supplies; ours stops at the insert's own last whole codon. Our protein is a contiguous
substring of theirs on **21/21** modules, so every score below is taken over exactly that
common region, and both encode the same protein there on 21/21.

**The reference does not reproduce under its declared seed.** Across five runs on
byte-identical inputs at `seed=42`, **21 of 21** shared modules received more than one
distinct sequence. That is the predicate the harness actually tests — it is not the stronger
claim that all five runs differ pairwise on every module, which was not checked. Nothing here
establishes *why* the seed does not fix the output; the per-run inputs and environment would
have to be preserved and compared before attributing a cause.

Our pilot is deterministic, which is asserted here only because it was measured: the pilot
was rerun under the same configuration and **all 21 modules' sequence hashes matched**, with
both runs preserved as `work/m3/inventory_pilot_run1.json` and `_run2.json`.

Scoring one draw of a stochastic optimiser against a deterministic result would report
sampling noise as a difference between the systems, so the comparison runs the reference five
times and reports the spread. (An earlier version of this document quoted a single draw —
0.6333, split 12–9 — which was not a sound measurement.)

CAI under our table, **one scorer applied to both** — neither system's self-reported score is
used, because their `codon_score` and our `cai_after` are not the same statistic:

| | mean CAI on the common region |
|---|---:|
| original deposited inventory | 0.2246 |
| **ours** (deterministic) | **0.6432** |
| **reference**, mean of 5 runs | **0.6394** |
| reference, per-run means | 0.6423, 0.6388, 0.6331, 0.6405, 0.6421 |

Per module, against the reference's per-module mean: reference higher on 7, ours higher on
14, none tied.

**Five runs do not pin the reference's mean.** An earlier five-run batch, on identical
inputs, gave 0.6324 with an 18-3 split; this one gives 0.6394 with 14-7. The batch mean moved
by 0.007 — larger than the 0.004 that separates it from our value here. Any statement that
depends on that gap is therefore not supported by five runs, and no larger batch was run.
The individual run inventories are preserved under `work/m3/reference/runs/` so the spread
can be rechecked rather than taken on trust.

**Reading.** Both systems substantially improve this CAI metric over the deposited inventory
under their respective constraints. **This comparison does not rank the optimisers.**

Our mean is above the reference's in both batches, and above it on a majority of modules in
both. That is as far as the evidence goes. No equivalence or superiority criterion was
declared in advance; the margin (0.004 in this batch, 0.011 in the previous one) is smaller
than the movement between batches; the reference's best single run reached 0.6423, within
0.001 of ours; and the constraint sets differ in ways that push in both directions. A
difference that small, without a pre-declared criterion and without a batch size that
stabilises it, establishes neither "equivalent" nor "better", and this document claims
neither.

The reproducibility difference is the more substantive observation, and it is reported as
what was observed in these runs, not as a judgement about the reference's design.

Runtime is reported but is not a claim: 42 modules in 35.1 s mean for the reference (0.84 s
each) against 21 modules in 1.6 s for us (0.08 s each). These are different amounts of work
on different inputs under different constraints, and the reference's own mean moved from
22.7 s to 35.1 s between batches on this machine, so the figure carries machine load as well.

Declared mismatches, stated rather than aligned:

| | reference | ours |
|---|---|---|
| GC band | 0.25–0.65 global, 0.15–0.85 windowed | 0.35–0.65 over 50 nt |
| max homopolymer | 3 | 4 |
| repeat *k* | 16 | 20 |
| scope | all 42 modules, from amino acids plus a coding mask | the 21 needed, from DNA, interfaces frozen |

A looser GC band and a stricter homopolymer limit push in opposite directions, so neither
system is simply the more constrained one.

What *matched* without adjustment: the forbidden-site set (SapI, BsaI, BpiI — BpiI being
BbsI's isoschizomer), genetic code 1, and the Pryor 2020 ligation data. Both inventories hold
`AATG` and `AGGT` variants, so `ppr_5prime_fusion_site` selects an assembly variant and does
not restrict the module-level comparison.

The reference reported QC PASS on 16/21 of the shared modules on its own terms, with hard
constraints satisfied on 21/21. Its own log records four modules where it took a codon below
its configured minimum relative adaptiveness to deplete a SapI site. **These are its observed
results under these inputs. They are not evidence about the reference generally, and they are
not evidence for us.**

One failure in this comparison was ours: the first attempt died on a `UnicodeEncodeError`
because the harness let a Windows subprocess inherit cp1252 while the reference printed an
arrow. That is recorded here as a harness fault and is not counted as a reference failure.

## 7. Cost

Pilot optimisation loop: **1.7 s of a 600 s cap**, 4 candidate seeds per module. That figure
is the loop only -- it excludes inventory loading, codon-table preparation, context
construction and serialisation, so it is not an end-to-end workflow runtime. 17 of the 21
modules are used by more than one target; that is a reuse count, not a planning-cache
measurement, and no cache instrumentation was recorded for this run.

Reference comparison: 5 runs averaging 22.7 s each, against a 3600 s cap per run.

## 8. Gate decision (E5)

| # | Criterion | Verdict |
|---|---|---|
| 1 | All six pilot targets compile with independently validated proteins and interfaces | **met** — 6/6, 0 failures, independent join |
| 2 | At least one reused module strictly improves the objective, no hidden mandatory-constraint regression | **met** — 17/17 reused modules improved; QC 21/21 PASS unchanged |
| 3 | Amount and distribution of improvement reported | **met** — §2 |
| 4 | Pilot completes within its cap | **met** — 1.7 s of 600 s |
| 5 | Reproducibility, M1 and release preparation complete enough to protect the internal freeze | **met** — `release_check.py` 10/10; see below |

E5 permits a limited user-facing route only if every criterion holds. All five now do, so
**the route is authorised.**

Criterion 5 closed on these: the 627-test suite passing with no skips unaccounted for; the
public example installing into a throwaway virtual environment and reproducing its documented
output verbatim with a fresh cache; the pinned 200-design corpus covering 9S/14S/19S with
zero failures; V5, V6 and V7 agreeing with the exports independently; and the M1 FASTA and
GenBank product exports plus the manifest now shipping inside the result package.

Two things it did **not** close, and which are not release gates: V1's matched-condition NEB
Ligase Fidelity Viewer panel, which needs a person at that web tool, and V2 vendor screening,
which was declared optional and subordinate to correctness.

Expansion toward all 42 modules is a separate bounded extension and is not authorised by this
result.

## 9. Commands

    python validation/experiments/m3_inventory_pilot.py
    python validation/experiments/m3_independent_qc.py
    python validation/experiments/m3_validate_inventory.py
    python validation/experiments/m3_reference_comparison.py --repeats 5

Artefacts: `work/m3/inventory_pilot.json`, `work/m3/independent_qc.json`,
`work/m3/inventory_validation.json`, `work/m3/reference/comparison.json`, the five
individual reference inventories under `work/m3/reference/runs/`, and the two repeated pilot
runs `work/m3/inventory_pilot_run1.json` and `_run2.json`.

A rerun of the reference comparison **archives** the previous directory rather than deleting
it; earlier it removed the evidence outright.
