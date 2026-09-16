# Full workflow — handoff for independent verification

Every claim names the command that produces it and the artefact it writes. Labels: **RAN**
(executed here), **READ** (inspected an artefact a previous run wrote), **ESTIMATE** (not
measured).

---

## 1. What is callable

Six tasks, each usable on its own, each writing persistent artefacts. `clippr.workflow`.

| # | Task | Function |
|---|---|---|
| 1 | design complete CDSs from RNA targets | `design_for_synthesis` |
| 2 | load an inventory and compile targets | `load_and_compile` |
| 3 | recode an inventory for a host, interfaces fixed | `recode_for_host` |
| 4 | optimise an inventory as a collection | `optimise_collection` |
| 5 | explore junction/codon trade-offs | `explore_interfaces` |
| 5b | commit to a front candidate, bound to its inventory | `select_interface` |
| 5c | build assembly-ready order substrates | `order_items_for` |
| 6 | check eligibility, plan pools, export | `plan_order` |

```python
from clippr import workflow as w
from clippr.ordering import CURRENT_OPOOL_50PMOL

# seeds=8 reproduces the saved inventory; seeds=4 gives a different one (see 12a)
recoded = w.recode_for_host("Table S1.xlsx", codon_table, "out/recode",
                            label="chlamy", seeds=8, wall_seconds=600)
front = w.explore_interfaces(recoded.artefacts["inventory"], targets,
                             codon_table, "out/interfaces", max_evaluations=100)

# Selection is bound to the inventory the front was measured against, and recomputes
# its objectives from the sequences it actually delivers.
selected = w.select_interface(recoded.artefacts["inventory"],
                              front.artefacts["front"], "out/selected")
recompiled = w.load_and_compile(selected.artefacts["inventory"], targets, "out/compile")

# Assembly-ready substrates, not bare inserts -- each digest-verified before export.
items = w.order_items_for(selected.artefacts["inventory"], "out/items")
order = w.plan_order(json.loads(Path(items.artefacts["order_items"]).read_text()),
                     CURRENT_OPOOL_50PMOL, "out/order")

w.write_package([recoded, front, selected, recompiled, items, order], "out/package")
```

## 2. Contracts, written before the code

- [`docs/workflow_contract.md`](docs/workflow_contract.md) — data contracts, coordinate
  conventions, failure states
- [`docs/objectives.md`](docs/objectives.md) — every objective's exact definition, including
  the awkward cases implementations differ on
- [`docs/methodology_sources.md`](docs/methodology_sources.md) — primary fact vs documented
  method vs engineering choice vs vendor rule

## 3. Gate results — all RAN

| Gate | Command | Result |
|---|---|---|
| A | (folded into the release check) | closed |
| B | `python validation/experiments/phaseb_gate.py` | **PASS** |
| C | `python validation/experiments/phasec_gate.py --max-proposals 168` | **PASS** |
| D | `python validation/experiments/phased_gate.py` | **PASS** |
| E | `python validation/experiments/phasee_gate.py` | **PASS** |
| F | `python validation/experiments/phasef_gate.py` | **PASS** |
| oracle | `python validation/experiments/matched_oracle_benchmark.py` | **PASS** — see §12b |
| level 1 | `python validation/experiments/level1_geometry.py` | **DERIVED AND FALSIFIED** — see §5d |
| §10 | `python validation/release_check.py --clean-install` | **10/10** |
| notebook | `python validation/run_notebook.py` | **PASS** — 21 cells, 0 failed, on the real input path |
| commands | `docs/verification_commands.md` | one invocation per operation |

Artefacts: `work/phase{b,c,d,e,f}/gate_*.json`, `work/oracle/matched_benchmark.json`,
`work/release/release_check.json`.

## 4. Coverage — RAN

**All 42 modules**, not 40. The 200-target corpus exercises 40; the two `AGGT` start variants
were unreachable because `parts._find` returned the first block match and they are identical
to their `AATG` twins in block, specificity residues and length. `select()` is now
variant-aware, and `witness_targets()` supplies a context for every module — **42/42**.

| | |
|---|---|
| modules loaded and validated | 42/42 |
| witness contexts | 42/42 |
| recoded realisations | 42/42 (39 improved, 3 unchanged) |
| products checked | 6/6 across 9S/14S/19S × both fusion sites |
| reload reproduces | version, sequences and recompiled product |

## 4a. What every result below is bound to

A number without its inputs is not reproducible, so the identities are stated once here and
referred to rather than repeated.

| input | identity |
|---|---|
| `data/grasp_supp/Table S1.xlsx` | sha256 `6324ede3e3863ced60f12cc822ebe68f530e9ebc6a5ebc45c490a8816930ef35` |
| `data/codon_tables/kazusa_3055.json` | sha256 `a9d8f6fc88e15694be1473123ff3dc9cb7e501777c462f3d244ec7368fb78ab3` |
| deposited inventory | version `7d550e628605e006` |
| recoded inventory (the shipped one) | version `3abdf2e9db7d805f`, `seeds=8` |
| reference checkout | `grasp-library-designer` v0.1.15, commit `8882759a` |
| block plasmids (§5d) | `9SDYW_level0.gb` sha256 `ec28074a288b3607cbb488d3ec938c35c17c8ef6d0aca499eaccf690c9860ad5` |
| level-1 products (§5d) | `9SrpoaDYW_variants_level1.gb` sha256 `0e91574ff0313eff6ae234ef93528a8f095380d2ffcadd15ed209187a06823c0` |

The last two were fetched 2026-09-15 from `github.com/farleykvdg/GRASP` at commit `b569a8c4`.
That repository declares **no licence**, so they are read for verification and never
redistributed; a derived overhang set is a fact and ships, the records themselves do not.

## 5. Optimiser results — RAN

**Fixed-interface recoding**, 42 modules: mean CAI **0.2079 → 0.6210** (39 improved, 3
unchanged), evaluated on the released fragment — see §5b.
The deposited modules were never optimised for this host, so this measures adding an
optimisation that was absent, not beating one.

**Collection optimisation**, matched budgets (168 proposals each, 3 seeds × both orders,
all 18 runs preserved). **Regenerated after the substrate-scope fix** — the earlier table on
this line was produced while the optimiser judged candidates on the bare insert, and is
superseded:

| mode | sharing Δ median | best | adaptation Δ median | accepted (median) |
|---|---:|---:|---:|---:|
| greedy | −2.591 | −2.774 | −0.00204 | 18.5 |
| anneal | **−2.647** | **−3.164** | −0.00274 | 93 |
| random (control) | −0.578 | −2.121 | −0.00067 | 82 |

Input `3abdf2e9db7d805f`, baseline adaptation 0.621004 and sharing/pair 19.8444.
Zero constraint violations and zero objective disagreements across all 18 runs — previously
16 of 18 returned an inventory the export step refused.

Both real methods still beat the control decisively. **What changed with the fix is the gap
between them**: annealing's advantage over greedy fell from 0.48 sharing points to **0.056**,
because the candidates it had been exploiting were ones the export step would have rejected.
The old claim that annealing pays "~2.6× more adaptation" for its edge no longer holds
either; the ratio is now 1.34×.

**Recommended default: greedy.** It reaches **97.87%** of annealing's median sharing
reduction with a third of the accepted changes and a smaller adaptation cost. (An earlier
draft said "within 2%", which is tighter than the measurement supports.) All three modes are
reproducible for a fixed seed, order and budget; none is seed-free deterministic.

**Joint junction/codon search** — *rebuilt after a reaction-model correction; see §5a.*
Level-0 reactions only, scored with **BbsI-HF**, 100 evaluations.

The **recommended** candidate raises predicted level-0 fidelity from **0.751683** to
**0.788894**, about **+3.72 percentage points**. The **highest fidelity observed** in the
bounded search is **0.794433**, about +4.28 points, at lower codon adaptation. The
recommendation follows the declared policy — retain candidates within 0.01 of the best
observed fidelity, then maximise adaptation — so it trades roughly 0.55 fidelity points for
adaptation. Observed front of 5, not degenerate. Level-1 reactions are now **scored too** — see §5d; level 0 still sets the minimum.

0.788894 is the recommended value, **not the best observed**. No overall assembly-yield claim
follows from either.

## 5a. A reaction-model correction — the headline changed

The orchestrator found that the previous Phase D result described a reaction that does not
happen, and was right.

**What was wrong.** `compile_target` assigned the join *between* two sub-assemblies to the
preceding reaction, and `joint_search` scored every reaction with **BsaI-HFv2**. Both are
wrong for this kit:

- I independently confirmed the cut geometry on all 42 primary `GRASP_-1.gb` records: **every
  one has two BbsI sites and no BsaI site.** A level-0 reaction is a BbsI reaction.
- The cut sequences are `CTCA/ACTC` (A modules), `ACTC/AAGA` (B), `AAGA/GCAC` (C),
  `GCAC/TGAA` (D), `TGAA/CGAG` (E). So every level-0 reaction has the same set:
  **CTCA, ACTC, AAGA, GCAC, TGAA, CGAG**.
- `CTTC`, `GTGA` and `CACG` are **not BbsI junctions at all**. They join assembled blocks,
  at level 1, cut by BsaI. Folding them into level-0 reactions manufactured the "19S third
  reaction bottleneck" at 0.596 that the old headline was built on.

**What changed.** `assembly_spec.ASSEMBLY_STAGES` now declares each stage's enzyme, matrix and
destination. Reactions carry their stage; level-0 reactions are scored with BbsI-HF; level-1
reactions were reported **unscorable** at the time, because their ends were thought to come
from a block plasmid this
package does not hold, and inventing a number for them would be the whole problem again. The
search varies level-0 junction classes only — four of them, exactly the internal junctions the
cut geometry predicts.

**The corrected result.** Baseline **0.751683**, matching the orchestrator's independently
computed value for that set. Recommended **0.788894** (+3.72 points); highest observed
**0.794433** (+4.28 points) at lower adaptation. The previous 0.596 → 0.711 figure is
superseded and should not be quoted.

The fixed-interface recoding and collection results are unaffected — they never depended on
the reaction model.

## 5b. Synthesis constraints are evaluated on what is ordered

A second scope correction, the same shape as §5a's and found the same way.

`pPR-1_19E_LN5N` passed the module contract at **GC 0.640** as a bare insert and reached
**0.660** as the released fragment — an E module's fragment is its insert plus the GC-rich
`CGAG` overhang. The order still reported success, because the constraint was evaluated before
the fragment existed.

**Repaired at source and at export.** `recode_inventory` now evaluates candidates on the
**released fragment**, so the violation is avoided rather than caught; `substrates.build`
additionally checks the wrapped sequence it hands out and records
`synthesis_problems`; and `order_items_for` fails when any item carries one, so a package
cannot report `ok=True` while exporting a sequence outside its own declared band.

The two intended BbsI sites are excluded by role — they are the construction. Everything else
applies to the whole ordered sequence.

**Cost of the tighter scope:** 39 modules improved instead of 40, 3 unchanged instead of 2,
and mean CAI **0.2079 → 0.6210** rather than 0.6342. The shipped inventory now has **zero**
substrates outside the contract.

## 5c. Search context travels with the front

A front records the scoring knobs it was measured under. `select_interface` was validating
every front with the default `k=20`, so a perfectly valid `k=10` search was refused. `k`,
genetic code, matrix and destination are now read from the front's binding and used for both
the check and the recomputation. `allow_mismatch=True` remains available but is not required
for a valid non-default context — suppressing the check was never the right answer.

## 5d. Level-1 block geometry — established; the reaction is still not scored

Two corrections in one section, in opposite directions. The original claim — that a level-1
reaction could not be scored "without block-plasmid and destination context this package does
not hold" — was wrong about the block plasmids. A later version of this section then said the
whole assembly was scored, which was wrong about the reaction. Both are recorded below, since
an intermediate wrong answer is part of how this one was reached.

**The block plasmids were already on disk**, in the provenance directory next to the module
records, fetched 2026-09-15 to verify the level-0 cut geometry and never reopened.
`9SDYW_level0.gb` holds 28 assembled block plasmids, each with exactly two BsaI sites and no
BbsI site — the mirror of the 42 module plasmids.

**And they are not needed — for the geometry.** Digesting all 28 gives 499 nt `AGGT → CTTC`
and 406 nt `CTTC → TTCG`. Joined at the shared overhang: **901 nt, exactly the compiled 9S
length**, and the joined sequence occurs verbatim inside the deposited level-1 product
`9SrpoADYW_pICH47802_lc`. The block vector contributes no bases, so the **block ends** are the
compiled product's own first and last interfaces — which Table S1 has always supplied.

That settles where the blocks are cut. It does **not** settle what else is in the tube with
them, which is what a fidelity needs.

Table S1 agrees on every interface, and was consulted only after the derivation was finished:
`1E → CTTC`, `14E → GTGA`, `19E → CACG`, `2E → TTCG`, `1A` entering at `AATG` or `AGGT`.

| reaction | participants | enzyme / matrix | figure |
|---|---|---|---:|
| level 0 (module inserts → block) | **complete** | BbsI-HF | **0.751683** — a fidelity |
| level 1 block subset, 9S | incomplete | BsaI-HFv2 | 0.998225 — a diagnostic |
| level 1 block subset, 14S | incomplete | BsaI-HFv2 | 0.998225 — a diagnostic |
| level 1 block subset, 19S | incomplete | BsaI-HFv2 | 0.996046 — a diagnostic |
| 19S subset **plus the known omitted ends** | still incomplete | BsaI-HFv2 | **0.742067** |

Only the first row is a reaction fidelity. The rest are scores of an overhang set that is not,
on its own, any tube's contents.

**No conclusion about which stage is the bottleneck follows from this, and the earlier claim
that one did is withdrawn.** Ligation fidelity is a property of the whole competing overhang
set in a tube, and the level-1 set here is incomplete. Adding only the ends we *know* are
omitted moves the 19S figure from 0.996046 to **0.742067** — below level 0's 0.751683. So the
one sensitivity calculation available points the opposite way to the claim it was used to
support, and the true participant list could move it further in either direction.

What the numbers do support: the level-0 reaction is scored at 0.751683 over a **complete**
overhang set, and the level-1 block subset would score 0.996–0.998 **if it were the whole
reaction, which it is not**. Those two are not comparable and are not compared here.

**How far the falsifier reaches.** Only the **9S** block plasmids are deposited, so 9S is
verified two independent ways — by digesting real plasmids and by reading Table S1 — and they
agree. The 14S and 19S sets follow from the same module interfaces in Table S1, which the 9S
case shows to be the right place to read them, but there are no 14S/19S block plasmids to
digest. They are one falsifier short of 9S, and that is stated rather than glossed.

**What a level-1 number does not cover.** The deposited construct co-assembles a P2L2S2 linker
(`TGTG → CAAC`) and the consensus DYW domain (`CAAC → GCTT`) in the same BsaI reaction, plus
at least one further part bridging `TTCG` to `TGTG` that is absent from the supplied records.
The figure above is for the **PPR block junctions only**; adding an editing domain adds
overhangs it does not include.

**The structured status, which is what a reader actually believes.** A caveat in prose beside
a `context_available: True` flag is not a caveat; the flag wins. So the two statuses are now
separate fields: `geometry_established: True` and `participants_established: False`.
`reaction_overhangs` is `null`, `scorable` is `False`, level 1 appears in
`unscorable_reactions`, and the PPR-only set lives in `block_subset_overhangs` labelled as a
diagnostic. Every scored reaction still records the set it was scored on and where its ends
came from (`ends_from`).

Derivation and both falsifiers: `validation/experiments/level1_geometry.py`, artefact
`work/level1/geometry.json`.

## 6. What the numbers do not mean

- **CAI is a predictor, not a measurement.** A 3× CAI gain is not a 3× expression gain.
- **Fidelity is a prediction from 2020 data, not a measured yield.** §7a shows CLIPPR, the
  GRASP reference and NEB's own viewer agreeing on the same sets, so the number is the one
  those tools compute — but all three predict from the same dataset, and none of them
  measured an assembly. The multi-reaction summary is the **minimum** across reactions,
  never a product.
- **The front is observed, not optimal.** Nondominated among candidates actually evaluated; a
  bounded beam search does not enumerate the feasible space.
- **The k-mer collection penalty is an engineering objective**, not a recombination
  probability. Lowering it does not necessarily lower the worst shared tract.
- **Nothing is biologically validated.** No binding, expression or assembly was measured.
- **No order was placed**, and local eligibility is not vendor approval.

## 7. Vendor support — RAN, with its limits stated

Prices carry a date or they are not offered. `HISTORICAL_OPOOL` holds this project's
long-standing 109.00 EUR / 1.63 EUR figures as **historical**: `price_status` blocks them from
surfacing as current and `priced_on` is absent, so `prices_usable` is False. Eligibility,
pooling and export work without any price.

Pooling matches the exhaustive optimum on every case n=2–10 with **zero gap** and refuses
nothing a valid partition covers — but that is under a flat per-pool tier, and it is not a
general optimality claim. Price arithmetic agrees with 5 hand-worked fixtures.

**Two profiles, because rule provenance and price provenance expire differently.**
`CURRENT_OPOOL_50PMOL` carries the published rules — 40–350 nt per oligo, 2–384 oligos per
pool at 50 pmol — with **no prices**, so cost reports unavailable while eligibility, pooling
and export work. `HISTORICAL_OPOOL` keeps the older figures, labelled historical.

This matters: the historical profile's 20 nt minimum **passed two 30-base oligos that the
current published minimum of 40 rejects** — a false positive in the direction that would reach
a vendor. Boundary cases at 39/40/350/351 now hold.

**Still not done:** the current rules were transcribed on 2026-09-15 rather than queried live,
and only the 50 pmol scale is modelled. Confirm against the current page before ordering.

## 7a. NEB Ligase Fidelity Viewer cross-check — RAN (by hand, 2026-09-15)

**CLIPPR and the reference agree on four specified BsaI cycling test sets. The recorded NEB
integer percentages are consistent with these scores under truncation. NEB's documented count
normalisation explains the expected scale factor.**

That sentence is the whole claim. It is deliberately narrower than the one this section
carried before.

| Set | CLIPPR | GRASP reference | NEB (BsaI-HFv2 37-16 cycling) |
|---|---:|---:|---:|
| 1 | 0.7812708969453874 | identical | 78% |
| 2 | 0.7826598240808516 | identical | 78% |
| 3 | 0.6556410942710312 | identical | 65% |
| 4 | 0.5962446032548091 | identical | 59% |

- CLIPPR and the reference agree **to floating-point precision**, not merely four decimals,
  by separate code paths — ours reads Pryor's `.xlsx`, theirs a `.csv` through dawdlib.
- NEB's help documents normalisation of total events to 100,000. Our workbook totals 203,364,
  so the expected count ratio is 2.03364 — which matches the observed 2.03–2.04.

### What this does *not* establish

**It does not establish the aggregation formula.** An alternative per-junction, two-strand
pooled formula gives 0.780582, 0.781971, 0.653014 and 0.594617 on the same four sets. Those
truncate to 78, 78, 65, 59 — **the same integers**. The panel cannot distinguish the two, and
NEB's current help does not document its display rounding rule. An earlier version of this
section claimed the match "established" our aggregation. It does not.

To settle it: a precision-resolved NEB output, a documented calculation or display rule, or a
discriminating set on which the two formulas differ by more than the rounding interval.

It also does not establish complete matrix identity, physical reaction coverage, or any
experimental validation. Three tools predicting from the same 2020 dataset and agreeing is
corroboration of implementation, not evidence about an assembly.

### The condition matters

First run against *BsaI-HFv2 37 static*, which put NEB 11–22 points above CLIPPR and looked
like a scorer defect. It was not. **Pryor 2020 Table S1 is the endpoint matrix for the
30-cycle 37↔16 °C reaction**, as the reference's own table name states. Against the matching
cycling row the discrepancy disappears. The error was in the instructions, not the code.

Raw data: `work/neb_observed.json`, `work/neb_vs_reference.json`. Procedure:
`validation/NEB_VIEWER_TASK.md`.

## 8. Defects the gates caught — all mine

1. **The optimiser constrained the wrong scope.** GC was enforced over the coding span handed
   to it, but the ordered fragment is the whole module and the frozen flanking bases sit
   outside that span. Two modules passed at 0.640 and hit 0.660 over the full module.
2. **Annealing was a random walk.** Temperature 0.02 against deltas of ~1e-3 gave ~90%
   acceptance, and it got one sweep where it needs many.
3. **The mode comparison was not budget-matched.** Greedy was capped at 42 proposals while
   annealing got 168. Matched, greedy gains substantially: median sharing Δ is −1.707 at
   84 proposals and −2.591 at 168.
4. **Pooling manufactured its own refusals.** First-fit turned 9 items into 4+4+1, and the
   pool of 1 fell below the minimum — though 3+3+3 is feasible.
5. **A starved search reported `complete`**, which reads as "no improvement exists" when
   nothing was tried.
6. **Absolute artefact paths** leaked the author's home directory into the result package and
   pinned it to one machine.
7. **Frame counts were transposed** in the write-up (4/8/28 → 4/28/8); the computation was
   always right.
8. **Table S1's hash was mislabelled** as the package source fingerprint.

Found by the orchestrator's full-workflow audit and fixed in this pass:

9. **Reaction membership and enzyme routing were physically wrong** — §5a. The headline
   changed as a result.
10. **Multi-pool ordering failed through the public workflow.** `plan_order` applied the
    per-pool maximum to the whole order, so nine oligos against a four-per-pool product were
    refused outright even though three pools of three is valid. Cost is now the sum over the
    planned pools.
11. **Selection was not connected to ordering.** `explore_interfaces` wrote a front but no
    inventory, and the notebook's order cell drew its sequences from the synthesis route's
    `lib.oligos()` — so a user could explore one inventory and order something unrelated.
    `select_interface` and `order_items_for` now close that path.
12. **The repeat objective measured the wrong scope.** Documented as order sequences, computed
    over assembled products. It is now measured over the distinct module inserts that are
    actually ordered, with the product figure kept as a named diagnostic.
13. **The inherited-problem exemption was far too broad.** `if problems and not
    module_constraint_problems(original.dna)` excused *every* problem in a module as soon as
    the original had any, so a module arriving with a homopolymer could acquire a new BsaI
    site and pass. Only the specific inherited problems are excused now.
14. **The result package named artefacts it did not contain**, emitting an unresolvable
    `[outside the package]` marker while still reporting `ok=True`. Artefacts are copied in,
    and a missing one prevents a complete-success verdict.
15. **"Deterministic" was wrong for the collection modes.** All three draw proposals from a
    seeded generator; they are reproducible for a fixed seed and settings, not seed-free
    deterministic. (The junction beam search genuinely uses no RNG.)

## 9. A finding about the deposited kit — RAN

Deposited module `pPR-1_D_LD5T` contains `TTTTT` and fails CLIPPR's max-homopolymer-4 rule:
**1 of 42 modules sits outside our synthesis profile.** That limit is our engineering choice,
not a published rule, so this is a statement about fit, not a defect in their kit. It surfaces
at load time through `inventories.synthesis_profile_notes()`. Recoding for *Chlamydomonas*
incidentally resolves it — the recoded inventory has none.

## 10. Tests — RAN

Counted from collection, not progress dots.

    python -m pytest --collect-only     # 883 collected
    python -m pytest                    # 883 passed

Skips are explicit: tests needing Supplementary Table S1 skip with a stated reason when it is
not supplied.

## 11. Supplied-input requirements — READ

| Input | Needed for | If absent |
|---|---|---|
| **Supplementary Table S1** | the deposited inventory | deliberately not shipped (`NOTICE.md`); `load_deposited` raises with a named reason |
| **a saved inventory** | tasks 2–5 | none — it carries its own sequences, proven in a clean install |
| **codon table** | recoding, optimisation | fetched once from Kazusa and cached, or supplied |
| **chloroplast genome** | host screening | fetched from NCBI on first use, or screening disabled with an honest status |

## 12. Relationship to the reference

GRASP Designer v0.1.15, commit `8882759ac267e79f`, is a **behavioural reference and
benchmark**. Its imports and execution live in an isolated harness; production and its tests
run without it installed. Reference implementation files were not read to guide coding, and no
functions, control flow, structures, parameter bundles or fixtures were copied or translated.

This project has had prior exposure to that repository. The process is independent
implementation using documented methods and primary data — **not** strict clean-room
development, and the engineering process alone does not settle every licensing question.

## 12a. Reproduction settings — RAN

Several figures reproduce only at a specific budget, and the defaults are smaller. Recording
them here because a result without its settings is not reproducible.

| Figure | Required setting | What the default gives |
|---|---|---|
| the saved 42-module inventory | `seeds=8` | `seeds=4` gives mean CAI 0.618830 and differs at ten modules |
| the three-mode comparison table | `--max-proposals 168` | the default 84 gives every mode a shorter budget |

Repeated four-attempt runs agree with each other; they are simply not the eight-attempt
result. All three collection modes are **reproducible for a fixed seed, order and budget** —
not seed-free deterministic. Record seed, processing order, budget and source identity beside
any saved result.

Full command list: `docs/verification_commands.md`.

## 12b. The matched oracle benchmark — RAN

    python validation/experiments/matched_oracle_benchmark.py

The acceptance requirement that stayed open through four review passes. It stayed open for a
reason worth stating: until the substrate scope was settled our deliverable was a bare insert
and theirs was a wrapped oligo, so a head-to-head would have compared two different kinds of
object. Both systems now emit a synthesis-ready sequence, so one checker can hold both.

**What this command does.** It reads a **saved** reference CSV and a **saved** CLIPPR
inventory and compares them. It executes neither optimiser, so it is a saved-output
cross-check, not a newly run head-to-head.

**What is compared.** One artefact per module per system — the sequence a vendor would be
asked to make. Theirs is `oligo_sequence_5to3` from their `optimized_library.csv`; ours is
`substrates.build(...).sequence`. All **42** modules join by name.

**Neither system is graded only on its home rules.** The checker is written inside the
benchmark, imports nothing from `clippr` or from the reference, and runs twice over both
systems' output. Reporting only "each passes its own rules" would be two self-assessments
printed side by side; the cross terms are the benchmark.

| clean ordered sequences, of 42 | under CLIPPR's rules | under the reference's rules |
|---|---:|---:|
| **ours** | **42** | 29 |
| **reference** | 6 | 41 |

Read this as a map of where the two rule sets differ, not a scoreboard:

- Our 13 failures under their regime are all **homopolymer runs of 4**. Their limit is 3;
  ours is 4. Adopting theirs would be adopting their policy, not fixing a defect.
- Their 36 failures under ours are all **GC windows above 0.65** — every one reported at
  0.660, which is the first offending window in each. Their windowed band is 0.15–0.85,
  ours 0.35–0.65.
- **Every cross-regime failure on both sides traces to a declared threshold difference.** That
  is the actual gate: an unexplained breach would be a real defect, and there are none.
- The one reference oligo that fails its **own** regime (`1A_5T_AGGT`, runs of 4) is flagged
  by their own QC as `WARNING` and emitted anyway. That is their policy, reported here, not
  something found against them.

**Codon adaptation, over the span both systems share.** Only **4 of 42** translated spans are
identical — in the other 38 the reference carries extra boundary residues — so the spans are
aligned before scoring and an ambiguous containment is refused rather than resolved:

| | CLIPPR | reference | CLIPPR ahead | reference ahead |
|---|---:|---:|---:|---:|
| over the shared span | **0.621004** | **0.658672** | **22** | **20** |
| scoring whole CDSs (a different region — not a comparison) | 0.621004 | 0.634781 | 29 | 13 |

The reference is ahead on the **mean**; CLIPPR is narrowly ahead on the **per-module split**,
22 to 20. Those point in different directions, so neither system is simply better here: the
reference's advantage is concentrated in the modules it wins. An earlier version of this
section quoted the second row, which compared different coding regions.

**Per preserved run**, since one output of a stochastic optimiser is an observation and not its
level: reference mean 0.654766 / 0.656756 / 0.656684 / 0.657501 / 0.658672 across run01–run05,
CLIPPR ahead on 23 / 22 / 24 / 24 / 22. The saved CSV behind the headline is run05, the most
favourable of the five.

**Determinism.** Ours is deterministic by construction. Across five preserved reference runs
at its declared seed, **41 of 42** modules received more than one CDS.

**Declared asymmetries**, stated rather than scored away:

- They emit vendor-acceptance and QC flags we do not; we emit pooled ordering they do not.
  **Correction, 2026-09-16:** an earlier version of this line said ligation fidelity was ours
  alone. It is not. The reference carries `ligation_fidelity.py`, wrapping GGAssembler/dawdlib
  `GGData.reaction_fidelity`, with per-stage calculators and protocol provenance (assay kind,
  cycling steps, ligase, DOIs, and a `grasp_status` marking each matrix a surrogate or proxy).
  It is absent from the reusable-library workflow's CSV, which is what was inspected; absence
  from one output was reported as absence from the tool.
- **The two wrappers are for different steps.** All 42 of their oligos release `ACAT … TTGT`,
  the **level −1** destination — they are the synthetic inserts for building the entry clone in
  pAGM1311. Ours is a level-0 BbsI cassette releasing the fragment an assembly consumes. A
  module plasmid is made at level −1 and used at level 0, so both are real and neither
  supersedes the other. Wrapper length and GC are reported, never ranked.
- The checker's enzyme exemption keys on the **recognition sequence**, not the name: BpiI and
  BbsI both read `GAAGAC`, and a name-based test reported our own two wrapper sites as
  unintended BpiI in all 42 modules. That was an artefact of the first checker, not a finding.

Artefact: `work/oracle/matched_benchmark.json`.

## 13. Runtime — RAN, on this machine

| Operation | Time |
|---|---|
| inventory import and validation (42 modules) | < 1 s |
| fixed-interface recoding, 42 modules | 3.9 s (cap 600 s) |
| collection optimisation, 168 proposals | ~4.6 s optimiser time per run |
| joint search, 100 evaluations | 0.37 s |
| compile one target | < 0.1 s |
| pooling and eligibility | < 0.1 s |
| full test suite | 4–11 min, machine-dependent (both observed) |
| release check with clean install | ~9 min |
| 200-design corpus regeneration | ~20 min |

## 14. The notebook

`notebooks/CLIPPR_designer.ipynb` is generated, not hand-edited: run
`python tools/build_notebook.py` after any API change. It now carries **38 cells** (21 code,
17 markdown) in two routes.

The original route designs a new CDS per target. The new **reusable-inventory route** wires in
the remaining five tasks as form-driven cells:

| Cell | Task |
|---|---|
| Inventory 1 | load the deposited kit, and report modules outside the synthesis profile |
| Inventory 2 | recode for your host, interfaces frozen |
| Inventory 3 | compile targets from the inventory |
| Inventory 4 | optimise the collection (greedy / anneal / random control) |
| Inventory 5 | explore junction trade-offs, with the interface-version warning |
| Inventory 6 | eligibility and pool plan |
| Inventory 7 | download the result package |

Those cells live in `tools/notebook_inventory_cells.py` rather than in the builder, which was
already long enough that appending a second complete workflow would have buried both.

Wiring this in exposed a gap in the notebook's own checker: it verified that imported names
exist, but handled neither `from X import Y as Z` aliasing nor submodule imports, so it
reported `clippr.inventories` and `clippr.workflow` as missing when both resolve. Fixed.

### The inventory cells had never executed — 2026-09-16

Three of them passed a variable named `table` as the codon table. Nothing binds `table` to a
codon table; the only assignment in the notebook is `table = result["oligos"][…]` in the
fragment-display cell, a pandas DataFrame. `recode_for_host`, `optimise_collection` and
`explore_interfaces` were each handed a DataFrame and raised as soon as a real input arrived.

Every check was green throughout. The static tests verify that imports resolve, dropdowns
resolve and regeneration is byte-stable — none of which touches whether a cell runs. The
execution harness did run the cells, but from the repository root, where `Table S1.xlsx` does
not exist, so every inventory cell took its `deposited is None` short-circuit and printed a
polite message. **The harness was exercising the absent-input branch and reporting PASS.**

Fixed in the generator: `host_codon_table` is resolved through the package's own resolver from
the form's `organism`, `codon_table_file` and `genetic_code_override`, under a name a display
cell cannot shadow. The notebook now runs **21 cells, 0 failed**, ending in a real package of
8 tasks with `ok=True`.

Two structural changes, because the branch being invisible matters more than the bug:

- `validation/run_notebook.py` **stages Supplementary Table S1 itself** and removes it
  afterwards, so the path a user takes is the path that runs. If the workbook is genuinely
  unavailable it says the check is weaker instead of passing the short-circuit.
- Release check 7 was a substring search over the notebook JSON for deferred experiment names.
  It is now **"notebook scope and execution"** and runs every cell. A notebook that parses is
  not a notebook that works, and the release suite could not previously tell them apart.

## 14a. The GC band is the dominant lever, and it is unsourced

Measured 2026-09-16, after the oracle comparison showed the reference ahead on mean codon
adaptation (0.658672 against 0.621004 over aligned spans).

Same solver, same seeds, same 12 modules, only the synthesis band changed:

| band | mean CAI |
|---|---:|
| ours — 0.35–0.65 per 50 nt, homopolymer ≤ 4 | 0.692434 |
| reference — 0.15–0.85 per 50 nt, homopolymer ≤ 3 | **0.854766** |

That first table changed **two** variables. Separated into a 2×2 over all 42 modules, the GC
effect survives — holding the cap at 4, first-12 reproduces 0.692434 → 0.854766 and all-42 with
full-substrate filtering goes **0.618830 → 0.819921**. Both systems share a codon component,
not a complete objective, so **no optimiser ranking follows.**

**`GC_BAND = (0.35, 0.65)` has no recorded source here** — but it is **not** tighter than
published guidance, as an earlier version of this section claimed. Twist publishes exactly
35–65% over 50 bp for codon optimisation. What is ours to answer for: we enforce the local half
and not the global 25–65% it is paired with, we enforce as hard what Twist frames as advisory,
and we order from IDT, whose oPools page states no GC rule at all.

**It is also not a one-line change.** `recoding.GC_BAND` is the validator's band; the solver's
`optimize_cds` keeps its own `gc_bounds` default. Patching one gives 0.657195, not 0.819921.
The real finding is that the same threshold lives in two places and neither knows about the
other.

Not changed in this pass. The full case, the vendor-source table and a six-step proposal are in
`docs/gc_band_and_constraint_model.md`. **The full case, with the vendor-profile comparison, the soft-versus-
hard constraint analysis and a five-step proposal, is `docs/gc_band_and_constraint_model.md`.**

## 15. Remaining limits

- IDT's published rules are now carried in `CURRENT_OPOOL_50PMOL` (§7), unpriced and
  transcribed rather than queried live; only the 50 pmol scale is modelled. Prices
  remain unavailable by design.
- V1's NEB panel is **done** — see §7a. The BbsI question it raised is **resolved**: level-0
  reactions are now scored with BbsI-HF (§5a).
- **What NEB's aggregation actually is remains unknown.** Two formulas fit the four observed
  integers equally well; a discriminating set or a precision-resolved output would settle it.
- **Level-1 reactions remain unscored**, for a narrower reason than before (§5d). The block
  release geometry is established and falsified two ways; the reaction's **participant list**
  is not. The deposited BsaI reaction co-assembles a linker and a DYW editing domain, and at
  least one further part is absent from the supplied records entirely. Adding only the known
  omitted ends moves 19S from 0.996046 to 0.742067, so a subset figure cannot stand in for the
  reaction.
- **The oracle benchmark covers two dimensions of plan §12, not all of them.** Sequence checks
  under both rule sets, and codon adaptation over aligned spans, are done. Repeat/sharing
  metrics, stage-specific fidelity, ordering and runtime under comparable conditions are
  **unperformed** — not passed. A PASS on the dimensions that ran is not parity.
- The matched whole-workflow benchmark against the newer reference is **done** — §12b. Its
  scope is the ordered sequence and codon adaptation. It does **not** compare assembly plans,
  because their wrapper is a level-1 cassette and ours a level-0 cassette; those are different
  physical reactions and remain uncompared.
- V2 vendor screening untouched; declared optional.
- M2's negative result stands unchanged. The Phase C and D searches have different search
  spaces and their own gates; they do not reopen it.
- No wet-lab validation of anything.
