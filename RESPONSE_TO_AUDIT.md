# Response to the constraint-enforcement work order, 2026-09-16

Four items, all executed, plus one finding of my own (§6) that retires the longest-standing
limit in the project. The previous pass's response (substrate scope, `k` propagation) is folded
into `FULL_WORKFLOW_HANDOFF.md` §5b and §5c and is not repeated here.

**State**

| | |
|---|---|
| tests | **883 collected, 883 passed** (was 846) |
| gates | B, C, D, E, F all **PASS** |
| level-1 geometry | **DERIVED AND FALSIFIED**; the reaction itself stays **unscored** — see §6 |
| matched oracle benchmark | **PASS** — the acceptance requirement that had never been run |
| release check | **10 of 10**, 0 failed, 0 needing a human (clean-install pass; check 7 now executes the notebook) |
| shipped inventory | `3abdf2e9db7d805f`, unchanged by this pass |

---

## 1. Constraint enforcement across optimisation

**Both findings reproduced exactly before anything was changed.**

`greedy seed42 fwd` on the recoded inventory returned `67cd5030521f9518`, in which
`pPR-1_19E_LN5N` carries a GC window of 0.660 at offset 55 — over the declared 0.65 ceiling,
in the sequence the export step would have ordered. Sixteen of the eighteen matched runs
returned an inventory with at least one such module.

The internal-BsaI probe `ACTC + ATGCGT×3 + GGTCTC + ATGCGT×5 + AAGA` returned
`synthesis_problems=()` while `Bio.Restriction` found one BsaI cut. `synthesis_problems`
scored GC and homopolymers and looked for no enzyme site at all.

### The repair is one function, not three patches

This was the **fourth** time in this project that a constraint was enforced on a sub-span of
what actually ships, and the previous three were each fixed in place. Fixing it a fourth time
the same way would have been the mistake.

`substrates.substrate_problems(insert, block)` is now **the** feasibility test. It builds the
released fragment, wraps it as it would be ordered, and judges that sequence. Three callers:

| caller | previously judged |
|---|---|
| `recoding.recode_inventory` | the released fragment (correct since the last pass) |
| `library_search._propose` | **the bare insert** |
| `joint_search.evaluate` | **the bare insert**, on both sides of the inherited/introduced comparison |
| `inventories.synthesis_profile_notes` | **the bare insert** — the function whose whole job is warning a user before they order |
| Gate D's front revalidation | **the bare insert** |

The last two I found myself while auditing the call sites after fixing the two you named, by
grepping for every caller of the old insert-level check. Both are the same defect; neither was
producing a wrong answer on the current inventory, but `synthesis_profile_notes` in particular
exists to tell a user what the order step will refuse, and it was looking one layer inside
that. Gate D still passes with the substrate-level check, front members unchanged.

Site checking is now part of it, role-aware: the wrapper's own two BbsI sites are excluded
**by position**, and every other occurrence of any enzyme in the active profile, on either
strand, is a breach.

### Result

| | before | after |
|---|---:|---:|
| matched runs returning an out-of-band substrate | **16 / 18** | **0 / 18** |
| objective disagreements with independent recomputation | 0 | 0 |

The adversarial probe now behaves correctly end to end:

    probe protein preserved:        True
    probe substrate problems:       ['GC 0.660 outside (0.35, 0.65) in the window at 55']
    joint_search.evaluate feasible: False
       reason: pPR-1_19E_LN5N: GC 0.660 outside (0.35, 0.65) in the window at 55

## 2. Regressions

New file `tests/test_substrate_scope.py` — eleven tests. They live in one file rather than in
the three modules they exercise because they are one defect.

| Test | Holds |
|---|---|
| `test_it_really_does_preserve_the_protein` | the probe is rejected for the right reason, not a length or translation change |
| `test_the_bare_insert_is_where_the_breach_hides` | worst GC is strictly higher on the fragment than on the insert |
| `test_the_shared_validator_refuses_it` | `substrate_problems` catches it |
| `test_the_joint_search_refuses_it_without_consulting_the_recoder` | injected past `apply_assignment`, so a correct recoder cannot mask the check |
| `test_the_clean_deposited_module_is_still_accepted` | the check rejects the probe, not everything reaching it |
| `test_an_internal_bsai_site_is_caught` | forward strand |
| `test_the_reverse_strand_counts_too` | bottom strand |
| `test_the_wrappers_own_bbsi_sites_are_not_counted` | success on the clean corpus |
| `test_greedy_seed_42_introduces_no_unorderable_module` | the named reproducer |
| `test_no_mode_introduces_one` | all three modes |
| `test_the_recoded_inventory_on_disk_is_orderable` | the artefact the rest of the workflow consumes |

**The probe is deterministic and independent of the recoder.** It is built by swapping every
codon of `pPR-1_19E_LN5N` for its most GC-rich synonym — arithmetic written in the test file —
so it cannot drift when the recoder changes.

**One judgement call, stated because it is arguable.** The optimiser tests assert that no
*new* offender is introduced, not that the output is absolutely clean. One deposited module
(`pPR-1_D_LD5T`) arrives carrying a run of five Ts and fails our homopolymer rule on arrival.
Demanding an absolutely clean output would make the deposited inventory unusable as a
baseline, and the incumbent has to be evaluable or there is nothing to compare against. If you
think that exemption is too generous, it is one line.

## 3. Gate coverage

**Phase C now validates saved full substrates.** It was calling `module_constraint_problems`
on `record.dna` — which is why it passed while 16 of 18 runs returned an inventory the next
step refused. It now calls `substrate_problems(r.dna, r.block)`, the same function the export
path uses.

**Gate F now executes collection optimisation and orders its output.** Previously it checked
that the API was present. An optimiser whose output the next step refuses passes that check.
It now runs `optimise_collection` in the clean installed environment and then
`order_items_for` on the result:

    collection optimisation: greedy: improved; 6 module(s) accepted, complete
    its output orders cleanly: True

Both new artefacts are in the Gate F package.

## 4. Regenerated comparisons

The whole Phase C table was stale. Re-run on `3abdf2e9db7d805f`, 168 proposals, 3 seeds ×
both orders, all 18 runs preserved:

| mode | sharing Δ median | best | adaptation Δ median | accepted (median) |
|---|---:|---:|---:|---:|
| greedy | −2.591 | −2.774 | −0.00204 | 18.5 |
| anneal | **−2.647** | **−3.164** | −0.00274 | 93 |
| random (control) | −0.578 | −2.121 | −0.00067 | 82 |

**The interesting part is what moved.** Annealing's advantage over greedy fell from 0.48
sharing points to **0.056** — because the candidates it had been exploiting were ones the
export step would have rejected. The "~2.6× more adaptation" language in the handoff is
withdrawn; the ratio is now 1.34×. Greedy remains the recommended default and the case for it
is now stronger, not weaker: **97.87%** of annealing's median sharing reduction, a third of the
accepted changes, smaller adaptation cost.

Also replaced:

- coverage "40 improved, 2 unchanged" → **39 improved, 3 unchanged**
- §12a `seeds=4` mean CAI 0.631980 → **0.618830** (still differs at ten modules)
- `docs/verification_commands.md` carried the same stale pair

Every results section is now preceded by §4a, which pins the run to exact source hashes:
Table S1 `6324ede3…`, codon table `a9d8f6fc…`, deposited `7d550e628605e006`, recoded
`3abdf2e9db7d805f`, reference v0.1.15 `8882759a`.

## 5. The matched oracle benchmark — run

`validation/experiments/matched_oracle_benchmark.py`, artefact
`work/oracle/matched_benchmark.json`, handoff §12b.

It compares the **ordered sequence** each system emits, per module, for all 42 — theirs
`oligo_sequence_5to3`, ours `substrates.build(...).sequence`. The checker is written inside
the benchmark, imports nothing from `clippr` or from the reference, and is run twice: once
under CLIPPR's rules, once under theirs. Grading each system only on its home rules would be
two self-assessments printed side by side.

| clean, of 42 | CLIPPR rules | reference rules |
|---|---:|---:|
| **ours** | **42** | 29 |
| **reference** | 6 | 41 |

- Our 13 failures under their regime are **all homopolymer runs of 4**. Their limit is 3.
- Their 36 failures under ours are **all GC windows above 0.65**, every one reported at 0.660
  — the first offending window in each. Their windowed band is 0.15–0.85.
- **Every cross-regime failure on both sides traces to a declared threshold difference.** That
  is the gate — an unexplained breach would be a real defect, and there are none. Requiring
  our output to be clean under *their* thresholds would make the benchmark a test of whether
  we adopted their policy.
- The one reference oligo failing its **own** regime (`1A_5T_AGGT`) is flagged `WARNING` by
  their own QC and emitted anyway. Their policy, reported, not a finding against them.

~~CAI under one scorer: ours 0.6210, reference 0.6348; ours ahead on 29.~~ **Superseded — see
§8.1.** Those figures compared different coding regions. Over the span both systems share:
ours **0.621004**, reference **0.658672**, ours ahead on **22** and behind on **20**.

Determinism: ours is deterministic by construction; across five preserved reference runs at
its declared seed, **41 of 42** modules received more than one CDS.

**A checker artefact I caught and am reporting rather than quietly fixing.** The first version
exempted intended sites by enzyme *name*. BpiI and BbsI are isoschizomers of `GAAGAC`, so it
reported our own two wrapper sites as unintended BpiI in all 42 modules — a 42/42 "failure"
that was entirely my checker. The exemption now keys on the recognition sequence.

---

## 6. Level-1 scoring — a standing limit that turned out to be a false claim

Not in your work order. I went looking because I had just told Zain that level-1 scoring was
blocked on data we did not have, and wanted to tell him where to get it.

**The data was already on disk**, fetched 2026-09-15 into the provenance directory alongside
the module records I had been using. `9SDYW_level0.gb` holds **28 assembled block plasmids**,
each carrying exactly two BsaI sites and no BbsI site — the exact mirror of the 42 module
plasmids' two BbsI and no BsaI. I had fetched them myself to verify the level-0 cut geometry
and then never opened them again.

**And they turn out not to be needed.** Digesting all 28 gives two fragment classes:

| released fragment | ends | n |
|---|---|---:|
| 499 nt | `AGGT → CTTC` | 14 |
| 406 nt | `CTTC → TTCG` | 14 |

499 + 406 − 4 = **901 nt**, which is exactly what this package compiles for 9S, and the joined
sequence occurs **verbatim** inside the deposited level-1 product `9SrpoADYW_pICH47802_lc`.
So the block vector contributes no bases: the released block fragment *is* the joined module
inserts, and a level-1 reaction's ends are the compiled product's own first and last
interfaces.

**Confirmed independently against Table S1**, consulted only after the derivation was
complete: `1E → CTTC`, `14E → GTGA`, `19E → CACG`, `2E → TTCG`, `1A` entering at `AATG` or
`AGGT` by variant. Every interface agrees.

So the standing limit — "its overhangs are set by the block plasmid and its destination, which
this package does not have" — was wrong on both clauses, and had been quoted unchallenged in
`assembly_spec`, the handoff, the verification commands and a passing test for weeks.

**What it is worth:**

| reaction | enzyme / matrix | predicted fidelity |
|---|---|---:|
| level 0 (module inserts → block) | BbsI-HF | **0.751683** |
| level 1, 9S | BsaI-HFv2 | 0.998225 |
| level 1, 14S | BsaI-HFv2 | 0.998225 |
| level 1, 19S | BsaI-HFv2 | 0.996046 |

**No conclusion about which stage is the bottleneck follows from this, and the earlier claim
that one did is withdrawn.** Ligation fidelity is a property of the whole competing overhang
set in a tube, and the level-1 set here is incomplete. Adding only the ends we *know* are
omitted moves the 19S figure from 0.996046 to **0.742067** — below level 0's 0.751683. So the
one sensitivity calculation available points the opposite way to the claim it was used to
support, and the true participant list could move it further in either direction.

What the numbers do support: the level-0 reaction is scored at 0.751683 over a **complete**
overhang set, and the level-1 block subset would score 0.996–0.998 **if it were the whole
reaction, which it is not**. Those two are not comparable and are not compared here.

**Wired in.** `ASSEMBLY_STAGES["level1"]` now declares `flanks_from_product`,
`_stage_summary` takes the product's terminal interfaces, and `joint_search` scores level-1
reactions against BsaI-HFv2. The headline minimum-across-reactions figure is unchanged,
because level 0 still sets it.

**The guard.** Every scored reaction records the overhang set it was scored on and where those
ends came from (`ends_from`). A fidelity without a set is the same defect, wearing a number
instead of a `None`, and there is a test for it.

**§8.2 supersedes the completeness claim in this section.** Level 1 is **not** scored: the
block geometry is established but the reaction's participant list is not, so
`unscorable_reactions` lists it again and the PPR-only figure is a labelled diagnostic.

**How far the falsifier reaches.** Only the **9S** block plasmids are deposited, so 9S is
verified two independent ways — by digesting real plasmids and by reading Table S1 — and they
agree. The 14S and 19S sets follow from the same module interfaces in Table S1, which the 9S
case shows to be the right place to read them, but there are no 14S/19S block plasmids to
digest. They are one falsifier short of 9S, and that is stated rather than glossed.

**Scope, stated because this is where I keep going wrong.** The deposited construct
co-assembles parts this package does not compile — a P2L2S2 linker (`TGTG → CAAC`) and the
consensus DYW domain (`CAAC → GCTT`) in the same BsaI reaction, plus at least one further
part bridging `TTCG` to `TGTG` absent from the supplied records. **The level-1 number above is
for the PPR block junctions only.** Adding an editing domain adds overhangs that are not in
it, and `LEVEL1_COASSEMBLED_PARTS_UNMODELLED` says so in the source.

Derivation and both falsifiers: `validation/experiments/level1_geometry.py`, artefact
`work/level1/geometry.json`.

---

## 7. A defect in the validation harness itself

Found while running the final release sequence, not by design.

`validation/independent/v5_digest_ligate.py` ended with:

    return 0 if ok == len(rows) and detected == 5 else 1

Two corrupt controls were added after that line was written. The script therefore reported
`corrupt_controls_detected: 7` of `corrupt_controls_total: 7` and `failures: []` in its JSON --
everything clean -- and **exited 1 anyway**, because `detected == 5` had become unsatisfiable.

It went unnoticed because the release check reads the JSON, not the exit code, so check 5 kept
passing on the evidence while the process kept saying it had failed. It only surfaced when a
shell loop checked the return status.

Fixed to compare against `len(corrupt_controls(rows[0]))` rather than a literal. Worth naming
because it is the same shape as the defects this pass was about: **a hardcoded expectation
that silently drifts away from what is actually measured.** A checker whose pass condition is
written down separately from the thing it counts will eventually disagree with itself, and the
disagreement is invisible as long as only one of the two is read.

---

## 8. Response to the independent review of 2026-09-16

Both findings reproduced before anything was changed, and both are accepted.

### 8.1 The CAI comparison used different coding regions

Reproduced exactly: **4 of 42** translated spans are identical, 38 differ, 0 are ambiguous.

| comparison | CLIPPR | reference | CLIPPR ahead | reference ahead |
|---|---:|---:|---:|---:|
| as submitted (unaligned) | 0.621004 | 0.634781 | 29 | 13 |
| **over the shared span** | **0.621004** | **0.658672** | **22** | **20** |

This is the same scope error as the previous four, relocated. Every earlier instance was a
*constraint* checked on a sub-span of what ships; this one is a *comparison* computed over a
sub-span of what both systems have in common. That is why four passes spent looking at
constraints did not find it.

**Fixed at the source, not in the write-up.** `matched_oracle_benchmark.py` now aligns the
translated spans before scoring, records each module's reference nucleotide bounds, protein
lengths and trimmed residue count, and **refuses** a module whose protein is absent or
ambiguously contained rather than choosing an offset. The whole-CDS figure is retained in the
output as `reference_mean_unaligned_do_not_quote` so the difference stays visible.

**Per preserved run, not one CSV.** As requested, the same alignment now runs over all five
preserved reference runs:

| run | reference mean | CLIPPR ahead, of 42 |
|---|---:|---:|
| run01 | 0.654766 | 23 |
| run02 | 0.656756 | 22 |
| run03 | 0.656684 | 24 |
| run04 | 0.657501 | 24 |
| run05 | 0.658672 | 22 |

Worth noting: the saved CSV used for the headline is run05 — the run where the reference scores
highest of the five. Not chosen for that reason, but the headline was resting on the most
favourable observation of a stochastic optimiser, which is exactly the hazard you named.

Your two fresh seed-42 runs (0.653995 and 0.653187, CLIPPR ahead on 25) sit just below this
range and are consistent with it.

**Regression:** `tests/test_oracle_alignment.py`, 9 tests. It asserts 4/42 identical, 42/42
unambiguous, the aligned figures 0.621004 / 0.658672 / 22, **and** the unaligned 0.634781 / 29,
with an explicit assertion that the two do not coincide — so a silent regression to unaligned
scoring fails here rather than being published.

### 8.2 Level-1 completeness contradicted its own caveat

Accepted, and the numeric sensitivity makes it sharper than "unsupported":

| set | BsaI-HFv2 |
|---|---:|
| 19S PPR block subset | 0.996046 |
| the same plus the **known** omitted ends `TGTG`, `CAAC`, `GCTT` | **0.742067** |
| level 0, participants complete | 0.751683 |

The known omission drops 19S *below* level 0. So "the bottleneck is unambiguously level 0" was
not merely beyond its evidence — the one sensitivity calculation available points the other
way. Withdrawn, along with "chosen poorly", from the handoff, this document and the derivation
script.

**The structured status now matches the prose**, which was the real defect: a caveat sitting
beside `context_available: True` is not a caveat, because a reader believes the field.

| field | now |
|---|---|
| `geometry_established` | `True` — the block release geometry, falsified two ways |
| `participants_established` | `False` |
| `context_available` | `False` |
| `reaction_overhangs` | `null` |
| `scorable` | `False` |
| `unscorable_reactions` | level 1 is listed again |
| `block_subset_overhangs` | the PPR-only set, in its own field, labelled a diagnostic |

Nothing was deleted: the verified geometry and the subset figure are both retained, in fields
that say what they are.

**Regressions:** `tests/test_level1.py` was rewritten — it had been passing while encoding the
incomplete-reaction assumption, exactly as you predicted. It now asserts the reaction is *not*
scorable, that the subset never appears as a fidelity, that the structured flags agree with the
prose, and the 0.996046 → 0.742067 → below-level-0 inequality that retires the headline.

### 8.3 Smaller corrections

- "within 2%" → greedy reaches **97.87%** of annealing's median sharing reduction. "About 98%"
  is what the measurement supports.
- The benchmark reads a saved reference CSV and a saved CLIPPR inventory. It is a
  **saved-output cross-check**, not a newly executed head-to-head optimisation, and is now
  described that way.
- Plan §12's benchmark protocol also covers repeat/sharing metrics, stage-specific fidelity,
  ordering and runtime under comparable conditions. Those dimensions are **unperformed**, not
  passed; a PASS on the sequence-check and CAI dimensions is not parity.

---

## 9. The notebook's inventory workflow had never executed

Not from any review. Found because Zain asked whether the Colab notebook works, and I checked
instead of answering from the test results.

`validation/run_notebook.py` reported **21 cells executed, 0 failed, PASS**. It was running the
wrong branch. The notebook looks for `Table S1.xlsx` by bare filename, as a Colab user would
after uploading it; the harness ran from the repository root where no such file exists, so
every inventory cell took its `deposited is None` short-circuit and printed a polite message.

Supplying the file:

    cell 16  Inventory 2 - recode for your host:
             AttributeError: 'int' object has no attribute 'replace'

**The cause.** Three inventory cells pass a variable named `table` as the codon table.
Nothing ever binds `table` to a codon table. The only assignment in the whole notebook is in
the fragment-display cell: `table = result["oligos"][...]` — a pandas DataFrame. So
`recode_for_host`, `optimise_collection` and `explore_interfaces` were each handed a DataFrame
where a codon table was expected, and raised the moment a real input reached them.

**The inventory half of the notebook had therefore never run.** It was written, and it was
statically tested — imports resolve, dropdowns resolve, regeneration is byte-stable, 19 tests
green — and none of that touches whether a cell runs.

**Fixed in the generator**, not the notebook: `tools/notebook_inventory_cells.py` now resolves
`host_codon_table` through the package's own resolver from the form's `organism`,
`codon_table_file` and `genetic_code_override`. The name is deliberately distinct so a display
cell cannot shadow it again. Two further undefined names surfaced on the way (`genetic_code`,
which the form calls `genetic_code_override`) and are fixed too.

Now, with the workbook present: **21 cells executed, 0 failed**, and the final cell writes a
real package — 8 tasks, `ok=True`, 0 failures.

**Two structural fixes, because the bug is less interesting than why it was invisible:**

1. `run_notebook.py` now **stages Supplementary Table S1 itself** and removes it afterwards, so
   the path a user takes is the path that runs. If the workbook is genuinely unavailable it
   says so and marks the check weaker, rather than reporting a pass for the absent branch.
2. Release check 7 was **"notebook scope"** — a substring search for deferred experiment names
   in the notebook JSON. It is now **"notebook scope and execution"** and actually runs every
   cell. A notebook that parses is not a notebook that works, and the release suite could not
   previously tell the difference.

This is the same shape as the substrate-scope defects: something was checked, the check was
green, and the check was of a narrower thing than the claim it was supporting. Here the
narrower thing was *a branch*.

---

## 10. Where the CAI gap actually comes from

Asked by Zain: the reference optimises codons better on the mean (0.658672 against 0.621004
over aligned spans) — why, and can we close it.

**It is not a better optimiser. It is a looser constraint**, and the effect is four times the
size of the gap it explains.

Controlled experiment: same solver (DNA Chisel, `CodonOptimize(method="use_best_codon")`), same
seed set, same 12 modules, same everything except the synthesis band.

| band | mean CAI |
|---|---:|
| ours — 0.35–0.65 per 50 nt, homopolymer ≤ 4 | 0.692434 |
| theirs — 0.15–0.85 per 50 nt, homopolymer ≤ 3 | **0.854766** |

The band is worth **~0.16 CAI**; the reported gap is 0.037. This is consistent with the matched
benchmark, where **36 of 42** of their ordered sequences breach our GC window, every one at
0.660. They buy codon adaptation with GC we forbid.

**The finding that follows.** `recoding.py:33` declares `GC_BAND = (0.35, 0.65)` with a careful
comment about *scope* and **no source**. The IDT oPools profile this band presumably serves
constrains oligo length and pool size and **says nothing about GC**.

**Corrected §12.1:** I also wrote that the band was tighter than any published vendor guideline.
It is not — Twist publishes exactly 35–65% over 50 bp. And §12.2: widening it is **not** one
line; the solver and the validator hold the threshold separately.

Recommended next step: **one versioned constraint profile consumed by every stage**, with
today's strict settings preserved as the named baseline. See §12.2 — the band is not a constant
to widen but a policy that has to reach the solver, the validator, export and the manifest
together.

A smaller lever noted in passing: `recode_inventory` calls `optimize_cds(...,
unique_kmer_size=None)`, disabling the repeat objective inside the solve, after which the
collection optimiser attacks the same axis separately. Two steps on one axis at different
times; worth testing whether one pass inside the solve dominates.

The full case — vendor-profile comparison, the soft-versus-hard constraint analysis, and a
five-step proposal — is `docs/gc_band_and_constraint_model.md`.

**Note on method.** The measurement above was made without reading the reference: its optimiser
settings are in the exported `config_as_run`, and the cause was settled by a controlled
experiment on CLIPPR's own code. Zain has since authorised reading the reference to learn from
(not to copy), and §11 records what that added — including two corrections to claims made here
before it was read. `docs/methodology_sources.md` §7 must be updated to match that change of
posture before either document is circulated.

---

## 11. What reading the reference added — and two corrections it forced

Zain authorised reading `grasp-library-designer` to learn from, not to copy, on 2026-09-16.
Recorded here because it changes two claims made earlier in this document and one design
question that was open.

### 11.1 Correction: the reference does compute ligation fidelity

§5 and the handoff said fidelity was absent from the reference. It is not.
`grasp_library/ligation_fidelity.py` — 442 lines — wraps GGAssembler/dawdlib's
`GGData.reaction_fidelity` with per-stage calculators and protocol provenance: assay kind,
cycling steps, ligase, buffer, Potapov and Pryor DOIs, and a `grasp_status` field marking each
matrix `"surrogate"` or `"proxy"` together with what it proxies for.

The error was inferential. I inspected the reusable-library workflow's 56-column CSV, found no
fidelity column, and reported absence from the *tool*. Absence from one output is not absence
from the package.

### 11.2 Correction: our level-1 flanking pair was incomplete, and their model says how

Their `grasp_level1_reaction_overhangs` composes the reaction as

    [left_final_cassette, left_ppr_outer] + block_joins + [right_ppr_outer, right_final_cassette]

— **four** flanking overhangs, with both pairs **required** and no defaults; it raises if either
is missing.

So `DESTINATION_OVERHANGS["level1"] = ("GGAG", "CGCT")` was never a competing answer to the
`AATG … TTCG` pair derived from the block plasmids. It is the *final cassette*; the PPR outer
pair sits inside it; the real reaction contains both.

| | our current set | under their composition |
|---|---:|---:|
| 9S | 0.998225 | **0.838330** |
| 19S | 0.996046 | **0.836500** |

This corroborates the review's finding from an independent direction, and supplies the
structure the participant list was missing. Two further points of agreement: their scored entry
point is named `grasp_first_stage_fidelity` — scoped to level 0, as ours now is — and their
Pareto search refuses to score partial sets, commenting that "partial beam states cannot yet
form a physical reaction." The discipline the review asked of us is the discipline they already
keep.

### 11.3 The architectural difference worth considering

They optimise a single weighted scalar under simulated annealing, in which GC, homopolymer and
repeat limits are **soft squared penalties** and only a forbidden restriction site is an
absolute veto. We give DNA Chisel the same limits as **hard constraints**.

Theirs therefore trades a hot GC window for a large codon gain and ships the result marked
`WARNING`; ours discards the candidate and returns the module *unchanged*, carrying its
original poorly adapted sequence. Our "unchanged" count is not an absence of opportunity — it
is the band vetoing every improvement, and it contributes to the CAI gap independently of the
band's width.

### 11.4 Three smaller ideas worth taking

- **`minimum_relative_adaptiveness = 0.20`** — the codon alphabet is pruned before the search
  runs, putting a floor under CAI and shrinking the space.
- **Rescue codons** — a scoped escape hatch below that floor, used *only* to repair a forbidden
  site, logged per codon, with an assertion that no unlogged position sits below the floor. We
  reject the whole candidate instead; theirs is better and no less honest.
- **Library-similarity penalty inside the per-module solve.** They pass other members' k-mers
  into the objective; we run a separate collection optimiser afterwards. Two passes on one axis
  where they use one.

### 11.5 What stays ours

Determinism; versioned inventories with content hashes and reload reproduction; the
released-fragment scope discipline; explicit participant-completeness status on every reaction;
collection optimisation with a random control; cross-regime benchmarking; one shared
feasibility validator on every acceptance path.

Their vendor modelling is more mature than ours. Our discipline about what a number means is
more developed than theirs. Neither observation is flattering to the other side of the ledger.

### 11.6 Documentation consequence

`docs/methodology_sources.md` §7 states that "reference implementation files are not read to
guide coding." That sentence is now false and must be rewritten before either document
circulates. The reference repository still declares **no licence**, so nothing is vendored,
translated or redistributed — derived facts ship, source files do not.

---

## 12. Corrections from the GC-policy audit

Four claims of mine were wrong. All four reproduced; I verified the first two myself rather
than accepting them.

### 12.1 "Tighter than any published vendor guideline" — **false**

Twist's gene-synthesis guidance states verbatim: *"We avoid fitted sequences that create global
GC% of less than 25% or more than 65% and local GC windows (50 bp) of less than 35% or more
than 65%."* That is CLIPPR's band and window **exactly**. Fetched and read directly.

I had compared against the reference implementation's *transcription* of Twist's rules, in a
document where I had myself written that the transcription should not be trusted for this
purpose. The narrower points survive: the band has no recorded source here, Twist must not be
retro-fitted as its origin, it is Twist's **codon-optimisation** guidance while we order from
IDT, we enforce the local half and not the global 25–65% it is paired with, and we enforce as
hard what that page frames as advisory.

### 12.2 "Widening the band is one line" — **false**

`recoding.GC_BAND` is the validator's; `codons.optimize_cds` carries its own
`gc_bounds=(0.35, 0.65)` default and `recode_inventory` never passes one. Patching the constant
alone gives 0.657195, not the fully broadened 0.819921. The solver keeps searching the narrow
space while the validator accepts a wider one.

That is a better finding than the one I claimed: **the same threshold lives in two places and
neither knows about the other.**

### 12.3 The experiment changed two variables

I varied the GC band *and* the homopolymer cap and reported it as "only the band changes". The
orchestrator's 2×2 over all 42 modules separates them, and the GC effect survives: holding the
cap at 4, first-12 reproduces 0.692434 → 0.854766, and all-42 with full-substrate filtering goes
0.618830 → 0.819921. The cap contributes little. Two things my 12-module cohort could not show:
coding-span CAI is not a deliverable, and in the cap-3 cells eight wide-band winners fail the
full-substrate homopolymer target — none of them in the first twelve.

### 12.4 "There is no artefact to carry forward" — **false**

`workflows.export_optimized_library` writes CSV, FASTA, Excel and annotated GenBank;
`workflows.compile_and_assemble_target` compiles a target against an already-optimised library.
I inferred absence from running one entry point — the same error as §11.1, made twice in the
same document. Whether their exports carry anything equivalent to our content-hashed version
and mismatch refusal is **unverified**, and cannot be inferred either way.

### What I am not claiming

No optimiser ranking. No whole-workflow superiority over GRASP. The 10/10 release check is a
statement about that checklist's scope. Hard-versus-soft at matched thresholds has not been
tested, and neither has a wider band's vendor acceptance.

---

## Still open, unchanged

- ~~**Level-1 reactions cannot be scored**~~ — **closed, and the claim was wrong.** See
  §6 below.
- **NEB's aggregation formula is still undetermined.** Two formulas fit the four observed
  integers equally well.
- **Assembly plans are not compared.** Their wrapper is a level-1 BsaI cassette, ours a
  level-0 BbsI cassette; those are different physical reactions and the benchmark does not
  pretend otherwise.
- The wrapper's spacer and padding composition are CLIPPR engineering choices, not published
  kit rules, and are stated as such in `substrates.py`.
- `repeat_burden` is 0 at k=20 over distinct inserts and over the wrapped sequences, so that
  axis remains degenerate.
- No wet-lab validation of anything.
