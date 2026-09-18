# Findings about the GRASP reference implementation

*Everything CLIPPR has established about `grasp-library-designer` — what it does, where it
agrees with us, where it differs, what it is better at, and what reading it corrected in our
own claims.*

**Subject.** GRASP Designer, `grasp-library-designer` v0.1.15, commit `8882759ac267e79f`;
21 modules, 12,578 lines under `grasp_library/`. Primary records and the earlier library
outputs come from `github.com/farleykvdg/GRASP` at commit `b569a8c4`.

**Standing.** It is a **behavioural reference and benchmark**, never implementation material.
The repository declares **no licence**: files are fetched for verification and never
redistributed, and derived facts (an overhang set, a threshold) may ship while source files may
not. Reading it to understand the method was authorised 2026-09-16; nothing is copied,
translated, or renamed into CLIPPR. `docs/methodology_sources.md` §7 carries the posture.

**How each finding was obtained** is marked, because it matters:

| mark | meaning |
|---|---|
| **[A]** | from exported artefacts — their CSV, `config_as_run`, preserved run JSON |
| **[E]** | from a controlled experiment on CLIPPR's own code |
| **[R]** | from reading their source, after authorisation |

---

## 1. What the tool is

A **codon optimiser with synthesis QC, Golden Gate junction search, and an oligo exporter**,
driven by a config dictionary and a Colab control panel. It takes amino-acid sequences plus a
per-codon coding mask and returns synthesis-ready oligos with a QC verdict. **[A][R]**

Its scope overlaps CLIPPR's design and recoding routes.

**Correction, 2026-09-16.** An earlier version of this line said it "has no concept of a
reusable versioned inventory: it redesigns all 42 modules from amino acids on every run." That
is **false**. `workflows.export_optimized_library` writes CSV, FASTA, Excel and annotated
GenBank, and `workflows.compile_and_assemble_target` takes an already-optimised library and
compiles a target against it without re-running optimisation. I inferred the absence from
running one entry point — the same error as §5.1, made twice in one document.

What a narrower comparison would still have to establish: whether their exported library carries
anything equivalent to CLIPPR's content-hashed inventory version, its reload-reproduction
guarantee, or its refusal to compile against a mismatched interface assignment. That comparison
has not been done, and the absence of such a contract **cannot** be inferred from the presence
or absence of an export function.

## 2. Where it independently agrees with us

These are the most valuable findings, because both systems reached them separately.

| finding | theirs | ours |
|---|---|---|
| level-0 destination overhangs | `GRASP_LEVEL0_EXTERNAL_OVERHANGS = ("CTCA", "CGAG")` **[R]** | `DESTINATION_OVERHANGS["level0"]` — identical |
| level-0 enzyme and matrix | "the Level 0 (BbsI-HF / BpiI) reaction score" **[R]** | BbsI-HF, after our reaction-model correction (§5a) |
| fidelity is per-reaction, never multiplied | "every PPR block uses the same junction set in a separate tube, so the objective does not multiply unrelated transformations" **[R]** | `docs/objectives.md` §5 — never an unlabelled product |
| partial overhang sets are not reactions | "partial beam states cannot yet form a physical reaction" **[R]** | participant-completeness status on every reaction |
| level 1 is not the scored stage | entry point is named `grasp_first_stage_fidelity` **[R]** | level 1 reported unscorable |
| synthesis QC is judged on the **full oligo** | prefix + CDS + suffix **[R]** | the released fragment / wrapped substrate |

The last four are the same disciplines an independent review demanded of CLIPPR this week. The
reference already kept them. That is corroboration of the review, not of us.

## 3. Where it works differently

### 3.1 Soft weighted penalties, not hard constraints **[R]**

Their objective is a single scalar maximised by simulated annealing:

    score = w_codon·codon_score
          − w_global_gc·gc_penalty²      − w_local_gc·local_gc_penalty
          − w_homopolymer·excess²        − w_internal_repeat·repeat_penalty
          − w_library_similarity·shared_kmer_penalty

Only a forbidden restriction site is an absolute veto (`−1e12`). Weights are configuration.

CLIPPR hands the same limits to DNA Chisel as **hard constraints**, with codon adaptation as
the objective. A candidate one window-point over the band is discarded.

**Consequence.** Theirs trades — a hot GC window for a large codon gain — and ships the result
marked `qc_status: WARNING`. Ours refuses and returns the module *unchanged*. Our "unchanged"
count is not an absence of opportunity; it is the band vetoing every improvement, and the three
modules concerned are named in `docs/synthesis_policy.md`.

**Measured 2026-09-17, and it vindicates their architecture.** An earlier version of this
paragraph claimed hard and soft enforcement produce identical output and that soft has no
independent value. That experiment never implemented soft enforcement -- `solver_bounds()`
returned the same arguments for `hard` and `target` -- and its conclusion was retracted.

Rerun against the repaired solver: enforcing CLIPPR's **strict** band as a target reaches mean
CAI 0.819921, identical module for module to widening the band to 0.15-0.85, against 0.618830
when the same thresholds are hard. So the soft-penalty approach the reference uses recovers the
full +0.201091 that a wider band buys, while still reporting all 40 out-of-band windows.

Their architecture is therefore a real advantage and not, as this document briefly claimed, a
compensation for a tight threshold. What it buys is the ability to keep a conservative
threshold as *guidance* without paying for it in codon adaptation.

### 3.2 Vendor profiles with a three-way rule split **[R]**

Five named profiles, each with a product URL:

| profile | global GC | window | homopolymer target | machine-hard |
|---|---|---|---|---|
| Twist · Express / Low complexity | 0.25–0.65 | 0.20–0.80 / 50 nt | 3 | ≤ 13 |
| Twist · Standard guidelines | 0.25–0.65 | 0.15–0.85 / 50 nt | 3 | ≤ 13 |
| Twist · Complex Genes tolerant | 0.25–0.65 | 0.10–0.90 / 50 nt | 3 | ≤ 30 |
| IDT · gBlocks / eBlocks conservative | 0.25–0.70 | 0.20–0.80 / 40 nt | 3 | — |
| Generic · conservative (default) | 0.30–… | 0.25–0.75 / 50 nt | 3 | — |
| **CLIPPR** | **none** | **0.35–0.65 / 50 nt** | 4 | — |

**That table is their transcription, not the vendors' pages**, and trusting it produced a false
claim in `docs/gc_band_and_constraint_model.md`. Twist's own guidance gives **35–65% over 50 bp**
for codon optimisation — exactly CLIPPR's band, and matching no row above. Read the primary
pages before adopting any number here.

And three distinct categories where we have one:

- `synthesis` — design targets, optimised against softly
- `machine_hard_constraints` — what is actually rejected (Twist's real homopolymer rule is
  < 14, not ≤ 3)
- `manual_review_rules` — vendor or human judgement, e.g. "no CcdB"

Our single band matches Twist's published codon-optimisation guidance, which no row above
reproduces. What is genuinely ours to answer for: we enforce the local half of Twist's pair and
**no global band at all**, and we enforce as hard what that page frames as optimisation
guidance.

### 3.3 A pruned codon alphabet, with a scoped escape hatch **[R]**

`build_allowed_codons(..., minimum_relative_adaptiveness=0.20)` removes rare codons from the
search space before annealing starts — a floor under CAI and a smaller space to search.

A second alphabet is built at floor `0.0` as `rescue_codons`, used **only** to repair a
forbidden restriction site, with every use logged per codon and an assertion that no unlogged
position sits below the floor. A hard constraint outranks a soft preference, the exception is
bounded, and it is visible in the output (`rescue_codon_count`, `rescue_codon_detail`,
`rescue_codon_reason`).

CLIPPR has neither the floor nor the repair: we discard the candidate.

### 3.4 Three Pareto objectives, one of them synthesis fitness **[R]**

`ObjectiveScores(ligation_fidelity, codon_optimality, synthesis)`, all maximised, with
`level_minus1_fidelity`, `level0_fidelity` and `level1_fidelity` **reported alongside but
excluded from the dominance test**.

CLIPPR's three are fidelity, adaptation and `repeat_burden` — and `repeat_burden` is
identically 0 at k=20 over distinct inserts, so our search is effectively two-objective. Theirs
uses a composite synthesis fitness in that slot, which does not degenerate.

They also handle the degeneracy problem directly: `prefer_ideal` applies a continuous
preference (GC → 0.5, a soft homopolymer ramp) "so Pareto ranking still moves under loose
vendor hard limits". A feasibility test gives no gradient once constraints are loose; a
continuous preference does.

### 3.5 Library similarity inside the solve **[R]**

`sequence_objective(..., external_kmers)` penalises k-mers shared with other library members
**during** each module's optimisation. CLIPPR recodes each module alone and then runs a
separate collection optimiser over the finished set — two passes on one axis where they use one.

### 3.6 Ligation matrices carry their protocol **[R]**

`ligation_fidelity.py` wraps GGAssembler/dawdlib `GGData.reaction_fidelity`, and each matrix is
annotated with `assay_kind` (`static_ligation` vs `golden_gate_cycling`), temperature and hours
or cycle steps, ligase, buffer, terminal step, `source_doi`, `source_url`, and a `grasp_status`
of `"surrogate"` or `"proxy"` together with `proxy_for`.

We hold the same Pryor and Potapov matrices with far less provenance attached — which is
exactly the axis where CLIPPR otherwise claims to be strict.

### 3.7 The level-1 reaction composition **[R]**

    grasp_level1_reaction_overhangs:
        [left_final_cassette, left_ppr_outer] + block_joins + [right_ppr_outer, right_final_cassette]

**Four** flanking overhangs, both pairs **required**, raising if either is absent.

This resolved an open CLIPPR question — see §5.2.

### 3.8 Pricing **[R]**

`idt_opools.py` carries a full IDT oPools list-price model: a €109.00 / 3,300-base pool floor,
per-base tiers at €0.038 / €0.025 / €0.013, €1.63 per oligo for phosphorylation, and all three
scales (1, 10 and 50 pmol) with their distinct count ranges.

CLIPPR models only the 50 pmol scale and **refuses to price at all**, on the grounds that an
undated number presented as current is worse than none. Their figures are undated too — but
they are *there*, and ours are not. Whether that is discipline or a gap is a judgement call
worth making explicitly rather than by default.

## 4. Where CLIPPR is ahead

- **Determinism.** Ours is deterministic by construction. Across five preserved runs at their
  declared seed 42, **41 of 42** modules received more than one CDS **[A]**; two further fresh
  runs reproduced the same 41/42 **[A]**. A user cannot reproduce their own result.
- **Versioned inventories.** Content-hashed, reload-reproducing, and refusing to compile against
  a mismatched interface assignment. They *do* export a reusable library and can compile a
  target against it (§1) — what is unverified is whether those exports carry an equivalent
  identity-and-mismatch contract. Claiming an advantage here needs that comparison; it is not
  established.
- **One shared feasibility validator.** Every acceptance path in CLIPPR calls
  `substrates.substrate_problems`. Their constraint logic is distributed across the objective,
  the QC pass and the vendor profile.
- **Explicit completeness status.** Every reaction reports whether its participants are known.
- **A random control** in collection optimisation — a method that cannot beat it has not been
  shown to work.
- **Cross-regime benchmarking.** Both systems graded under both rule sets, with every
  cross-regime failure attributed to a declared threshold difference.

Their vendor modelling is more mature. Our discipline about what a number means is more
developed. Neither observation flatters the other column.

## 5. Corrections this reading forced

### 5.1 "The reference computes no ligation fidelity" — **wrong**

Stated in the handoff and in `RESPONSE_TO_AUDIT.md` §5 on the basis that its reusable-library
CSV has no fidelity column among 56, and its `config_as_run` names no matrix. Both true; the
conclusion was not. `ligation_fidelity.py` is 442 lines of exactly that, with per-stage
calculators.

**Absence from one workflow's output was reported as absence from the tool.** The failure was
inference presented as fact, and it is the same error shape as every scope defect this project
has recorded: a narrower observation carrying a wider claim.

### 5.2 Our level-1 flanking pair was incomplete

CLIPPR derived the PPR block ends as `AATG`/`AGGT … TTCG` from the 28 deposited block plasmids,
and treated `DESTINATION_OVERHANGS["level1"] = ("GGAG", "CGCT")` as a competing, unverified
alternative. Their composition (§3.7) shows it is neither competing nor wrong: it is the
**final cassette** pair, and the real reaction contains both pairs.

| | CLIPPR's current set | under their composition **[E]** |
|---|---:|---:|
| 9S | 0.998225 | **0.838330** |
| 19S | 0.996046 | **0.836500** |

This corroborates the independent review's finding — that a PPR-only subset is not the
reaction's fidelity — from a direction neither the review nor CLIPPR had used.

### 5.3 "Their wrapper is a level-1 BsaI cassette" — **wrong** **[A]**

Digesting their oligos shows all 42 release `ACAT … TTGT`, which is the **level −1**
destination. Their oligo is the synthetic insert for building the entry clone in pAGM1311;
CLIPPR's substrate is the level-0 cassette releasing the fragment an assembly consumes. A
module plasmid is *made* at level −1 and *used* at level 0 — both artefacts are real, for
different steps, and neither supersedes the other.

## 6. The measured comparison

| | CLIPPR | reference |
|---|---|---|
| mean CAI, aligned spans **[A]** | 0.621004 | **0.658672** |
| modules ahead **[A]** | **22** | 20 |
| ordered sequences clean under CLIPPR's rules **[A]** | **42 / 42** | 6 / 42 |
| ordered sequences clean under the reference's rules **[A]** | 29 / 42 | **41 / 42** |
| CDS reproducibility at one declared seed **[A]** | deterministic | 41 / 42 differ |

Every cross-regime failure on both sides traces to a **declared threshold difference** — ours
are all homopolymer runs of 4 against their limit of 3; theirs are all GC windows at 0.660
against our ceiling of 0.65. Neither system breaches the other's rules in a way the rule sets
do not already predict.

**The CAI difference is dominated by a constraint difference. [E]** My first experiment changed
the band *and* the homopolymer cap while reporting only the band; separated into a 2×2 over all
42 modules, the GC effect survives — holding the cap at 4, the first-12 figures are 0.692434 →
0.854766, and over all 42 with full-substrate filtering, 0.618830 → 0.819921.

Both systems share a mean-log-relative-adaptiveness *component*; theirs also carries weighted
GC, homopolymer, repeat and similarity penalties in the same scalar. That is not a shared
objective, and no optimiser ranking follows from any of this. Full argument and its limits:
`docs/gc_band_and_constraint_model.md`.

## 7. What is worth adopting, ranked

Each is a method to reimplement from the documented approach, not code to take.

1. **Vendor profiles with the target / machine-hard / manual-review split** (§3.2). Fixes an
   unsourced constant, makes IDT-versus-Twist a user choice, and converts an invisible ceiling
   into a declared trade. Any threshold adopted must be re-sourced from the vendor's own page —
   their table is *their* transcription.
2. **A minimum-relative-adaptiveness floor on the codon alphabet** (§3.3). Cheap, and puts a
   floor under CAI.
3. **Scoped rescue codons with a logged exception and an assertion** (§3.3), replacing our
   discard-the-candidate behaviour.
4. **The level-1 cassette ends** (§3.7), to complete the participant list — with the
   co-assembled linker and DYW domain still outstanding.
5. **A non-degenerate third Pareto axis** (§3.4): a composite synthesis fitness in place of
   `repeat_burden`, plus a continuous `prefer_ideal` preference so ranking still moves when
   constraints are loose.
6. **Library similarity inside the per-module solve** (§3.5), rather than a second pass.
7. **Protocol provenance on the ligation matrices** (§3.6) — assay kind, cycling steps, DOI,
   and an explicit surrogate/proxy status.

None of these is implemented. Items 1 and 4 change published numbers and should be decisions,
not edits.

## 8. What this reading does **not** establish

- **That their output manufactures.** No sequence from either system has been synthesised.
  Twist and IDT score submissions with proprietary models; a published guideline is not an
  acceptance guarantee.
- **That a higher CAI is a better protein.** It is a sequence statistic. `docs/objectives.md`
  §8 states the limit and it applies to both systems.
- **That their thresholds are correct.** They are their transcription of vendor documentation,
  read here and not yet verified against the primary pages.
- **That their level-1 composition is complete.** It omits the co-assembled P2L2S2 linker and
  DYW domain that the deposited construct's own records show in the same BsaI reaction.
- **Anything about the quality of their code.** Not assessed, not the point, and not ours to
  judge.

## 9. Module map, for orientation

| module | lines | what it holds |
|---|---:|---|
| `optimizer.py` | 1124 | annealing, weighted objective, allowed-codon construction, site repair, synthesis QC, library driver |
| `gga_split.py` | 1125 | Golden Gate splitting |
| `import_grasp.py` | 961 | reading the deposited records |
| `workflows.py` | 967 | orchestration |
| `control_panel.py` | 984 | Colab UI |
| `genbank_export.py` | 920 | export |
| `synthesis_vendors.py` | 560 | vendor profiles, ligation protocol metadata, enzyme sets |
| `oneshot.py` | 564 | single-target route |
| `objectives.py` | 465 | the three Pareto objectives |
| `ligation_fidelity.py` | 442 | dawdlib/GGAssembler wrapper, per-stage reaction composition |
| `pareto.py` | 304 | dominance, front, knee point, overhang search |
| `restriction_sites.py` | 209 | site handling |
| `idt_opools.py` | 122 | IDT oPools list-price model |

Read for this document: `synthesis_vendors.py`, `optimizer.py`, `objectives.py`, `pareto.py`,
`ligation_fidelity.py`, `idt_opools.py`. The rest are mapped but not examined.

---

## Cross-module sharing — measured against them, and we come off worse

§12 asked for matched-condition comparisons. Breaches (under both rule sets) and codon
adaptation (over the aligned span) were done. Repetition and cross-module sharing had only ever
been measured for CLIPPR against itself, which is a self-assessment.

Measured 2026-09-17 by `validation/experiments/matched_oracle_benchmark.py`, one implementation
written in that file and applied to both systems' coding spans — deliberately not either side's
own scorer, because the reference carries a similarity penalty in its objective and CLIPPR
carries `collection_sharing` in `synthesis_fitness`, so using either would grade one system with
the term it was optimising for. k = 12, folded on reverse complement:

| | internal repetition | cross-module sharing | distinct 12-mers |
|---|---:|---:|---:|
| CLIPPR | 0.000000 | **0.964147** | 552 |
| GRASP Designer | 0.000000 | **0.911273** | **1005** |

Per module, CLIPPR shares more on **33 of 42**, less on 2, tied on 7.

**Neither system produces a module that repeats internally** — no duplicated 12-mer inside any
single sequence, on either side. That axis is settled and it is a tie.

It also corroborates an earlier result by a different route. `repeat_burden`, CLIPPR's original
third Pareto axis, was identically zero at k = 20, 16, 12 and 10 and registered only at k = 8,
which is why it was replaced: a repeated 10-mer inside a ~100 nt module is genuinely rare. This
measurement reaches the same conclusion from the reference's sequences as well as ours, so the
zero is a property of the molecules rather than of our scorer.

**On sharing across the collection the reference is ahead**, and the gap is larger than the
means suggest: their 42 designs draw on **1005 distinct 12-mers against our 552**. Nearly twice
the sequence diversity over the same 42 proteins.

**Why, mechanically — and I had this wrong the first time I wrote it down.**

My first version of this section said nothing in CLIPPR's per-module path sees collection
sharing. That is false. `library_search.optimise_library` carries a `sharing` term weighted
**1.0, equal to adaptation** — not an afterthought. What has no sharing term is
`recoding.recode_inventory`, and that is what produced `work/phaseb/inventory_recoded.json`, the
inventory this benchmark compares. So the comparison above pitted our *no-sharing* artefact
against their *with-similarity* one.

Running the fair version:

| inventory | cross-module sharing | distinct 12-mers |
|---|---:|---:|
| `inventory_recoded.json` — recoded per module, no sharing term | 0.964147 | 552 |
| `inventory_greedy.json` — our collection optimiser, sharing weighted 1.0 | 0.953602 | 621 |
| GRASP Designer | **0.911273** | **1005** |

Our collection optimiser helps and does not close the gap: it recovers about **a fifth** of the
distance on sharing (0.964 → 0.954 against their 0.911) and adds 69 distinct 12-mers against
their lead of 453.

Two differences plausibly explain the rest, neither of them measured here. Their similarity
penalty sits **inside each module's solve**, so every candidate codon choice is scored against
the rest of the collection as it is made; ours ranks whole proposals *after* a per-module
optimiser has already produced them, over a bounded proposal budget. And CLIPPR recodes each
module toward the same codon optimum by default, so identical amino-acid stretches receive
identical codons unless something actively pushes them apart.

**What it is worth, honestly.** Both figures are high because these modules are tandem repeats
of one template — the architecture, not either design. A 5-point difference on a quantity that
starts at 0.91 is not obviously a large practical change, and no synthesis outcome was measured:
this is a sequence property, and whether it changes an oligo pool's behaviour at the bench is
exactly the sort of claim this project does not make. But it is a real, matched, reproducible
difference on an axis where we had reported only our own number.

**The improvement it suggests**, not yet made: move the similarity term into the per-module
solve rather than scoring finished proposals, and give `recode_inventory` one at all. Both trade
codon adaptation for diversity at a rate nobody here has measured, so each needs its own
experiment — closing a gap is not by itself an improvement.
