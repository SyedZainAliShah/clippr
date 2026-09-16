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

Its scope overlaps CLIPPR's design and recoding routes. It has no concept of a *reusable
versioned inventory*: it redesigns all 42 modules from amino acids on every run.

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
marked `qc_status: WARNING`. Ours refuses and returns the module *unchanged*, carrying its
original poorly adapted sequence. Our "unchanged" count is not an absence of opportunity; it is
the band vetoing every improvement.

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

And three distinct categories where we have one:

- `synthesis` — design targets, optimised against softly
- `machine_hard_constraints` — what is actually rejected (Twist's real homopolymer rule is
  < 14, not ≤ 3)
- `manual_review_rules` — vendor or human judgement, e.g. "no CcdB"

Our single band is tighter than every row above, including both "conservative" profiles, and we
enforce **no global band at all** — the opposite emphasis to every vendor document cited.

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
- **Versioned inventories.** Content-hashed, reload-reproducing, recompilable. They redesign
  from amino acids every run; there is no artefact to carry forward.
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

**The CAI difference is a constraint difference, not an optimiser difference. [E]** Same
solver, same seeds, same 12 modules, only the band changed:

| band | mean CAI |
|---|---:|
| CLIPPR — 0.35–0.65 / 50 nt, homopolymer ≤ 4 | 0.692434 |
| reference — 0.15–0.85 / 50 nt, homopolymer ≤ 3 | **0.854766** |

The band is worth ≈ 0.16 CAI; the observed gap is 0.037. Both systems maximise the same
quantity and neither is cleverer at it. Full argument: `docs/gc_band_and_constraint_model.md`.

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
