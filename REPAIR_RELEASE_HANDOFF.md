# CLIPPR repair release — executor handoff, revision 4

**Executor: Claude. Orchestrator: Codex. Revised 10 October 2026.**

Response to `repair_plan_20261009` and to three independent reviews: `repair_review_20261009`,
`repair_review_20261010_revision2` and `repair_review_20261010_revision3`. Nothing has been
committed, pushed, published, ordered or sent. All changes are local and uncommitted on
`479690b53663c397926aa431e72f4965dfbe4998`.

The third review accepted the export predicate, the real notebook execution, the
missing-input refusal, the missing-data verifier exit and the state-table repair. It found
three narrower issues and was right about all three. They are closed below.

| Revision 3 claimed | Actually |
|---|---|
| "All five findings closed" for the display finding | The **API** boundary was closed. Every surface a person reads still led with a green coding-span `PASS`, and the audit never mentioned the rejection |
| check 13 "binds the release" | It required fields to be non-empty without recomputing them. Six altered-record controls passed, and a bundle of four files passed with no B–E or V5–V7 evidence |
| "R10 delivered" | The dated successor exists. Four of its claims exceeded the evidence, and the R10 requirement set had no deliverable map |

---

## 1. Current measured state

Evidence bundle: `work/repair_release_20261010-232032/`. All three earlier bundles are
preserved unchanged, including the `before/` snapshot of the pre-repair records.

| | |
|---|---|
| package source | `0f4d63540cd2c291e6c85185b8d2c58232b0d300f2df37060e9f465135bb14ed` |
| wheel | `db2a77070c628ff9…`, built from that source |
| full suite | **1108 passed of 1108 collected, 0 failed**, exit 0 (was 1062) |
| focused suite | **132 passed, 0 failures, 0 errors, 0 skipped**, 40.733 s, JUnit in the bundle |
| aggregate release check | **13 of 13 passed, 0 failed, 0 needing a human**, exit 0 |
| clean-environment install | **PASS**, reproducing the **new** summary verbatim |
| Gate F | **PASS**, 12 artefacts re-verified after relocation |
| Gates B, C, D, E | **PASS** |
| V5, V6, V7 | **PASS**, records in the bundle |
| designer notebook | 14 ran, 0 failed, **8 postconditions**, **2 inputs staged and hashed** |
| inventory notebook | 11 ran, 0 failed, **8 postconditions**, **2 inputs staged and hashed** |
| corpus contract | 200 designs, 0 accepted, **200 rejected**, 0 not assessed, 200 comparisons, **0 disagreements** |
| fragments | 1050 checked, **277 with a hard finding**, 0 unusable |
| frozen record | `release_check.json` in the bundle, sha256 `d6e40803d2ce45b3…` |

**The bundle is now self-contained and immutable.** Eleven records plus a frozen copy of the
aggregate with its hash beside it. `work/release/` and `work/independent/` remain as mutable
conveniences, which is what the renderer incident exposed.

---

## 2. The three findings of the third review

### Finding 1 — status propagation incomplete in the designer display. **Closed.**

The card's badge text and colour came from `result["qc"]["status"]`, so a refused design
rendered a green `PASS`; the audit printed `constraints satisfied` and `QC PASS` with no final
verdict; and `result["summary"]` — the README quickstart's entire printed output — ended
`QC PASS.`

The audit now reads:

```
  SYNTHESIS      NOT ORDERABLE
                 - the final-artifact status is 'fail', not one of pass, warning
                 - 1 hard policy breach(es) on the ordered sequences: AAAAUGUGG_F1:
                   GC 0.660 outside (0.35, 0.65) in the window at 1
  final artifacts    fail  (the sequences that would be ordered)
  optimiser      constraints satisfied  (the coding span)
  coding-span QC PASS  (the CDS only, not the ordered sequences)
```

The card's badge is `NOT ORDERABLE` in the refusal colour, with a panel naming each reason
and stating that the files are not an order; `coding-span QC` and `final artifacts` are
separate scoped statistics; GC and repeats say `coding span`. The fragment section no longer
says "fragments to order" or "orderable piece" and states that a row is not an order. The
download table labels diagnostic keys `DIAGNOSTIC ONLY - rejected …`.

The summary, and therefore the README quickstart, is now:

```
AAAAUGUGG (9S) -> 302 aa, 906 nt, 4 fragments. Fidelity 0.828.
Coding-span QC PASS. NOT ORDERABLE (final artifacts: fail).
```

Release check 3 parses that comment from the README as its expected clean-install output, so
it self-updated and **passed on the new string**.

The retained CDS pass is kept, not deleted. It is true, and scoping it is the honest fix.

**The display-level falsifier is asserted twice.** `tests/test_display_verdict.py` (14
controls) executes the generated card cell against a real refused design and checks the badge,
the colour, the breach figure, the diagnostic wording and the scoped CDS statistics, with a
positive case proving the card is not hardwired to refuse. And the designer notebook's own
postconditions now include `'NOT ORDERABLE' in card` and
`'NOT ORDERABLE' in result['audit'].report()`, so the release gate fails if the display
regresses. Neither is `bool(audit.report())`.

### Finding 2 — check 13 checked presence, not identity. **Closed.**

All six of the review's counterexamples now fail, as do five more:

| Control | Was | Now |
|---|---|---|
| `notebook_sha256` = 64 zeros | `ok=True` | FAIL, "records notebook 000000000000 but … hashes to …" |
| staged Table S1 hash = zeros | `ok=True` | FAIL, "its source now hashes to …" |
| Gate F source + wheel hash = zeros | `ok=True` | FAIL, "gate_f built 000000000000, but this checkout is …" |
| `corpus_sha256` = zeros | `ok=True` | FAIL, corpus hash mismatch |
| `expectations=[]`, `postconditions_declared=8` | `ok=True` | FAIL, "claims 8 postcondition(s) but lists 0" |
| 200 rows, all copies of row 1 | `ok=True` | FAIL, "not the 200 declared panel rows" |
| four-record directory, no B–E or V5–V7 | `ok=True` | FAIL, "missing from …" |
| `--expect "1"` as the only postcondition | would have passed | FAIL, "declares no postcondition covering 'deposited'" |
| a tampered validator record | not inspected | FAIL, "a different version of v6_constraints.py" |
| a wheel edited after the record | not inspected | FAIL, "the file hashes to …" |
| a wheel path that does not exist | not inspected | FAIL, "no retained wheel to verify its hash against" |

`REQUIRED_RECORDS` names eleven records; every stated identity is recomputed from the file it
names; `NOTEBOOK_CONTRACT` names the result each notebook must have asserted; and the corpus
panel is compared as a **multiset**, so a duplicated panel fails while a legitimately repeated
target does not.

The corpus report and V5–V7 now carry `provenance` with the validator's own SHA-256, because
the package fingerprint does not cover `validation/*.py`. The three validators gained `--out`
so the bundle holds its own copies, and `check_independent_agreement` reads the bundle when one
is given.

**The review's criticism of my own positive fixture is accepted.** It built a Gate F record
holding only `{"gate": "PASS"}`, a zeroed notebook hash and arbitrary corpus rows, and called
that the passing case — it blessed a fixture the gate should reject. The positive control is
now **the actual release bundle**, and every negative control is a *copy* of it with one thing
changed. 34 controls in `tests/test_release_evidence_binding.py`, including one that re-hashes
the live bundle to prove no control wrote into it.

Adverse controls run only in `tmp_path`. The bundle is copied, never modified.

### Finding 3 — successor claims and the R10 map. **Closed, with one item open.**

The four overstatements are corrected, and §10 of `RELEASE_STATUS_2026-10-10.md` states what
revision 1 claimed beside what is correct, so anyone who read it knows which claims to drop:

1. **"every current design is rejected"** → "All 200 saved corpus designs were rejected on
   revalidation under this policy, and the demonstrated default one-shot example is also
   rejected." The corpus is the preserved historical panel, **not** newly generated.
2. **V6's scope** → "200/200 satisfy the V6 coding-sequence constraints; 200/200 fail the
   separate final-oligo policy check." Both are true. Release check 5's own printed sentence
   was also unscoped and now says `the V6 coding-sequence constraints`.
3. **"every figure is from the evidence directory"** → the bundle holds its own copy of every
   record, including a frozen `release_check.json` with its hash; the shared paths are named
   as mutable conveniences.
4. **the declared data route** → the review is right that the designer record staged nothing
   while the run read a codon table and a chloroplast genome from a per-user cache. The
   harness gained `--env`, both notebooks now stage every substantive input from a
   repository-relative path, and the hashes are in the records.

**R10 deliverable map:** §9 of the successor maps each original requirement to a current
section or file and marks its state. One item is **open and not closed by this pass**: the
attribution and input-rights statement for Table S1 and the deposited kit. Verification of the
current official publication requirements is also explicitly **not** certified. R10 is
therefore not wholly closed, and the successor says so rather than claiming it is.

---

## 3. The GC rejection: unchanged, and still the open decision

200 of 200 saved designs fail the strict final-oligo rule. 1050 fragments checked, **277 with
a breach**; 268 worst windows at 0.66 and 9 at 0.68 against a 0.65 ceiling. These agree with
all three independent counts.

Cause, traced: the one-shot optimiser enforces GC on the **coding span**. The ordered sequence
is the destination's entry overhang, plus the CDS slice, plus the enzyme wrapper. On the
shipped example the CDS 50-nt windows run 0.540–0.640, inside the band, while the ordered
oligo reaches 0.660 at offset 1. The payload begins `CTCA`, which no CDS window covered.

The threshold was not weakened and no optimiser default was changed. **77 of the failing
fragments are F2–F7** (F2=5, F3=22, F4=24, F5=2, F6=12, F7=12), so an outer-gene-only fix
cannot settle it.

### The pilot has now been run, and it succeeded

Both arms, 3 targets × 3 seeds, declared before execution. Full report:
`PILOT_WRAPPER_AWARE_2026-10-10.md`.

| Arm | Accepted | Independent re-read of the written export |
|---|---|---|
| wrapper-aware constraint applied | **9 of 9** | exit 0 — 9 accepted, 0 rejected, 0 disagreements |
| constraint omitted (negative control) | **0 of 9** | exit 0 — 0 accepted, 9 rejected, 0 disagreements |

The arms differ in one thing: whether the solver saw the ordered fragment or only the coding
span. Same targets, seeds, frozen geometry, policy and codon table. **The negative control is
what makes the result attributable** — without it, 9 of 9 on a first attempt would be
indistinguishable from a lucky re-roll.

Worst CAI cost **0.005**; two attempts improved. Protein never changed, fragment boundaries
and construction DNA byte-identical, 0 invariant violations across 9 attempts. 339 s of a
declared 5400 s budget.

**The broader panel has now been run too.** All 200 pinned designs, both arms —
`PANEL_WRAPPER_AWARE_2026-10-11.md`:

| | treatment | control |
|---|---:|---:|
| accepted | **192 / 200** | **0 / 200** |
| hard findings | 277 → **8** | 277 → 277 |
| reproduced the corpus bit-for-bit | 0 | **200 / 200** |
| independent re-read | 192 acc, 8 rej, 0 disagreements | 0 acc, 200 rej, 0 disagreements |

9S 100/100, 19S 48/50, 14S 44/50 — the 19-mers were expected to be hardest and were not; that
is unexplained. CAI cost over the accepted set: median −0.003, worst −0.032, 10 improved.

The eight failures are a **new** mode: every one breaches the 0.35 **floor** at 0.340, where
every baseline breach was the 0.65 ceiling. Pushing GC down to clear the ceiling overshoots
below the floor on one interior fragment. All eight carry `constraints_ok=False`, gave up in
≤0.08 s, and wrote nothing orderable.

**This still changes nothing about the shipped default.** `optimize_cds(extra_constraints=...)`
is empty for every shipped route and a test asserts no `src/` module other than its definition
references it. Switching it on would move every published number and still refuse eight named
targets; that decision is not taken here and is not a dry-lab result.

One seed throughout (42, matching the corpus baselines, so the comparison is like-for-like).
Whether another seed rescues any of the eight is untested.

---

## 4. The demonstration

Both routes are verified inside the notebooks, on the wheel in §1, with every input staged and
hashed.

The one-shot route is the **rejection** example, and the notebook now says so where a reader
looks first:

```
result card badge              : NOT ORDERABLE
AAAAUGUGG final-artifact status : fail
synthesis_ready                : False
hard  : AAAAUGUGG_F1: GC 0.660 outside (0.35, 0.65) in the window at 1
paths : diagnostic_oligo_csv present, oligo_csv absent
```

The inventory route is the usable demonstration:

```
42 assembly-ready substrates from selected-recommended@a6c3074572de546c
order plan : 42 oligos, eligible, 1 pool(s)
package    : 8 tasks, ok=True, 0 failure(s)
vendor approval: not obtained; local checks against a transcribed vendor profile,
                 and no order was placed
```

The wording for a preliminary advice session is in §4 and §5 of
`RELEASE_STATUS_2026-10-10.md`. The review's instruction is satisfied without a manual
annotation: the card itself now carries the rejection, so the one-shot example can be shown as
a diagnostic without someone having to explain the badge away.

---

## 5. What remains unassessed

- **Whether any vendor would refuse these sequences.** A failed engineering-policy check is
  not a vendor refusal. No order was submitted and no quote obtained. The vendor profile is a
  local transcription dated 2026-09-15 modelling only the 50 pmol scale.
- **Whether the wrapper-aware route should be the default.** A positive final-export control
  now exists over the whole declared panel — 192 of 200, independently rechecked — but the
  **shipped** route still produces nothing orderable, eight targets resist the approach
  entirely, and only one seed was tried. Switching the default is an open decision, stated with
  its three options in `PANEL_WRAPPER_AWARE_2026-10-11.md`.
- **Browser or Colab execution.** Only installed local notebook execution was verified.
- **Level-1 assembly context.** Participants are not established from a flat role list, so no
  level-1 score is reported.
- **The wiki-facing attribution text.** `NOTICE.md` already carries the software's
  third-party material and input rights, including the deliberate non-redistribution of the
  GRASP insert sequences and the unlicensed PPRmatcher table; `data/` is gitignored and absent
  from the wheel. An earlier draft of this handoff said the statement was unwritten. That was
  wrong, and it is corrected in §9 of the successor. A wiki page remains a separate
  deliverable.
- **Current official publication requirements.** Not verified.
- **The exact historical bytes behind the `c76dc3ca…` fingerprint.** The line-ending mechanism
  is demonstrated; the historical mix is not recovered.
- **Any biological property.** Nothing here measures binding, expression, regulation,
  specificity or assembly yield.

## 6. Reproduction

§8 of `RELEASE_STATUS_2026-10-10.md`, which is self-contained, writes into a fresh directory,
and stages every substantive input from a repository-relative path. The acceptance gate is the
last command there: `validation\state_table.py` run on its own must exit 0. It does.

The full suite was **not** run manually and then again inside the aggregate. Check 1 is the
measurement; the separate focused run exists for its JUnit record and covers a different, named
subset.

## 7. Scope of this pass

Changed: project code (`design.py`, `audit.py`, `library.py`), the notebook generator, both
generated notebooks, the README quickstart, three independent validators, the notebook harness,
the release check, and four documents. Added three test modules.

Did **not**: commit, push, publish, submit an order, or message anyone. 39 files are modified
or untracked locally.

Also corrected here: the three inputs the notebooks stage are **not** tracked in the
repository. `data/` is gitignored, so a fresh clone must fetch the codon table and genome and
download Table S1 before the documented commands run. Their hashes are now listed in §8 of the
successor so a reproducer can confirm identical bytes.
