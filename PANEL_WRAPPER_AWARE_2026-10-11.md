# Wrapper-aware optimisation over the full 200-design panel

**11 October 2026.** The predeclared broader comparison the pilot called for. Both arms, every
design reported. Nothing switched on.

## The question

The pilot recovered 9 of 9 on three targets. That is a feasibility answer, not a panel answer.
This runs all 200 pinned designs:

> With the geometry frozen and only synonymous codon choices changed, how many of the 200
> designs can satisfy the declared local GC band on **every** ordered fragment?

## Answer

**192 of 200.** Eight fail, for one specific and diagnosable reason.

| | treatment | control (constraint omitted) |
|---|---:|---:|
| accepted | **192 / 200** | **0 / 200** |
| rejected | 8 | 200 |
| hard findings across all fragments | 277 → **8** | 277 → 277 |
| reproduced the corpus sequence bit-for-bit | 0 | **200 / 200** |
| CAI delta | median −0.003 | **exactly 0.0, all 200** |
| invariant violations | **0** | 0 |
| independent re-read of the written export | 192 accepted, 8 rejected, 0 disagreements | 0 accepted, 200 rejected, 0 disagreements |
| elapsed | 12.9 min | 3.7 min |

**The control is also the setup check, and it passed exactly.** With the constraint removed,
the same seed, protein, codon table and locked overhangs returned the corpus coding sequence
**bit-for-bit on all 200 designs**, with a CAI delta of exactly zero and the identical 277 hard
findings. So this harness is indistinguishable from whatever generated the corpus, and the 192
recoveries are attributable to the constraint and to nothing else.

## By architecture

| | accepted | failed |
|---|---:|---:|
| 9S (100 designs) | **100** | 0 |
| 14S (50 designs) | 44 | 6 |
| 19S (50 designs) | 48 | 2 |

The 19-mers were expected to be hardest — most fragments, most locked junctions, least codon
freedom. **They were not.** The 14-mers failed three times as often. That is unexplained and is
not claimed to be understood.

## Cost, over the 192 that worked

```
CAI delta    median -0.003    worst -0.032    best +0.006
             10 improved, 181 slightly worse, only 4 worse than 0.02
solve time   median 1.9 s     p95 15.2 s      max 50.1 s
```

Predicted ligation fidelity is **unchanged by construction** — the geometry is frozen, so the
overhangs and therefore the fidelity cannot move. It is not reported as a gain.

## The eight failures

Every one breaches the **floor**, not the ceiling:

```
ACAGUGCGACAUUA_F4 : GC 0.340 outside (0.35, 0.65)
GGCCAAUUCGCAUG_F4 : GC 0.340 outside (0.35, 0.65)
GCAAACAAAAGCAU_F4 : GC 0.340 outside (0.35, 0.65)
GGCCCGGCCGCGAC_F6 : GC 0.340 outside (0.35, 0.65)
UUCUCUCACGUUGU_F4 : GC 0.340 outside (0.35, 0.65)
CGGCAAUUAUGAAA_F4 : GC 0.340 outside (0.35, 0.65)
UUAAUCCAUAAAACACUAG_F5 : GC 0.340 outside (0.35, 0.65)
ACUCGGUAUCGAGCCUAAA_F7 : GC 0.340 outside (0.35, 0.65)
```

**This is a new failure mode, not the old one.** Every breach in the baseline corpus was the
ceiling — 277 fragments above 0.65, none near the floor. Pushing GC down to clear the ceiling
overshoots below 0.35 on one interior fragment. On these eight targets there is a window where
no synonymous choice satisfies both ends of the band at once.

Three properties make this a clean limitation rather than a hazard:

- **The solver reported it.** All eight carry `constraints_ok=False`.
- **It gave up immediately** — 0.08 s at worst, against up to 50 s on designs it solved. An
  instant contradiction, not an exhausted search.
- **Nothing orderable was written.** The final check rejected all eight and the independent
  checker agreed, with zero disagreements.

The failing fragments are F4, F5, F6 and F7 — all interior, never F1. The pilot's baseline
failures were F1. Different position entirely, so the pilot gave no warning of this.

The CAI of those eight collapses to about 0.18, which is near-random. That is the **unsolved**
sequence being handed back, not a trade-off anyone chose. The first version of this report's
summary folded those into the CAI distribution and printed a worst-case cost of −0.41; that
number described a solver failure, not a price any accepted design paid. The summary is
recomputed from the same untouched per-design records, with CAI scoped to the accepted set and
failures reported separately. Both files carry a `summary_recomputed` note saying so.

## What this does and does not support

**Supports:** showing the solver the ordered fragment instead of the coding span recovers 96%
of a panel that is currently 0% orderable, at a negligible codon cost, with no rule relaxed and
no invariant lost. There is now a positive final-export control under the strict policy, over
the declared panel, independently verified.

**Does not support:**

- **It is still not switched on.** `optimize_cds(extra_constraints=...)` is empty for every
  shipped route and a test asserts no `src/` module other than its definition references it.
  Turning it on is a decision with consequences for every published number and is not taken
  here.
- **8 of 200 remain unsolvable** by this approach. A route that refuses 4% of targets is a
  different product claim from one that refuses none, and the eight are named above rather
  than rounded away.
- **One seed.** Seed 42 throughout, matching the corpus baselines so the comparison is
  like-for-like. Whether another seed rescues some of the eight is untested.
- **Nothing biological.** No binding, expression, regulation or assembly yield was measured. No
  stage score was combined into a yield.
- **No vendor involvement.** A passing engineering-policy check is not vendor acceptance. No
  order was submitted.
- **The band itself still has no recorded provenance.** These designs satisfy a constant whose
  origin this repository has never recorded. See `docs/gc_band_and_constraint_model.md`.

## Evidence

In `work/repair_release_20261010-232032/panel/`:

| File | What it is |
|---|---|
| `treatment.json` | all 200 designs, baseline and candidate, per-design timings and findings |
| `treatment_candidates.jsonl` | the 200 exported candidate fragment sets, as written |
| `treatment_independent.json` | the independent checker's verdict on those bytes |
| `control.json`, `control_candidates.jsonl`, `control_independent.json` | the constraint-omitted arm |

## Reproduction

```powershell
Set-Location -LiteralPath 'C:\CodePractice\Drylab\clippr'
$Py = '.\.venv\Scripts\python.exe'
$Out = "work\panel_$(Get-Date -Format yyyyMMdd-HHmmss)"

& $Py validation\experiments\wrapper_aware_panel.py --out "$Out\treatment.json"
& $Py validation\experiments\wrapper_aware_panel.py --out "$Out\control.json" --without-constraint
```

The control must report **0 accepted and 200 reproducing the corpus bit-for-bit.** If it
reproduces fewer, the harness differs from whatever generated the corpus and the treatment arm
proves nothing.

## The decision this leaves

Three options, and the evidence now exists to argue about them:

1. **Keep the strict band, keep the current default.** The de novo route stays unusable and the
   rejection is published. No new risk; no new capability.
2. **Switch the wrapper-aware constraint on.** 192 of 200 become orderable at a negligible
   codon cost. Costs: eight named targets still refuse, every published number moves, and the
   constraint needs its own regression coverage in the shipped path.
3. **Keep the band and investigate the eight.** Whether a different seed, a different split or
   a relaxed interior geometry recovers them is unknown and is a smaller experiment than
   option 2.

This document does not choose. Option 2 is the one that changes a default, and a default change
is not a dry-lab result.
