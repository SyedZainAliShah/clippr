# Wrapper-aware optimisation pilot — result

**10 October 2026.** Predeclared experiment, run once, both arms. Nothing committed to a
default. **RAN** = executed here.

## The question

The one-shot route enforces the local GC band on the **coding span**. What a vendor receives
is the destination's entry overhang, plus the coding slice, plus the enzyme wrapper. All 200
saved corpus designs pass their coding-span checks and are refused on the ordered fragment.

> With the geometry frozen and only synonymous codon choices changed, can every ordered
> fragment satisfy the band the project already declares?

## Answer

**Yes, on this pilot.** 9 of 9 attempts produced an ordered fragment set that both the shipped
checker and the independent checker accept, with no rule relaxed.

| Arm | Accepted | Independent re-read |
|---|---|---|
| **Constraint applied** | **9 of 9** | exit 0 — 9 accepted, 0 rejected, 0 disagreements |
| **Constraint omitted** (negative control) | **0 of 9** | exit 0 — 0 accepted, 9 rejected, 0 disagreements |

The two arms differ in one thing only: whether the optimiser was shown the ordered fragment
or just the coding span. Same targets, same seeds, same frozen geometry, same policy, same
codon table. **This is why the result is attributable to the constraint** and not to a lucky
re-roll of the solver.

## What was declared before running

```
targets  AAAAUGUGG (9S), UUACACGUGCGUAC (14S), CUAUCACAUCACAUAAGCG (19S)
seeds    1, 42, 1234
budget   600 s per attempt, 5400 s total
policy   clippr-strict-legacy — local GC 0.35-0.65 over 50 nt, homopolymer <= 4
frozen   protein, genetic code, codon table, policy, fragment boundaries,
         overhang assignment, construction DNA
changed  synonymous codon choices only
```

Used: **339 s of 5400**. Per-attempt solve time 0.6–11.0 s.

## Cost

| | |
|---|---|
| worst CAI loss | **0.005** (0.603 → 0.598) |
| best CAI change | **+0.002** — two attempts improved |
| protein changed | **never** — translation re-checked against the baseline protein |
| geometry changed | **never** — flanks byte-identical, fragment boundaries identical |
| invariant violations across 9 attempts | **0** |

Codon adaptation is essentially unaffected. The breach was a narrow one — 0.660 against a
0.650 ceiling — and fixing it costs a handful of codon choices.

## How it works

A DNA Chisel constraint that evaluates GC windows over `left_flank + coding_slice +
right_flank` for every fragment, with breaching windows reported back in coding-span
coordinates clipped to the mutable region. The flanks are **derived from the baseline's own
oligo table**, not rebuilt from the wrapper rules, so the pilot cannot disagree with
production about where the payload sits; a slice that does not occur exactly once in the
ordered string is refused rather than guessed.

It is evaluated globally, not localised: one codon change moves a window in a fragment whose
flanks lie outside any local view.

## What this does not establish

- **It is not a default change.** The constraint is passed through
  `optimize_cds(extra_constraints=...)`, which is empty for every shipped route. Promoting it
  needs the predeclared broader panel, reported in full — all targets, not only recoveries.
- **Three targets and three seeds is a feasibility pilot**, not a performance estimate over
  arbitrary targets. The 19-mer took 11 s at its worst here; that does not bound the general
  case.
- **No vendor approved anything.** A passing engineering-policy check is not a vendor
  acceptance. No order was submitted.
- **Nothing biological was measured.** Not binding, not expression, not assembly yield. No
  stage score was combined into a yield figure.
- **The saved corpus is unchanged.** Its 200 designs remain rejected and were not
  regenerated; this pilot designed fresh candidates for three targets.

## Evidence

In `work/repair_release_20261010-232032/pilot/`:

| File | What it is |
|---|---|
| `wrapper_aware_pilot.json` | treatment arm, all nine attempts with baseline, candidate, flanks, timings |
| `wrapper_aware_pilot_candidates.jsonl` | the exported candidate fragments, as written |
| `wrapper_aware_pilot_independent.json` | the independent checker's verdict on those bytes |
| `negative_control.json` | the constraint-omitted arm |
| `negative_control_candidates.jsonl`, `negative_control_independent.json` | and its independent verdict |

## Reproduction

```powershell
Set-Location -LiteralPath 'C:\CodePractice\Drylab\clippr'
$Py = '.\.venv\Scripts\python.exe'
$Out = "work\pilot_$(Get-Date -Format yyyyMMdd-HHmmss)"

& $Py validation\experiments\wrapper_aware_pilot.py --out "$Out\treatment.json"
& $Py validation\experiments\wrapper_aware_pilot.py --out "$Out\control.json" --without-constraint
```

The control must report **0 accepted**. If it does not, the treatment arm proves nothing and
the result is void.

## Next, if this is taken further

The predeclared broader comparison: all 200 corpus targets under the same frozen-geometry
rules, reporting every target rather than only the recoveries, with runtime measured. Only
then is a default change a question worth asking.
