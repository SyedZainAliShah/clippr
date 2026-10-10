# Wrapper-aware optimisation over the full 200-design panel

**11 October 2026, revision 2.** The predeclared broader comparison the pilot called for.
Both arms, every design reported. Nothing switched on.

> **Revision 1 of this document reported 192 of 200 and described a new failure mode in
> which eight targets could not satisfy both ends of the GC band.** That was wrong. The
> eight failures were a defect in the experiment's own constraint, not a property of the
> targets. Corrected result and the defect are in §5. Revision 1 should not be quoted.

## The question

The pilot recovered 9 of 9 on three targets — a feasibility answer, not a panel answer. This
runs all 200 pinned designs:

> With the geometry frozen and only synonymous codon choices changed, how many of the 200
> designs can satisfy the declared local GC band on **every** ordered fragment?

## 1. Answer

**200 of 200.**

| | treatment | control (constraint omitted) |
|---|---:|---:|
| accepted | **200 / 200** | **0 / 200** |
| rejected | 0 | 200 |
| hard findings across all fragments | 277 → **0** | 277 → 277 |
| reproduced the corpus sequence bit-for-bit | 0 | **200 / 200** |
| CAI delta | median −0.003 | **exactly 0.0, all 200** |
| invariant violations | **0** | 0 |
| independent re-read of the written export | **200 accepted, 0 rejected, 0 disagreements** | 0 accepted, 200 rejected, 0 disagreements |
| elapsed | 13.5 min | 4.0 min |

Every architecture, complete: **9S 100/100, 14S 50/50, 19S 50/50.**

## 2. Why the control is the load-bearing part

With the constraint removed — same seed, protein, codon table and locked overhangs — the
harness returned the corpus coding sequence **bit-for-bit on all 200 designs**, with a CAI
delta of exactly zero and the identical 277 hard findings.

So this harness is indistinguishable from whatever generated the corpus, and the recoveries
are attributable to the constraint and to nothing else. Without that check, 200 of 200 on a
first attempt would be indistinguishable from a harness that quietly differed from the one
under test.

## 3. Cost

```
CAI delta    median -0.003    worst -0.032    best +0.006
             13 improved, only 4 worse than 0.02
solve time   median 2.1 s     p95 14.3 s      max 51.7 s
```

Predicted ligation fidelity is **unchanged by construction** — the geometry is frozen, so the
overhangs and therefore the fidelity cannot move. It is not reported as a gain.

## 4. How it works

A DNA Chisel constraint that evaluates GC windows over `left_flank + coding_slice +
right_flank` for every fragment, with breaching windows reported back in coding-span
coordinates clipped to the mutable region. The flanks are derived from each design's own
stored oligo table, not rebuilt from the wrapper rules, so the experiment cannot disagree with
production about where the payload sits; a slice that does not occur exactly once in the
ordered string is refused rather than guessed.

It is evaluated globally, not localised: one codon change moves a window in a fragment whose
flanks lie outside any local view.

## 5. The defect in revision 1, and what it fabricated

Revision 1 reported 8 failures and built a narrative on them: that pushing GC down to clear
the 0.65 ceiling overshoots below the 0.35 floor, that all eight breached at 0.340 on interior
fragments F4–F7, and that no synonymous choice could satisfy both ends on those targets.

**None of that was real.** The constraint's `localized` method was declared as

```python
def localized(self, location, problem=None):
```

while some of DNA Chisel's constraint-resolution paths call it with a further keyword:

```
TypeError: EnforceGCInOrderedFragments.localized()
           got an unexpected keyword argument 'with_righthand'
```

On those eight designs the solver took that path, the call raised, the solver reported
`constraints_ok=False` and handed back the **unoptimised** reverse-translated sequence. The
"0.340 floor breach" was that raw sequence's GC. The CAI collapse to ~0.18 was the same
artifact. The clue was in the record the whole time: those eight gave up in ≤0.08 s while
solved designs took up to 51 s — an instant crash, not an exhausted search.

Adding `with_righthand=True` to the signature fixed all eight. Re-tested individually first
(8 of 8 accepted), then the whole panel re-run, because the same broken path was reachable for
every design and not only the eight that surfaced it.

Two wrong claims were published from this experiment before it was right:

| Claimed | Actually |
|---|---|
| "worst CAI cost −0.41" | solver failures averaged into the cost distribution; true worst is **−0.032** |
| "8 of 200 fail; a new floor-breach failure mode" | an exception in the experiment's own code; **200 of 200** succeed |

Both are recorded here rather than quietly replaced. The per-design records in the bundle are
the raw evidence for the corrected numbers.

## 6. What this does and does not support

**Supports:** showing the solver the ordered fragment instead of the coding span recovers the
**entire** declared panel, at a negligible codon cost, with no rule relaxed, no invariant lost
and no geometry changed. There is now a positive final-export control under the strict policy
over all 200 designs, independently verified from the written bytes.

**Does not support:**

- **It is still not switched on.** `optimize_cds(extra_constraints=...)` is empty for every
  shipped route and a test asserts no `src/` module other than its definition references it.
  Switching it on moves every published number and needs its own regression coverage in the
  shipped path. That decision is not taken here.
- **One seed.** Seed 42 throughout, matching the corpus baselines so the comparison is
  like-for-like. Behaviour at other seeds is untested.
- **One host and one policy.** *Chlamydomonas* nuclear codon table, `clippr-strict-legacy`,
  `igem_rfc1000` enzymes, level-0 destination, BsaI. Nothing about other hosts or profiles.
- **Random targets.** The 200 are uniformly random RNA strings generated by this project, not
  biologically chosen targets.
- **Nothing biological.** No binding, expression, regulation or assembly yield was measured.
  No stage score was combined into a yield.
- **No vendor involvement.** A passing engineering-policy check is not vendor acceptance. No
  order was submitted. Twist's own published acceptance criteria are substantially looser than
  the band satisfied here — see `docs/synthesis_policy.md`.

## 7. Evidence

In `work/repair_release_20261010-232032/panel/`:

| File | What it is |
|---|---|
| `treatment_v2.json` | the corrected run: all 200 designs, baseline and candidate, per-design timings and findings |
| `treatment_v2_candidates.jsonl` | the 200 exported candidate fragment sets, as written |
| `treatment_v2_independent.json` | the independent checker's verdict on those bytes |
| `control_v2.json`, `control_v2_candidates.jsonl`, `control_v2_independent.json` | the constraint-omitted arm |
| `treatment.json`, `control.json` | **revision 1, retained.** The run with the defect, kept as the record of what was wrong |

## 8. Reproduction

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

## 9. The decision this leaves

1. **Keep the strict band and the current default.** The de novo route stays unusable and the
   rejection is published. No new risk, no new capability.
2. **Switch the wrapper-aware constraint on.** All 200 become orderable at a negligible codon
   cost. Costs: every published number moves, and the constraint needs regression coverage in
   the shipped path rather than only in an experiment.
3. **Revisit the band itself**, now that its provenance is settled: it is Twist's published
   *codon-optimisation guidance*, while Twist's own stated acceptance criteria for a 50 bp
   window are 10–90%. That is a separate decision from this experiment and should not be
   folded into it.

This document does not choose. Option 2 changes a default, and a default change is not a
dry-lab result.
