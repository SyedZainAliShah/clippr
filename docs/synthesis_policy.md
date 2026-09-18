# The synthesis policy, and why it is an object rather than a constant

*Implemented 2026-09-16, in response to an independent review. The strict profile is the
default and reproduces every published number; nothing here changes a shipped result.*

---

## The two failures this exists to prevent

### A threshold living in two places

`recoding.GC_BAND` was the **validator's** band. `codons.optimize_cds` carried its own
`gc_bounds=(0.35, 0.65)` default, and `recode_inventory` never passed one. The two numbers were
equal, so nothing looked wrong — until someone tried to change one.

I claimed in an earlier draft that widening the band was "a one-line change". A review tested
it. Patching the constant alone gives a four-seed mean of **0.657195**, not the **0.819921** a
genuinely broadened policy reaches, because the solver goes on searching the narrow space while
the validator accepts a wider one. The claim was false, and the reason it was false is the real
finding: **a threshold that appears twice is not a constant, it is an undocumented policy.**

A profile now supplies `solver_bounds()` to the optimiser *and* the limits to the validator, so
the two cannot drift apart again. The regression asserts both numbers, and the broad figure is
the one that matters — the validator alone cannot produce it.

### An advisory rule enforced as a rejection

Twist publishes, for codon optimisation: *"global GC% of less than 25% or more than 65% and
local GC windows (50 bp) of less than 35% or more than 65%."* That is CLIPPR's band exactly.
But Twist frames it as guidance applied **during optimisation**, and separates it from the
sequence-acceptance rules that decide whether an order is taken.

CLIPPR enforced it as a hard rejection. Measured cost: about **0.16 CAI**, and 3 of 42 modules
returned unchanged rather than improved. The inverse error is just as bad — a soft warning read
as vendor approval — so the two can never be allowed to share a return value.

---

## What a profile carries

A number cannot be audited. `0.65` does not say where it came from, which product it applies
to, when it was read, or whether breaching it stops an order or costs a point.

| field | why |
|---|---|
| `source`, `product`, `read_on` | a vendor *name* is not a policy; the product and the date are |
| `local_gc`, `window`, `global_gc`, `max_homopolymer` | the thresholds |
| `enforcement` | per rule: `hard`, `target`, or `vendor_judgement` |
| `unresolved` | what this profile does **not** decide, so it is not mistaken for checked |
| `version` | content hash over every verdict-changing field |

**The three-way split is the point.**

    hard              a breach is infeasible; the candidate is refused
    target            a breach is reported and delivered; an optimiser steers away
    vendor_judgement  not machine-checkable here; recorded so it is not mistaken for checked

An unlisted rule is **hard**. Silence must not relax a constraint.

**No profile can soften a forbidden site.** An enzyme cutting where it should not is never a
matter of preference, and a profile that declares otherwise is ignored on that rule alone.
There is a test.

---

## The shipped profiles

| | `clippr-strict-legacy` | `twist-codon-optimisation-guidance` | `broad-experimental` |
|---|---|---|---|
| local GC / 50 nt | 0.35–0.65 **hard** | 0.35–0.65 **target** | 0.15–0.85 **hard** |
| global GC | none | 0.25–0.65 **target** | none |
| homopolymer | ≤ 4 **hard** | ≤ 4 **target** | ≤ 4 **hard** |
| source | **unrecorded** | Twist gene resources, read 2026-09-16 | CLIPPR experiment |
| for ordering? | yes — the default | yes, with warnings | **no** |
| mean CAI, 42 modules, 4 seeds | **0.618830** | — | **0.819921** |

`STRICT_LEGACY` is the default. Its band matches Twist's published guidance exactly — and its
`source` still says *unrecorded*, because this repository never recorded that as its origin.
Matching a guideline after the fact is not being derived from it, and back-dating a citation
would be the same class of error the rest of this project spent a week correcting.

`BROAD_EXPERIMENTAL` is for measuring what the strict band costs. It is **not** a vendor rule
and no vendor has been shown to accept sequences designed under it.

---

## The four-arm benchmark — RAN

`validation/experiments/policy_benchmark.py`, all 42 modules, four seeds, matched budgets,
every arm judged on its full ordered substrate under **its own** policy. Artefact:
`work/policy/policy_benchmark.json`.

| arm | mean CAI | improved | worse | hard breaches | advisory |
|---|---:|---:|---:|---:|---:|
| `strict_hard` (today's default) | 0.618830 | 39 | 0 | 0 | 0 |
| `broad_hard` | **0.819921** | 42 | 0 | 0 | 0 |
| `strict_soft` | **0.819921** | 42 | 0 | 0 | **40** |
| `broad_soft` | 0.819921 | 42 | 0 | 0 | 0 |

Rerun 2026-09-17. The first run of the two soft arms did not implement soft enforcement and
reported `strict_soft` at 0.657195 with 3 advisories; see section B.

### A — hard versus hard, different thresholds

Broad ahead on **40 of 42**, strict ahead on **0**, two tied. Mean **+0.201091**, median gain
+0.165763, largest single gain +0.670973, **largest loss +0.000000**. No module regresses, and
both arms deliver an inventory their own policy accepts.

The two tied modules are `pPR-1_2E_LD` and `pPR-1_2E_LN` — already at their ceiling under both.

### B - same thresholds, hard versus soft

**Retracted once, now measured.** The first version of this experiment did not implement soft
enforcement at all: `solver_bounds()` returned identical arguments for `hard` and `target`, so
both "soft" arms were a hard solver read by a permissive validator. That is recorded here
rather than deleted, because the wrong version of this table was quoted in three documents and
a reader needs to be able to tell which one they have.

Rerun 2026-09-17 against the repaired solver, where a rule named in `soft_rules` becomes a DNA
Chisel *objective* rather than a constraint:

| arm | mean CAI | improved | advisory warnings | seconds |
|---|---:|---:|---:|---:|
| `strict_hard` (today's default) | 0.618830 | 39 / 42 | 0 | 2.4 |
| `broad_hard` | 0.819921 | 42 / 42 | 0 | 2.7 |
| **`strict_soft`** | **0.819921** | **42 / 42** | **40** | 6.2 |
| `broad_soft` | 0.819921 | 42 / 42 | 0 | 2.2 |

**Soft enforcement at the strict thresholds recovers the entire gain of widening the band.**
Not approximately: `strict_soft` and `broad_hard` agree on every one of the 42 modules to six
decimal places, and their inventory means are bit-identical. The strict band costs 0.201091
mean CAI when it is a wall and **nothing at all** when it is a gradient.

The earlier, broken version of this arm scored 0.657195 and recovered only three modules for
+0.038. That number was a permissive validator reading a narrow solver, which is why it
coincided exactly with the one-line `GC_BAND` patch. The real figure is +0.201091 over 40
modules.

`broad_soft` ties `broad_hard` on all 42 with zero warnings, which is the control: at
thresholds that do not bind there is nothing to soften, so the mechanism costs nothing when it
is not needed.

**What this changes, and what it does not.** On these two axes `strict_soft` dominates both
alternatives: it matches the broad band's CAI while keeping the 40 advisories the broad band
discards. Each advisory names its window, for example `pPR-1_14A_LD5N` reporting
*GC 0.660 outside (0.35, 0.65) in the window at 30 (worst 0.720 at 60)*.

It is not free, and the cost is not only the 2.6x runtime. A soft policy **delivers sequences
outside Twist's published local-GC guidance** and records that it did; a hard policy refuses to
produce them. That is a trade of a guarantee for a documented warning, and which one a team
wants depends on whether anyone reads the warnings. This experiment measures the trade; it does
not make the choice. **No default is changed here.**
### The 0.660 windows, and what each policy does with them

Forty of the 42 modules have at least one 50 nt window outside the strict band once the solver
is free to improve them; the recurring figure is **0.660**, the same number that has run through
this investigation from the first substrate-scope defect onward. What separates the arms is what
they do about it:

| policy | those 40 modules | what you are told |
|---|---|---|
| `strict_hard` | 39 improved as far as the band allows, 3 left at the deposited sequence | nothing; the band is invisible in the output |
| `broad_hard` | all 42 improved freely | nothing; 0.660 is inside the declared band |
| `strict_soft` | all 42 improved freely | **40 advisories**, each naming its window |

The three modules `strict_hard` cannot improve at all — `pPR-1_14A_LD5T`, `pPR-1_14E_LN5N`,
`pPR-1_19E_LN5N` — are the extreme case of the same effect: every candidate that improves them
lands over the band, so a hard band leaves them at 0.188–0.199 while any other arm reaches
0.840–0.860. They are the largest single gains in the table (+0.671, +0.643, +0.643).

### Why the retracted number coincided with the one-line patch

The retracted `strict_soft` scored **0.657195**, exactly what the orchestrator's audit obtained
by patching `recoding.GC_BAND` alone. At the time this was reported as independent
corroboration. It was the opposite: both configurations were the *same* thing — a narrow solver
read by a permissive validator — so the agreement was a duplicate, not a replication. It is kept
here because an exact coincidence between two routes is the sort of evidence that feels
strongest and, in this case, was diagnostic of the bug rather than of the result.

## What this does not settle

- **Which profile should be the default.** The benchmark above is the input, and it now points
  somewhere specific: `strict_soft` matches the broad band's CAI exactly while keeping all 40
  advisories, so the choice is no longer "conservative or optimised" but "refuse or record". It
  still does **not** decide the question, because none of it is evidence of vendor acceptance —
  the only thing that matters for an order — and a soft policy ships the out-of-band sequences a
  hard one refuses. The decision needs a human willing to own that risk.
- **Whether a wider band manufactures.** Twist and IDT both score submissions with proprietary
  models. A published guideline is not an acceptance guarantee and a looser one is not a
  licence.
- ~~Hard versus soft at matched thresholds.~~ **Measured 2026-09-17**, after the first attempt
  was retracted for not implementing soft enforcement. Soft at the strict thresholds reaches the
  broad band's CAI exactly, on every module, while reporting 40 advisories. See the benchmark
  above; the question is closed and the answer is the opposite of the retracted one.
- **That CAI predicts expression.** It does not. `docs/objectives.md` §8 holds.
- **IDT's rules.** We order from IDT, whose oPools page states no numeric GC threshold at all.
  The gBlocks FAQ is a different product. No profile here is an IDT acceptance model.

---

## Using one

    from clippr.synthesis_profile import TWIST_CODON_GUIDANCE
    from clippr.recoding import recode_inventory

    report = recode_inventory(kit, table, label="chlamy", profile=TWIST_CODON_GUIDANCE)
    report.as_dict()["profile"]["version"]      # bound to the result

By name, for a config file or a form: `"clippr-strict-legacy"`,
`"twist-codon-optimisation-guidance"`, `"broad-experimental"`. An unknown name raises rather
than falling back to a default — a typo must not silently select a different policy.

Feasibility and quality are answered separately:

    substrate_problems(insert, block, profile)   # hard breaches — refuses an order
    substrate_warnings(insert, block, profile)   # target breaches — reported, still shipped
    substrate_findings(insert, block, profile)   # both, split

---

## Related

- `docs/gc_band_and_constraint_model.md` — the measurement, the vendor sources and their limits
- `docs/objectives.md` §0 — evaluate the contract on what leaves the building
- `docs/reference_implementation_findings.md` §3.2 — where the three-way split came from
