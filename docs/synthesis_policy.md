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

## What this does not settle

- **Which profile should be the default.** The measured case supports making the policy
  explicit, not widening it. Changing the default changes published results and needs its own
  benchmark: all 42 modules, fixed seeds and budgets, the same full-substrate checks, and
  per-module differences — not a mean alone.
- **Whether a wider band manufactures.** Twist and IDT both score submissions with proprietary
  models. A published guideline is not an acceptance guarantee and a looser one is not a
  licence.
- **Hard versus soft at matched thresholds.** Every experiment so far varies hard bounds
  against hard bounds. Whether soft penalties beat hard constraints *at the same numbers* is a
  separate question and has not been tested.
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
