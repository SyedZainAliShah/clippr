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
| `strict_soft` | 0.657195 | 42 | 0 | 0 | 3 |
| `broad_soft` | 0.819921 | 42 | 0 | 0 | 0 |

### A — hard versus hard, different thresholds

Broad ahead on **40 of 42**, strict ahead on **0**, two tied. Mean **+0.201091**, median gain
+0.165763, largest single gain +0.670973, **largest loss +0.000000**. No module regresses, and
both arms deliver an inventory their own policy accepts.

The two tied modules are `pPR-1_2E_LD` and `pPR-1_2E_LN` — already at their ceiling under both.

### B — same thresholds, hard versus soft

**This is the experiment nothing had run, and it changes the picture.**

| thresholds | soft ahead | hard ahead | tied | mean difference |
|---|---:|---:|---:|---:|
| strict (0.35–0.65) | 3 | 0 | 39 | +0.038365 |
| broad (0.15–0.85) | 0 | 0 | **42** | **+0.000000** |

At broad thresholds, hard and soft enforcement are **exactly identical** — every module, to
twelve decimal places. Soft enforcement has no independent value. It matters only when a
threshold *binds*, and what it does then is recover the candidates that threshold vetoes.

So the soft-penalty architecture is not the advantage it looked like. It is a **compensation
for a threshold that is too tight**, and the honest comparison is between thresholds, not
between enforcement models. An earlier draft of the reference findings implied otherwise.

### The three unchanged modules, fully explained

`strict_soft`'s three advisory warnings are **exactly** the three modules `strict_hard` leaves
unchanged:

    pPR-1_14A_LD5T   GC 0.660 outside (0.35, 0.65) at window 14
    pPR-1_14E_LN5N   GC 0.660 outside (0.35, 0.65) at window 55
    pPR-1_19E_LN5N   GC 0.660 outside (0.35, 0.65) at window 55

All three at **0.660** — the same figure that has run through this entire investigation, from
the first substrate-scope defect onward. Those modules are not unimprovable; every candidate
that improves them lands 0.010 over the band.

Note the scale difference: softening recovers three modules for +0.038, while widening lifts
**forty** for +0.201. The band does not merely veto three modules outright — it constrains the
search for nearly all of them.

### An independent arrival at the same number

`strict_soft` scores **0.657195**, which is exactly what the orchestrator's audit obtained by
patching `recoding.GC_BAND` alone. That is not a coincidence: both configurations are a narrow
solver with a permissive validator. Two different routes to the same value, which corroborates
both — and identifies precisely what the "one-line" patch was doing.

## What this does not settle

- **Which profile should be the default.** The benchmark above is now that input, and it is
  favourable: +0.201 mean, 40 of 42 improved, nothing worse, and every arm delivering an
  inventory its own policy accepts. It still does **not** decide the question, because none of
  it is evidence of vendor acceptance — which is the only thing that matters for an order.
  The decision needs a human who is willing to own that risk.
- **Whether a wider band manufactures.** Twist and IDT both score submissions with proprietary
  models. A published guideline is not an acceptance guarantee and a looser one is not a
  licence.
- ~~Hard versus soft at matched thresholds.~~ **Now tested** — see the benchmark above. At
  non-binding thresholds the two are identical to twelve decimal places; soft only recovers
  what a binding threshold vetoes. That question is closed.
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
