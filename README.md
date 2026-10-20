# PPR-opt

**Design tooling for synthetic PPR regulators.** iGEM Marburg 2026.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SyedZainAliShah/clippr/blob/main/notebooks/PPROPT_designer.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-1a7f5a.svg)](LICENSE)


Give it a target RNA sequence and it returns a PPR protein **designed to recognise that
sequence under the PPR code**, a synthesisable coding sequence, a Golden Gate assembly plan,
and the DNA fragments to order. No binding was measured; the recognition claim is what the
code predicts, not an experimental result.

PPR proteins are built from tandem ~31-residue repeats, and each repeat reads exactly one
RNA base through two specificity residues — the **PPR code** (`TN`→A, `NN`→C, `TD`→G,
`ND`→U). The protein is therefore a deterministic function of the target, which is why this
package generates designs rather than searching a parts catalogue.

Built for *Chlamydomonas reinhardtii*, which is the organism this project's wet lab works
in. Any host can be used by supplying a codon table.

---

## Start here

Two notebooks, both form-driven and needing no code:

| | what it does | you need |
|---|---|---|
| [**Designer**](https://colab.research.google.com/github/SyedZainAliShah/clippr/blob/main/notebooks/PPROPT_designer.ipynb) | designs a new PPR for any target you choose | nothing |
| [**Inventory**](https://colab.research.google.com/github/SyedZainAliShah/clippr/blob/main/notebooks/PPROPT_inventory.ipynb) | builds from the deposited 42-module GRASP kit | Supplementary Table S1 |

Each notebook ends with a plain-language glossary of every term it uses.

## Install

**As a package**, Python 3.10+:

```bash
pip install git+https://github.com/SyedZainAliShah/clippr.git
```

**To reproduce the figures below exactly**, use the pinned set rather than the floors in
`pyproject.toml`:

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # Linux/macOS: .venv/bin/python
.venv/Scripts/python -m pip install -e .
.venv/Scripts/python -m pytest                            # expect 1226 passed, 0 failed
```

## A worked example

```python
import json
from ppropt import design_oneshot

codon_table = json.load(open("kazusa_3055.json", encoding="utf-8"))
r = design_oneshot("AAAAUGUGG", codon_table=codon_table, genetic_code=1,
                   check_offtarget=False, outdir="out")
print(r["summary"])
# AAAAUGUGG (9S) -> 302 aa, 906 nt, 4 fragments. Fidelity 0.828. Coding-span QC PASS. Orderable under the declared software checks.
```

**Two verdicts, different scopes, both always reported.** `Coding-span QC` checks the gene
itself. The synthesis verdict checks the sequences you would actually send to a vendor —
the entry overhang, the coding slice and the enzyme wrapper together. One passing says
nothing about the other.

A design that fails the second is never written as an order: its files go out under
`diagnostic_` keys with `REJECTED` in every filename, so they can be inspected but no
downstream step can mistake them for an order.

> **Changed 2026-10-12.** Until then the de novo route produced nothing orderable. The GC
> band was enforced on the coding span, while the ordered fragment also carries the wrapper
> — so a window straddling that boundary was never shown to the optimiser, and all 200
> corpus designs broke it. The optimiser now sees the whole ordered fragment and re-spells
> codons against it. The band itself did not change. `wrapper_aware=False` reproduces the
> earlier behaviour exactly, and is the negative control.

## Two routes

| | synthesise | use parts you own |
|---|---|---|
| designs | any target | only what the kit can spell |
| you supply | a target sequence | the GRASP kit, and Table S1 |
| output | fragments to order | which modules to use, and which wells to pick |
| new DNA | the whole gene | recoded modules only |

## What this does not claim

The limits matter more than the features here, so they are not buried:

- **Nothing has been validated at the bench.** Every number is a software check or a
  prediction from a published model.
- **No vendor has accepted or refused a sequence.** Nothing was ever submitted. The policy
  that refuses fragments is this project's own, and is stricter than any vendor rule we
  could verify first-hand — see [`docs/vendor_rule_comparison.md`](docs/vendor_rule_comparison.md).
- **Predicted fidelity** comes from published ligation-count matrices, not measured
  assembly efficiency.
- **QC** is a sequence-complexity check, not calibrated against vendor outcomes.
- **Predicted PPR affinity is an unvalidated annotation**, from a P-type table applied to an
  S-type scaffold. It is reported beside sequence separation and never allowed to override it.
- **`orthogonal.py` is a capability, not a validated result** — unit tests only, because no
  ground truth exists.

**Computable properties and model-based annotations are kept apart deliberately.** Sequence
uniqueness, restriction-site absence, translation, codon constraints and assembly geometry
are directly computable: the package either satisfies them or reports that it could not.
Predicted fidelity, QC verdicts and host-occurrence warnings are annotations — useful for
ranking and screening, not guarantees of binding, expression or assembly.

## Reading further

This README is the overview. The detail lives in `docs/`:

| | |
|---|---|
| [`objectives.md`](docs/objectives.md) | what is optimised, what is constrained, and what none of it measures |
| [`gc_band_and_constraint_model.md`](docs/gc_band_and_constraint_model.md) | where the GC band came from, and what it does not establish |
| [`vendor_rule_comparison.md`](docs/vendor_rule_comparison.md) | what each vendor actually publishes, verified first-hand |
| [`synthesis_policy.md`](docs/synthesis_policy.md) | the named policies and their sources |
| [`methodology_sources.md`](docs/methodology_sources.md) | the published source behind every scientific constant |
| [`workflow_contract.md`](docs/workflow_contract.md) | what each step promises its caller |
| [`verification_commands.md`](docs/verification_commands.md) | one invocation per claim in this README |

## Verification

Thirteen release checks, each either run or explicitly marked as read:

```bash
.venv/Scripts/python validation/release_check.py
```

They recompute rather than re-read: every recorded result states the artefacts it was
produced from, and the checker re-derives those identities instead of trusting the record.
A record naming no source, or whose stated hash no longer matches, fails.

Independent validators (V5 digest-and-ligate, V6 constraints, V7 export semantics) are
written without importing this package, so they cannot inherit its mistakes.

## Licence

MIT — see [LICENSE](LICENSE). [NOTICE.md](NOTICE.md) records the third-party data that
ships with the package, its attribution, and the published sources every scientific
constant derives from.
