# CLIPPR

**Design tooling for synthetic PPR regulators.** iGEM Marburg 2026.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SyedZainAliShah/clippr/blob/main/notebooks/CLIPPR_designer.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-1a7f5a.svg)](LICENSE)


Give it a target RNA sequence and it returns a PPR protein that binds that sequence, a
synthesisable coding sequence, a Golden Gate assembly plan, and the DNA fragments to order.

PPR proteins are built from tandem ~31-residue repeats, and each repeat reads exactly one
RNA base through two specificity residues — the **PPR code** (`TN`→A, `NN`→C, `TD`→G,
`ND`→U). The protein is therefore a deterministic function of the target, which is why this
package generates designs rather than searching a parts catalogue.

Built for *Chlamydomonas reinhardtii*, which is the organism this project's wet lab works
in. Any host can be used by supplying a codon table.

---

## Install

```bash
pip install git+https://github.com/SyedZainAliShah/clippr.git
```

Python 3.10+. Or skip the install entirely and **[open the notebook in Colab](https://colab.research.google.com/github/SyedZainAliShah/clippr/blob/main/notebooks/CLIPPR_designer.ipynb)** — a form-driven version that needs no code.

## A worked example

```python
from clippr import design_oneshot

r = design_oneshot("AAAAUGUGG", outdir="out")
print(r["summary"])
# AAAAUGUGG (9S) -> 302 aa, 906 nt, 4 fragments. Fidelity 0.830. QC PASS. 109.00 EUR.
```

The target's length picks the architecture: 9, 14 or 19 bases give 9S, 14S or 19S. Four
files land in `out/` — an order CSV, the oligos as FASTA, the assembled gene, and an
annotated GenBank record with every repeat labelled by the base it reads.

Every design also carries its own decision record:

```python
print(r["audit"].report())
```

```
design audit — AAAAUGUGG (9S)
  host           c_reinhardtii_nuclear, genetic code 1
  enzyme profile igem_rfc1000 (BsaI, BbsI, SapI)
  cuts           76, 151, 226
  fidelity       0.830 predicted  (ceiling 0.830, set by the destination pair)
  constraints    satisfied
  QC             PASS

  selected overhangs
    CGAC  selected   cut   76  highest predicted fidelity among feasible  [2 local realizations]
    AATG  selected   cut  151  highest predicted fidelity among feasible  [3 local realizations]
    AGCA  selected   cut  226  highest predicted fidelity among feasible  [6 local realizations]
```

Rejected candidates appear there too, each with the reason it was ruled out. A surprising
design can be interrogated rather than taken on trust.

---

## What this adds over the published workflow

**1. It optimises for the right host.** The published GRASP workflow targets
*E. coli*. Its coding sequences score CAI 0.411 against the *Chlamydomonas* nuclear table.
CLIPPR uses the host's own codon usage. (A codon-usage statement, not an expression claim —
no expression was measured.)

**2. It removes the repeats that tandem-repeat proteins produce.** "Use the best codon
everywhere" emits the same DNA once per repeat. Measured across 9S–19S targets:

| setting | duplicated 20-mers | longest exact repeat | CAI |
|---|---|---|---|
| best codon only | 425 | **95 nt** | 0.843 |
| unique 20-mers *(default)* | 0 | 19 nt | 0.793 |
| E. coli-optimised baseline | 129 | 44 nt | 0.411 |

Costs about 6% codon adaptation. Across 200 designs, 197 come out with no repeated 20-mer.
Uniquifying is an *objective*, not a constraint, so it cannot be promised.

**3. It checks the host, in the tier that matters.** A PPR binds **RNA**, so the question
is not whether a sequence appears in the genome but whether it appears in a transcript, in
the sense orientation. Measured against the *Chlamydomonas* chloroplast (203,828 bp,
34.5% GC; 109 annotated transcripts covering 43.5% of it), 200 random targets per length:

| target length | in genomic DNA | **in a transcript** |
|---|---|---|
| 9 nt | 97 of 200 (48%) | **32 of 200 (16%)** |
| 14 nt | 0 | 0 |
| 19 nt | 0 | 0 |

Counting genomic DNA on both strands overstates the risk roughly threefold: a match in a
non-transcribed region is not an RNA off-target, and neither is a reverse-complement match
in DNA. Every design is checked and reports both tiers separately, naming the gene when a
transcript is hit.

Occurrence is a *necessary* condition for an off-target interaction, never a sufficient
one. No binding affinity is predicted.

**4. Its synthesis QC actually discriminates.** The baseline flagged WARNING on 200 of 200
designs — a verdict that never varies carries no information. The same sequences here split
12 PASS / 83 WARNING / 105 FAIL, with every threshold documented by the percentile it sits
at.

## Designing a whole library

Fifty regulators is not fifty designs in fifty folders — it is one order sheet.

```python
from clippr import design_library

lib = design_library(my_targets, outdir="library")
print(lib.summary())
```

```
library of 50 designs
  fragments      212 across all designs, 54,900 bases
  pooled cost    109.00 EUR as one pool, versus 5450.00 separately — 5341.00 EUR saved
  QC             48 PASS, 2 WARNING
  closest targets UUACACGUG/ACGUACGUA differ at 5
```

An oligo pool is priced per pool, essentially flat across these sizes, so ordering a
library as one pool rather than fifty is the single largest cost decision in the workflow.
`write_library` emits a combined order sheet, a vendor-ready oPool CSV, a QC table, one
GenBank per design, and the cross-talk matrix.

### Cross-talk: two tiers, one gate

Cross-talk is a property of the set, not of any member: a design that is perfect alone is
useless if another target in the library sits one base away. `crosstalk.py` reports it in
two tiers that are deliberately not mixed.

| tier | what it is | status |
|---|---|---|
| **A — Hamming distance** | two targets differ in *k* of *n* positions | **the hard criterion, and the only thing that gates** |
| **B — predicted affinity** | score the PPR designed for A against B, relative to its own target | an annotation; gates nothing |

Tier A makes no biological claim. It says two sequences differ in *k* places, which is
geometry, and stays true whatever anyone later learns about PPR binding.

Tier B needs a PPR specificity table, and comes with a caveat that cannot be argued away:
the available table (Yan et al., as distributed with PPRmatcher) was derived from **P-type**
PPR motifs, while the GRASP scaffold here is **S-type**. All four code pairs GRASP uses do
score their cognate base highest in that table, which is reassuring, but a P-type model has
not been shown to apply to an S-type scaffold. So a high score means *look*, never *fail*.

**The table is not distributed with CLIPPR.** PPRmatcher carries no licence, so vendoring
it would be a rights problem however useful it is. Supply your own path; with no table the
module reports tier A alone and says so.

```python
from clippr import compare_tiers, load_ppr_scores

scores = load_ppr_scores("Yan.tsv")          # obtained separately
print(lib.crosstalk(scores=scores))
compare_tiers(lib.targets, scores)["tiers_disagree"]
```

The useful question about a model whose applicability is unproven is not "is it right" but
**"does it point anywhere Hamming does not"** — which is what `compare_tiers` answers. Its
`tiers_disagree` flag never means the gate moved; it means a pair is worth a human look.

### DNA shared between members — the other library risk

Cross-talk asks whether two PPRs could bind each other's target. **Homology asks whether two
*genes* share enough identical DNA to recombine**, which is a different question with a
different answer. Every member of a PPR library carries the same scaffold, so it is never
trivially no. Measured on five 9S designs:

| | before | after |
|---|---|---|
| longest shared stretch | 84 nt | **47 nt** |
| pairs sharing ≥ 50 nt | 10 of 10 | **0 of 10** |

Every shared stretch began at position 0 in both members and decoded to `MQGGNSEEPRKSFDERPER…`
— the fixed 23-residue N-terminal scaffold. It is not a codon-diversification failure: across
369 nt of *identical protein* the longest shared DNA run was only 59 nt, so the repeat body is
already well separated. The scaffold is simply protein-identical by construction and receives
the same codons every time.

The fix is constructive. Each member gets its own synonymous encoding of the scaffold, locked in
place — deterministic in the member index, so a library stays reproducible.

```python
from clippr import diversify_library

res = diversify_library(my_targets)
print(res["max_before"], "->", res["max_after"])   # 84 -> 47
```

Two simpler approaches were tried first and are documented in `homology.py` so they are not
retried: forbidding the shared stretch outright only forces a one-base change, and forbidding
every k-mer inside it is **not monotone** — at one ban width the worst stretch went from 84 nt
to 86. A diversified member is now accepted only when it is no worse than the baseline it
replaces.

**What this does not claim.** There is no evidence here for a 50 nt danger threshold in
*Chlamydomonas* — that default is a rule of thumb from general practice and nothing in this
repository derives one. A shared stretch is a *necessary substrate* for homologous recombination,
never a prediction that it will occur, and the risk depends on the physical library architecture:
one construct per strain is a different situation from many constructs entering the same nuclear
genome, or from a pooled DNA mixture. The defensible sentence is that **independently designed
members acquired substantial unintended DNA identity through a shared scaffold, and synonymous
redesign reduced the longest shared tract from 84 to 47 nt.** Everything beyond that needs the
bench.

## Is the chosen design good relative to the alternatives?

`design_oneshot` ranks assembly plans by predicted fidelity and keeps the first that works,
discarding the rest unexamined — so it could not answer that question. `search.py` carries
several plans all the way through codon optimisation and QC, then applies a **declared priority
hierarchy** rather than arbitrary weights:

1. every constraint satisfied and QC not FAIL
2. predicted fidelity within a tolerance of the best feasible value
3. prefer QC PASS, then the higher optimiser score
4. tie-break on fidelity

```python
from clippr import design_searched

r = design_searched("AAAAUGUGG", budget=8)
print(r["certificate"])
```

It returns **one answer with the alternatives shown** — what was considered, how far the choice
sits from the best available value on each objective, whether anything dominates it, and what
the nearest alternatives would cost. Not a Pareto front: arbitrating a trade-off surface is not
work a wet lab asked for.

**Measured caveat, stated because it matters.** Across all three architectures, predicted
fidelity was **identical for every candidate** and QC passed for every candidate — fidelity is
capped by the destination overhang pair, and `safe_overhangs` has already removed the designs
that would have failed QC. So the hierarchy is currently *single-objective in practice*, and
this is not a multi-objective optimiser. What it does deliver is a better design than the
first-feasible plan (it changed the answer for all three architectures) and an evidenced answer
to the question above.

## Constraints are declared, not assumed

Which Type IIS sites are excluded is a *policy*, and the reasons differ in kind:

| profile | excludes | why |
|---|---|---|
| `assembly` | BsaI, BbsI | this assembly's own chemistry |
| `igem_rfc1000` *(default)* | + SapI | iGEM RFC[1000] requires BsaI and SapI absent |
| `moclo_compat` | + BsmBI | keeps later MoClo-family levels open — a preference |

```python
design_oneshot("AAAAUGUGG", enzyme_profile="moclo_compat")
```

The package keeps published fact, assembly mechanics and project choices in separate
modules (`biology.py`, `assembly_spec.py`, `policy.py`) so a reader can tell what GRASP
requires from what CLIPPR chose.

## Verification

Checked against a fixed 200-design reference corpus:

| check | result |
|---|---|
| Protein sequence | **200/200 exact** |
| Cut geometry | **750/750 exact** |
| Oligos, given the same CDS and cuts | **200/200 byte-identical** |
| Codon constraints, independently verified | **200/200** |
| End-to-end pipeline | **200/200, zero exceptions, seed-reproducible** |

Plus 429 unit tests, including an exhaustive comparison of the overhang feasibility filter
against an independently written brute-force oracle.

```bash
pytest                                    # unit tests
python validation/report.py --quick       # every validation check, one table
```

Byte-identical oligos are the load-bearing result: for equivalent inputs this package
reproduces the corpus construction exactly, so the differences above are attributable to
architecture and added constraints rather than a silently changed construct definition.

## Provenance

Every constant traces to a primary published source, not to any other implementation:

- scaffold sequences and the PPR code are derived from Farley et al. (2025) and its
  supplementary tables by `tools/derive_scaffold.py`, which assembles the deposited
  modules per the published recipe;
- ligation matrices come from the Pryor et al. (2020) supplement directly, and ship with
  the package (see [NOTICE.md](NOTICE.md) for their attribution);
- `biology.py`, `assembly_spec.py` and `policy.py` separate published fact from assembly
  mechanics from this project's own choices, so a reader can tell which is which.

## Caveats

- **Predicted fidelity** comes from published ligation-count matrices, not measured
  assembly efficiency.
- **QC** is a sequence-complexity and feasibility check, not calibrated against vendor
  outcomes.
- **Cost** is an IDT oPools list price, not a quote, and is flat across pool sizes in this
  range — it cannot rank designs.
- **`orthogonal.py` is a capability, not a validated result.** Every other module is checked
  against a 200-design oracle; this one has unit tests only, because no ground truth exists.
- **Predicted PPR affinity is an unvalidated annotation**, from a P-type table applied to an
  S-type scaffold. It is reported beside sequence separation and never allowed to override
  it — see [Cross-talk: two tiers, one gate](#cross-talk-two-tiers-one-gate).
- **Nothing here has been validated at the bench.**

### Computable properties versus model-based annotations

CLIPPR keeps these apart deliberately. Sequence uniqueness, restriction-site absence,
translation, codon constraints and assembly geometry are **directly computable** — the
package either satisfies them or reports that it could not. Predicted ligation fidelity,
synthesis QC verdicts and host-occurrence warnings are **annotations**, useful for ranking
and screening but not guarantees of binding, expression or assembly success. Experimental
validation is required for any biological claim.

## Licence

MIT — see [LICENSE](LICENSE). [NOTICE.md](NOTICE.md) records the third-party data that
ships with the package, its attribution, and the published sources every scientific
constant derives from.
