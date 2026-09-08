# CLIPPR

**Design tooling for synthetic PPR regulators.** iGEM Marburg 2026.

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

Python 3.10+. Runs anywhere, including a fresh Colab runtime — see
[`notebooks/CLIPPR_designer.ipynb`](notebooks/CLIPPR_designer.ipynb) for a form-driven
version that needs no code.

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

**3. Its synthesis QC actually discriminates.** The baseline flagged WARNING on 200 of 200
designs — a verdict that never varies carries no information. The same sequences here split
12 PASS / 83 WARNING / 105 FAIL, with every threshold documented by the percentile it sits
at.

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

Plus 283 unit tests, including an exhaustive comparison of the overhang feasibility filter
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
- **Nothing here has been validated at the bench.**

## Licence

MIT — see [LICENSE](LICENSE). [NOTICE.md](NOTICE.md) records the third-party data that
ships with the package, its attribution, and the published sources every scientific
constant derives from.
