# Third-party material and attribution

CLIPPR itself is MIT-licensed; see [LICENSE](LICENSE). This file records material that
travels with the package but originates elsewhere, and the sources every scientific
constant is derived from.

## Data redistributed with the package

**`src/clippr/data/matrices/BsaI-HFv2.xlsx`** and **`src/clippr/data/matrices/BbsI-HF.xlsx`**

Supplementary Tables S1 and S4 of:

> Pryor JM, Potapov V, Kucera RB, Bilotti K, Cantor EJ, Lohman GJS (2020).
> *Enabling one-pot Golden Gate assemblies of unprecedented complexity using
> data-optimized assembly design.* PLoS ONE 15(9): e0238592.
> https://doi.org/10.1371/journal.pone.0238592

PLoS ONE publishes under **CC BY 4.0**, so these are redistributed here with attribution.
They are required at runtime — `overhangs.set_fidelity` cannot score an overhang set
without them, which is why they ship inside the package rather than being downloaded.

## Data deliberately *not* redistributed

**PPR specificity scores (`Yan.tsv`, as distributed with
[PPRmatcher](https://github.com/ian-small/PPRmatcher))** are supported by
`src/clippr/crosstalk.py` but are **not** included in this package.

That repository declares no licence, so no redistribution permission is granted and none
can be inferred. `load_ppr_scores` therefore takes a path the user supplies; with no table,
`crosstalk` reports sequence separation alone and names the source. The scores read as a
plain TSV, so any table in that layout works — nothing in the code is specific to that file
beyond its format.

Note also that the table derives from **P-type** PPR experiments while the GRASP scaffold
is **S-type**, which is why the module treats its output as an annotation and never as a
criterion.

## Sources of derived constants

**Scaffold sequences, the PPR code and the repeat template** in `src/clippr/biology.py`
and `src/clippr/scaffold.json` are derived from:

> Dennis M, Low SY, Viljoen A, Pullakhandam A, Colas des Francs-Small C, Campbell-Clause L, Bond CS, Small I, Kwok van der Giezen FM (2025). GRASP: a modular toolkit for building synthetic pentatricopeptide repeat RNA-binding proteins. Nucleic Acids Research.
> https://doi.org/10.1093/nar/gkaf1169

`tools/derive_scaffold.py` assembles them from the deposited modules in the paper's
supplementary tables, following the published recipe. Re-run that script rather than
editing the values by hand.

**The GRASP module inventory** in `src/clippr/parts.json` is derived from Supplementary Table S1
of the same paper by `tools/derive_parts.py`: plasmid identity, Golden Gate overhangs, plate
positions and insert *lengths*. **Insert sequences are deliberately not included** — selecting
modules does not need them, since the point of the route is that the user already holds the
physical plasmids, and omitting them keeps this a derived index rather than a republication of
the table. Anyone needing the sequences has the paper and the Addgene kit entry.

The *approach* of compiling a target RNA into an ordered part list was taken from the reference
implementation, which we read. Our implementation shares no source with it and takes its data
from the published table, but the idea is credited here rather than presented as our own.

**The iGEM RFC[1000] requirement** encoded in `policy.ENZYME_PROFILES` — that BsaI and
SapI recognition sites be absent from participating parts — comes from iGEM's Type IIS
assembly documentation.

## What is *not* derived from a published source

`src/clippr/policy.py` holds this project's own choices: which enzyme sites to exclude,
which host to optimise for, where QC thresholds sit. Each carries its reason in the file.
The three-layer split across `biology.py`, `assembly_spec.py` and `policy.py` exists so
that published fact stays distinguishable from project policy.
