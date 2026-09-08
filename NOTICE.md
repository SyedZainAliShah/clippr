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

## Sources of derived constants

**Scaffold sequences, the PPR code and the repeat template** in `src/clippr/biology.py`
and `src/clippr/scaffold.json` are derived from:

> Farley KV, et al. (2025). Nucleic Acids Research.
> https://doi.org/10.1093/nar/gkaf1169

`tools/derive_scaffold.py` assembles them from the deposited modules in the paper's
supplementary tables, following the published recipe. Re-run that script rather than
editing the values by hand.

**The iGEM RFC[1000] requirement** encoded in `policy.ENZYME_PROFILES` — that BsaI and
SapI recognition sites be absent from participating parts — comes from iGEM's Type IIS
assembly documentation.

## What is *not* derived from a published source

`src/clippr/policy.py` holds this project's own choices: which enzyme sites to exclude,
which host to optimise for, where QC thresholds sit. Each carries its reason in the file.
The three-layer split across `biology.py`, `assembly_spec.py` and `policy.py` exists so
that published fact stays distinguishable from project policy.
