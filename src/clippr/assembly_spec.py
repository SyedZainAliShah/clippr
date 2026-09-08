"""Layer 2 of 3 — assembly mechanics. How Golden Gate works for this construct.

Facts about the deposited vectors and the cut geometry: true of the assembly, but not
biology in the sense `biology.py` means, and not a choice in the sense `policy.py` means.

    biology.py         published fact
    assembly_spec.py   assembly mechanics      <- you are here
    policy.py          this project's choices
"""
from __future__ import annotations

#: Deposited destination overhangs as *coding sites*, 5' then 3'.
#:
#: The 3' entry is the coding site, not the overhang: the enzyme leaves its reverse
#: complement. Level 0 is stored ("CTCA", "CGAG"), but the strands that actually meet in
#: the tube are CTCA and CTCG -- which differ by a single nucleotide and score 0.794
#: predicted fidelity, where the stored pair scores a misleading 1.000.
#: Build reaction sets with `overhangs.reaction_overhangs()`, which applies the flip,
#: rather than passing these values to a fidelity function directly.
DESTINATION_OVERHANGS: dict[str, tuple[str, str]] = {
    "level_minus1": ("ACAT", "TTGT"),
    "level0": ("CTCA", "CGAG"),
    "level1": ("GGAG", "CGCT"),
}

#: Residue offsets from an ARELF motif at which the corpus designs place cuts.
#: Measured across all 200 stored designs by `validation/compare_cuts.py`. Recorded as an
#: observation about the corpus, not, not as a constraint CLIPPR imposes -- cuts
#: are *not* restricted to these positions; see `arelf.tunable_cuts`.
ARELF_OFFSETS_OBSERVED: tuple[int, ...] = (3, 5, 7, 8, 13, 16, 17, 18, 21, 22, 28, 29)

#: Kept for backwards compatibility; prefer ARELF_OFFSETS_OBSERVED.
ARELF_OFFSETS: tuple[int, ...] = tuple(range(12))

#: A codon-aligned cut at residue `c` splits the CDS at nucleotide `3c`, and the overhang
#: it leaves is `cds[3c - overhang_len : 3c]` -- the wobble base of codon `c-2` plus all of
#: codon `c-1`. Established by measurement: 750/750 cut contexts across 200 designs.
#: The overhang length itself is enzyme-specific and read from Biopython, not fixed here.
CUT_IS_CODON_ALIGNED: bool = True
