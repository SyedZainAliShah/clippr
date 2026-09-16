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
#: Two different things are called "level 1" in this package and they are not the same set.
#: These are the **acceptor vector's own** ends, which the one-shot design path assembles a
#: fresh construct into. A GRASP *block* assembly's level-1 reaction is a different reaction
#: with different ends -- the compiled product's own interfaces -- see `LEVEL1_BLOCK_ENDS`.
#: Neither figure contradicts the other; they describe different tubes.
DESTINATION_OVERHANGS: dict[str, tuple[str, str]] = {
    "level_minus1": ("ACAT", "TTGT"),
    "level0": ("CTCA", "CGAG"),
    "level1": ("GGAG", "CGCT"),
}

#: Which enzyme cuts at each assembly stage, and which measured matrix scores it.
#:
#: **This is not a free choice.** Every one of the 42 deposited GRASP level -1 plasmids
#: carries exactly two BbsI sites and no BsaI site, so a level-0 reaction -- the one that
#: assembles module inserts into a block -- is cut by BbsI and must be scored against the
#: BbsI-HF matrix. Scoring it with BsaI-HFv2 describes a reaction that does not happen.
#:
#: The level-1 reaction, which joins assembled blocks, is a BsaI reaction. It was reported
#: unscorable here until 2026-09-16 on the grounds that its ends "are set by the block plasmid
#: and its destination, which this package does not have". That was wrong twice over: the
#: block plasmids were in the provenance directory all along, and they turn out not to be
#: needed, because the released block fragment is exactly the joined module inserts. See
#: `LEVEL1_BLOCK_ENDS`.
ASSEMBLY_STAGES: dict[str, dict] = {
    "level0": {
        "enzyme": "BbsI",
        "matrix": "BbsI-HF",
        "role": "assemble module inserts into one block",
        "destination": ("CTCA", "CGAG"),
        "context_available": True,
    },
    "level1": {
        "enzyme": "BsaI",
        "matrix": "BsaI-HFv2",
        "role": "assemble blocks into the final construct",
        # The flanking ends are the compiled product's own termini, not a fixed vector pair:
        # the 5' end is whichever `1A` fusion-site variant the product uses, so a constant
        # here would be wrong for half the inventory. `("GGAG", "CGCT")` stood here as an
        # unverified MoClo default and never matched the deposited construct.
        "destination": None,
        "flanks_from_product": True,
        # The block *geometry* is established (LEVEL1_BLOCK_ENDS). The reaction's full
        # participant list is not: the deposited construct co-assembles parts this package
        # does not compile. Ligation fidelity depends on the whole competing overhang set, so
        # a PPR-only score is a diagnostic, not this reaction's fidelity -- see
        # LEVEL1_PARTICIPANTS_ESTABLISHED.
        "geometry_established": True,
        "context_available": False,
    },
}

#: The level-0 internal junctions the deposited kit uses, read off the primary cut geometry
#: rather than inferred from the insert sequences. Together with the level-0 destination they
#: are the complete level-0 reaction set: CTCA, ACTC, AAGA, GCAC, TGAA, CGAG.
LEVEL0_INTERNAL_JUNCTIONS: tuple[str, ...] = ("ACTC", "AAGA", "GCAC", "TGAA")

#: Joins *between* assembled blocks. These appear in a compiled product's sequence, but they
#: are not BbsI junctions and do not belong to any level-0 reaction. Treating them as level-0
#: ends produced a spurious "third reaction bottleneck" in the 19S architecture.
BLOCK_JOIN_OVERHANGS: tuple[str, ...] = ("CTTC", "GTGA", "CACG")

#: The level-1 cut geometry, derived from the 28 deposited block plasmids in
#: `9SDYW_level0.gb` and independently confirmed against Table S1's module interfaces.
#:
#: Each block plasmid carries exactly two BsaI sites and no BbsI site -- the mirror of the 42
#: module plasmids. Digesting them releases two fragment classes, 499 nt `AGGT -> CTTC` and
#: 406 nt `CTTC -> TTCG`. Joined at the shared overhang they give 901 nt, which is exactly the
#: length this package compiles for 9S, and the joined sequence occurs verbatim inside the
#: deposited level-1 product `9SrpoADYW_pICH47802_lc`.
#:
#: The consequence that makes level 1 scorable without any vector: **the released block
#: fragment is the joined module inserts and nothing else.** The block vector contributes no
#: bases, so a reaction's ends are the compiled product's own first and last interfaces,
#: which Table S1 supplies. The internal joins agree exactly with the module table --
#: `1E -> CTTC`, `14E -> GTGA`, `19E -> CACG`, `2E -> TTCG`, and `1A` entering at `AATG` or
#: `AGGT` depending on the variant chosen.
LEVEL1_BLOCK_ENDS: dict[str, object] = {
    "entry": ("AATG", "AGGT"),          # the two 1A fusion-site variants
    "exit": "TTCG",                     # 2E, in every architecture
    "internal": ("CTTC", "GTGA", "CACG"),
    "source": "9SDYW_level0.gb @ b569a8c4 (28 block plasmids) + Supplementary Table S1",
    "verified_by": ("fragment lengths 499 + 406 - 4 = 901 nt, the compiled 9S length; "
                    "the joined fragment occurs verbatim in 9SrpoADYW_pICH47802_lc"),
}

#: Whether the **complete** level-1 reaction can be scored. It cannot, and the reason is
#: narrower than the one this flag used to carry. The block release geometry is established
#: (`LEVEL1_BLOCK_ENDS`); what is missing is the reaction's *participant list*. The deposited
#: BsaI reaction also contains a P2L2S2 linker and the consensus DYW domain, and at least one
#: further part bridging `TTCG` to `TGTG` that the supplied records do not contain.
#:
#: This matters numerically, not just formally. Adding only the *known* omitted ends to the
#: 19S set moves the predicted fidelity from 0.996046 to **0.742067** -- below the level-0
#: figure of 0.751683. So a PPR-only score cannot be used to say which stage is the
#: bottleneck; on the one sensitivity calculation available it points the other way.
LEVEL1_CONTEXT_AVAILABLE = False
LEVEL1_PARTICIPANTS_ESTABLISHED = False

#: The parts the deposited level-1 reaction is known to also contain, read from the same
#: primary file as the block plasmids. Not a complete list -- `TTCG` to `TGTG` is unbridged --
#: which is exactly why the reaction stays unscored.
LEVEL1_KNOWN_COASSEMBLED_ENDS: tuple[str, ...] = ("TGTG", "CAAC", "GCTT")

#: What a level-1 score does **not** cover. The deposited construct co-assembles parts this
#: package does not compile -- a P2L2S2 linker (`TGTG -> CAAC`) and the consensus DYW domain
#: (`CAAC -> GCTT`), both in the same BsaI reaction, and at least one further part bridging
#: `TTCG` to `TGTG` that is absent from the supplied records. A level-1 fidelity reported here
#: is therefore for the **PPR block junctions only**, and adding an editing domain adds
#: overhangs that are not in it.
LEVEL1_COASSEMBLED_PARTS_UNMODELLED = True


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
