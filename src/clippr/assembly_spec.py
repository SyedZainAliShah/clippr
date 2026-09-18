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


#: The reference implementation's composition rule for a level-1 reaction, recorded because it
#: is the structure our own derivation was missing -- not because it is complete. Reading
#: `grasp_level1_reaction_overhangs` shows the set assembled as:
#:
#:     [final_cassette_5, ppr_outer_5] + block_joins + [ppr_outer_3, final_cassette_3]
#:
#: with both flanking pairs required and no defaults. Their `final_cassette` is the MoClo
#: level-1 acceptor (`GGAG`/`CGCT`) and their `ppr_outer` is `AGGT`/`TTCG`.
#:
#: **They omit the linker and DYW domain too**, and report a `level1_fidelity` from the
#: incomplete set regardless. The composition rule is worth having; their completeness claim
#: is not.
LEVEL1_REFERENCE_COMPOSITION = {
    "rule": "final_cassette_5, ppr_outer_5, *block_joins, ppr_outer_3, final_cassette_3",
    "final_cassette": ("GGAG", "CGCT"),
    "ppr_outer": ("AGGT", "TTCG"),
    "source": "grasp-library-designer ligation_fidelity.grasp_level1_reaction_overhangs",
    "omits": ("P2L2S2 linker TGTG->CAAC", "consensus DYW domain CAAC->GCTT",
              "at least one part bridging TTCG to TGTG"),
}


#: Roles a complete level-1 reaction must account for. Derived from the deposited construct:
#: the PPR blocks sit inside a final cassette, and the same BsaI reaction carries a P2L2S2
#: linker and the consensus DYW domain. The bridge between `TTCG` and `TGTG` is demanded
#: explicitly because the supplied records do not contain it -- so completeness is currently
#: unreachable, and that is the honest state rather than a hidden one.
LEVEL1_REQUIRED_ROLES: tuple[str, ...] = (
    "ppr_blocks", "final_cassette", "linker", "editing_domain", "bridge",
)


def level1_participants(block_subset, *, roles=None, evidence: str = "") -> dict:
    """Assemble a level-1 participant list, and say whether it is **complete**.

    `roles` maps a role from `LEVEL1_REQUIRED_ROLES` to the overhangs it contributes. A
    reaction is complete only when every required role is present, every contributed end is a
    four-base overhang, and the ends form a single closed chain -- each internal overhang
    appearing exactly twice, and exactly two ends appearing once.

    `evidence` is required but **not sufficient**. An earlier version set completeness from
    `bool(extra) and bool(evidence)`, so one overhang and a sentence saying "DYW and bridging
    parts are missing" returned `participants_established=True`. A citation field is not a
    completeness check, and prose cannot discharge a contract.

    Nothing here is scored. The caller gets the set, an explicit verdict and, when incomplete,
    the roles it is missing.
    """
    supplied = {str(role): [str(o).upper().replace("U", "T") for o in (ends or ())]
                for role, ends in (roles or {}).items()}
    subset = [str(o).upper().replace("U", "T") for o in block_subset]
    supplied.setdefault("ppr_blocks", subset)

    missing = [role for role in LEVEL1_REQUIRED_ROLES if not supplied.get(role)]
    malformed = sorted({o for ends in supplied.values() for o in ends
                        if len(o) != 4 or set(o) - set("ACGT")})

    overhangs: list[str] = []
    for role in LEVEL1_REQUIRED_ROLES:
        overhangs.extend(supplied.get(role, []))
    for role, ends in sorted(supplied.items()):
        if role not in LEVEL1_REQUIRED_ROLES:
            overhangs.extend(ends)

    # A reaction is a chain: internal junctions are shared by two parts, the two outermost
    # ends by one each. Anything else is not a set of fragments that can assemble in order.
    counts: dict[str, int] = {}
    for overhang in overhangs:
        counts[overhang] = counts.get(overhang, 0) + 1
    terminal = sorted(o for o, n in counts.items() if n == 1)
    overused = sorted(o for o, n in counts.items() if n > 2)
    chain_ok = len(terminal) == 2 and not overused

    complete = bool(evidence) and not missing and not malformed and chain_ok
    reasons = []
    if not evidence:
        reasons.append("no evidence supplied naming where these participants came from")
    if missing:
        reasons.append(f"no participant supplied for {', '.join(missing)}")
    if malformed:
        reasons.append(f"not four-base overhangs: {', '.join(malformed)}")
    if not chain_ok and not missing:
        reasons.append(
            f"the ends do not form one chain: {len(terminal)} appear once "
            f"(expected 2)" + (f", and {', '.join(overused)} appear more than twice"
                               if overused else ""))

    return {
        "overhangs": overhangs,
        "block_subset": subset,
        "roles_supplied": {role: list(ends) for role, ends in sorted(supplied.items())},
        "roles_required": list(LEVEL1_REQUIRED_ROLES),
        "roles_missing": missing,
        "participants_established": complete,
        "evidence": evidence or None,
        "why_incomplete": None if complete else "; ".join(reasons),
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
