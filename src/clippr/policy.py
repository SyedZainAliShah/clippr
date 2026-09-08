"""Layer 3 of 3 — what CLIPPR *chose*. None of this is required by GRASP or by biology.

Everything here is a project decision: which enzyme sites to exclude, which host to
optimise for, where to put a QC threshold. Each is defensible and each is documented with
its reason, but a future contributor must be able to tell them apart from published fact.
That is the whole point of the split:

    biology.py         published fact
    assembly_spec.py   assembly mechanics
    policy.py          this project's choices   <- you are here

**Why this file exists.** A constraint invented beyond every available source needs its own
justification. Excluding BsmBI was exactly that -- it appears in the paper, in the recovered
reference corpus, and in the iGEM standard precisely nowhere -- and it made 13 of 50
19S designs infeasible in our own pipeline until `arelf.safe_overhangs` was added. Nothing
in this file should be mistaken for something GRASP biologically requires.
"""
from __future__ import annotations

#: Why each Type IIS site is excluded. These are different kinds of claim and collapsing
#: them into one blacklist hides an assumption.
#:
#: - ``assembly``      the enzyme this assembly actually cuts with. An internal site would
#:   be cut during the reaction, so its absence is chemistry, not preference.
#: - ``igem_rfc1000``  required by iGEM's Type IIS assembly standard RFC[1000], which
#:   mandates that BsaI and SapI sites be absent from participating parts. A real
#:   requirement, but one that follows from a declared standard rather than the reaction.
#: - ``downstream``    kept clear only so later MoClo levels stay available. A preference
#:   about hypothetical future use, and the weakest of the three.
ENZYME_ROLES: dict[str, str] = {
    "BsaI": "assembly",
    "BbsI": "assembly",
    "SapI": "igem_rfc1000",
    "BsmBI": "downstream",
}

#: Named constraint profiles, narrowest to widest. Each is cumulative over the last, but
#: the reason each enzyme joins differs, so read the roles above rather than the tuple.
#:
#: ``assembly`` -- BsaI and BbsI, and both belong on evidence, not assumption. BsaI wraps
#: every fragment (it is the `wrap_enzyme` on all 950 stored oligos) and BbsI is the level 0
#: ligation enzyme, which is why this project scores the level 0 destination pair against
#: the BbsI-HF matrix. Two enzymes because the assembly is hierarchical, not because one
#: standard demands both.
#:
#: ``igem_rfc1000`` -- adds SapI. RFC[1000]'s own formal requirement is that **BsaI and
#: SapI** sites be absent; it does not make BbsI illegal. BbsI is in this profile because
#: CLIPPR's assembly uses it, not because the standard asks for it. The profile is the
#: union of two independent requirements, not an inheritance.
#:
#: ``moclo_compat`` -- adds BsmBI, purely to keep later MoClo-family levels available.
#: Opt-in, and the name says *a* cross-compatibility choice: MoClo and GoldenBraid
#: families use different Type IIS enzymes, so there is no single "MoClo set".
#:
#: **Keep the three sources of evidence separate:**
#:
#:   * the **paper** establishes GRASP's assembly architecture;
#:   * the **reference corpus** records what its generator actually
#:     excluded, namely "BsaI, SapI, BpiI" (BpiI is an isoschizomer of BbsI, same GAAGAC
#:     site) -- an artifact of a lawfully retained copy, not a published standard;
#:   * **iGEM documentation** establishes the RFC[1000] requirement, naming BsaI and SapI.
#:
#: BsmBI appears in none of the three. Choose the profile from the actual assembly
#: specification once the backbone is fixed.
ENZYME_PROFILES: dict[str, tuple[str, ...]] = {
    "none": (),
    "assembly": ("BsaI", "BbsI"),
    "igem_rfc1000": ("BsaI", "BbsI", "SapI"),
    "moclo_compat": ("BsaI", "BbsI", "SapI", "BsmBI"),
}

DEFAULT_ENZYME_PROFILE: str = "igem_rfc1000"


def enzymes_for(profile: str | tuple[str, ...] = DEFAULT_ENZYME_PROFILE) -> tuple[str, ...]:
    """Resolve a profile name, or pass an explicit tuple of enzyme names through."""
    if not isinstance(profile, str):
        return tuple(profile)
    try:
        return ENZYME_PROFILES[profile]
    except KeyError as e:
        raise ValueError(
            f"unknown enzyme profile {profile!r}; have {sorted(ENZYME_PROFILES)}"
        ) from e


#: The host this project designs for. Chlamydomonas reinhardtii, NCBI taxonomy id.
#: A project choice: the wet lab's organism, not something GRASP requires.
CHLAMYDOMONAS_TAXID: int = 3055

#: Named hosts -> (Kazusa taxid, NCBI genetic code). Only the verified entry is here; a
#: guessed chloroplast taxid would silently produce wrong DNA.
ORGANISMS: dict[str, tuple[int, int]] = {
    "c_reinhardtii_nuclear": (CHLAMYDOMONAS_TAXID, 1),
}

#: Longest and shortest fragment to order, in residues. 90 aa = 270 nt of coding plus
#: overhangs and arms. Reproduces the corpus splits: 302 aa -> 4, 612 aa -> 7.
MAX_FRAGMENT_AA: int = 90
MIN_FRAGMENT_AA: int = 40

#: k for the sequence-uniqueness objective and the QC repeat metric. 20 because it
#: dominated 15 on measurement: same duplicate count, less codon adaptation surrendered.
REPEAT_K: int = 20
