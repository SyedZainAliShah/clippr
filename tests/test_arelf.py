"""Unit tests for cut placement and overhang reachability."""
from __future__ import annotations

import pytest

from clippr import constants as C
from clippr.arelf import (
    achievable_overhangs,
    candidate_cuts,
    cut_context,
    find_motifs,
    offsets_from_motif,
    tunable_cuts,
)
from clippr.ppr import describe

#: Residue offsets from the preceding ARELF at which the corpus designs cut.
#: Measured across all 200 stored designs by validation/compare_cuts.py.
REFERENCE_OFFSETS = [3, 5, 7, 8, 13, 16, 17, 18, 21, 22, 28, 29]


@pytest.fixture(scope="module")
def protein():
    return describe("AAAAUGUGG")["aa_sequence"]


class TestMotifs:
    def test_one_motif_per_repeat(self, protein):
        assert len(find_motifs(protein)) == 9

    def test_motifs_are_evenly_spaced(self, protein):
        m = find_motifs(protein)
        assert {b - a for a, b in zip(m, m[1:])} == {len(C.REPEAT_TEMPLATE.format(
            fifth="X", last="Y"))}

    def test_each_hit_really_is_the_motif(self, protein):
        assert all(protein[i:i + 5] == "ARELF" for i in find_motifs(protein))

    def test_absent_motif_returns_empty(self):
        assert find_motifs("MMMMMMMM") == []


class TestCutContext:
    def test_is_the_four_bases_before_the_cut(self):
        cds = "AAACCCGGGTTTAAACCC"
        assert cut_context(cds, 4) == cds[8:12]

    def test_offset_shifts_by_codons(self):
        cds = "AAACCCGGGTTTAAACCC"
        assert cut_context(cds, 2, 2) == cut_context(cds, 4)

    def test_spans_a_codon_boundary(self):
        """The overhang is the wobble base of one codon plus the whole of the next."""
        cds = "AAA" "CCC" "GGG" "TTT"
        assert cut_context(cds, 3) == "C" + "GGG"

    @pytest.mark.parametrize("cut", [0, 1, 100])
    def test_out_of_range_rejected(self, cut):
        with pytest.raises(ValueError):
            cut_context("AAACCCGGGTTT", cut)


class TestAchievable:
    def test_contains_the_actual_context(self, protein):
        """Whatever a cut currently reads must be one of the options it can reach."""
        from clippr.arelf import cut_context as ctx
        cds = "".join(C.SYNONYMOUS_CODONS[a][0] for a in protein)
        for cut in (50, 88, 171, 238):
            assert ctx(cds, cut) in achievable_overhangs(protein, cut)

    def test_counts_follow_the_genetic_code(self, protein):
        """Options = (distinct wobble bases of residue c-2) x (codons of residue c-1)."""
        for cut in (73, 78, 88, 98):
            a, b = protein[cut - 2], protein[cut - 1]
            expect = len({c[2] for c in C.SYNONYMOUS_CODONS[a]}) * len(C.SYNONYMOUS_CODONS[b])
            assert len(achievable_overhangs(protein, cut)) == expect

    def test_methionine_is_the_narrow_case(self, protein):
        """M has one codon, so a cut just after it offers only wobble variation."""
        cut = next(c for c in range(2, len(protein)) if protein[c - 1] == "M")
        assert len(achievable_overhangs(protein, cut)) <= 4

    def test_all_options_are_well_formed(self, protein):
        for o in achievable_overhangs(protein, 88):
            assert len(o) == 4 and set(o) <= set("ACGT")

    def test_needs_two_preceding_residues(self, protein):
        with pytest.raises(ValueError):
            achievable_overhangs(protein, 1)


def _brute_force_feasible(protein, cut_aa, overhang, enzymes):
    """Independent oracle: exhaustively enumerate the local synonymous space.

    Deliberately written without reusing any of safe_overhangs' logic, so agreement
    between the two is evidence rather than tautology.
    """
    from Bio import Restriction
    rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    sites = []
    for name in C.enzymes_for(enzymes):
        s = str(getattr(Restriction, name).site).upper()
        sites += [s] if rc(s) == s else [s, rc(s)]

    left_res = protein[cut_aa - 2]
    right_res = protein[cut_aa] if cut_aa < len(protein) else None
    lefts = [c for c in C.SYNONYMOUS_CODONS[left_res] if c[2] == overhang[0]]
    rights = list(C.SYNONYMOUS_CODONS[right_res]) if right_res else [""]
    if not lefts:
        return False
    # the middle codon is pinned by the overhang; it must code the right residue
    from Bio.Seq import Seq
    if str(Seq(overhang[1:]).translate()) != protein[cut_aa - 1]:
        return False
    for L in lefts:
        for R in rights:
            window = L + overhang[1:] + R
            if not any(s in window for s in sites):
                return True
    return False


class TestSafeOverhangsProperty:
    """Soundness and completeness against brute-force local enumeration.

    This is the test that pins the bug class, not one instance: every accepted overhang
    must have a demonstrable realization, and every rejected one must have none.
    """

    @pytest.mark.parametrize("profile", ["none", "assembly", "igem_rfc1000", "moclo_compat"])
    def test_sound_and_complete_on_every_residue_pair(self, profile):
        from clippr.arelf import safe_overhangs
        residues = sorted(C.SYNONYMOUS_CODONS)
        checked = 0
        for a in residues:
            for b in residues:
                prot = "M" + a + b + "A"
                cut = 3
                accepted = set(safe_overhangs(prot, cut, enzymes=profile))
                for oh in achievable_overhangs(prot, cut):
                    checked += 1
                    want = _brute_force_feasible(prot, cut, oh, profile)
                    assert (oh in accepted) is want, (
                        f"{profile}: residues {a}{b}, overhang {oh}: "
                        f"filter says {oh in accepted}, brute force says {want}")
        assert checked > 2000, f"only {checked} combinations exercised"

    def test_agrees_on_a_real_protein(self):
        from clippr.arelf import safe_overhangs
        prot = describe("AAAGCGGCACUUGUGAAGU")["aa_sequence"]
        for cut in range(2, len(prot), 7):
            accepted = set(safe_overhangs(prot, cut))
            for oh in achievable_overhangs(prot, cut):
                assert (oh in accepted) is _brute_force_feasible(
                    prot, cut, oh, C.DEFAULT_ENZYME_PROFILE)

    def test_accepted_is_always_a_subset_of_achievable(self):
        from clippr.arelf import safe_overhangs
        prot = describe("AAAAUGUGG")["aa_sequence"]
        for cut in range(2, len(prot) + 1):
            assert set(safe_overhangs(prot, cut)) <= set(achievable_overhangs(prot, cut))

    def test_a_missing_wobble_codon_rejects(self):
        """The empty-lefts case: 'impossible' must never become 'unconstrained'.

        Methionine has one codon, ATG, so no M codon ends in C. An overhang claiming a
        C wobble from an M at cut-2 has no realization and must be rejected even though
        the surrounding sequence is entirely free of restriction sites.
        """
        from clippr.arelf import safe_overhangs
        prot = "AMKAAA"
        cut = 3                                   # residues cut-2, cut-1 are M, K
        assert prot[cut - 2] == "M" and prot[cut - 1] == "K"
        assert not [c for c in C.SYNONYMOUS_CODONS["M"] if c[2] == "C"]
        assert _brute_force_feasible(prot, cut, "CAAA", "none") is False
        assert "CAAA" not in safe_overhangs(prot, cut, enzymes="none")


class TestSafeOverhangs:
    """Regression cover for the failure that broke 13 of 50 19S designs."""

    def test_subset_of_achievable(self, protein):
        from clippr.arelf import safe_overhangs
        for cut in (73, 88, 175, 262):
            assert set(safe_overhangs(protein, cut)) <= set(
                achievable_overhangs(protein, cut))

    def test_rejects_an_overhang_that_forces_a_site(self):
        """The real corpus case: AGAC at Arg-Asp, with the site on the reverse strand.

        AGAC fixes the Asp codon to GAC and forces the Arg codon to end in A -- both
        options (AGA, CGA) give ...GAGAC. When the next residue's codons all begin with C
        or G, as alanine's do, the result is GAGACC or GAGACG: BsaI or BsmBI reverse
        complement. A forward-only scan misses this entire class.
        """
        from clippr.arelf import safe_overhangs
        prot = "MAAAAARDAAAAA"
        cut = prot.index("D") + 1
        assert prot[cut - 2] == "R" and prot[cut - 1] == "D"
        assert "AGAC" in achievable_overhangs(prot, cut)
        # only BsmBI can fire here -- alanine's codons all begin with G, so GAGACG is
        # reachable but GAGACC is not. The exclusion therefore needs the strict profile.
        assert "AGAC" not in safe_overhangs(prot, cut, enzymes="moclo_compat")

    def test_the_forced_site_really_exists_on_a_built_sequence(self):
        """Verify by construction with Biopython, not by reasoning about the code."""
        from Bio import Restriction
        from Bio.Seq import Seq
        prot = "MAAAAARDAAAAA"
        cut = prot.index("D") + 1                      # residues cut-2, cut-1 are R, D
        cds = ("".join(C.SYNONYMOUS_CODONS[a][0] for a in prot[:cut - 2])
               + "CGA" + "GAC"
               + "".join(C.SYNONYMOUS_CODONS[a][0] for a in prot[cut:]))
        assert str(Seq(cds).translate()) == prot
        assert cds[3 * cut - 4:3 * cut] == "AGAC"
        assert Restriction.BsmBI.search(Seq(cds)), "the forced site is not actually there"

    def test_an_escape_keeps_the_overhang_available(self):
        """Same overhang, but a next residue whose codons can start with A or T."""
        from clippr.arelf import safe_overhangs
        prot = "MAAAAARDMAAAA"          # methionine follows: ATG starts with A
        cut = prot.index("D") + 1
        assert "AGAC" in safe_overhangs(prot, cut)

    def test_keeps_overhangs_with_an_escape(self, protein):
        from clippr.arelf import safe_overhangs
        assert safe_overhangs(protein, 88), "no safe overhang at a routine cut"

    @pytest.mark.parametrize("profile,expect", [
        ("moclo_compat", False),   # BsmBI excluded, so AGAC is infeasible here
        ("igem_rfc1000", True),    # the default; BsmBI allowed, and alanine's codons all
                                   # begin with G, so only GAGACG could fire -- AGAC stays
        ("assembly", True),
        ("none", True),
    ])
    def test_profile_changes_what_is_feasible(self, profile, expect):
        """Which sites are excluded is a declared policy, not a universal fact.

        This case is exactly the trade to decide deliberately: BsmBI is excluded only to
        keep later MoClo levels open, and excluding it is what makes AGAC infeasible here.
        """
        from clippr.arelf import safe_overhangs
        prot = "MAAAAARDAAAAA"
        cut = prot.index("D") + 1
        assert ("AGAC" in safe_overhangs(prot, cut, enzymes=profile)) is expect

    def test_profiles_are_ordered_by_strictness(self):
        p = C.ENZYME_PROFILES
        assert (set(p["none"]) < set(p["assembly"]) < set(p["igem_rfc1000"])
                < set(p["moclo_compat"]))

    def test_every_enzyme_has_a_declared_role(self):
        assert set(C.ENZYME_ROLES) == set(C.ENZYME_PROFILES["moclo_compat"])

    def test_default_is_what_the_project_knows_it_needs(self):
        """Not the strictest: BsmBI exclusion is speculative and is opt-in."""
        assert C.DEFAULT_ENZYME_PROFILE == "igem_rfc1000"
        assert "BsmBI" not in C.enzymes_for()
        assert C.ENZYME_ROLES["BsmBI"] == "downstream"

    def test_unknown_profile_rejected(self):
        with pytest.raises(ValueError, match="unknown enzyme profile"):
            C.enzymes_for("whatever")


class TestTunableCuts:
    def test_spans_the_reference_offsets(self, protein):
        """Our criterion must not exclude positions the corpus designs actually use."""
        assert set(REFERENCE_OFFSETS) <= set(offsets_from_motif(protein))

    def test_default_only_excludes_dead_positions(self, protein):
        cuts = tunable_cuts(protein)
        assert all(len(achievable_overhangs(protein, c)) >= 2 for c in cuts)

    def test_raising_the_threshold_narrows(self, protein):
        assert len(tunable_cuts(protein, 12)) < len(tunable_cuts(protein, 2))


class TestCandidateCuts:
    def test_fragments_respect_the_bounds(self, protein):
        n = len(protein)
        for cuts in candidate_cuts(protein, 4, 65, 90)[:50]:
            bounds = [0] + cuts + [n]
            assert all(65 <= b - a <= 90 for a, b in zip(bounds, bounds[1:]))

    def test_cuts_are_ascending_and_distinct(self, protein):
        for cuts in candidate_cuts(protein, 4, 65, 90)[:50]:
            assert cuts == sorted(set(cuts))

    def test_single_fragment_needs_no_cuts(self, protein):
        assert candidate_cuts(protein, 1, 0, 10_000) == [[]]

    def test_impossible_split_rejected(self, protein):
        with pytest.raises(ValueError, match="cannot split"):
            candidate_cuts(protein, 4, 200, 250)

    def test_refuses_to_truncate_silently(self, protein):
        """A truncated list looks like a sample but is biased towards early cuts."""
        with pytest.raises(ValueError, match="narrow the fragment bounds"):
            candidate_cuts(protein, 4, 40, 120, max_results=100)

    def test_prunes_instead_of_filtering(self, protein):
        """Tight bounds must return quickly, not walk every combination first."""
        import time
        t0 = time.time()
        out = candidate_cuts(protein, 4, 74, 78)
        assert time.time() - t0 < 5.0
        assert out


class TestBalancedCuts:
    def test_returns_most_balanced_first(self, protein):
        from clippr.arelf import balanced_cuts
        splits = balanced_cuts(protein, 4, 40, 120, limit=20)
        n = len(protein)

        def spread(cuts):
            b = [0] + cuts + [n]
            return max(y - x for x, y in zip(b, b[1:]))

        assert spread(splits[0]) <= spread(splits[-1])

    def test_respects_the_limit(self, protein):
        from clippr.arelf import balanced_cuts
        assert len(balanced_cuts(protein, 4, 40, 120, limit=7)) <= 7

    def test_all_splits_are_legal(self, protein):
        from clippr.arelf import balanced_cuts
        n = len(protein)
        for cuts in balanced_cuts(protein, 4, 60, 95, limit=20):
            b = [0] + cuts + [n]
            assert all(60 <= y - x <= 95 for x, y in zip(b, b[1:]))
            assert cuts == sorted(set(cuts))

    def test_scales_to_the_largest_architecture(self):
        """19S split seven ways is where exhaustive enumeration became infeasible."""
        import time
        from clippr.arelf import balanced_cuts
        big = describe("AAAGCGGCACUUGUGAAGU")["aa_sequence"]
        t0 = time.time()
        splits = balanced_cuts(big, 7, 40, 90, limit=40)
        assert time.time() - t0 < 5.0
        assert splits
