"""Unit tests for ligation-fidelity scoring and overhang selection."""
from __future__ import annotations

import pytest

from clippr.overhangs import (
    best_set,
    enumerate_candidates,
    fidelity_components,
    load_matrix,
    palindromic,
    reaction_overhangs,
    reverse_complement,
    set_fidelity,
    valid_set,
)

#: Totals of the published count tables. A parser that silently drops a row or reads a
#: header as data will not hit these.
MATRIX_SUMS = {"BsaI-HFv2": 203364, "BbsI-HF": 189282}


class TestMatrix:
    @pytest.mark.parametrize("name,total", MATRIX_SUMS.items())
    def test_shape_and_total(self, name, total):
        counts, labels = load_matrix(name)
        assert counts.shape == (256, 256)
        assert len(labels) == 256
        assert counts.sum() == total

    def test_labels_are_every_four_mer(self):
        _, labels = load_matrix("BsaI-HFv2")
        assert len(set(labels)) == 256
        assert all(len(o) == 4 and set(o) <= set("ACGT") for o in labels)

    def test_correct_pairing_is_the_diagonal(self):
        """Row i is labelled revcomp(labels[i]), so on-target counts lie on the diagonal."""
        counts, labels = load_matrix("BsaI-HFv2")
        for o in ("ACAT", "CTCA", "GGAG"):
            i = labels.index(o)
            assert counts[i, i] == max(counts[:, i]), f"{o} is not its own best partner"


class TestSetFidelity:
    @pytest.mark.parametrize("ohs,matrix,want", [
        (["ACAT", "ACAA"], "BsaI-HFv2", 0.938296),   # level -1 destination pair
        (["CTCA", "CTCG"], "BbsI-HF", 0.794433),     # level 0 destination pair
        (["GGAG", "AGCG"], "BsaI-HFv2", 1.000000),   # level 1, MoClo standard sites
    ])
    def test_published_pairs(self, ohs, matrix, want):
        assert set_fidelity(ohs, matrix) == pytest.approx(want, abs=1e-6)

    def test_orientation_invariant(self):
        s = ["ACAT", "CTCA", "GGAG"]
        flipped = [reverse_complement(o) for o in s]
        assert set_fidelity(s) == pytest.approx(set_fidelity(flipped), abs=1e-12)
        one = [reverse_complement(s[0])] + s[1:]
        assert set_fidelity(s) == pytest.approx(set_fidelity(one), abs=1e-12)

    def test_bounded_in_unit_interval(self):
        for s in (["GGAG"], ["GGAG", "AGCG"], ["ACAT", "ACAA", "CTCA", "GGAG"]):
            assert 0.0 <= set_fidelity(s) <= 1.0

    def test_geometric_mean_of_components(self):
        import math
        s = ["ACAT", "ACAA"]
        fwd, rev = fidelity_components(s)
        assert set_fidelity(s) == pytest.approx(math.sqrt(fwd * rev), abs=1e-12)

    def test_more_overhangs_cannot_help(self):
        """Adding a competitor can only take probability away."""
        base = set_fidelity(["GGAG", "AGCG"])
        assert set_fidelity(["GGAG", "AGCG", "GGAA"]) <= base + 1e-12

    def test_matrices_disagree(self):
        """The two enzymes are different measurements; a shared cache would hide that."""
        assert set_fidelity(["CTCA", "CTCG"], "BsaI-HFv2") != set_fidelity(
            ["CTCA", "CTCG"], "BbsI-HF")

    @pytest.mark.parametrize("bad", [[], ["ACG"], ["ACGTA"], ["ACGN"], ["AC GT"]])
    def test_rejects_malformed(self, bad):
        with pytest.raises(ValueError):
            set_fidelity(bad)

    def test_accepts_rna_spelling(self):
        assert set_fidelity(["ACAU", "ACAA"]) == set_fidelity(["ACAT", "ACAA"])

    def test_unknown_matrix_rejected(self):
        with pytest.raises(ValueError, match="unknown matrix"):
            set_fidelity(["GGAG", "AGCG"], "NotAnEnzyme")


class TestValidity:
    def test_palindrome_detected(self):
        assert palindromic("AATT") and palindromic("GGCC")
        assert not palindromic("GGAG")

    @pytest.mark.parametrize("bad", [
        ["AATT", "GGAG"],      # palindrome
        ["GGAG", "GGAG"],      # repeat
        ["GGAG", "CTCC"],      # a member's reverse complement
    ])
    def test_invalid_sets(self, bad):
        assert not valid_set(bad)

    def test_valid_set(self):
        assert valid_set(["GGAG", "AGCG", "AATG"])


class TestReactionOverhangs:
    def test_flips_the_three_prime_coding_site(self):
        """The stored level 0 pair is ("CTCA","CGAG"); the strands present are CTCA/CTCG."""
        assert reaction_overhangs([], "level0") == ["CTCA", "CTCG"]

    def test_the_flip_changes_the_verdict(self):
        """Guards the specific false all-clear this helper exists to prevent."""
        from clippr import constants as C
        stored = list(C.DESTINATION_OVERHANGS["level0"])
        assert set_fidelity(stored, "BbsI-HF") == pytest.approx(1.0, abs=1e-9)
        assert set_fidelity(reaction_overhangs([], "level0"), "BbsI-HF") < 0.8

    def test_junctions_are_kept_in_order(self):
        got = reaction_overhangs(["CGCT", "AATG", "TGTT"], "level0")
        assert got == ["CTCA", "CGCT", "AATG", "TGTT", "CTCG"]

    def test_unknown_level_rejected(self):
        with pytest.raises(ValueError, match="unknown level"):
            reaction_overhangs([], "level7")


class TestSelection:
    def test_enumerate_drops_palindromes(self):
        assert enumerate_candidates(["AATT,GGAG,GGCC"]) == [["GGAG"]]

    def test_enumerate_deduplicates(self):
        assert enumerate_candidates(["GGAG,GGAG,AATG"]) == [["GGAG", "AATG"]]

    def test_enumerate_accepts_lists(self):
        assert enumerate_candidates([["GGAG", "AATG"]]) == [["GGAG", "AATG"]]

    def test_enumerate_rejects_all_palindromic(self):
        with pytest.raises(ValueError, match="palindromic"):
            enumerate_candidates(["AATT,GGCC"])

    def test_best_set_picks_the_maximum(self):
        candidates = [["GGAA", "AGCG"], ["AATG"]]
        chosen, score = best_set(candidates, fixed=["GGAG"])
        assert chosen == ["AGCG", "AATG"]
        assert score == pytest.approx(set_fidelity(["AGCG", "AATG", "GGAG"]), abs=1e-12)

    def test_best_set_honours_fixed_overhangs(self):
        """A candidate identical to a fixed overhang is not a legal choice."""
        chosen, _ = best_set([["GGAG", "AGCG"]], fixed=["GGAG"])
        assert chosen == ["AGCG"]

    def test_best_set_refuses_to_guess_when_huge(self):
        with pytest.raises(ValueError, match="max_combinations"):
            best_set([["ACGT"] * 40] * 5, max_combinations=100)

    def test_best_set_reports_no_valid_combination(self):
        with pytest.raises(ValueError, match="no valid overhang combination"):
            best_set([["GGAG"], ["GGAG"]])
