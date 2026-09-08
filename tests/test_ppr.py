"""Unit tests for target -> PPR code -> protein."""
from __future__ import annotations

import pytest

from clippr import constants as C
from clippr.ppr import (
    architecture_of,
    code_to_protein,
    code_to_rna,
    describe,
    normalize_target,
    rna_to_code,
)


class TestNormalize:
    def test_accepts_dna_and_rna(self):
        assert normalize_target("ttacacgtg") == "UUACACGUG"
        assert normalize_target("UUACACGUG") == "UUACACGUG"

    def test_strips_whitespace(self):
        assert normalize_target(" UUA CAC GUG ") == "UUACACGUG"

    @pytest.mark.parametrize("bad", ["", "UUACACGU", "UUACACGUGA", "N" * 9])
    def test_rejects_bad_input(self, bad):
        with pytest.raises(ValueError):
            normalize_target(bad)

    def test_rejects_non_acgu(self):
        with pytest.raises(ValueError, match="non-ACGU"):
            normalize_target("UUACACGUX")


class TestCode:
    def test_round_trip(self):
        for seq in ["UUACACGUG", "ACGUACGUA", "A" * 14, "GCUAAAGACUUGCAGGCUA"]:
            assert code_to_rna(rna_to_code(seq)) == normalize_target(seq)

    def test_code_pairs_are_known(self):
        for code in rna_to_code("ACGUACGUA"):
            assert code in C.CODE_TO_BASE

    def test_published_mapping(self):
        # Farley et al. Fig 1B: Thr/Asn->A, Asn/Asn->C, Thr/Asp->G, Asn/Asp->U
        assert C.CODE_TO_BASE == {"TN": "A", "NN": "C", "TD": "G", "ND": "T"}


class TestProtein:
    def test_one_repeat_per_base(self):
        for seq, n in [("UUACACGUG", 9), ("A" * 14, 14), ("C" * 19, 19)]:
            assert describe(seq)["n_arelf"] == n

    def test_length_is_linear_in_target(self):
        a, b = describe("A" * 9)["aa_length"], describe("A" * 14)["aa_length"]
        assert b - a == 5 * len(C.REPEAT_TEMPLATE.format(fifth="X", last="Y"))

    def test_starts_with_methionine(self):
        assert describe("UUACACGUG")["aa_sequence"].startswith("M")

    def test_start_codon_can_be_omitted(self):
        with_m = describe("UUACACGUG", include_start_codon=True)["aa_sequence"]
        without = describe("UUACACGUG", include_start_codon=False)["aa_sequence"]
        assert with_m == "M" + without

    def test_architecture_detection(self):
        assert architecture_of("UUACACGUG") == "9S"
        assert architecture_of("A" * 14) == "14S"
        assert architecture_of("A" * 19) == "19S"

    def test_unknown_code_rejected(self):
        with pytest.raises(ValueError, match="unknown PPR code"):
            code_to_protein(["ZZ"])

    def test_deterministic(self):
        assert describe("GCUAAAGAC") == describe("GCUAAAGAC")
