"""Unit tests for two-tier cross-talk assessment.

The load-bearing property is not that the affinity model is right -- its applicability to
this scaffold is unestablished -- but that it can never gate a decision. Tier A alone
decides; tier B may only add "review".
"""
from __future__ import annotations

import pytest

from clippr.crosstalk import (
    PairRisk,
    assess,
    compare_tiers,
    load_ppr_scores,
    predicted_affinity,
    report,
)

#: A minimal table in the PPRmatcher layout, covering the four codes GRASP uses. Written
#: inline so the suite needs no third-party file -- the real table is not distributed.
TABLE = "\n".join([
    "For motif types: P",
    "5th/last\tA\tC\tG\tU",
    "TN\t0.70\t-0.51\t0.53\t-0.22",
    "NN\t-0.44\t0.77\t-0.07\t0.62",
    "TD\t0.01\t-1.0\t0.84\t-0.02",
    "ND\t-0.53\t0.66\t0.55\t0.83",
])


@pytest.fixture(scope="module")
def scores(tmp_path_factory):
    p = tmp_path_factory.mktemp("tbl") / "scores.tsv"
    p.write_text(TABLE, encoding="utf-8")
    return load_ppr_scores(p)


class TestLoading:
    def test_parses_the_documented_layout(self, scores):
        assert set(scores) == {"TN", "NN", "TD", "ND"}
        assert set(scores["TN"]) == set("ACGT")

    def test_u_is_stored_as_t(self, scores):
        assert "T" in scores["ND"] and "U" not in scores["ND"]

    def test_cognate_base_scores_highest_for_every_grasp_code(self, scores):
        """Independent corroboration that the table is sane for these four codes.

        The PPR code says TN reads A, NN reads C, TD reads G, ND reads U. A scoring table
        derived from separate experiments agreeing on all four is worth asserting.
        """
        for code, base in (("TN", "A"), ("NN", "C"), ("TD", "G"), ("ND", "T")):
            assert max(scores[code], key=scores[code].get) == base

    def test_missing_header_rejected(self, tmp_path):
        p = tmp_path / "bad.tsv"
        p.write_text("nothing useful here\n", encoding="utf-8")
        with pytest.raises(ValueError, match="no header row"):
            load_ppr_scores(p)

    def test_wrong_bases_rejected(self, tmp_path):
        p = tmp_path / "bad.tsv"
        p.write_text("x\n5th/last\tA\tC\tG\tX\nTN\t1\t1\t1\t1\n", encoding="utf-8")
        with pytest.raises(ValueError, match="expected the four bases"):
            load_ppr_scores(p)


class TestAffinity:
    def test_a_target_scores_itself_at_one(self, scores):
        assert predicted_affinity("AAAAUGUGG", "AAAAUGUGG", scores) == pytest.approx(1.0)

    def test_a_different_target_scores_lower(self, scores):
        assert predicted_affinity("AAAAUGUGG", "CCCCCCCCC", scores) < 1.0

    def test_is_directional(self, scores):
        """The PPR for A against B is not the same question as the PPR for B against A."""
        a, b = "AAAAUGUGG", "GCUAAAGAC"
        assert predicted_affinity(a, b, scores) != predicted_affinity(b, a, scores)

    def test_length_mismatch_rejected(self, scores):
        with pytest.raises(ValueError, match="cannot score"):
            predicted_affinity("AAAAUGUGG", "AAAAUGUGGAAAAA", scores)

    def test_missing_code_rejected(self, tmp_path):
        p = tmp_path / "part.tsv"
        p.write_text("x\n5th/last\tA\tC\tG\tU\nTN\t1\t0\t0\t0\n", encoding="utf-8")
        with pytest.raises(ValueError, match="no entry for code pair"):
            predicted_affinity("AAAAUGUGG", "AAAAUGUGG", load_ppr_scores(p))


class TestTierSeparation:
    """Tier B must never gate. This is the property that makes the model safe to include."""

    def test_hamming_alone_decides_too_close(self, scores):
        pair = PairRisk("AAAAUGUGG", "AAAAUGUGA", 1, 0.1, 0.1)
        assert pair.flag == "TOO CLOSE", "low affinity must not rescue a close pair"

    def test_affinity_can_only_add_review(self, scores):
        pair = PairRisk("AAAAUGUGG", "CCCCCCCCC", 9, 0.99, 0.99)
        assert pair.flag == "review", "high affinity must not fail a well-separated pair"

    def test_a_clean_pair_is_ok(self):
        assert PairRisk("AAAAUGUGG", "CCCCCCCCC", 9, 0.1, 0.1).flag == "ok"

    def test_tier_a_works_without_any_table(self):
        pairs = assess(["AAAAUGUGG", "CCCCCCCCC", "GGGGGGGGG"])
        assert pairs and all(p.worst_affinity is None for p in pairs)

    def test_report_says_when_tier_b_is_absent(self):
        text = report(["AAAAUGUGG", "CCCCCCCCC"])
        assert "no scoring table supplied" in text

    def test_report_carries_the_applicability_caveat(self, scores):
        text = report(["AAAAUGUGG", "CCCCCCCCC"], scores)
        assert "P-type" in text and "S-type" in text
        assert "gates nothing" in text


class TestAssess:
    def test_orders_closest_first(self, scores):
        pairs = assess(["AAAAUGUGG", "AAAAUGUGA", "CCCCCCCCC"], scores)
        assert [p.hamming for p in pairs] == sorted(p.hamming for p in pairs)

    def test_skips_pairs_of_unequal_length(self, scores):
        pairs = assess(["AAAAUGUGG", "GCUAAAGACUUGCA"], scores)
        assert pairs == []

    def test_every_same_length_pair_appears(self, scores):
        targets = ["AAAAUGUGG", "CCCCCCCCC", "GGGGGGGGG", "UUUUUUUUU"]
        assert len(assess(targets, scores)) == 6


class TestCompareTiers:
    def test_reports_whether_the_model_changes_anything(self, scores):
        """The honest use of a model whose applicability is unproven."""
        out = compare_tiers(["AAAAUGUGG", "CCCCCCCCC", "GGGGGGGGG"], scores)
        for key in ("hamming_flagged", "affinity_flagged", "agree", "affinity_only",
                    "tiers_disagree"):
            assert key in out

    def test_disagreement_is_reported_without_moving_the_gate(self, scores):
        """The case that motivates the whole two-tier split.

        AAAAUGUGG and GGGGGGGGG differ in 6 of 9 positions, so tier A accepts them. But
        the PPR designed for the first scores the second at 0.82 of its own target,
        because TN and ND both score G positively in the table. The model would have you
        look; Hamming would not. Both facts are reported and the gate does not move.
        """
        out = compare_tiers(["AAAAUGUGG", "CCCCCCCCC", "GGGGGGGGG"], scores, min_hamming=5)
        assert out["hamming_flagged"] == [], "tier A accepts this library"
        assert out["affinity_only"] == [("AAAAUGUGG", "GGGGGGGGG")]
        assert out["tiers_disagree"] is True

    def test_a_crowded_library_is_caught_by_tier_a(self, scores):
        out = compare_tiers(["AAAAUGUGG", "AAAAUGUGA", "AAAAUGUAA"], scores, min_hamming=5)
        assert out["hamming_flagged"], "tier A must catch these on its own"
