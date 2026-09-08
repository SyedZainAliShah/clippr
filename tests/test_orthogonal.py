"""Unit tests for orthogonal target-set design."""
from __future__ import annotations

import pytest

from clippr.orthogonal import (
    crosstalk_report,
    design_orthogonal_set,
    distance,
)
from clippr.ppr import describe


class TestDistance:
    def test_identical_is_zero(self):
        assert distance("ACGUACGUA", "ACGUACGUA", [1.0] * 9) == 0

    def test_counts_mismatches(self):
        assert distance("AAAAAAAAA", "AAAAAAAAC", [1.0] * 9) == 1
        assert distance("AAAAAAAAA", "CCCCCCCCC", [1.0] * 9) == 9

    def test_respects_weights(self):
        w = [0.0] * 8 + [5.0]
        assert distance("AAAAAAAAA", "CAAAAAAAA", w) == 0      # unweighted position
        assert distance("AAAAAAAAA", "AAAAAAAAC", w) == 5      # weighted position

    def test_length_mismatch_rejected(self):
        with pytest.raises(ValueError):
            distance("ACGU", "ACGUA", [1.0] * 4)


class TestDesign:
    def test_returns_requested_count(self):
        ts = design_orthogonal_set(6, length=9, seed=0, iterations=500)
        assert len(ts.targets) == 6

    def test_targets_are_distinct(self):
        ts = design_orthogonal_set(8, length=9, seed=0, iterations=500)
        assert len(set(ts.targets)) == 8

    def test_targets_are_valid_rna(self):
        ts = design_orthogonal_set(5, length=9, seed=0, iterations=500)
        for t in ts.targets:
            assert set(t) <= set("ACGU")
            assert len(t) == 9

    def test_every_target_designs_a_binder(self):
        ts = design_orthogonal_set(5, length=9, seed=0, iterations=500)
        for t in ts.targets:
            d = describe(t)
            assert d["n_arelf"] == 9
            assert d["aa_sequence"].startswith("M")

    def test_beats_random_separation(self):
        """Optimised sets should separate better than an arbitrary draw."""
        import random
        ts = design_orthogonal_set(8, length=9, seed=1, iterations=4000)
        rng = random.Random(1)
        rand = ["".join(rng.choice("ACGU") for _ in range(9)) for _ in range(8)]
        w = [1.0] * 9
        rand_min = min(distance(a, b, w)
                       for i, a in enumerate(rand) for b in rand[i + 1:])
        assert ts.min_distance > rand_min

    def test_deterministic_under_seed(self):
        a = design_orthogonal_set(6, length=9, seed=42, iterations=800)
        b = design_orthogonal_set(6, length=9, seed=42, iterations=800)
        assert a.targets == b.targets

    def test_forbid_is_respected(self):
        ts = design_orthogonal_set(4, length=9, seed=0, iterations=300,
                                   forbid=["AAAAAAAAA", "CCCCCCCCC"])
        assert "AAAAAAAAA" not in ts.targets
        assert "CCCCCCCCC" not in ts.targets

    @pytest.mark.parametrize("metric", ["uniform", "weighted"])
    def test_both_metrics_work(self, metric):
        ts = design_orthogonal_set(6, length=9, metric=metric, seed=0, iterations=500)
        assert ts.min_distance > 0
        assert ts.metric == metric

    def test_unknown_metric_rejected(self):
        with pytest.raises(ValueError, match="unknown metric"):
            design_orthogonal_set(4, length=9, metric="nonsense", iterations=10)

    def test_impossible_request_rejected(self):
        with pytest.raises(ValueError):
            design_orthogonal_set(1, length=9)
        with pytest.raises(ValueError):
            design_orthogonal_set(5, length=1)      # only 4 distinct 1-mers


class TestReport:
    def test_report_mentions_every_target(self):
        ts = design_orthogonal_set(5, length=9, seed=0, iterations=300)
        text = ts.report()
        for t in ts.targets:
            assert t in text

    def test_crosstalk_matrix_is_symmetric(self):
        ts = design_orthogonal_set(5, length=9, seed=0, iterations=300)
        for i in range(5):
            for j in range(5):
                assert ts.distances[i][j] == ts.distances[j][i]

    def test_crosstalk_report_renders(self):
        ts = design_orthogonal_set(4, length=9, seed=0, iterations=300)
        assert crosstalk_report(ts.targets).count("[") >= 4
