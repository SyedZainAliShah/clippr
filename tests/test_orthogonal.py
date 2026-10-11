"""Unit tests for orthogonal target-set design."""
from __future__ import annotations

import pytest

from ppropt.orthogonal import (
    crosstalk_report,
    design_orthogonal_set,
    distance,
)
from ppropt.ppr import describe


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


class TestWeightedMetricDomain:
    """The optional weighted metric covers positions 1-10 and nothing else.

    `POSITION_WEIGHTS.get(i + 1, 0.0)` gave every later position weight zero, so two
    14-mers differing only after position 10 scored weighted distance 0 against uniform
    distance 4, and the 19-mer pair 0 against 9. A reassuring zero for targets nobody
    measured is worse than a refusal.
    """

    def test_fourteen_mers_differing_only_past_position_ten(self):
        a, b = "A" * 14, "A" * 10 + "U" * 4
        from ppropt.orthogonal import _weights, distance
        assert distance(a, b, _weights(14, "uniform")) == 4
        with pytest.raises(ValueError, match="positions 1-10"):
            _weights(14, "weighted")

    def test_nineteen_mers_differing_only_past_position_ten(self):
        a, b = "A" * 19, "A" * 10 + "U" * 9
        from ppropt.orthogonal import _weights, distance
        assert distance(a, b, _weights(19, "uniform")) == 9
        with pytest.raises(ValueError, match="positions 1-10"):
            _weights(19, "weighted")

    def test_a_refusal_is_not_a_silent_downgrade_to_uniform(self):
        from ppropt.orthogonal import _weights
        with pytest.raises(ValueError):
            _weights(14, "weighted")

    def test_nine_positions_still_work(self):
        from ppropt.orthogonal import _weights
        w = _weights(9, "weighted")
        assert len(w) == 9 and sum(w) == pytest.approx(9.0)

    def test_ten_positions_are_the_boundary_and_are_allowed(self):
        from ppropt.orthogonal import _weights, WEIGHTED_METRIC_MAX_LENGTH
        assert WEIGHTED_METRIC_MAX_LENGTH == 10
        assert len(_weights(10, "weighted")) == 10

    def test_a_short_weight_vector_is_refused_rather_than_zipped(self):
        """`zip` stops at the shortest argument, so this compared a prefix silently."""
        from ppropt.orthogonal import distance
        with pytest.raises(ValueError, match="weight vector"):
            distance("AAAAAAAAA", "AAAAAAAAU", [1.0] * 4)

    def test_negative_and_non_finite_weights_are_refused(self):
        from ppropt.orthogonal import distance
        for bad in ([1.0, -1.0, 1.0], [1.0, float("inf"), 1.0], [1.0, float("nan"), 1.0]):
            with pytest.raises(ValueError, match="finite and non-negative"):
                distance("AAA", "AUA", bad)
