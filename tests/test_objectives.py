"""The objective definitions, checked against hand-worked values.

Expected values here are computed by hand or from the declared formula, never by calling the
function under test and recording what it said.
"""
from __future__ import annotations

import math

import pytest

from clippr.objectives import (codon_adaptation, collection_penalty, distinct_repeated_words,
                               dominates, duplicated_kmers, nondominated,
                               relative_adaptiveness, worst_shared_tract)

#: Two amino acids, two codons each, deliberately lopsided so weights are easy to verify.
TOY = {"K": {"AAA": 0.8, "AAG": 0.2}, "F": {"TTT": 0.5, "TTC": 0.5},
       "M": {"ATG": 1.0}, "W": {"TGG": 1.0}}


class TestRelativeAdaptiveness:
    def test_weights_are_relative_to_the_best_synonym(self):
        w = relative_adaptiveness(TOY)
        assert w["AAA"] == 1.0
        assert w["AAG"] == pytest.approx(0.25)
        assert w["TTT"] == w["TTC"] == 1.0


class TestCodonAdaptation:
    def test_all_optimal_codons_give_one(self):
        assert codon_adaptation("AAATTT", TOY)["cai"] == pytest.approx(1.0)

    def test_the_geometric_mean_is_the_declared_formula(self):
        # one AAG (0.25) and one TTT (1.0): exp(mean(log w)) = sqrt(0.25)
        got = codon_adaptation("AAGTTT", TOY)
        assert got["cai"] == pytest.approx(math.sqrt(0.25))
        assert got["codons"] == 2

    def test_single_codon_amino_acids_are_excluded(self):
        """M and W have weight 1 by construction; including them would dilute the score."""
        with_them = codon_adaptation("AAGATGTGG", TOY)
        without = codon_adaptation("AAG", TOY)
        assert with_them["cai"] == pytest.approx(without["cai"])
        assert with_them["codons"] == 1

    def test_a_codon_absent_from_the_table_is_counted_not_assumed(self):
        got = codon_adaptation("AAGCCC", TOY)
        assert got["missing_weights"] == 1
        assert got["codons"] == 1

    def test_a_zero_weight_gives_zero_not_negative_infinity(self):
        table = {"K": {"AAA": 1.0, "AAG": 0.0}}
        got = codon_adaptation("AAG", table)
        assert got["cai"] == 0.0
        assert got["zero_weight_codons"] == ["AAG"]

    def test_the_scored_interval_is_reported(self):
        got = codon_adaptation("AAATTTAAA", TOY)
        assert got["nt"] == 9 and got["codons"] == 3


class TestRepetition:
    def test_duplicated_counts_occurrences_beyond_the_first(self):
        assert duplicated_kmers("AAAA", k=2) == 2          # 3 windows, 1 distinct

    def test_distinct_words_counts_words_not_occurrences(self):
        assert distinct_repeated_words("AAAA", k=2) == 1

    def test_they_differ_and_are_not_interchangeable(self):
        seq = "ATATATAT"
        assert duplicated_kmers(seq, k=2) != distinct_repeated_words(seq, k=2)

    def test_a_sequence_shorter_than_k_has_no_repeats(self):
        assert duplicated_kmers("AT", k=20) == 0


class TestCollectionSharing:
    def test_identical_records_share_every_window(self):
        got = collection_penalty({"a": "ACGTACGT", "b": "ACGTACGT"}, k=4)
        # 5 windows each, all shared
        assert got["penalty"] == 5
        assert got["k"] == 4

    def test_disjoint_records_share_nothing(self):
        assert collection_penalty({"a": "AAAAAAAA", "b": "CCCCCCCC"}, k=4)["penalty"] == 0

    def test_a_single_record_has_no_pairs(self):
        assert collection_penalty({"a": "ACGTACGT"}, k=4)["penalty"] == 0

    def test_worst_tract_names_the_pair(self):
        got = worst_shared_tract({"a": "TTTTACGTACGTTTTT", "b": "GGGACGTACGTGGG"})
        assert got["nt"] == len("ACGTACGT")
        assert set(got["pair"]) == {"a", "b"}


class TestDominance:
    DIRECTIONS = {"gain": "max", "cost": "min"}

    def test_better_on_both_dominates(self):
        assert dominates({"gain": 2, "cost": 1}, {"gain": 1, "cost": 2}, self.DIRECTIONS)

    def test_equal_does_not_dominate(self):
        assert not dominates({"gain": 1, "cost": 1}, {"gain": 1, "cost": 1},
                             self.DIRECTIONS)

    def test_a_trade_off_dominates_neither_way(self):
        a, b = {"gain": 2, "cost": 2}, {"gain": 1, "cost": 1}
        assert not dominates(a, b, self.DIRECTIONS)
        assert not dominates(b, a, self.DIRECTIONS)

    def test_a_difference_within_tolerance_is_not_dominance(self):
        """Otherwise floating-point noise manufactures a front."""
        assert not dominates({"gain": 1 + 1e-12, "cost": 1}, {"gain": 1, "cost": 1},
                             self.DIRECTIONS)

    def test_an_unknown_direction_is_refused(self):
        with pytest.raises(ValueError, match="expected 'max' or 'min'"):
            dominates({"gain": 1}, {"gain": 0}, {"gain": "sideways"})


class TestNondominated:
    DIRECTIONS = {"gain": "max", "cost": "min"}

    def test_a_dominated_point_is_excluded(self):
        candidates = [{"gain": 3, "cost": 1}, {"gain": 1, "cost": 5}]
        assert nondominated(candidates, self.DIRECTIONS) == [0]

    def test_a_trade_off_keeps_both(self):
        candidates = [{"gain": 3, "cost": 3}, {"gain": 1, "cost": 1}]
        assert nondominated(candidates, self.DIRECTIONS) == [0, 1]

    def test_exact_duplicates_are_both_kept(self):
        """Neither strictly beats the other, so neither may be dropped silently."""
        candidates = [{"gain": 2, "cost": 2}, {"gain": 2, "cost": 2}]
        assert nondominated(candidates, self.DIRECTIONS) == [0, 1]

    def test_an_empty_candidate_set_gives_an_empty_front(self):
        assert nondominated([], self.DIRECTIONS) == []

    def test_a_single_candidate_is_its_own_front(self):
        assert nondominated([{"gain": 1, "cost": 1}], self.DIRECTIONS) == [0]
