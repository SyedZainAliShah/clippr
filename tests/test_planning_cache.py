"""The planning cache must be invisible except in how long planning takes.

Everything here uses hand-built candidate pools rather than real proteins: the property
under test is the cache key, and exercising it through a 19S planner would turn a unit test
into a five-minute search. Equivalence on real corpus proteins is
`validation/check_planning_cache.py`.
"""
from __future__ import annotations

import pytest

from clippr import overhangs
from clippr.overhangs import best_set, cache_stats, clear_cache, matrix_fingerprint

POOLS = [["AGCG", "ATTC"], ["AATG", "GGTA"], ["CGCT", "ACTC"]]


@pytest.fixture(autouse=True)
def _clean_cache():
    clear_cache()
    yield
    clear_cache()


class TestTransparency:
    def test_a_hit_returns_the_same_answer_as_a_miss(self):
        cold = best_set(POOLS)
        warm = best_set(POOLS)
        assert cold == warm
        assert cache_stats()["hits"] == 1

    def test_a_hit_returns_a_fresh_list(self):
        """A shared list would let one caller's sort corrupt every later hit."""
        first, _ = best_set(POOLS)
        first.append("MUTATED")
        second, _ = best_set(POOLS)
        assert "MUTATED" not in second

    def test_clearing_resets_hits_and_entries(self):
        best_set(POOLS)
        clear_cache()
        assert cache_stats() == {"hits": 0, "misses": 0, "entries": 0}


class TestKeyCoverage:
    """Each input that can change the answer must change the key."""

    def test_the_matrix_is_part_of_the_key(self):
        best_set(POOLS, matrix="BsaI-HFv2")
        best_set(POOLS, matrix="BbsI-HF")
        assert cache_stats() == {"hits": 0, "misses": 2, "entries": 2}

    def test_the_fixed_overhangs_are_part_of_the_key(self):
        best_set(POOLS, fixed=["CTCA", "CTCG"])
        best_set(POOLS, fixed=["GGAG", "CGCT"])
        assert cache_stats()["misses"] == 2

    def test_the_candidate_pools_are_part_of_the_key(self):
        best_set(POOLS)
        best_set([POOLS[0], POOLS[1], ["CGCT", "TTCG"]])
        assert cache_stats()["misses"] == 2

    def test_pool_order_is_part_of_the_key(self):
        """Junction order decides which overhang is locked where, so it is not a set."""
        best_set(POOLS)
        best_set(list(reversed(POOLS)))
        assert cache_stats()["misses"] == 2

    def test_a_scorer_change_invalidates_the_cache(self, monkeypatch):
        """The failure this guards: a stale entry outliving a change to the formula.

        Version 1 chose different junctions from version 2 on 5 of 5 targets tested, so an
        entry that survived that change would resurrect the superseded selection.
        """
        best_set(POOLS)
        monkeypatch.setattr(overhangs, "SCORER_VERSION", 99)
        best_set(POOLS)
        assert cache_stats() == {"hits": 0, "misses": 2, "entries": 2}

    def test_the_search_budget_is_part_of_the_key(self):
        """Conservative on purpose: a raised limit can only be searched, never assumed.

        Both limits happen to admit this pool, so the answers agree -- but a key that
        ignored the budget would also serve a truncated search's answer to a raised one.
        """
        a = best_set(POOLS, max_combinations=500_000)
        b = best_set(POOLS, max_combinations=10)
        assert cache_stats() == {"hits": 0, "misses": 2, "entries": 2}
        assert a == b

    def test_the_table_contents_are_part_of_the_key(self):
        """A fingerprint that never moves is not protecting anything."""
        assert matrix_fingerprint("BsaI-HFv2") != matrix_fingerprint("BbsI-HF")


class TestFailuresAreNotCached:
    def test_an_impossible_pool_still_raises_on_the_second_call(self):
        """Caching only successes keeps the retry path intact for the caller above."""
        impossible = [["AATT"], ["AATT"]]          # palindromic, so no valid set exists
        with pytest.raises(ValueError):
            best_set(impossible)
        with pytest.raises(ValueError):
            best_set(impossible)
        assert cache_stats()["entries"] == 0
