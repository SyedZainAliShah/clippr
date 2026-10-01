"""Order-product eligibility, pricing and pooling — Gate E's boundary and arithmetic cases.

Expected costs here are hand-calculated from the declared profile, never taken from a run of
the function under test. A pricing test that records what the code produced is a change
detector, not a check.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from clippr.ordering import (HISTORICAL_OPOOL, ProductProfile, check_eligibility,
                             load_profile, plan_pools, write_order_files)


def profile(**kw) -> ProductProfile:
    base = dict(name="test", vendor="V", region="EU",
                min_oligo_nt=20, max_oligo_nt=60, min_oligos=2, max_oligos=4,
                scales=("50 pmol",))
    base.update(kw)
    return ProductProfile(**base)


def seq(n: int) -> str:
    return ("ACGT" * (n // 4 + 1))[:n]


class TestLengthBoundaries:
    """Just below, at, and just above each limit."""

    def test_one_below_the_minimum_length_is_a_violation(self):
        got = check_eligibility([seq(19), seq(30)], profile())
        assert not got["eligible"]
        assert any("below the minimum" in v for v in got["violations"])

    def test_exactly_the_minimum_length_is_allowed(self):
        assert check_eligibility([seq(20), seq(30)], profile())["eligible"]

    def test_exactly_the_maximum_length_is_allowed(self):
        assert check_eligibility([seq(60), seq(30)], profile())["eligible"]

    def test_one_above_the_maximum_length_is_a_violation(self):
        got = check_eligibility([seq(61), seq(30)], profile())
        assert not got["eligible"]
        assert any("above the maximum" in v for v in got["violations"])


class TestCountBoundaries:
    """Count limits are **per pool**. An order larger than one pool gets split, not refused."""

    def test_too_few_for_even_one_pool_is_a_violation(self):
        got = check_eligibility([seq(30)], profile())
        assert not got["eligible"]
        assert any("cannot fill even one pool" in v for v in got["violations"])

    def test_the_minimum_count_message_says_no_filler_is_added(self):
        """Padding an order to reach a threshold would ship something unasked for."""
        got = check_eligibility([seq(30)], profile())
        assert any("no filler" in v for v in got["violations"])

    def test_exactly_the_minimum_count_is_allowed(self):
        assert check_eligibility([seq(30)] * 2, profile())["eligible"]

    def test_exactly_the_maximum_count_is_allowed(self):
        assert check_eligibility([seq(30)] * 4, profile())["eligible"]

    def test_more_than_one_pool_is_not_an_order_wide_violation(self):
        """Applying the per-pool maximum to the whole order refused every split order."""
        assert check_eligibility([seq(30)] * 9, profile())["eligible"]

    def test_the_per_pool_check_still_enforces_the_maximum(self):
        got = check_eligibility([seq(30)] * 5, profile(), per_pool=True)
        assert not got["eligible"]
        assert any("per pool" in v for v in got["violations"])

    def test_the_per_pool_check_still_enforces_the_minimum(self):
        got = check_eligibility([seq(30)], profile(), per_pool=True)
        assert not got["eligible"]
        assert any("per pool" in v for v in got["violations"])


class TestOtherEligibility:
    def test_empty_input_is_refused_not_crashed(self):
        got = check_eligibility([], profile())
        assert not got["eligible"] and got["oligos"] == 0

    def test_invalid_bases_are_named(self):
        got = check_eligibility([seq(28) + "NNXX", seq(30)], profile())
        assert any("outside" in v for v in got["violations"])

    def test_duplicates_warn_but_do_not_block(self):
        got = check_eligibility([seq(30), seq(30)], profile())
        assert got["eligible"]
        assert any("duplicate" in w for w in got["warnings"])

    def test_violations_warnings_and_unresolved_are_separate(self):
        got = check_eligibility([seq(30)] * 2, profile(unresolved_rules=("vendor terms",)))
        assert got["violations"] == [] and got["unresolved"] == ["vendor terms"]

    def test_scales_are_offered_only_when_eligible(self):
        assert check_eligibility([seq(30)] * 2, profile())["eligible_scales"]
        assert check_eligibility([seq(30)], profile())["eligible_scales"] == []


class TestHistoricalProfile:
    def test_the_historical_profile_still_checks_eligibility(self):
        """A stale rule set must still answer the eligibility question."""
        assert check_eligibility([seq(30)] * 3, HISTORICAL_OPOOL)["eligible"]


class TestPooling:
    def _items(self, n, length=30):
        return [{"sequence": seq(length), "design": f"d{i}", "module": f"m{i}",
                 "quantity": 1} for i in range(n)]

    def test_items_are_split_to_respect_the_maximum(self):
        got = plan_pools(self._items(9), profile())     # max 4 per pool
        assert got["feasible"] and got["pool_count"] == 3

    def test_every_item_lands_in_exactly_one_pool(self):
        got = plan_pools(self._items(9), profile())
        placed = [m for p in got["pools"] for m in p["members"]]
        assert len(placed) == 9
        assert len({m["design"] for m in placed}) == 9

    def test_source_mappings_and_quantities_survive(self):
        """A pool plan that loses which design an oligo came from cannot be ordered against."""
        got = plan_pools(self._items(6), profile())
        for pool in got["pools"]:
            for member in pool["members"]:
                assert member["design"] and member["module"] and member["quantity"] == 1

    def test_duplicate_sequences_are_not_merged_away(self):
        items = [{"sequence": seq(30), "design": "a", "quantity": 1},
                 {"sequence": seq(30), "design": "b", "quantity": 1}]
        got = plan_pools(items, profile())
        placed = [m for p in got["pools"] for m in p["members"]]
        assert {m["design"] for m in placed} == {"a", "b"}

    def test_unused_capacity_is_reported(self):
        got = plan_pools(self._items(5), profile())
        assert sum(p["unused_capacity"] for p in got["pools"]) == 3    # 8 slots, 5 used

    def test_an_oversized_oligo_makes_the_plan_infeasible(self):
        items = self._items(3) + [{"sequence": seq(999), "design": "x"}]
        got = plan_pools(items, profile())
        assert not got["feasible"] and "exceed" in got["reason"]

    def test_a_pool_below_the_minimum_is_refused_without_filler(self):
        got = plan_pools(self._items(1), profile())
        assert not got["feasible"] and "no filler" in got["reason"]

    def test_empty_input_is_refused(self):
        assert not plan_pools([], profile())["feasible"]

    def test_it_does_not_claim_to_be_optimal(self):
        got = plan_pools(self._items(6), profile())
        assert "not a proven cost optimum" in got["optimality"]

    def test_it_is_deterministic(self):
        a = plan_pools(self._items(9), profile())
        b = plan_pools(self._items(9), profile())
        assert a["pools"] == b["pools"]


class TestOrderFiles:
    def test_the_files_round_trip(self, tmp_path):
        items = [{"sequence": seq(30), "design": f"d{i}", "module": f"m{i}", "quantity": 2}
                 for i in range(4)]
        plan = plan_pools(items, profile())
        paths = write_order_files(plan, profile(), tmp_path)

        with open(paths["order_sequences"], encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 4 and set(rows[0]) == {"Pool name", "Sequence"}

        with open(paths["pool_membership"], encoding="utf-8", newline="") as handle:
            members = list(csv.DictReader(handle))
        assert len(members) == 4
        assert {m["design"] for m in members} == {"d0", "d1", "d2", "d3"}
        assert all(int(m["quantity"]) == 2 for m in members)
        assert all(int(m["length"]) == len(m["sequence"]) for m in members)

    def test_the_report_carries_the_profile_and_its_limits(self, tmp_path):
        plan = plan_pools([{"sequence": seq(30), "design": "a"},
                           {"sequence": seq(30), "design": "b"}], profile())
        paths = write_order_files(plan, profile(), tmp_path)
        report = json.loads(Path(paths["order_plan"]).read_text(encoding="utf-8"))
        assert report["profile"]["limits"]["oligos"] == [2, 4]


class TestProfileFile:
    def test_a_supplied_profile_loads(self, tmp_path):
        path = tmp_path / "p.json"
        path.write_text(json.dumps({
            "name": "supplied", "vendor": "V", "region": "EU",
            "min_oligo_nt": 20, "max_oligo_nt": 60, "min_oligos": 2, "max_oligos": 8,
            "source_url": "https://example.invalid",
        }), encoding="utf-8")
        loaded = load_profile(path)
        assert loaded.name == "supplied"
        assert loaded.max_oligos == 8
        assert check_eligibility([seq(30)] * 3, loaded)["eligible"]
