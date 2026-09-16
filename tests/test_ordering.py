"""Order-product eligibility, pricing and pooling — Gate E's boundary and arithmetic cases.

Expected costs here are hand-calculated from the declared profile, never taken from a run of
the function under test. A pricing test that records what the code produced is a change
detector, not a check.
"""
from __future__ import annotations

import csv
import json
from decimal import Decimal
from pathlib import Path

from clippr.ordering import (HISTORICAL_OPOOL, PriceTier, ProductProfile, check_eligibility,
                             estimate_cost, load_profile, plan_pools, write_order_files)


def profile(**kw) -> ProductProfile:
    base = dict(name="test", vendor="V", region="EU", currency="EUR",
                min_oligo_nt=20, max_oligo_nt=60, min_oligos=2, max_oligos=4,
                scales=("50 pmol",),
                price_tiers=(PriceTier(1, 4, Decimal("100.00")),
                             PriceTier(5, 10, Decimal("150.00")),
                             PriceTier(11, None, Decimal("10.00"), per="oligo")),
                priced_on="2026-09-15", price_status="current",
                modifications={"phos": {"per_oligo": "1.63"}})
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


class TestPricingArithmetic:
    """Every expected total below is worked out by hand from the declared tiers."""

    def test_a_pool_tier_does_not_scale_with_count(self):
        got = estimate_cost([seq(30)] * 3, profile())
        assert got["total"] == "100.00"

    def test_the_tier_transition_is_at_the_declared_boundary(self):
        # tier 1 covers 1-4 at 100.00; tier 2 covers 5-10 at 150.00
        four = estimate_cost([seq(30)] * 4, profile(max_oligos=20))
        five = estimate_cost([seq(30)] * 5, profile(max_oligos=20))
        assert four["total"] == "100.00"
        assert five["total"] == "150.00"

    def test_a_per_oligo_tier_multiplies(self):
        # tier 3: 11+ oligos at 10.00 each -> 12 * 10.00 = 120.00
        got = estimate_cost([seq(30)] * 12, profile(max_oligos=20))
        assert got["total"] == "120.00"

    def test_modifications_are_priced_per_oligo(self):
        # 100.00 + 1.63 * 3 = 104.89
        got = estimate_cost([seq(30)] * 3, profile(), modifications=("phos",))
        assert got["total"] == "104.89"

    def test_tax_and_shipping_compose_in_the_declared_order(self):
        # 100.00 + 12.50 shipping = 112.50; 19% tax = 21.375 -> 21.38; total 133.88
        got = estimate_cost([seq(30)] * 3, profile(),
                            tax_rate=Decimal("0.19"), shipping=Decimal("12.50"))
        assert got["total"] == "133.88"
        assert got["tax_included"] and got["shipping_included"]

    def test_tax_and_shipping_are_never_assumed(self):
        got = estimate_cost([seq(30)] * 3, profile())
        assert not got["tax_included"] and not got["shipping_included"]

    def test_rounding_is_half_up_on_the_cent(self):
        # 100.00 * 0.125 = 12.50 exactly; use a rate that lands on a half cent
        got = estimate_cost([seq(30)] * 3, profile(), tax_rate=Decimal("0.0725"))
        # 100.00 * 0.0725 = 7.25 -> total 107.25
        assert got["total"] == "107.25"

    def test_an_unknown_modification_is_refused_not_ignored(self):
        got = estimate_cost([seq(30)] * 3, profile(), modifications=("gold plating",))
        assert not got["available"] and "no rule for modification" in got["reason"]

    def test_a_count_outside_every_tier_is_refused(self):
        narrow = profile(price_tiers=(PriceTier(5, 6, Decimal("100.00")),), max_oligos=20)
        got = estimate_cost([seq(30)] * 2, narrow)
        assert not got["available"] and "no price tier" in got["reason"]

    def test_an_estimate_never_calls_itself_a_quote(self):
        assert estimate_cost([seq(30)] * 3, profile())["is_a_quote"] is False


class TestPriceProvenance:
    def test_the_historical_profile_refuses_to_price(self):
        """Historical constants must never surface as a current number."""
        got = estimate_cost([seq(30)] * 3, HISTORICAL_OPOOL)
        assert not got["available"]
        assert got["price_status"] == "historical"

    def test_the_historical_profile_still_checks_eligibility(self):
        """The useful half must not depend on a price that has gone stale."""
        assert check_eligibility([seq(30)] * 3, HISTORICAL_OPOOL)["eligible"]

    def test_an_undated_profile_is_not_usable_for_pricing(self):
        assert not profile(priced_on=None).prices_usable

    def test_a_dated_estimate_carries_its_date_and_source(self):
        got = estimate_cost([seq(30)] * 3, profile(source_url="https://example.invalid"))
        assert got["priced_on"] == "2026-09-15"
        assert got["source_url"] == "https://example.invalid"


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
        report = json.loads(Path(paths["estimate_report"]).read_text(encoding="utf-8"))
        assert report["profile"]["limits"]["oligos"] == [2, 4]
        assert report["profile"]["price_status"] == "current"


class TestProfileFile:
    def test_a_supplied_dated_profile_loads_and_prices(self, tmp_path):
        path = tmp_path / "p.json"
        path.write_text(json.dumps({
            "name": "supplied", "vendor": "V", "region": "EU", "currency": "EUR",
            "min_oligo_nt": 20, "max_oligo_nt": 60, "min_oligos": 2, "max_oligos": 8,
            "price_tiers": [{"min_oligos": 1, "max_oligos": None, "price": "77.77"}],
            "priced_on": "2026-01-01", "price_status": "current",
            "source_url": "https://example.invalid",
        }), encoding="utf-8")
        loaded = load_profile(path)
        assert loaded.prices_usable
        assert estimate_cost([seq(30)] * 3, loaded)["total"] == "77.77"
