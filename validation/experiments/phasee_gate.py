"""Phase E gate — pooling against exhaustive optima, and pricing against hand arithmetic.

Gate E, from the full-workflow plan: boundaries just below, at and above every limit; tier
transitions; empty input; ineligible counts; modifications; invalid bases; quantity-preserving
grouping; CSV round trips; and price arithmetic against independently hand-calculated
fixtures. Most of those live in `tests/test_ordering.py`, which is where they belong.

What needs a script is the part a unit test cannot express: **how far the pooling heuristic is
from optimal**. The plan forbids calling a heuristic globally cost-optimal, so this enumerates
every valid partition on small cases and reports the gap instead of asserting there is none.

No order is placed, and nothing here contacts a vendor.

    python validation/experiments/phasee_gate.py
"""
from __future__ import annotations

import argparse
import json
import sys
from decimal import Decimal
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from clippr.ordering import (HISTORICAL_OPOOL, PriceTier, ProductProfile,  # noqa: E402
                             check_eligibility, estimate_cost, plan_pools)


def partitions(items: list[int], min_size: int, max_size: int):
    """Every way to split `items` into groups within [min_size, max_size].

    Brute force and obviously correct, because it is the reference the heuristic is measured
    against. Only usable on small inputs, which is exactly what it is for.
    """
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for size in range(min_size - 1, min(max_size, len(items))):
        for companions in combinations(rest, size):
            group = [first, *companions]
            remainder = [x for x in rest if x not in companions]
            for tail in partitions(remainder, min_size, max_size):
                yield [group, *tail]


def pool_cost(groups, profile) -> Decimal | None:
    """Total cost of a partition under the profile, or None if any group is unpriceable."""
    total = Decimal("0")
    for group in groups:
        tier = next((t for t in profile.price_tiers
                     if t.min_oligos <= len(group)
                     and (t.max_oligos is None or len(group) <= t.max_oligos)), None)
        if tier is None:
            return None
        total += tier.cost(len(group), 0)
    return total


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "work" / "phasee"))
    args = ap.parse_args()

    profile = ProductProfile(
        name="synthetic", vendor="V", region="EU", currency="EUR",
        min_oligo_nt=10, max_oligo_nt=100, min_oligos=2, max_oligos=4,
        price_tiers=(PriceTier(1, 4, Decimal("100.00")),),
        priced_on="2026-09-15", price_status="current")

    print("E1  pooling heuristic against exhaustive optima")
    rows, worst_gap = [], Decimal("0")
    for n in range(2, 11):
        items = [{"sequence": "ACGT" * 5, "design": f"d{i}"} for i in range(n)]
        plan = plan_pools(items, profile)

        best = None
        for groups in partitions(list(range(n)), profile.min_oligos, profile.max_oligos):
            cost = pool_cost(groups, profile)
            if cost is not None and (best is None or cost < best):
                best = cost

        if not plan["feasible"]:
            rows.append({"items": n, "heuristic": None, "optimal": str(best) if best else None,
                         "reason": plan["reason"]})
            print(f"      n={n:>2}  heuristic infeasible; exhaustive optimum "
                  f"{best if best is not None else 'none'}")
            continue

        got = pool_cost([[m["design"] for m in p["members"]] for p in plan["pools"]], profile)
        gap = (got - best) if best is not None else Decimal("0")
        worst_gap = max(worst_gap, gap)
        rows.append({"items": n, "pools": plan["pool_count"], "heuristic": str(got),
                     "optimal": str(best), "gap": str(gap)})
        print(f"      n={n:>2}  {plan['pool_count']} pool(s)  heuristic {got}  "
              f"optimum {best}  gap {gap}")

    # The heuristic must never be infeasible where a valid partition exists: that would be
    # the heuristic manufacturing the refusal, which is what plain first-fit did at n=9.
    manufactured = [r for r in rows if r["heuristic"] is None and r["optimal"] is not None]
    print(f"    cases refused by the heuristic where a valid partition exists: "
          f"{len(manufactured)}")
    print(f"    worst cost gap from optimal: {worst_gap}")

    print("\nE2  price arithmetic against hand-calculated fixtures")
    fixtures = [
        # (oligos, modifications, tax, shipping, hand-worked total)
        (3, (), None, None, "100.00"),
        (3, ("phos",), None, None, "104.89"),          # 100.00 + 1.63*3
        (3, (), Decimal("0.19"), None, "119.00"),      # 100.00 * 1.19
        (3, (), None, Decimal("12.50"), "112.50"),
        (3, (), Decimal("0.19"), Decimal("12.50"), "133.88"),   # 112.50 * 1.19
    ]
    priced = ProductProfile(**{**profile.__dict__,
                               "modifications": {"phos": {"per_oligo": "1.63"}}})
    arithmetic_ok = True
    for count, mods, tax, ship, expected in fixtures:
        got = estimate_cost(["ACGT" * 5] * count, priced, modifications=mods,
                            tax_rate=tax, shipping=ship)
        agrees = got.get("total") == expected
        arithmetic_ok &= agrees
        print(f"      {count} oligos mods={mods or '-'} tax={tax or '-'} ship={ship or '-'}"
              f"  expected {expected}  got {got.get('total')}  "
              f"{'ok' if agrees else 'MISMATCH'}")

    print("\nE3  price provenance")
    historical = estimate_cost(["ACGT" * 5] * 3, HISTORICAL_OPOOL)
    refuses = not historical["available"] and historical["price_status"] == "historical"
    still_checks = check_eligibility(["ACGT" * 5] * 3, HISTORICAL_OPOOL)["eligible"]
    print(f"      historical constants refuse to price: {refuses}")
    print(f"      eligibility still works without a price: {still_checks}")
    print(f"      unresolved rules declared: {len(HISTORICAL_OPOOL.unresolved_rules)}")

    passed = not manufactured and arithmetic_ok and refuses and still_checks
    print(f"\nGATE E: {'PASS' if passed else 'FAIL'}")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "gate_e.json").write_text(json.dumps({
        "pooling": rows,
        "worst_cost_gap_from_optimal": str(worst_gap),
        "heuristic_manufactured_refusals": manufactured,
        "price_arithmetic_matches_hand_fixtures": arithmetic_ok,
        "historical_constants_refuse_to_price": refuses,
        "eligibility_works_without_price": still_checks,
        "no_order_placed": True,
        "scope": ("local implementation checks only; current vendor documentation remains "
                  "the authority and no order was placed or submitted"),
        "gate": "PASS" if passed else "FAIL",
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out / 'gate_e.json'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
