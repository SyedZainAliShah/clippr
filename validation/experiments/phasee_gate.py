"""Phase E gate — pooling against exhaustive optima.

Gate E, from the full-workflow plan: boundaries just below, at and above every limit; tier
transitions; empty input; ineligible counts; invalid bases; quantity-preserving
grouping; and CSV round trips. Most of those live in `tests/test_ordering.py`, which is where they belong.

What needs a script is the part a unit test cannot express: **how far the pooling heuristic is
from optimal**. The plan forbids calling a heuristic globally optimal, so this enumerates
every valid partition on small cases and reports the gap instead of asserting there is none.

No order is placed, and nothing here contacts a vendor.

    python validation/experiments/phasee_gate.py
"""
from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from clippr.ordering import (HISTORICAL_OPOOL, ProductProfile,  # noqa: E402
                             check_eligibility, plan_pools)


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



def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "work" / "phasee"))
    args = ap.parse_args()

    profile = ProductProfile(
        name="synthetic", vendor="V", region="EU",
        min_oligo_nt=10, max_oligo_nt=100, min_oligos=2, max_oligos=4)

    print("E1  pooling heuristic against exhaustive optima")
    rows, worst_gap = [], 0
    for n in range(2, 11):
        items = [{"sequence": "ACGT" * 5, "design": f"d{i}"} for i in range(n)]
        plan = plan_pools(items, profile)

        # A flat per-pool tier made cost a restatement of pool count, so the gap is now
        # measured on pools directly -- the quantity the heuristic actually controls.
        best = None
        for groups in partitions(list(range(n)), profile.min_oligos, profile.max_oligos):
            if best is None or len(groups) < best:
                best = len(groups)

        if not plan["feasible"]:
            rows.append({"items": n, "heuristic": None, "optimal": best,
                         "reason": plan["reason"]})
            print(f"      n={n:>2}  heuristic infeasible; fewest possible pools {best}")
            continue

        got = plan["pool_count"]
        gap = (got - best) if best is not None else 0
        worst_gap = max(worst_gap, gap)
        rows.append({"items": n, "pools": got, "heuristic": got,
                     "optimal": best, "gap": gap})
        print(f"      n={n:>2}  heuristic {got} pool(s)  fewest possible {best}  "
              f"gap {gap}")

    # The heuristic must never be infeasible where a valid partition exists: that would be
    # the heuristic manufacturing the refusal, which is what plain first-fit did at n=9.
    manufactured = [r for r in rows if r["heuristic"] is None and r["optimal"] is not None]
    print(f"    cases refused by the heuristic where a valid partition exists: "
          f"{len(manufactured)}")
    print(f"    worst pool-count gap from optimal: {worst_gap}")

    print("\nE2  a stale rule set still answers eligibility")
    still_checks = check_eligibility(["ACGT" * 5] * 3, HISTORICAL_OPOOL)["eligible"]
    print(f"      eligibility works on the legacy profile: {still_checks}")
    print(f"      unresolved rules declared: {len(HISTORICAL_OPOOL.unresolved_rules)}")

    passed = not manufactured and still_checks
    print(f"\nGATE E: {'PASS' if passed else 'FAIL'}")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "gate_e.json").write_text(json.dumps({
        "pooling": rows,
        "worst_pool_count_gap_from_optimal": worst_gap,
        "heuristic_manufactured_refusals": manufactured,
        "eligibility_works_on_legacy_profile": still_checks,
        "no_order_placed": True,
        "scope": ("local implementation checks only; current vendor documentation remains "
                  "the authority and no order was placed or submitted"),
        "gate": "PASS" if passed else "FAIL",
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out / 'gate_e.json'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
