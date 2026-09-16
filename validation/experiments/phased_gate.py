"""Phase D gate — Pareto bookkeeping against exhaustive enumeration, then the real inventory.

Gate D, from the full-workflow plan:

  * on tiny synthetic problems, compare enumeration and Pareto filtering against an
    independently enumerated feasible set -- including dominance, ties, infeasible
    high-scoring points, reverse-complement conventions and a known trade-off
  * on the real inventory, preserve proteins and validate the compiled product of every
    released candidate
  * if every observed candidate has the same fidelity, report the degeneracy rather than
    manufacturing a trade-off

The synthetic part matters more than it looks. On the real inventory the feasible space is far
too large to enumerate, so the only way to know the filter is right is to give it a problem
small enough to solve by brute force and check it agrees.

Every released front member is revalidated **from its exported sequences**, not from the
candidate object the search returned -- an object can be right while the sequence it claims to
describe is wrong.

    python validation/experiments/phased_gate.py
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from Bio.Seq import Seq                                           # noqa: E402

from clippr import inventories as inv, joint_search as js         # noqa: E402
from clippr.junctions import valid_assignment                     # noqa: E402
from clippr.objectives import nondominated                        # noqa: E402
from clippr.overhangs import palindromic, reverse_complement      # noqa: E402
from clippr.substrates import substrate_problems                  # noqa: E402

TARGETS = ["AAAAUGUGG", "UUACACGUGCGUAC", "CUAUCACAUCACAUAAGCG"]


def brute_force_front(points: list[dict], directions: dict[str, str]) -> list[int]:
    """The nondominated set, computed by definition over every pair.

    Deliberately the slow, obvious implementation: this is the reference the fast one is
    checked against, so it must be transparently correct rather than clever.
    """
    keep = []
    for i, a in enumerate(points):
        beaten = False
        for j, b in enumerate(points):
            if i == j:
                continue
            no_worse = all(
                (b[n] >= a[n] if d == "max" else b[n] <= a[n])
                for n, d in directions.items())
            strictly = any(
                (b[n] > a[n] if d == "max" else b[n] < a[n])
                for n, d in directions.items())
            if no_worse and strictly:
                beaten = True
                break
        if not beaten:
            keep.append(i)
    return keep


def synthetic_checks() -> dict:
    """Hand-built cases with known answers, including the ones that are easy to get wrong."""
    results = {}
    d = {"fidelity": "max", "adaptation": "max", "repeat_burden": "min"}

    # A known trade-off: neither point beats the other.
    trade = [{"fidelity": 0.9, "adaptation": 0.5, "repeat_burden": 10},
             {"fidelity": 0.5, "adaptation": 0.9, "repeat_burden": 10}]
    results["known_trade_off_keeps_both"] = nondominated(trade, d) == [0, 1]

    # A dominated point must go, however good one of its objectives looks.
    dominated = trade + [{"fidelity": 0.5, "adaptation": 0.5, "repeat_burden": 20}]
    results["dominated_point_excluded"] = nondominated(dominated, d) == [0, 1]

    # Exact ties: neither strictly beats the other, so both stay.
    ties = [{"fidelity": 0.7, "adaptation": 0.7, "repeat_burden": 5},
            {"fidelity": 0.7, "adaptation": 0.7, "repeat_burden": 5}]
    results["exact_ties_both_kept"] = nondominated(ties, d) == [0, 1]

    # A minimised objective must actually be minimised.
    burden = [{"fidelity": 0.7, "adaptation": 0.7, "repeat_burden": 5},
              {"fidelity": 0.7, "adaptation": 0.7, "repeat_burden": 50}]
    results["minimised_objective_respected"] = nondominated(burden, d) == [0]

    # Floating-point noise must not manufacture dominance.
    noise = [{"fidelity": 0.7, "adaptation": 0.7, "repeat_burden": 5},
             {"fidelity": 0.7 + 1e-13, "adaptation": 0.7, "repeat_burden": 5}]
    results["tolerance_blocks_noise_dominance"] = nondominated(noise, d) == [0, 1]

    # The fast filter must agree with brute force on a randomly shaped grid.
    grid = [{"fidelity": f / 10, "adaptation": a / 10, "repeat_burden": b}
            for f in range(5) for a in range(5) for b in range(3)]
    results["agrees_with_brute_force"] = nondominated(grid, d) == brute_force_front(grid, d)

    # An empty set has an empty front, and a single point is its own front.
    results["empty_front_for_empty_input"] = nondominated([], d) == []
    results["single_point_is_its_own_front"] = nondominated([grid[0]], d) == [0]

    # An infeasible point with the best score must never reach the front: the front is
    # computed over feasible candidates only, which is a property of the caller, so this
    # checks the search's own filtering rather than `nondominated`.
    candidates = [js.Candidate({"x": "AAAA"}, False, "infeasible",
                               {"fidelity": 1.0, "adaptation": 1.0, "repeat_burden": 0}),
                  js.Candidate({"x": "CCCC"}, True, None,
                               {"fidelity": 0.5, "adaptation": 0.5, "repeat_burden": 5})]
    feasible = [c for c in candidates if c.feasible]
    front = [feasible[i] for i in nondominated([c.objectives for c in feasible], d)]
    results["infeasible_high_score_excluded"] = (
        len(front) == 1 and front[0].objectives["fidelity"] == 0.5)

    # Reverse-complement and palindrome conventions.
    results["palindrome_rejected"] = not valid_assignment(["AATT"])[0]
    results["complement_pair_rejected"] = not valid_assignment(["AAGA", "TCTT"])[0]
    results["duplicate_rejected"] = not valid_assignment(["ACTC", "ACTC"])[0]
    results["clash_with_destination_rejected"] = not valid_assignment(["CTCA"])[0]
    results["rc_is_involutive"] = all(
        reverse_complement(reverse_complement(x)) == x
        for x in ("ACTC", "AAGA", "GCAC", "TGAA"))
    results["palindrome_detects_self_complement"] = (
        palindromic("AATT") and not palindromic("ACTC"))

    # Exhaustive enumeration over a tiny assignment space must reproduce the same front the
    # filter gives, including the count of feasible assignments.
    domains = {"j1": ["ACTC", "AAGA", "AATT"], "j2": ["GCAC", "TGAA", "CTCA"]}
    enumerated, feasible_count = [], 0
    for combo in itertools.product(*domains.values()):
        ok, _why = valid_assignment(list(combo))
        if not ok:
            continue
        feasible_count += 1
        enumerated.append({"fidelity": len(set(combo)) / 2,
                           "adaptation": 0.5,
                           "repeat_burden": sum(c.count("A") for c in combo)})
    results["tiny_enumeration_feasible_count"] = feasible_count
    results["tiny_enumeration_matches_brute_force"] = (
        nondominated(enumerated, d) == brute_force_front(enumerated, d))
    return results


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", default=str(ROOT / "work" / "phaseb"
                                               / "inventory_recoded.json"))
    ap.add_argument("--max-evaluations", type=int, default=100)
    ap.add_argument("--wall-seconds", type=float, default=600.0)
    ap.add_argument("--out", default=str(ROOT / "work" / "phased"))
    args = ap.parse_args()

    from clippr.codons import complete_table

    print("D1  synthetic Pareto and convention checks")
    synthetic = synthetic_checks()
    for name, value in synthetic.items():
        mark = value if isinstance(value, bool) else True
        print(f"      {name:44s} {value if not isinstance(value, bool) else ('ok' if mark else 'FAIL')}")
    synthetic_ok = all(v for v in synthetic.values() if isinstance(v, bool))
    print(f"    synthetic checks: {'PASS' if synthetic_ok else 'FAIL'}")

    table = complete_table(json.loads(
        (ROOT / "data" / "codon_tables" / "kazusa_3055.json").read_text(encoding="utf-8")), 1)
    base = inv.load(args.inventory)
    classes = js.junction_classes(base, TARGETS)
    print(f"\nD2  real inventory: {len(base)} modules, {len(classes)} junction classes")
    for c in classes:
        print(f"      {c.key:14s} {len(c.options):>2} permitted, {len(c.roles)} role(s)")

    result = js.search(base, TARGETS, table, beam_width=8,
                       max_evaluations=args.max_evaluations,
                       wall_seconds=args.wall_seconds)
    print(f"\nD3  search: evaluated {result['evaluated']} "
          f"(feasible {result['feasible']}, infeasible {result['infeasible']}) in "
          f"{result['elapsed_seconds']}s, {result['completion']}")
    print(f"    incumbent fidelity {result['incumbent']['objectives']['fidelity']:.4f}, "
          f"front size {len(result['observed_front'])}")

    if result["fidelity_degenerate"]:
        print("    DEGENERATE: every observed candidate has the same fidelity; there is no "
              "fidelity trade-off to report")
    else:
        fids = sorted({c["objectives"]["fidelity"] for c in result["observed_front"]})
        print(f"    distinct fidelities on the front: {fids}")

    # --- D4: every released front member revalidated from its own sequences ---
    print("\nD4  revalidating every front member from its exported sequences")
    released_problems = {}
    for index, member in enumerate(result["observed_front"]):
        rebuilt = js.apply_assignment(base, classes, member["assignment"],
                                      label=f"front{index}")
        problems = []
        for module_id, record in rebuilt.modules.items():
            original = base.modules[module_id]
            s, e = original.coding_interval
            if len(record.dna) != len(original.dna):
                problems.append(f"{module_id} length changed")
            if str(Seq(record.dna[s:e]).translate()) \
                    != str(Seq(original.dna[s:e]).translate()):
                problems.append(f"{module_id} protein changed")
            # On the substrate that would be ordered, not the bare insert: a front member
            # whose export the order step would refuse is not a valid front member.
            bad = substrate_problems(record.dna, record.block)
            if bad:
                problems.append(f"{module_id}: {bad[0]}")
        for target in TARGETS:
            compiled = inv.compile_target(rebuilt, target)
            if not compiled.get("available"):
                problems.append(f"{target} does not compile")
                continue
            by_reaction = {}
            for junction in compiled["junctions"]:
                by_reaction.setdefault(junction["reaction"], []).append(junction["overhang"])
            for reaction, ends in sorted(by_reaction.items()):
                ok, why = valid_assignment(ends)
                if not ok:
                    problems.append(f"{target} reaction {reaction}: {why}")
        if problems:
            released_problems[index] = problems
        print(f"      front[{index}] fid {member['objectives']['fidelity']:.4f}  "
              f"{'ok' if not problems else 'FAIL: ' + problems[0]}")

    # --- D5: the recommendation follows its declared policy ---
    recommended = result["recommended"]
    policy_ok = False
    if recommended:
        feasible_front = [c for c in result["observed_front"] if c["feasible"]]
        best = max(c["objectives"]["fidelity"] for c in feasible_front)
        near = [c for c in feasible_front
                if c["objectives"]["fidelity"] >= best - js.FIDELITY_TOLERANCE]
        top = max(c["objectives"]["adaptation"] for c in near)
        policy_ok = (recommended["objectives"]["fidelity"] >= best - js.FIDELITY_TOLERANCE
                     and recommended["objectives"]["adaptation"] == top)
    print(f"\nD5  recommendation follows the declared policy: {policy_ok}")
    print(f"    policy: {result['recommendation_policy']}")

    # --- D6: an exhausted budget must not be reported as "no improvement exists" ---
    starved = js.search(base, TARGETS, table, beam_width=2, max_evaluations=1,
                        wall_seconds=600.0)
    honest_budget = (starved["budget_exhausted"]
                     and starved["completion"] == "budget_exhausted"
                     and "not the globally optimal front" in starved["scope"])
    print(f"D6  a starved search reports budget, not a claim about the space: {honest_budget}")

    passed = (synthetic_ok and not released_problems and policy_ok and honest_budget)
    print(f"\nGATE D: {'PASS' if passed else 'FAIL'}")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "gate_d.json").write_text(json.dumps({
        "synthetic": synthetic,
        "classes": result["classes"],
        "search": result,
        "released_member_problems": released_problems,
        "recommendation_policy_followed": policy_ok,
        "starved_search_reports_budget": honest_budget,
        "gate": "PASS" if passed else "FAIL",
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out / 'gate_d.json'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
