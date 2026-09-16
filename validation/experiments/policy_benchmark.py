"""Benchmark a synthesis policy change before anyone considers moving the default.

This is item 5 of the constraint-profile work order, and it exists because the measurements
that prompted the profile work were not good enough to justify a default change: they used 12
of 42 modules, changed two variables at once, and scored the best *coding* candidate rather
than a deliverable inventory.

**Two experiments, kept apart on purpose.**

    A  hard versus hard    strict (0.35-0.65) against broad (0.15-0.85), both enforced hard.
                           This is the one already run; it is repeated here over all 42 modules
                           with full-substrate filtering, matched seeds and matched budgets.

    B  hard versus soft    the *same* thresholds, enforced hard in one arm and as targets in
                           the other. Nothing measured so far says anything about this, and it
                           is the question that decides whether the reference's soft-penalty
                           architecture is worth adopting. Confounding it with A is what made
                           the earlier numbers unusable.

**Per-module differences, not a mean.** A mean hides which modules moved and by how much, and
a policy that improves 40 modules slightly while destroying two is not the same as one that
improves all 42. Every module's before/after is written out.

**Feasibility is judged under each arm's own profile**, on the full ordered substrate. An arm
that produces a higher CAI while returning substrates its own policy refuses has not won
anything, and the earlier 12-module experiment could not see that.

    python validation/experiments/policy_benchmark.py [--seeds 4] [--out work/policy]
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

TABLE_S1 = ROOT / "data" / "grasp_supp" / "Table S1.xlsx"
CODON_TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"


def measure(deposited, raw_table, table, profile, seeds: int, label: str) -> dict:
    """Recode the whole kit under one policy, and report what it actually delivers."""
    from clippr import substrates as sub
    from clippr.objectives import codon_adaptation
    from clippr.recoding import recode_inventory

    started = time.perf_counter()
    report = recode_inventory(deposited, raw_table, label=label, seeds=seeds, profile=profile)
    elapsed = time.perf_counter() - started

    per_module, hard, soft = {}, {}, {}
    for module_id, record in sorted(report.inventory.modules.items()):
        start, end = record.coding_interval
        after = codon_adaptation(record.dna[start:end], table)["cai"]
        original = deposited.modules[module_id]
        before = codon_adaptation(
            original.dna[slice(*original.coding_interval)], table)["cai"]
        per_module[module_id] = {"before": round(before, 6), "after": round(after, 6),
                                 "delta": round(after - before, 6),
                                 "changed": record.dna != original.dna}
        found = sub.substrate_findings(record.dna, record.block, profile)
        if found["hard"]:
            hard[module_id] = found["hard"]
        if found["target"]:
            soft[module_id] = found["target"]

    deltas = [row["delta"] for row in per_module.values()]
    return {
        "profile": profile.as_dict(),
        "seconds": round(elapsed, 2),
        "inventory_version": report.inventory.version,
        "improved": len(report.improved), "unchanged": len(report.unchanged),
        "mean_cai": round(statistics.mean(r["after"] for r in per_module.values()), 6),
        "median_delta": round(statistics.median(deltas), 6),
        "worst_delta": round(min(deltas), 6),
        "modules_worse": sum(1 for d in deltas if d < -1e-9),
        # The acceptance question: does this arm deliver what its own policy permits?
        "hard_breaches_under_own_profile": hard,
        "target_warnings_under_own_profile": soft,
        "per_module": per_module,
    }


def compare(left: dict, right: dict) -> dict:
    """Per-module difference between two arms, so a mean cannot hide a regression."""
    modules = sorted(set(left["per_module"]) & set(right["per_module"]))
    rows = []
    for module_id in modules:
        a, b = left["per_module"][module_id], right["per_module"][module_id]
        rows.append({"module": module_id, "left": a["after"], "right": b["after"],
                     "difference": round(b["after"] - a["after"], 6)})
    diffs = [row["difference"] for row in rows]
    return {
        "modules": len(rows),
        "right_ahead": sum(1 for d in diffs if d > 1e-9),
        "left_ahead": sum(1 for d in diffs if d < -1e-9),
        "tied": sum(1 for d in diffs if abs(d) <= 1e-9),
        "mean_difference": round(statistics.mean(diffs), 6) if diffs else None,
        "median_difference": round(statistics.median(diffs), 6) if diffs else None,
        "largest_gain": round(max(diffs), 6) if diffs else None,
        "largest_loss": round(min(diffs), 6) if diffs else None,
        "rows": rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--out", default=str(ROOT / "work" / "policy"))
    args = ap.parse_args()

    if not TABLE_S1.is_file():
        print(f"missing {TABLE_S1}; this benchmark needs the deposited kit")
        return 2

    from clippr import inventories as inv
    from clippr import synthesis_profile as sp
    from clippr.codons import complete_table

    raw = json.loads(CODON_TABLE.read_text(encoding="utf-8"))
    table = complete_table(raw, 1)
    deposited = inv.load_deposited(TABLE_S1)
    print(f"{len(deposited)} modules, {args.seeds} seeds per module, matched budgets\n")

    # --- the four arms ---
    # A: hard vs hard, different thresholds.  B: same thresholds, hard vs soft.
    soft_strict = sp.STRICT_LEGACY.narrowed(
        name="strict-thresholds-as-targets",
        enforcement={"local_gc": sp.TARGET, "homopolymer": sp.TARGET,
                     "forbidden_sites": sp.HARD},
        note="experiment B: the strict numbers, enforced softly")
    soft_broad = sp.BROAD_EXPERIMENTAL.narrowed(
        name="broad-thresholds-as-targets",
        enforcement={"local_gc": sp.TARGET, "homopolymer": sp.TARGET,
                     "forbidden_sites": sp.HARD},
        note="experiment B: the broad numbers, enforced softly")

    arms = {}
    for key, profile in (("strict_hard", sp.STRICT_LEGACY),
                         ("broad_hard", sp.BROAD_EXPERIMENTAL),
                         ("strict_soft", soft_strict),
                         ("broad_soft", soft_broad)):
        arms[key] = measure(deposited, raw, table, profile, args.seeds, f"bench-{key}")
        a = arms[key]
        print(f"  {key:12s} mean CAI {a['mean_cai']:.6f}  improved {a['improved']:2d}  "
              f"worse {a['modules_worse']:2d}  hard breaches {len(a['hard_breaches_under_own_profile'])}  "
              f"warnings {len(a['target_warnings_under_own_profile'])}  {a['seconds']:.1f}s")

    experiment_a = compare(arms["strict_hard"], arms["broad_hard"])
    experiment_b_strict = compare(arms["strict_hard"], arms["strict_soft"])
    experiment_b_broad = compare(arms["broad_hard"], arms["broad_soft"])

    print("\nA  hard vs hard, different thresholds (strict -> broad):")
    print(f"     broad ahead on {experiment_a['right_ahead']}, strict ahead on "
          f"{experiment_a['left_ahead']}, tied {experiment_a['tied']}; "
          f"mean {experiment_a['mean_difference']:+.6f}, "
          f"worst single module {experiment_a['largest_loss']:+.6f}")
    print("\nB  same thresholds, hard vs soft:")
    for name, got in (("strict", experiment_b_strict), ("broad", experiment_b_broad)):
        print(f"     {name:6s} soft ahead on {got['right_ahead']}, hard ahead on "
              f"{got['left_ahead']}, tied {got['tied']}; "
              f"mean {got['mean_difference']:+.6f}")

    # --- the acceptance question, stated rather than inferred ---
    verdict = {
        key: {"delivers_under_own_policy": not arm["hard_breaches_under_own_profile"],
              "warnings": len(arm["target_warnings_under_own_profile"])}
        for key, arm in arms.items()
    }
    print("\nevery arm must deliver an inventory its own policy accepts:")
    for key, got in verdict.items():
        print(f"     {key:12s} {got['delivers_under_own_policy']}  "
              f"({got['warnings']} advisory warning(s))")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "seeds": args.seeds,
        "modules": len(deposited),
        "deposited_version": deposited.version,
        "arms": arms,
        "experiment_a_hard_vs_hard": experiment_a,
        "experiment_b_hard_vs_soft": {"strict": experiment_b_strict,
                                      "broad": experiment_b_broad},
        "delivers_under_own_policy": verdict,
        "not_established": [
            "vendor acceptance of any arm; no sequence here has been synthesised",
            "that a higher CAI is a better protein",
            "anything about runtime beyond this recoding stage on this machine",
            "which profile should be the default; that is a decision, and this is its input",
        ],
    }
    (out / "policy_benchmark.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {out / 'policy_benchmark.json'}")
    print("\nThis is input to a decision, not the decision. No default is changed here.")
    return 0 if all(v["delivers_under_own_policy"] for v in verdict.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
