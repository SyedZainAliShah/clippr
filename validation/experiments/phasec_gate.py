"""Phase C gate — the collection optimiser, compared against its own controls.

Gate C, from the full-workflow plan:

  * output inventories satisfy every hard constraint
  * objective deltas agree with a full independent recomputation
  * RNG reproducibility holds, and reload/recompile reproduces
  * no improvement returns the incumbent

Three modes on the same input, budgets and contexts, across three fixed seeds and both
processing orders. **Every run is preserved.** Picking the best seed after seeing the results
would turn a search into a lottery ticket chosen in hindsight, so all of them are reported and
the recommendation is made from the distribution.

`random` is a control, not a competitor: same budget, proposals accepted without regard to the
objective. An optimiser that cannot beat it has not been shown to work.

Objective deltas are recomputed here from the returned sequences, not read from the search's
own report -- a search grading its own improvement has measured nothing.

    python validation/experiments/phasec_gate.py [--seeds 3]
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from clippr import inventories as inv                             # noqa: E402
from clippr.library_search import MODES, optimise_library         # noqa: E402
from clippr.objectives import codon_adaptation, collection_penalty  # noqa: E402
from clippr.substrates import substrate_problems                  # noqa: E402

K = 20


def independent_objectives(library, table) -> dict:
    """Recompute the collection's objectives from the inventory's own sequences.

    Deliberately not `library_search.collection_objectives`: that is the function the search
    optimises against, and using it to check the search would only prove it is consistent
    with itself.
    """
    cais = []
    for record in library.modules.values():
        start, end = record.coding_interval
        cais.append(codon_adaptation(record.dna[start:end], table)["cai"])
    records = {m: r.dna for m, r in library.modules.items()}
    penalty = collection_penalty(records, k=K)["penalty"]
    pairs = max(1, len(records) * (len(records) - 1) // 2)
    return {"adaptation": round(sum(cais) / len(cais), 6),
            "sharing_penalty": penalty,
            "sharing_per_pair": round(penalty / pairs, 4)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", default=str(ROOT / "work" / "phaseb"
                                               / "inventory_recoded.json"))
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--max-proposals", type=int, default=84)
    ap.add_argument("--wall-seconds", type=float, default=600.0)
    ap.add_argument("--out", default=str(ROOT / "work" / "phasec"))
    args = ap.parse_args()

    from clippr.codons import complete_table

    raw_table = json.loads((ROOT / "data" / "codon_tables" / "kazusa_3055.json")
                           .read_text(encoding="utf-8"))
    table = complete_table(raw_table, 1)
    base = inv.load(args.inventory)
    baseline = independent_objectives(base, table)
    print(f"input: {len(base)} modules, version {base.version}")
    print(f"baseline (recomputed independently): adaptation {baseline['adaptation']}, "
          f"sharing/pair {baseline['sharing_per_pair']}\n")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    runs, violations, disagreements = [], {}, []

    seeds = [42 + i for i in range(args.seeds)]
    for mode in MODES:
        for seed in seeds:
            for reverse in (False, True):
                result = optimise_library(
                    base, raw_table, mode=mode, seed=seed, k=K,
                    wall_seconds=args.wall_seconds, max_proposals=args.max_proposals,
                    reverse_order=reverse)
                got = independent_objectives(result.inventory, table)

                # Every delivered module must satisfy the contract on the substrate that
                # would actually be ordered. Checking `r.dna` let this gate pass while 16 of
                # 18 runs returned an inventory the next step refused to export.
                bad = {m: substrate_problems(r.dna, r.block)
                       for m, r in result.inventory.modules.items()
                       if substrate_problems(r.dna, r.block)}
                if bad:
                    violations[f"{mode}/s{seed}/{'rev' if reverse else 'fwd'}"] = bad

                # The search's own report must agree with an independent recomputation.
                claimed = result.final_objectives["adaptation"]
                if abs(claimed - got["adaptation"]) > 5e-4:
                    disagreements.append(
                        {"run": f"{mode}/s{seed}", "claimed": claimed,
                         "independent": got["adaptation"]})

                record = {**result.as_dict(), "independent": got,
                          "reverse_order": reverse,
                          "sharing_delta": round(got["sharing_per_pair"]
                                                 - baseline["sharing_per_pair"], 4),
                          "adaptation_delta": round(got["adaptation"]
                                                    - baseline["adaptation"], 6)}
                runs.append(record)
                print(f"  {mode:7s} seed {seed} {'rev' if reverse else 'fwd'}  "
                      f"accepted {record['accepted']:>2}  "
                      f"adaptation {got['adaptation']:.4f} "
                      f"({record['adaptation_delta']:+.4f})  "
                      f"sharing/pair {got['sharing_per_pair']:.2f} "
                      f"({record['sharing_delta']:+.2f})  {record['completion']}")

    # --- reproducibility: the same seed must give the same sequences ---
    a = optimise_library(base, raw_table, mode="anneal", seed=42, k=K,
                         wall_seconds=args.wall_seconds, max_proposals=args.max_proposals)
    b = optimise_library(base, raw_table, mode="anneal", seed=42, k=K,
                         wall_seconds=args.wall_seconds, max_proposals=args.max_proposals)
    reproducible = a.inventory.version == b.inventory.version
    different_seed = optimise_library(base, raw_table, mode="anneal", seed=99, k=K,
                                      wall_seconds=args.wall_seconds,
                                      max_proposals=args.max_proposals)
    seed_matters = different_seed.inventory.version != a.inventory.version
    print(f"\nreproducible under the same seed: {reproducible}")
    print(f"a different seed explores differently: {seed_matters}")

    # --- no improvement must return the incumbent, not a worse inventory ---
    starved = optimise_library(base, raw_table, mode="greedy", seed=42, k=K,
                               wall_seconds=0.0, max_proposals=0)
    returns_incumbent = (starved.inventory.version == base.version
                         and not starved.improved and starved.budget_exhausted)
    print(f"a starved run returns the incumbent and says so: {returns_incumbent}")

    # --- reload and recompile ---
    saved = out / "inventory_greedy.json"
    inv.save(a.inventory, saved)
    reloaded = inv.load(saved)
    roundtrip = (reloaded.version == a.inventory.version
                 and inv.compile_target(reloaded, "AAAAUGUGG")["product"]
                 == inv.compile_target(a.inventory, "AAAAUGUGG")["product"])
    print(f"reload reproduces and recompiles: {roundtrip}")

    # --- mode comparison, over every preserved run ---
    print("\nby mode (all runs, none discarded):")
    summary = {}
    for mode in MODES:
        subset = [r for r in runs if r["mode"] == mode]
        shares = [r["sharing_delta"] for r in subset]
        adapts = [r["adaptation_delta"] for r in subset]
        summary[mode] = {
            "runs": len(subset),
            "sharing_delta_median": round(statistics.median(shares), 4),
            "sharing_delta_best": round(min(shares), 4),
            "adaptation_delta_median": round(statistics.median(adapts), 6),
            "accepted_median": statistics.median(r["accepted"] for r in subset),
            "seconds_median": round(statistics.median(r["elapsed_seconds"] for r in subset), 2),
        }
        print(f"  {mode:7s} n={len(subset)}  sharing delta median "
              f"{summary[mode]['sharing_delta_median']:+.3f} (best "
              f"{summary[mode]['sharing_delta_best']:+.3f})  adaptation delta median "
              f"{summary[mode]['adaptation_delta_median']:+.5f}")

    beats_control = (summary["anneal"]["sharing_delta_median"]
                     < summary["random"]["sharing_delta_median"])
    print(f"\nanneal beats the random control on median sharing: {beats_control}")
    if not beats_control:
        print("  -- reported as measured. An optimiser that cannot beat its control on "
              "these\n     settings has not been shown to work, and that is the finding.")

    passed = (not violations and not disagreements and reproducible and seed_matters
              and returns_incumbent and roundtrip)
    print(f"\nGATE C: {'PASS' if passed else 'FAIL'}")
    if violations:
        print(f"  constraint violations in {len(violations)} run(s)")
    if disagreements:
        print(f"  {len(disagreements)} objective disagreements with independent recomputation")

    (out / "gate_c.json").write_text(json.dumps({
        "input_version": base.version, "baseline": baseline,
        "runs": runs, "by_mode": summary,
        "reproducible_same_seed": reproducible, "different_seed_differs": seed_matters,
        "starved_returns_incumbent": returns_incumbent, "reload_roundtrip": roundtrip,
        "constraint_violations": violations,
        "objective_disagreements": disagreements,
        "anneal_beats_random_control": beats_control,
        "note": ("objectives recomputed independently of library_search.collection_objectives;"
                 " all runs preserved, no seed selected after the fact"),
        "gate": "PASS" if passed else "FAIL",
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out / 'gate_c.json'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
