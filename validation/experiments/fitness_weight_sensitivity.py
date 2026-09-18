"""Does the Pareto front survive the weights I chose?

`synthesis_fitness` combines four terms under weights (0.45 / 0.25 / 0.15 / 0.15) that are a
declared engineering choice and **not a measurement**. A front that only exists at those
numbers would be an artefact of my judgement rather than a property of the candidates, and
nobody could tell the difference by reading the code.

So this re-runs the same search under several weightings and reports how much the front and the
recommendation move. Three outcomes, all worth knowing:

  * the recommendation is stable       the weights are not load-bearing; say so and move on
  * the front membership shifts a lot  the axis is real but the weights need justifying
  * a weighting degenerates the axis   that weighting is unusable and the flag should catch it

The equal-weight and single-term arms are deliberately extreme. A sensitivity sweep that only
tries nearby values cannot tell a robust result from a lucky one.

    python validation/experiments/fitness_weight_sensitivity.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

RECODED = ROOT / "work" / "phaseb" / "inventory_recoded.json"
CODON_TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"
TARGETS = ["AAAAUGUGG", "UUACACGUGCGUAC", "CUAUCACAUCACAUAAGCG"]

#: Weightings to try. The shipped one, an even split, and one arm per term taken to dominance.
WEIGHTINGS = {
    "shipped": {"gc_centrality": 0.45, "homopolymer_headroom": 0.25,
                "internal_repetition": 0.15, "collection_sharing": 0.15},
    "equal": {"gc_centrality": 0.25, "homopolymer_headroom": 0.25,
              "internal_repetition": 0.25, "collection_sharing": 0.25},
    "gc_only": {"gc_centrality": 1.0, "homopolymer_headroom": 0.0,
                "internal_repetition": 0.0, "collection_sharing": 0.0},
    "homopolymer_only": {"gc_centrality": 0.0, "homopolymer_headroom": 1.0,
                         "internal_repetition": 0.0, "collection_sharing": 0.0},
    "sharing_only": {"gc_centrality": 0.0, "homopolymer_headroom": 0.0,
                     "internal_repetition": 0.0, "collection_sharing": 1.0},
    "repetition_only": {"gc_centrality": 0.0, "homopolymer_headroom": 0.0,
                        "internal_repetition": 1.0, "collection_sharing": 0.0},
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--evaluations", type=int, default=100)
    ap.add_argument("--out", default=str(ROOT / "work" / "fitness" / "weight_sensitivity.json"))
    args = ap.parse_args()

    if not RECODED.is_file():
        print(f"missing {RECODED}; run Phase B first")
        return 2

    from clippr import inventories as inv
    from clippr import joint_search as js
    from clippr import synthesis_fitness as sf
    from clippr.codons import complete_table

    table = complete_table(json.loads(CODON_TABLE.read_text(encoding="utf-8")), 1)
    library = inv.load(RECODED)

    original = dict(sf.WEIGHTS)
    arms = {}
    try:
        for name, weights in WEIGHTINGS.items():
            sf.WEIGHTS.clear()
            sf.WEIGHTS.update(weights)
            spread = sf.inventory_fitness(library)
            result = js.search(library, TARGETS, table,
                               max_evaluations=args.evaluations, wall_seconds=600)
            front = result["observed_front"]
            recommended = result["recommended"]
            arms[name] = {
                "weights": dict(weights),
                "axis_spread": spread["spread"],
                "axis_degenerate": spread["degenerate"],
                "front_size": len(front),
                "front_assignments": sorted(
                    tuple(sorted(c["assignment"].items())) for c in front),
                "recommended_assignment": dict(sorted(recommended["assignment"].items())),
                "recommended_fidelity": recommended["objectives"]["fidelity"],
                "recommended_adaptation": recommended["objectives"]["adaptation"],
            }
            print(f"  {name:18s} spread {spread['spread']:.6f}  front {len(front):2d}  "
                  f"recommended fidelity {recommended['objectives']['fidelity']:.6f}  "
                  f"degenerate={spread['degenerate']}")
    finally:
        sf.WEIGHTS.clear()
        sf.WEIGHTS.update(original)

    shipped = arms["shipped"]
    same_recommendation = [n for n, a in arms.items()
                           if a["recommended_assignment"] == shipped["recommended_assignment"]]
    degenerate = [n for n, a in arms.items() if a["axis_degenerate"]]
    front_sizes = {n: a["front_size"] for n, a in arms.items()}

    print(f"\n  recommendation identical to shipped under: {sorted(same_recommendation)}")
    print(f"  front sizes: {front_sizes}")
    if degenerate:
        print(f"  WEIGHTINGS THAT DEGENERATE THE AXIS: {degenerate}")
        print("  -- those are unusable, and the degenerate flag caught them.")

    verdict = ("the recommendation does not depend on the weights"
               if len(same_recommendation) == len(arms)
               else "the recommendation moves with the weights; they are load-bearing "
                    "and need justifying")
    print(f"\n  {verdict}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "inventory_version": library.version,
        "evaluations": args.evaluations,
        "arms": arms,
        "recommendation_stable_under": sorted(same_recommendation),
        "degenerate_weightings": degenerate,
        "verdict": verdict,
        "not_established": [
            "that any weighting is correct; this measures sensitivity, not truth",
            "anything about inventories other than this one",
        ],
    }, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
