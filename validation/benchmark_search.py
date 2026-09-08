"""Does exploring the design space beat simply not exploring it?

**Why this exists.** `search.py` evaluates several assembly plans in full and selects between
finished designs. That is more work than `design_oneshot`, which keeps the first plan that
satisfies its constraints. More work is not the same as a better answer, and the measured
degeneracy of the fidelity and QC objectives means the search has only one axis left to improve.
So the honest question is whether it beats cheaper alternatives on that axis, and the honest
outcome if it does not is to delete the module.

Four strategies, each given the same candidate pool for the same target:

    greedy      the first feasible plan -- what design_oneshot does
    random      one feasible plan chosen at random
    top2        the better of the two highest-fidelity plans
    search      the full selection over the whole pool

Reported on the optimiser score, which is the only objective measured to vary. Higher is better;
values are comparable within a target and meaningless across targets, so the comparison is by
within-target gap, never by averaging raw scores.

**Read the result the right way round.** `search` scoring a gap of 0.00 is partly definitional:
with fidelity and QC degenerate, `select` reduces to the argmax of the optimiser score, so it
cannot fail to find the best candidate in the pool it was given. The informative quantities are
therefore the *other* rows — how much the cheaper strategies leave behind on candidates that were
already computed and discarded. That is a statement about how costly the first-feasible heuristic
is, not proof that the selection rule is clever.

The gap is also in DNA Chisel's own objective units, so it is an internal measure of codon and
k-mer-uniqueness quality. It is not a validated biological improvement, and the relative figure
below is the honest way to size it.

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\benchmark_search.py [--targets N]
"""
from __future__ import annotations

import argparse
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from clippr.design import design_oneshot            # noqa: E402
from clippr.ppr import describe                     # noqa: E402
from clippr.search import explore, select           # noqa: E402
from integration import corpus                      # noqa: E402

BUDGET = 6


def strategies(cands, rng):
    """Each strategy's chosen candidate, from the same evaluated pool."""
    feasible = [c for c in cands if c.feasible]
    if not feasible:
        return {}
    out = {"greedy": feasible[0], "random": rng.choice(feasible)}
    out["top2"] = max(feasible[:2], key=lambda c: c.objectives_score)
    out["search"] = select(cands)[0]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=int, default=12)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    targets = corpus()[:args.targets]
    names = ["greedy", "random", "top2", "search"]
    gaps: dict[str, list[float]] = {n: [] for n in names}
    wins: dict[str, int] = {n: 0 for n in names}
    scale: list[float] = []          # typical |score| per target, to size the gap relatively
    n_done = 0

    for i, rna in enumerate(targets, 1):
        protein = describe(rna)["aa_sequence"]
        # One baseline design resolves the codon table, genetic code, effective enzyme profile
        # and fragment count exactly as design_oneshot would, so the pool the strategies choose
        # from is the pool design_oneshot itself would have seen.
        base = design_oneshot(rna, check_offtarget=False)
        cands = explore(protein, base["codon_table"],
                        n_fragments=base["n_fragments"], destination=None,
                        matrix="BsaI-HFv2",
                        genetic_code=base["genetic_code"],
                        enzyme_profile=base["enzyme_profile_effective"],
                        budget=BUDGET)
        picks = strategies(cands, rng)
        if not picks:
            continue
        n_done += 1
        best = max(c.objectives_score for c in cands if c.feasible)
        scale.append(abs(best))
        for n in names:
            gaps[n].append(best - picks[n].objectives_score)
        top = min(names, key=lambda n: best - picks[n].objectives_score)
        wins[top] += 1
        print(f"  {i}/{len(targets)} {rna}: "
              + "  ".join(f"{n} {best - picks[n].objectives_score:+.1f}" for n in names),
              flush=True)

    print(f"\n{n_done} targets, budget {BUDGET}")
    print(f"{'strategy':10s}{'mean gap to best':>18s}{'median':>10s}{'exact hits':>12s}")
    print("-" * 52)
    for n in names:
        g = gaps[n]
        if not g:
            continue
        print(f"{n:10s}{statistics.mean(g):>18.2f}{statistics.median(g):>10.2f}"
              f"{sum(1 for x in g if x < 1e-9):>8d}/{len(g)}")

    print("\nGap is distance from the best feasible optimiser score for that target, so 0.00")
    print("means the strategy found the best design in the pool. Lower is better.")
    print("\nsearch reaching 0.00 is definitional -- select() is the argmax on the only axis")
    print("that varies. The informative rows are the others: what the cheap strategies leave")
    print("behind on candidates that were already computed.")

    g_search = statistics.mean(gaps["search"]) if gaps["search"] else float("inf")
    g_greedy = statistics.mean(gaps["greedy"]) if gaps["greedy"] else float("inf")
    if scale:
        pct = 100.0 * (g_greedy - g_search) / statistics.mean(scale)
        print(f"\ngreedy leaves {g_greedy - g_search:.2f} optimiser-score points behind on "
              f"average, about {pct:.1f}% of the")
        print("typical score magnitude for these targets. Modest, real, and in DNA Chisel's")
        print("own units -- not a validated biological gain.")
    if g_search < g_greedy:
        print("\nExploring earns its runtime on the one objective that varies.")
    else:
        print("\nExploring does NOT beat the first-feasible plan. On this evidence the module")
        print("should be deleted rather than kept for the sake of the machinery.")
    sys.exit(0 if g_search <= g_greedy else 1)


if __name__ == "__main__":
    main()
