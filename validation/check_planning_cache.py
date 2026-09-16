"""Prove the planning cache changes nothing except how long planning takes.

The risk this checks is not that the cache is slow but that it is *wrong*: an entry reused
across targets whose inputs differ only in something the key omits would return a confident
plan built for another problem. So the test is cross-target, not repeat-call. For each
distinct candidate pool:

    1. clear the cache and plan the representative target      -> the cold answer
    2. clear, plan a *different* target that shares the pool   -> fills the cache
    3. plan the representative again                           -> must equal the cold answer

A repeat call in one process would pass even if the key were nonsense. Step 2 is what makes
step 3 a real hit.

Plans are compared whole -- every split's cuts, overhangs and score, in order -- because
order and fallback alternatives are load-bearing: `design_oneshot` walks the first
`OPTIMISE_ATTEMPTS` plans and takes the first whose constraints hold.

    python validation/check_planning_cache.py [--architectures 9 14 19] [--full]

Exits non-zero on any mismatch.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation"))

from integration import corpus                                        # noqa: E402
from clippr import constants as C                                     # noqa: E402
from clippr.arelf import balanced_cuts, safe_overhangs                # noqa: E402
from clippr.design import (MAX_FRAGMENT_AA, MIN_FRAGMENT_AA, SPLIT_BUDGET,  # noqa: E402
                           _plan_fragments)
from clippr.overhangs import (SCORER_VERSION, cache_stats, clear_cache,  # noqa: E402
                              enumerate_candidates, matrix_fingerprint)
from clippr.ppr import describe                                       # noqa: E402

MATRIX = "BsaI-HFv2"


def _n_fragments(protein: str) -> int:
    return max(1, -(-len(protein) // MAX_FRAGMENT_AA))


def pool_signature(protein: str, profile=C.DEFAULT_ENZYME_PROFILE) -> str:
    """Hash the ordered candidate pools a planner would search. Cheap by design.

    Measured at 0.005-0.107 s per protein against 12 s to search three of those pools at
    19S, which is the whole reason a cache is worth having.
    """
    n = _n_fragments(protein)
    splits = []
    for cuts in balanced_cuts(protein, n, MIN_FRAGMENT_AA, MAX_FRAGMENT_AA,
                              limit=SPLIT_BUDGET):
        try:
            options = enumerate_candidates(
                [",".join(safe_overhangs(protein, c, enzymes=profile)) for c in cuts])
        except ValueError as exc:
            options = {"error": str(exc)}
        splits.append([list(cuts), options])
    payload = {"n_fragments": n, "matrix": MATRIX, "profile": list(C.enzymes_for(profile)),
               "scorer": SCORER_VERSION, "table": matrix_fingerprint(MATRIX),
               "splits": splits}
    return hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()


def plan(protein: str):
    """The full ordered plan list, normalised for comparison."""
    plans, ceiling = _plan_fragments(protein, _n_fragments(protein), None, MATRIX,
                                     C.DEFAULT_ENZYME_PROFILE)
    return [[list(cuts), list(ohs), round(score, 12)] for cuts, ohs, score in plans], ceiling


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--architectures", nargs="*", type=int, default=[9, 14, 19])
    ap.add_argument("--full", action="store_true",
                    help="check every target, not one per distinct pool (hours at 19S)")
    args = ap.parse_args()

    targets = [t for t in corpus() if len(t) in args.architectures]
    proteins = {t: describe(t)["aa_sequence"] for t in targets}

    started = time.perf_counter()
    groups: dict[str, list[str]] = defaultdict(list)
    for t in targets:
        groups[pool_signature(proteins[t])].append(t)
    print(f"{len(targets)} targets -> {len(groups)} distinct candidate pools "
          f"({time.perf_counter() - started:.2f}s to classify)")
    for sig, members in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        lengths = sorted({len(m) for m in members})
        print(f"   {sig[:12]}  {len(members):>3} targets  lengths {lengths}")

    failures = []
    for sig, members in sorted(groups.items()):
        checked = members if args.full else members[:1]
        for target in checked:
            protein = proteins[target]

            clear_cache()
            t0 = time.perf_counter()
            cold, cold_ceiling = plan(protein)
            cold_seconds = time.perf_counter() - t0

            # Fill the cache from a *different* target sharing this pool, so the next call
            # is a genuine cross-target hit rather than a repeat of its own result.
            other = next((m for m in members if m != target), None)
            clear_cache()
            if other is not None:
                plan(proteins[other])
            t0 = time.perf_counter()
            warm, warm_ceiling = plan(protein)
            warm_seconds = time.perf_counter() - t0

            same = cold == warm and cold_ceiling == warm_ceiling
            if not same:
                failures.append((target, cold, warm))
            speedup = cold_seconds / warm_seconds if warm_seconds else float("inf")
            print(f"   {target:22s} plans={len(cold):>2} cold={cold_seconds:8.3f}s "
                  f"warm={warm_seconds:7.3f}s  x{speedup:6.1f}  "
                  f"{'MATCH' if same else '*** MISMATCH ***'}"
                  f"{'' if other else '  (only member; repeat-call only)'}")

    print("\nchecking that changed inputs are not served a stale plan")
    protein = proteins[targets[0]]
    clear_cache()
    base, _ = plan(protein)
    variants = {
        "different destination": lambda: _plan_fragments(
            protein, _n_fragments(protein), ("GGAG", "CGCT"), MATRIX,
            C.DEFAULT_ENZYME_PROFILE),
        "different matrix": lambda: _plan_fragments(
            protein, _n_fragments(protein), None, "BbsI-HF", C.DEFAULT_ENZYME_PROFILE),
        "different enzyme profile": lambda: _plan_fragments(
            protein, _n_fragments(protein), None, MATRIX, ("BsaI",)),
    }
    for name, call in variants.items():
        plans, _ = call()
        got = [[list(c), list(o), round(s, 12)] for c, o, s in plans]
        note = "differs from the cached plan" if got != base else "identical plan"
        print(f"   {name:26s} {note}")

    print(f"\ncache: {cache_stats()}")
    if failures:
        print(f"\nFAIL: {len(failures)} target(s) planned differently warm than cold")
        return 1
    print("PASS: every checked plan list is identical cold and warm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
