"""Per-stage runtime, measured rather than asserted.

Section 12 of the agreed plan asks for runtime under comparable conditions. Neither CLIPPR nor
the GRASP reference had any timing instrumentation: grepping the reference's 21 modules for
`perf_counter`, `elapsed` or equivalent finds nothing, and CLIPPR only ever recorded a single
figure for the 200-design corpus. "How long does the workflow take" therefore had no answer
from either side, and a whole-workflow speed comparison would have been invented.

**What this measures.** Each user-visible stage separately, in the order a user runs them, with
the machine and interpreter recorded. Stages are timed in isolation and again as a chain, since
a sum of isolated timings is not a workflow time -- the later stages consume the earlier ones'
output and the cost of producing it is real.

**What it deliberately does not do.**

  * It does not time the reference. Running their optimiser under our budgets would compare two
    different amounts of work and call it speed; a matched-budget comparison is a separate
    experiment with its own design.
  * It does not extrapolate. One machine, one interpreter, one inventory size; the numbers are
    an observation about this run, recorded with enough context to be repeated.
  * It reports every repeat, not a best-of. A minimum hides variance, and variance is what a
    user actually experiences.

    python validation/experiments/runtime_profile.py [--repeats 3]
"""
from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

TABLE_S1 = ROOT / "data" / "grasp_supp" / "Table S1.xlsx"
CODON_TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"
TARGETS = ["AAAAUGUGG", "UUACACGUGCGUAC", "CUAUCACAUCACAUAAGCG"]


class Stopwatch:
    """Wall-clock per stage. `perf_counter` because it is monotonic and not wall-time-adjusted."""

    def __init__(self) -> None:
        self.timings: dict[str, list[float]] = {}

    def time(self, name: str, work):
        started = time.perf_counter()
        result = work()
        self.timings.setdefault(name, []).append(time.perf_counter() - started)
        return result

    def summary(self) -> dict:
        return {
            name: {
                "repeats": len(values),
                "seconds": [round(v, 4) for v in values],
                "median": round(statistics.median(values), 4),
                "min": round(min(values), 4),
                "max": round(max(values), 4),
                # Reported because a user experiences the spread, not the best case.
                "spread": round(max(values) - min(values), 4),
            }
            for name, values in sorted(self.timings.items())
        }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--out", default=str(ROOT / "work" / "runtime" / "runtime_profile.json"))
    args = ap.parse_args()

    if not TABLE_S1.is_file():
        print(f"missing {TABLE_S1}; the deposited kit is needed to time the inventory route")
        return 2

    from clippr import inventories as inv
    from clippr import joint_search as js
    from clippr import library_search as ls
    from clippr import ordering, substrates
    from clippr.codons import complete_table
    from clippr.recoding import recode_inventory

    raw = json.loads(CODON_TABLE.read_text(encoding="utf-8"))
    table = complete_table(raw, 1)
    watch = Stopwatch()

    print(f"{args.repeats} repeats, {args.seeds} recoding seeds, "
          f"{platform.python_version()} on {platform.machine()}\n")

    recoded = None
    for _ in range(args.repeats):
        deposited = watch.time("01 load deposited kit", lambda: inv.load_deposited(TABLE_S1))
        report = watch.time(
            "02 recode for host",
            lambda: recode_inventory(deposited, raw, label="runtime", seeds=args.seeds))
        recoded = report.inventory

        watch.time("03 compile 3 targets",
                   lambda: [inv.compile_target(recoded, t) for t in TARGETS])
        watch.time(
            "04 optimise collection (greedy, 42)",
            lambda: ls.optimise_library(recoded, raw, mode="greedy", seed=42,
                                        max_proposals=42, wall_seconds=600))
        classes = watch.time("05 build junction classes",
                             lambda: js.junction_classes(recoded, TARGETS))
        watch.time(
            "06 interface search (100 evals)",
            lambda: js.search(recoded, TARGETS, table, max_evaluations=100,
                              wall_seconds=600))
        watch.time("07 build 42 order substrates",
                   lambda: substrates.substrates_for(recoded))
        watch.time(
            "08 check order eligibility",
            lambda: ordering.check_eligibility(
                [s.sequence for s in substrates.substrates_for(recoded)],
                ordering.CURRENT_OPOOL_50PMOL))
        del classes

    summary = watch.summary()
    for name, got in summary.items():
        print(f"  {name:36s} median {got['median']:7.3f}s   "
              f"range {got['min']:.3f}-{got['max']:.3f}")

    chain_median = sum(got["median"] for got in summary.values())
    print(f"\n  {'sum of stage medians':36s}        {chain_median:7.3f}s")
    print("  A sum of medians is not a workflow time -- no single run took exactly this.")

    payload = {
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor() or "unreported",
        },
        "settings": {"repeats": args.repeats, "recoding_seeds": args.seeds,
                     "targets": TARGETS, "modules": len(recoded) if recoded else None},
        "stages": summary,
        "sum_of_stage_medians_seconds": round(chain_median, 3),
        "not_established": [
            "any comparison with the reference implementation; it is not timed here and "
            "running it under our budgets would compare different amounts of work",
            "behaviour on other machines, interpreters or inventory sizes",
            "notebook execution time, which the notebook harness reports separately",
        ],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
