"""Per-stage runtime, measured rather than asserted.

Section 12 of the agreed plan asks for runtime under comparable conditions. Neither CLIPPR nor
the GRASP reference had any timing instrumentation: grepping the reference's 21 modules for
`perf_counter`, `elapsed` or equivalent finds nothing, and CLIPPR only ever recorded a single
figure for the 200-design corpus. "How long does the workflow take" therefore had no answer
from either side, and a whole-workflow speed comparison would have been invented.

**What this measures, in two passes that answer different questions.**

  *Isolated* — each stage run from the same fixed input, repeated. This is the right way to
  compare stages against each other, because every one starts from identical work.

  *Chained* — one honest end-to-end route where every stage consumes the previous stage's
  **actual output**: the collection optimiser's inventory is what the interface search explores,
  the selected inventory is what gets exported. Timed as a whole and per stage.

The distinction is the point, and an earlier version of this script got it wrong in a way worth
recording. Its docstring claimed stages were timed "in isolation and again as a chain" and no
chain existed: stage 04 called `optimise_library` and threw the returned inventory away, stage
06 ran a search whose result never reached an export, and the reported total was a **sum of
medians** — a number no single run ever took. A sum of isolated stage timings systematically
understates a workflow, because the later stages are cheap on the deposited kit and expensive
on the recoded, optimised one they actually receive.

The chain runs twice, under the default policy and under an explicit non-default one, because
a policy changes how much work the solver does and "the runtime" of a configurable pipeline is
not one number.

**What this deliberately does not do.**

  * It does not time the reference. Running their optimiser under our budgets would compare two
    different amounts of work and call it speed; a matched-budget comparison is a separate
    experiment with its own design.
  * It does not extrapolate. One machine, one interpreter, one inventory size; the numbers are
    an observation about this run, recorded with enough context to be repeated.
  * It reports every repeat, not a best-of. A minimum hides variance, and variance is what a
    user actually experiences.

Both user-visible routes are covered. The inventory route (recode a kit, compile from it)
runs by default; the one-shot synthesis route (design a CDS from an RNA target, no kit) is
behind `--oneshot` because a single 19S design runs into minutes and would triple the default
cost. That asymmetry is itself the most useful runtime fact this script produces.

    python validation/experiments/runtime_profile.py [--repeats 3] [--oneshot]
"""
from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import tempfile
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


def timed_isolated(args, raw, table) -> dict:
    """Every stage from the same fixed input, so the stages are comparable with each other."""
    from clippr import inventories as inv
    from clippr import joint_search as js
    from clippr import library_search as ls
    from clippr import ordering, substrates
    from clippr.recoding import recode_inventory

    watch = Stopwatch()
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
        watch.time("05 build junction classes", lambda: js.junction_classes(recoded, TARGETS))
        watch.time(
            "06 interface search (100 evals)",
            lambda: js.search(recoded, TARGETS, table, max_evaluations=100, wall_seconds=600))
        watch.time("07 build 42 order substrates", lambda: substrates.substrates_for(recoded))
        watch.time(
            "08 check order eligibility",
            lambda: ordering.check_eligibility(
                [s.sequence for s in substrates.substrates_for(recoded)],
                ordering.CURRENT_OPOOL_50PMOL))
    return watch.summary()


def timed_chain(args, raw, profile, label: str) -> dict:
    """One route, start to finish, each stage consuming what the previous one produced.

    Uses the workflow entry points rather than the engine, because those are what a user runs
    and what the published artefacts come from. Every inventory path handed to a stage is the
    path the previous stage wrote.
    """
    from clippr import workflow as w

    watch = Stopwatch()
    totals: list[float] = []
    outcome: dict = {}

    for repeat in range(args.repeats):
        out = Path(tempfile.mkdtemp(prefix=f"chain_{label}_{repeat}_"))
        started = time.perf_counter()

        recoded = watch.time("01 recode for host", lambda: w.recode_for_host(
            TABLE_S1, raw, out / "recode", label="runtime", profile=profile))
        inventory = recoded.artefacts["inventory"]

        optimised = watch.time("02 optimise collection", lambda: w.optimise_collection(
            inventory, raw, out / "optimise", mode="greedy", profile=profile,
            seed=42, max_proposals=42, wall_seconds=600))
        # The output, not the input. The previous version of this script timed this stage and
        # then explored the *pre-optimisation* inventory, so nothing downstream saw its work.
        inventory = optimised.artefacts["inventory"]

        front = watch.time("03 explore interfaces", lambda: w.explore_interfaces(
            inventory, TARGETS, raw, out / "explore", profile=profile,
            max_evaluations=100, wall_seconds=600))

        selected = watch.time("04 select interface", lambda: w.select_interface(
            inventory, front.artefacts["front"], out / "select", codon_table=raw))
        if selected.ok:
            inventory = selected.artefacts["inventory"]

        export = watch.time("05 export order items", lambda: w.order_items_for(
            inventory, out / "order", profile=profile))

        totals.append(time.perf_counter() - started)
        outcome = {
            "recoded_ok": recoded.ok,
            "optimiser_improved": optimised.data.get("improved"),
            "front_size": front.data.get("front_size", front.summary),
            "selection_ok": selected.ok,
            "selection_reason": (selected.failures or [{}])[0].get("reason") if not selected.ok
                                else None,
            "export_ok": export.ok,
            "export_failures": len(export.failures),
            "synthesis_profile": export.data.get("synthesis_profile")
                                 or (profile.name if profile is not None else "default"),
        }

    return {
        "stages": watch.summary(),
        "total_seconds": [round(t, 4) for t in totals],
        "total_median": round(statistics.median(totals), 3),
        # The number that matters, and the one the old script could not produce: a real
        # end-to-end time, not a sum of stages measured against a cheaper input.
        "sum_of_stage_medians": round(
            sum(s["median"] for s in watch.summary().values()), 3),
        "outcome_of_last_run": outcome,
    }


def timed_oneshot(args, raw) -> dict:
    """The synthesis route: one CDS designed from an RNA target, no inventory involved.

    Scaling is the point here, so all three architectures are timed rather than one. The cost
    is dominated by constraint solving over a CDS that grows with the target: a 19S design is a
    ~2 kb coding sequence carrying a local GC band, a homopolymer cap and three forbidden
    recognition sites at once.
    """
    from clippr.design import design_oneshot

    watch = Stopwatch()
    sizes = {}
    for target, arch in (("AAAAUGUGG", "9S"),
                         ("GCUAAAGACUUGCA", "14S"),
                         ("AAAGCGGCACUUGUGAAGU", "19S")):
        for _ in range(args.repeats):
            # Off-target screening is switched off here, and that is not a shortcut: on a
            # cache miss it fetches a genome from NCBI under a 90 s timeout, and a timing
            # script that can wait on a remote host is not measuring runtime. Screening cost
            # is a separate question and deserves its own measurement, with the fetch and the
            # scan reported apart rather than folded into a design time.
            #
            # (With screening on, a run here sat for over two minutes at 0.03 s of CPU -- that
            # is a stall on something, not computation. The chloroplast genome is cached on
            # disk, so "it was fetching" is not established; the cause was not diagnosed and
            # is not claimed. What is established is that the screening path can block, which
            # is reason enough to keep it out of a timing measurement.)
            got = watch.time(f"design {arch} ({len(target)} nt target)",
                             lambda t=target: design_oneshot(t, codon_table=raw,
                                                             check_offtarget=False))
        sizes[arch] = {"target_nt": len(target), "cds_nt": len(got["cds"]),
                       "protein_aa": len(got["protein"])}
    return {"stages": watch.summary(), "designs": sizes,
            "scope": ("codon design and fragmenting only; off-target screening is disabled "
                      "because it fetches a genome and would make this a network measurement")}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--out", default=str(ROOT / "work" / "runtime" / "runtime_profile.json"))
    ap.add_argument("--oneshot", action="store_true",
                    help="also time the synthesis route (all three architectures; slow -- a "
                         "19S design alone runs into minutes)")
    args = ap.parse_args()

    if not TABLE_S1.is_file():
        print(f"missing {TABLE_S1}; the deposited kit is needed to time the inventory route")
        return 2

    from clippr import synthesis_profile as sp
    from clippr.codons import complete_table

    raw = json.loads(CODON_TABLE.read_text(encoding="utf-8"))
    table = complete_table(raw, 1)

    print(f"{args.repeats} repeats, {args.seeds} recoding seeds, "
          f"{platform.python_version()} on {platform.machine()}\n")

    print("isolated stages (each from the same fixed input)")
    isolated = timed_isolated(args, raw, table)
    for name, got in isolated.items():
        print(f"  {name:36s} median {got['median']:7.3f}s   "
              f"range {got['min']:.3f}-{got['max']:.3f}")
    isolated_sum = sum(got["median"] for got in isolated.values())
    print(f"  {'sum of isolated medians':36s}        {isolated_sum:7.3f}s")
    print("  Not a workflow time: no run took this, and each stage saw a cheaper input than\n"
          "  the chain gives it.\n")

    routes = {}
    for label, profile in (("default", None), ("broad-experimental", sp.BROAD_EXPERIMENTAL)):
        print(f"chained route: {label}")
        got = timed_chain(args, raw, profile, label)
        routes[label] = got
        for name, stage in got["stages"].items():
            print(f"  {name:36s} median {stage['median']:7.3f}s   "
                  f"range {stage['min']:.3f}-{stage['max']:.3f}")
        print(f"  {'END TO END (measured, median)':36s}        {got['total_median']:7.3f}s")
        print(f"  outcome: {json.dumps(got['outcome_of_last_run'])}\n")

    oneshot = {}
    if args.oneshot:
        print("one-shot synthesis route (no inventory; this is the slow one)")
        oneshot = timed_oneshot(args, raw)
        for name, stage in oneshot["stages"].items():
            print(f"  {name:36s} median {stage['median']:8.3f}s   "
                  f"range {stage['min']:.3f}-{stage['max']:.3f}")
        print("  Cost is dominated by constraint solving and grows steeply with the target:")
        print("  the CDS lengthens and every window constraint applies across all of it.")
        print()

    payload = {
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor() or "unreported",
        },
        "settings": {"repeats": args.repeats, "recoding_seeds": args.seeds,
                     "targets": TARGETS},
        "isolated_stages": isolated,
        "sum_of_isolated_medians_seconds": round(isolated_sum, 3),
        "chained_routes": routes,
        "oneshot_synthesis_route": oneshot or {
            "measured": False,
            "reason": "not requested; pass --oneshot (it is slow)"},
        "how_to_read_this": (
            "The chained totals are measured end to end and are the runtime a user "
            "experiences. The isolated figures compare stages against each other and must "
            "not be summed into a workflow time -- each stage there receives a cheaper input "
            "than the chain hands it. The one-shot figures, when present, are a different "
            "route entirely -- designing a CDS from scratch rather than compiling one from a "
            "kit -- and are not comparable with the inventory stages."),
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
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
