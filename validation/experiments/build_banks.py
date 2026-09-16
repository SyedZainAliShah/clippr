"""D2 — generate the per-target candidate banks the M2 comparison will select from.

One row per (target, seed). Every candidate is validated before it is banked: translation
back to the expected protein, the locked junction overhangs still present at their cuts, the
optimiser's own constraint verdict, and synthesis QC. A candidate that fails any of those is
recorded as a failure rather than dropped, so the denominator stays honest.

Candidates are deduplicated by sequence and the four counts the plan asks for are reported
separately: requested, generated, valid, distinct. If a target yields fewer than four
distinct valid candidates that is the result, not a prompt to keep searching -- searching
until the bank is full would make the bank size a function of how hard we looked.

Resumable and bounded through `clippr.experiment`, keyed on the frozen configuration's hash
so rows from a different configuration can never be merged into this bank.

    python validation/experiments/build_banks.py [--panels 9S 14S 19S] [--wall-seconds N]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from Bio.Data import CodonTable                                   # noqa: E402
from Bio.Seq import Seq                                           # noqa: E402

from clippr import experiment, manifest                           # noqa: E402
from clippr.design import design_oneshot                          # noqa: E402
from clippr.overhangs import SCORER_VERSION, cache_stats          # noqa: E402
from clippr.ppr import describe                                   # noqa: E402

CONFIG = Path(__file__).resolve().parent / "m2_joint_selection.json"


def cai(cds: str, table) -> float:
    """CAI over the codons whose amino acid has a synonymous choice, as in the corpus run."""
    forward = CodonTable.unambiguous_dna_by_id[1].forward_table
    weights = {c.replace("U", "T"): v / max(d.values())
               for aa, d in table.items() for c, v in d.items() if max(d.values()) > 0}
    used = [weights[cds[i:i + 3]] for i in range(0, len(cds), 3)
            if forward.get(cds[i:i + 3]) not in ("M", "W", None)]
    return math.exp(sum(math.log(w) for w in used) / len(used)) if used and min(used) > 0 else 0.0


def duplicated_20mers(cds: str) -> int:
    windows = [cds[i:i + 20] for i in range(len(cds) - 19)]
    return len(windows) - len(set(windows))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panels", nargs="*", default=["9S", "14S", "19S"])
    ap.add_argument("--wall-seconds", type=float, default=None,
                    help="overrides the budget frozen in the configuration")
    ap.add_argument("--outdir", default=str(ROOT / "work" / "m2"))
    args = ap.parse_args()

    raw = CONFIG.read_bytes()
    cfg = json.loads(raw)
    config_sha = hashlib.sha256(raw).hexdigest()
    shared = cfg["shared_settings"]
    table = json.loads((ROOT / shared["codon_table"]).read_text(encoding="utf-8"))

    # Identity must name the source that produced the rows, not only the configuration. A
    # journal keyed on config alone cannot tell two edits of the same dirty file apart, and
    # this checkout is uncommitted.
    identity = manifest.run_identity(
        table, shared["genetic_code"], matrix=shared["matrix"],
        experiment=cfg["experiment"], config_version=cfg["version"],
        config_sha256=config_sha, enzyme_profile=shared["enzyme_profile"],
        destination=shared["destination"], seeds=list(cfg["candidates"]["seeds"]))

    # The registered budget applies unless deliberately overridden on the command line.
    budget = (args.wall_seconds if args.wall_seconds is not None
              else cfg["budget"]["wall_seconds_total"])

    items = [(panel, target, seed)
             for panel in args.panels
             for target in cfg["panels"][panel]
             for seed in cfg["candidates"]["seeds"]]
    print(f"config {config_sha[:16]}  |  {len(items)} candidate slots requested")
    print(f"identity: {json.dumps(identity, sort_keys=True)[:110]}...\n")

    def build(item):
        panel, target, seed = item
        d = design_oneshot(target, codon_table=table,
                           genetic_code=shared["genetic_code"],
                           enzyme_profile=shared["enzyme_profile"],
                           matrix=shared["matrix"],
                           check_offtarget=shared["check_offtarget"], seed=seed)
        cds = d["cds"]
        expected = describe(target)["aa_sequence"]
        problems = []
        if str(Seq(cds).translate()) != expected:
            problems.append("translation does not match the expected protein")
        for cut, oh in zip(d["cuts"], d["junction_overhangs"]):
            if cds[3 * cut - 4:3 * cut] != oh:
                problems.append(f"locked overhang {oh} absent at cut {cut}")
        if not d["constraints_ok"]:
            problems.append("optimiser could not satisfy every constraint")
        if d["qc"]["status"] != "PASS":
            problems.append(f"synthesis QC {d['qc']['status']}")
        return {
            "panel": panel, "target": target, "seed": seed,
            "cds": cds, "cds_sha256": hashlib.sha256(cds.encode()).hexdigest()[:16],
            "cuts": list(d["cuts"]), "overhangs": list(d["junction_overhangs"]),
            "fidelity": d["fidelity"], "qc": d["qc"]["status"],
            "constraints_ok": bool(d["constraints_ok"]),
            "longest_repeat": d["qc"]["longest_repeat"],
            "duplicated_20mers": duplicated_20mers(cds),
            "cai": cai(cds, d["codon_table"]),
            "gc": round((cds.count("G") + cds.count("C")) / len(cds), 4),
            "valid": not problems, "problems": problems,
            "timings": d["timings"],
        }

    out = Path(args.outdir)
    result = experiment.run(
        items, build, path=out / "candidate_bank.jsonl", identity=identity,
        key=lambda it: f"{it[1]}:seed{it[2]}", wall_seconds=budget,
        on_progress=lambda i, n, k: print(f"  [{i:>3}/{n}] {k}", flush=True))

    rows = result["rows"]
    print(f"\nstatus {result['status']}  |  requested (this run) {result['requested']}"
          f"  new {result['completed_this_run']}  resumed {result['resumed']}"
          f"  rows total {result['completed']}  failures {len(result['failures'])}"
          f"  duplicates dropped {result['duplicates_dropped']}")
    print(f"budget {budget}s wall  |  spent {result['wall_seconds_all_runs']:.1f}s across all "
          f"runs  |  exceeded {result['budget_exceeded']}")
    print(f"planning cache: {cache_stats()}")

    print(f"\n{'target':22s} {'panel':>5} {'req':>4} {'gen':>4} {'valid':>6} {'distinct':>9}")
    summary = {}
    for target in dict.fromkeys(r["target"] for r in rows):
        mine = [r for r in rows if r["target"] == target]
        valid = [r for r in mine if r["valid"]]
        distinct = len({r["cds"] for r in valid})
        summary[target] = {"panel": mine[0]["panel"],
                           "requested": cfg["candidates"]["per_target"],
                           "generated": len(mine), "valid": len(valid), "distinct": distinct}
        flag = "" if distinct >= cfg["candidates"]["per_target"] else "   << short"
        print(f"{target:22s} {mine[0]['panel']:>5} {summary[target]['requested']:>4} "
              f"{len(mine):>4} {len(valid):>6} {distinct:>9}{flag}")

    short = [t for t, s in summary.items() if s["distinct"] < cfg["candidates"]["per_target"]]
    problems = Counter(p for r in rows for p in r["problems"])
    payload = {"identity": identity, "status": result["status"],
               "requested_this_run": result["requested"],
               "rows_total": result["completed"],
               "completed_this_run": result["completed_this_run"],
               "resumed": result["resumed"],
               "duplicates_dropped": result["duplicates_dropped"],
               "budget_wall_seconds": budget,
               "wall_seconds_all_runs": result["wall_seconds_all_runs"],
               "budget_exceeded": result["budget_exceeded"],
               "failures": result["failures"], "per_target": summary,
               "targets_short_of_four_distinct": short,
               "validation_problems": dict(problems),
               "planning_cache": cache_stats()}
    (out / "bank_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                                           encoding="utf-8")
    if short:
        print(f"\n{len(short)} target(s) yielded fewer than "
              f"{cfg['candidates']['per_target']} distinct valid candidates: {short}")
        print("recorded as the result; not searched further")
    if problems:
        print(f"\nvalidation problems: {dict(problems)}")
    print(f"\nwrote {out / 'candidate_bank.jsonl'} and bank_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
