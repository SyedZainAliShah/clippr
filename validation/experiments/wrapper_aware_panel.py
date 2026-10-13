"""Wrapper-aware optimisation over the whole pinned 200-design panel.

The pilot (`PILOT_WRAPPER_AWARE_2026-10-10.md`) answered the feasibility question on three
targets: showing the solver the ordered fragment rather than the coding span recovered 9 of 9,
against 0 of 9 with the constraint removed. Three targets is not an answer about the panel,
and promoting the constraint to a default needs one.

This runs every target in the pinned corpus and **reports all of them**, not only the
recoveries. The corpus rows are the baselines: they were generated at seed 42 under a recorded
configuration, so re-optimising at the same seed on the same frozen geometry is a like-for-like
comparison and no baseline is recomputed.

Two arms:

  treatment   the wrapper-aware constraint applied
  control     the constraint omitted, everything else identical

The control is also a **setup check**. With the constraint off, the same seed, protein, table
and locked overhangs should return the corpus coding sequence bit for bit. If it does not, this
harness differs from whatever produced the corpus and the treatment arm proves nothing; the
mismatch is counted and reported rather than smoothed over.

Nothing here measures a biological property. Fidelity is reported to show it did **not**
change — the geometry is frozen, so it cannot — and no stage score is combined into a yield.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation" / "experiments"))

from ppropt import constants as C                                     # noqa: E402
from ppropt.assembly import build_oligos                              # noqa: E402
from ppropt.codons import optimize_cds, table_from_kazusa             # noqa: E402
from ppropt.design import final_artifact_report                       # noqa: E402
from ppropt.objectives import codon_adaptation                        # noqa: E402
from ppropt.synthesis_profile import resolve                          # noqa: E402
from wrapper_aware_pilot import (frozen_geometry, make_constraint,    # noqa: E402
                                 make_structure_constraint, fold_dg,
                                 STRUCTURE_DG_LIMIT, STRUCTURE_WINDOW_NT)

#: Declared before execution.
BUDGET = {
    # The pre-change corpus: this experiment measures against what the software produced
    # before the ordered-fragment constraint, so its baseline must not follow the
    # canonical name when that name moves to the current corpus.
    "panel": "work/independent/corpus200_prewrapper.jsonl",
    "designs": 200,
    "seed": 42,
    "wall_seconds_per_design": 300,
    "wall_seconds_total": 10800,
}

PROFILE = "ppropt-strict-legacy"
ENZYME_PROFILE = "igem_rfc1000"


def _rows(path: Path) -> list[dict]:
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        entry = json.loads(line)
        if "__run__" in entry or "__run_identity__" in entry:
            continue
        out.append(entry)
    return out


def repeat_metrics(seq: str, min_k: int = 8) -> dict:
    """The two repeat figures IDT's complexity report names.

    IDT flags exact repeats at **13 nt** and a repeat burden above **40%** of the sequence.
    Our uniquifier works at k=20, so repeats of 13-19 nt pass us and concern them. These
    measure that gap rather than arguing about it.
    """
    n = len(seq)
    index: dict[str, list[int]] = {}
    for i in range(n - min_k + 1):
        index.setdefault(seq[i:i + min_k], []).append(i)
    covered = bytearray(n)
    for positions in index.values():
        if len(positions) > 1:
            for p in positions:
                covered[p:p + min_k] = b"" * min_k

    longest = 0
    for length in range(min(120, n), min_k - 1, -1):
        seen = set()
        for i in range(n - length + 1):
            chunk = seq[i:i + length]
            if chunk in seen:
                longest = length
                break
            seen.add(chunk)
        if longest:
            break

    return {"longest_exact_repeat_nt": longest,
            "bases_in_repeats_pct": round(100.0 * sum(covered) / max(1, n), 1),
            "min_k": min_k}


def _baseline_frame(row: dict):
    """The corpus row's stored oligo table, as a frame, so geometry comes from the record."""
    import pandas as pd

    return pd.DataFrame(row["oligos"])


def one(row: dict, table, policy, with_constraint: bool,
        with_structure: bool = False, kmer: int | None = None) -> dict:
    """One corpus design: re-optimise on its frozen geometry and re-judge."""
    target = row.get("target") or row.get("target_rna")
    rec: dict = {"target": target, "architecture": row.get("architecture")}
    t0 = time.perf_counter()

    baseline = {"cds": row["cds"], "oligos": _baseline_frame(row),
                "protein": row["protein"], "cuts": row["cuts"],
                "junction_overhangs": row["junction_overhangs"]}
    enzyme = str(baseline["oligos"].iloc[0].get("wrap_enzyme") or "BsaI")

    base_report = final_artifact_report(baseline["oligos"], enzyme=enzyme,
                                        enzymes=C.enzymes_for(ENZYME_PROFILE))
    rec["baseline"] = {
        "status": base_report["status"],
        "hard": list(base_report.get("hard") or []),
        "hard_count": len(base_report.get("hard") or []),
        "cai": round(codon_adaptation(row["cds"], table)["cai"], 4),
        "fidelity": row.get("fidelity"),
        "cds_sha256": hashlib.sha256(row["cds"].encode()).hexdigest()[:16],
    }

    try:
        geometry = frozen_geometry(baseline)
    except ValueError as exc:
        rec["candidate"] = {"outcome": f"geometry not identifiable: {exc}"}
        rec["wall_seconds"] = round(time.perf_counter() - t0, 2)
        return rec

    locked = {3 * c - 4: o for c, o in zip(baseline["cuts"],
                                           baseline["junction_overhangs"])}
    t1 = time.perf_counter()
    try:
        opt = optimize_cds(
            baseline["protein"], locked_sites=locked, codon_table=table,
            genetic_code=1, seed=BUDGET["seed"],
            enzymes=C.enzymes_for(ENZYME_PROFILE),
            **({"unique_kmer_size": kmer} if kmer else {}),
            extra_constraints=(
                ([make_constraint(geometry, policy)] if with_constraint else [])
                + ([make_structure_constraint(geometry)] if with_structure else [])))
    except Exception as exc:                        # noqa: BLE001 - reported, not hidden
        rec["candidate"] = {"outcome": f"solver raised {type(exc).__name__}: {exc}"}
        rec["wall_seconds"] = round(time.perf_counter() - t0, 2)
        return rec
    solve = time.perf_counter() - t1

    cds = opt["cds"]
    oligos = build_oligos(cds, baseline["cuts"], baseline["junction_overhangs"],
                          None, enzyme, prefix=target)
    report = final_artifact_report(oligos, enzyme=enzyme,
                                   enzymes=C.enzymes_for(ENZYME_PROFILE))

    moved = [f["fragment_id"] for f, (_, r) in zip(geometry, oligos.iterrows())
             if (int(r["cds_start"]), int(r["cds_end"])) != (f["cds_start"], f["cds_end"])]
    reflanked = [f["fragment_id"] for f, (_, r) in zip(geometry, oligos.iterrows())
                 if not (str(r["oligo_sequence_5to3"]).startswith(f["left"])
                         and str(r["oligo_sequence_5to3"]).endswith(f["right"]))]

    from Bio.Seq import Seq
    translates = str(Seq(cds).translate()).rstrip("*") == baseline["protein"].rstrip("*")

    # Measured in both arms whether or not it was constrained, so the cost of adding the
    # term and the effect of omitting it are both visible.
    leading = [fold_dg(str(r["oligo_sequence_5to3"])[:STRUCTURE_WINDOW_NT])
               for _, r in oligos.iterrows()]
    base_leading = [fold_dg(str(r["oligo_sequence_5to3"])[:STRUCTURE_WINDOW_NT])
                    for _, r in baseline["oligos"].iterrows()]

    rec["candidate"] = {
        "status": report["status"],
        "hard": list(report.get("hard") or []),
        "hard_count": len(report.get("hard") or []),
        "advisory_count": len(report.get("advisory") or []),
        "constraints_ok": opt["constraints_ok"],
        "cai": round(codon_adaptation(cds, table)["cai"], 4),
        "protein_unchanged": translates,
        "geometry_unchanged": not moved and not reflanked,
        "moved_fragments": moved,
        "reflanked_fragments": reflanked,
        "cds_sha256": hashlib.sha256(cds.encode()).hexdigest()[:16],
        # The control's setup check: with the constraint off this must be True.
        "reproduces_corpus_cds": cds == row["cds"],
        "solve_seconds": round(solve, 2),
        "optimiser_summary": (opt["summary"] or "")[-300:],
        "repeats": repeat_metrics(cds),
        "leading_window_dg": [round(x, 2) for x in leading],
        "leading_window_breaches": sum(1 for x in leading if x < STRUCTURE_DG_LIMIT),
    }
    rec["baseline"]["repeats"] = repeat_metrics(row["cds"])
    rec["baseline"]["leading_window_dg"] = [round(x, 2) for x in base_leading]
    rec["baseline"]["leading_window_breaches"] = sum(
        1 for x in base_leading if x < STRUCTURE_DG_LIMIT)
    rec["candidate"]["accepted"] = bool(
        report["status"] in ("pass", "warning") and opt["constraints_ok"]
        and translates and not moved and not reflanked)
    rec["export_row"] = {"target": target, "cds": cds, "protein": baseline["protein"],
                         "architecture": row.get("architecture"),
                         "oligos": json.loads(oligos.to_json(orient="records"))}
    rec["wall_seconds"] = round(time.perf_counter() - t0, 2)
    return rec


def _median(xs):
    return round(statistics.median(xs), 1) if xs else None


def summarise(attempts: list[dict]) -> dict:
    """Every target, grouped. A recovery rate is meaningless without the denominator."""
    by_arch: dict[str, dict] = {}
    for a in attempts:
        arch = a.get("architecture") or "?"
        slot = by_arch.setdefault(arch, {"designs": 0, "accepted": 0, "rejected": 0,
                                         "not_attempted": 0})
        slot["designs"] += 1
        cand = a.get("candidate") or {}
        if "status" not in cand:
            slot["not_attempted"] += 1
        elif cand.get("accepted"):
            slot["accepted"] += 1
        else:
            slot["rejected"] += 1

    judged = [a for a in attempts if "status" in (a.get("candidate") or {})]
    accepted = [a for a in judged if a["candidate"]["accepted"]]
    failed = [a for a in judged if not a["candidate"]["accepted"]]
    # CAI over the ACCEPTED designs only. A design whose solver gave up returns the
    # unsolved sequence, whose CAI is near-random -- folding that into the distribution
    # reported a 0.41 "cost" that nobody traded for and that no accepted design paid.
    deltas = sorted(round(a["candidate"]["cai"] - a["baseline"]["cai"], 4)
                    for a in accepted)
    solves = sorted(a["candidate"]["solve_seconds"] for a in judged)
    remaining = sum(a["candidate"]["hard_count"] for a in judged)
    baseline_hard = sum(a["baseline"]["hard_count"] for a in attempts)

    def pct(xs, q):
        return xs[min(len(xs) - 1, int(q * len(xs)))] if xs else None

    return {
        "by_architecture": by_arch,
        "designs": len(attempts),
        "accepted": sum(1 for a in judged if a["candidate"]["accepted"]),
        "rejected": sum(1 for a in judged if not a["candidate"]["accepted"]),
        "not_attempted": len(attempts) - len(judged),
        "baseline_hard_findings_total": baseline_hard,
        "candidate_hard_findings_total": remaining,
        "cai_delta_accepted_only": {
            "scope": "the accepted designs; a solver failure is reported separately",
            "min": deltas[0] if deltas else None,
            "median": pct(deltas, 0.5),
            "max": deltas[-1] if deltas else None,
            "improved": sum(1 for d in deltas if d > 0),
            "worsened": sum(1 for d in deltas if d < 0),
            "worse_than_0_02": sum(1 for d in deltas if d < -0.02)},
        # Named, counted and diagnosed rather than averaged away.
        "failures": {
            "designs": len(failed),
            "all_reported_by_the_solver": all(
                not a["candidate"]["constraints_ok"] for a in failed) if failed else None,
            "by_architecture": {k: sum(1 for a in failed
                                       if a.get("architecture") == k)
                                for k in sorted({a.get("architecture") for a in failed})},
            "breached_floor_only": all(
                "0.340" in h or "0.34 " in h
                for a in failed for h in a["candidate"]["hard"]) if failed else None,
            "max_solve_seconds": max((a["candidate"]["solve_seconds"] for a in failed),
                                     default=None),
            "targets": [a["target"] for a in failed],
            "means": ("the solver reported constraints_ok=False and the final check "
                      "rejected the result; nothing orderable was written"),
        },
        "solve_seconds": {"min": solves[0] if solves else None,
                          "median": pct(solves, 0.5), "p95": pct(solves, 0.95),
                          "max": solves[-1] if solves else None},
        "invariant_violations": sum(
            1 for a in judged
            if not (a["candidate"]["protein_unchanged"]
                    and a["candidate"]["geometry_unchanged"])),
        "fidelity_note": ("the geometry is frozen, so the predicted ligation fidelity is "
                          "unchanged by construction; it is not re-reported as a gain"),
        "reproduces_corpus_cds": sum(1 for a in judged
                                     if a["candidate"]["reproduces_corpus_cds"]),
        "repeats": {
            "baseline_longest_nt": max((a["baseline"]["repeats"]["longest_exact_repeat_nt"]
                                        for a in judged), default=None),
            "candidate_longest_nt": max((a["candidate"]["repeats"]["longest_exact_repeat_nt"]
                                         for a in judged), default=None),
            "baseline_burden_pct_median": _median(
                [a["baseline"]["repeats"]["bases_in_repeats_pct"] for a in judged]),
            "candidate_burden_pct_median": _median(
                [a["candidate"]["repeats"]["bases_in_repeats_pct"] for a in judged]),
            "candidates_over_idt_13nt": sum(
                1 for a in judged
                if a["candidate"]["repeats"]["longest_exact_repeat_nt"] >= 13),
            "idt_limits": "exact repeats < 13 nt; bases in repeats < 40%",
        },
        "leading_window_breaches": {
            "baseline": sum(a["baseline"].get("leading_window_breaches", 0)
                            for a in judged),
            "candidate": sum(a["candidate"].get("leading_window_breaches", 0)
                             for a in judged),
            "limit": "dG < -8.0 kcal/mol in the first 48 nt",
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--without-constraint", action="store_true",
                    help="the control arm, and the setup check described in the docstring")
    ap.add_argument("--with-structure", action="store_true",
                    help=("also constrain the leading-window fold, the quantity GC is a "
                          "proxy for"))
    ap.add_argument("--kmer", type=int, default=None, metavar="K",
                    help=("uniquify all K-mers instead of the shipped default of 20. "
                          "IDT flags exact repeats from 13 nt, so repeats of 13-19 nt "
                          "currently pass us and concern them."))
    ap.add_argument("--limit", type=int, default=None, help="first N designs only")
    args = ap.parse_args()

    panel = ROOT / BUDGET["panel"]
    rows = _rows(panel)
    if args.limit:
        rows = rows[:args.limit]

    policy = resolve(PROFILE)
    table = table_from_kazusa()
    with_constraint = not args.without_constraint

    result = {
        "experiment": "wrapper_aware_panel",
        "arm": "treatment" if with_constraint else "control (constraint omitted)",
        "wrapper_aware_constraint_applied": with_constraint,
        "structure_constraint_applied": args.with_structure,
        "unique_kmer_size": args.kmer or "shipped default (20)",
        "structure_rule": {"window_nt": STRUCTURE_WINDOW_NT,
                           "dg_limit": STRUCTURE_DG_LIMIT,
                           "source": "Twist codon-optimisation FAQ, gene-synthesis guidance; "
                                     "no oPools structure rule exists to use"},
        "question": ("over the whole pinned panel, how many designs can satisfy the "
                     "declared local GC band on every ordered fragment using synonymous "
                     "codon choices alone?"),
        "declared_budget": BUDGET,
        "panel_sha256": hashlib.sha256(panel.read_bytes()).hexdigest(),
        "designs_in_panel": len(rows),
        "policy": {"name": policy.name, "local_gc": list(policy.local_gc),
                   "window": policy.window, "max_homopolymer": policy.max_homopolymer},
        "frozen": ["protein", "genetic_code", "codon_table", "policy",
                   "fragment_boundaries", "overhang_assignment", "construction_dna"],
        "changed": ["synonymous_codon_choices"],
        "reports": "every design in the panel, accepted or not",
        "attempts": [],
    }

    started = time.perf_counter()
    for n, row in enumerate(rows, 1):
        if time.perf_counter() - started > BUDGET["wall_seconds_total"]:
            result["stopped_early"] = (
                f"total wall budget of {BUDGET['wall_seconds_total']}s reached after "
                f"{n - 1} designs")
            break
        rec = one(row, table, policy, with_constraint,
                  with_structure=args.with_structure, kmer=args.kmer)
        result["attempts"].append(rec)
        cand = rec.get("candidate", {})
        mark = ("ok " if cand.get("accepted") else
                "   " if "status" in cand else "???")
        if n % 10 == 0 or not cand.get("accepted"):
            print(f"  [{mark}] {n:3d}/{len(rows)} {rec['target']:20s} "
                  f"{rec['architecture'] or '?':4s} "
                  f"base {rec['baseline']['status']:12s} -> "
                  f"cand {cand.get('status', cand.get('outcome', '?')):12s} "
                  f"{cand.get('solve_seconds', 0):6.1f}s")

    result["elapsed_seconds"] = round(time.perf_counter() - started, 2)
    result["summary"] = summarise(result["attempts"])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    exports = [a.pop("export_row") for a in result["attempts"] if "export_row" in a]
    if exports:
        jsonl = out.with_name(out.stem + "_candidates.jsonl")
        jsonl.write_text("\n".join(json.dumps(r) for r in exports) + "\n",
                         encoding="utf-8")
        verdict = out.with_name(out.stem + "_independent.json")
        proc = subprocess.run(
            [sys.executable,
             str(ROOT / "validation" / "independent" / "final_artifact_contract.py"),
             "--jsonl", str(jsonl), "--out", str(verdict),
             "--expect-designs", str(len(exports))],
            capture_output=True, text=True, cwd=str(ROOT))
        independent = (json.loads(verdict.read_text(encoding="utf-8"))
                       if verdict.is_file() else {})
        result["independent_recheck"] = {
            "candidates_written": str(jsonl), "report": str(verdict),
            "exit_code": proc.returncode,
            "accepted": independent.get("accepted"),
            "rejected": independent.get("rejected"),
            "not_assessed": independent.get("not_assessed"),
            "comparisons_performed": independent.get("comparisons_performed"),
            "disagreements": len(independent.get("disagreements_with_production") or []),
        }

    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    s = result["summary"]
    print(f"\n=== {result['arm']} ===")
    print(f"  designs            {s['designs']}")
    print(f"  accepted           {s['accepted']}")
    print(f"  rejected           {s['rejected']}")
    print(f"  not attempted      {s['not_attempted']}")
    print(f"  hard findings      {s['baseline_hard_findings_total']} baseline -> "
          f"{s['candidate_hard_findings_total']} candidate")
    ca = s["cai_delta_accepted_only"]
    print(f"  CAI (accepted)     min {ca['min']}  median {ca['median']}  max {ca['max']}  "
          f"({ca['improved']} improved, {ca['worsened']} worse, "
          f"{ca['worse_than_0_02']} worse than 0.02)")
    fl = s["failures"]
    if fl["designs"]:
        print(f"  failures           {fl['designs']}, all solver-reported="
              f"{fl['all_reported_by_the_solver']}, floor-only="
              f"{fl['breached_floor_only']}, max solve {fl['max_solve_seconds']}s")
    print(f"  solve seconds      median {s['solve_seconds']['median']}  "
          f"p95 {s['solve_seconds']['p95']}  max {s['solve_seconds']['max']}")
    print(f"  invariant breaches {s['invariant_violations']}")
    print(f"  reproduced corpus  {s['reproduces_corpus_cds']} of {s['designs']}")
    rp = s["repeats"]
    print(f"  longest repeat     {rp['baseline_longest_nt']} nt -> "
          f"{rp['candidate_longest_nt']} nt   (IDT flags >= 13)")
    print(f"  repeat burden      {rp['baseline_burden_pct_median']}% -> "
          f"{rp['candidate_burden_pct_median']}%   (IDT limit < 40%)")
    print(f"  designs >= 13 nt   {rp['candidates_over_idt_13nt']} of {s['designs']}")
    lw = s["leading_window_breaches"]
    print(f"  leading-window dG  {lw['baseline']} breaching -> {lw['candidate']}")
    for arch, v in sorted(s["by_architecture"].items()):
        print(f"    {arch:5s} {v['accepted']:3d} accepted / {v['designs']:3d}")
    if "independent_recheck" in result:
        ir = result["independent_recheck"]
        print(f"  independent        exit {ir['exit_code']}, {ir['accepted']} accepted, "
              f"{ir['rejected']} rejected, {ir['disagreements']} disagreement(s)")
    print(f"  elapsed            {result['elapsed_seconds']}s of a declared "
          f"{BUDGET['wall_seconds_total']}s")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
