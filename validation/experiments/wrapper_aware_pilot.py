"""Wrapper-aware optimisation pilot — can the ordered fragment satisfy the strict band?

The shipped one-shot route constrains GC on the **coding span**. What a vendor receives is
the destination's entry overhang, plus the coding slice, plus the enzyme wrapper. On the
saved corpus every design passes its coding-span checks and 277 of 1050 ordered fragments
breach the 0.35-0.65 local band, so all 200 designs are refused. 77 of those fragments are
F2-F7, so constraining only the outer gene boundaries cannot settle it.

This asks one question: with the geometry frozen and only synonymous codon choices changed,
can every ordered fragment satisfy the band the project already declares?

What is frozen, so that a difference is attributable to the codon choice alone:
protein sequence, genetic code, codon table, policy, fragment boundaries, overhang
assignment, and the fixed construction DNA. Nothing is relaxed. No shipped default is
changed: the constraint is passed through `optimize_cds(extra_constraints=...)`, which is
empty for every shipped route.

Budgets are declared in BUDGET below and recorded in the result before anything runs. An
exhausted budget is reported as "no feasible candidate found within this budget", which is
not the same statement as "infeasible" -- that would need a proof, and this is a search.

Nothing here measures a biological property, and no stage score is multiplied into a yield.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation" / "independent"))

from clippr import constants as C                                     # noqa: E402
from clippr.assembly import build_oligos                              # noqa: E402
from clippr.codons import optimize_cds, table_from_kazusa             # noqa: E402
from clippr.design import design_oneshot, final_artifact_report       # noqa: E402
from clippr.objectives import codon_adaptation                        # noqa: E402
from clippr.synthesis_profile import resolve                          # noqa: E402

#: Declared before execution. A run that exceeds these stops and says so.
BUDGET = {
    "targets": ["AAAAUGUGG", "UUACACGUGCGUAC", "CUAUCACAUCACAUAAGCG"],
    "seeds": [1, 42, 1234],
    "wall_seconds_per_attempt": 600,
    "wall_seconds_total": 5400,
}

#: The pilot's own policy reading, taken from the shipped profile rather than retyped.
PROFILE = "clippr-strict-legacy"


def _gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / max(1, len(seq))


def frozen_geometry(baseline: dict) -> list[dict]:
    """Each fragment's coding slice and the fixed DNA that flanks it when ordered.

    Derived from the baseline's own oligo table rather than rebuilt from the wrapper rules,
    so the pilot cannot disagree with production about where the payload sits. The slice
    must occur exactly once in the ordered string or the mapping is ambiguous and we refuse.
    """
    cds = baseline["cds"]
    out = []
    for _, row in baseline["oligos"].iterrows():
        full = str(row["oligo_sequence_5to3"])
        a, b = int(row["cds_start"]), int(row["cds_end"])
        chunk = cds[a:b]
        if full.count(chunk) != 1:
            raise ValueError(
                f"{row['fragment_id']}: its coding slice occurs {full.count(chunk)} times "
                f"in the ordered sequence, so the flanks cannot be identified")
        i = full.index(chunk)
        out.append({"fragment_id": str(row["fragment_id"]), "cds_start": a, "cds_end": b,
                    "left": full[:i], "right": full[i + len(chunk):]})
    return out


def make_constraint(fragments, policy):
    """A DNA Chisel constraint that scores the ordered fragment, not the coding span.

    Built here rather than in `src/` because this is an experiment. Promoting it would be a
    separate decision with its own evidence, and one passing example would not justify it.
    """
    import dnachisel as dc

    lo, hi = policy.local_gc
    window = policy.window

    class EnforceGCInOrderedFragments(dc.Specification):
        """Every GC window of every ordered fragment, flanks included."""

        def __init__(self, boost=1.0):
            self.boost = boost

        def evaluate(self, problem):
            seq = str(problem.sequence)
            breaches, worst = [], None
            for f in fragments:
                a, b = f["cds_start"], f["cds_end"]
                full = f["left"] + seq[a:b] + f["right"]
                offset = a - len(f["left"])        # full[i] is CDS base offset + i
                for i in range(0, max(1, len(full) - window + 1)):
                    gc = _gc(full[i:i + window])
                    if lo <= gc <= hi:
                        continue
                    if worst is None or abs(gc - (lo + hi) / 2) > abs(worst - (lo + hi) / 2):
                        worst = gc
                    # Clipped to the mutable coding span: the flanks cannot be changed, so
                    # pointing the solver at them would ask for the impossible.
                    start, end = max(a, offset + i), min(b, offset + i + window)
                    if end > start:
                        breaches.append((start, end))
            locations = [dc.Location(s, e) for s, e in sorted(set(breaches))]
            return dc.SpecEvaluation(
                self, problem, score=-len(locations), locations=locations,
                message=(f"{len(locations)} ordered-fragment window(s) outside "
                         f"({lo}, {hi})" + (f", worst {worst:.3f}" if worst else "")))

        def localized(self, location, problem=None, with_righthand=True):
            # Evaluated globally: one codon change moves a window in a fragment whose
            # flanks sit outside this location, so a localised view would miss it.
            #
            # `with_righthand` is passed by some of DNA Chisel's resolution paths. Omitting
            # it from the signature raised TypeError inside the solver, which reported
            # constraints_ok=False and handed back the unoptimised sequence -- read as
            # eight targets where the band was unsatisfiable. It was this.
            return self

        def label_parameters(self):
            return [("window", str(window)), ("gc", f"{lo}-{hi}")]

    return EnforceGCInOrderedFragments()


def attempt(target: str, seed: int, table, policy, budget_left: float,
            with_constraint: bool = True) -> dict:
    """One target at one seed: baseline, then the wrapper-aware candidate."""
    record = {"target": target, "seed": seed}
    t0 = time.perf_counter()

    baseline = design_oneshot(target, codon_table=table, check_offtarget=False, seed=seed)
    base_report = baseline["final_artifacts"]
    record["baseline"] = {
        "status": base_report["status"],
        "synthesis_ready": baseline["synthesis_ready"],
        "hard": list(base_report.get("hard") or []),
        "constraints_ok": baseline["constraints_ok"],
        "cai": round(codon_adaptation(baseline["cds"], table)["cai"], 4),
        "fidelity": baseline["fidelity"],
        "cds_sha256_prefix": baseline["cds"][:12],
    }

    geometry = frozen_geometry(baseline)
    record["frozen"] = {
        "cuts": list(baseline["cuts"]),
        "junction_overhangs": list(baseline["junction_overhangs"]),
        "fragments": [{"fragment_id": f["fragment_id"], "cds_start": f["cds_start"],
                       "cds_end": f["cds_end"], "left": f["left"], "right": f["right"],
                       "left_gc": round(_gc(f["left"]), 3),
                       "right_gc": round(_gc(f["right"]), 3)} for f in geometry],
    }

    # Exactly the baseline's inputs, plus the one new constraint.
    locked = {3 * c - 4: o for c, o in zip(baseline["cuts"],
                                           baseline["junction_overhangs"])}
    spent = time.perf_counter() - t0
    allowance = min(BUDGET["wall_seconds_per_attempt"], budget_left - spent)
    if allowance <= 0:
        record["candidate"] = {"outcome": "not attempted: total budget exhausted"}
        record["wall_seconds"] = round(time.perf_counter() - t0, 2)
        return record

    t1 = time.perf_counter()
    try:
        opt = optimize_cds(
            baseline["protein"], locked_sites=locked, codon_table=table,
            genetic_code=baseline["genetic_code"], seed=seed,
            enzymes=C.enzymes_for("igem_rfc1000"),
            extra_constraints=([make_constraint(geometry, policy)]
                               if with_constraint else []))
    except Exception as exc:                        # noqa: BLE001 - reported, not hidden
        record["candidate"] = {"outcome": f"solver raised {type(exc).__name__}: {exc}"}
        record["wall_seconds"] = round(time.perf_counter() - t0, 2)
        return record
    solve_seconds = time.perf_counter() - t1

    cds = opt["cds"]
    if cds == baseline["cds"]:
        record["candidate"] = {"outcome": "the solver returned the baseline unchanged"}

    # Rebuild on the frozen geometry and re-judge with production, not with our own copy
    # of the rule. The point is whether the shipped checker accepts the result.
    oligos = build_oligos(cds, baseline["cuts"], baseline["junction_overhangs"],
                          baseline.get("destination"), baseline["oligos"].iloc[0]
                          ["wrap_enzyme"], prefix=target)
    report = final_artifact_report(oligos, enzyme=baseline["oligos"].iloc[0]["wrap_enzyme"],
                                   enzymes=C.enzymes_for("igem_rfc1000"))

    # Geometry must still be the frozen one, or the comparison is not like for like.
    moved = [f["fragment_id"] for f, (_, row) in zip(geometry, oligos.iterrows())
             if (int(row["cds_start"]), int(row["cds_end"]))
             != (f["cds_start"], f["cds_end"])]
    # The construction DNA must be the same bytes, not merely the same coordinates. If a
    # flank moved, a GC win could have come from the wrapper rather than from the codons.
    reflanked = []
    for f, (_, row) in zip(geometry, oligos.iterrows()):
        full = str(row["oligo_sequence_5to3"])
        if not (full.startswith(f["left"]) and full.endswith(f["right"])):
            reflanked.append(f["fragment_id"])

    record["candidate"] = {
        **record.get("candidate", {}),
        "status": report["status"],
        "hard": list(report.get("hard") or []),
        "advisory": len(report.get("advisory") or []),
        "constraints_ok": opt["constraints_ok"],
        "optimiser_summary": (opt["summary"] or "")[-400:],
        "cai": round(codon_adaptation(cds, table)["cai"], 4),
        "protein_unchanged": _translates_to(cds, baseline["protein"],
                                            baseline["genetic_code"]),
        "geometry_unchanged": not moved and not reflanked,
        "moved_fragments": moved,
        "reflanked_fragments": reflanked,
        "cds_changed_from_baseline": cds != baseline["cds"],
        "solve_seconds": round(solve_seconds, 2),
    }
    record["export_row"] = {
        "target": target, "cds": cds, "protein": baseline["protein"],
        "oligos": json.loads(oligos.to_json(orient="records")),
    }
    record["candidate"]["accepted"] = bool(
        report["status"] in ("pass", "warning") and opt["constraints_ok"]
        and record["candidate"]["protein_unchanged"]
        and record["candidate"]["geometry_unchanged"])
    record["wall_seconds"] = round(time.perf_counter() - t0, 2)
    return record


def _translates_to(cds: str, protein: str, genetic_code: int) -> bool:
    """The invariant that must not be traded away for a GC win."""
    from Bio.Seq import Seq

    got = str(Seq(cds).translate(table=genetic_code))
    return got.rstrip("*") == protein.rstrip("*")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="JSON record path")
    ap.add_argument("--targets", default=None,
                    help="comma-separated override of the declared pilot targets")
    ap.add_argument("--seeds", default=None, help="comma-separated override of the seeds")
    ap.add_argument("--without-constraint", action="store_true",
                    help=("the negative control: re-optimise on the frozen geometry with "
                          "the wrapper-aware constraint omitted. If candidates still pass, "
                          "the constraint is not what fixed them and the result is void."))
    args = ap.parse_args()

    targets = ([t.strip() for t in args.targets.split(",") if t.strip()]
               if args.targets else list(BUDGET["targets"]))
    seeds = ([int(s) for s in args.seeds.split(",") if s.strip()]
             if args.seeds else list(BUDGET["seeds"]))

    policy = resolve(PROFILE)
    table = table_from_kazusa()

    result = {
        "experiment": "wrapper_aware_pilot",
        "question": ("with the geometry frozen and only synonymous codon choices changed, "
                     "can every ordered fragment satisfy the declared local GC band?"),
        "declared_budget": BUDGET,
        "targets_run": targets,
        "seeds_run": seeds,
        "policy": {"name": policy.name, "local_gc": list(policy.local_gc),
                   "window": policy.window,
                   "max_homopolymer": policy.max_homopolymer},
        "frozen": ["protein", "genetic_code", "codon_table", "policy",
                   "fragment_boundaries", "overhang_assignment", "construction_dna"],
        "changed": ["synonymous_codon_choices"],
        "wrapper_aware_constraint_applied": not args.without_constraint,
        "arm": ("treatment" if not args.without_constraint
                else "negative control: constraint omitted, everything else identical"),
        "scope": ("a search under a declared budget; nothing here is a biological "
                  "measurement and no stage score is combined into a yield"),
        "attempts": [],
    }

    started = time.perf_counter()
    for target in targets:
        for seed in seeds:
            left = BUDGET["wall_seconds_total"] - (time.perf_counter() - started)
            print(f"\n=== {target} seed {seed}  ({left:.0f}s of total budget left) ===")
            rec = attempt(target, seed, table, policy, left,
                          with_constraint=not args.without_constraint)
            result["attempts"].append(rec)
            base, cand = rec["baseline"], rec.get("candidate", {})
            print(f"    baseline  {base['status']:12s} CAI {base['cai']:.3f}"
                  f"  {len(base['hard'])} hard")
            if "status" in cand:
                print(f"    candidate {cand['status']:12s} CAI {cand['cai']:.3f}"
                      f"  {len(cand['hard'])} hard"
                      f"  accepted={cand['accepted']}"
                      f"  {cand['solve_seconds']:.1f}s")
                for h in cand["hard"][:2]:
                    print(f"      {h}")
            else:
                print(f"    candidate {cand.get('outcome')}")

    accepted = [a for a in result["attempts"]
                if a.get("candidate", {}).get("accepted")]
    elapsed = time.perf_counter() - started
    result["elapsed_seconds"] = round(elapsed, 2)
    result["accepted_attempts"] = len(accepted)
    result["total_attempts"] = len(result["attempts"])
    result["budget_exhausted"] = elapsed >= BUDGET["wall_seconds_total"]
    result["outcome"] = (
        f"{len(accepted)} of {len(result['attempts'])} attempts produced an ordered "
        f"fragment set the shipped checker accepts"
        if accepted else
        "no feasible candidate found within this budget; this is a search result, not a "
        "proof that none exists")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # The positive control the work order asked for: re-read the written export with the
    # independent checker, not with the production function that produced the verdict
    # above. A production "pass" confirmed by the same code is not independent.
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
        disagreements = len(independent.get("disagreements_with_production") or [])
        result["independent_recheck"] = {
            "candidates_written": str(jsonl),
            "report": str(verdict),
            "exit_code": proc.returncode,
            "accepted": independent.get("accepted"),
            "rejected": independent.get("rejected"),
            "not_assessed": independent.get("not_assessed"),
            "comparisons_performed": independent.get("comparisons_performed"),
            "disagreements": disagreements,
            "means": ("an independent implementation re-read the exported bytes; exit 0 "
                      "with 0 disagreements is the positive final-export control"),
        }
        print(f"\nindependent recheck: exit {proc.returncode}, "
              f"{independent.get('accepted')} accepted, "
              f"{independent.get('rejected')} rejected, "
              f"{disagreements} disagreement(s)")

    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print(f"\n{result['outcome']}")
    print(f"elapsed {elapsed:.0f}s of a declared {BUDGET['wall_seconds_total']}s budget")
    print(f"wrote {out}")
    if accepted:
        print("A positive control exists. It does NOT justify changing a default: that "
              "needs the predeclared broader panel, reported in full.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
