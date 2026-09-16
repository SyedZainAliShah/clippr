"""D3/D4 — compare three ways of choosing a library's DNA, under the frozen configuration.

    independent            the current default: one design per target at the base seed
    scaffold_diversified   homology.diversify_library with the corrected encoding re-roll
    joint_selected         choose jointly across targets from the banked candidates

The primary objective, tie-breaks, flag threshold and decision gate were fixed in
`m2_joint_selection.json` before any of this ran; its hash is recorded with every result.

**On precomputing the pair matrix.** The orchestrator's caution against it was sized for 50
targets and 8 candidates -- 78,400 exact comparisons, hours of work. These panels are 6
targets and 4 candidates: 15 target pairs x 16 candidate pairs = 240 comparisons, seconds.
With the matrix in hand the exhaustive optimum over 4^6 = 4096 assignments is table lookups,
so both a greedy and an exact answer are affordable and both are reported. The exhaustive
result is exact *for this panel size only*; it is not a claim that the approach scales, and
greedy is reported alongside precisely so its gap is visible.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from clippr import experiment, manifest                                  # noqa: E402
from clippr.homology import diversify_library, longest_shared, shared_kmers  # noqa: E402

CONFIG = Path(__file__).resolve().parent / "m2_joint_selection.json"


def profile(sequences: dict[str, str], threshold: int) -> dict:
    """The declared primary metric plus the declared secondary reports, for one library."""
    pairs = []
    for a, b in itertools.combinations(sorted(sequences), 2):
        length, start_a, start_b = longest_shared(sequences[a], sequences[b])[:3]
        pairs.append({"a": a, "b": b, "shared_nt": length,
                      "start_a": start_a, "start_b": start_b,
                      "shared_20mers": shared_kmers(sequences[a], sequences[b])})
    return {
        "worst_shared_tract_nt": max((p["shared_nt"] for p in pairs), default=0),
        "pairs_at_or_over_threshold": sum(1 for p in pairs if p["shared_nt"] >= threshold),
        "total_shared_20mers": sum(p["shared_20mers"] for p in pairs),
        "n_pairs": len(pairs),
        "pairs": pairs,
    }


def _rank(profile_dict, assignment) -> tuple:
    """The declared objective and tie-breaks, in order. Lower is better on every element."""
    return (profile_dict["worst_shared_tract_nt"],
            profile_dict["pairs_at_or_over_threshold"],
            profile_dict["total_shared_20mers"],
            tuple(assignment))


def greedy(targets, bank, distance, threshold) -> tuple[list[int], dict]:
    """Pick each target's candidate against those already fixed, in the declared order."""
    chosen: list[int] = []
    for position, target in enumerate(targets):
        best = None
        for index in range(len(bank[target])):
            worst = max((distance[(targets[k], chosen[k], target, index)]
                         for k in range(position)), default=0)
            key = (worst, index)
            if best is None or key < best[0]:
                best = (key, index)
        chosen.append(best[1])
    sequences = {t: bank[t][i] for t, i in zip(targets, chosen)}
    return chosen, profile(sequences, threshold)


def exhaustive(targets, bank, distance, threshold) -> tuple[list[int], dict]:
    """Exact optimum by the declared ranking. Affordable only at this panel size."""
    # Two passes. The first ranks every assignment using only the precomputed distances --
    # the earlier version re-profiled inside the loop, redoing exact comparisons it had
    # already paid for, which is why it cost minutes rather than lookups. The second
    # profiles only the assignments tying for the best primary value, which is the smallest
    # set on which the declared tie-breaks can be applied.
    scored, lowest = [], None
    for assignment in itertools.product(*(range(len(bank[t])) for t in targets)):
        worst = 0
        for i, j in itertools.combinations(range(len(targets)), 2):
            worst = max(worst, distance[(targets[i], assignment[i],
                                         targets[j], assignment[j])])
        if lowest is None or worst < lowest:
            lowest, scored = worst, [assignment]
        elif worst == lowest:
            scored.append(assignment)

    best = None
    for assignment in scored:
        prof = profile({t: bank[t][a] for t, a in zip(targets, assignment)}, threshold)
        key = _rank(prof, assignment)
        if best is None or key < best[0]:
            best = (key, list(assignment), prof)
    return best[1], best[2]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panels", nargs="*", default=["9S", "14S", "19S"])
    ap.add_argument("--bank", default=str(ROOT / "work" / "m2" / "candidate_bank.jsonl"))
    ap.add_argument("--outdir", default=str(ROOT / "work" / "m2"))
    ap.add_argument("--historical", action="store_true",
                    help="compare a bank whose identity differs, labelled as unverified")
    args = ap.parse_args()

    raw = CONFIG.read_bytes()
    cfg = json.loads(raw)
    config_sha = hashlib.sha256(raw).hexdigest()
    threshold = cfg["flag_threshold_nt"]
    shared_settings = cfg["shared_settings"]
    table = json.loads((ROOT / shared_settings["codon_table"]).read_text(encoding="utf-8"))

    # Verify the bank's identity header rather than skipping it. Discarding line 1 and
    # trusting the rest would let rows generated under a different configuration, scorer or
    # source be compared as though they were one experiment.
    lines = [l for l in Path(args.bank).read_text(encoding="utf-8").splitlines() if l.strip()]
    header = json.loads(lines[0]).get(experiment.IDENTITY_KEY)
    if header is None:
        raise SystemExit(f"{args.bank} has no identity header; it was not written by "
                         f"clippr.experiment and cannot be trusted as one experiment")
    if header.get("config_sha256") != config_sha:
        raise SystemExit(
            f"bank was built under configuration {str(header.get('config_sha256'))[:16]} "
            f"but this is {config_sha[:16]}; rebuild the bank rather than comparing across "
            f"configurations")
    rows, journal_failures = [], []
    for line in lines[1:]:
        record = json.loads(line)
        if experiment.FAILURE_KEY in record:
            journal_failures.append(record[experiment.FAILURE_KEY])
        elif experiment.RUN_KEY not in record:
            rows.append(record)
    # Printing the remaining identity fields without checking them let a bank declaring
    # scorer_version 999 and a different matrix fingerprint be compared as a matched run.
    current = manifest.run_identity(table, shared_settings["genetic_code"],
                                    matrix=shared_settings["matrix"])
    # `source_sha256` belongs in this set. Leaving it out meant a bank whose scorer, matrix
    # and codon table all matched was compared as a current run even when it declared a
    # source that never existed -- the one field that identifies the code which built it.
    checked = ("source_sha256", "scorer_version", "matrix_fingerprint",
               "codon_table_effective_sha256")
    mismatched = {k: (header.get(k), current[k]) for k in checked
                  if header.get(k) is not None and header[k] != current[k]}
    unpinned = [k for k in checked if header.get(k) is None]
    if mismatched and not args.historical:
        drift = ("The package source has changed since this bank was built, so it is not a "
                 "matched run even if every other field agrees. "
                 if "source_sha256" in mismatched else "")
        raise SystemExit(
            "bank identity does not match this source: "
            + "; ".join(f"{k} bank={b} now={n}" for k, (b, n) in mismatched.items())
            + ". " + drift
            + "Rerun the bank, or pass --historical to compare it as an explicitly "
              "unverified historical artefact.")
    print(f"bank identity: source {str(header.get('source_sha256'))[:16]} "
          f"scorer {header.get('scorer_version')} matrix {header.get('matrix_fingerprint')}")
    if unpinned:
        print(f"  NOT PINNED by this bank's header: {unpinned} -- it predates full run "
              f"identity and is not retroactively relabelled")
    if mismatched:
        print(f"  HISTORICAL MODE: identity differs on {sorted(mismatched)}; "
              f"results are not a matched current run")
    if journal_failures:
        print(f"bank journal records {len(journal_failures)} failure(s): "
              f"{[f['key'] for f in journal_failures][:5]}")
    print(f"config {config_sha[:16]}  |  bank rows {len(rows)}  |  threshold {threshold} nt\n")

    results = {}
    for panel in args.panels:
        targets = cfg["panels"][panel]
        bank = {}
        for t in targets:
            seqs = [r["cds"] for r in rows
                    if r["target"] == t and r["valid"] and r["panel"] == panel]
            bank[t] = list(dict.fromkeys(seqs))
        expected = cfg["candidates"]["per_target"]
        short = {t: len(v) for t, v in bank.items() if len(v) < expected}
        if any(not v for v in bank.values()):
            print(f"{panel}: bank incomplete, skipping"); continue
        if short:
            print(f"  coverage: {len(short)} target(s) below {expected} candidates: {short}")

        print(f"===== {panel}  ({len(targets)} targets, "
              f"{[len(bank[t]) for t in targets]} distinct candidates) =====")

        # 1 - independent: the first banked candidate is the base-seed design
        t0 = time.perf_counter()
        independent = profile({t: bank[t][0] for t in targets}, threshold)
        t_independent = time.perf_counter() - t0

        # 2 - scaffold diversification, generating its own designs as the method does
        t0 = time.perf_counter()
        div = diversify_library(targets, codon_table=table,
                                genetic_code=shared_settings["genetic_code"],
                                enzyme_profile=shared_settings["enzyme_profile"],
                                check_offtarget=shared_settings["check_offtarget"])
        diversified = profile(div["cds"], threshold)
        t_diversified = time.perf_counter() - t0

        # 3 - joint selection over the banked candidates
        t0 = time.perf_counter()
        distance = {}
        for a, b in itertools.combinations(targets, 2):
            for i, sa in enumerate(bank[a]):
                for j, sb in enumerate(bank[b]):
                    d = longest_shared(sa, sb)[0]
                    distance[(a, i, b, j)] = d
                    distance[(b, j, a, i)] = d
        t_matrix = time.perf_counter() - t0
        comparisons = len(distance) // 2

        t0 = time.perf_counter()
        g_pick, g_prof = greedy(targets, bank, distance, threshold)
        t_greedy = time.perf_counter() - t0
        r_pick, r_prof = greedy(list(reversed(targets)), bank, distance, threshold)
        t0 = time.perf_counter()
        e_pick, e_prof = exhaustive(targets, bank, distance, threshold)
        t_exhaustive = time.perf_counter() - t0

        incumbent_name, incumbent = min(
            [("independent", independent), ("scaffold_diversified", diversified)],
            key=lambda kv: _rank(kv[1], ()))
        # The gate is the preregistered method -- deterministic greedy -- and turns on a
        # strict improvement in the declared primary metric alone. Comparing the whole
        # ranking tuple let a secondary-only difference count as a primary win. Exhaustive
        # is a diagnostic at this panel size, not the registered production method.
        primary = "worst_shared_tract_nt"
        improved = g_prof[primary] < incumbent[primary]
        exhaustive_improved = e_prof[primary] < incumbent[primary]
        shipped_name, shipped = (("greedy", g_prof) if improved
                                 else (incumbent_name, incumbent))

        for name, prof, secs in (("independent", independent, t_independent),
                                 ("scaffold_diversified", diversified, t_diversified),
                                 ("joint greedy", g_prof, t_greedy),
                                 ("joint exhaustive", e_prof, t_exhaustive)):
            print(f"  {name:22s} worst {prof['worst_shared_tract_nt']:>4} nt   "
                  f"flagged {prof['pairs_at_or_over_threshold']:>2}/{prof['n_pairs']}   "
                  f"20mers {prof['total_shared_20mers']:>5}   {secs:7.2f}s")
        print(f"  {'greedy order dependence':22s} forward {g_pick} -> "
              f"{g_prof['worst_shared_tract_nt']} nt | reversed {list(reversed(r_pick))} -> "
              f"{r_prof['worst_shared_tract_nt']} nt")
        print(f"  pair matrix: {comparisons} comparisons in {t_matrix:.2f}s")
        print(f"  incumbent {incumbent_name} ({incumbent['worst_shared_tract_nt']} nt); "
              f"registered greedy {'IMPROVED' if improved else 'did not improve'} -> "
              f"shipping {shipped_name} ({shipped['worst_shared_tract_nt']} nt)")
        print(f"  exhaustive diagnostic "
              f"{'would improve' if exhaustive_improved else 'would not improve'} "
              f"({e_prof['worst_shared_tract_nt']} nt); not the registered method\n")

        results[panel] = {
            "targets": targets,
            "candidates_per_target": [len(bank[t]) for t in targets],
            "independent": independent, "scaffold_diversified": diversified,
            "joint_greedy": g_prof, "joint_exhaustive": e_prof,
            "greedy_assignment": g_pick, "greedy_reversed_assignment": list(reversed(r_pick)),
            "greedy_order_dependent": g_prof["worst_shared_tract_nt"] !=
                                      r_prof["worst_shared_tract_nt"],
            "exhaustive_assignment": e_pick,
            "incumbent": incumbent_name, "improved": improved,
            "exhaustive_would_improve": exhaustive_improved, "shipped": shipped_name,
            "gate_method": "greedy (preregistered); exhaustive is a diagnostic",
            "selected_cds": {
                "independent": {t: bank[t][0] for t in targets},
                "scaffold_diversified": div["cds"],
                "joint_greedy": {t: bank[t][i] for t, i in zip(targets, g_pick)},
                "joint_exhaustive": {t: bank[t][i] for t, i in zip(targets, e_pick)},
            },
            "seconds": {"independent": t_independent, "scaffold_diversified": t_diversified,
                        "pair_matrix": t_matrix, "greedy": t_greedy,
                        "exhaustive": t_exhaustive},
            "pair_comparisons": comparisons,
        }

    improved_panels = [p for p, r in results.items() if r["improved"]]
    gate = cfg["decision_gate"]
    # The gate is defined over all three panels. Judging it on a subset would let a partial
    # run read as a settled verdict, so say so rather than quietly reporting one.
    complete = len(results) == len(cfg["panels"])
    verdict = (("KEEP" if len(improved_panels) >= 2 else "DO NOT KEEP") if complete
               else f"PROVISIONAL ({len(results)} of {len(cfg['panels'])} panels run)")
    print(f"decision gate: improves on {len(improved_panels)} of {len(results)} panel(s) "
          f"{improved_panels} -> {verdict}")
    print(f"  rule: {gate['keep_the_feature_if']}")
    print(f"  otherwise: {gate['otherwise']}")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "method_comparison.json").write_text(
        json.dumps({"config_sha256": config_sha, "threshold_nt": threshold,
                    "panels": results, "improved_panels": improved_panels,
                    "panels_run": sorted(results), "gate_complete": complete,
                    "verdict": verdict}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {out / 'method_comparison.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
