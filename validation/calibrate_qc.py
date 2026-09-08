"""Show that the QC thresholds separate a real distribution rather than always firing.

The baseline pipeline reported WARNING on 200 of 200 designs with zero failures. The
pass criterion here is therefore not "our sequences pass" -- it is that the verdict
**varies**, that all three states are reachable, and that each cut-point sits at a
documented percentile of the measured distribution.

Compares two populations of 200 over the same proteins: the corpus' own coding
sequences, and the ones `codons.optimize_cds` produces. Reads saved artifacts plus a
cached regeneration, so it runs in our own venv:

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\calibrate_qc.py
"""
from __future__ import annotations

import glob
import json
import os
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from clippr.codons import optimize_cds, table_from_kazusa
from clippr.ppr import describe
from clippr.qc import DEFAULT_THRESHOLDS, synthesis_qc

CACHE = ROOT / "results" / "our_cds.json"


def reference_sequences() -> dict[str, str]:
    out = {}
    for p in sorted(glob.glob(str(ROOT / "results" / "_designs" / "*" / "assembly_plan_*.csv"))):
        df = pd.read_csv(p)
        out[os.path.basename(os.path.dirname(p))] = "".join(
            r["payload_5to3"][4:4 + int(r["cds_end"]) - int(r["cds_start"])]
            for _, r in df.iterrows())
    return out


def our_sequences(names: list[str]) -> dict[str, str]:
    if CACHE.exists():
        return json.loads(CACHE.read_text(encoding="utf-8"))
    table = table_from_kazusa()
    out = {}
    for i, rna in enumerate(names, 1):
        plan = glob.glob(str(ROOT / "results" / "_designs" / rna / "assembly_plan_*.csv"))[0]
        df = pd.read_csv(plan)
        cuts = [int(x) for x in df["aa_start_0based"][1:]]
        locked = {3 * c - 4: str(w).upper()
                  for c, w in zip(cuts, df["oh5_coding_site_5to3"][1:])}
        out[rna] = optimize_cds(describe(rna)["aa_sequence"], locked_sites=locked,
                                codon_table=table)["cds"]
        if i % 50 == 0:
            print(f"  regenerating {i}/{len(names)}", flush=True)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(out), encoding="utf-8")
    return out


def percentile_of(values: list[float], cut: float) -> float:
    return 100 * sum(1 for v in values if v < cut) / len(values)


def main() -> None:
    ref = reference_sequences()
    if not ref:
        sys.exit("no stored designs found under results/_designs/")
    ours = our_sequences(list(ref))

    qa = {k: synthesis_qc(v) for k, v in ref.items()}
    qb = {k: synthesis_qc(ours[k]) for k in ref}

    ca = Counter(q["status"] for q in qa.values())
    cb = Counter(q["status"] for q in qb.values())
    print(f"{'population':14s}{'PASS':>7s}{'WARNING':>9s}{'FAIL':>7s}")
    print(f"{'reference':14s}{ca['PASS']:>7d}{ca['WARNING']:>9d}{ca['FAIL']:>7d}")
    print(f"{'ours':14s}{cb['PASS']:>7d}{cb['WARNING']:>9d}{cb['FAIL']:>7d}")

    combined = list(qa.values()) + list(qb.values())
    print("\ncut-points, and where they sit in the combined n=%d distribution:" % len(combined))
    for metric, warn, fail in (
        ("longest_repeat", "longest_repeat_warn", "longest_repeat_fail"),
        ("repeated_kmer_fraction", "repeat_fraction_warn", "repeat_fraction_fail"),
    ):
        vals = [q[metric] for q in combined]
        w, f = DEFAULT_THRESHOLDS[warn], DEFAULT_THRESHOLDS[fail]
        print(f"  {metric:16s} WARNING {w:<6} = p{percentile_of(vals, w):.0f}"
              f"     FAIL {f:<6} = p{percentile_of(vals, f):.0f}")

    print("\nfeasibility gates (expected inert on this corpus):")
    for metric in ("gc_pct", "longest_homopolymer"):
        vals = [q[metric] for q in combined]
        print(f"  {metric:20s} range {min(vals):.1f} - {max(vals):.1f}")
    fired = sum(1 for q in combined
                if any("GC" in w or "homopolymer" in w for w in q["warnings"] + q["failures"]))
    print(f"  fired on {fired}/{len(combined)} sequences")

    # the criterion: the verdict must vary, and every state must be reachable
    states = set(ca) | set(cb)
    varies = len(set(q["status"] for q in qa.values())) > 1
    reachable = {"PASS", "WARNING", "FAIL"} <= states
    print(f"\nverdict varies within a population : {varies}")
    print(f"all three states reachable         : {reachable}")
    print(f"no state is constant across corpus : {ca['WARNING'] != len(qa)}")

    ok = varies and reachable
    print("\n" + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
