"""Constraint-satisfaction check for codon optimisation, over the 200-design corpus.

Verified by constraint satisfaction, not by output matching: the corpus was optimised by
simulated annealing and we use DNA Chisel, so the two reach different valid answers and
comparing sequences would only measure that. What must hold is that every sequence we
emit is legal.

Reads only saved artifacts under `results/_designs/`, so it imports nothing from the
reference and runs in our own venv:

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\compare_codons.py [--limit N]

For every design it re-optimises the protein under that design's own cut positions and
asserts: translation preserved, locked overhangs verbatim, no blacklisted enzyme site,
no homopolymer over 4, GC inside the window. Failures are counted and reported, never
raised -- a protein the solver cannot satisfy is a result worth knowing.
"""
from __future__ import annotations

import argparse
import glob
import os
import statistics
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from clippr.codons import BLACKLIST_ENZYMES, optimize_cds, table_from_kazusa
from clippr.ppr import describe

MAX_HOMOPOLYMER = 4
GC_BOUNDS, GC_WINDOW = (0.35, 0.65), 50


def enzyme_sites() -> dict[str, list[str]]:
    """Recognition sites to search for, each with its reverse complement."""
    import dnachisel as dc

    comp = str.maketrans("ACGT", "TGCA")
    out = {}
    for e in BLACKLIST_ENZYMES:
        site = dc.EnzymeSitePattern(e).sequence.upper()
        out[e] = sorted({site, site.translate(comp)[::-1]})
    return out


def violations(cds: str, protein: str, locked: dict[int, str], sites) -> list[str]:
    bad = []
    if str(Seq(cds).translate()) != protein:
        bad.append("translation")
    for start, want in locked.items():
        if cds[start:start + len(want)] != want:
            bad.append(f"locked@{start}")
    for enzyme, patterns in sites.items():
        if any(p in cds for p in patterns):
            bad.append(f"site:{enzyme}")
    run = longest = 1
    for a, b in zip(cds, cds[1:]):
        run = run + 1 if a == b else 1
        longest = max(longest, run)
    if longest > MAX_HOMOPOLYMER:
        bad.append(f"homopolymer:{longest}")
    for i in range(0, max(1, len(cds) - GC_WINDOW + 1)):
        w = cds[i:i + GC_WINDOW]
        gc = sum(1 for b in w if b in "GC") / len(w)
        if not GC_BOUNDS[0] <= gc <= GC_BOUNDS[1]:
            bad.append(f"gc_window@{i}:{gc:.2f}")
            break
    return bad


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="only the first N designs")
    args = ap.parse_args()

    plans = sorted(glob.glob(str(ROOT / "results" / "_designs" / "*" / "assembly_plan_*.csv")))
    if args.limit:
        plans = plans[:args.limit]
    if not plans:
        sys.exit("no stored designs found under results/_designs/")

    table = table_from_kazusa()
    reasons: Counter[str] = Counter()
    n_ok = n_solver_ok = 0
    gcs, longest_repeats = [], []

    for i, path in enumerate(plans, 1):
        rna = os.path.basename(os.path.dirname(path))
        df = pd.read_csv(path)
        protein = describe(rna)["aa_sequence"]
        cuts = [int(x) for x in df["aa_start_0based"][1:]]
        locked = {3 * c - 4: str(w).upper()
                  for c, w in zip(cuts, df["oh5_coding_site_5to3"][1:])}

        res = optimize_cds(protein, locked_sites=locked, codon_table=table,
                           genetic_code=1, seed=42)
        n_solver_ok += bool(res["constraints_ok"])
        cds = res["cds"]
        bad = violations(cds, protein, locked, enzyme_sites())
        if bad:
            for b in bad:
                reasons[b.split("@")[0].split(":")[0]] += 1
        else:
            n_ok += 1
        gcs.append(sum(1 for b in cds if b in "GC") / len(cds))
        seen, longest = set(), 0
        for j in range(len(cds) - 19):
            if cds[j:j + 20] in seen:
                longest = 20
                break
            seen.add(cds[j:j + 20])
        longest_repeats.append(longest)
        if i % 25 == 0:
            print(f"  ... {i}/{len(plans)}", flush=True)

    n = len(plans)
    print(f"\ndesigns optimised        : {n}")
    print(f"solver reported all-pass : {n_solver_ok}/{n}")
    print(f"independently verified   : {n_ok}/{n}")
    print(f"GC content               : {min(gcs):.3f}-{max(gcs):.3f} "
          f"(mean {statistics.mean(gcs):.3f})")
    print(f"designs with a repeated 20-mer : {sum(1 for x in longest_repeats if x)}/{n}")
    if reasons:
        print("\nviolations by kind:")
        for k, v in reasons.most_common():
            print(f"  {k:16s} {v}")
    ok = n_ok == n
    print("\n" + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
