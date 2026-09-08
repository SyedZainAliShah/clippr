"""Oracle check: our cut geometry against 200 stored reference designs.

Runs entirely off saved artifacts in `results/_designs/`, so it imports nothing from the
third-party code and can run in our own venv:

    PYTHONUTF8=1 python validation/compare_cuts.py

For every fragment boundary in every stored assembly plan it asserts three things:

  1. the concatenated fragment protein equals what `clippr.ppr.describe` builds
  2. `arelf.cut_context(cds, cut_aa)` equals the overhang the corpus actually uses
  3. that overhang is reachable by synonymous substitution from the protein alone,
     i.e. `arelf.achievable_overhangs` contains it

(3) is the load-bearing one: it shows our codon model spans the space the corpus
searched, which is what lets us optimise cuts without reproducing its heuristic.
"""
from __future__ import annotations

import glob
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from clippr.arelf import achievable_overhangs, cut_context, find_motifs, offsets_from_motif
from clippr.ppr import describe


def cds_of(df: pd.DataFrame) -> str:
    """Rebuild the coding sequence from the fragment payloads.

    payload = [5' overhang] + cds[start:end], and the final fragment carries the
    destination 3' overhang after its coding span.
    """
    parts = []
    for _, r in df.iterrows():
        span = int(r["cds_end"]) - int(r["cds_start"])
        parts.append(r["payload_5to3"][4:4 + span])
    return "".join(parts)


def main() -> None:
    plans = sorted(glob.glob(str(ROOT / "results" / "_designs" / "*" / "assembly_plan_*.csv")))
    if not plans:
        sys.exit("no stored designs found under results/_designs/")

    n_prot = n_prot_ok = 0
    n_cut = n_ctx_ok = n_reach_ok = 0
    failures: list[str] = []
    all_offsets: set[int] = set()

    for path in plans:
        rna = os.path.basename(os.path.dirname(path))
        df = pd.read_csv(path)
        protein = "".join(df["aa_sequence"])
        cds = cds_of(df)

        n_prot += 1
        if protein == describe(rna)["aa_sequence"]:
            n_prot_ok += 1
        else:
            failures.append(f"{rna}: protein differs from clippr.ppr.describe")

        if len(cds) != 3 * len(protein):
            failures.append(f"{rna}: cds {len(cds)} nt != 3 x {len(protein)} aa")
            continue

        cuts = [int(x) for x in df["aa_start_0based"][1:]]
        for cut, want in zip(cuts, df["oh5_coding_site_5to3"][1:]):
            n_cut += 1
            got = cut_context(cds, cut)
            if got == want:
                n_ctx_ok += 1
            elif len(failures) < 5:
                failures.append(f"{rna} cut {cut}: cut_context {got} != reference {want}")
            if want in achievable_overhangs(protein, cut):
                n_reach_ok += 1
            elif len(failures) < 5:
                failures.append(f"{rna} cut {cut}: {want} unreachable synonymously")

        all_offsets |= set(offsets_from_motif(protein, cuts))

    print(f"designs compared        : {n_prot}")
    print(f"protein matches ours    : {n_prot_ok}/{n_prot}")
    print(f"cuts compared           : {n_cut}")
    print(f"cut_context exact       : {n_ctx_ok}/{n_cut}")
    print(f"overhang reachable      : {n_reach_ok}/{n_cut}")
    print(f"ARELF-relative offsets  : {sorted(all_offsets)}")

    if failures:
        print("\nfailures:")
        for f in failures:
            print("  " + f)
    ok = n_prot_ok == n_prot and n_ctx_ok == n_cut and n_reach_ok == n_cut
    print("\n" + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
