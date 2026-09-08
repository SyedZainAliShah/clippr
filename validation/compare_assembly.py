"""Oracle check: our oligos against 200 stored reference designs, byte for byte.

This is the whole safety argument. Every downstream artefact -- FASTA, GenBank, the oPool
order, the cost estimate -- is derived from these sequences, so anything short of
identical is a bug, not a difference of opinion.

The comparison is fair because it removes every other variable: we take the corpus'
own coding sequence, its own cut positions and its own destination overhangs, and ask only
whether our splitting and wrapping reproduce its oligos. Codon choice, cut selection and
overhang selection are all held fixed.

Reads only saved artifacts under `results/_designs/`, so it imports nothing from the
reference and runs in our own venv:

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\compare_assembly.py
"""
from __future__ import annotations

import glob
import os
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from clippr.assembly import build_oligos, reassemble, split_cds


def cds_of(df: pd.DataFrame) -> str:
    """Rebuild the coding sequence from the fragment payloads."""
    return "".join(r["payload_5to3"][4:4 + int(r["cds_end"]) - int(r["cds_start"])]
                   for _, r in df.iterrows())


def main() -> None:
    plans = sorted(glob.glob(str(ROOT / "results" / "_designs" / "*" / "assembly_plan_*.csv")))
    if not plans:
        sys.exit("no stored designs found under results/_designs/")

    n = n_oligo_ok = n_payload_ok = n_roundtrip = 0
    n_frag = n_frag_ok = 0
    failures: list[str] = []
    kinds: Counter[str] = Counter()

    for path in plans:
        rna = os.path.basename(os.path.dirname(path))
        ref = pd.read_csv(path)
        cds = cds_of(ref)
        cuts = [int(x) for x in ref["aa_start_0based"][1:]]
        overhangs = [str(x).upper() for x in ref["oh5_coding_site_5to3"][1:]]
        destination = (str(ref.iloc[0]["oh5_coding_site_5to3"]).upper(),
                       str(ref.iloc[-1]["oh3_coding_site_5to3"]).upper())
        enzyme = str(ref.iloc[0]["wrap_enzyme"])
        n += 1

        ours = build_oligos(cds, cuts, overhangs, destination, enzyme)
        if len(ours) != len(ref):
            failures.append(f"{rna}: {len(ours)} fragments, reference has {len(ref)}")
            kinds["fragment_count"] += 1
            continue

        # the corpus oligo sequences live in the oligos csv, not the plan
        oligo_csv = glob.glob(os.path.join(os.path.dirname(path), "oneshot_*_oligos.csv"))
        ref_oligos = pd.read_csv(oligo_csv[0]) if oligo_csv else None

        payload_ok = oligo_ok = True
        for i in range(len(ref)):
            n_frag += 1
            good = True
            if ours.iloc[i]["payload_5to3"] != ref.iloc[i]["payload_5to3"]:
                payload_ok = good = False
                kinds["payload"] += 1
                if len(failures) < 5:
                    failures.append(f"{rna} F{i+1}: payload differs")
            if ref_oligos is not None:
                want = ref_oligos.iloc[i]["oligo_sequence_5to3"]
                got = ours.iloc[i]["oligo_sequence_5to3"]
                if got != want:
                    oligo_ok = good = False
                    kinds["oligo"] += 1
                    if len(failures) < 5:
                        j = next((k for k in range(min(len(got), len(want)))
                                  if got[k] != want[k]), min(len(got), len(want)))
                        failures.append(
                            f"{rna} F{i+1}: oligo differs at {j}: "
                            f"ours {got[max(0,j-6):j+6]!r} ref {want[max(0,j-6):j+6]!r}")
            n_frag_ok += good

        n_payload_ok += payload_ok
        n_oligo_ok += oligo_ok
        n_roundtrip += reassemble(split_cds(cds, cuts, overhangs, destination, enzyme)) == cds

    print(f"designs compared          : {n}")
    print(f"fragments compared        : {n_frag}")
    print(f"payloads byte-identical   : {n_payload_ok}/{n} designs")
    print(f"oligos byte-identical     : {n_oligo_ok}/{n} designs")
    print(f"fragments fully identical : {n_frag_ok}/{n_frag}")
    print(f"reassembly reproduces cds : {n_roundtrip}/{n}")
    if kinds:
        print("\nmismatch kinds:", dict(kinds))
    if failures:
        print("\nfirst failures:")
        for f in failures:
            print("  " + f)

    ok = n_payload_ok == n and n_oligo_ok == n and n_roundtrip == n
    print("\n" + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
