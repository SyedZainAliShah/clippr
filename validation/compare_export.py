"""Oracle check: our exported files against 200 stored reference designs.

Three things are checked, each for every design:

  1. `opool_quote` reproduces all seven numeric fields of the corpus quote
  2. the GenBank we write re-parses, and its sequence and CDS translation round-trip
  3. the oligo FASTA re-parses and every record matches the order CSV byte for byte

Reads only saved artifacts under `results/_designs/`, so it runs in our own venv:

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\compare_export.py
"""
from __future__ import annotations

import glob
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from clippr.assembly import build_oligos, split_cds
from clippr.export import opool_quote, write_fasta, write_genbank
from clippr.ppr import describe

QUOTE_FIELDS = ("n_oligos", "total_bases", "min_oligo_nt", "max_oligo_nt",
                "dna_eur", "phospho_eur", "total_eur")


def main() -> None:
    plans = sorted(glob.glob(str(ROOT / "results" / "_designs" / "*" / "assembly_plan_*.csv")))
    if not plans:
        sys.exit("no stored designs found under results/_designs/")

    tmp = Path(tempfile.mkdtemp())
    n = n_quote = n_gb = n_fa = 0
    bad: Counter[str] = Counter()
    failures: list[str] = []

    for path in plans:
        rna = os.path.basename(os.path.dirname(path))
        ref = pd.read_csv(path)
        cds = "".join(r["payload_5to3"][4:4 + int(r["cds_end"]) - int(r["cds_start"])]
                      for _, r in ref.iterrows())
        cuts = [int(x) for x in ref["aa_start_0based"][1:]]
        dest = (str(ref.iloc[0]["oh5_coding_site_5to3"]),
                str(ref.iloc[-1]["oh3_coding_site_5to3"]))
        n += 1

        frags = split_cds(cds, cuts, destination=dest)
        oligos = build_oligos(cds, cuts, destination=dest, prefix=rna)

        # 1. the quote
        q = opool_quote(oligos)
        qcsv = glob.glob(os.path.join(os.path.dirname(path), "oneshot_*_idt_opool.csv"))
        if qcsv:
            refq = pd.read_csv(qcsv[0]).iloc[0]
            diffs = [f for f in QUOTE_FIELDS if q[f] != refq[f]]
            if diffs:
                bad["quote"] += 1
                if len(failures) < 5:
                    f0 = diffs[0]
                    failures.append(f"{rna}: quote {f0} ours={q[f0]} ref={refq[f0]}")
            else:
                n_quote += 1

        # 2. GenBank round-trip
        gb = write_genbank(cds, frags, rna, tmp / f"{rna}.gb", destination=dest)
        rec = SeqIO.read(gb, "genbank")
        cds_feats = [f for f in rec.features if f.type == "CDS"]
        ok = (str(rec.seq) == cds and len(cds_feats) == 1
              and cds_feats[0].qualifiers["translation"][0] == str(Seq(cds).translate())
              and int(cds_feats[0].location.end) == len(cds))
        # every feature must lie inside the record
        ok = ok and all(0 <= int(f.location.start) <= int(f.location.end) <= len(cds)
                        for f in rec.features)
        n_gb += ok
        if not ok:
            bad["genbank"] += 1
            if len(failures) < 5:
                failures.append(f"{rna}: genbank round-trip failed")

        # 3. oligo FASTA
        fa = write_fasta(oligos, tmp / f"{rna}.fasta")
        recs = list(SeqIO.parse(fa, "fasta"))
        ok = (len(recs) == len(oligos)
              and all(str(r.seq) == o for r, o in
                      zip(recs, oligos["oligo_sequence_5to3"])))
        n_fa += ok
        if not ok:
            bad["fasta"] += 1

        # the annotation must agree with the PPR code, not just be well-formed
        d = describe(rna)
        labels = [f.qualifiers.get("label", [""])[0] for f in rec.features]
        for i, code in enumerate(d["code_pairs"], start=1):
            want = f"PPR{i} 5th {code[0]}"
            if want not in labels:
                bad["annotation"] += 1
                if len(failures) < 5:
                    failures.append(f"{rna}: missing feature {want!r}")
                break

    print(f"designs compared        : {n}")
    print(f"oPool quote matches     : {n_quote}/{n}")
    print(f"GenBank round-trips     : {n_gb}/{n}")
    print(f"oligo FASTA round-trips : {n_fa}/{n}")
    if bad:
        print("\nfailures by kind:", dict(bad))
    for f in failures:
        print("  " + f)

    ok = n_quote == n and n_gb == n and n_fa == n and not bad
    print("\n" + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
