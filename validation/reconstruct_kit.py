"""Reconstruct the deposited-kit product for the six reference targets, and export it.

Produces the artefact the kit-versus-synthesis comparison has twice lacked: the actual
joined sequence, its stage products, every junction with coordinates, the reading frame and
what the frame implies about scope. A reviewer can re-derive all of it from the same
Supplementary Table S1 without running this script.

    python validation/reconstruct_kit.py --table "data/grasp_supp/Table S1.xlsx" \
        --outdir work/kit_products

Table S1 is not shipped with the package (see NOTICE.md); supply your own copy. Without it
this exits with a clear message rather than a traceback, and module selection still works.
"""
from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from clippr.homology import longest_shared                          # noqa: E402
from clippr.parts import select                                     # noqa: E402
from clippr.products import (InsertRecordsInvalid, InsertsUnavailable,  # noqa: E402
                             load_inserts, reconstruct, translate, write)

TARGETS = ["AAAAUGUGG", "GCUAAAGAC", "UUACACGUG", "CGUACGUAC", "AUCGAUCGA", "GGCCAAUUG"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", default=str(ROOT / "data" / "grasp_supp" / "Table S1.xlsx"))
    ap.add_argument("--outdir", default=str(ROOT / "work" / "kit_products"))
    ap.add_argument("--targets", nargs="*", default=TARGETS)
    args = ap.parse_args()

    try:
        sequences, source = load_inserts(args.table)
    except (InsertsUnavailable, InsertRecordsInvalid) as exc:
        print(f"cannot reconstruct: {exc}")
        return 2

    out = Path(args.outdir)
    print(f"insert records : {source.n_records} from {source.path}")
    print(f"sha256         : {source.sha256}\n")

    products, rejected = {}, {}
    for target in args.targets:
        plan = select(target)
        if not plan.available:
            rejected[target] = plan.reason
            print(f"{target:12s} REJECTED  {plan.reason}")
            continue
        rec = reconstruct(plan, sequences, source)
        tr = translate(rec, plan.modules)
        write(rec, out / f"{target}.json", translation=tr)
        products[target] = rec.final.sequence
        print(f"{target:12s} {len(rec.final.sequence):5d} nt  "
              f"{len(plan.modules):2d} modules  {len(rec.stages)} stages  "
              f"frame {tr.frame}  {len(tr.protein)} aa  {tr.stop_codons} stops")

    # The number two earlier attempts got wrong, now taken from exported sequence.
    worst = max(((longest_shared(products[a], products[b])[0], a, b)
                 for a, b in combinations(sorted(products), 2)), default=(0, "", ""))
    print(f"\nlongest shared tract between two kit products: {worst[0]} nt "
          f"({worst[1]} / {worst[2]})")
    print("measured on the exported sequences with homology.longest_shared, the same "
          "function the synthesis route uses")

    summary = {"source": source.as_dict(), "targets": args.targets,
               "reconstructed": sorted(products), "rejected": rejected,
               "worst_shared_nt": worst[0], "worst_pair": [worst[1], worst[2]]}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                                      encoding="utf-8")
    print(f"\nwrote {len(products)} product file(s) and summary.json to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
