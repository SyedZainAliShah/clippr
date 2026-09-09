"""Derive the GRASP module inventory from the published supplementary tables.

PRIMARY SOURCE ONLY. Reads Farley et al. (2025) Supplementary Table S1 and emits
`src/clippr/parts.json`. Nothing is taken from the reference implementation, which ships the
same modules as its own Geneious plasmid exports; this reads the paper's published table
instead, and `validation/compare_parts.py` checks the two agree.

**Insert DNA is deliberately not emitted.** Selecting modules needs their identity, their
Golden Gate overhangs and their plate position -- not their sequence, because the point of the
route is that the user already holds the physical plasmids. Leaving the sequences out keeps
this a derived index rather than a republication of the paper's table, and anyone who wants
the sequences has the paper and Addgene.

Table S1 columns: Plasmid ID, Plasmid Vector, Plate ID, Insert Sequence, 5' overhang,
3' overhang.

Module naming, derived from the IDs themselves and confirmed against Table S2:

    pPR-1_<block>_L<last>5<fifth>        e.g. pPR-1_B_LD5N   -> block B, last D, fifth N
    pPR-1_1A_5<fifth>_<overhang>         the first block carries no last residue
    pPR-1_2E_L<last>                     the terminal block carries no fifth residue

Run once; re-run only if the source table changes.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SUP = ROOT / "data" / "grasp_supp"
OUT = ROOT / "src" / "clippr" / "parts.json"

#: 5th residue + last residue -> the base that pair reads. Same mapping as biology.py, repeated
#: here because this script must run standalone against the paper.
CODE_TO_BASE = {"TN": "A", "NN": "C", "TD": "G", "ND": "U"}


def parse_id(plasmid_id: str) -> dict:
    """block, last and fifth residues from a plasmid ID."""
    body = plasmid_id.replace("pPR-1_", "")
    # the first block: 1A_5N_AATG -- a fifth residue and an explicit 5' overhang, no last
    m = re.fullmatch(r"(1A)_5([NT])_([ACGT]{4})", body)
    if m:
        return {"block": m.group(1), "last": None, "fifth": m.group(2)}
    # the terminal block: 2E_LD -- a last residue only
    m = re.fullmatch(r"(2E)_L([DN])", body)
    if m:
        return {"block": m.group(1), "last": m.group(2), "fifth": None}
    # every other block: B_LD5N, 14E_LN5T, ...
    m = re.fullmatch(r"([0-9]*[A-Z])_L([DN])5([NT])", body)
    if m:
        return {"block": m.group(1), "last": m.group(2), "fifth": m.group(3)}
    raise ValueError(f"unrecognised module id: {plasmid_id}")


def main() -> None:
    s1 = pd.read_excel(SUP / "Table S1.xlsx")
    modules = []
    for _, row in s1.iterrows():
        pid = str(row["Plasmid ID"]).strip()
        info = parse_id(pid)
        insert = re.sub(r"[^ACGT]", "", str(row["Insert Sequence"]).upper())
        modules.append({
            "plasmid_id": pid,
            "vector": str(row["Plasmid Vector"]).strip(),
            "plate": str(row["Plate ID"]).strip(),
            "five_overhang": str(row["5ʹ overhang"]).strip().upper(),
            "three_overhang": str(row["3ʹ overhang"]).strip().upper(),
            "insert_length": len(insert),          # length only, never the sequence
            **info,
        })

    blocks = sorted({m["block"] for m in modules})
    graph: dict[str, list[str]] = {}
    for m in modules:
        graph.setdefault(m["five_overhang"], [])
        if m["three_overhang"] not in graph[m["five_overhang"]]:
            graph[m["five_overhang"]].append(m["three_overhang"])

    data = {
        "provenance": (
            "Derived from Farley et al. 2025 Supplementary Table S1 by tools/derive_parts.py. "
            "Module identity, overhangs and plate positions only; insert sequences are "
            "deliberately omitted -- see the paper or Addgene for those. No part of the "
            "reference implementation was used."
        ),
        "source": "Farley KV et al. (2025) Nucleic Acids Research, 10.1093/nar/gkaf1169",
        "kit": "GRASP Cloning Kit, 42 plasmids in pAGM1311",
        "code_to_base": CODE_TO_BASE,
        "code_note": ("A module named L<last>5<fifth> contributes the 5th residue for its own "
                      "base and the last residue for the preceding one; the base at position k "
                      "is read by (fifth of module k, last of module k+1). Derived from Table "
                      "S2 and validated on 28 of its 31 published variants."),
        "blocks": blocks,
        "overhang_graph": graph,
        "modules": modules,
    }
    OUT.write_text(json.dumps(data, indent=1, sort_keys=False) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}: {len(modules)} modules, {len(blocks)} blocks")
    print(f"  blocks: {', '.join(blocks)}")
    print(f"  overhang graph: {len(graph)} nodes")


if __name__ == "__main__":
    main()
