"""Derive GRASP scaffold constants from the published supplementary tables.

PRIMARY SOURCE ONLY. Reads Farley et al. (2025) Supplementary Tables S1 and S2 and
reconstructs the wild-type protein by assembling the deposited modules per the published
recipe. Nothing here is copied from any third-party implementation.

  Table S1  42 deposited modules: insert DNA, 5' and 3' overhangs
  Table S2  per protein variant, which modules compose CDS1 and CDS2

Module ids encode their variable residues: "<block>_L<last>5<fifth>", e.g. B_LD5N has
last=D, fifth=N. Table S2 refers to them as "<block>_<last><fifth>", e.g. B_DN.

Emits src/clippr/scaffold.json. Run once; re-run only if the source tables change.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
SUP = ROOT / "data" / "grasp_supp"
OUT = ROOT / "src" / "clippr" / "scaffold.json"


def load_modules() -> dict[str, dict]:
    s1 = pd.read_excel(SUP / "Table S1.xlsx")
    mods = {}
    for _, r in s1.iterrows():
        pid = str(r["Plasmid ID"]).strip()
        key = pid.split("_", 1)[1] if "_" in pid else pid
        mods[key] = {
            "dna": str(r["Insert Sequence"]).strip().upper(),
            "oh5": str(r["5ʹ overhang"]).strip().upper(),
            "oh3": str(r["3ʹ overhang"]).strip().upper(),
        }
    return mods


def resolve(token: str, mods: dict) -> str:
    """Table S2 token -> Table S1 module key.

    Three naming shapes occur in the deposited set:
      B_DN   -> B_LD5N        both residues specified
      1A_N   -> 1A_5N_AATG    only a fifth residue, plus a fusion-site suffix
      2E_D   -> 2E_LD         only a last residue (terminal blocks carry no fifth)
    """
    block, _, code = token.partition("_")
    same_block = [k for k in mods if k.split("_")[0] == block]

    if len(code) == 2:                                   # last + fifth
        want = f"{block}_L{code[0]}5{code[1]}"
        if want in mods:
            return want

    if len(code) == 1:
        for pattern in (f"{block}_L{code}", f"{block}_5{code}"):
            exact = [k for k in same_block if k == pattern]
            if exact:
                return exact[0]
        prefixed = [k for k in same_block if k.startswith(f"{block}_5{code}_")]
        if prefixed:
            return sorted(prefixed)[0]
        loose = [k for k in same_block if f"5{code}" in k or f"L{code}" in k]
        if loose:
            return sorted(loose)[0]

    raise KeyError(f"cannot resolve {token!r}; block {block} has {sorted(same_block)}")


def assemble(tokens: list[str], mods: dict) -> str:
    """Golden Gate: adjacent modules share a 4-nt overhang, so trim the duplicate."""
    dna = ""
    for i, tok in enumerate(tokens):
        m = mods[resolve(tok, mods)]
        seq = m["dna"]
        if i and seq[:4] == prev_oh3:
            seq = seq[4:]
        dna += seq
        prev_oh3 = m["oh3"]
    return dna


def find_repeat(protein: str) -> tuple[str, int]:
    """Strongest periodicity, and the consensus repeat starting at the ARELF motif."""
    best = max(range(28, 40),
               key=lambda p: sum(protein[i] == protein[i + p]
                                 for i in range(len(protein) - p)) / (len(protein) - p))
    return best


def main() -> None:
    mods = load_modules()
    s2 = pd.read_excel(SUP / "Table S2.xlsx")
    wt = s2[s2.sPPR.astype(str).str.strip() == "p0"].iloc[0]

    cds1 = assemble(str(wt["CDS1 Modules"]).split(), mods)
    cds2 = assemble(str(wt["CDS2 Modules"]).split(), mods)

    prot = {}
    for name, dna in (("CDS1", cds1), ("CDS2", cds2)):
        for frame in range(3):
            trimmed = dna[frame: len(dna) - ((len(dna) - frame) % 3)]
            aa = str(Seq(trimmed).translate())
            if "*" not in aa[:-1] and "ARELF" in aa:
                prot[name] = aa
                print(f"{name}: frame {frame}, {len(aa)} aa")
                print(f"   {aa}")
                break
        else:
            raise SystemExit(f"{name}: no clean ARELF-containing frame")

    joined = prot["CDS1"] + prot["CDS2"]
    period = find_repeat(joined)
    print(f"\nrepeat periodicity: {period} aa")

    # the ARELF motif anchors each repeat; positions of every occurrence
    arelf = [m.start() for m in re.finditer("ARELF", joined)]
    print(f"ARELF occurrences at: {arelf}")
    spacings = [b - a for a, b in zip(arelf, arelf[1:])]
    print(f"spacings: {spacings}")

    data = {
        "provenance": (
            "Derived from Farley et al. 2025 Supplementary Tables S1 and S2 by assembling "
            "the deposited modules per the published recipe for variant p0. No third-party "
            "implementation was consulted."
        ),
        "code_to_base": {"TN": "A", "NN": "C", "TD": "G", "ND": "T"},
        "code_to_base_source": "Farley et al. 2025, Fig 1B: Thr/Asn->A, Asn/Asn->C, Thr/Asp->G, Asn/Asp->U",
        "cds1_protein": prot["CDS1"],
        "cds2_protein": prot["CDS2"],
        "repeat_period_aa": period,
        "arelf_positions_wt": arelf,
        "arelf_spacings": spacings,
        "architectures": {"9S": 9, "14S": 14, "19S": 19},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
