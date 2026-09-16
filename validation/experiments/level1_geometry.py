"""Derive the level-1 reaction geometry from the deposited block plasmids, and falsify it.

`assembly_spec` reported a level-1 reaction unscorable until 2026-09-16, on the grounds that
"its overhangs are set by the block plasmid and its destination, which this package does not
have". That was wrong twice: the block plasmids were in the provenance directory the whole
time, and they turn out not to be needed.

**The derivation.** Each of the 28 assembled block plasmids in `9SDYW_level0.gb` carries
exactly two BsaI sites and no BbsI site -- the mirror of the 42 module plasmids, which carry
two BbsI and no BsaI. BsaI is `GGTCTC(1/5)`: a forward site at *i* leaves the 4-nt 5' overhang
at `[i+7, i+11)`, a reverse site at *i* leaves one at `[i-5, i-1)`. The released fragment runs
from the first base of its 5' overhang to the last base of its 3' overhang, carrying both --
the same span rule the level-0 derivation needed, and the same place that derivation first got
it wrong.

**What this establishes and what it does not.** It establishes the *block release geometry*.
It does **not** make the level-1 reaction scorable: that reaction co-assembles parts this
package does not compile, and ligation fidelity is a property of the whole competing overhang
set. Adding only the ends we know are omitted moves the 19S figure from 0.996046 to 0.742067,
below level 0 -- so a subset score cannot stand in for the reaction, and no bottleneck
conclusion follows from it.

**The falsifier, and the geometry it establishes.** The two fragment classes
are 499 nt `AGGT -> CTTC` and 406 nt `CTTC -> TTCG`. Joined at the shared overhang they are
901 nt -- *exactly* what this package compiles for 9S -- and the joined sequence occurs
verbatim inside the deposited level-1 product. So the block vector contributes no bases: the
released block fragment is the joined module inserts, and a level-1 reaction's ends are the
compiled product's own first and last interfaces, which Table S1 already supplies.

**Confirmed independently against Table S1**, which was never consulted by the derivation:
`1E -> CTTC`, `14E -> GTGA`, `19E -> CACG`, `2E -> TTCG`, with `1A` entering at `AATG` or
`AGGT` according to the variant.

**The participants that are missing.** A P2L2S2 linker (`TGTG -> CAAC`) and the consensus DYW
domain (`CAAC -> GCTT`) are in the same BsaI reaction, and at least one further part bridging
`TTCG` to `TGTG` is absent from the supplied records entirely. Until that list is complete the
reaction stays unscored, and the subset figure is reported as a named diagnostic.

    python validation/experiments/level1_geometry.py
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

#: Primary GenBank records, read for verification and never redistributed (`NOTICE.md`).
#: Set CLIPPR_GRASP_PRIMARY to point elsewhere; an absolute path here would pin this script
#: to one machine.
PRIMARY = Path(os.environ.get("CLIPPR_GRASP_PRIMARY")
               or ROOT.parent / "work" / "provenance" / "grasp_primary")
LEVEL0_BLOCKS = PRIMARY / "9SDYW_level0.gb"
LEVEL1_PRODUCTS = PRIMARY / "9SrpoaDYW_variants_level1.gb"
TABLE_S1 = ROOT / "data" / "grasp_supp" / "Table S1.xlsx"


def is_block(record) -> bool:
    """An assembled block plasmid, as opposed to a single part deposited alongside them."""
    return (record.description or "").lower().startswith(("pagm9121", "ppr"))


def bsai_release(seq: str) -> tuple[str, str, str]:
    """(fragment, 5' overhang, 3' overhang) from a circular plasmid with one site each way."""
    s = seq.upper()
    n = len(s)
    doubled = s + s
    forward = [i for i in range(n) if doubled.startswith("GGTCTC", i)]
    reverse = [i for i in range(n) if doubled.startswith("GAGACC", i)]
    if len(forward) != 1 or len(reverse) != 1:
        raise ValueError(f"expected one BsaI site each way, found "
                         f"{len(forward)} forward and {len(reverse)} reverse")
    start = forward[0] + 7
    end = reverse[0] - 1
    if end <= start:
        end += n                                  # the fragment crosses the origin
    fragment = doubled[start:end]
    return fragment, fragment[:4], fragment[-4:]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "work" / "level1" / "geometry.json"))
    args = ap.parse_args()

    if not LEVEL0_BLOCKS.is_file():
        print(f"missing {LEVEL0_BLOCKS}; this derivation needs the primary GenBank records")
        return 2

    from Bio import SeqIO

    records = list(SeqIO.parse(LEVEL0_BLOCKS, "genbank"))
    blocks = [r for r in records if is_block(r)]
    parts = [r for r in records if not is_block(r)]
    print(f"{len(records)} records in {LEVEL0_BLOCKS.name}: "
          f"{len(blocks)} assembled blocks, {len(parts)} single parts\n")

    classes: dict[tuple[str, str], list[str]] = collections.defaultdict(list)
    fragments: dict[str, str] = {}
    cuts_ok = True
    for record in blocks:
        try:
            fragment, five, three = bsai_release(str(record.seq))
        except ValueError as exc:
            print(f"  {record.name}: {exc}")
            cuts_ok = False
            continue
        classes[(five, three)].append(record.name)
        fragments[record.name] = fragment

    print("released block fragments, by exposed ends:")
    sizes = {}
    for (five, three), names in sorted(classes.items()):
        lengths = sorted({len(fragments[n]) for n in names})
        sizes[f"{five}->{three}"] = lengths
        print(f"   {five} -> {three}   n={len(names):2d}   fragment {lengths} nt")

    fives = {f for f, _ in classes}
    threes = {t for _, t in classes}
    entry, exit_ = sorted(fives - threes), sorted(threes - fives)
    internal = sorted(fives & threes)
    print(f"\n   entry {entry}   internal {internal}   exit {exit_}")

    # --- falsifier 1: joined blocks must reproduce a deposited level-1 product ---
    joined_length = None
    verbatim = []
    if LEVEL1_PRODUCTS.is_file() and internal:
        products = list(SeqIO.parse(LEVEL1_PRODUCTS, "genbank"))
        left = [n for n in fragments if fragments[n].endswith(internal[0])]
        right = [n for n in fragments if fragments[n].startswith(internal[0])]
        for a in left[:4]:
            for b in right[:4]:
                joined = fragments[a] + fragments[b][4:]
                joined_length = len(joined)
                for product in products:
                    doubled = (str(product.seq) + str(product.seq)).upper()
                    if joined in doubled:
                        verbatim.append({"blocks": [a, b], "nt": len(joined),
                                         "found_in": product.name})
                        break
                if verbatim:
                    break
            if verbatim:
                break
        print(f"\n   joined block pair: {joined_length} nt")
        for hit in verbatim:
            print(f"   found verbatim in {hit['found_in']}")
        if not verbatim:
            print("   NOT found in any level-1 product -- the derivation is wrong")

    # --- falsifier 2: the module table must say the same thing, independently ---
    table_agrees, table_ends = None, {}
    if TABLE_S1.is_file():
        from clippr import inventories as inv

        deposited = inv.load_deposited(TABLE_S1)
        for record in deposited.modules.values():
            if record.block in ("1E", "14E", "19E", "2E"):
                table_ends.setdefault(record.block, set()).add(record.three_interface)
            if record.block == "1A":
                table_ends.setdefault("1A entry", set()).add(record.five_interface)
        table_ends = {k: sorted(v) for k, v in sorted(table_ends.items())}
        print("\n   Table S1, consulted only now:")
        for block, ends in table_ends.items():
            print(f"      {block:10s} {ends}")
        table_agrees = (table_ends.get("1E") == ["CTTC"]
                        and table_ends.get("14E") == ["GTGA"]
                        and table_ends.get("19E") == ["CACG"]
                        and table_ends.get("2E") == ["TTCG"])
        print(f"   agrees with the block-plasmid derivation: {table_agrees}")

    # --- what it is worth: level-1 versus level-0 fidelity ---
    from clippr.overhangs import set_fidelity

    def score(overhangs, matrix):
        got = set_fidelity(list(overhangs), matrix=matrix)
        return round(got["fidelity"] if isinstance(got, dict) else got, 6)

    # The PPR block subset of the level-1 reaction. NOT the reaction: the deposited BsaI
    # reaction co-assembles parts this package does not compile, and ligation fidelity is a
    # property of the whole competing set.
    level1 = {
        "9S": ["AATG", "CTTC", "TTCG"],
        "14S": ["AATG", "CTTC", "GTGA", "TTCG"],
        "19S": ["AATG", "CTTC", "GTGA", "CACG", "TTCG"],
    }
    known_extra = ["TGTG", "CAAC", "GCTT"]      # linker + DYW, from this same primary file
    level0 = ["CTCA", "ACTC", "AAGA", "GCAC", "TGAA", "CGAG"]
    scores = {name: score(s, "BsaI-HFv2") for name, s in level1.items()}
    with_known = {name: score(s + known_extra, "BsaI-HFv2") for name, s in level1.items()}
    level0_score = score(level0, "BbsI-HF")
    print(chr(10) + "   level 0, participants complete (BbsI-HF): "
          f"{level0_score:.6f}  -- a reaction fidelity")
    print("   level-1 PPR block subsets (BsaI-HFv2) -- diagnostics, not reaction "
          "fidelities:")
    for name, value in scores.items():
        print(f"      {name:4s} subset {value:.6f}   subset + known omitted ends "
              f"{with_known[name]:.6f}")
    print(chr(10) + "   Adding only the ends we KNOW are omitted moves 19S from "
          f"{scores[chr(39)+chr(39)] if False else scores['19S']:.6f} to "
          f"{with_known['19S']:.6f},")
    print(f"   below level 0's {level0_score:.6f}. No conclusion about which stage is the "
          "bottleneck")
    print("   follows from these numbers, and the true participant list is still unknown.")

    passed = bool(cuts_ok and verbatim and table_agrees)
    print(f"\nLEVEL-1 GEOMETRY: {'DERIVED AND FALSIFIED' if passed else 'NOT ESTABLISHED'}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "source": str(LEVEL0_BLOCKS),
        "blocks": len(blocks), "single_parts": len(parts),
        "fragment_classes": sizes,
        "entry": entry, "internal": internal, "exit": exit_,
        "joined_nt": joined_length,
        "found_verbatim_in_level1_product": verbatim,
        "table_s1_interfaces": table_ends,
        "table_s1_agrees": table_agrees,
        "level0_fidelity_bbsi_hf": level0_score,
        "level1_block_subset_diagnostic": {
            "scores": scores,
            "plus_known_omitted_ends": with_known,
            "known_omitted_ends": known_extra,
            "is": ("a score of the PPR block subset, not a reaction fidelity: the deposited "
                   "BsaI reaction co-assembles parts this package does not compile"),
        },
        "not_covered": ("co-assembled linker TGTG->CAAC and DYW domain CAAC->GCTT are in the "
                        "same BsaI reaction and are not compiled here; at least one further "
                        "part bridging TTCG to TGTG is absent from the supplied records"),
        "result": "derived and falsified" if passed else "not established",
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
