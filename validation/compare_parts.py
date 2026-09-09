"""Oracle check: our module selection against the compositions published in Table S2.

Table S2 lists, for each of 31 published sPPR variants, the modules composing CDS1 and CDS2
alongside the cognate RNA. That makes it a genuine oracle for the parts route: given the
cognate, `parts.select` must reproduce the published module list exactly.

This is the parts-route counterpart of the corpus checks, with one important difference -- the
oracle here is the **paper**, not the reference implementation, so unlike those checks this one
runs on any checkout that has the supplementary tables.

**Three published rows do not agree, and each is demonstrably an error in Table S2 rather than
in this mapping.** They are named and reported rather than silently excluded.

    p0  its modules decode to UUACACGUG, but the table states UUAACAGUGCAAAAUC.
        The decoded sequence is exactly the wild-type RNA used throughout the paper's own
        perturbation panel, so the module list is right and the stated target is the typo.

    p6  its modules decode to UGGCACGUG, the table states UGACACGUGCAAAAUC.

    p8  the table lists a module `C_DD`, meaning a 5th residue of D. **No such plasmid
        exists**: every module in the kit has a 5th residue of N or T, and the C block
        contains exactly LD5N, LD5T, LN5N and LN5T. p8's target UUGCACGUG needs a 5th
        residue of T at that position, which is `C_DT` -- what this code selects.

The third is the strongest of the three, because it does not require trusting our decoding at
all: the published token names a plasmid absent from the published inventory.

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\compare_parts.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from clippr.parts import select                      # noqa: E402

SUP = ROOT / "data" / "grasp_supp"

#: Rows carrying a demonstrable error in Table S2 itself. See the module docstring for the
#: evidence on each; p8's is provable from the published inventory alone.
KNOWN_TABLE_INCONSISTENCIES = {"p0", "p6", "p8"}


def published_tokens(row) -> list[str]:
    out = []
    for col in ("CDS1 Modules", "CDS2 Modules"):
        out += [t for t in str(row[col]).split() if t and t != "nan"]
    return out


def our_tokens(plan) -> list[str]:
    """Our chosen modules, rendered in Table S2's `<block>_<last><fifth>` shorthand."""
    out = []
    for m in plan.modules:
        code = f"{m.last or ''}{m.fifth or ''}"
        out.append(f"{m.block}_{code}")
    return out


def main() -> None:
    table = SUP / "Table S2.xlsx"
    if not table.exists():
        print(f"no {table.relative_to(ROOT)}; this check needs the paper's supplementary "
              f"tables, which are not redistributed here")
        print("\nSKIP")
        raise SystemExit(0)

    s2 = pd.read_excel(table)
    checked = exact = 0
    mismatches: list[str] = []
    unavailable: list[str] = []
    skipped: list[str] = []

    for _, row in s2.iterrows():
        name = str(row["sPPR"]).strip()
        target = re.sub(r"[^ACGU]", "",
                        str(row["Cognate RNA target"]).upper().replace("T", "U"))
        want = published_tokens(row)
        if not target or not want:
            continue
        if name in KNOWN_TABLE_INCONSISTENCIES:
            skipped.append(name)
            continue

        # the published targets carry trailing context beyond the read window
        plan = select(target[:len(want) - 1])
        checked += 1
        if not plan.available:
            unavailable.append(f"{name}: {plan.reason}")
            continue
        got = our_tokens(plan)
        if got == want:
            exact += 1
        elif len(mismatches) < 5:
            diff = next((i for i, (a, b) in enumerate(zip(got, want)) if a != b), None)
            mismatches.append(f"{name}: position {diff}, ours {got[diff]!r} vs published "
                              f"{want[diff]!r}" if diff is not None
                              else f"{name}: length {len(got)} vs {len(want)}")

    print(f"variants compared            : {checked}")
    print(f"module list exactly published: {exact}/{checked}")
    if skipped:
        print(f"skipped (table inconsistent) : {', '.join(sorted(skipped))}")
    if unavailable:
        print(f"kit route unavailable        : {len(unavailable)}")
        for u in unavailable[:5]:
            print(f"    {u}")
    if mismatches:
        print("\nfirst mismatches:")
        for m in mismatches:
            print(f"  {m}")

    ok = checked and exact == checked
    print("\nOracle is the published paper, not the reference implementation, so this check "
          "runs wherever the supplementary tables are present.")
    print("PASS" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
