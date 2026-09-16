"""M3 — recompute every promised module invariant from the DNA, trusting nothing the pilot said.

The pilot calls the production optimiser and records whatever `constraints_ok` it returned.
That is the optimiser grading its own work, and the saved records do not contain the complete
before/after QC the M3 result document claims. This recomputes each invariant from the
original and recoded DNA directly.

Nothing here imports `clippr.qc` or `clippr.codons`. The constraint set, the genetic code
handling and the CAI definition are written out below, so agreement with the pilot is
corroboration rather than an echo.

**The invariants, one per promise M3 makes:**

  * length unchanged -- a recoded module must drop into the same assembly
  * both four-base interfaces byte-identical to the original
  * in-frame protein identical to the original's, in the module's own declared frame
  * no forbidden recognition site, either strand
  * GC within the declared band, in every window
  * no homopolymer run over the limit
  * CAI recomputed here, and the claimed improvement reproduced

A module failing any of these is a defect regardless of what the pilot recorded. All 21
intended modules stay in the denominator: a module the pilot skipped is a module not
delivered, not a module excluded from the count.

    python validation/experiments/m3_independent_qc.py
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from Bio.Data import CodonTable                                   # noqa: E402
from Bio.Seq import Seq                                           # noqa: E402

from clippr.parts import select                                   # noqa: E402
from clippr.products import load_inserts                          # noqa: E402

SIX = ["AAAAUGUGG", "GCUAAAGAC", "UUACACGUG", "CGUACGUAC", "AUCGAUCGA", "GGCCAAUUG"]
INTERFACE = 4

#: The declared contract, written out rather than imported from `policy`/`constants`.
BLACKLIST = {"BsaI": "GGTCTC", "BbsI": "GAAGAC", "SapI": "GCTCTTC"}
GC_BAND = (0.35, 0.65)
GC_WINDOW = 50
MAX_HOMOPOLYMER = 4
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def cai(dna: str, table: dict) -> float:
    """Geometric mean of relative adaptiveness over the codons that carry a choice.

    Methionine and tryptophan have one codon each, so including them would dilute the score
    with values that are 1.0 by construction and cannot be optimised.
    """
    forward = CodonTable.unambiguous_dna_by_id[1].forward_table
    weights = {c.replace("U", "T"): v / max(d.values())
               for _aa, d in table.items() for c, v in d.items() if max(d.values()) > 0}
    used = [weights[dna[i:i + 3]] for i in range(0, len(dna) - 2, 3)
            if dna[i:i + 3] in weights and forward.get(dna[i:i + 3]) not in ("M", "W", None)]
    if not used or min(used) <= 0:
        return 0.0
    return math.exp(sum(math.log(w) for w in used) / len(used))


def module_frames() -> dict[str, int]:
    """Each module's reading frame, derived here from the kit layout, not read from the pilot.

    The product's coding frame is 1, set by the AATG fusion site, so codon starts sit at
    product positions congruent to 1 mod 3. The first whole codon at or after a module's
    offset therefore begins `(1 - offset) % 3` bases into that module.
    """
    frames: dict[str, int] = {}
    sequences, _ = load_inserts(str(ROOT / "data" / "grasp_supp" / "Table S1.xlsx"))
    for target in SIX:
        offset = 0
        for module in select(target).modules:
            frame = (1 - offset) % 3
            prior = frames.get(module.plasmid_id)
            if prior is not None and prior != frame:
                raise SystemExit(f"{module.plasmid_id} occupies frames {prior} and {frame}")
            frames[module.plasmid_id] = frame
            offset += len(sequences[module.plasmid_id]) - INTERFACE
    return frames


def coding_span(seq: str, frame: int) -> str:
    n = (len(seq) - frame) // 3
    return seq[frame:frame + 3 * n]


def invariants(original: str, recoded: str, frame: int, table: dict) -> dict:
    """Every promise, checked. `problems` empty means the module is deliverable."""
    problems = []

    if len(recoded) != len(original):
        problems.append(f"length {len(recoded)} != original {len(original)}")
    if recoded[:INTERFACE] != original[:INTERFACE]:
        problems.append(f"5' interface {recoded[:INTERFACE]} != {original[:INTERFACE]}")
    if recoded[-INTERFACE:] != original[-INTERFACE:]:
        problems.append(f"3' interface {recoded[-INTERFACE:]} != {original[-INTERFACE:]}")

    old_coding, new_coding = coding_span(original, frame), coding_span(recoded, frame)
    old_protein = str(Seq(old_coding).translate())
    new_protein = str(Seq(new_coding).translate())
    if new_protein != old_protein:
        first = next((i for i, (a, b) in enumerate(zip(new_protein, old_protein)) if a != b),
                     min(len(new_protein), len(old_protein)))
        problems.append(f"protein differs at residue {first}")
    if "*" in new_protein:
        problems.append(f"stop codon at residue {new_protein.index('*')}")

    for name, site in BLACKLIST.items():
        for strand, pattern in (("forward", site), ("reverse", rc(site))):
            at = recoded.find(pattern)
            if at >= 0:
                problems.append(f"{name} site on the {strand} strand at {at}")

    for start in range(0, max(1, len(recoded) - GC_WINDOW + 1)):
        window = recoded[start:start + GC_WINDOW]
        if len(window) < GC_WINDOW:
            break
        gc = (window.count("G") + window.count("C")) / len(window)
        if not GC_BAND[0] <= gc <= GC_BAND[1]:
            problems.append(f"GC {gc:.3f} outside {GC_BAND} in the window at {start}")
            break

    run, prev = 1, ""
    for base in recoded:
        run = run + 1 if base == prev else 1
        prev = base
        if run > MAX_HOMOPOLYMER:
            problems.append(f"homopolymer run longer than {MAX_HOMOPOLYMER}")
            break

    return {"problems": problems,
            "cai_before": round(cai(old_coding, table), 4),
            "cai_after": round(cai(new_coding, table), 4),
            "protein_aa": len(new_protein),
            "coding_nt": len(new_coding)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pilot", default=str(ROOT / "work" / "m3" / "inventory_pilot.json"))
    ap.add_argument("--table", default=str(ROOT / "data" / "grasp_supp" / "Table S1.xlsx"))
    args = ap.parse_args()

    pilot = json.loads(Path(args.pilot).read_text(encoding="utf-8"))
    table = json.loads((ROOT / "data" / "codon_tables" / "kazusa_3055.json")
                       .read_text(encoding="utf-8"))
    originals, source = load_inserts(args.table)
    frames = module_frames()

    # Every module the six targets require, whether or not the pilot produced one.
    intended = sorted(frames)
    delivered = pilot["modules"]
    print(f"{len(intended)} modules required by the six reference targets; pilot delivered "
          f"{len(delivered)}")
    print(f"  source {source.sha256[:16]}  |  contract {sorted(BLACKLIST)}, GC {GC_BAND} "
          f"over {GC_WINDOW} nt, homopolymer <= {MAX_HOMOPOLYMER}\n")

    rows, failing, disagreements, missing = [], [], [], []
    for pid in intended:
        if pid not in delivered:
            missing.append(pid)
            continue
        record = delivered[pid]
        checked = invariants(originals[pid], record["sequence"], frames[pid], table)
        checked["module"] = pid
        checked["frame"] = frames[pid]
        checked["frame_agrees_with_pilot"] = frames[pid] == record["frame"]
        if not checked["frame_agrees_with_pilot"]:
            checked["problems"].append(
                f"frame {frames[pid]} derived here, pilot recorded {record['frame']}")

        # The pilot's own numbers, compared rather than trusted.
        claimed_gain = record["cai_after"] > record["cai_before"]
        actual_gain = checked["cai_after"] > checked["cai_before"]
        if claimed_gain != actual_gain:
            disagreements.append((pid, record["cai_after"], checked["cai_after"]))
        checked["pilot_cai_after"] = record["cai_after"]
        checked["cai_matches_pilot"] = abs(checked["cai_after"]
                                           - record["cai_after"]) < 5e-4

        rows.append(checked)
        if checked["problems"]:
            failing.append(pid)
        status = "PASS" if not checked["problems"] else "FAIL"
        print(f"  {pid:22s} frame {checked['frame']}  {checked['coding_nt']:>4} nt  "
              f"CAI {checked['cai_before']:.4f} -> {checked['cai_after']:.4f}  "
              f"{'=' if checked['cai_matches_pilot'] else '≠'}pilot  {status}")
        for problem in checked["problems"]:
            print(f"      {problem}")

    passing = len(rows) - len(failing)
    print(f"\nindependently PASS all invariants: {passing}/{len(intended)} intended modules")
    if missing:
        print(f"  not delivered by the pilot: {len(missing)} {missing}")
    print(f"  CAI reproduced within 5e-4 of the pilot: "
          f"{sum(1 for r in rows if r['cai_matches_pilot'])}/{len(rows)}")
    print(f"  improvement direction disagreements with the pilot: {len(disagreements)}")
    print(f"  frames derived here agree with the pilot's: "
          f"{sum(1 for r in rows if r['frame_agrees_with_pilot'])}/{len(rows)}")

    out = ROOT / "work" / "m3" / "independent_qc.json"
    out.write_text(json.dumps(
        {"intended_modules": len(intended), "delivered": len(rows),
         "not_delivered": missing, "passing": passing, "failing": failing,
         "disagreements": [{"module": m, "pilot": p, "independent": i}
                           for m, p, i in disagreements],
         "contract": {"blacklist": BLACKLIST, "gc_band": GC_BAND, "gc_window": GC_WINDOW,
                      "max_homopolymer": MAX_HOMOPOLYMER},
         "independence": ("constraint set, frames and CAI recomputed here; imports no "
                          "clippr.qc and no clippr.codons"),
         "modules": rows}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0 if passing == len(intended) else 1


if __name__ == "__main__":
    raise SystemExit(main())
