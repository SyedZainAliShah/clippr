"""V1 preflight — state exactly what CLIPPR's fidelity number is, before comparing it.

The orchestrator's point: a discrepancy against NEB's viewer is uninterpretable until both
the *data* and the *formula* are pinned on each side. CLIPPR implements one aggregation
convention (A) and documents a second (B) that can reverse rankings between sets, so
"we disagree with NEB" could mean a wrong table, a wrong condition, a different aggregation,
a different treatment of reverse complements, or an actual defect.

This writes the card that makes the comparison decidable: source workbook and sheet, table
orientation, counts, the exact numerator and denominator, the aggregation, whether
destination ends are included, and both conventions evaluated on the registered panel. It
compares nothing and needs no network. Running it before touching the viewer is what makes
a later disagreement diagnosable rather than merely alarming.

    python validation/independent/v1_scorer_card.py
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from clippr.overhangs import (SCORER_VERSION, _clean, _index,  # noqa: E402
                              fidelity_components, matrix_fingerprint,
                              reaction_overhangs, reverse_complement, set_fidelity,
                              valid_set)

MATRICES = {"BsaI-HFv2": "S1 Table. BsaI-HFv2", "BbsI-HF": "Table S4. BbsI-HF"}

#: Registered before evaluation. Grouped so the physical-model comparison is never mixed
#: with input-handling checks, which the viewer may reject outright rather than score.
PANEL = {
    "destination_pair_level0": ["CTCA", "CTCG"],
    "destination_pair_level_minus1": ["ACAT", "ACAA"],
    "level1_moclo_sites": ["GGAG", "AGCG"],
    "junctions_9S_selected": ["ATTC", "AATG", "AGCG"],
    "full_reaction_9S": None,                       # filled from reaction_overhangs below
    # Was ["GGAG","AGCG","AATG","CGCT"] until this card screened it: CGCT is the reverse
    # complement of AGCG, so that set is not a reaction -- it mixed a stored coding site
    # with an overhang, the exact confusion `reaction_overhangs()` exists to remove. It
    # scored 0.250 under A and 1.000 under B, which is what an undefined input looks like.
    "ordinary_valid_set": ["GGAG", "AGCG", "AATG", "TTCG"],
    "mismatch_prone_control": ["AAAG", "AAAC", "AAAT"],
    "aggregation_sensitive_A_over_B": ["GCAC", "GTGA", "TTGA", "CATA"],
    "aggregation_sensitive_B_over_A": ["GGTT", "AAAC", "AGGA", "ATCC"],
}
#: Separate from the physical panel: these probe input handling, not the model.
INPUT_HANDLING = {
    "palindrome": ["AATT", "GGAG"],
    "reverse_complement_pair": ["CTCA", "TGAG"],
    "duplicate_member": ["GGAG", "GGAG"],
}


def convention_b(overhangs, matrix: str) -> float:
    """The documented alternative: one product of per-junction ratios, both strands summed.

    Kept here rather than in the package because it is not what CLIPPR ships; it exists so
    the card can report both and so neither is chosen after seeing an external number.
    """
    ohs = _clean(overhangs)
    counts, idx = _index(matrix)
    pool = set(ohs) | {reverse_complement(o) for o in ohs}
    product = 1.0
    for o in ohs:
        numerator = counts[idx[o], idx[o]] + counts[idx[reverse_complement(o)],
                                                    idx[reverse_complement(o)]]
        denominator = sum(counts[idx[reverse_complement(q)], idx[o]]
                          + counts[idx[reverse_complement(q)], idx[reverse_complement(o)]]
                          for q in pool)
        if denominator == 0:
            return 0.0
        product *= numerator / denominator
    return product


def rounding_interval(decimal_places: int) -> float:
    """Half a unit in the last displayed place, on a 0-1 scale, for a percentage display."""
    return 0.5 * 10 ** (-decimal_places) / 100


def main() -> int:
    PANEL["full_reaction_9S"] = reaction_overhangs(["ATTC", "AATG", "AGCG"], "level0")

    card = {
        "purpose": "pin CLIPPR's fidelity model before any external comparison",
        "scorer_version": SCORER_VERSION,
        "data": {},
        "model": {
            "numerator": "counts[idx[O], idx[O]] -- O ligating its own Watson-Crick partner",
            "denominator": ("sum over every strand present of counts[idx[P], idx[O]] + "
                            "counts[idx[rc(P)], idx[O]], i.e. the set together with its "
                            "reverse complements"),
            "denominator_source": ("Pryor et al., Materials and Methods: N_total counts O's "
                                   "ligations 'to any overhangs in the set and its WC pair'"),
            "aggregation_implemented": ("A: geometric mean of the forward and reverse "
                                        "directional products"),
            "aggregation_alternative": ("B: single product of per-junction ratios with "
                                        "numerator and denominator each summed over both "
                                        "strands; reported here, not shipped"),
            "reverse_complements": ("pooled into the competitor set by the denominator; "
                                    "the score is invariant to reverse-complementing any "
                                    "member or the whole set"),
            "destination_ends": ("excluded unless the caller passes them; "
                                 "reaction_overhangs() adds them explicitly"),
            "competitor_multiplicity": "one per declared junction; no weighting applied",
        },
        "registered_panel": {},
        "input_handling": {},
        "rounding_intervals_on_0_1_scale": {
            f"{d}_decimal_percent": rounding_interval(d) for d in (0, 1, 2)},
        "not_established_by_agreement": [
            "observed assembly yield",
            "that any published experiment is defective",
            "that convention A is correct and B is not",
        ],
    }

    import openpyxl
    for matrix, sheet in MATRICES.items():
        path = ROOT / "src" / "clippr" / "data" / "matrices" / f"{matrix}.xlsx"
        counts, labels = _index(matrix)[0], None
        wb = openpyxl.load_workbook(path, read_only=True)
        card["data"][matrix] = {
            "workbook": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "sheet": wb.sheetnames[0],
            "expected_sheet": sheet,
            "shape": list(counts.shape),
            "total_counts": int(counts.sum()),
            "orientation": ("columns are the 256 overhangs in alphabetical order; rows are "
                            "their reverse complements, checked at load rather than assumed"),
            "content_fingerprint": matrix_fingerprint(matrix),
            "condition_metadata_in_file": "none beyond the sheet title",
            "condition_to_confirm_against_viewer": (
                "Pryor et al. 2020 supplementary S1 (BsaI-HFv2) / S4 (BbsI-HF); the ligation "
                "condition is not recorded in the workbook and must be read from the paper "
                "before selecting a dataset in the viewer"),
        }
        wb.close()

    invalid = {n: o for n, o in PANEL.items() if not valid_set(o)}
    if invalid:
        raise SystemExit(
            f"registered panel contains sets that are not reactions: {sorted(invalid)}. "
            f"Invalid inputs belong in INPUT_HANDLING, never in the physical-model panel -- "
            f"conventions diverge arbitrarily on them and the comparison means nothing.")

    for name, ohs in PANEL.items():
        entry = {"overhangs": list(ohs), "valid_set": True,
                 "includes_destination_ends": name == "full_reaction_9S"}
        for matrix in MATRICES:
            fwd, rev = fidelity_components(ohs, matrix)
            entry[matrix] = {
                "A_implemented": set_fidelity(ohs, matrix),
                "B_alternative": convention_b(ohs, matrix),
                "forward_product": fwd,
                "reverse_product": rev,
            }
            entry[matrix]["A_minus_B"] = (entry[matrix]["A_implemented"]
                                          - entry[matrix]["B_alternative"])
        card["registered_panel"][name] = entry

    for name, ohs in INPUT_HANDLING.items():
        try:
            value = set_fidelity(ohs, "BsaI-HFv2")
            card["input_handling"][name] = {"overhangs": ohs, "scored": value,
                                            "note": "scored; validity is a separate check"}
        except Exception as exc:                       # noqa: BLE001 - recorded as behaviour
            card["input_handling"][name] = {"overhangs": ohs,
                                            "raised": f"{type(exc).__name__}: {exc}"}

    out = ROOT / "work" / "independent" / "v1_scorer_card.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"scorer version {SCORER_VERSION}")
    for matrix, d in card["data"].items():
        print(f"  {matrix}: sheet {d['sheet']!r}  counts {d['total_counts']:,}  "
              f"fingerprint {d['content_fingerprint']}")
    print(f"\n{'registered set':32s} {'A (shipped)':>13} {'B (alt)':>13} {'A-B':>11}")
    for name, entry in card["registered_panel"].items():
        e = entry["BsaI-HFv2"]
        print(f"  {name:30s} {e['A_implemented']:13.9f} {e['B_alternative']:13.9f} "
              f"{e['A_minus_B']:+11.9f}")
    worst = max(abs(v["BsaI-HFv2"]["A_minus_B"]) for v in card["registered_panel"].values())
    print(f"\nlargest |A-B| on the registered panel: {worst:.9f}")
    print(f"rounding interval if the viewer shows 1 decimal percent: "
          f"{rounding_interval(1):.7f} on a 0-1 scale")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
