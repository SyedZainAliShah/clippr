"""V6 — check CLIPPR's own designs against the declared hard constraints, independently.

The existing codon comparison re-optimises under the *reference's* cut positions and locked
overhangs, so it never exercises the plans CLIPPR chooses for itself. This checks those
plans, and it derives every expectation from the declared contract rather than from the
design result or its QC verdict -- a checker that reads `constraints_ok` and agrees with it
has measured nothing.

Nothing here imports `clippr.qc` or `clippr.codons`. Enzyme recognition sites, the GC band,
the window size and the homopolymer limit are written out below.

**One stated dependency.** The full-protein check rebuilds the expected sequence from
`biology.N_TERMINAL` and `REPEAT_TEMPLATE`, which are production constants. They are derived
from the paper and were checked against the authors' deposited construct in V3, so this is
stronger than the decode-only check it replaced -- but it is not a wholly independent source,
and it is not `describe()`, which is the function under test. Decoding the specificity
residues alone accepted a Q->F change at residue 1.

**Hard constraints** (a violation is a defect) are kept apart from **objectives** (CAI,
repeat reduction), which are reported but never used to pass or fail a design.

    python validation/independent/v6_constraints.py [--limit N]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from journal import load_designs, provenance                       # noqa: E402

#: Declared in `policy.ENZYME_PROFILES` as the igem_rfc1000 profile. Written out rather than
#: imported so this is a second statement of the contract, not a restatement of the code.
BLACKLIST = {"BsaI": "GGTCTC", "BbsI": "GAAGAC", "SapI": "GCTCTTC"}
GC_BAND = (0.35, 0.65)
GC_WINDOW = 50
MAX_HOMOPOLYMER = 4
#: The published PPR code, transcribed from the paper.
CODE = {("T", "N"): "A", ("N", "N"): "C", ("T", "D"): "G", ("N", "D"): "U"}
REPEAT_LEN = 31
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

TRANSLATION = {}


def _codon_table():
    if not TRANSLATION:
        from Bio.Data import CodonTable
        TRANSLATION.update(CodonTable.unambiguous_dna_by_id[1].forward_table)
        for stop in CodonTable.unambiguous_dna_by_id[1].stop_codons:
            TRANSLATION[stop] = "*"
    return TRANSLATION


def translate(dna: str) -> str:
    table = _codon_table()
    return "".join(table.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


def rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def decoded_target(protein: str) -> str | None:
    """Read the target back out of the protein using the published code.

    Independent of `describe()`: it locates the repeats by their invariant spacing from the
    end of the N-terminal region and reads the fifth and last residue of each.
    """
    anchor = protein.find("GVVS")
    if anchor < 0:
        return None
    base = anchor + 4
    n = (len(protein) - base) // REPEAT_LEN
    if n < 1 or (len(protein) - base) % REPEAT_LEN:
        return None
    pairs = [(protein[base + REPEAT_LEN * i + 1], protein[base + REPEAT_LEN * i + 27])
             for i in range(n)]
    return "".join(CODE.get(p, "?") for p in pairs)


def expected_protein(target: str) -> str | None:
    """Rebuild the protein a target should give, from the template and the published code.

    Uses the scaffold constants, which derive from the paper and were checked against the
    authors' deposited construct in V3 -- not from `describe()`, which is the function whose
    output this is meant to test.
    """
    from clippr.biology import N_TERMINAL, REPEAT_TEMPLATE

    inverse = {base: pair for pair, base in CODE.items()}
    repeats = []
    for base in target.upper().replace("T", "U"):
        if base not in inverse:
            return None
        fifth, last = inverse[base]
        repeats.append(REPEAT_TEMPLATE.format(fifth=fifth, last=last))
    return N_TERMINAL + "".join(repeats)


def violations(row: dict) -> list[str]:
    """Every hard-constraint breach in one design, named."""
    cds, target = row["cds"], row["target"]
    out = []

    if len(cds) % 3:
        out.append(f"CDS length {len(cds)} is not a whole number of codons")
    protein = translate(cds)
    if "*" in protein[:-1]:
        out.append(f"internal stop codon at residue {protein.index('*')}")
    if "X" in protein:
        out.append("untranslatable codon")

    read_back = decoded_target(protein)
    if read_back != target:
        out.append(f"protein decodes to {read_back!r}, not the requested {target!r}")

    # Decoding reads only the fifth and last residue of each repeat, so a change anywhere
    # else in the scaffold passed silently -- Q->F at residue 1 was accepted. Rebuild the
    # whole expected protein from the template and compare it in full.
    expected = expected_protein(target)
    if expected is not None and protein != expected:
        first = next((i for i, (a, b) in enumerate(zip(protein, expected)) if a != b),
                     min(len(protein), len(expected)))
        out.append(f"protein differs from the template-built expectation at residue {first} "
                   f"({protein[first:first + 1]!r} vs {expected[first:first + 1]!r})")

    if len(row["junction_overhangs"]) != len(row["cuts"]):
        out.append(f"{len(row['junction_overhangs'])} junction overhangs recorded for "
                   f"{len(row['cuts'])} cuts; absent metadata cannot be checked")
    for cut, overhang in zip(row["cuts"], row["junction_overhangs"]):
        found = cds[3 * cut - 4:3 * cut]
        if found != overhang:
            out.append(f"locked overhang {overhang} absent at cut {cut} (found {found})")

    for name, site in BLACKLIST.items():
        for strand, pattern in (("forward", site), ("reverse", rc(site))):
            at = cds.find(pattern)
            if at >= 0:
                out.append(f"{name} site on the {strand} strand at {at}")

    for start in range(0, max(1, len(cds) - GC_WINDOW + 1)):
        window = cds[start:start + GC_WINDOW]
        if len(window) < GC_WINDOW:
            break
        gc = (window.count("G") + window.count("C")) / len(window)
        if not GC_BAND[0] <= gc <= GC_BAND[1]:
            out.append(f"GC {gc:.3f} outside {GC_BAND} in the window at {start}")
            break

    run, prev = 1, ""
    for base in cds:
        run = run + 1 if base == prev else 1
        prev = base
        if run > MAX_HOMOPOLYMER:
            out.append(f"homopolymer run longer than {MAX_HOMOPOLYMER}")
            break
    return out


def objectives(row: dict) -> dict:
    """Reported, never used to pass or fail. These are preferences, not contract terms."""
    cds = row["cds"]
    windows = [cds[i:i + 20] for i in range(len(cds) - 19)]
    return {"duplicated_20mers": len(windows) - len(set(windows)),
            "gc_overall": round((cds.count("G") + cds.count("C")) / len(cds), 4)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jsonl", default=str(ROOT / "work" / "independent"
                                        / "corpus200_pinned.jsonl"))
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    rows = load_designs(args.jsonl)
    if args.limit:
        rows = rows[:args.limit]

    print(f"{len(rows)} designs, checked against the contract rather than their own verdict")
    print(f"  blacklist {sorted(BLACKLIST)} on both strands | GC {GC_BAND} over {GC_WINDOW} nt"
          f" | homopolymer <= {MAX_HOMOPOLYMER}")

    clean, offenders, disagreements, obj = 0, [], [], []
    for row in rows:
        bad = violations(row)
        obj.append(objectives(row))
        if bad:
            offenders.append((row["target"], bad))
        else:
            clean += 1
        # The design's own claim, compared with ours -- reported, not used as the expectation
        claimed = bool(row.get("constraints_ok")) and row.get("qc", {}).get("status") == "PASS"
        if claimed == bool(bad):
            disagreements.append((row["target"], claimed, bad))

    print(f"\nsatisfying every hard constraint: {clean}/{len(rows)}")
    for target, bad in offenders[:10]:
        print(f"   {target}: {'; '.join(bad[:3])}")
    print(f"disagreements with the design's own constraints_ok/QC: {len(disagreements)}")
    for target, claimed, bad in disagreements[:5]:
        print(f"   {target}: claimed_ok={claimed} but independent check found {bad[:2]}")

    dup = [o["duplicated_20mers"] for o in obj]
    gc = sorted(o["gc_overall"] for o in obj)
    print(f"\nobjectives (reported, not gating): designs with zero duplicated 20-mers "
          f"{sum(1 for d in dup if d == 0)}/{len(dup)}; GC min {gc[0]} median "
          f"{gc[len(gc) // 2]} max {gc[-1]}")

    out = ROOT / "work" / "independent" / "v6_constraints.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(
        {"provenance": provenance(__file__, args.jsonl),
         "designs": len(rows), "clean": clean,
         "offenders": [{"target": t, "violations": v} for t, v in offenders],
         "disagreements_with_claimed_status": len(disagreements),
         "contract": {"blacklist": BLACKLIST, "gc_band": GC_BAND, "gc_window": GC_WINDOW,
                      "max_homopolymer": MAX_HOMOPOLYMER},
         "independence": ("expectations derived from the declared contract and the published "
                          "PPR code; imports no clippr.qc or clippr.codons")},
        indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0 if clean == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
