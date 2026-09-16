"""M3 pilot — recode the kit modules for one host while keeping every interface intact.

Bounded exactly as the orchestrator specified: the 21 modules the six reference targets
need, one declared expression context, at most four candidates per module, and a wall-clock
cap for the whole pilot. It does not touch the other 21 modules, does not search hosts, and
does not redesign interfaces.

**Context is the hard part, and it is measured rather than assumed.** A module's insert is
not a standalone coding sequence: it joins its neighbours through a shared four-base
overhang, and codons run across those joins. Two facts settle how to handle that, both
measured across all 200 corpus targets:

  * every target has a kit route, using 40 distinct modules;
  * **no module ever appears at more than one reading frame** -- 4 modules sit at frame 0,
    28 at frame 1, 8 at frame 2. Not because contributions are codon multiples (only 18 of
    42 are) but because the kit's architecture fixes each module's position class.
    (An earlier draft transposed the frame 1 and frame 2 counts here. The computation was
    always the one above; only the prose was wrong. Over all 42 modules, including the two
    AGGT start variants the corpus never selects, the distribution is 4 / 30 / 8.)

So each module has one frame, and can be recoded in it. And because the leading and trailing
partial codons are always inside the four-base interfaces, freezing those interfaces also
freezes every boundary-spanning base -- the neighbour's codons cannot be disturbed.

**What this can and cannot show.** A redesigned inventory can improve codon-table agreement
in a declared host while preserving proteins and interfaces. It cannot reduce DNA shared
between products that reuse the same records: recoding a shared module changes both products
identically. M3 is a reusable-inventory capability, not a route to library diversification.

    python validation/experiments/m3_inventory_pilot.py [--wall-seconds 600]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from Bio.Data import CodonTable                                   # noqa: E402
from Bio.Seq import Seq                                           # noqa: E402

from clippr import constants as C                                 # noqa: E402
from clippr.codons import complete_table, optimize_cds            # noqa: E402
from clippr.parts import select                                   # noqa: E402
from clippr.products import load_inserts                          # noqa: E402

SIX = ["AAAAUGUGG", "GCUAAAGAC", "UUACACGUG", "CGUACGUAC", "AUCGAUCGA", "GGCCAAUUG"]
INTERFACE = 4
TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"


def module_contexts(sequences: dict[str, str]) -> dict[str, dict]:
    """Each pilot module with the frame it occupies and the interfaces it must keep."""
    contexts: dict[str, dict] = {}
    for target in SIX:
        plan = select(target)
        offset = 0
        for module in plan.modules:
            seq = sequences[module.plasmid_id]
            # The product's coding frame is 1 (set by the AATG fusion site), so codon starts
            # sit at product positions congruent to 1 mod 3. The first one at or after this
            # module's offset is therefore (1 - offset) % 3 bases into the module. Getting
            # this sign backwards translates the wrong frame, which the stop-codon guard
            # below caught on 16 of 21 modules.
            frame = (1 - offset) % 3
            prior = contexts.get(module.plasmid_id)
            if prior and prior["frame"] != frame:
                raise SystemExit(f"{module.plasmid_id} appears at frames {prior['frame']} "
                                 f"and {frame}; the pilot assumes one frame per module")
            contexts[module.plasmid_id] = {
                "sequence": seq, "frame": frame, "length": len(seq),
                "five_interface": seq[:INTERFACE], "three_interface": seq[-INTERFACE:],
                "used_by": sorted(set((prior or {}).get("used_by", [])) | {target}),
            }
            offset += len(seq) - INTERFACE
    return contexts


def coding_span(context: dict) -> tuple[int, int]:
    """First and last index of the whole codons this module carries in its own frame."""
    start = context["frame"]
    n = (context["length"] - start) // 3
    return start, start + 3 * n


def cai(dna: str, table) -> float:
    forward = CodonTable.unambiguous_dna_by_id[1].forward_table
    weights = {c.replace("U", "T"): v / max(d.values())
               for aa, d in table.items() for c, v in d.items() if max(d.values()) > 0}
    used = [weights[dna[i:i + 3]] for i in range(0, len(dna) - 2, 3)
            if forward.get(dna[i:i + 3]) not in ("M", "W", None)]
    return math.exp(sum(math.log(w) for w in used) / len(used)) if used and min(used) > 0 else 0.0


def redesign(context: dict, table, seed: int, per_module_seconds: float) -> dict | None:
    """One recoded candidate, or None if the optimiser could not satisfy the constraints."""
    seq = context["sequence"]
    start, end = coding_span(context)
    protein = str(Seq(seq[start:end]).translate())
    if "*" in protein:
        return {"skipped": "module carries a stop codon in its own frame"}

    # Freeze both interfaces, but only the part of each that lies inside this module's
    # coding span -- bases outside it are carried over verbatim and never reach the
    # optimiser. Locking the whole four bases would address positions past the end of the
    # coding sequence being optimised.
    length = context["length"]
    locked = {}
    head_end = min(INTERFACE, end)
    if head_end > start:
        locked[0] = seq[start:head_end]
    tail_from = max(length - INTERFACE, start)
    if end > tail_from:
        locked[tail_from - start] = seq[tail_from:end]

    began = time.perf_counter()
    result = optimize_cds(protein, locked_sites=locked,
                          codon_table=table, genetic_code=1, seed=seed,
                          enzymes=C.DEFAULT_ENZYME_PROFILE, unique_kmer_size=None)
    elapsed = time.perf_counter() - began
    if elapsed > per_module_seconds:
        return {"timed_out": elapsed}
    return {"cds": result["cds"], "constraints_ok": bool(result["constraints_ok"]),
            "seconds": elapsed}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--wall-seconds", type=float, default=600.0)
    ap.add_argument("--per-module-seconds", type=float, default=10.0)
    ap.add_argument("--candidates", type=int, default=4)
    ap.add_argument("--table", default=str(ROOT / "data" / "grasp_supp" / "Table S1.xlsx"))
    args = ap.parse_args()

    sequences, source = load_inserts(args.table)
    table = complete_table(json.loads(TABLE.read_text(encoding="utf-8")), 1)
    contexts = module_contexts(sequences)
    print(f"pilot inventory: {len(contexts)} modules required by the six reference targets")
    print(f"  frames present: {sorted({c['frame'] for c in contexts.values()})}")
    print(f"  source: {source.sha256[:16]}  |  context: Chlamydomonas nuclear, code 1")
    print(f"  budget: {args.wall_seconds}s total, {args.per_module_seconds}s per module, "
          f"{args.candidates} candidates each\n")

    started = time.perf_counter()
    records, timed_out, unchanged = {}, [], []
    for pid in sorted(contexts):
        if time.perf_counter() - started >= args.wall_seconds:
            timed_out.append(pid)
            continue
        ctx = contexts[pid]
        start, end = coding_span(ctx)
        original = ctx["sequence"]
        before = cai(original[start:end], table)

        best = None
        for seed in range(42, 42 + args.candidates):
            got = redesign(ctx, table, seed, args.per_module_seconds)
            if not got or "cds" not in got or not got["constraints_ok"]:
                continue
            candidate = original[:start] + got["cds"] + original[end:]
            if candidate[:INTERFACE] != ctx["five_interface"] or \
               candidate[-INTERFACE:] != ctx["three_interface"]:
                continue
            if str(Seq(candidate[start:end]).translate()) != \
               str(Seq(original[start:end]).translate()):
                continue
            score = cai(candidate[start:end], table)
            if best is None or score > best["cai_after"]:
                best = {"sequence": candidate, "cai_after": score, "seed": seed,
                        "seconds": got["seconds"]}
        if best is None:
            unchanged.append(pid)
            continue
        improved = best["cai_after"] > before
        records[pid] = {
            "original_sha256": hashlib.sha256(original.encode()).hexdigest()[:16],
            "redesigned_sha256": hashlib.sha256(best["sequence"].encode()).hexdigest()[:16],
            "cai_before": round(before, 4), "cai_after": round(best["cai_after"], 4),
            "improved": improved, "identical": best["sequence"] == original,
            "frame": ctx["frame"], "used_by": ctx["used_by"], "seed": best["seed"],
            "sequence": best["sequence"] if improved else original,
        }
        flag = "improved" if improved else "kept original"
        print(f"  {pid:22s} frame {ctx['frame']}  CAI {before:.4f} -> "
              f"{best['cai_after']:.4f}  {flag}")

    elapsed = time.perf_counter() - started
    improved = [p for p, r in records.items() if r["improved"]]
    reused = [p for p, r in records.items() if len(r["used_by"]) > 1]
    print(f"\n{len(records)}/{len(contexts)} modules redesigned in {elapsed:.1f}s "
          f"(cap {args.wall_seconds}s)")
    print(f"  improved on CAI            : {len(improved)}")
    print(f"  kept original (no gain)    : {len(records) - len(improved)}")
    print(f"  no valid candidate         : {len(unchanged)} {unchanged[:4]}")
    print(f"  not attempted (budget)     : {len(timed_out)}")
    print(f"  reused by >1 target        : {len(reused)} of {len(records)}")

    out = ROOT / "work" / "m3"
    out.mkdir(parents=True, exist_ok=True)
    (out / "inventory_pilot.json").write_text(json.dumps(
        {"context": {"host": "Chlamydomonas nuclear", "genetic_code": 1,
                     "table": str(TABLE), "source_sha256": source.sha256},
         "budget": {"wall_seconds": args.wall_seconds,
                    "per_module_seconds": args.per_module_seconds,
                    "candidates": args.candidates, "elapsed_seconds": elapsed},
         "modules": records, "no_valid_candidate": unchanged,
         "not_attempted": timed_out,
         "limits": ["recoding shared modules cannot reduce DNA shared between products "
                    "that reuse them", "one host context; interfaces unchanged"]},
        indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {out / 'inventory_pilot.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
