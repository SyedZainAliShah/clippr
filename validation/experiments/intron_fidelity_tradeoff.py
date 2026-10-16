"""Best achievable fidelity against the number of introns, per architecture and level.

The decision-relevant shape. Exhaustive search over (site x overhang) is intractable for
19S -- choosing 7 positions from ~90 candidates is billions of combinations -- so positions
are fixed by even spacing and `overhangs.best_set` then solves the overhang choice, which is
the part it was written for and provably optimal over that space.
"""
import json, sys, time
from pathlib import Path
from ppropt import DESTINATION_OVERHANGS, design_oneshot, introns, substrates as sub
from ppropt.arelf import safe_overhangs
from ppropt.assembly import reverse_complement, wrap_fragment, wrapper_intended_sites
from ppropt.overhangs import best_set
from ppropt.synthesis_profile import STRICT_LEGACY as P

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
RES = OUT / "tradeoff.json"
rows = json.loads(RES.read_text(encoding="utf-8")) if RES.is_file() else []
done = {r["arch"] for r in rows}
MIN_EX, MAX_EX = introns.MIN_EXON_NT, introns.MAX_EXON_NT
I = introns.RBCS2I1


def options(cds, prot, s):
    """Overhangs usable at codon-aligned nt position s: NG/GN, and a clean intron part."""
    if s % 3 or s < MIN_EX or s > len(cds) - MIN_EX or cds[s] != "G":
        return []
    out = []
    for o in safe_overhangs(prot, s // 3, "igem_rfc1000"):
        if not o.endswith("G"):
            continue
        pay = o + I + cds[s:s + 4]
        w = wrap_fragment(pay, introns.PART_ENZYME)
        f = sub.synthesis_findings(w, P, intended_sites=wrapper_intended_sites(pay, introns.PART_ENZYME),
                                   enzymes=("BsaI", "BbsI", "SapI"))
        if not f["hard"]:
            out.append(o)
    return out


for target, arch in (("AAAAUGUGG", "9S"), ("UUACACGUGCGUAC", "14S"),
                     ("CUAUCACAUCACAUAAGCG", "19S")):
    if arch in done:
        continue
    t0 = time.time()
    r = design_oneshot(target, organism="c_reinhardtii_nuclear", codon_table=None,
                       genetic_code=None, enzyme_profile="igem_rfc1000", enzyme="BsaI",
                       matrix="BsaI-HFv2", destination=DESTINATION_OVERHANGS["level0"],
                       seed=42, check_offtarget=False, outdir=None)
    cds, prot = r["cds"], r["protein"]
    usable = {s: o for s in range(0, len(cds), 3) if (o := options(cds, prot, s))}

    curve = []
    for k in range(1, 9):
        if (k + 1) * MIN_EX > len(cds):
            break
        # Even spacing, then the nearest usable site to each ideal position.
        ideals = [round((i + 1) * len(cds) / (k + 1)) for i in range(k)]
        picked, taken = [], set()
        for want in ideals:
            cand = sorted((abs(s - want), s) for s in usable if s not in taken)
            if not cand:
                break
            taken.add(cand[0][1]); picked.append(cand[0][1])
        picked.sort()
        if len(picked) != k:
            continue
        b = [0] + picked + [len(cds)]
        if not all(MIN_EX <= y - x <= MAX_EX for x, y in zip(b, b[1:])):
            continue
        entry = {"k": k, "positions": picked,
                 "exons": [y - x for x, y in zip(b, b[1:])]}
        for level in ("level0", "level1"):
            d = DESTINATION_OVERHANGS[level]
            try:
                chosen, fid = best_set([usable[s] for s in picked],
                                       fixed=[d[0], reverse_complement(d[1])],
                                       matrix="BsaI-HFv2")
                entry[level] = {"fidelity": round(fid, 4), "overhangs": chosen}
            except Exception as exc:
                entry[level] = {"error": f"{type(exc).__name__}"}
        curve.append(entry)

    row = {"arch": arch, "cds_nt": len(cds), "usable_sites": len(usable),
           "baseline_fidelity": round(r["fidelity"], 4), "curve": curve,
           "seconds": round(time.time() - t0, 1)}
    rows.append(row)
    RES.write_text(json.dumps(rows, indent=1), encoding="utf-8")
    print(f"{arch}  CDS {len(cds)} nt, {len(usable)} usable sites, baseline {row['baseline_fidelity']}")
    for e in curve:
        l0 = e["level0"].get("fidelity", "-"); l1 = e["level1"].get("fidelity", "-")
        print(f"   {e['k']} intron(s): level0 {l0}  level1 {l1}   exons {e['exons']}")
    print(flush=True)
print("done")
