"""Phase B gate — all 42 modules load, recode, compile and reload, checked independently.

Gate B, from the full-workflow plan:

  * all 42 records load and have a validated context
  * every module has a valid unchanged or recoded realisation
  * every tested product preserves its intended protein and interfaces
  * reload and recompile reproduce the same sequences
  * missing, swapped or incompatible records fail

Adaptation gains are reported **separately** from the capability pass. A gate that mixed them
would let a large CAI improvement paper over a module that never compiled.

Nothing here imports `clippr.qc` or `clippr.codons`. The constraint set, the frames and the
product join are written out or recomputed, so agreement with the production path is
corroboration rather than an echo of it.

    python validation/experiments/phaseb_gate.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from Bio.Seq import Seq                                           # noqa: E402

from clippr import inventories as inv                             # noqa: E402
from clippr.parts import FUSION_SITES, inventory as kit_index     # noqa: E402

INTERFACE = 4
BLACKLIST = {"BsaI": "GGTCTC", "BbsI": "GAAGAC", "SapI": "GCTCTTC"}
GC_BAND = (0.35, 0.65)
GC_WINDOW = 50
MAX_HOMOPOLYMER = 4
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

#: One target of each buildable architecture, plus the witness targets that reach the two
#: AGGT start variants the 200-target corpus never selects.
ARCHITECTURES = ["AAAAUGUGG", "UUACACGUGCGUAC", "CUAUCACAUCACAUAAGCG"]


def rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def module_problems(original, recoded) -> list[str]:
    """Every invariant a recoded module must hold, recomputed from the DNA."""
    problems = []
    if len(recoded.dna) != len(original.dna):
        problems.append(f"length {len(recoded.dna)} != {len(original.dna)}")
    if recoded.dna[:INTERFACE] != original.dna[:INTERFACE]:
        problems.append("5' interface changed")
    if recoded.dna[-INTERFACE:] != original.dna[-INTERFACE:]:
        problems.append("3' interface changed")
    if recoded.frame != original.frame:
        problems.append(f"frame {recoded.frame} != {original.frame}")

    start, end = original.coding_interval
    was = str(Seq(original.dna[start:end]).translate())
    now = str(Seq(recoded.dna[start:end]).translate())
    if now != was:
        problems.append("protein changed")
    if "*" in now:
        problems.append("stop codon in frame")

    for name, site in BLACKLIST.items():
        for strand, pattern in (("forward", site), ("reverse", rc(site))):
            at = recoded.dna.find(pattern)
            if at >= 0:
                problems.append(f"{name} {strand} site at {at}")

    for i in range(0, max(1, len(recoded.dna) - GC_WINDOW + 1)):
        window = recoded.dna[i:i + GC_WINDOW]
        if len(window) < GC_WINDOW:
            break
        gc = (window.count("G") + window.count("C")) / len(window)
        if not GC_BAND[0] <= gc <= GC_BAND[1]:
            problems.append(f"GC {gc:.3f} outside {GC_BAND} at {i}")
            break

    run, prev = 1, ""
    for base in recoded.dna:
        run = run + 1 if base == prev else 1
        prev = base
        if run > MAX_HOMOPOLYMER:
            problems.append(f"homopolymer run > {MAX_HOMOPOLYMER}")
            break
    return problems


def join_independently(library, module_ids) -> tuple[str, list[str]]:
    """Splice at the four-base overlaps here, not through `compile_target`."""
    problems = []
    product = library.modules[module_ids[0]].dna
    for left, right in zip(module_ids, module_ids[1:]):
        a, b = library.modules[left].dna, library.modules[right].dna
        if a[-INTERFACE:] != b[:INTERFACE]:
            problems.append(f"{left}/{right} overlap {a[-INTERFACE:]} vs {b[:INTERFACE]}")
        product += b[INTERFACE:]
    return product, problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--deposited", default=str(ROOT / "data" / "grasp_supp" / "Table S1.xlsx"))
    ap.add_argument("--recoded", default=str(ROOT / "work" / "phaseb"
                                             / "inventory_recoded.json"))
    ap.add_argument("--out", default=str(ROOT / "work" / "phaseb" / "gate_b.json"),
                    help="JSON file to write (not a directory)")
    ap.add_argument("--scratch", default=None,
                    help="where temporary probe files go; defaults beside --out, so an "
                         "independent review never writes into the executor's evidence")
    args = ap.parse_args()

    deposited = inv.load_deposited(args.deposited)
    recoded = inv.load(args.recoded)
    indexed = {m.plasmid_id for m in kit_index()}
    print(f"kit index {len(indexed)} | deposited {len(deposited)} | "
          f"recoded {len(recoded)}")

    # --- B1: every indexed module is present in both inventories ---
    absent_dep = sorted(indexed - set(deposited.modules))
    absent_rec = sorted(indexed - set(recoded.modules))
    print(f"\nB1 every indexed module present: deposited missing {absent_dep or 'none'}, "
          f"recoded missing {absent_rec or 'none'}")

    # --- B2: every module has a witness context, so its frame is derived not guessed ---
    witnesses = inv.witness_targets()
    unwitnessed = sorted(indexed - set(witnesses))
    print(f"B2 witness context for every module: {len(witnesses)}/{len(indexed)}, "
          f"missing {unwitnessed or 'none'}")

    # --- B3: every module has a valid realisation ---
    failing, gains = {}, []
    for module_id in sorted(indexed):
        problems = module_problems(deposited.modules[module_id], recoded.modules[module_id])
        if problems:
            failing[module_id] = problems
        gains.append(module_id)
    print(f"B3 valid realisation for every module: {len(gains) - len(failing)}/{len(gains)}")
    for module_id, problems in list(failing.items())[:5]:
        print(f"      {module_id}: {'; '.join(problems[:2])}")

    # --- B4: products preserve protein and interfaces, under both fusion sites ---
    compiled, product_problems = 0, {}
    for target in ARCHITECTURES:
        for fusion_site in FUSION_SITES:
            old = inv.compile_target(deposited, target, fusion_site=fusion_site)
            new = inv.compile_target(recoded, target, fusion_site=fusion_site)
            ids = [m for m, _v in new["modules"]]
            independent, join_problems = join_independently(recoded, ids)

            problems = list(join_problems)
            if independent != new["product"]:
                problems.append("independent join differs from compile_target")
            if len(new["product"]) != len(old["product"]):
                problems.append(f"product {len(new['product'])} vs {len(old['product'])} nt")
            frame = new["coding_interval"][0]
            old_coding = old["product"][frame:]
            new_coding = new["product"][frame:]
            old_p = str(Seq(old_coding[:len(old_coding) // 3 * 3]).translate())
            new_p = str(Seq(new_coding[:len(new_coding) // 3 * 3]).translate())
            if old_p != new_p:
                problems.append("assembled protein differs from the deposited build")
            if "*" in new_p:
                problems.append("stop codon in the assembled product")
            for name, site in BLACKLIST.items():
                for strand, pattern in (("forward", site), ("reverse", rc(site))):
                    if pattern in new["product"]:
                        problems.append(f"{name} {strand} site spans a junction")

            key = f"{target}/{fusion_site}"
            compiled += 1
            if problems:
                product_problems[key] = problems
            print(f"      {key:26s} {len(new['product']):>5} nt  {len(new_p):>4} aa  "
                  f"{'ok' if not problems else 'FAIL'}")
    print(f"B4 products valid: {compiled - len(product_problems)}/{compiled}")

    # --- B5: reload and recompile reproduce the same sequences ---
    scratch = Path(args.scratch) if args.scratch else Path(args.out).parent
    scratch.mkdir(parents=True, exist_ok=True)
    roundtrip = scratch / "_roundtrip.json"
    inv.save(recoded, roundtrip)
    reloaded = inv.load(roundtrip)
    same_version = reloaded.version == recoded.version
    same_modules = all(reloaded.modules[m].dna == recoded.modules[m].dna
                       for m in recoded.modules)
    recompiled = inv.compile_target(reloaded, ARCHITECTURES[0])
    original = inv.compile_target(recoded, ARCHITECTURES[0])
    same_product = recompiled["product"] == original["product"]
    roundtrip.unlink(missing_ok=True)
    print(f"B5 reload reproduces: version {same_version}, sequences {same_modules}, "
          f"recompiled product {same_product}")

    # --- B6: incompatible and missing records must fail, not compile ---
    rejected = {}
    broken = inv.derive(recoded, "broken", {})
    # Delete a module the target actually needs. Removing an arbitrary one proved nothing:
    # the first probe deleted a module this target never selects, so compiling succeeded and
    # the probe reported "not rejected" for a case it had not actually created.
    needed = [m for m, _v in inv.compile_target(recoded, ARCHITECTURES[0])["modules"]]
    victim = needed[0]
    del broken.modules[victim]
    try:
        inv.compile_target(broken, ARCHITECTURES[0])
        rejected["missing module"] = False
    except inv.IncompatibleInventory:
        rejected["missing module"] = True

    swapped = inv.derive(recoded, "swapped", {})
    target_ids = [m for m, _v in inv.compile_target(recoded, ARCHITECTURES[0])["modules"]]
    second = swapped.modules[target_ids[1]]
    swapped.modules[target_ids[1]] = inv.ModuleRecord(
        **{**second.__dict__, "dna": "TTTT" + second.dna[INTERFACE:],
           "five_interface": "TTTT"})
    try:
        inv.compile_target(swapped, ARCHITECTURES[0])
        rejected["incompatible interface"] = False
    except inv.IncompatibleInventory:
        rejected["incompatible interface"] = True

    edited = json.loads(Path(args.recoded).read_text(encoding="utf-8"))
    # Guarantee a real change. "A" + dna[1:] is a no-op on a sequence already starting with
    # A, and the probe then reported "not rejected" for an edit it had not actually made.
    original_dna = edited["modules"][victim]["dna"]
    edited["modules"][victim]["dna"] = ("T" if original_dna[0] != "T" else "G")         + original_dna[1:]
    tampered = scratch / "_tampered.json"
    tampered.write_text(json.dumps(edited), encoding="utf-8")
    try:
        inv.load(tampered)
        rejected["edited saved file"] = False
    except inv.InventoryInvalid:
        rejected["edited saved file"] = True
    tampered.unlink(missing_ok=True)

    print("B6 corruption probes (each must be rejected):")
    for name, caught in rejected.items():
        print(f"      {name:26s} {'rejected' if caught else 'NOT REJECTED'}")

    # --- adaptation gains, reported apart from the capability gate ---
    report_path = ROOT / "work" / "phaseb" / "recoding_report.json"
    gains_summary = {}
    if report_path.is_file():
        report = json.loads(report_path.read_text(encoding="utf-8"))
        scored = [o for o in report["objectives"].values() if "cai_after" in o]
        if scored:
            gains_summary = {
                "modules_scored": len(scored),
                "improved": sum(1 for o in scored if o["cai_after"] > o["cai_before"]),
                "mean_before": round(sum(o["cai_before"] for o in scored) / len(scored), 4),
                "mean_after": round(sum(o["cai_after"] for o in scored) / len(scored), 4),
            }
            print(f"\nadaptation gains (reported, not part of the gate): "
                  f"{gains_summary['improved']}/{gains_summary['modules_scored']} improved, "
                  f"mean CAI {gains_summary['mean_before']} -> {gains_summary['mean_after']}")

    passed = (not absent_dep and not absent_rec and not unwitnessed and not failing
              and not product_problems and same_version and same_modules and same_product
              and all(rejected.values()))
    print(f"\nGATE B: {'PASS' if passed else 'FAIL'}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "indexed_modules": len(indexed),
        "deposited_version": deposited.version, "recoded_version": recoded.version,
        "absent_from_deposited": absent_dep, "absent_from_recoded": absent_rec,
        "modules_without_witness": unwitnessed,
        "module_failures": failing,
        "products_checked": compiled, "product_failures": product_problems,
        "reload": {"version": same_version, "sequences": same_modules,
                   "recompiled_product": same_product},
        "corruption_probes": rejected,
        "adaptation_gains_reported_separately": gains_summary,
        "independence": ("constraints, frames and the product join recomputed here; "
                         "imports no clippr.qc and no clippr.codons"),
        "gate": "PASS" if passed else "FAIL",
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
