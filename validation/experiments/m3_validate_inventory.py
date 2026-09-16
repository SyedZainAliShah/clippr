"""M3 / E4 — assemble from the redesigned inventory and check what per-module work cannot.

Recoding each module in isolation proves nothing about the assembly. Two modules can each
satisfy every constraint and still create a forbidden site *across the join between them*,
because neither optimiser ever saw the other's sequence. This joins the redesigned records,
verifies every overlap, translates the product and checks the constraints on the assembled
molecule -- which is the only place a junction-spanning defect can appear.

The join is written out here rather than taken from `products.reconstruct`, so a fault in the
production reconstructor cannot hide a fault in the inventory.

Also separated deliberately, per the orchestrator: **synthesis constraints on individual
modules** are not the same measurement as **repeats in the assembled product**. Reusing a
module necessarily repeats its DNA in every product that contains it, so a product-level
repeat count is not a defect in the module.

    python validation/experiments/m3_validate_inventory.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from Bio.Seq import Seq                                           # noqa: E402

from clippr.parts import select                                   # noqa: E402
from clippr.products import load_inserts                          # noqa: E402

SIX = ["AAAAUGUGG", "GCUAAAGAC", "UUACACGUG", "CGUACGUAC", "AUCGAUCGA", "GGCCAAUUG"]
BLACKLIST = {"BsaI": "GGTCTC", "BbsI": "GAAGAC", "SapI": "GCTCTTC"}
INTERFACE = 4
PRODUCT_FRAME = 1
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def join(module_ids: list[str], inventory: dict[str, str]) -> tuple[str, list[str]]:
    """Splice modules at their four-base overlaps, reporting every join that fails."""
    problems = []
    product = inventory[module_ids[0]]
    for left, right in zip(module_ids, module_ids[1:]):
        a, b = inventory[left], inventory[right]
        if a[-INTERFACE:] != b[:INTERFACE]:
            problems.append(f"{left} ends {a[-INTERFACE:]} but {right} starts "
                            f"{b[:INTERFACE]}")
        product += b[INTERFACE:]
    return product, problems


def sites_in(product: str) -> list[str]:
    found = []
    for name, site in BLACKLIST.items():
        for strand, pattern in (("forward", site), ("reverse", rc(site))):
            at = product.find(pattern)
            if at >= 0:
                found.append(f"{name} {strand} at {at}")
    return found


def interfaces_vs_originals(inventory: dict[str, str],
                            originals: dict[str, str]) -> list[str]:
    """Every interface that has drifted from the deposited original, named.

    Compatible overlaps are not the promise. The promise is that each module keeps the
    *original's* four bases, so a recoded module drops into an assembly built from unrecoded
    neighbours. Two modules could agree perfectly with each other and both have drifted from
    the kit; an overlap check alone would not notice.
    """
    drift = []
    for pid, seq in sorted(inventory.items()):
        original = originals[pid]
        if seq[:INTERFACE] != original[:INTERFACE]:
            drift.append(f"{pid} 5' {seq[:INTERFACE]} != {original[:INTERFACE]}")
        if seq[-INTERFACE:] != original[-INTERFACE:]:
            drift.append(f"{pid} 3' {seq[-INTERFACE:]} != {original[-INTERFACE:]}")
    return drift


def duplicated_20mers(seq: str) -> int:
    windows = [seq[i:i + 20] for i in range(len(seq) - 19)]
    return len(windows) - len(set(windows))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pilot", default=str(ROOT / "work" / "m3" / "inventory_pilot.json"))
    ap.add_argument("--table", default=str(ROOT / "data" / "grasp_supp" / "Table S1.xlsx"))
    args = ap.parse_args()

    pilot = json.loads(Path(args.pilot).read_text(encoding="utf-8"))
    redesigned = {pid: rec["sequence"] for pid, rec in pilot["modules"].items()}
    originals, source = load_inserts(args.table)

    print(f"redesigned inventory: {len(redesigned)} modules  |  context "
          f"{pilot['context']['host']}, code {pilot['context']['genetic_code']}")

    # One versioned sequence per identity: every product must use the same record.
    identities = {pid: hashlib.sha256(seq.encode()).hexdigest()[:16]
                  for pid, seq in redesigned.items()}
    if len(set(redesigned)) != len(redesigned):
        print("  FAIL: a module identity maps to more than one sequence")
        return 1
    same_length = [pid for pid, seq in redesigned.items()
                   if len(seq) != len(originals[pid])]
    print(f"  one sequence per identity: yes  |  length preserved: "
          f"{len(redesigned) - len(same_length)}/{len(redesigned)}")

    # Compatible overlaps are not the promise. The promise is that each module keeps the
    # *original's* four-base interfaces, so a recoded module drops into an assembly built
    # from unrecoded neighbours. A pair of modules could agree with each other perfectly and
    # both have drifted from the kit, and an overlap check alone would not notice.
    interface_drift = interfaces_vs_originals(redesigned, originals)
    print(f"  interfaces identical to the deposited originals: "
          f"{len(redesigned) - len({d.split()[0] for d in interface_drift})}"
          f"/{len(redesigned)}")
    for drift in interface_drift[:5]:
        print(f"      {drift}")

    failures, rows = [], []
    for target in SIX:
        plan = select(target)
        ids = [m.plasmid_id for m in plan.modules]
        missing = [i for i in ids if i not in redesigned]
        if missing:
            failures.append((target, f"inventory lacks {missing}"))
            continue

        new_product, join_problems = join(ids, redesigned)
        old_product, _ = join(ids, originals)
        coding_new = new_product[PRODUCT_FRAME:]
        coding_new = coding_new[:len(coding_new) // 3 * 3]
        coding_old = old_product[PRODUCT_FRAME:]
        coding_old = coding_old[:len(coding_old) // 3 * 3]
        protein_new = str(Seq(coding_new).translate())
        protein_old = str(Seq(coding_old).translate())

        problems = list(join_problems)
        if len(new_product) != len(old_product):
            problems.append(f"product length {len(new_product)} vs {len(old_product)}")
        if protein_new != protein_old:
            problems.append("assembled protein differs from the original inventory's")
        if "*" in protein_new:
            problems.append(f"stop codon in the assembled product at {protein_new.index('*')}")
        junction_sites = sites_in(new_product)
        if junction_sites:
            problems.append(f"forbidden site in the assembled product: {junction_sites}")

        changed = sum(1 for a, b in zip(new_product, old_product) if a != b)
        rows.append({
            "target": target, "modules": len(ids),
            "product_nt": len(new_product), "protein_aa": len(protein_new),
            "bases_changed": changed,
            "percent_changed": round(100 * changed / len(old_product), 1),
            "duplicated_20mers_new": duplicated_20mers(new_product),
            "duplicated_20mers_old": duplicated_20mers(old_product),
            "problems": problems,
        })
        if problems:
            failures.append((target, "; ".join(problems[:2])))
        status = "ok" if not problems else "FAIL"
        print(f"  {target:12s} {len(new_product):>5} nt  {len(protein_new):>4} aa  "
              f"{changed:>4} bases changed ({rows[-1]['percent_changed']:>4.1f}%)  "
              f"dup20mers {rows[-1]['duplicated_20mers_old']} -> "
              f"{rows[-1]['duplicated_20mers_new']}  {status}")

    print(f"\nassembled products valid: {len(rows) - len(failures)}/{len(SIX)}")
    for target, why in failures:
        print(f"   {target}: {why}")

    # Product-level sharing is a property of reuse, not of the recoding.
    if len(rows) == len(SIX):
        new_products = {t: join([m.plasmid_id for m in select(t).modules], redesigned)[0]
                        for t in SIX}
        old_products = {t: join([m.plasmid_id for m in select(t).modules], originals)[0]
                        for t in SIX}
        from clippr.homology import longest_shared
        worst_new = max(longest_shared(new_products[a], new_products[b])[0]
                        for a, b in combinations(SIX, 2))
        worst_old = max(longest_shared(old_products[a], old_products[b])[0]
                        for a, b in combinations(SIX, 2))
        print(f"\nlongest shared tract between two kit products: {worst_old} nt before, "
              f"{worst_new} nt after recoding")
        print("   recoding shared modules cannot reduce this: both products contain the same"
              " records, so they change identically")

    # Corruption probes. A validator that only ever sees correct input reports a number, not
    # a verdict: these establish that the checks above would actually object.
    print("\ncorruption probes (each must be rejected):")
    probes, caught = [], 0
    plan_ids = [m.plasmid_id for m in select(SIX[0]).modules]
    first = plan_ids[0]

    # Judged by `interfaces_vs_originals`, the same function the real inventory goes
    # through -- a probe scored by its own private test would prove nothing about the gate.
    altered_entry = dict(redesigned)
    head = redesigned[first]
    flipped = "T" if head[0] != "T" else "G"
    altered_entry[first] = flipped + head[1:]
    probes.append(("external entry base altered on the first module",
                   bool(interfaces_vs_originals(altered_entry, originals))))

    broken_join = dict(redesigned)
    broken_join[plan_ids[1]] = "TTTT" + redesigned[plan_ids[1]][INTERFACE:]
    _seq, join_problems = join(plan_ids, broken_join)
    probes.append(("junction overhang mutated on the second module", bool(join_problems)))

    planted = dict(redesigned)
    mid = len(redesigned[plan_ids[1]]) // 2
    planted[plan_ids[1]] = (redesigned[plan_ids[1]][:mid] + "GGTCTC"
                            + redesigned[plan_ids[1]][mid + 6:])
    product, _ = join(plan_ids, planted)
    probes.append(("BsaI site planted inside a module", bool(sites_in(product))))

    for name, detected in probes:
        caught += bool(detected)
        print(f"   {name:52s} {'rejected' if detected else 'NOT REJECTED'}")

    out = ROOT / "work" / "m3"
    (out / "inventory_validation.json").write_text(json.dumps(
        {"modules": len(redesigned), "identities": identities, "products": rows,
         "interfaces_match_originals": not interface_drift,
         "interface_drift": interface_drift,
         "corruption_probes": {name: bool(d) for name, d in probes},
         "failures": [{"target": t, "why": w} for t, w in failures],
         "shared_tract_nt": {"before": worst_old, "after": worst_new}
         if len(rows) == len(SIX) else None,
         "independence": "join written out here; production reconstructor not used",
         "separation": ("module synthesis constraints and product repeat counts are "
                        "reported apart; reuse necessarily repeats module DNA")},
        indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {out / 'inventory_validation.json'}")
    return 1 if (failures or interface_drift or caught != len(probes)) else 0


if __name__ == "__main__":
    raise SystemExit(main())
