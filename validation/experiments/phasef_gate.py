"""Phase F gate — the documented workflow, executed from a clean installed package.

Gate F, from the full-workflow plan: run the documented API paths from a clean install;
verify a saved inventory reopens and compiles **without the reference repository or hidden
author paths**; use distributable synthetic fixtures where primary sequences cannot be
bundled; and state the supplied-input requirements plainly.

The saved inventory is the interesting case. It carries its own sequences, so a clean
environment can load and compile it with no access to Supplementary Table S1 — which is the
file this package deliberately does not ship. That is what makes a saved inventory a
distributable artefact rather than a pointer to one.

The clean environment is built with `PYTHONPATH`, `PYTHONHOME` and `PYTHONSTARTUP` dropped, so
an inherited path cannot let it import the source checkout and pass on the strength of code
that was never installed.

    python validation/experiments/phasef_gate.py
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PYTHON = ROOT / ".venv" / "Scripts" / "python.exe"
if not PYTHON.is_file():
    PYTHON = Path(sys.executable)

#: Runs inside the clean environment. Everything it needs is either installed with the
#: package or handed to it as a file; it reads nothing from the developer's checkout.
CLEAN_SCRIPT = '''
import json, sys
from pathlib import Path

area = Path(sys.argv[1])
out = {}

import clippr
from clippr import inventories as inv, workflow as w
from clippr.ordering import HISTORICAL_OPOOL

out["import_origin"] = clippr.__file__
out["api_present"] = sorted(
    n for n in ("design_for_synthesis", "load_and_compile", "recode_for_host",
                "optimise_collection", "explore_interfaces", "select_interface",
                "order_items_for", "plan_order", "write_package")
    if hasattr(w, n))

# A saved inventory carries its own sequences: no Table S1 required.
library = inv.load(area / "inventory.json")
out["inventory"] = {"modules": len(library), "version": library.version,
                    "label": library.label}
out["validate"] = inv.validate(library)

compiled = inv.compile_target(library, "AAAAUGUGG")
out["compiled"] = {"available": compiled["available"],
                   "product_nt": compiled["product_nt"],
                   "modules": len(compiled["modules"]),
                   "reactions": compiled["reactions"]}

targets = ["AAAAUGUGG", "UUACACGUGCGUAC"]
result = w.load_and_compile(area / "inventory.json", targets, area / "task2")
out["task2"] = {"ok": result.ok, "summary": result.summary}

# The whole connected path, in the clean environment. Loading and compiling alone cannot
# establish selection or assembly-ready ordering, which is what this gate now claims.
table = json.loads((area / "codon_table.json").read_text(encoding="utf-8"))
front = w.explore_interfaces(area / "inventory.json", targets, table, area / "front",
                             max_evaluations=12, wall_seconds=600)
out["front"] = {"ok": front.ok, "summary": front.summary}

selected = w.select_interface(area / "inventory.json", front.artefacts["front"],
                              area / "selected", codon_table=table)
out["selected"] = {"ok": selected.ok, "summary": selected.summary,
                   "recomputed": selected.data.get("objectives_recomputed"),
                   "mismatches": selected.data.get("binding_mismatches")}

recompiled = w.load_and_compile(selected.artefacts["inventory"], targets,
                                area / "recompiled")
out["recompiled"] = {"ok": recompiled.ok,
                     "version": recompiled.data["inventory_version"]}

# Collection optimisation, executed rather than merely present. Checking that the API
# exists could not catch an optimiser whose output the next step refuses.
optimised = w.optimise_collection(area / "inventory.json", table, area / "collection",
                                  mode="greedy", max_proposals=42, wall_seconds=600)
optimised_items = w.order_items_for(optimised.artefacts["inventory"],
                                    area / "collection_items")
out["collection"] = {"ok": optimised.ok, "summary": optimised.summary,
                     "order_ok": optimised_items.ok,
                     "order_failures": optimised_items.failures[:3]}

items_result = w.order_items_for(selected.artefacts["inventory"], area / "items")
items = json.loads(Path(items_result.artefacts["order_items"]).read_text(
    encoding="utf-8"))
out["order_items"] = {"ok": items_result.ok, "count": len(items),
                      "form": items_result.data["form"],
                      "ends": [list(e) for e in items_result.data["exposed_ends"]]}

# A mismatched front must be refused, in the clean environment too.
mismatch = w.select_interface(area / "other_inventory.json", front.artefacts["front"],
                              area / "mismatch", codon_table=table)
out["mismatch_refused"] = not mismatch.ok

# Task 6 needs no supplied input at all.
order = w.plan_order(items, HISTORICAL_OPOOL, area / "task6")
out["task6"] = {"ok": order.ok, "summary": order.summary,
                "cost_available": order.data["cost"]["available"]}

package = w.write_package(
    [result, front, selected, recompiled, optimised, optimised_items, items_result,
     order], area / "package",
    inputs={"inventory": "inventory.json"})
out["package"] = str(package)
out["package_ok"] = json.loads(Path(package).read_text(encoding="utf-8"))["ok"]

print("JSON_BEGIN" + json.dumps(out) + "JSON_END")
'''


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", default=str(ROOT / "work" / "phaseb"
                                               / "inventory_recoded.json"))
    ap.add_argument("--out", default=str(ROOT / "work" / "phasef"))
    args = ap.parse_args()

    source_inventory = Path(args.inventory)
    if not source_inventory.is_file():
        print(f"no saved inventory at {source_inventory}; run the Phase B recoding first")
        return 2

    area = Path(tempfile.mkdtemp(prefix="clippr-phasef-"))
    try:
        print("F1  building a wheel and installing it into a throwaway environment")
        build = subprocess.run([str(PYTHON), "-m", "build", "--wheel", "--outdir",
                                str(area / "dist")], cwd=ROOT, capture_output=True,
                               text=True)
        wheels = list((area / "dist").glob("*.whl"))
        if build.returncode or not wheels:
            print(f"    wheel build failed: {build.stderr.strip()[-200:]}")
            return 1
        print(f"    built {wheels[0].name}")

        venv = area / "env"
        subprocess.run([str(PYTHON), "-m", "venv", str(venv)], check=True,
                       capture_output=True)
        python = venv / "Scripts" / "python.exe"
        if not python.is_file():
            python = venv / "bin" / "python"
        install = subprocess.run([str(python), "-m", "pip", "install", "--quiet",
                                  str(wheels[0])], capture_output=True, text=True)
        if install.returncode:
            print(f"    install failed: {install.stderr.strip()[-300:]}")
            return 1
        print("    installed")

        # A saved inventory carries its own sequences. A second, different inventory and a
        # codon table are handed across so the clean run can exercise selection, ordering and
        # the mismatched-front refusal.
        shutil.copy(source_inventory, area / "inventory.json")
        shutil.copy(ROOT / "data" / "codon_tables" / "kazusa_3055.json",
                    area / "codon_table.json")
        other = ROOT / "work" / "phasec" / "inventory_greedy.json"
        shutil.copy(other if other.is_file() else source_inventory,
                    area / "other_inventory.json")
        script = area / "run_clean.py"
        script.write_text(CLEAN_SCRIPT, encoding="utf-8")

        env = {k: v for k, v in os.environ.items()
               if k not in ("PYTHONPATH", "PYTHONHOME", "PYTHONSTARTUP")}
        env.update({"CLIPPR_CACHE_DIR": str(area / "cache"), "PYTHONIOENCODING": "utf-8",
                    "PYTHONNOUSERSITE": "1"})

        print("\nF2  running the documented paths in the clean environment")
        got = subprocess.run([str(python), str(script), str(area)], cwd=area,
                             capture_output=True, text=True, env=env, timeout=1800)
        if "JSON_BEGIN" not in got.stdout:
            print(f"    clean run produced no result (exit {got.returncode})")
            print(f"    stderr: {(got.stderr or '').strip()[-400:]}")
            return 1
        payload = json.loads(got.stdout.split("JSON_BEGIN")[1].split("JSON_END")[0])

        inside = str(venv.resolve()) in str(Path(payload["import_origin"]).resolve())
        print(f"    imported from {payload['import_origin']}")
        print(f"    resolves inside the clean environment: {inside}")
        print(f"    workflow API present: {len(payload['api_present'])}/9")
        print(f"    inventory loaded without Table S1: {payload['inventory']['modules']} "
              f"modules, version {payload['inventory']['version']}")
        print(f"    validate: {payload['validate'] or 'clean'}")
        print(f"    compiled: {payload['compiled']}")
        print(f"    task 2: {payload['task2']['summary']}")
        print(f"    task 6: {payload['task6']['summary']}")
        print(f"    package written and ok: {payload['package_ok']}")

        # --- F3: no author path leaked into the artefacts ---
        print("\nF3  checking artefacts for developer-specific paths")
        leaked = []
        for path in (area / "package").rglob("*.json"):
            text = path.read_text(encoding="utf-8", errors="replace")
            for marker in (str(ROOT), "syedz", "CodePractice", "grasp-library-designer"):
                if marker in text:
                    leaked.append(f"{path.name}: {marker}")
        print(f"    developer paths found in the package: {leaked or 'none'}")

        print(f"    front: {payload['front']['summary']}")
        print(f"    selected: {payload['selected']['summary']}")
        print(f"    recompiled from {payload['recompiled']['version']}")
        print(f"    order items: {payload['order_items']['count']} "
              f"{payload['order_items']['form']}, ends "
              f"{payload['order_items']['ends']}")
        print(f"    mismatched front refused: {payload['mismatch_refused']}")
        print(f"    collection optimisation: {payload['collection']['summary']}")
        print(f"    its output orders cleanly: {payload['collection']['order_ok']}")

        passed = (inside and len(payload["api_present"]) == 9
                  and not payload["validate"]
                  and payload["compiled"]["available"]
                  and payload["task2"]["ok"] and payload["task6"]["ok"]
                  and payload["front"]["ok"] and payload["selected"]["ok"]
                  and payload["selected"]["recomputed"]
                  and payload["recompiled"]["ok"]
                  and payload["order_items"]["form"] == "assembly_ready"
                  and payload["mismatch_refused"]
                  and payload["collection"]["ok"]
                  and payload["collection"]["order_ok"]
                  and payload["package_ok"] and not leaked)
        print(f"\nGATE F: {'PASS' if passed else 'FAIL'}")

        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        (out / "gate_f.json").write_text(json.dumps({
            "wheel": wheels[0].name,
            "clean_environment": payload,
            "import_resolves_inside_clean_env": inside,
            "developer_paths_in_artefacts": leaked,
            "supplied_input_requirements": {
                "saved inventory": "none -- it carries its own sequences",
                "deposited inventory": ("Supplementary Table S1 must be supplied; it is "
                                        "deliberately not shipped, see NOTICE.md"),
                "codon table": ("supplied by the caller, or fetched once from Kazusa and "
                                "cached"),
                "host screening": "genome fetched from NCBI on first use, or disabled",
            },
            "gate": "PASS" if passed else "FAIL",
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote {out / 'gate_f.json'}")
        return 0 if passed else 1
    finally:
        shutil.rmtree(area, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
