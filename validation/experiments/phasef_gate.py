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
import hashlib
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
import hashlib, json, sys
from pathlib import Path

area = Path(sys.argv[1])
out = {"stages": []}

def record(stage, ok, detail=None):
    out["stages"].append({"stage": stage, "ok": bool(ok), "detail": detail})
    return ok

def finish(code=0):
    print("JSON_BEGIN" + json.dumps(out) + "JSON_END")
    sys.exit(code)

import clippr
from clippr import inventories as inv, workflow as w
from clippr.ordering import HISTORICAL_OPOOL

out["import_origin"] = clippr.__file__
try:
    from clippr.manifest import source_fingerprint
    out["installed_source_fingerprint"] = source_fingerprint()["sha256"]
except Exception as exc:
    out["installed_source_fingerprint"] = f"unavailable: {exc}"
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
record("inventory", not out["validate"], out["validate"] or None)

compiled = inv.compile_target(library, "AAAAUGUGG")
out["compiled"] = {"available": compiled["available"],
                   "product_nt": compiled["product_nt"],
                   "modules": len(compiled["modules"]),
                   "reactions": compiled["reactions"]}
record("compile_target", compiled["available"])

targets = ["AAAAUGUGG", "UUACACGUGCGUAC"]
result = w.load_and_compile(area / "inventory.json", targets, area / "task2")
out["task2"] = {"ok": result.ok, "summary": result.summary}
if not record("load_and_compile", result.ok, result.failures[:3]):
    finish(1)

table = json.loads((area / "codon_table.json").read_text(encoding="utf-8"))
front = w.explore_interfaces(area / "inventory.json", targets, table, area / "front",
                             max_evaluations=12, wall_seconds=600)
out["front"] = {"ok": front.ok, "summary": front.summary}
if not record("explore_interfaces", front.ok, front.failures[:3]):
    finish(1)
front_path = front.artefacts["front"]

# --- R1 controls, in the installed environment ---------------------------------------
# Control 2: an explicitly supplied table.
selected = w.select_interface(area / "inventory.json", front_path, area / "selected",
                              codon_table=table)
out["selected"] = {"ok": selected.ok, "summary": selected.summary,
                   "recomputed": selected.data.get("objectives_recomputed"),
                   "mismatches": selected.data.get("binding_mismatches")}
if not record("select_interface_explicit", selected.ok, selected.failures[:3]):
    finish(1)
explicit_version = json.loads(
    Path(selected.artefacts["inventory"]).read_text(encoding="utf-8"))["version"]

# Control 1: the same selection with no table argument, resolving through the cache. This
# is the path that failed from an installed wheel; an explicit-table gate cannot see it.
default_sel = w.select_interface(area / "inventory.json", front_path,
                                 area / "selected_default")
default_version = None
if default_sel.ok:
    default_version = json.loads(
        Path(default_sel.artefacts["inventory"]).read_text(encoding="utf-8"))["version"]
out["select_default_table"] = {
    "ok": default_sel.ok, "summary": default_sel.summary,
    "version": default_version, "explicit_version": explicit_version,
    "same_as_explicit": default_version == explicit_version}
record("select_interface_cached_default",
       default_sel.ok and default_version == explicit_version,
       default_sel.failures[:3] or None)

# Control 3: a different effective table must be refused, publishing nothing.
altered = json.loads(json.dumps(table))
for _aa, _codons in altered.items():
    if isinstance(_codons, dict) and len(_codons) > 1:
        _c = sorted(_codons)[0]
        _codons[_c] = round(1.0 - float(_codons[_c]), 6)
        break
altered_sel = w.select_interface(area / "inventory.json", front_path, area / "altered",
                                 codon_table=altered)
out["altered_table_refused"] = {
    "refused": not altered_sel.ok,
    "published_nothing": not (area / "altered" / "inventory.json").exists(),
    "reasons": [f.get("reason", "")[:120] for f in (altered_sel.failures or [])][:2]}
record("altered_table_refused",
       not altered_sel.ok and not (area / "altered" / "inventory.json").exists())

recompiled = w.load_and_compile(selected.artefacts["inventory"], targets,
                                area / "recompiled")
out["recompiled"] = {"ok": recompiled.ok,
                     "version": recompiled.data.get("inventory_version")}
if not record("recompile_from_selected", recompiled.ok, recompiled.failures[:3]):
    finish(1)

# Collection optimisation, executed rather than merely present.
optimised = w.optimise_collection(area / "inventory.json", table, area / "collection",
                                  mode="greedy", max_proposals=42, wall_seconds=600)
if not record("optimise_collection", optimised.ok, optimised.failures[:3]):
    out["collection"] = {"ok": False, "summary": optimised.summary}
    finish(1)
optimised_items = w.order_items_for(optimised.artefacts["inventory"],
                                    area / "collection_items")
out["collection"] = {"ok": optimised.ok, "summary": optimised.summary,
                     "order_ok": optimised_items.ok,
                     "order_failures": optimised_items.failures[:3]}
record("collection_orders_cleanly", optimised_items.ok)

items_result = w.order_items_for(selected.artefacts["inventory"], area / "items")
if not record("order_items_for", items_result.ok, items_result.failures[:3]):
    finish(1)
items = json.loads(Path(items_result.artefacts["order_items"]).read_text(
    encoding="utf-8"))
out["order_items"] = {"ok": items_result.ok, "count": len(items),
                      "form": items_result.data["form"],
                      "ends": [list(e) for e in items_result.data["exposed_ends"]]}

# A mismatched front must be refused. The control is the optimised inventory produced
# above: a known different, valid inventory whose version is asserted to differ, rather
# than a copy of the same file when a historical Phase C run happens to be absent.
other_path = Path(optimised.artefacts["inventory"])
other_version = json.loads(other_path.read_text(encoding="utf-8"))["version"]
out["mismatch_control"] = {"other_version": other_version,
                           "source_version": library.version,
                           "genuinely_different": other_version != library.version}
mismatch = w.select_interface(other_path, front_path, area / "mismatch",
                              codon_table=table)
out["mismatch_refused"] = (not mismatch.ok) and other_version != library.version
record("mismatched_front_refused", out["mismatch_refused"])

order = w.plan_order(items, HISTORICAL_OPOOL, area / "task6")
out["task6"] = {"ok": order.ok, "summary": order.summary,
                "eligible": order.data["eligibility"]["eligible"],
                "per_pool_eligible": [bool(r.get("eligible"))
                                      for r in order.data["per_pool_eligibility"]],
                "pool_count": order.data["pool_count"],
                "optimality": order.data["optimality"],
                "vendor_approval": order.data["vendor_approval"]}
if not record("plan_order", order.ok, order.failures[:3]):
    finish(1)

package = w.write_package(
    [result, front, selected, recompiled, optimised, optimised_items, items_result,
     order], area / "package",
    inputs={"inventory": "inventory.json"})
out["package"] = str(package)
manifest = json.loads(Path(package).read_text(encoding="utf-8"))
out["package_ok"] = manifest["ok"]
# Hash every artefact the manifest names, so the driver can re-verify them after the
# package has been copied out of the temporary area.
digests = {}
for task in manifest.get("tasks", []):
    for name, rel in (task.get("artefacts") or {}).items():
        f = Path(package).parent / rel
        if f.is_file():
            digests[rel] = hashlib.sha256(f.read_bytes()).hexdigest()
out["package_artefact_sha256"] = digests
record("write_package", manifest["ok"])

finish(0)
'''


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", default=str(ROOT / "work" / "phaseb"
                                               / "inventory_recoded.json"))
    ap.add_argument("--out", default=str(ROOT / "work" / "phasef"),
                    help="directory (not a file) for gate_f.json and the kept evidence")
    ap.add_argument("--keep-environment", action="store_true",
                    help=("build the wheel and validation environment under --out and keep "
                          "them, so the installed notebook checks can use the same wheel"))
    args = ap.parse_args()

    source_inventory = Path(args.inventory)
    if not source_inventory.is_file():
        print(f"no saved inventory at {source_inventory}; run the Phase B recoding first")
        return 2

    # Absolute: the clean run executes with cwd=area, so a relative script path would
    # resolve against area and double the prefix.
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    # Created where it will live. Moving an already-created virtual environment breaks its
    # recorded interpreter paths, so --keep-environment must not relocate one afterwards.
    if args.keep_environment:
        area = out / "environment"
        if area.exists():
            # Never delete a previous run's retained environment: it is the evidence a
            # notebook result was bound to. Refuse and let the caller choose a new --out.
            print(f"{area} already exists. It holds an earlier run's wheel and "
                  f"environment, which later results may cite, so it is not removed. "
                  f"Point --out at a new directory.")
            return 2
        area.mkdir(parents=True)
        cleanup = False
    else:
        area = Path(tempfile.mkdtemp(prefix="clippr-phasef-"))
        cleanup = True

    record: dict = {"kept_environment": args.keep_environment}
    try:
        print("F1  building a wheel and installing it into a throwaway environment")
        build = subprocess.run([str(PYTHON), "-m", "build", "--wheel", "--outdir",
                                str(area / "dist")], cwd=ROOT, capture_output=True,
                               text=True)
        (out / "wheel_build.log").write_text((build.stdout or "") + (build.stderr or ""),
                                             encoding="utf-8")
        wheels = list((area / "dist").glob("*.whl"))
        if build.returncode or not wheels:
            print(f"    wheel build failed: {build.stderr.strip()[-200:]}")
            _write_record(out, {**record, "gate": "FAIL",
                                "stage": "wheel_build",
                                "detail": build.stderr.strip()[-2000:]})
            return 1
        wheel = wheels[0]
        print(f"    built {wheel.name}")

        venv = area / "env"
        subprocess.run([str(PYTHON), "-m", "venv", str(venv)], check=True,
                       capture_output=True)
        python = venv / "Scripts" / "python.exe"
        if not python.is_file():
            python = venv / "bin" / "python"
        install = subprocess.run([str(python), "-m", "pip", "install", "--quiet",
                                  str(wheel)], capture_output=True, text=True)
        if install.returncode:
            print(f"    install failed: {install.stderr.strip()[-300:]}")
            _write_record(out, {**record, "gate": "FAIL", "stage": "install",
                                "detail": install.stderr.strip()[-2000:]})
            return 1
        print("    installed")

        freeze = subprocess.run([str(python), "-m", "pip", "freeze"],
                                capture_output=True, text=True)
        record["environment_dependencies"] = sorted(
            l for l in (freeze.stdout or "").splitlines() if l.strip())

        # A saved inventory carries its own sequences. The codon table is handed across and
        # also placed in the cache, so the cached-default control resolves offline.
        shutil.copy(source_inventory, area / "inventory.json")
        table_src = ROOT / "data" / "codon_tables" / "kazusa_3055.json"
        shutil.copy(table_src, area / "codon_table.json")
        cache_tables = area / "cache" / "codon_tables"
        cache_tables.mkdir(parents=True, exist_ok=True)
        shutil.copy(table_src, cache_tables / "kazusa_3055.json")
        script = area / "run_clean.py"
        script.write_text(CLEAN_SCRIPT, encoding="utf-8")

        record["inputs"] = {"inventory": {"path": str(source_inventory),
                                          "sha256": _sha256(source_inventory)},
                            "codon_table": {"path": str(table_src),
                                            "sha256": _sha256(table_src)}}
        record["wheel"] = {"name": wheel.name, "sha256": _sha256(wheel),
                           "path": str(wheel.resolve()) if args.keep_environment else None}

        env = {k: v for k, v in os.environ.items()
               if k not in ("PYTHONPATH", "PYTHONHOME", "PYTHONSTARTUP")}
        env.update({"CLIPPR_CACHE_DIR": str(area / "cache"), "PYTHONIOENCODING": "utf-8",
                    "PYTHONNOUSERSITE": "1"})

        print("\nF2  running the documented paths in the clean environment")
        got = subprocess.run([str(python), str(script), str(area)], cwd=area,
                             capture_output=True, text=True, env=env, timeout=1800)
        # Written before anything can go wrong, and before any cleanup.
        (out / "clean_run.stdout.log").write_text(got.stdout or "", encoding="utf-8")
        (out / "clean_run.stderr.log").write_text(got.stderr or "", encoding="utf-8")

        if "JSON_BEGIN" not in got.stdout:
            print(f"    clean run produced no result (exit {got.returncode})")
            print(f"    stderr: {(got.stderr or '').strip()[-400:]}")
            print(f"    logs kept in {out}")
            _write_record(out, {**record, "gate": "FAIL", "stage": "clean_run",
                                "exit_code": got.returncode,
                                "detail": (got.stderr or "").strip()[-2000:]})
            return 1
        payload = json.loads(got.stdout.split("JSON_BEGIN")[1].split("JSON_END")[0])

        inside = str(venv.resolve()) in str(Path(payload["import_origin"]).resolve())
        print(f"    imported from {payload['import_origin']}")
        print(f"    resolves inside the clean environment: {inside}")
        print(f"    workflow API present: {len(payload['api_present'])}/9")
        print(f"    inventory loaded without Table S1: {payload['inventory']['modules']} "
              f"modules, version {payload['inventory']['version']}")
        print(f"    validate: {payload['validate'] or 'clean'}")
        print(f"    task 2: {payload['task2']['summary']}")
        failed_stages = [s["stage"] for s in payload.get("stages", []) if not s["ok"]]
        print(f"    stages failed: {failed_stages or 'none'}")

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
        print(f"    selected (explicit table): {payload['selected']['summary']}")
        sd = payload.get("select_default_table", {})
        print(f"    selected (cached default): ok={sd.get('ok')} "
              f"same version as explicit={sd.get('same_as_explicit')}")
        print(f"    altered table refused: {payload['altered_table_refused']['refused']}, "
              f"published nothing: "
              f"{payload['altered_table_refused']['published_nothing']}")
        print(f"    recompiled from {payload['recompiled']['version']}")
        print(f"    order items: {payload['order_items']['count']} "
              f"{payload['order_items']['form']}")
        mc = payload["mismatch_control"]
        print(f"    mismatch control is a different inventory: "
              f"{mc['genuinely_different']} ({mc['other_version']} vs "
              f"{mc['source_version']})")
        print(f"    mismatched front refused: {payload['mismatch_refused']}")
        print(f"    order plan: {payload['task6']['summary']}; "
              f"optimality={payload['task6']['optimality']}")
        print(f"    vendor approval: {payload['task6']['vendor_approval']}")

        # --- F4: the package must survive leaving the temporary area ---
        print("\nF4  relocating the result package and re-verifying its artefacts")
        kept_package = out / "package"
        if kept_package.exists():
            print(f"{kept_package} already exists; not overwriting earlier evidence")
            return 2
        shutil.copytree(area / "package", kept_package)
        expected = payload.get("package_artefact_sha256", {})
        relocated_problems = []
        for rel, digest in sorted(expected.items()):
            f = kept_package / rel
            if not f.is_file():
                relocated_problems.append(f"{rel}: missing after relocation")
            elif _sha256(f) != digest:
                relocated_problems.append(f"{rel}: content changed after relocation")
        print(f"    {len(expected)} artefact(s) checked after relocation; "
              f"problems: {relocated_problems or 'none'}")

        passed = (inside and len(payload["api_present"]) == 9
                  and not payload["validate"]
                  and not failed_stages
                  and payload["compiled"]["available"]
                  and payload["task2"]["ok"] and payload["task6"]["ok"]
                  and payload["task6"]["eligible"]
                  and all(payload["task6"]["per_pool_eligible"])
                  and payload["front"]["ok"] and payload["selected"]["ok"]
                  and payload["selected"]["recomputed"]
                  and sd.get("ok") and sd.get("same_as_explicit")
                  and payload["altered_table_refused"]["refused"]
                  and payload["altered_table_refused"]["published_nothing"]
                  and payload["recompiled"]["ok"]
                  and payload["order_items"]["form"] == "assembly_ready"
                  and payload["mismatch_control"]["genuinely_different"]
                  and payload["mismatch_refused"]
                  and payload["collection"]["ok"]
                  and payload["collection"]["order_ok"]
                  and payload["package_ok"] and not leaked
                  and expected and not relocated_problems)
        print(f"\nGATE F: {'PASS' if passed else 'FAIL'}")

        record.update({
            "clean_environment": payload,
            "import_resolves_inside_clean_env": inside,
            "installed_import_origin": payload["import_origin"],
            "installed_source_fingerprint": payload.get("installed_source_fingerprint"),
            "developer_paths_in_artefacts": leaked,
            "failed_stages": failed_stages,
            "relocated_package": str(kept_package),
            "relocated_artefact_problems": relocated_problems,
            "supplied_input_requirements": {
                "saved inventory": "none -- it carries its own sequences",
                "deposited inventory": ("Supplementary Table S1 must be supplied; it is "
                                        "deliberately not shipped, see NOTICE.md"),
                "codon table": ("supplied by the caller, or fetched once from Kazusa and "
                                "cached"),
                "host screening": "genome fetched from NCBI on first use, or disabled",
            },
            "gate": "PASS" if passed else "FAIL",
        })
        if args.keep_environment:
            # Added only after the gate's own clean run, so the clean-install claim still
            # covers the wheel and its declared dependencies alone. The notebooks call
            # display(), which Colab provides and a bare wheel environment does not.
            extra = subprocess.run([str(python), "-m", "pip", "install", "--quiet",
                                    "ipython"], capture_output=True, text=True)
            record["notebook_runtime_dependency"] = {
                "package": "ipython",
                "installed_after_the_gate_run": True,
                "ok": extra.returncode == 0,
                "why": ("notebook cells call display(); it is not a CLIPPR dependency and "
                        "is not part of the clean-install claim"),
            }
            if extra.returncode:
                print(f"    warning: could not add ipython: "
                      f"{extra.stderr.strip()[-200:]}")
            record["validation_python"] = str(python.resolve())
            record["wheel_path"] = str(wheel.resolve())
            record["wheel_sha256"] = record["wheel"]["sha256"]
        _write_record(out, record)
        return 0 if passed else 1
    finally:
        if cleanup:
            shutil.rmtree(area, ignore_errors=True)


def _write_record(out: Path, record: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    path = out / "gate_f.json"
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {path}")


if __name__ == "__main__":
    raise SystemExit(main())
