"""Run a notebook's code cells against an installed build, from a scratch directory.

`validation/run_notebook.py` puts the working tree's `src/` first on `sys.path`, so it proves
the notebook agrees with the code being edited. It cannot catch the failure that matters to a
user: the notebook depending on something the **installed** package does not have.

This runs the notebook under the interpreter of an environment where CLIPPR is installed, and
refuses to proceed if the import does not resolve inside that environment. Earlier versions
deleted `ROOT/src` from `sys.path` and then `chdir`-ed into the repository, which left every
file in the checkout reachable -- a notebook could read an input nobody had declared and the
run would still look clean. The working directory is now a scratch directory and inputs are
staged into it explicitly, with their hashes recorded.

    <installed-env>/python validation/notebook_on_published_build.py \\
        --notebook notebooks/CLIPPR_designer.ipynb \\
        --workdir work/run/notebook_designer \\
        --out work/run/notebook_designer.json \\
        --input data/codon_tables/kazusa_3055.json=data/codon_tables/kazusa_3055.json

Cells whose first line is a `%` or `!` magic are the Colab install and upload cells; they are
skipped and named. Any other cell that cannot run is an error, not a skip.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

#: Written into a working directory this runner created, so a later run can tell its
#: own scratch space from a directory someone else owns.
OWNERSHIP_MARKER = ".clippr-run"

#: Lines a run may drop before executing a cell: dependency installation and extension
#: loading, which mean nothing outside IPython. Dropping the LINE is all this permits.
SKIPPABLE_LINES = ("%pip", "!pip", "!git", "%load_ext")

#: Cell magics this runner knows how to honour by executing the body underneath them.
#: `%%capture` only suppresses display, so running the body is the faithful reading.
SUPPORTED_CELL_MAGICS = ("%%capture",)

#: Kept for callers that import it. `%%capture` and `%load_ext` were once here and let a
#: whole cell vanish: a required cell beginning `%%capture` and then raising reported
#: "1 skipped" and PASS.
INSTALL_MAGICS = SKIPPABLE_LINES


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def classify_cell(src: str) -> tuple[list[str], list[str], list[str]]:
    """Split a cell into code to run, magic lines to drop, and magic we will not excuse.

    Line by line, deliberately. Skipping the whole cell whenever it held a permitted
    magic is what let required work disappear: naming `%time` unsupported left the
    mechanism in place, so `%%capture` followed by a raise still reported PASS.
    """
    keep, dropped, unsupported = [], [], []
    for ln in src.splitlines():
        bare = ln.strip()
        if not ln.lstrip().startswith(("%", "!")):
            keep.append(ln)
        elif bare.startswith(SUPPORTED_CELL_MAGICS):
            dropped.append(bare)            # honoured by running the body underneath
        elif bare.startswith(SKIPPABLE_LINES):
            dropped.append(bare)
        else:
            unsupported.append(bare)
    return keep, dropped, unsupported


def _has_statement(source: str) -> bool:
    """Whether source holds work, as opposed to blank lines and comments."""
    import ast

    try:
        return bool(ast.parse(source).body)
    except SyntaxError:
        return True         # not ours to excuse; let the executor report it


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--notebook", required=True)
    ap.add_argument("--workdir", required=True,
                    help="scratch directory to run in; the repository is not used as cwd")
    ap.add_argument("--out", required=True, help="JSON report path")
    ap.add_argument("--input", action="append", default=[], metavar="DEST=SRC",
                    help="stage SRC into the working directory at DEST; repeatable")
    ap.add_argument("--env", action="append", default=[], metavar="NAME=VALUE",
                    help=("set an environment variable before the build is imported. A "
                          "relative VALUE resolves inside the working directory, which is "
                          "how CLIPPR_CACHE_DIR points at a staged cache so the bytes the "
                          "run consumes are declared and hashed. Repeatable."))
    ap.add_argument("--require-input", action="append", default=[], metavar="DEST",
                    help=("a declared input this notebook needs; the run refuses to "
                          "start unless it was staged. Repeatable."))
    ap.add_argument("--expect", action="append", default=[], metavar="EXPR",
                    help=("a postcondition evaluated in the notebook namespace after the "
                          "run; all must hold. Without one, a notebook that guarded "
                          "every step on a missing input reported ten executed cells, "
                          "no workflow and PASS. Repeatable."))
    args = ap.parse_args()

    notebook = Path(args.notebook).resolve()
    workdir = Path(args.workdir).resolve()
    out = Path(args.out).resolve()

    # A run owns its working directory or refuses to use it. The previous version called
    # shutil.rmtree on any existing --workdir before validating anything, so one mistyped
    # path destroyed whatever was there.
    marker = workdir / OWNERSHIP_MARKER
    if workdir.exists():
        if not workdir.is_dir():
            print(f"--workdir exists and is not a directory: {workdir}")
            return 2
        contents = [q for q in workdir.iterdir() if q.name != OWNERSHIP_MARKER]
        if contents and not marker.is_file():
            print(f"refusing to use {workdir}: it is not empty and carries no "
                  f"{OWNERSHIP_MARKER} marker, so this run does not own it. Point "
                  f"--workdir at a new path.")
            return 2
        for q in contents:                      # ours, by the marker; safe to clear
            shutil.rmtree(q) if q.is_dir() else q.unlink()
    workdir.mkdir(parents=True, exist_ok=True)
    marker.write_text("created by validation/notebook_on_published_build.py\n",
                      encoding="utf-8")

    staged = {}
    for pair in args.input:
        dest, _, src = pair.partition("=")
        source = Path(src).resolve()
        if not source.is_file():
            print(f"declared input missing: {source}")
            return 2
        target = (workdir / dest).resolve()
        if workdir not in target.parents:
            # `..` escapes and absolute destinations would stage outside the sandbox.
            print(f"staged input destination escapes the working directory: {dest}")
            return 2
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(source, target)
        staged[dest] = {"from": str(source), "sha256": _sha256(source)}

    # Before the build is imported, so a cache override is in place when clippr resolves
    # its paths. Without this the installed wheel read a user-level cache whose bytes no
    # record named: the designer notebook consumed a codon table and a genome from
    # %LOCALAPPDATA%\clippr and declared no inputs at all.
    environment = {}
    for pair in args.env:
        name, _, value = pair.partition("=")
        if not name:
            print(f"--env needs NAME=VALUE, got {pair!r}")
            return 2
        resolved = Path(value)
        if not resolved.is_absolute():
            resolved = (workdir / value).resolve()
            if workdir not in resolved.parents and resolved != workdir:
                print(f"--env {name} resolves outside the working directory: {value}")
                return 2
            value = str(resolved)
        os.environ[name] = value
        environment[name] = value

    unstaged = [d for d in args.require_input if d not in staged]
    if unstaged:
        print(f"required input not staged: {', '.join(unstaged)}. Pass it with "
              f"--input DEST=SRC; a run without its declared inputs is not a run.")
        return 2

    # The working tree must not be importable, and the repository must not be the cwd.
    for entry in list(sys.path):
        if entry and Path(entry).resolve() in (REPO / "src", REPO):
            sys.path.remove(entry)
    os.chdir(workdir)

    try:
        import clippr
    except Exception as exc:                        # noqa: BLE001 - reported, not hidden
        print(f"cannot import clippr: {type(exc).__name__}: {exc}")
        return 1
    origin = Path(clippr.__file__).resolve()
    prefix = Path(sys.prefix).resolve()
    inside = prefix == origin or prefix in origin.parents
    print(f"clippr imported from {origin}")
    print(f"resolves inside this environment ({sys.prefix}): {inside}")
    if not inside:
        print("refusing to report on a build that is not the installed one")
        return 1

    # The installed package must be the source under review. A green notebook run against a
    # wheel built before the repairs proves nothing about the repairs -- which is exactly
    # what happened once: nine files differed and the run still passed.
    installed_fp = source_fp = None
    try:
        from clippr.manifest import source_fingerprint
        installed_fp = source_fingerprint()["sha256"]
        checkout = REPO / "src" / "clippr" / "manifest.py"
        if checkout.is_file():
            import importlib.util
            spec = importlib.util.spec_from_file_location("_checkout_manifest", checkout)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            source_fp = module.source_fingerprint()["sha256"]
    except Exception as exc:                        # noqa: BLE001 - reported, not hidden
        print(f"cannot compare installed and source fingerprints: "
              f"{type(exc).__name__}: {exc}")
        return 1
    if source_fp and installed_fp != source_fp:
        print(f"installed package {installed_fp[:12]} does not match the source "
              f"{source_fp[:12]}; rebuild the wheel before running notebooks")
        return 1
    print(f"installed package fingerprint matches the source: {str(installed_fp)[:12]}")

    cells = json.loads(notebook.read_text(encoding="utf-8"))["cells"]
    env: dict = {"__name__": "__main__"}
    ran = failed = 0
    records = []

    for i, cell in enumerate(cells):
        if cell["cell_type"] != "code":
            continue
        src = "".join(cell["source"])
        title = next((ln for ln in src.splitlines() if ln.startswith("#@title")), "")
        label = title.replace("#@title ", "").split("{")[0].strip() or f"cell {i}"
        # Classify line by line. Skipping the whole cell whenever it held a permitted
        # magic is what let required work disappear: the fix for `%time` left the
        # mechanism in place, so `%%capture` plus a raise still reported PASS.
        keep, dropped, unsupported = classify_cell(src)
        if unsupported:
            failed += 1
            records.append({"cell": i, "label": label, "status": "failed",
                            "error": ("unsupported magic in a required cell: "
                                      + unsupported[0][:60])})
            print(f"  [FAIL] {label}")
            print(f"         unsupported magic: {unsupported[0][:60]}")
            continue
        body = "\n".join(keep)
        if dropped and not _has_statement(body):
            # Nothing but permitted magic lines: a genuine install cell.
            records.append({"cell": i, "label": label, "status": "skipped",
                            "reason": f"install line only: {dropped[0][:60]}"})
            print(f"  [skip] {label}  ({dropped[0][:40]})")
            continue
        try:
            exec(compile(body, f"<cell {i}>", "exec"), env)
            ran += 1
            records.append({"cell": i, "label": label, "status": "ok",
                            **({"dropped_lines": dropped} if dropped else {})})
            print(f"  [ ok ] {label}")
        except Exception as exc:                    # noqa: BLE001 - reported, not hidden
            failed += 1
            detail = "".join(traceback.format_exception_only(type(exc), exc)).strip()
            records.append({"cell": i, "label": label, "status": "failed",
                            "error": detail})
            print(f"  [FAIL] {label}\n         {detail}")

    # Postconditions: did the notebook do the work, or only execute cleanly? Evaluated
    # in the notebook's own namespace, so an expectation naming a result the notebook
    # never built fails on the NameError rather than being quietly absent.
    expectations = []
    for expr in args.expect:
        if failed:
            expectations.append({"expect": expr, "ok": False,
                                 "error": "not evaluated: a cell had already failed"})
            continue
        try:
            value = eval(expr, env)                 # noqa: S307 - operator-supplied
            held = bool(value)
            entry = {"expect": expr, "ok": held, "value": repr(value)[:200]}
        except Exception as exc:                    # noqa: BLE001 - reported, not hidden
            held = False
            entry = {"expect": expr, "ok": False,
                     "error": f"{type(exc).__name__}: {exc}"}
        expectations.append(entry)
        print(f"  [{'  ok ' if held else 'FAIL'}] expect: {expr[:70]}")
        if not held and "error" in entry:
            print(f"         {entry['error'][:100]}")
    unmet = [e for e in expectations if not e["ok"]]
    if unmet:
        failed += len(unmet)

    skipped = [r for r in records if r["status"] == "skipped"]
    report = {
        "notebook": str(notebook),
        "notebook_sha256": _sha256(notebook),
        "interpreter": sys.executable,
        "environment_prefix": sys.prefix,
        "import_origin": str(origin),
        "import_resolves_inside_environment": inside,
        "installed_source_fingerprint": installed_fp,
        "checkout_source_fingerprint": source_fp,
        "workdir": str(workdir),
        "staged_inputs": staged,
        "required_inputs": list(args.require_input),
        "environment_overrides": environment,
        "expectations": expectations,
        "expected_code_cells": ran + failed + len(skipped),
        "executed": ran,
        "failed": failed,
        "skipped": skipped,
        "cells": records,
        # Executing cleanly is not the same as doing the work. A notebook with no
        # declared postcondition can still only report what it is: cells that ran.
        "ok": failed == 0 and ran > 0 and not unmet,
        "postconditions_declared": len(expectations),
        "scope": ("installed local notebook execution; no browser or Colab run is implied"),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    if ran == 0:
        # A notebook that executed nothing has not been shown to work.
        failed = failed or 1
        print("no code cell executed; this is not a passing run")
    print(f"\n{ran} ran, {failed} failed, {len(skipped)} skipped (named above)")
    print(f"wrote {out}")
    print("PASS -- the installed build runs this notebook" if not failed
          else "FAIL -- this notebook does not run on the installed build")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
