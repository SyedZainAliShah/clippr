"""Execute every code cell of the Colab notebook, in order, as a user would.

`tests/test_notebook.py` checks the notebook statically -- that it parses, that its dropdowns
resolve, that its imports still exist. That is fast and runs on every commit. This is the other
half: actually running the thing, which catches what static analysis cannot, such as a call that
type-checks but raises, or a cell that depends on a variable an earlier cell no longer defines.

The install cell is skipped, since the package is already importable when this runs. Every other
cell executes for real against live code, in one shared namespace, exactly as "Run all" would.

Two dependencies are Colab-provided rather than ours: IPython (for `display` and the file
download helper) and jinja2 (which pandas' `.style` needs). Both are absent from a bare
environment. The notebook is written to degrade gracefully without jinja2, and this script
reports rather than fails when IPython is missing, so a checkout without them still gets a
useful answer instead of a misleading red.

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\run_notebook.py

A notebook and its form values can both be chosen, which is what makes this a check of
the route a user actually takes rather than of the defaults alone. `--set` rewrites a
`#@param` assignment before the cell runs, exactly as typing in the form does:

    ... validation\\run_notebook.py --notebook notebooks\\PPROPT_inventory.ipynb \\
        --set organism=e_coli --set label=ecoli-demo
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

NOTEBOOK = ROOT / "notebooks" / "PPROPT_designer.ipynb"

#: The inventory cells short-circuit when this file is absent, and the absent branch is the
#: cheap one. Running only that branch is what let the inventory half of the notebook ship
#: having never executed: `table` was passed as a codon table and is a display DataFrame, so
#: three cells raised the moment a real input reached them, while every test stayed green.
#: The harness now stages the file itself, so the path a user takes is the path that runs.
SUPPLIED_TABLE_S1 = ROOT / "data" / "grasp_supp" / "Table S1.xlsx"


def _literal(value: str):
    """`24` as an int, `True` as a bool, anything else as the string it came in as.

    A form field declared `{type:"integer"}` hands the cell an int, so quoting everything
    makes `max_evaluations` a string and the first comparison against it raises
    `'>=' not supported between instances of 'int' and 'str'` -- a harness artefact that
    looks exactly like a notebook bug.
    """
    import ast

    try:
        parsed = ast.literal_eval(value)
    except (ValueError, SyntaxError):
        return value
    return parsed if isinstance(parsed, (int, float, bool)) else value


def _override(src: str, settings: dict[str, str]) -> str:
    """Rewrite `name = <literal>` to the requested value, keeping the `#@param` comment.

    Colab's form widgets do exactly this: the comment declares the widget, the assignment
    holds the value. Setting the variable after the cell has run would be too late, because
    the cell uses it in the same breath.
    """
    for name, value in settings.items():
        pattern = re.compile("^" + re.escape(name) + " = [^" + chr(10) + "]*$", re.MULTILINE)
        found = pattern.search(src)
        if not found:
            continue
        comment = ""
        if "#@param" in found.group(0):
            comment = "  #@param" + found.group(0).split("#@param", 1)[1]
        src = pattern.sub(lambda _m: f"{name} = {_literal(value)!r}{comment}", src, count=1)
    return src


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Execute every code cell of a PPROPT notebook, as a user would.")
    parser.add_argument("--notebook", default=str(NOTEBOOK),
                        help="notebook to execute (default: the designer notebook)")
    parser.add_argument("--set", action="append", default=[], metavar="NAME=VALUE",
                        help="override a #@param value, as typing in the form would")
    args = parser.parse_args()
    settings = dict(s.split("=", 1) for s in args.set)

    try:
        import IPython  # noqa: F401
    except ImportError:
        print("IPython is not installed here, so cells that call display() or the Colab")
        print("download helper cannot run. Colab always provides it. Install it to run this")
        print("check locally:  pip install ipython")
        print("\nSKIP")
        raise SystemExit(0)

    # The notebook looks for the workbook by bare filename, as a Colab user would after
    # uploading it. Stage it in the working directory and take it away afterwards.
    staged = Path.cwd() / SUPPLIED_TABLE_S1.name
    staged_by_us = False
    if SUPPLIED_TABLE_S1.is_file() and not staged.exists():
        staged.write_bytes(SUPPLIED_TABLE_S1.read_bytes())
        staged_by_us = True
        print(f"staged {staged.name} so the inventory cells run for real")
    elif not SUPPLIED_TABLE_S1.is_file():
        print("Supplementary Table S1 is not present, so the inventory cells will only")
        print("exercise their absent-input branch. That is a weaker check, and it is the")
        print("one that previously hid a crash on the real path.")

    try:
        _run(staged, Path(args.notebook), settings)
    finally:
        if staged_by_us:
            staged.unlink(missing_ok=True)


def _run(staged: Path, notebook: Path, settings: dict[str, str]) -> None:
    nb = json.loads(notebook.read_text(encoding="utf-8"))
    cells = [c for c in nb["cells"] if c["cell_type"] == "code"]
    print(f"{len(cells)} code cells in {notebook.name}")
    if settings:
        print("form values: " + ", ".join(f"{k}={v}" for k, v in settings.items()))
    print()

    env: dict = {"__name__": "__main__"}
    failures: list[tuple[int, str, str]] = []
    ran = 0

    for i, cell in enumerate(cells, 1):
        src = _override("".join(cell["source"]), settings)
        title = next((ln for ln in src.splitlines() if ln.startswith("#@title")), "")
        label = title.replace("#@title", "").split("{")[0].strip() or f"cell {i}"
        if any(ln.lstrip().startswith(("%", "!")) for ln in src.splitlines()):
            print(f"[{i:2}/{len(cells)}] skip  {label}  (installs the package)")
            continue
        try:
            exec(compile(src, f"<notebook cell {i}>", "exec"), env)
            ran += 1
            print(f"[{i:2}/{len(cells)}] ok    {label}", flush=True)
        except Exception as exc:                       # any failure is a broken notebook
            failures.append((i, label, f"{type(exc).__name__}: {exc}"))
            print(f"[{i:2}/{len(cells)}] FAIL  {label}", flush=True)
            traceback.print_exc(limit=4)

    print(f"\n{ran} cells executed, {len(failures)} failed")
    for i, label, err in failures:
        print(f"  cell {i}  {label}: {err}")
    print("\n" + ("PASS" if not failures else "FAIL"))
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
