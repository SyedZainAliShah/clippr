"""Run the notebook against an *installed* CLIPPR, the way Colab actually does.

`validation/run_notebook.py` puts the working tree's `src/` first on `sys.path`, so it proves
the notebook agrees with the code being edited. It cannot catch the failure that matters to a
user: the notebook depending on something the **published** package does not have.

That failure happened. A notebook cell was added that imported `clippr.synthesis_profile` and
passed `profile=` to four workflow calls. The module exists on a feature branch; the notebook's
own Setup cell installs from the repository's default branch, where it does not. Every Colab
session broke, and `run_notebook.py` stayed green throughout, because it was testing the
working tree the cell had been written against.

So this imports nothing from the working tree. Run it with the interpreter of a virtualenv that
has CLIPPR pip-installed from GitHub, and it executes the notebook's code cells against that
build:

    python -m venv /tmp/sim
    /tmp/sim/bin/python -m pip install ipython "git+https://github.com/SyedZainAliShah/clippr.git"
    /tmp/sim/bin/python validation/notebook_on_published_build.py . notebooks/CLIPPR_designer.ipynb

`ipython` is installed because Colab has it and some cells use `display()`; without it the
harness reports failures that a real Colab would not have.

Cells beginning with a `%` or `!` line are skipped, as the in-tree harness does -- those are
the install and upload magics, which have no meaning outside a notebook.

Exits non-zero if any cell fails.
"""
from __future__ import annotations

import json
import os
import sys
import traceback
from pathlib import Path

ROOT = Path(sys.argv[1]).resolve()
NOTEBOOK = Path(sys.argv[2]).resolve()

# The whole point: the working tree must not be importable.
for entry in list(sys.path):
    if entry and Path(entry).resolve() == (ROOT / "src"):
        sys.path.remove(entry)
os.chdir(ROOT)                      # data files resolve relative to the repo

import clippr                        # noqa: E402

print(f"clippr from {Path(clippr.__file__).parent}")
print(f"working tree src on path: "
      f"{any(Path(p).resolve() == ROOT / 'src' for p in sys.path if p)}")
print()

cells = json.loads(NOTEBOOK.read_text(encoding="utf-8"))["cells"]
env: dict = {"__name__": "__main__"}
ran = failed = skipped = 0
failures = []

for i, cell in enumerate(cells):
    if cell["cell_type"] != "code":
        continue
    src = "".join(cell["source"])
    title = next((ln for ln in src.splitlines() if ln.startswith("#@title")), "")
    label = title.replace("#@title ", "").split("{")[0].strip() or f"cell {i}"

    if any(ln.lstrip().startswith(("%", "!")) for ln in src.splitlines()):
        skipped += 1
        print(f"  [skip] {label}  (notebook magic)")
        continue

    try:
        exec(compile(src, f"<cell {i}>", "exec"), env)
        ran += 1
        print(f"  [ ok ] {label}")
    except Exception as exc:                       # noqa: BLE001 - reported, not hidden
        failed += 1
        failures.append((i, label, exc))
        print(f"  [FAIL] {label}")
        print("         " + "".join(
            traceback.format_exception_only(type(exc), exc)).strip())

print()
print(f"{ran} ran, {failed} failed, {skipped} skipped (magics)")
if failures:
    print("\nFAILURES:")
    for i, label, exc in failures:
        print(f"  cell {i}  {label}: {type(exc).__name__}: {exc}")
print("\n" + ("PASS -- the published build runs this notebook"
              if not failed else "FAIL -- this notebook does not run on the published build"))
sys.exit(1 if failed else 0)
