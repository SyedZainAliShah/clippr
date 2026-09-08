"""Run every validation script and print one table.

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\report.py [--quick]

Checks read saved artifacts and run in this package's own environment.

Some validation scripts are kept locally rather than distributed -- they depend on an
external oracle that not every checkout has. They are picked up automatically when
present and simply absent otherwise, so this report never silently drops a check it
could have run. Set `CLIPPR_ORACLE_PYTHON` to the interpreter those scripts need.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUR_PY = ROOT / ".venv" / "Scripts" / "python.exe"

#: Interpreter for checks that need an external oracle, if this checkout has one.
ORACLE_PY = Path(os.environ["CLIPPR_ORACLE_PYTHON"]) if os.environ.get(
    "CLIPPR_ORACLE_PYTHON") else None

#: Checks kept outside the distributed repository, run only when the file is present.
OPTIONAL = {"compare_ppr.py": ("protein sequence vs an external oracle", True, [])}

#: (script, what it proves, needs the oracle interpreter, args when --quick)
CHECKS = [
    ("compare_cuts.py", "cut geometry over 200 designs", False, []),
    ("compare_assembly.py", "oligos byte-identical", False, []),
    ("compare_export.py", "quotes, GenBank, FASTA", False, []),
    ("compare_codons.py", "constraints satisfied", False, ["--limit", "20"]),
    ("calibrate_qc.py", "QC verdict separates", False, []),
    ("integration.py", "end-to-end, 200 targets", False, ["--limit", "20"]),
]


def all_checks() -> list[tuple[str, str, bool, list[str]]]:
    """The distributed checks, plus any optional local ones this checkout happens to have."""
    found = [(name, *spec) for name, spec in OPTIONAL.items()
             if (ROOT / "validation" / name).exists()]
    return found + CHECKS


def run(script: str, python: Path | None, args: list[str]) -> tuple[str, float, str]:
    if python is None:
        return "SKIP", 0.0, "set CLIPPR_ORACLE_PYTHON to run this check"
    if not python.exists():
        return "SKIP", 0.0, f"{python.name} not found"
    t0 = time.time()
    p = subprocess.run([str(python), str(ROOT / "validation" / script), *args],
                       capture_output=True, text=True, env={"PYTHONUTF8": "1", **_env()})
    dt = time.time() - t0
    tail = [l for l in (p.stdout or "").strip().splitlines() if l.strip()]
    detail = tail[-1] if tail else (p.stderr or "").strip().splitlines()[-1:] or [""]
    if isinstance(detail, list):
        detail = detail[0] if detail else ""
    return ("PASS" if p.returncode == 0 else "FAIL"), dt, detail


def _env() -> dict:

    return {k: v for k, v in os.environ.items()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="shorten the slow checks; not a substitute for a full run")
    args = ap.parse_args()

    print(f"{'check':24s}{'proves':34s}{'result':>8s}{'sec':>8s}")
    print("-" * 74)
    results = []
    for script, proves, needs_oracle, quick_args in all_checks():
        python = ORACLE_PY if needs_oracle else OUR_PY
        status, dt, detail = run(script, python, quick_args if args.quick else [])
        results.append((script, status))
        print(f"{script:24s}{proves:34s}{status:>8s}{dt:>8.1f}")
        if status == "FAIL":
            print(f"    {detail}")
        elif status == "SKIP":
            print(f"    skipped: {detail}")

    failed = [s for s, st in results if st == "FAIL"]
    skipped = [s for s, st in results if st == "SKIP"]
    print("-" * 74)
    print(f"{len(results) - len(failed) - len(skipped)} passed, "
          f"{len(failed)} failed, {len(skipped)} skipped")
    if args.quick:
        print("--quick shortened compare_codons and integration; run without it before "
              "calling the build done.")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
