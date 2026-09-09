"""Run every validation script and print one table, in two clearly separated tiers.

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\report.py [--quick]

**The two tiers are not the same kind of evidence, and the split is the point.**

*Reproducible anywhere* checks run on a bare checkout with no extra data. Their inputs are
generated from fixed seeds, so anyone who clones this repository can re-run them and get the
same answer. This is the evidence a reader can actually verify.

*Development record* checks compare our output against a 200-design reference corpus. That
corpus is third-party program output, 38 MB, and is not redistributed, so these checks skip
on any checkout that lacks it. They are the stronger claim -- byte-level parity -- but a
reader cannot confirm them, and this report must never let a machine that happens to hold
the corpus imply otherwise.

Set `CLIPPR_ORACLE_PYTHON` for the few checks that additionally need an external
interpreter; they are absent from most checkouts and are skipped rather than dropped.
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

#: The reference corpus. Absent from every checkout but the ones it was built on.
CORPUS = ROOT / "results" / "_designs"

#: Checks kept outside the distributed repository, run only when the file is present.
OPTIONAL = {"compare_ppr.py": ("protein sequence vs an external oracle", "oracle", [])}

#: (script, what it proves, what it needs, args when --quick). `needs` is None for checks
#: that run on a bare checkout, "corpus" for the parity checks, "oracle" for an external
#: interpreter.
CHECKS = [
    ("check_regression.py", "our own output has not drifted", None, ["--limit", "20"]),
    ("integration.py", "end-to-end, 200 targets", None, ["--limit", "20"]),
    ("benchmark_search.py", "exploring beats not exploring", None, ["--targets", "6"]),
    ("run_notebook.py", "every notebook cell executes", None, []),
    ("compare_cuts.py", "cut geometry over 200 designs", "corpus", []),
    ("compare_assembly.py", "oligos byte-identical", "corpus", []),
    ("compare_export.py", "quotes, GenBank, FASTA", "corpus", []),
    ("compare_codons.py", "constraints satisfied", "corpus", ["--limit", "20"]),
    ("calibrate_qc.py", "QC verdict separates", "corpus", []),
]

TIERS = {None: "reproducible on any checkout",
         "corpus": "development record — needs the reference corpus",
         "oracle": "development record — needs an external oracle"}


def all_checks() -> list[tuple[str, str, str | None, list[str]]]:
    """The distributed checks, plus any optional local ones this checkout happens to have."""
    found = [(name, *spec) for name, spec in OPTIONAL.items()
             if (ROOT / "validation" / name).exists()]
    return CHECKS + found


def run(script: str, python: Path | None, args: list[str],
        needs: str | None = None) -> tuple[str, float, str]:
    if needs == "corpus" and not CORPUS.exists():
        return "SKIP", 0.0, "reference corpus not present in this checkout"
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

    results = []
    for tier, heading in TIERS.items():
        group = [c for c in all_checks() if c[2] == tier]
        if not group:
            continue
        print(f"\n{heading.upper()}")
        print(f"{'check':24s}{'proves':34s}{'result':>8s}{'sec':>8s}")
        print("-" * 74)
        for script, proves, needs, quick_args in group:
            python = ORACLE_PY if needs == "oracle" else OUR_PY
            status, dt, detail = run(script, python,
                                     quick_args if args.quick else [], needs)
            results.append((script, status, tier))
            print(f"{script:24s}{proves:34s}{status:>8s}{dt:>8.1f}")
            if status == "FAIL":
                print(f"    {detail}")
            elif status == "SKIP":
                print(f"    skipped: {detail}")

    failed = [s for s, st, _ in results if st == "FAIL"]
    skipped = [s for s, st, _ in results if st == "SKIP"]
    here = [s for s, st, t in results if st == "PASS" and t is None]
    parity_ran = [s for s, st, t in results if st == "PASS" and t is not None]
    print("\n" + "=" * 74)
    print(f"{len(results) - len(failed) - len(skipped)} passed, "
          f"{len(failed)} failed, {len(skipped)} skipped")
    print(f"verifiable on a bare checkout: {len(here)} of "
          f"{sum(1 for c in all_checks() if c[2] is None)}")
    if parity_ran:
        print(f"{len(parity_ran)} parity checks ran because this machine holds the reference "
              "data. They are a development record, not evidence a reader can confirm — "
              "cite them as such.")
    if args.quick:
        print("--quick shortened compare_codons, integration and check_regression; run "
              "without it before calling the build done.")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
