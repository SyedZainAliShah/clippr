"""The gate. Deterministic, non-negotiable, and not a model.

**Why this exists.** Every mistake this project has shipped was a mistake of *judgement*
under pressure, not of capability: a suite run on a subset and reported as the whole, a
`file:line` citation that pointed at the wrong line, a heredoc that mangled an escape, a
notebook whose builder changed without the notebook being rebuilt. Each one was caught, late,
by a human reading output. None was caught by prose telling someone to be careful.

The lesson this repository already learned, in `CLAUDE.md`, is that **every rule framed as
"reject if X" was followed and every rule framed as "remember to do X" was ignored.** So this
file contains no advice. It runs checks and it exits non-zero.

**On tamper-resistance, honestly.** No file on a disk an agent can write to is immutable, and
claiming otherwise would be the sort of overstatement this gate exists to prevent. What is
achievable, and what this provides, is two things:

1. **Determinism.** The verdict is an exit code computed from measurements. There is no
   wording to soften and no reviewer to persuade. An agent cannot talk a count of failures
   down to zero.
2. **Tamper-evidence.** `--self` prints this file's SHA-256. Record it before a task and
   compare after; `git diff validation/gate.py` shows any edit. Agents are instructed never
   to modify this file, and a change to it during a task is itself the finding.

Usage:

    python validation/gate.py                 # everything, including the full suite
    python validation/gate.py --deep          # also execute both notebooks, cell by cell
    python validation/gate.py --quick         # everything except the suite, for iteration
    python validation/gate.py --self          # print this file's SHA-256 and exit
    python validation/gate.py --json out.json # also write a machine-readable verdict

Exit code is 0 only when every check passes. `--quick` can never exit 0 on its own terms:
it reports `INCOMPLETE`, because a gate that has not run the suite has not gated anything.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYTHON = Path(sys.executable)

#: Phrases this repository has banned, with the reason, so a failure explains itself.
#: Sourced from CLAUDE.md. Each is matched case-insensitively against changed markdown.
BANNED = {
    "orderable piece": "a fragment is not an order; say 'ordered fragment' or 'oligo'",
    "the capacity is 84": "say '84 encodings were obtained'; capacity is not what was measured",
    "vendor approved": "no vendor has approved anything; we never submitted a sequence",
    "vendor accepted": "no vendor has accepted anything; we never submitted a sequence",
}

#: Directories whose every .py file must parse. Catches heredoc escape mangling, which has
#: silently corrupted a generator twice.
PARSE_DIRS = ("src", "tools", "validation", "tests")

#: The notebooks, and the generator that must reproduce them byte for byte.
NOTEBOOKS = {
    "notebooks/PPROPT_designer.ipynb": "tools/build_notebook.py",
    "notebooks/PPROPT_inventory.ipynb": "tools/build_inventory_notebook.py",
}


class Check:
    """One verdict. `ok is None` means the check could not run, which is not a pass."""

    def __init__(self, name: str, ok: bool | None, detail: str, data: dict | None = None):
        self.name, self.ok, self.detail, self.data = name, ok, detail, data or {}

    def line(self) -> str:
        mark = {True: "PASS", False: "FAIL", None: "SKIP"}[self.ok]
        return f"  [{mark}] {self.name}: {self.detail}"

    def as_dict(self) -> dict:
        return {"check": self.name, "ok": self.ok, "detail": self.detail, **self.data}


def _run(args: list[str], **kw) -> subprocess.CompletedProcess:
    env = {**os.environ, "PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8"}
    return subprocess.run(args, cwd=str(ROOT), capture_output=True, text=True,
                          encoding="utf-8", errors="replace", env=env, **kw)


# --------------------------------------------------------------------------- checks


def check_suite() -> Check:
    """The full test suite, counted from JUnit XML rather than stdout.

    Not from stdout on purpose. pytest's summary line does not reliably reach a redirected
    stdout in this environment, which is exactly how a run that never finished was once
    reported as passing. The XML is written by pytest itself and carries the counts as
    attributes, so there is nothing to parse out of prose and nothing to mis-read.
    """
    with tempfile.TemporaryDirectory() as tmp:
        xml = Path(tmp) / "junit.xml"
        started = time.time()
        proc = _run([str(PYTHON), "-m", "pytest", "tests/", "-q",
                     f"--junitxml={xml}"])
        elapsed = time.time() - started
        if not xml.is_file():
            return Check("suite", False,
                         f"pytest wrote no XML after {elapsed:.0f}s; it did not complete "
                         f"(exit {proc.returncode})")
        root = ET.parse(xml).getroot()
        s = root[0] if root.tag == "testsuites" else root
        a = s.attrib
        total = int(a["tests"])
        failures, errors = int(a.get("failures", 0)), int(a.get("errors", 0))
        skipped = int(a.get("skipped", 0))
        passed = total - skipped - failures - errors
        broken = [f'{tc.get("classname")}::{tc.get("name")}'
                  for tc in s for k in ("failure", "error") if tc.findall(k)]
        detail = (f"{total} tests, {passed} passed, {skipped} skipped, "
                  f"{failures} failures, {errors} errors, {elapsed:.0f}s")
        if failures or errors:
            detail += " -- " + ", ".join(broken[:5])
            if len(broken) > 5:
                detail += f" (+{len(broken) - 5} more)"
        return Check("suite", failures == 0 and errors == 0, detail,
                     {"tests": total, "passed": passed, "skipped": skipped,
                      "failures": failures, "errors": errors, "broken": broken})


def check_python_parses() -> Check:
    """Every tracked .py compiles. Catches a generator corrupted by shell escaping."""
    bad = []
    for d in PARSE_DIRS:
        for path in sorted((ROOT / d).rglob("*.py")):
            if "__pycache__" in path.parts or ".venv" in path.parts:
                continue
            try:
                ast.parse(path.read_text(encoding="utf-8"))
            except (SyntaxError, UnicodeDecodeError) as exc:
                bad.append(f"{path.relative_to(ROOT)}: {exc}")
    return Check("python parses", not bad,
                 "every tracked .py parses" if not bad else "; ".join(bad[:3]),
                 {"unparseable": bad})


def check_notebooks_current() -> Check:
    """Each notebook is exactly what its generator produces, right now.

    A generator edited without a rebuild ships a notebook that disagrees with the code that
    claims to produce it, and no test of either half notices.
    """
    stale, missing = [], []
    for nb, gen in NOTEBOOKS.items():
        nb_path, gen_path = ROOT / nb, ROOT / gen
        if not nb_path.is_file() or not gen_path.is_file():
            missing.append(nb if not nb_path.is_file() else gen)
            continue
        before = nb_path.read_bytes()
        proc = _run([str(PYTHON), str(gen_path)])
        after = nb_path.read_bytes()
        if proc.returncode != 0:
            stale.append(f"{gen} failed to run: {(proc.stderr or '')[-160:]}")
        elif before != after:
            stale.append(f"{nb} was not current; {gen} rebuilt it")
    if missing:
        return Check("notebooks current", False, f"missing: {missing}")
    return Check("notebooks current", not stale,
                 "both notebooks match their generators" if not stale
                 else "; ".join(stale), {"stale": stale})


def check_notebook_line_endings() -> Check:
    """Notebooks are LF. `.gitattributes` marks *.ipynb -text, so CRLF breaks CI."""
    bad = [nb for nb in NOTEBOOKS if b"\r\n" in (ROOT / nb).read_bytes()]
    return Check("notebook line endings", not bad,
                 "LF throughout" if not bad else f"CRLF found in {bad}", {"crlf": bad})


def _changed_markdown() -> list[Path]:
    """Markdown touched relative to origin/main, falling back to the working tree."""
    for args in (["git", "diff", "--name-only", "origin/main...HEAD"],
                 ["git", "diff", "--name-only", "HEAD"]):
        proc = _run(args)
        if proc.returncode == 0 and proc.stdout.strip():
            names = [n for n in proc.stdout.split("\n") if n.strip().endswith(".md")]
            if names:
                return [ROOT / n for n in names if (ROOT / n).is_file()]
    return []


#: `path/to/file.py:123` inside prose or a markdown link.
_CITATION = re.compile(r"\b((?:[\w.\-]+/)*[\w.\-]+\.(?:py|json|yml|yaml|ipynb)):(\d+)\b")


def check_citations() -> Check:
    """Every `file:line` in changed markdown points at a line that exists.

    This catches the mechanical half of a bad citation: a file that is not there, or a line
    number past the end of it. It cannot check that the line *says* what the prose claims --
    that is a reader's job, and the researcher agent's. Three of thirteen citations in one
    document were off by a line or two, so the mechanical half is worth automating.
    """
    problems, checked = [], 0
    for md in _changed_markdown():
        text = md.read_text(encoding="utf-8", errors="replace")
        for raw, lineno in _CITATION.findall(text):
            target = ROOT / raw
            if not target.is_file():
                for base in ("src/ppropt", "src", "tools", "validation", "tests"):
                    if (ROOT / base / raw).is_file():
                        target = ROOT / base / raw
                        break
            if not target.is_file():
                continue  # not a repo path; prose may legitimately name a foreign file
            checked += 1
            n = len(target.read_text(encoding="utf-8", errors="replace").splitlines())
            if int(lineno) > n or int(lineno) < 1:
                problems.append(f"{md.name} cites {raw}:{lineno} but it has {n} lines")
    return Check("citations resolve", not problems,
                 f"{checked} file:line citations resolve" if not problems
                 else "; ".join(problems[:4]), {"checked": checked, "problems": problems})


def check_banned_phrases() -> Check:
    """Phrases CLAUDE.md forbids, checked in changed markdown."""
    hits = []
    for md in _changed_markdown():
        lowered = md.read_text(encoding="utf-8", errors="replace").lower()
        for phrase, why in BANNED.items():
            if phrase in lowered:
                hits.append(f"{md.name}: '{phrase}' -- {why}")
    return Check("banned phrases", not hits,
                 "none present" if not hits else "; ".join(hits), {"hits": hits})


def check_notebooks_execute() -> Check:
    """Both notebooks run every cell in order against live code.

    Catches a cell that type-checks and raises: the inventory notebook once shipped having
    never executed, because a display DataFrame was passed where a codon table was expected
    and three cells raised the moment a real input reached them, while every static test
    stayed green.

    **What it does not catch**, stated so nobody relies on it for this: a name defined and
    used within one cell is fine here but can still break a *slice* of that cell. The
    `chosen_table_source` regression was exactly that shape -- defined beside the design
    call, used in the result card, both in cell 7 -- so running cells in order never saw it.
    `tests/test_display_verdict.py` caught it, by rendering the card from the cell's tail in
    a narrower namespace. That test is the guard for that class, not this check.

    Slow -- several minutes per notebook -- so it runs only under `--deep`.
    """
    broken, ran = [], []
    for nb in NOTEBOOKS:
        proc = _run([str(PYTHON), str(ROOT / "validation" / "run_notebook.py"),
                     "--notebook", str(ROOT / nb)])
        tail = (proc.stdout or "").strip().splitlines()
        verdict = next((ln for ln in reversed(tail)
                        if ln.strip() in ("PASS", "FAIL", "SKIP")), "")
        counts = next((ln for ln in reversed(tail) if "cells executed" in ln), "")
        name = Path(nb).name
        if verdict == "PASS":
            ran.append(f"{name} ({counts})")
        elif verdict == "SKIP":
            broken.append(f"{name}: skipped -- {tail[0] if tail else 'no reason given'}")
        else:
            broken.append(f"{name}: {counts or (proc.stderr or '')[-160:]}")
    return Check("notebooks execute", not broken,
                 "; ".join(ran) if not broken else "; ".join(broken),
                 {"ran": ran, "broken": broken})


def check_tree_state() -> Check:
    """Reports what is uncommitted. Advisory: a dirty tree is normal mid-task."""
    proc = _run(["git", "status", "--porcelain"])
    lines = [ln for ln in proc.stdout.split("\n") if ln.strip()]
    tracked = [ln for ln in lines if not ln.startswith("??")]
    return Check("working tree", True,
                 "clean" if not lines else
                 f"{len(tracked)} tracked change(s), {len(lines) - len(tracked)} untracked",
                 {"tracked": tracked})


# --------------------------------------------------------------------------- driver


def self_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true",
                    help="skip the suite; reports INCOMPLETE and never exits 0")
    ap.add_argument("--deep", action="store_true",
                    help="also execute both notebooks cell by cell; adds several minutes")
    ap.add_argument("--self", action="store_true", help="print this file's SHA-256 and exit")
    ap.add_argument("--json", metavar="PATH", help="also write the verdict as JSON")
    args = ap.parse_args()

    if getattr(args, "self"):
        print(self_sha256())
        return

    print(f"gate.py {self_sha256()[:16]}  root {ROOT}")
    print()

    checks = [
        check_python_parses(),
        check_notebooks_current(),
        check_notebook_line_endings(),
        check_citations(),
        check_banned_phrases(),
        check_tree_state(),
    ]
    if args.deep:
        print("  executing both notebooks ...", flush=True)
        checks.append(check_notebooks_execute())
    if args.quick:
        checks.append(Check("suite", None, "skipped by --quick"))
    else:
        print("  running the full suite; this takes several minutes ...", flush=True)
        checks.append(check_suite())

    for c in checks:
        print(c.line())

    failed = [c for c in checks if c.ok is False]
    incomplete = [c for c in checks if c.ok is None]
    verdict = "FAIL" if failed else ("INCOMPLETE" if incomplete else "PASS")

    print()
    print(f"VERDICT: {verdict}")
    if failed:
        print(f"  {len(failed)} check(s) failed: {', '.join(c.name for c in failed)}")
    if incomplete:
        print(f"  {len(incomplete)} check(s) did not run: "
              f"{', '.join(c.name for c in incomplete)}")
        print("  INCOMPLETE is not a pass. Nothing ships on it.")

    if args.json:
        Path(args.json).write_text(json.dumps({
            "verdict": verdict,
            "gate_sha256": self_sha256(),
            "checks": [c.as_dict() for c in checks],
        }, indent=1) + "\n", encoding="utf-8")

    raise SystemExit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
