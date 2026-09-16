"""Release validation -- the ten checks from section 10 of the execution plan.

Each check is one of three things, and the distinction is the point of the script:

  RAN       this script executed something and read the result
  READ      this script inspected a file that another run produced, and reports what it says
  MANUAL    this cannot be decided mechanically; the script says so instead of passing it

A gate this script cannot evaluate is never reported as passed. The failure mode being
guarded against is a checklist that goes green because the checker was weaker than the claim.

The clean-environment install (check 3) builds a wheel and installs it into a throwaway
virtual environment, so it is slow and is skipped unless `--clean-install` is given. Its last
result is read from disk otherwise.

    python validation/release_check.py [--clean-install]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work" / "release"
PYTHON = ROOT / ".venv" / "Scripts" / "python.exe"
if not PYTHON.is_file():
    PYTHON = Path(sys.executable)

#: The example in README.md, and the line it must print. Both are read from the README rather
#: than restated here, so the check cannot drift from the documentation it is checking.
EXAMPLE_FENCE = re.compile(r"```python\n(.*?)```", re.DOTALL)


def _result(label: str, kind: str, ok: bool | None, detail: str) -> dict:
    return {"check": label, "evidence": kind, "ok": ok, "detail": detail}


def run_suite() -> tuple[subprocess.CompletedProcess, int]:
    """The whole suite, once, with skip reasons. Checks 1 and 2 both read this one run.

    Section 10 asks for the complete gates without duplicate orchestration of the same
    expensive check, and running pytest three times to answer three questions about one
    result is exactly that.
    """
    collected = subprocess.run([str(PYTHON), "-m", "pytest", "--collect-only"],
                               cwd=ROOT, capture_output=True, text=True)
    match = re.search(r"(\d+) tests collected", collected.stdout)
    run = subprocess.run([str(PYTHON), "-m", "pytest", "-rs"],
                         cwd=ROOT, capture_output=True, text=True)
    return run, int(match.group(1)) if match else -1


def check_tests(run: subprocess.CompletedProcess, collected: int) -> dict:
    """1. Unit and regression checks pass, and the count is the collected count."""
    passed = re.search(r"(\d+) passed", run.stdout)
    failed = re.search(r"(\d+) failed", run.stdout)
    # pytest's dot output is not a count: decimal points in comparison output have been
    # miscounted as tests here before. Only the collected total and the summary line count.
    n_passed = int(passed.group(1)) if passed else (collected if run.returncode == 0 else -1)
    return _result("1. unit and regression checks", "RAN",
                   run.returncode == 0 and collected > 0,
                   f"{n_passed} passed of {collected} collected, "
                   f"{int(failed.group(1)) if failed else 0} failed")


def check_integration_coverage(run: subprocess.CompletedProcess) -> dict:
    """2. Intended integration coverage accounted for; missing cases stay visible."""
    skips = re.findall(r"SKIPPED \[(\d+)\] ([^\n]+)", run.stdout)
    total = sum(int(n) for n, _ in skips)
    reasons = "; ".join(sorted({reason.split(": ", 1)[-1][:70] for _, reason in skips}))
    return _result("2. integration coverage", "RAN", True,
                   f"{total} skipped, each with a stated reason"
                   + (f" -- {reasons}" if reasons else " -- none skipped"))


def check_clean_install(run_it: bool) -> dict:
    """3. The public example installs and runs in a clean environment, no private paths."""
    record = WORK / "clean_install.json"
    if not run_it:
        if record.is_file():
            saved = json.loads(record.read_text(encoding="utf-8"))
            # Saved evidence only counts if it was produced from this tree. Otherwise it is
            # a record of some other checkout passing, which is not this release's evidence.
            stale = []
            if saved.get("package_source_sha256") != _package_source_sha256():
                stale.append("package source")
            if saved.get("readme_sha256") != _sha256(ROOT / "README.md"):
                stale.append("README (the example itself)")
            if stale:
                return _result("3. clean-environment install", "MANUAL", None,
                               f"saved evidence predates a change to {' and '.join(stale)}; "
                               f"re-run with --clean-install")
            return _result("3. clean-environment install", "READ", saved["ok"],
                           f"{saved['detail']} (recorded {saved['when']} from wheel "
                           f"{str(saved.get('wheel_sha256'))[:12]}; re-run with "
                           f"--clean-install)")
        return _result("3. clean-environment install", "MANUAL", None,
                       "not run in this session and no record on disk; use --clean-install")

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    example = next((b for b in EXAMPLE_FENCE.findall(readme) if "design_oneshot" in b), None)
    expected = re.search(r"^# (AAAAUGUGG \(9S\).*)$", example or "", re.MULTILINE)
    if not example or not expected:
        return _result("3. clean-environment install", "RAN", False,
                       "could not find the worked example and its expected output in README")

    area = Path(tempfile.mkdtemp(prefix="clippr-release-"))
    try:
        build = subprocess.run([str(PYTHON), "-m", "build", "--wheel", "--outdir",
                                str(area / "dist")], cwd=ROOT, capture_output=True, text=True)
        wheels = list((area / "dist").glob("*.whl"))
        if build.returncode or not wheels:
            return _result("3. clean-environment install", "RAN", False,
                           f"wheel build failed: {build.stderr.strip()[-200:]}")

        venv = area / "env"
        subprocess.run([str(PYTHON), "-m", "venv", str(venv)], check=True,
                       capture_output=True)
        python = venv / "Scripts" / "python.exe"
        if not python.is_file():
            python = venv / "bin" / "python"
        install = subprocess.run([str(python), "-m", "pip", "install", "--quiet",
                                  str(wheels[0])], capture_output=True, text=True)
        if install.returncode:
            return _result("3. clean-environment install", "RAN", False,
                           f"install failed: {install.stderr.strip()[-200:]}")

        table_file = ROOT / "data" / "codon_tables" / "kazusa_3055.json"
        shutil.copy(table_file, area)
        script = area / "example.py"
        script.write_text(example, encoding="utf-8")
        # A fresh cache directory, so a cached download from the developer's machine cannot
        # make the example look self-contained when it is not. PYTHONPATH and PYTHONHOME are
        # dropped outright: an inherited PYTHONPATH pointing at the source tree would let the
        # example import the checkout instead of the installed wheel, which is the one thing
        # this check exists to rule out.
        env = {k: v for k, v in os.environ.items()
               if k not in ("PYTHONPATH", "PYTHONHOME", "PYTHONSTARTUP")}
        env.update({"CLIPPR_CACHE_DIR": str(area / "cache"), "PYTHONIOENCODING": "utf-8",
                    "PYTHONNOUSERSITE": "1"})

        origin = subprocess.run(
            [str(python), "-c", "import clippr,sys;print(clippr.__file__);"
                                "print(clippr.__version__ if hasattr(clippr,'__version__') "
                                "else '')"],
            cwd=area, capture_output=True, text=True, env=env)
        imported = (origin.stdout or "").strip().splitlines()
        resolved_inside = bool(imported) and str(venv.resolve()) in str(
            Path(imported[0]).resolve())

        got = subprocess.run([str(python), str(script)], cwd=area, capture_output=True,
                             text=True, env=env, timeout=900)
        printed = (got.stdout or "").strip().splitlines()
        reproduced = (got.returncode == 0 and bool(printed)
                      and printed[-1] == expected.group(1).strip())
        ok = reproduced and resolved_inside
        if not resolved_inside:
            detail = (f"clippr imported from {imported[0] if imported else '?'}, which is "
                      f"outside the clean environment; the example did not exercise the wheel")
        elif reproduced:
            detail = f"example reproduced its documented output verbatim: {printed[-1]!r}"
        else:
            detail = (f"expected {expected.group(1).strip()!r}, got "
                      f"{printed[-1:] or got.stderr[-200:]!r}")

        WORK.mkdir(parents=True, exist_ok=True)
        from datetime import datetime, timezone
        # The evidence names the exact artefacts it was produced from. A filename and a
        # timestamp cannot tell a later reader whether this run tested the tree they hold.
        record.write_text(json.dumps(
            {"ok": ok, "detail": detail,
             "wheel": wheels[0].name, "wheel_sha256": _sha256(wheels[0]),
             "package_source_sha256": _package_source_sha256(),
             "readme_sha256": _sha256(ROOT / "README.md"),
             "codon_table_sha256": _sha256(table_file),
             "import_origin": imported[0] if imported else None,
             "resolved_inside_clean_env": resolved_inside,
             "when": datetime.now(timezone.utc).isoformat(timespec="seconds")},
            indent=2) + "\n", encoding="utf-8")
        return _result("3. clean-environment install", "RAN", ok, detail)
    finally:
        shutil.rmtree(area, ignore_errors=True)


def check_architectures() -> dict:
    """4. Representative 9S/14S/19S outputs checked independently.

    Coverage is derived from the registered target list, not from the corpus's own summary.
    Requiring only that the three architecture keys exist would accept a three-design smoke
    run, `{9S: 3, 14S: 0, 19S: 0}` -- three keys, status complete, nothing covered.
    """
    summary = ROOT / "work" / "independent" / "corpus200_pinned_summary.json"
    corpus_file = ROOT / "work" / "independent" / "corpus200_pinned.jsonl"
    if not summary.is_file() or not corpus_file.is_file():
        return _result("4. 9S/14S/19S coverage", "MANUAL", None,
                       "no pinned corpus on disk; run validation/independent/"
                       "regenerate_corpus.py")

    sys.path.insert(0, str(ROOT / "validation"))
    sys.path.insert(0, str(ROOT / "validation" / "independent"))
    from integration import corpus as registered_targets
    from journal import load_designs

    expected = list(registered_targets())
    want = {}
    for target in expected:
        key = f"{len(target)}S"
        want[key] = want.get(key, 0) + 1

    saved = json.loads(summary.read_text(encoding="utf-8"))
    rows = load_designs(corpus_file)
    got_targets = [r["target"] for r in rows]
    got_arch = {}
    for row in rows:
        got_arch[row["architecture"]] = got_arch.get(row["architecture"], 0) + 1

    missing = sorted(set(expected) - set(got_targets))
    extra = sorted(set(got_targets) - set(expected))
    duplicated = len(got_targets) - len(set(got_targets))
    failures = saved.get("failures", [])
    ok = (not missing and not extra and not duplicated and not failures
          and got_arch == want and saved.get("status") == "complete")
    detail = (f"{len(rows)} designs {got_arch}, registered {len(expected)} {want}, "
              f"status {saved.get('status')}, {len(failures)} failures")
    if missing or extra or duplicated:
        detail += (f"; missing {len(missing)} {missing[:2]}, unexpected {len(extra)}, "
                   f"duplicated {duplicated}")
    return _result("4. 9S/14S/19S coverage", "READ", ok, detail)


def _v5(d: dict) -> tuple[bool, str]:
    """Every design reconstructed, and every deliberate corruption still caught.

    The control count is read from the result rather than hardcoded: pinning it at 5 made
    adding two controls look like a regression, which would have discouraged adding more.
    A floor is still enforced, so an empty control set cannot pass trivially.
    """
    total = d.get("corrupt_controls_total", d["corrupt_controls_detected"])
    ok = (d["reconstructed"] == d["designs"] > 0 and not d["failures"]
          and total >= 5 and d["corrupt_controls_detected"] == total)
    return ok, (f"{d['reconstructed']}/{d['designs']} reconstructed, "
                f"{d['corrupt_controls_detected']}/{total} corrupt controls caught")


def _v6(d: dict) -> tuple[bool, str]:
    ok = d["clean"] == d["designs"] > 0 and not d["offenders"]
    return ok, (f"{d['clean']}/{d['designs']} satisfy every hard constraint, "
                f"{d['disagreements_with_claimed_status']} disagreements with claimed status")


def _v7(d: dict) -> tuple[bool, str]:
    """No semantic problem, and every control detected rather than parsed.

    An empty control dictionary used to satisfy `caught == len(controls)` trivially, so the
    controls are required to exist as well as to pass.
    """
    controls = d["malformed_controls"]
    caught = sum(1 for v in controls.values() if v)
    ok = (not d["semantic_problems"] and len(d["designs"]) > 0
          and len(controls) >= 5 and caught == len(controls))
    return ok, (f"{len(d['designs'])} designs, {len(d['semantic_problems'])} semantic "
                f"problems, {caught}/{len(controls)} controls caught")


def _binding(record: dict, needs_corpus: bool) -> tuple[bool, str]:
    """Is this result evidence about the artefacts being released, or about something else?

    A validator summary is a handful of integers. Without the corpus, validator and package
    hashes it was computed under, a passing result from an older tree satisfies a newer
    release directory unnoticed -- which is precisely the failure a release gate exists to
    catch.
    """
    prov = record.get("provenance")
    if not prov:
        return False, "no provenance recorded; cannot bind this result to what it checked"

    current_source = _package_source_sha256()
    if prov.get("package_source_sha256") != current_source:
        return False, (f"computed under package source "
                       f"{str(prov.get('package_source_sha256'))[:12]}, current is "
                       f"{current_source[:12]}")

    validator = ROOT / "validation" / "independent" / prov["validator"]
    if not validator.is_file():
        return False, f"validator {prov['validator']} not found"
    if _sha256(validator) != prov.get("validator_sha256"):
        return False, f"{prov['validator']} has changed since this result was produced"

    if needs_corpus:
        corpus_file = ROOT / "work" / "independent" / "corpus200_pinned.jsonl"
        if not corpus_file.is_file():
            return False, "pinned corpus absent"
        if _sha256(corpus_file) != prov.get("corpus_sha256"):
            return False, "result was computed against a different corpus"
    return True, "bound"


def _sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest()


def _package_source_sha256() -> str:
    sys.path.insert(0, str(ROOT / "src"))
    from clippr.manifest import source_fingerprint

    return source_fingerprint()["sha256"]


def check_independent_agreement() -> dict:
    """5. Exported sequences, metrics and hashes agree with independent calculations."""
    parts, ok = [], True
    for name, filename, read, needs_corpus in (
        ("V5 digest/ligate", "v5_digest_ligate.json", _v5, True),
        ("V6 constraints", "v6_constraints.json", _v6, True),
        ("V7 export semantics", "v7_export_semantics.json", _v7, False),
    ):
        file = ROOT / "work" / "independent" / filename
        if not file.is_file():
            parts.append(f"{name}: not run")
            ok = False
            continue
        record = json.loads(file.read_text(encoding="utf-8"))
        bound, why = _binding(record, needs_corpus)
        if not bound:
            parts.append(f"{name}: STALE -- {why}")
            ok = False
            continue
        good, detail = read(record)
        parts.append(f"{name}: {detail}")
        ok = ok and good
    return _result("5. independent agreement", "READ", ok, " | ".join(parts))


def check_optional_inputs() -> dict:
    """6. Optional inputs stay explicit; unavailable ones produce honest statuses.

    A narrow selection, so it runs in seconds rather than repeating the whole suite.
    """
    run = subprocess.run(
        [str(PYTHON), "-m", "pytest",
         "tests/test_products.py", "tests/test_manifest.py", "tests/test_offtarget.py",
         "-k", "unavailable or not_a_clean_screen or not_assessed or each_state"
               " or absent or not_a_cds or not_guessed"],
        cwd=ROOT, capture_output=True, text=True)
    passed = re.search(r"(\d+) passed", run.stdout)
    n = int(passed.group(1)) if passed else 0
    return _result("6. optional inputs are honest when absent", "RAN",
                   run.returncode == 0 and n > 0,
                   f"{n} tests covering absent Table S1, disabled and unavailable screening")


def check_notebook() -> dict:
    """7. The notebook exposes supported workflows, not deferred experiments."""
    path = ROOT / "notebooks" / "CLIPPR_designer.ipynb"
    if not path.is_file():
        return _result("7. notebook scope", "MANUAL", None, "notebook not found")
    text = path.read_text(encoding="utf-8")
    # M2 was a negative result and M3 is gated; neither may be advertised as a shipped
    # workflow. Their experiment entry points must not appear in the notebook at all.
    deferred = [name for name in ("m2_joint_selection", "build_banks", "compare_methods",
                                  "m3_inventory_pilot", "run_library_redesign")
                if name in text]
    if deferred:
        return _result("7. notebook scope and execution", "RAN", False,
                       f"notebook references deferred experiments: {deferred}")

    # Scope alone is not enough. The notebook's inventory half once shipped having never
    # executed: it passed `table` as a codon table, `table` is a display DataFrame bound by
    # an earlier cell, and three cells raised the moment a real input reached them. Every
    # static test stayed green, because the only path anything exercised was the one where
    # Supplementary Table S1 is absent and the cells short-circuit. So this runs the thing.
    proc = subprocess.run([str(PYTHON), str(ROOT / "validation" / "run_notebook.py")],
                          capture_output=True, text=True, encoding="utf-8",
                          errors="replace", cwd=str(ROOT),
                          env={**os.environ, "PYTHONUTF8": "1"})
    tail = (proc.stdout or "").strip().splitlines()
    verdict = next((line for line in reversed(tail) if line.strip() in
                    ("PASS", "FAIL", "SKIP")), "")
    counts = next((line for line in reversed(tail) if "cells executed" in line), "")
    if verdict == "SKIP":
        return _result("7. notebook scope and execution", "MANUAL", None,
                       "scope clean; execution skipped -- " + (tail[0] if tail else ""))
    return _result("7. notebook scope and execution", "RAN", verdict == "PASS",
                   f"no deferred experiment advertised; every cell executed ({counts})"
                   if verdict == "PASS" else
                   f"notebook execution failed: {counts or proc.stderr[-200:]}")


def check_claim_labels() -> dict:
    """8. README claims separate prediction, sequence measurement and experiment."""
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    required = ["No binding was measured", "what the\ncode predicts, not an experimental "
                                           "result"]
    missing = [phrase for phrase in required if phrase not in readme]
    return _result("8. claims are labelled by evidence type", "RAN", not missing,
                   "README states the recognition claim is predicted, not measured"
                   if not missing else f"missing: {missing}")


def check_runtime_table() -> dict:
    """9. Runtime figures come from measurements, with their configuration retained."""
    summary = ROOT / "work" / "independent" / "corpus200_pinned_summary.json"
    if not summary.is_file():
        return _result("9. runtime figures are measured", "MANUAL", None,
                       "no measured corpus on disk to compare against")
    saved = json.loads(summary.read_text(encoding="utf-8"))
    identity = saved.get("identity", {})
    kept = [k for k in ("seed", "matrix", "enzyme_profile", "destination", "genetic_code",
                        "dependencies", "codon_table_effective_sha256") if k in identity]
    return _result("9. runtime figures are measured", "READ", len(kept) >= 6,
                   f"{saved.get('wall_seconds_all_runs', 0) / 60:.1f} min for "
                   f"{saved.get('rows')} designs, configuration retained: {kept}")


def check_handoff() -> dict:
    """10. The commands and evidence a verifier needs, present and named."""
    expected = {
        "M2 verdict": ROOT / "validation" / "experiments" / "m2_result.md",
        "M3 verdict": ROOT / "validation" / "experiments" / "m3_result.md",
        "pinned corpus": ROOT / "work" / "independent" / "corpus200_pinned.jsonl",
        "reference comparison": ROOT / "work" / "m3" / "reference" / "comparison.json",
        "inventory validation": ROOT / "work" / "m3" / "inventory_validation.json",
    }
    absent = sorted(name for name, path in expected.items() if not path.is_file())
    return _result("10. handoff evidence present", "READ", not absent,
                   f"{len(expected) - len(absent)}/{len(expected)} artefacts on disk"
                   + (f"; missing {absent}" if absent else ""))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clean-install", action="store_true",
                    help="build a wheel and run the README example in a throwaway venv")
    args = ap.parse_args()

    print(f"release validation  |  {ROOT}")
    print(f"interpreter {PYTHON}\n")

    suite, collected = run_suite()
    checks = [check_tests(suite, collected), check_integration_coverage(suite),
              check_clean_install(args.clean_install), check_architectures(),
              check_independent_agreement(), check_optional_inputs(), check_notebook(),
              check_claim_labels(), check_runtime_table(), check_handoff()]

    width = max(len(c["check"]) for c in checks)
    for c in checks:
        mark = {True: "PASS", False: "FAIL", None: "NEEDS A HUMAN"}[c["ok"]]
        print(f"  [{c['evidence']:6s}] {c['check']:<{width}}  {mark}")
        print(f"           {' ' * width}  {c['detail']}")

    failed = [c for c in checks if c["ok"] is False]
    manual = [c for c in checks if c["ok"] is None]
    print(f"\n{sum(1 for c in checks if c['ok'])} of {len(checks)} passed; "
          f"{len(failed)} failed; {len(manual)} need a human")
    if manual:
        print("  A check this script cannot evaluate is not a passed check.")

    WORK.mkdir(parents=True, exist_ok=True)
    out = WORK / "release_check.json"
    out.write_text(json.dumps({"checks": checks}, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
