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
import hashlib
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

    area = Path(tempfile.mkdtemp(prefix="ppropt-release-"))
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
        env.update({"PPROPT_CACHE_DIR": str(area / "cache"), "PYTHONIOENCODING": "utf-8",
                    "PYTHONNOUSERSITE": "1"})

        origin = subprocess.run(
            [str(python), "-c", "import ppropt,sys;print(ppropt.__file__);"
                                "print(ppropt.__version__ if hasattr(ppropt,'__version__') "
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
            detail = (f"ppropt imported from {imported[0] if imported else '?'}, which is "
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
    # Scoped. V6 checks the coding sequence; the final wrapped oligo is a separate check
    # that rejects this same corpus. "every hard constraint" read as both.
    return ok, (f"{d['clean']}/{d['designs']} satisfy the V6 coding-sequence constraints, "
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


def _inputs_fingerprint() -> str:
    sys.path.insert(0, str(ROOT / "validation"))
    from state_table import inputs_fingerprint

    return inputs_fingerprint()


def _git_head() -> str:
    """HEAD, or a stated absence. Never raises: this runs inside the record write, after the
    whole suite, so an exception here would discard the entire run's evidence."""
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(ROOT),
                             capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError) as exc:
        return f"unavailable: {type(exc).__name__}"
    if out.returncode != 0:
        return "unavailable: not a git checkout"
    return out.stdout.strip() or "unavailable: empty rev-parse"


def _package_source_sha256() -> str:
    sys.path.insert(0, str(ROOT / "src"))
    from ppropt.manifest import source_fingerprint

    return source_fingerprint()["sha256"]


def check_independent_agreement(evidence_dir: str | None = None) -> dict:
    """5. Exported sequences, metrics and hashes agree with independent calculations.

    Reads the release bundle when one is given. The shared `work/independent` path stays
    as the default convenience, but a release must be judged on its own records, not on
    whatever a later ad-hoc run happened to leave in a shared directory.
    """
    parts, ok = [], True
    for name, filename, read, needs_corpus in (
        ("V5 digest/ligate", "v5_digest_ligate.json", _v5, True),
        ("V6 constraints", "v6_constraints.json", _v6, True),
        ("V7 export semantics", "v7_export_semantics.json", _v7, False),
    ):
        file = ROOT / "work" / "independent" / filename
        if evidence_dir and (Path(evidence_dir) / filename).is_file():
            file = Path(evidence_dir) / filename
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


#: Both shipped notebooks. The inventory route used to be the second half of the designer
#: notebook, and this gate ran it by running that one file. It is its own notebook now, so
#: naming only the designer here would leave half the product never executed -- which is the
#: exact failure the comment below records.
NOTEBOOKS = ("PPROPT_designer.ipynb", "PPROPT_inventory.ipynb")


def check_notebook() -> dict:
    """7. The notebooks expose supported workflows, not deferred experiments."""
    paths = [ROOT / "notebooks" / name for name in NOTEBOOKS]
    absent = [p.name for p in paths if not p.is_file()]
    if absent:
        return _result("7. notebook scope", "MANUAL", None, f"not found: {absent}")
    text = "".join(p.read_text(encoding="utf-8") for p in paths)
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
    executed, failed, skipped = [], [], []
    for path in paths:
        proc = subprocess.run([str(PYTHON), str(ROOT / "validation" / "run_notebook.py"),
                               "--notebook", str(path)],
                              capture_output=True, text=True, encoding="utf-8",
                              errors="replace", cwd=str(ROOT),
                              env={**os.environ, "PYTHONUTF8": "1"})
        tail = (proc.stdout or "").strip().splitlines()
        verdict = next((line for line in reversed(tail) if line.strip() in
                        ("PASS", "FAIL", "SKIP")), "")
        counts = next((line for line in reversed(tail) if "cells executed" in line), "")
        if verdict == "SKIP":
            skipped.append(f"{path.name}: {tail[0] if tail else 'skipped'}")
        elif verdict == "PASS":
            executed.append(f"{path.name} ({counts})")
        else:
            failed.append(f"{path.name}: {counts or proc.stderr[-200:]}")

    if failed:
        return _result("7. notebook scope and execution", "RAN", False,
                       "notebook execution failed -- " + "; ".join(failed))
    if skipped and not executed:
        return _result("7. notebook scope and execution", "MANUAL", None,
                       "scope clean; execution skipped -- " + "; ".join(skipped))
    detail = "no deferred experiment advertised; every cell executed in " +              ", ".join(executed)
    if skipped:
        detail += "; skipped " + "; ".join(skipped)
    return _result("7. notebook scope and execution", "RAN", True, detail)


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


def check_badge_count(run: subprocess.CompletedProcess, collected: int) -> dict:
    """12. The notebook's tests badge equals the count the suite actually reports.

    The badge is a committed constant so notebook generation stays reproducible from the
    repository alone. That makes it capable of going stale, which is what this check is for:
    it read 283 against a real 949, then 927 against a measured 974, both times for long
    enough to reach users.
    """
    import re as _re

    source = (ROOT / "tools" / "build_notebook.py").read_text(encoding="utf-8")
    declared = _re.search(r"^TESTS = (\d+)$", source, _re.M)
    if not declared:
        return _result("12. notebook tests badge is current", "RAN", False,
                       "tools/build_notebook.py has no `TESTS = <int>` line to check")
    passed = _re.search(r"(\d+) passed", run.stdout)
    if not passed or run.returncode != 0:
        return _result("12. notebook tests badge is current", "MANUAL", None,
                       "the suite reported no usable count, so the badge cannot be checked")
    shown, measured = int(declared.group(1)), int(passed.group(1))
    return _result("12. notebook tests badge is current", "RAN", shown == measured,
                   f"badge says {shown}, suite reports {measured}"
                   + ("" if shown == measured
                      else f" -- set TESTS = {measured} in tools/build_notebook.py and rebuild"))


#: Every record a candidate release must carry, and where in the bundle it lives. B-E and
#: V5-V7 used to be absent from this list, so a directory holding four files passed while
#: the rest of the evidence could be anything, anywhere, from any run.
REQUIRED_RECORDS = {
    "gate_b": ("gate_b.json",),
    "gate_c": ("phasec", "gate_c.json"),
    "gate_d": ("phased", "gate_d.json"),
    "gate_e": ("phasee", "gate_e.json"),
    "gate_f": ("phasef", "gate_f.json"),
    "notebook_designer": ("notebook_designer.json",),
    "notebook_inventory": ("notebook_inventory.json",),
    "final_artifact_contract": ("final_artifact_contract.json",),
    "v5_digest_ligate": ("v5_digest_ligate.json",),
    "v6_constraints": ("v6_constraints.json",),
    "v7_export_semantics": ("v7_export_semantics.json",),
}

#: What each notebook must actually have demonstrated. A positive postcondition count is
#: not a contract: `--expect "1"` satisfied the old gate. Each entry is a substring that
#: must appear in some declared expectation for that notebook.
NOTEBOOK_CONTRACT = {
    "notebook_designer": ("final_artifacts", "synthesis_ready", "oligos"),
    "notebook_inventory": ("deposited", "order_items", "order.ok", "summary["),
}

#: The inputs each notebook must have been given, staged and hashed. The staged-input loop
#: can only re-hash what a record declares, so a record declaring none passed every check
#: while the run had in fact read its codon table and genome from wherever the machine
#: happened to keep them. Naming them here makes "declared nothing" a failure.
NOTEBOOK_REQUIRED_INPUTS = {
    "notebook_designer": ("codon_tables/kazusa_3055.json", "genomes/NC_005353.1.gb"),
    "notebook_inventory": ("Table S1.xlsx", "codon_tables/kazusa_3055.json"),
}

#: The panel the corpus contract must cover.
EXPECTED_DESIGNS = 200


def _under_root(path) -> Path:
    """Anchor a recorded path to the repository, not to whatever directory we ran from.

    Records store paths as they were typed, which is usually repo-relative. Resolving
    those against the process CWD made the verdict depend on where the check was invoked:
    the same bundle verified from the repo root and failed from anywhere else.
    """
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def _file_sha256(path) -> str | None:
    try:
        return hashlib.sha256(_under_root(path).read_bytes()).hexdigest()
    except OSError:
        return None


def _check_notebook_record(name, record, current, problems):
    """A notebook record, verified against the files it names."""
    if not record.get("ok"):
        problems.append(f"{name} records ok={record.get('ok')!r}")

    # Recompute, do not read. A zeroed hash used to pass because the field was non-empty.
    stated = record.get("notebook_sha256")
    notebook = record.get("notebook")
    if not stated:
        problems.append(f"{name} carries no notebook_sha256")
    elif not notebook or not _under_root(notebook).is_file():
        problems.append(f"{name} names no readable notebook to verify its hash against")
    else:
        actual = _file_sha256(notebook)
        if actual != stated:
            problems.append(f"{name} records notebook {str(stated)[:12]} but "
                            f"{Path(notebook).name} hashes to {str(actual)[:12]}")

    for field in ("installed_source_fingerprint", "checkout_source_fingerprint",
                  "interpreter", "environment_prefix"):
        if not record.get(field):
            problems.append(f"{name} carries no {field}")
    installed = record.get("installed_source_fingerprint")
    checkout = record.get("checkout_source_fingerprint")
    if installed and checkout and installed != checkout:
        problems.append(f"{name} ran against a wheel that is not this source")
    if installed and current and installed != current:
        problems.append(f"{name} ran against {str(installed)[:12]}, but this checkout "
                        f"is {current[:12]}")

    # The loop below re-hashes what the record declares, so it can only ever be as strong
    # as the declaration. Require the inputs first; otherwise declaring none passes.
    staged_inputs = record.get("staged_inputs") or {}
    for required in NOTEBOOK_REQUIRED_INPUTS.get(name, ()):
        if not any(required in str(dest) for dest in staged_inputs):
            problems.append(f"{name} staged no {required}, so the run cannot be shown "
                            f"to have used the inputs this release names")

    # Every staged input, re-hashed from the file it came from.
    for dest, staged in staged_inputs.items():
        source = staged.get("from")
        stated = staged.get("sha256")
        if not stated:
            problems.append(f"{name} staged {dest} with no hash")
        elif not source or not _under_root(source).is_file():
            problems.append(f"{name} staged {dest} from a source that is not readable")
        elif _file_sha256(source) != stated:
            problems.append(f"{name} staged {dest} as {str(stated)[:12]}, but its "
                            f"source now hashes to {str(_file_sha256(source))[:12]}")

    # Cells that ran are not work done.
    declared = record.get("postconditions_declared")
    entries = record.get("expectations") or []
    if not declared:
        problems.append(f"{name} declared no postcondition, so it cannot show the "
                        f"notebook did its work")
    elif declared != len(entries):
        problems.append(f"{name} claims {declared} postcondition(s) but lists "
                        f"{len(entries)}")
    unmet = [e for e in entries if not e.get("ok")]
    if unmet:
        problems.append(f"{name} postcondition not met: "
                        f"{unmet[0].get('expect', '?')[:60]}")

    # And the postconditions must be the declared ones, not any true statement.
    text = " ".join(str(e.get("expect", "")) for e in entries)
    for token in NOTEBOOK_CONTRACT.get(name, ()):
        if token not in text:
            problems.append(f"{name} declares no postcondition covering {token!r}, so "
                            f"its result contract is not established")


def _check_gate_f_record(record, current, problems):
    """Gate F names the wheel it built. Verify the wheel, not the claim."""
    if record.get("gate") != "PASS":
        problems.append(f"gate_f records {record.get('gate')!r}")
    installed = record.get("installed_source_fingerprint")
    if not installed:
        problems.append("gate_f carries no installed_source_fingerprint")
    elif current and installed != current:
        problems.append(f"gate_f built {str(installed)[:12]}, but this checkout is "
                        f"{current[:12]}")
    stated = (record.get("wheel") or {}).get("sha256") or record.get("wheel_sha256")
    path = (record.get("wheel") or {}).get("path") or record.get("wheel_path")
    if not stated:
        problems.append("gate_f carries no wheel hash")
    elif not path or not Path(path).is_file():
        problems.append("gate_f names no retained wheel to verify its hash against")
    elif _file_sha256(path) != stated:
        problems.append(f"gate_f records wheel {str(stated)[:12]} but the file hashes to "
                        f"{str(_file_sha256(path))[:12]}")


def _check_corpus_record(record, current, problems):
    """The corpus report, against the corpus file and the checker that produced it."""
    rows = record.get("rows") or []
    checked = len(rows)
    compared = record.get("comparisons_performed", 0)
    if compared <= 0:
        problems.append("the corpus checker performed no comparison")
    if record.get("disagreements_with_production"):
        problems.append("the corpus checker disagrees with production")
    if compared != checked:
        problems.append(f"the corpus checker compared {compared} of {checked} rows")
    if checked != EXPECTED_DESIGNS:
        problems.append(f"the corpus checker covered {checked} designs, not the declared "
                        f"{EXPECTED_DESIGNS}")
    if record.get("not_assessed"):
        problems.append(f"{record['not_assessed']} design(s) could not be assessed; "
                        f"required data is missing")

    prov = record.get("provenance") or {}
    # The package fingerprint does not cover validation/*.py, so the checker's own bytes
    # are part of what this result is evidence about.
    validator = ROOT / "validation" / "independent" / "final_artifact_contract.py"
    stated = prov.get("validator_sha256")
    if not stated:
        problems.append("the corpus report carries no validator_sha256")
    elif _file_sha256(validator) != stated:
        problems.append("the corpus report was produced by a different version of "
                        "final_artifact_contract.py than the one on disk")
    # "Compare only if both exist" lets a record with no stated package pass while being
    # bound to no source at all. A missing identity is a failure, not a skip.
    if not prov.get("package_source_sha256"):
        problems.append("the corpus report carries no package_source_sha256, so it is "
                        "not bound to any source tree")
    elif current and prov["package_source_sha256"] != current:
        problems.append(f"the corpus report ran against "
                        f"{prov['package_source_sha256'][:12]}, not {current[:12]}")

    # The panel must be the declared source rows, repeats included. A report whose rows
    # were one result copied two hundred times passed on its row count alone.
    corpus = prov.get("corpus_path")
    stated_corpus = prov.get("corpus_sha256")
    if not stated_corpus:
        problems.append("the corpus report carries no corpus_sha256")
    elif not corpus or not _under_root(corpus).is_file():
        problems.append("the corpus report names no readable corpus to verify")
    elif _file_sha256(corpus) != stated_corpus:
        problems.append(f"the corpus report records {str(stated_corpus)[:12]} but "
                        f"{Path(corpus).name} hashes to "
                        f"{str(_file_sha256(corpus))[:12]}")
    else:
        declared = []
        for line in _under_root(corpus).read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                entry = json.loads(line)
            except ValueError:
                continue
            if "__run__" in entry or "__run_identity__" in entry:
                continue
            # The same resolution the checker uses, so the two lists are comparable.
            declared.append(entry.get("target") or entry.get("target_rna"))
        reported = [r.get("target") for r in rows]
        if sorted(filter(None, reported)) != sorted(filter(None, declared)):
            problems.append(f"the corpus report's {len(reported)} rows are not the "
                            f"{len(declared)} declared panel rows")


def _check_validator_record(name, record, current, problems):
    """V5-V7: produced by the checker on disk, against this package."""
    prov = record.get("provenance") or {}
    stated = prov.get("validator_sha256")
    script = ROOT / "validation" / "independent" / f"{name}.py"
    if not stated:
        problems.append(f"{name} carries no validator_sha256")
    elif _file_sha256(script) != stated:
        problems.append(f"{name} was produced by a different version of {name}.py than "
                        f"the one on disk")
    if not prov.get("package_source_sha256"):
        problems.append(f"{name} carries no package_source_sha256, so it is not bound "
                        f"to any source tree")
    elif current and prov["package_source_sha256"] != current:
        problems.append(f"{name} ran against {prov['package_source_sha256'][:12]}, not "
                        f"{current[:12]}")


def check_evidence_dir(path: str | None) -> dict | None:
    """13. This run's evidence is present, complete and bound to this release.

    Returns None when no directory was given, so an ordinary run is unchanged.

    Every required record must be in the bundle, and every identity a record states is
    recomputed from the file it names. The earlier version compared two fingerprints "if
    both exist" and otherwise checked that fields were non-empty, so a record whose
    notebook hash, input hash, wheel hash or corpus hash had been replaced by zeros passed,
    as did a report whose two hundred rows were one result copied two hundred times.
    """
    if not path:
        return None
    where = Path(path)
    required = {name: where.joinpath(*parts)
                for name, parts in REQUIRED_RECORDS.items()}
    missing = sorted(name for name, f in required.items() if not f.is_file())
    if missing:
        return _result("13. this run's evidence is present", "RAN", False,
                       f"missing from {where}: {', '.join(missing)}")

    current = _package_source_sha256()
    problems: list[str] = []
    for name, f in required.items():
        try:
            record = json.loads(f.read_text(encoding="utf-8"))
        except ValueError as exc:
            problems.append(f"{name}: unreadable ({exc})")
            continue
        if name == "gate_f":
            _check_gate_f_record(record, current, problems)
        elif name.startswith("gate_"):
            # `not in ("PASS", None)` admitted a record with no verdict at all: a gate file
            # truncated to {} read as a pass. Every writer emits PASS or FAIL, so a missing
            # or null verdict has no legitimate source.
            if record.get("gate") is None:
                problems.append(f"{name} records no gate verdict")
            elif record["gate"] != "PASS":
                problems.append(f"{name} records {record['gate']!r}")
        elif name.startswith("notebook"):
            _check_notebook_record(name, record, current, problems)
        elif name == "final_artifact_contract":
            _check_corpus_record(record, current, problems)
        elif name.startswith("v"):
            _check_validator_record(name, record, current, problems)

    return _result("13. this run's evidence is present", "RAN", not problems,
                   "; ".join(problems) if problems
                   else f"{len(required)} records, every stated identity recomputed and "
                        f"matched, all from {where}")


def check_state_table() -> dict:
    """11. The handover's State table is rendered from measurements, not typed from memory."""
    script = ROOT / "validation" / "state_table.py"
    if not script.is_file():
        return _result("11. handover State table is measured", "MANUAL", None,
                       "validation/state_table.py is absent")
    proc = subprocess.run([str(PYTHON), str(script)], capture_output=True, text=True,
                          encoding="utf-8", errors="replace", cwd=str(ROOT))
    first = (proc.stdout or proc.stderr or "").strip().splitlines()
    return _result("11. handover State table is measured", "RAN", proc.returncode == 0,
                   first[0] if first else "state_table.py produced no output")


def main() -> int:
    from datetime import datetime, timezone

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clean-install", action="store_true",
                    help="build a wheel and run the README example in a throwaway venv")
    ap.add_argument("--evidence-dir", default=None,
                    help=("read this run's gate, notebook and artifact records from here "
                          "instead of the historical fixed paths; a missing required "
                          "record fails rather than falling back"))
    ap.add_argument("--refresh-state", action="store_true",
                    help=("after the measured record is written, re-render the handover's "
                          "State table and verify it; default behaviour does not mutate"))
    args = ap.parse_args()

    print(f"release validation  |  {ROOT}")
    print(f"interpreter {PYTHON}\n")

    suite, collected = run_suite()
    checks = [check_tests(suite, collected), check_integration_coverage(suite),
              check_clean_install(args.clean_install), check_architectures(),
              check_independent_agreement(args.evidence_dir), check_optional_inputs(), check_notebook(),
              check_claim_labels(), check_runtime_table(), check_handoff(),
              check_badge_count(suite, collected)]

    # Check 11 reads the record this run writes, so it runs after the write -- but it is
    # appended to `checks` and the record rewritten, so the tally, the record and the
    # generated row all describe eleven checks. Reporting "10 of 10" above a FAIL line was
    # the generator committing the sin it exists to prevent.
    width = max(len(c["check"]) for c in checks)
    for c in checks:
        mark = {True: "PASS", False: "FAIL", None: "NEEDS A HUMAN"}[c["ok"]]
        print(f"  [{c['evidence']:6s}] {c['check']:<{width}}  {mark}")
        print(f"           {' ' * width}  {c['detail']}")

    # Every check that does not read the record must be in the record before it is
    # written. Appending check 13 after the write meant the state table was rendered from
    # an interim 11-check file and the final file immediately invalidated it: the table
    # said "11 of 11" while the record said 12 of 12.
    evidence = check_evidence_dir(args.evidence_dir)
    if evidence is not None:
        checks.append(evidence)
        mark = {True: "PASS", False: "FAIL", None: "NEEDS A HUMAN"}[evidence["ok"]]
        print(f"  [{evidence['evidence']:6s}] {evidence['check']:<{width}}  {mark}")
        print(f"           {' ' * width}  {evidence['detail']}")

    WORK.mkdir(parents=True, exist_ok=True)
    out = WORK / "release_check.json"
    # Bind the record to the tree. Without these, anything reading it later cannot tell
    # whether the numbers still describe this commit.
    record = {"when": datetime.now(timezone.utc).isoformat(timespec="seconds"),
              "head": _git_head(),
              "package_source_sha256": _package_source_sha256(),
              "inputs_sha256": _inputs_fingerprint(),
              "checks": checks}
    out.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")

    # Only check 11 runs after the write, because it is the one check that reads the
    # record. Its own result is the single row the final rewrite adds, and the state
    # renderer excludes it, so the rendered tally stays true.
    if args.refresh_state:
        # Render from the record this run wrote, then verify. A renderer failure is a
        # release failure; it is never converted into success.
        script = ROOT / "validation" / "state_table.py"
        written = subprocess.run([str(PYTHON), str(script), "--write"],
                                 capture_output=True, text=True, cwd=str(ROOT))
        tail = (written.stdout or written.stderr).strip().splitlines()
        print(f"\n  state table --write: {(tail[-1] if tail else '(no output)')[:100]}")
        if written.returncode != 0:
            # A renderer that failed is a release failure. Its exit code was discarded,
            # so a failed rewrite was followed by a verification of whatever the file
            # already said.
            checks.append(_result("11. handover State table is measured", "RAN", False,
                                  f"state_table.py --write exited {written.returncode}: "
                                  f"{(written.stderr or '').strip()[:120]}"))
            record["checks"] = checks
            out.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
            print("\nthe State table could not be rendered; not reporting a tally")
            return 1

    state = check_state_table()
    mark = {True: "PASS", False: "FAIL", None: "NEEDS A HUMAN"}[state["ok"]]
    print(f"  [{state['evidence']:6s}] {state['check']:<{width}}  {mark}")
    print(f"           {' ' * width}  {state['detail']}")

    checks.append(state)
    failed = [c for c in checks if c["ok"] is False]
    manual = [c for c in checks if c["ok"] is None]
    print(f"\n{sum(1 for c in checks if c['ok'] is True)} of {len(checks)} passed; "
          f"{len(failed)} failed; {len(manual)} need a human")
    if manual:
        print("  A check this script cannot evaluate is not a passed check.")

    # Rewritten so the record names eleven checks. This does not re-stale the block: the
    # binding is the inputs digest, and `_release_row` is rendered from this final list.
    record["checks"] = checks
    out.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    # A copy inside the bundle, so every figure in a release document resolves to a file
    # later runs do not rewrite. `work/release/release_check.json` stays as the convenience
    # path; it is mutable, and a control run once overwrote it.
    if args.evidence_dir:
        frozen = Path(args.evidence_dir) / "release_check.json"
        frozen.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        digest = hashlib.sha256(frozen.read_bytes()).hexdigest()
        (Path(args.evidence_dir) / "release_check.sha256").write_text(
            f"{digest}  release_check.json\n", encoding="utf-8")
        print(f"\nfroze {frozen}\n  sha256 {digest}")
    return 1 if failed or manual else 0


if __name__ == "__main__":
    raise SystemExit(main())
