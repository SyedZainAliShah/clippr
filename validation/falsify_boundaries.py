"""Prove the boundary suite is not vacuous.

A regression that passes on the broken code tests nothing, and this project has shipped that
mistake before: six defects crossed a green test suite because every test exercised a path the
defect did not touch. So each repair is reverted here, one at a time, in a scratch worktree,
and the suite is run against the broken code. A regression that stays green under the defect it
claims to cover is reported as MISSED.

It has already earned this. On its first run it found three regressions that did not
discriminate what they claimed to: a permissive-export test that passed even when the profile
was ignored (both saved inventories satisfy the default, so only a *pair* of disagreeing calls
shows the argument arrived), a completeness test whose role map did not form a valid chain (so
the missing citation was never the deciding factor), and a revert this script had invented
rather than one the code ever had.

    python validation/falsify_boundaries.py

Exits non-zero if any reverted defect leaves the suite green.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PY = Path(sys.executable)
WORK = Path(tempfile.gettempdir()) / "clippr_falsify_tree"

#: Each entry reverts one repair by patching the *current* source back to its broken form, and
#: names the tests that must go red. `old` is what the repaired code says; `new` is the defect.
DEFECTS = [
    dict(
        name="F3 target collapses back to hard",
        file="src/clippr/synthesis_profile.py",
        old='''                "global_gc_bounds": tuple(self.global_gc) if self.global_gc else None,
                "soft_rules": self.soft_rules}''',
        new='''                "global_gc_bounds": tuple(self.global_gc) if self.global_gc else None,
                "soft_rules": ()}''',
        expect=["test_the_two_enforcements_send_different_arguments",
                "test_an_unsatisfiable_target_still_returns_a_sequence"],
    ),
    # The two halves of the F5 repair are not independent. Once __post_init__ rebuilds the map
    # into a fresh proxy, sharing it in narrowed() is harmless, because the child copies it on
    # the way in -- reverting that line alone leaves the suite green, and the standalone revert
    # was a defect this falsifier invented rather than one the code ever had. The state the
    # code was actually in is both at once: a plain dict, handed to the child by reference.
    dict(
        name="F5 as it was: a plain dict, shared with the child",
        file="src/clippr/synthesis_profile.py",
        old='''        object.__setattr__(self, "enforcement", MappingProxyType(dict(self.enforcement)))''',
        new='''        object.__setattr__(self, "enforcement", self.enforcement)''',
        expect=["test_the_enforcement_map_refuses_mutation",
                "test_a_derived_policy_does_not_share_its_parents_map"],
    ),
    dict(
        name="an enzyme site may be declared soft",
        file="src/clippr/synthesis_profile.py",
        old='''            if rule in NEVER_SOFT and how != HARD:''',
        new='''            if False:''',
        expect=["test_an_enzyme_site_can_never_be_declared_a_preference"],
    ),
    dict(
        name="F1 export ignores the profile it was handed",
        file="src/clippr/workflow.py",
        old='''                 for s in sub.substrates_for(library, profile)]''',
        new='''                 for s in sub.substrates_for(library)]''',
        expect=["test_two_policies_over_one_inventory_give_two_verdicts",
                "test_a_tighter_policy_fails_on_its_own_violations"],
    ),
    dict(
        name="F1 substrates ignore the profile they were built with",
        file="src/clippr/substrates.py",
        old='''def substrate_problems(insert: str, block: str, profile=None) -> list[str]:''',
        new='''def substrate_problems(insert: str, block: str, profile=None) -> list[str]:
    profile = None''',
        expect=["test_the_same_module_is_judged_differently_under_each",
                "test_a_tighter_policy_fails_on_its_own_violations"],
    ),
    dict(
        name="F4 the collection solver drops the bounds",
        file="src/clippr/library_search.py",
        old='''    bounds = _resolve(profile).solver_bounds()''',
        new='''    bounds = {}''',
        expect=["test_a_proposal_is_solved_under_the_profile_not_the_default"],
    ),
    dict(
        name="F2 the front records no policy",
        file="src/clippr/joint_search.py",
        old='''"synthesis_profile": _resolve_profile(profile).as_dict()''',
        new='''"synthesis_profile": None''',
        expect=["test_a_front_records_the_policy_it_was_measured_under",
                "test_every_point_on_the_front_selects_under_that_policy"],
    ),
    dict(
        name="F2 selection publishes after a failed verification",
        file="src/clippr/workflow.py",
        old='''    if not verified.feasible:
        # Never publish an artefact the verifier just rejected.''',
        new='''    if False:
        # Never publish an artefact the verifier just rejected.''',
        expect=["test_a_candidate_that_fails_verification_publishes_nothing"],
    ),
    dict(
        name="F6 completeness is prose plus one overhang",
        file="src/clippr/assembly_spec.py",
        old='''    complete = bool(evidence) and not missing and not malformed and chain_ok''',
        new='''    complete = bool(evidence) and bool(subset)''',
        expect=["test_prose_and_one_overhang_do_not_establish_completeness",
                "test_evidence_alone_is_not_enough"],
    ),
    dict(
        name="F6 completeness stops requiring a citation",
        file="src/clippr/assembly_spec.py",
        old='''    complete = bool(evidence) and not missing and not malformed and chain_ok''',
        new='''    complete = not missing and not malformed and chain_ok''',
        expect=["test_an_otherwise_complete_set_is_refused_for_the_citation_alone"],
    ),
]

TESTS = "tests/test_policy_boundaries.py"


def run(cmd, cwd=REPO, check=True):
    got = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True, shell=False)
    if check and got.returncode:
        raise SystemExit(f"{' '.join(map(str, cmd))} failed:\n{got.stdout}\n{got.stderr}")
    return got


def failing_tests(tree: Path) -> set[str]:
    """Names of the tests that fail in this tree."""
    got = subprocess.run([str(PY), "-m", "pytest", TESTS, "-q", "--no-header", "-p",
                          "no:cacheprovider", "--tb=no"],
                         cwd=str(tree), capture_output=True, text=True)
    out = got.stdout + got.stderr
    names = set(re.findall(r"::(\w+)::(\w+)", out))
    collected = {test for _cls, test in names}
    # An import-time failure fails everything; say so rather than reporting a clean sweep.
    if "error" in out.lower() and not collected:
        return {"<collection error>"}
    return collected


def main() -> int:
    if WORK.exists():
        run(["git", "worktree", "remove", "--force", str(WORK)], check=False)
    base = run(["git", "merge-base", "HEAD", "main"]).stdout.strip()
    print(f"worktree at merge-base {base[:10]}")
    run(["git", "worktree", "add", "--detach", str(WORK), base])

    # The tests and the repaired sources come from the working tree; only one file at a time
    # is reverted, so everything else is the current code.
    for path in ("src/clippr", "tests", "pyproject.toml", "data", "work"):
        source, target = REPO / path, WORK / path
        if source.is_dir():
            shutil.rmtree(target, ignore_errors=True)
            shutil.copytree(source, target)
        elif source.is_file():
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read_bytes())

    clean = failing_tests(WORK)
    if clean:
        print(f"!! the suite is not green in the worktree: {sorted(clean)}")
        return 2
    print("worktree reproduces a green suite; reverting repairs one at a time\n")

    verdicts = []
    for defect in DEFECTS:
        target = WORK / defect["file"]
        original = target.read_text(encoding="utf-8")
        if defect["old"] not in original:
            verdicts.append((defect["name"], "PATCH DID NOT APPLY", []))
            print(f"[skip] {defect['name']}: anchor not found in {defect['file']}")
            continue
        target.write_text(original.replace(defect["old"], defect["new"], 1), encoding="utf-8")
        red = failing_tests(WORK)
        target.write_text(original, encoding="utf-8")

        caught = sorted(set(defect["expect"]) & red)
        missed = sorted(set(defect["expect"]) - red)
        verdicts.append((defect["name"], "caught" if caught and not missed else
                         ("partial" if caught else "MISSED"), sorted(red)))
        mark = "ok " if caught and not missed else ("~~ " if caught else "!! ")
        print(f"{mark}{defect['name']}")
        print(f"     expected red: {', '.join(defect['expect'])}")
        print(f"     actually red: {', '.join(sorted(red)) or '(none -- the suite stayed green)'}")
        if missed:
            print(f"     NOT CAUGHT  : {', '.join(missed)}")
        print()

    run(["git", "worktree", "remove", "--force", str(WORK)], check=False)

    bad = [name for name, verdict, _ in verdicts if verdict not in ("caught",)]
    print("=" * 78)
    for name, verdict, _ in verdicts:
        print(f"  {verdict:<18} {name}")
    print("=" * 78)
    if bad:
        print(f"\n{len(bad)} defect(s) not fully caught by the regressions.")
        return 1
    print(f"\nAll {len(DEFECTS)} reverted defects are caught by the boundary suite.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
