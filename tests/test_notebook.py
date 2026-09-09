"""Static checks on the Colab notebook.

**Why these exist.** The notebook is the only interface most users will ever touch, and until
now nothing in the suite covered it: a renamed function, a dropdown option that no longer
resolves, or a hand-edit destined to be overwritten would all leave 439 tests passing and the
notebook broken. These are the failure modes that are silent everywhere else.

They are deliberately all *static* -- parse, compare, cross-reference -- so they cost
milliseconds and run on every commit. Actually executing the cells is a different and much
slower job, and lives in `validation/run_notebook.py`.
"""
from __future__ import annotations

import ast
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "CLIPPR_designer.ipynb"
GENERATOR = ROOT / "tools" / "build_notebook.py"


@pytest.fixture(scope="module")
def nb():
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def code_cells(nb):
    return ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]


@pytest.fixture(scope="module")
def all_code(code_cells):
    return "\n".join(code_cells)


def as_python(source: str) -> str:
    """Notebook source with IPython magics blanked, so `ast` can parse it.

    `%pip install ...` and `!command` are valid notebook lines and invalid Python. Blanking
    rather than deleting them keeps line numbers aligned, so a reported syntax error still
    points at the right line of the cell.
    """
    return "\n".join("" if line.lstrip().startswith(("%", "!")) else line
                     for line in source.splitlines())


def params(source: str) -> list[tuple[str, str, str]]:
    """(variable, literal default, the #@param spec) for every Colab form field."""
    out = []
    for line in source.splitlines():
        if "#@param" not in line:
            continue
        assign, spec = line.split("#@param", 1)
        name, _, default = assign.partition("=")
        out.append((name.strip(), default.strip(), spec.strip()))
    return out


class TestStructure:
    def test_is_valid_notebook_json(self, nb):
        assert nb["nbformat"] == 4
        assert nb["cells"]

    def test_every_code_cell_parses(self, code_cells):
        for i, src in enumerate(code_cells, 1):
            try:
                ast.parse(as_python(src))
            except SyntaxError as exc:
                pytest.fail(f"code cell {i} does not parse: {exc}")

    def test_magics_are_only_used_for_installation(self, code_cells):
        """A stray magic elsewhere would run in Colab and fail everywhere else."""
        for i, src in enumerate(code_cells, 1):
            for line in src.splitlines():
                if line.lstrip().startswith(("%", "!")):
                    assert "pip" in line, f"cell {i} uses a non-pip magic: {line.strip()}"

    def test_no_outputs_are_committed(self, nb):
        """Stored outputs bloat the diff and show users stale results from another machine."""
        for i, c in enumerate(nb["cells"], 1):
            if c["cell_type"] == "code":
                assert not c.get("outputs"), f"cell {i} has stored output"
                assert not c.get("execution_count"), f"cell {i} has an execution count"

    def test_cell_titles_are_unique(self, code_cells):
        """Colab shows the title when a form cell is collapsed; duplicates are unnavigable."""
        titles = [m.group(1).strip() for src in code_cells
                  if (m := re.search(r"#@title (.+)", src))]
        assert len(titles) == len(set(titles)), "duplicate cell titles"

    def test_the_referenced_diagram_exists(self, nb):
        """The header markdown points at a file that has to ship with the notebook."""
        md = "\n".join("".join(c["source"]) for c in nb["cells"]
                       if c["cell_type"] == "markdown")
        for name in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", md):
            if name.startswith("http"):
                continue
            assert (NOTEBOOK.parent / name).exists(), f"missing asset {name}"


class TestNamesResolve:
    """Every clippr name the notebook imports must still exist. This is the silent one."""

    def test_imported_names_exist(self, all_code):
        import clippr

        missing = []
        for module, names in re.findall(r"from (clippr[.\w]*) import ([^\n(]+)", all_code):
            mod = clippr
            if module != "clippr":
                mod = __import__(module, fromlist=["_"])
            for name in (n.strip() for n in names.split(",")):
                if name and not hasattr(mod, name):
                    missing.append(f"{module}.{name}")
        assert not missing, f"the notebook imports names that no longer exist: {missing}"

    def test_design_oneshot_accepts_every_keyword_the_notebook_passes(self, all_code):
        """A renamed or removed parameter breaks 'Run all' while the suite stays green."""
        import inspect

        from clippr import design_oneshot

        accepted = set(inspect.signature(design_oneshot).parameters)
        call = re.search(r"design_oneshot\(\s*(.*?)\)\s*$", all_code,
                         re.DOTALL | re.MULTILINE)
        assert call, "expected a design_oneshot call in the notebook"
        passed = set(re.findall(r"(\w+)\s*=", call.group(1)))
        assert passed <= accepted, f"unknown keywords: {sorted(passed - accepted)}"


class TestFormDefaults:
    """Colab's defaults are what a stranger runs unedited, so they must all be valid."""

    def test_every_param_has_a_spec(self, code_cells):
        for src in code_cells:
            for name, default, spec in params(src):
                assert spec.startswith("{") or spec.startswith("["), \
                    f"{name} has a malformed #@param spec: {spec}"
                assert default, f"{name} has no default"

    def test_default_target_is_a_valid_architecture(self, all_code):
        from clippr import describe

        default = next(d for n, d, _ in params(all_code) if n == "target_rna")
        assert describe(default.strip('"'))["architecture"] in ("9S", "14S", "19S")

    def test_organism_options_all_resolve(self, all_code):
        from clippr import ORGANISMS

        spec = next(s for n, _, s in params(all_code) if n == "organism")
        for option in re.findall(r'"([^"]+)"', spec):
            assert option in ORGANISMS, f"organism dropdown offers unknown {option!r}"

    def test_enzyme_profile_options_all_resolve(self, all_code):
        from clippr import ENZYME_PROFILES

        spec = next(s for n, _, s in params(all_code) if n == "enzyme_profile")
        for option in re.findall(r'"([^"]+)"', spec):
            assert option in ENZYME_PROFILES, f"unknown enzyme profile {option!r}"

    def test_destination_options_all_resolve(self, all_code):
        from clippr import DESTINATION_OVERHANGS

        spec = next(s for n, _, s in params(all_code) if n == "destination_level")
        for option in re.findall(r'"([^"]+)"', spec):
            assert option in DESTINATION_OVERHANGS, f"unknown destination {option!r}"

    def test_ligation_table_options_exist_as_files(self, all_code):
        spec = next(s for n, _, s in params(all_code) if n == "ligation_table")
        matrices = ROOT / "src" / "clippr" / "data" / "matrices"
        for option in re.findall(r'"([^"]+)"', spec):
            assert (matrices / f"{option}.xlsx").exists(), f"no matrix file for {option!r}"

    def test_assembly_enzyme_options_are_known_to_biopython(self, all_code):
        from Bio import Restriction

        spec = next(s for n, _, s in params(all_code) if n == "assembly_enzyme")
        for option in re.findall(r'"([^"]+)"', spec):
            assert getattr(Restriction, option, None) is not None, f"unknown enzyme {option!r}"

    def test_every_default_is_among_its_own_options(self, code_cells):
        """A default outside its dropdown makes the form open in an invalid state."""
        for src in code_cells:
            for name, default, spec in params(src):
                if not spec.startswith("["):
                    continue
                options = re.findall(r'"([^"]+)"', spec)
                assert default.strip('"') in options, \
                    f"{name} defaults to {default} which is not among {options}"


class TestGeneratedNotebookIsCurrent:
    """The notebook is generated, so a hand-edit is silently lost on the next build."""

    def test_regenerating_would_not_change_it(self, tmp_path):
        import subprocess
        import sys

        before = NOTEBOOK.read_bytes()
        result = subprocess.run([sys.executable, str(GENERATOR)], cwd=ROOT,
                                capture_output=True, text=True,
                                env={"PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8",
                                     "PATH": ""} | dict(__import__("os").environ))
        assert result.returncode == 0, f"the generator failed: {result.stderr[-400:]}"
        after = NOTEBOOK.read_bytes()
        if before != after:
            NOTEBOOK.write_bytes(before)      # leave the tree as we found it
            pytest.fail("notebooks/CLIPPR_designer.ipynb differs from what "
                        "tools/build_notebook.py generates — edit the generator, not the "
                        "notebook, then rebuild")


class TestColabLinks:
    def test_badge_and_install_point_at_the_same_repository(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        badge = re.search(r"colab\.research\.google\.com/github/([\w-]+/[\w-]+)/", readme)
        install = re.search(r"pip install git\+https://github\.com/([\w-]+/[\w-]+)", readme)
        assert badge and install, "README should carry both a Colab badge and an install line"
        assert badge.group(1) == install.group(1), (
            f"the Colab badge points at {badge.group(1)} but install uses "
            f"{install.group(1)}")

    def test_badge_points_at_the_notebook_that_exists(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        m = re.search(r"colab\.research\.google\.com/github/[\w-]+/[\w-]+/blob/[\w-]+/(\S+?)\)",
                      readme)
        assert m, "no Colab notebook path found in the README"
        assert (ROOT / m.group(1)).exists(), f"badge points at missing {m.group(1)}"
