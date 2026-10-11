"""Controls for what the notebook harness is allowed to skip.

The harness once skipped any cell containing a permitted magic. Naming `%time`
unsupported fixed one symptom and left the mechanism: a required cell beginning
`%%capture` and then raising reported "1 skipped" and PASS. Skipping is now a property
of a line, not of a cell.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

HARNESS = (Path(__file__).resolve().parents[1] / "validation"
           / "notebook_on_published_build.py")


@pytest.fixture(scope="module")
def harness():
    spec = importlib.util.spec_from_file_location("_nb_harness", HARNESS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestCellClassification:
    def test_a_pure_install_cell_is_skippable(self, harness):
        keep, dropped, unsupported = harness.classify_cell("%pip install ppropt\n")
        assert dropped and not unsupported
        assert not harness._has_statement("\n".join(keep))

    def test_capture_does_not_hide_the_body(self, harness):
        """The exact control: a required cell whose first line is `%%capture`."""
        src = "%%capture\nraise RuntimeError('this cell must not be skipped')\n"
        keep, dropped, unsupported = harness.classify_cell(src)

        assert not unsupported
        assert dropped == ["%%capture"]
        body = "\n".join(keep)
        assert harness._has_statement(body), "the body must still be work to run"
        with pytest.raises(RuntimeError):
            exec(compile(body, "<control>", "exec"), {})

    def test_load_ext_does_not_hide_the_body(self, harness):
        src = "%load_ext autoreload\nresult = 1 / 0\n"
        keep, dropped, unsupported = harness.classify_cell(src)
        assert dropped == ["%load_ext autoreload"] and not unsupported
        assert harness._has_statement("\n".join(keep))

    def test_an_install_line_beside_real_work_keeps_the_work(self, harness):
        keep, dropped, _ = harness.classify_cell("%pip install ppropt\nx = 2\n")
        assert dropped and harness._has_statement("\n".join(keep))

    @pytest.mark.parametrize("magic", ["%time", "%%time", "%matplotlib inline",
                                       "!rm -rf out", "%cd .."])
    def test_unrecognised_magic_is_refused_not_skipped(self, harness, magic):
        _, _, unsupported = harness.classify_cell(f"{magic}\nx = 1\n")
        assert unsupported == [magic]

    def test_capture_is_not_a_skippable_line(self, harness):
        """Regression on the constant itself, which is how the defect was introduced."""
        assert "%%capture" not in harness.SKIPPABLE_LINES
        assert "%load_ext" in harness.SKIPPABLE_LINES

    def test_comments_alone_are_not_work(self, harness):
        assert not harness._has_statement("# just a note\n\n")
        assert harness._has_statement("x = 1\n")

    def test_broken_source_is_not_excused(self, harness):
        """A cell we cannot parse is the executor's to report, not ours to skip."""
        assert harness._has_statement("def (:\n")
