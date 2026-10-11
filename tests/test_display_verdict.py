"""Controls for what a person actually reads: the card, the audit and the summary.

The API refused to export a rejected design while every user-facing surface still led
with a green `PASS` taken from the coding-span QC. These tests assert the rejection is
visible, and that the retained CDS result is labelled with its scope rather than deleted.

The card is built by `tools/build_notebook.py` and lives in the generated notebook, so it
is exercised here by running that cell's source against a real design.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from ppropt.design import design_oneshot

from test_design import TOY_TABLE

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "PPROPT_designer.ipynb"


@pytest.fixture(scope="module")
def rejected():
    """The shipped example, which the final-artifact check genuinely refuses."""
    out = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE, check_offtarget=False)
    assert out["synthesis_ready"] is False, "fixture assumes this example is refused"
    assert out["final_artifacts"]["status"] == "fail"
    return out


def _card_source():
    """The generated notebook's result cell, minus the design call and the display call."""
    cells = json.loads(NOTEBOOK.read_text(encoding="utf-8"))["cells"]
    for cell in cells:
        src = "".join(cell["source"])
        if cell["cell_type"] == "code" and 'qc = result["qc"]' in src:
            body = src[src.index('qc = result["qc"]'):]
            return body.replace("display(HTML(card))", "")
    raise AssertionError("the result card cell is not in the generated notebook")


@pytest.fixture(scope="module")
def card(rejected):
    env = {"result": rejected, "HTML": lambda x: x, "display": lambda x: x}
    exec(compile(_card_source(), "<card>", "exec"), env)
    return env["card"]


class TestResultCard:
    def test_the_badge_is_the_synthesis_verdict(self, card):
        """It read `qc["status"]`, so a refused design showed a green PASS."""
        assert "NOT ORDERABLE" in card
        assert not re.search(r'border-radius:99px">PASS<', card)

    def test_the_badge_colour_is_the_refusal_colour(self, card):
        """The tone came from the CDS status, so a refused design rendered green."""
        head = card[:card.index("border-radius:99px")]
        assert "#a8402c" in head, "the badge should carry the FAIL colour"
        assert "#1a7f5a" not in head, "the badge must not be the PASS green"

    def test_the_card_states_the_breach(self, card):
        assert "must not be ordered" in card
        assert "0.660" in card and "0.35, 0.65" in card

    def test_the_card_says_the_files_are_not_an_order(self, card):
        assert "diagnostic_" in card and "REJECTED" in card
        assert "not an order" in card

    def test_the_cds_result_is_kept_and_scoped(self, card):
        """The coding span genuinely passes. That stays, with its scope named."""
        assert "coding-span QC" in card
        assert "the CDS only" in card
        assert "GC, coding span" in card
        assert "final artifacts" in card
        assert "the sequences that would be ordered" in card

    def test_an_accepted_design_reads_as_accepted(self, rejected):
        """The positive side: the card must not say NOT ORDERABLE unconditionally."""
        accepted = dict(rejected)
        accepted["synthesis_ready"] = True
        accepted["not_ready_because"] = []
        accepted["final_artifacts"] = dict(rejected["final_artifacts"], status="pass")
        env = {"result": accepted, "HTML": lambda x: x, "display": lambda x: x}
        exec(compile(_card_source(), "<card>", "exec"), env)

        assert "READY TO ORDER" in env["card"]
        assert "must not be ordered" not in env["card"]


class TestAuditReport:
    def test_the_report_leads_with_the_refusal(self, rejected):
        text = rejected["audit"].report()
        assert "NOT ORDERABLE" in text
        assert text.index("NOT ORDERABLE") < text.index("coding-span QC")

    def test_the_report_gives_the_reason(self, rejected):
        text = rejected["audit"].report()
        assert "0.660" in text
        assert "final artifacts" in text

    def test_the_report_scopes_what_passed(self, rejected):
        """"constraints satisfied" and "QC PASS" are both true here and both narrow."""
        text = rejected["audit"].report()
        assert "coding-span QC PASS" in text
        assert "the CDS only, not the ordered sequences" in text
        assert "constraints satisfied  (the coding span)" in text

    def test_an_unscoped_qc_line_is_gone(self, rejected):
        text = rejected["audit"].report()
        assert not re.search(r"^  QC +PASS$", text, re.M)


class TestSummaryLine:
    def test_the_summary_states_the_verdict(self, rejected):
        """`r["summary"]` is the README quickstart's whole printed output."""
        assert "NOT ORDERABLE" in rejected["summary"]
        assert "Coding-span QC PASS" in rejected["summary"]

    def test_the_readme_quickstart_matches_production(self, rejected):
        """Check 3 parses this comment as the expected clean-install output."""
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        assert f"# {rejected['summary']}" in readme, (
            "the README quickstart output no longer matches what design_oneshot prints")


class TestNotebookWording:
    def test_a_fragment_row_is_not_called_an_order(self):
        text = NOTEBOOK.read_text(encoding="utf-8")
        assert "orderable piece" not in text
        assert "fragments to order" not in text

    def test_the_fragment_section_says_a_row_is_not_an_order(self):
        text = NOTEBOOK.read_text(encoding="utf-8")
        assert "A row here is not an order" in text
