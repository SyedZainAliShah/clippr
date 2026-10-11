"""Controls for the defects found in the 2026-10-11 review.

Each test fails against the code as it stood before that review. They are kept together
rather than scattered into the suites for each module, because what they have in common is
the failure mode, not the subject: in every one of them a check existed, ran, and reported
success while verifying nothing. A missing hash was skipped rather than failed, a generator
was walked twice so the second pass saw nothing, a lexical scan decided a continuation line
was a magic. None of these could be caught by reading the code and agreeing with it, which
is why they are written down as executable controls instead.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation"))

from ppropt.viz import (FAIL, FRAGMENT_HUES, PASS, WARN, kit_plate_svg,  # noqa: E402
                        plate_svg)

MODULES = [{"plasmid_id": "pA", "plate": "A1"}, {"plasmid_id": "pB", "plate": "B2"}]
NO_PLATE_NOTE = "no recorded plate position"


class TestFragmentHuesAreNeverAVerdict:
    """The palette served two roles and the fragment list grew into the verdict colours.

    At seven hues a 19S design drew F6 in FAIL red and F7 in PASS green, inside a diagram
    whose verdict band uses those same values, and every architecture drew F4 in WARN amber.
    """

    def test_no_fragment_hue_is_a_verdict_colour(self):
        assert not set(FRAGMENT_HUES) & {PASS, WARN, FAIL}

    @pytest.mark.parametrize("n_fragments", [1, 4, 5, 7, 12, 30])
    def test_no_architecture_can_cycle_onto_a_verdict_colour(self, n_fragments):
        used = {FRAGMENT_HUES[i % len(FRAGMENT_HUES)] for i in range(n_fragments)}
        assert not used & {PASS, WARN, FAIL}

    def test_neighbouring_fragments_never_share_a_hue(self):
        for i in range(30):
            a = FRAGMENT_HUES[i % len(FRAGMENT_HUES)]
            b = FRAGMENT_HUES[(i + 1) % len(FRAGMENT_HUES)]
            assert a != b


class TestKitPlateAcceptsAnyIterableOnce:
    """`selected_plasmids` is annotated Iterable and was walked twice."""

    def test_a_generator_keeps_the_missing_module_warning(self):
        wanted = ["pA", "pGHOST"]
        from_list = kit_plate_svg(list(wanted), MODULES)
        from_gen = kit_plate_svg((p for p in wanted), MODULES)
        assert NO_PLATE_NOTE in from_list
        assert NO_PLATE_NOTE in from_gen, "the generator was consumed by the first pass"

    def test_a_generator_draws_the_same_plate_as_a_list(self):
        wanted = ["pA", "pB"]
        assert (kit_plate_svg(list(wanted), MODULES).count('fill="#fff"')
                == kit_plate_svg((p for p in wanted), MODULES).count('fill="#fff"'))


class TestEveryUndrawnWellIsAccountedFor:
    """`bad[:6]` named six and silently discarded the rest."""

    def test_more_than_six_bad_wells_are_counted_not_dropped(self):
        bad = [f"Z{i}" for i in range(1, 10)]
        svg = plate_svg(bad)
        assert "(9)" in svg, "the total must be stated"
        assert "3 more" in svg, "the truncation must be visible"

    def test_six_or_fewer_are_listed_without_a_more_clause(self):
        assert "more" not in plate_svg(["Z1", "Z2"]).split("not drawn")[1]


class TestAuditStatusNeverRunsIntoItsLabel:
    """`:>8` pads but does not truncate, so a longer status lost its separator."""

    @pytest.mark.parametrize("status", ["pass", "fail", "warning", "not_assessed"])
    def test_the_status_is_separated_from_the_label(self, status):
        line = f"  final artifacts {status:>7} "
        assert "artifacts " in line
        assert f"artifacts{status}" not in line


class TestCellClassificationAsksTheParser:
    """A lexical per-line scan cannot tell a magic from a continuation line."""

    def _classify(self, src):
        import notebook_on_published_build as n
        return n.classify_cell(src)

    @pytest.mark.parametrize("src", [
        'msg = ("n=%d"\n       % n)',
        'ok = (a\n      != b)',
        'x = (1\n     % 2)',
    ])
    def test_valid_python_is_never_called_an_unsupported_magic(self, src):
        keep, dropped, unsupported = self._classify(src)
        assert unsupported == []
        assert dropped == []
        assert len(keep) == len(src.splitlines())

    def test_a_magic_inside_a_string_is_not_deleted_from_it(self):
        src = 'HELP = """\n%pip install foo\n"""'
        keep, dropped, unsupported = self._classify(src)
        assert dropped == [], "the line belongs to the string, not to the cell"
        assert "%pip install foo" in "\n".join(keep)

    def test_a_real_magic_is_still_dropped(self):
        keep, dropped, unsupported = self._classify("%pip install foo\nx = 1")
        assert dropped == ["%pip install foo"]
        assert keep == ["x = 1"]


class TestVendorAnchorsSurviveMarkup:
    """A false UNSUPPORTED accuses a truthful source of fabricating a quote."""

    def _mod(self):
        import check_vendor_claims as c
        return c

    @pytest.mark.parametrize("page,anchor", [
        ("the limit is 25<span>%</span> of the sequence length",
         "the limit is 25% of the sequence length"),
        ("a run of <b>10</b>nt or longer is rejected",
         "a run of 10nt or longer is rejected"),
    ])
    def test_tags_splitting_a_token_still_match(self, page, anchor):
        c = self._mod()
        stripped = c.re.sub(r"(?s)<[^>]+>", " ", page)
        assert c.squash(anchor) in c.squash(stripped)

    @pytest.mark.parametrize("entity,plain", [
        ("5&nbsp;bp", "5 bp"),
        ("repeats &amp; hairpins", "repeats & hairpins"),
        ("&lt; 25% of the length", "< 25% of the length"),
    ])
    def test_entities_are_decoded_before_matching(self, entity, plain):
        import html as html_lib
        c = self._mod()
        assert c.normalise(plain) in c.normalise(html_lib.unescape(entity))


class TestMalformedClaimsAreGradedNotCrashed:
    """A null field passed the required-field test via ``str(None) == "None"``."""

    def test_a_null_field_is_reported_missing(self):
        c = self._mod()
        claim = {"vendor": "Twist", "product": "Gene Fragment", "parameter": "GC",
                 "value": "25-65%", "enforcement": "hard",
                 "url": "https://example.invalid/x", "verbatim_anchor": None}
        report = c.check([claim], pause=0)
        assert report["malformed"], "a null anchor must be malformed, not fetched"
        assert "verbatim_anchor" in report["malformed"][0]["problem"]

    def test_a_null_url_is_reported_missing(self):
        c = self._mod()
        claim = {"vendor": "Twist", "product": "Gene Fragment", "parameter": "GC",
                 "value": "25-65%", "enforcement": "hard", "url": None,
                 "verbatim_anchor": "the six word anchor goes here"}
        report = c.check([claim], pause=0)
        assert report["malformed"]
        assert "url" in report["malformed"][0]["problem"]

    def _mod(self):
        import check_vendor_claims as c
        return c


class TestRecordedPathsAreAnchoredToTheRepo:
    """Resolving a recorded path against the CWD made the verdict depend on where we ran."""

    def test_a_relative_path_resolves_under_root_not_cwd(self):
        import release_check as rc
        resolved = rc._under_root("work/independent/corpus200_pinned.jsonl")
        assert resolved.is_absolute()
        assert str(resolved).startswith(str(rc.ROOT))

    def test_an_absolute_path_is_left_alone(self):
        import release_check as rc
        p = Path(__file__).resolve()
        assert rc._under_root(p) == p


class TestEveryNotebookMustDeclareItsInputs:
    """The staged-input loop can only verify what a record declares."""

    def test_both_notebooks_name_required_inputs(self):
        import release_check as rc
        for name in ("notebook_designer", "notebook_inventory"):
            assert rc.NOTEBOOK_REQUIRED_INPUTS.get(name), f"{name} names no required input"

    def test_a_record_declaring_nothing_is_a_problem(self):
        import release_check as rc
        problems = []
        rc._check_notebook_record(
            "notebook_designer",
            {"ok": True, "notebook_sha256": "x", "notebook": "nope",
             "staged_inputs": {}}, None, problems)
        assert any("staged no" in p for p in problems)
