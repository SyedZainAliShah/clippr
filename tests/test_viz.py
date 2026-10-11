"""Controls for the inline SVG: the plate map and the assembly diagram.

These draw things a bench scientist acts on. A wrong well is a wrong pipette; a wrong
overhang label is a failed ligation. So the tests check the content of the drawing, not that
a drawing was produced.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from clippr.viz import (FAIL, FILL_OPACITY, MUTED, PASS, WELL_R,
                        assembly_svg, kit_plate_svg, plate_svg)

PARTS = Path(__file__).resolve().parents[1] / "src" / "clippr" / "parts.json"


@pytest.fixture(scope="module")
def modules():
    return json.loads(PARTS.read_text(encoding="utf-8"))["modules"]


def well_circles(svg: str) -> list[str]:
    """Only the wells. The legend draws swatches too, at a different radius."""
    return re.findall(rf'<circle cx="[\d.]+" cy="[\d.]+" r="{WELL_R}"', svg)


def wells_drawn_filled(svg: str) -> set[str]:
    """The wells whose label is drawn inside a filled circle."""
    return set(re.findall(r">([A-H]\d{1,2})</text>", svg))


class TestPlateMap:
    def test_every_well_of_a_96_plate_is_drawn(self):
        assert len(well_circles(plate_svg(["A1"]))) == 96

    def test_selected_wells_are_labelled_in_place(self):
        svg = plate_svg(["C1", "B5", "F12"])
        assert wells_drawn_filled(svg) == {"C1", "B5", "F12"}

    def test_the_label_is_inside_the_well_not_a_legend(self):
        """A printed page has no hover. The well must say which well it is."""
        svg = plate_svg(["C1"])
        assert ">C1</text>" in svg

    def test_an_unselected_well_carries_no_label(self):
        svg = plate_svg(["C1"])
        assert "B5" not in wells_drawn_filled(svg)

    def test_an_invalid_address_is_reported_not_dropped(self):
        """Silently omitting a well would be a missing pipette step."""
        svg = plate_svg(["C1", "ZZ9", "Q13", "A99"])
        assert "not a 96-well address" in svg
        for bad in ("ZZ9", "Q13", "A99"):
            assert bad in svg
        assert wells_drawn_filled(svg) == {"C1"}

    def test_lowercase_and_padding_are_accepted(self):
        assert wells_drawn_filled(plate_svg([" c1 ", "b5"])) == {"C1", "B5"}

    def test_the_rest_of_the_kit_is_drawn_faintly(self, modules):
        """The reader is holding a populated plate, not a sparse scatter of dots."""
        svg = kit_plate_svg([m["plasmid_id"] for m in modules[:5]], modules)
        assert wells_drawn_filled(svg) == {m["plate"] for m in modules[:5]}
        faint = svg.count(f'r="{WELL_R}" fill="{MUTED}" fill-opacity="{FILL_OPACITY}"')
        assert faint == len({m["plate"] for m in modules}) - 5

    def test_the_legend_names_each_kind_of_well(self):
        """Three fills mean three things; prose in a caption is easy to skip."""
        svg = plate_svg(["C1"], present=["B5"])
        for label in ("to pick", "other kit module", "empty"):
            assert f">{label}</text>" in svg

    def test_the_legend_can_be_turned_off(self):
        """The swatches go; the caption, which also says "to pick", stays."""
        assert 'r="6"' not in plate_svg(["C1"], legend=False)
        assert 'r="6"' in plate_svg(["C1"], legend=True)

    def test_a_module_with_no_plate_position_is_named(self, modules):
        svg = kit_plate_svg(["pPR-NOT-IN-THE-KIT"], modules)
        assert "no recorded plate position" in svg

    def test_greys_follow_the_page_rather_than_being_fixed(self):
        """A fixed grey cannot be dark on white and light on black. currentColor can."""
        svg = plate_svg(["C1"], present=["B5"], title="t")
        assert "currentColor" in svg
        assert "rgba(" not in svg, "no baked-in grey should survive"
        assert "#ffffff" not in svg.lower().replace('fill="#fff"', "")


@pytest.fixture(scope="module")
def design():
    from clippr.codons import table_from_kazusa
    from clippr.design import design_oneshot

    return design_oneshot("AAAAUGUGG", codon_table=table_from_kazusa(),
                          check_offtarget=False)


class TestAssemblyDiagram:

    def test_one_rect_per_fragment(self, design):
        svg = assembly_svg(design["oligos"])
        assert svg.count("<rect") == len(design["oligos"])

    def test_the_junction_label_is_the_coding_strand_identity(self, design):
        """The physical 3' end of one fragment is the reverse complement of the next
        fragment's 5' end -- F1 ends GAAT where F2 begins ATTC. Labelling the physical end
        gives one junction two different words, and `junction_overhangs` reports neither.
        """
        svg = assembly_svg(design["oligos"])
        for overhang in design["junction_overhangs"]:
            assert overhang in svg, f"junction {overhang} is not drawn"

    def test_the_caption_says_which_strand_is_labelled(self, design):
        """Ambiguity here is a mis-ordered oligo, so it is spelled out rather than implied."""
        svg = assembly_svg(design["oligos"])
        assert "coding strand" in svg
        assert "not as the physical single-stranded end" in svg

    def test_a_tab_is_drawn_at_every_internal_junction(self, design):
        svg = assembly_svg(design["oligos"])
        assert svg.count("<path") == len(design["oligos"]) - 1

    def test_the_final_verdict_is_shown_when_given(self, design):
        svg = assembly_svg(design["oligos"], orderable=False, verdict="NOT ORDERABLE")
        assert "NOT ORDERABLE" in svg
        assert FAIL in svg, "a refusal must not be drawn in the pass colour"
        assert PASS not in svg

    def test_an_orderable_design_is_not_drawn_as_a_refusal(self, design):
        svg = assembly_svg(design["oligos"], orderable=True, verdict="Orderable")
        assert FAIL not in svg

    def test_lengths_are_the_ordered_oligo_not_the_coding_slice(self, design):
        """The wrapper is part of what gets synthesised and paid for."""
        svg = assembly_svg(design["oligos"])
        for length in design["oligos"]["oligo_length"]:
            assert f">{length} nt</text>" in svg

    def test_an_empty_table_says_so_rather_than_drawing_an_empty_axis(self):
        assert "no fragments to draw" in assembly_svg([])

    def test_it_accepts_plain_dictionaries(self, design):
        """The notebooks hold records, not always a frame."""
        rows = design["oligos"].to_dict("records")
        assert assembly_svg(rows).count("<rect") == len(rows)

    def test_colours_inherit(self, design):
        svg = assembly_svg(design["oligos"], title="t")
        assert "currentColor" in svg
