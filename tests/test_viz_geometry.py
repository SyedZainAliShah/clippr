"""Geometry controls: the drawing is well-formed, symmetric, and survives its edge cases.

Visual bugs are the ones that get shipped, because a picture always looks like a picture.
The first version of this module had a left margin of 32 against a right margin of 48 and
frame insets of 12/16/10/16, and it looked fine. So symmetry is arithmetic here, not
judgement, and every coordinate is checked against the viewBox rather than eyeballed.
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET

import pytest

from ppropt.viz import (FRAME_INSET, PAD, PITCH, TRACK, WELL_R, assembly_svg,
                        plate_svg)

SVG_NS = "{http://www.w3.org/2000/svg}"


def tree(svg: str):
    """Parsed, which also proves the markup is well-formed."""
    return ET.fromstring(svg)


def viewbox(svg: str) -> tuple[float, float, float, float]:
    return tuple(float(v) for v in
                 re.search(r'viewBox="([^"]+)"', svg).group(1).split())


def circles(svg: str, radius: float | None = None):
    out = []
    for el in tree(svg).iter(f"{SVG_NS}circle"):
        r = float(el.get("r"))
        if radius is None or r == radius:
            out.append((float(el.get("cx")), float(el.get("cy")), r))
    return out


def rects(svg: str):
    return [(float(e.get("x")), float(e.get("y")),
             float(e.get("width")), float(e.get("height")))
            for e in tree(svg).iter(f"{SVG_NS}rect")]


# --------------------------------------------------------------------- plate

PLATE_CASES = [
    pytest.param([], id="nothing selected"),
    pytest.param(["A1"], id="one corner"),
    pytest.param(["H12"], id="other corner"),
    pytest.param(["A1", "A12", "H1", "H12"], id="four corners"),
    pytest.param([f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)], id="all 96"),
    pytest.param(["C1", "C1", "C1"], id="duplicates"),
    pytest.param([" c1 ", "B5"], id="padding and case"),
]


class TestPlateGeometry:
    @pytest.mark.parametrize("selected", PLATE_CASES)
    def test_it_is_well_formed_xml(self, selected):
        tree(plate_svg(selected, title="t"))

    @pytest.mark.parametrize("selected", PLATE_CASES)
    def test_the_lattice_is_always_96_wells_on_an_exact_grid(self, selected):
        wells = circles(plate_svg(selected), radius=WELL_R)
        assert len(wells) == 96
        xs, ys = sorted({w[0] for w in wells}), sorted({w[1] for w in wells})
        assert len(xs) == 12 and len(ys) == 8
        assert all(round(b - a, 6) == PITCH for a, b in zip(xs, xs[1:]))
        assert all(round(b - a, 6) == PITCH for a, b in zip(ys, ys[1:]))

    @pytest.mark.parametrize("selected", PLATE_CASES)
    def test_nothing_is_drawn_outside_the_viewbox(self, selected):
        svg = plate_svg(selected, title="a title", note="a note")
        _, _, w, h = viewbox(svg)
        for cx, cy, r in circles(svg):
            assert 0 <= cx - r and cx + r <= w, f"circle escapes horizontally at {cx}"
            assert 0 <= cy - r and cy + r <= h, f"circle escapes vertically at {cy}"
        for x, y, rw, rh in rects(svg):
            assert 0 <= x and x + rw <= w
            assert 0 <= y and y + rh <= h

    def test_the_outer_margins_are_equal_left_and_right(self):
        svg = plate_svg(["A1"], title="t")
        _, _, w, _ = viewbox(svg)
        frame = rects(svg)[0]
        left = frame[0] - PAD                       # row letters live in this band
        right = w - (frame[0] + frame[2])
        assert right == PAD, f"right margin {right} should be the pad {PAD}"
        assert left >= 0

    def test_the_frame_insets_are_equal_on_all_four_sides(self):
        """12/16/10/16 is what this looked like before it was derived."""
        svg = plate_svg(["A1"])
        fx, fy, fw, fh = rects(svg)[0]
        wells = circles(svg, radius=WELL_R)
        xs, ys = [c[0] for c in wells], [c[1] for c in wells]
        assert min(xs) - WELL_R - fx == FRAME_INSET
        assert (fx + fw) - (max(xs) + WELL_R) == FRAME_INSET
        assert min(ys) - WELL_R - fy == FRAME_INSET
        assert (fy + fh) - (max(ys) + WELL_R) == FRAME_INSET

    def test_the_grid_sits_centred_inside_its_frame(self):
        svg = plate_svg(["A1"])
        fx, fy, fw, fh = rects(svg)[0]
        wells = circles(svg, radius=WELL_R)
        xs, ys = [c[0] for c in wells], [c[1] for c in wells]
        assert abs(((min(xs) + max(xs)) / 2) - (fx + fw / 2)) < 1e-6
        assert abs(((min(ys) + max(ys)) / 2) - (fy + fh / 2)) < 1e-6

    def test_a_title_moves_the_plate_down_and_changes_nothing_else(self):
        plain, titled = plate_svg(["A1"]), plate_svg(["A1"], title="t")
        assert viewbox(plain)[2] == viewbox(titled)[2]
        shift = viewbox(titled)[3] - viewbox(plain)[3]
        assert shift > 0
        a = sorted(c[1] for c in circles(plain, WELL_R))
        b = sorted(c[1] for c in circles(titled, WELL_R))
        assert all(round(y2 - y1, 6) == shift for y1, y2 in zip(a, b))

    def test_every_selected_well_is_filled_and_no_other(self):
        """Read the attribute rather than matching the serialised text."""
        svg = plate_svg(["C1", "B5"])
        solid = [el for el in tree(svg).iter(f"{SVG_NS}circle")
                 if float(el.get("r")) == WELL_R
                 and (el.get("fill") or "").startswith("#")]
        assert len(solid) == 2

    def test_all_96_selected_still_fits(self):
        svg = plate_svg([f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)])
        assert len(re.findall(r'fill="#fff"', svg)) == 96

    def test_duplicates_are_drawn_once(self):
        assert plate_svg(["C1", "C1", "C1"]).count(">C1</text>") == 1


# ------------------------------------------------------------------ assembly

def fake_oligos(n: int, lengths: list[int] | None = None):
    """n fragments tiling a coding sequence, with junction overhangs."""
    lengths = lengths or [300] * n
    rows, start = [], 0
    for i, span in enumerate(lengths):
        rows.append({"fragment_id": f"T_F{i + 1}", "cds_start": start,
                     "cds_end": start + span, "oligo_length": span + 32,
                     "oh5_coding_site_5to3": "AAAA" if i else "CTCA",
                     "oh3_coding_site_5to3": "GGGG" if i < n - 1 else "CGAG",
                     "wrap_enzyme": "BsaI"})
        start += span
    return rows


ASSEMBLY_CASES = [
    pytest.param(fake_oligos(1), id="single fragment"),
    pytest.param(fake_oligos(4), id="9S, four fragments"),
    pytest.param(fake_oligos(7), id="19S, seven fragments"),
    pytest.param(fake_oligos(3, [10, 1500, 10]), id="wildly uneven"),
    pytest.param(fake_oligos(12), id="more than the hue palette"),
]


class TestAssemblyGeometry:
    @pytest.mark.parametrize("rows", ASSEMBLY_CASES)
    def test_it_is_well_formed_xml(self, rows):
        tree(assembly_svg(rows, title="t", verdict="v", orderable=False))

    @pytest.mark.parametrize("rows", ASSEMBLY_CASES)
    def test_the_fragments_tile_the_track_exactly(self, rows):
        bars = sorted(rects(assembly_svg(rows)), key=lambda r: r[0])
        assert round(bars[0][0], 6) == PAD
        assert round(bars[-1][0] + bars[-1][2], 6) == PAD + TRACK
        for a, b in zip(bars, bars[1:]):
            assert abs(b[0] - (a[0] + a[2])) < 0.2, "fragments must not gap or overlap"

    @pytest.mark.parametrize("rows", ASSEMBLY_CASES)
    def test_the_left_and_right_margins_match(self, rows):
        svg = assembly_svg(rows)
        _, _, w, _ = viewbox(svg)
        bars = sorted(rects(svg), key=lambda r: r[0])
        assert round(bars[0][0], 6) == round(w - (bars[-1][0] + bars[-1][2]), 6)

    @pytest.mark.parametrize("rows", ASSEMBLY_CASES)
    def test_nothing_escapes_the_viewbox(self, rows):
        svg = assembly_svg(rows, title="t", verdict="a fairly long verdict line")
        _, _, w, h = viewbox(svg)
        for x, y, rw, rh in rects(svg):
            assert 0 <= x and x + rw <= w + 1e-6
            assert 0 <= y and y + rh <= h

    @pytest.mark.parametrize("rows", ASSEMBLY_CASES)
    def test_every_fragment_keeps_its_identifier(self, rows):
        """A thin fragment loses its inside label; it must gain one underneath."""
        svg = assembly_svg(rows)
        for i in range(len(rows)):
            assert f">F{i + 1}</text>" in svg, f"F{i + 1} is unlabelled"

    @pytest.mark.parametrize("rows", ASSEMBLY_CASES)
    def test_one_tab_per_internal_junction(self, rows):
        svg = assembly_svg(rows)
        assert len(list(tree(svg).iter(f"{SVG_NS}path"))) == len(rows) - 1

    def test_a_single_fragment_has_no_junction(self):
        assert "<path" not in assembly_svg(fake_oligos(1))

    def test_hues_cycle_rather_than_running_out(self):
        svg = assembly_svg(fake_oligos(12))
        assert len(rects(svg)) == 12

    def test_a_zero_length_span_still_draws_something_clickable(self):
        rows = [{"fragment_id": "F1", "cds_start": 0, "cds_end": 0, "oligo_length": 0}]
        assert rects(assembly_svg(rows))[0][2] >= 2.0

    def test_markup_in_a_fragment_name_is_escaped(self):
        rows = fake_oligos(1)
        rows[0]["fragment_id"] = "<script>alert(1)</script>"
        svg = assembly_svg(rows)
        assert "<script>" not in svg
        assert "&lt;script&gt;" in svg

    def test_a_verdict_band_shifts_the_bar_and_nothing_else(self):
        plain = assembly_svg(fake_oligos(4))
        with_v = assembly_svg(fake_oligos(4), verdict="NOT ORDERABLE", orderable=False)
        assert viewbox(plain)[2] == viewbox(with_v)[2]
        assert viewbox(with_v)[3] > viewbox(plain)[3]
