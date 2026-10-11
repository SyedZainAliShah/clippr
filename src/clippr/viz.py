"""Inline SVG for the two things a bench scientist reads off a screen and acts on.

A plate map and an assembly diagram are not decoration. A list of wells is transcribed by
hand into a pipetting order, and a list of fragment coordinates is mentally reassembled to
check the overhangs line up. Both are error-prone in the same way: the information is present
and the reader has to do the work.

Two rules this module keeps, both learned the hard way elsewhere in this package:

**Nothing important is hover-only.** Tooltips do not survive export to PDF or static HTML, so
every fact a reader needs is drawn.

**Colours inherit.** Notebooks run on light and dark backgrounds. Text uses `currentColor` and
fills are translucent greys, so nothing renders as a white box on a dark page. The accent
colours are mid-tone on purpose: they read against either.
"""
from __future__ import annotations

from typing import Iterable, Mapping, Sequence

#: Mid-tone and chosen to stay legible on a light or a dark background.
PASS = "#1a7f5a"
WARN = "#9a6b1f"
FAIL = "#a8402c"
MUTED = "rgba(128,145,138,.45)"
FILL = "rgba(128,145,138,.16)"

#: Fragment hues, cycled. Distinct at a glance without being a traffic light -- these say
#: "different piece", not "good" or "bad", and must not be read as a verdict.
FRAGMENT_HUES = ("#3b6ea5", "#6b4fa0", "#1a7f5a", "#9a6b1f", "#8a4f6d", "#2d7d8a", "#7a6a3a")

ROWS = "ABCDEFGH"
COLUMNS = range(1, 13)


def _esc(text) -> str:
    return (str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _parse_well(well: str) -> tuple[str, int] | None:
    """`"C1"` -> `("C", 1)`. Returns None for anything that is not a 96-well address."""
    well = str(well).strip().upper()
    if len(well) < 2 or well[0] not in ROWS:
        return None
    try:
        column = int(well[1:])
    except ValueError:
        return None
    return (well[0], column) if column in COLUMNS else None


def plate_svg(selected: Iterable[str], *, present: Iterable[str] | None = None,
              title: str | None = None, note: str | None = None,
              accent: str = PASS) -> str:
    """A 96-well map with `selected` wells filled.

    `present` is the rest of the plate's occupied wells -- the kit's other modules -- drawn
    faintly so the reader can see the plate they are holding rather than a sparse scatter of
    dots. Wells in `selected` that are not valid 96-well addresses are listed under the plate
    instead of being dropped silently.
    """
    chosen, bad = {}, []
    for well in selected:
        parsed = _parse_well(well)
        (bad.append(str(well)) if parsed is None
         else chosen.setdefault(parsed, str(well).strip().upper()))
    occupied = {p for p in (_parse_well(w) for w in (present or ())) if p}

    pitch, radius = 38, 14
    left, top = 46, 52 if title else 34
    width = left + pitch * 12 + 24
    height = top + pitch * 8 + (44 if bad else 20)

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
             f'width="{width}" height="{height}" role="img" '
             f'style="max-width:100%;font-family:ui-sans-serif,system-ui,sans-serif">']
    if title:
        parts.append(f'<text x="{left - 24}" y="24" font-size="14" font-weight="600" '
                     f'fill="currentColor">{_esc(title)}</text>')

    # The plate outline, so it reads as a plate and not a grid of dots.
    parts.append(f'<rect x="{left - 26}" y="{top - 24}" width="{pitch * 12 + 18}" '
                 f'height="{pitch * 8 + 16}" rx="10" fill="none" stroke="{MUTED}"/>')

    for c in COLUMNS:
        parts.append(f'<text x="{left + pitch * (c - 1)}" y="{top - 6}" font-size="10" '
                     f'text-anchor="middle" fill="currentColor" opacity=".6">{c}</text>')
    for r, row in enumerate(ROWS):
        parts.append(f'<text x="{left - 20}" y="{top + pitch * r + 4}" font-size="10" '
                     f'text-anchor="middle" fill="currentColor" opacity=".6">{row}</text>')

    for r, row in enumerate(ROWS):
        for c in COLUMNS:
            cx, cy = left + pitch * (c - 1), top + pitch * r
            key = (row, c)
            if key in chosen:
                parts.append(f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="{accent}"/>')
                # The label goes inside the well, so a printed page still says which.
                parts.append(f'<text x="{cx}" y="{cy + 4}" font-size="9.5" font-weight="700" '
                             f'text-anchor="middle" fill="#fff">{row}{c}</text>')
            elif key in occupied:
                parts.append(f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="{FILL}" '
                             f'stroke="{MUTED}"/>')
            else:
                parts.append(f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" '
                             f'stroke="{MUTED}" stroke-dasharray="2 3"/>')

    caption = note or (f"{len(chosen)} well(s) to pick"
                       + (f", {len(occupied - set(chosen))} other kit modules shown faintly"
                          if occupied else ""))
    parts.append(f'<text x="{left - 26}" y="{top + pitch * 8 + 4}" font-size="11" '
                 f'fill="currentColor" opacity=".7">{_esc(caption)}</text>')
    if bad:
        parts.append(f'<text x="{left - 26}" y="{top + pitch * 8 + 22}" font-size="11" '
                     f'fill="{FAIL}">not a 96-well address, not drawn: '
                     f'{_esc(", ".join(bad[:6]))}</text>')
    parts.append("</svg>")
    return "".join(parts)


def kit_plate_svg(selected_plasmids: Iterable[str], modules: Sequence[Mapping],
                  **kwargs) -> str:
    """`plate_svg` driven by the GRASP module records, which carry `plasmid_id` and `plate`."""
    where = {m["plasmid_id"]: m["plate"] for m in modules if m.get("plate")}
    wanted = [where[p] for p in selected_plasmids if p in where]
    missing = [p for p in selected_plasmids if p not in where]
    note = kwargs.pop("note", None)
    if missing:
        note = (note or f"{len(wanted)} well(s) to pick") + \
            f"; {len(missing)} selected module(s) have no recorded plate position"
    return plate_svg(wanted, present=where.values(), note=note, **kwargs)


def assembly_svg(oligos, *, title: str | None = None, orderable: bool | None = None,
                 verdict: str | None = None) -> str:
    """The fragments end to end, with each junction's overhang drawn between them.

    `oligos` is the design's oligo table. `orderable` and `verdict` describe the **final
    artifact** check, not the coding-span QC: a design can pass the latter and fail the
    former, and drawing one verdict where the reader expects the other is a mistake this
    package has already made once.
    """
    rows = list(oligos.to_dict("records")) if hasattr(oligos, "to_dict") else list(oligos)
    if not rows:
        return ('<div style="font:13px ui-sans-serif,system-ui;opacity:.7">'
                'no fragments to draw</div>')

    spans = [(int(r.get("cds_start", 0)), int(r.get("cds_end", 0))) for r in rows]
    total = max(e for _, e in spans) or 1
    track, left, top = 760, 16, 62 if title else 40
    height = top + 108

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" '
             f'viewBox="0 0 {track + left * 2} {height}" width="{track + left * 2}" '
             f'height="{height}" role="img" '
             f'style="max-width:100%;font-family:ui-sans-serif,system-ui,sans-serif">']
    if title:
        parts.append(f'<text x="{left}" y="22" font-size="14" font-weight="600" '
                     f'fill="currentColor">{_esc(title)}</text>')
    if verdict:
        tone = FAIL if orderable is False else PASS if orderable else MUTED
        parts.append(f'<text x="{left}" y="40" font-size="11.5" font-weight="700" '
                     f'fill="{tone}">{_esc(verdict)}</text>')

    bar_y, bar_h = top + 16, 34
    for i, (row, (start, end)) in enumerate(zip(rows, spans)):
        x = left + track * start / total
        w = max(2.0, track * (end - start) / total)
        hue = FRAGMENT_HUES[i % len(FRAGMENT_HUES)]
        parts.append(f'<rect x="{x:.1f}" y="{bar_y}" width="{w:.1f}" height="{bar_h}" '
                     f'rx="4" fill="{hue}" opacity=".85"/>')

        ident = str(row.get("fragment_id") or row.get("order_fragment_id") or f"F{i + 1}")
        if w > 44:
            parts.append(f'<text x="{x + w / 2:.1f}" y="{bar_y + 22}" font-size="11.5" '
                         f'font-weight="700" text-anchor="middle" fill="#fff">'
                         f'{_esc(ident.split("_")[-1])}</text>')
        length = row.get("oligo_length") or (end - start)
        parts.append(f'<text x="{x + w / 2:.1f}" y="{bar_y + bar_h + 15}" font-size="10" '
                     f'text-anchor="middle" fill="currentColor" opacity=".65">'
                     f'{length} nt</text>')

        # The junction as both fragments see it on the coding strand. The physical 3' end
        # of one fragment is the reverse complement of the next fragment's 5' end -- F1 ends
        # GAAT where F2 begins ATTC -- so drawing the physical end labels one junction with
        # two different words. `oh3_coding_site_5to3` is the shared identity, and the one
        # `junction_overhangs` reports.
        exit_oh = row.get("oh3_coding_site_5to3") or row.get("three_prime_end_overhang")
        if exit_oh and i < len(rows) - 1:
            jx = x + w
            parts.append(f'<path d="M {jx - 7:.1f} {bar_y + 5} '
                         f'q 7 {bar_h / 2 - 5} 0 {bar_h - 10} '
                         f'l 14 0 q -7 -{bar_h / 2 - 5} 0 -{bar_h - 10} z" '
                         f'fill="currentColor" opacity=".28"/>')
            parts.append(f'<text x="{jx:.1f}" y="{bar_y - 7}" font-size="10.5" '
                         f'font-family="ui-monospace,monospace" font-weight="700" '
                         f'text-anchor="middle" fill="currentColor">{_esc(exit_oh)}</text>')

    ends = (rows[0].get("oh5_coding_site_5to3") or rows[0].get("five_prime_end_overhang"),
            rows[-1].get("oh3_coding_site_5to3") or rows[-1].get("three_prime_end_overhang"))
    if ends[0]:
        parts.append(f'<text x="{left}" y="{bar_y - 7}" font-size="10.5" '
                     f'font-family="ui-monospace,monospace" text-anchor="start" '
                     f'fill="currentColor" opacity=".75">{_esc(ends[0])}</text>')
    if ends[1]:
        parts.append(f'<text x="{left + track}" y="{bar_y - 7}" font-size="10.5" '
                     f'font-family="ui-monospace,monospace" text-anchor="end" '
                     f'fill="currentColor" opacity=".75">{_esc(ends[1])}</text>')

    enzyme = rows[0].get("wrap_enzyme") or "BsaI"
    parts.append(f'<text x="{left}" y="{height - 8}" font-size="11" fill="currentColor" '
                 f'opacity=".65">{len(rows)} fragments, {total} nt coding sequence. '
                 f'Four-base {_esc(enzyme)} junctions are labelled on the coding strand, '
                 f'not as the physical single-stranded end. Lengths are the ordered '
                 f'oligos, wrapper included.</text>')
    parts.append("</svg>")
    return "".join(parts)
