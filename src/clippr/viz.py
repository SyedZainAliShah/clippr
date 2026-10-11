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

#: Fully saturated, and still legible on both backgrounds -- the two are usually traded
#: against each other, so the palette was solved rather than chosen. Each hue is the most
#: saturated sRGB colour at its angle that keeps white label text at WCAG AA (>= 4.5:1) and
#: stays visible on a dark page (>= 4.0:1). The previous muted set managed 4.7-6.5 on white
#: but only 2.9-4.1 on dark, which is why it looked drab and washed out at once.
#:
#: One palette serves both roles. These work as a fill under white text and as text on either
#: background, so a verdict and a well use the same red and mean the same thing.
PASS = "#008554"
WARN = "#ad6500"
FAIL = "#eb0046"

#: Greys follow the page instead of fighting it: `currentColor` inherits the notebook's text
#: colour, so the same markup gives dark lines on white and light lines on black. A fixed
#: grey cannot do both.
MUTED = "currentColor"
MUTED_OPACITY = ".32"
FILL_OPACITY = ".14"

#: Fragment hues, cycled. Distinct at a glance without being a traffic light -- these say
#: "different piece", not "good" or "bad", and must not be read as a verdict. Deliberately
#: no green or red among the first few, so a fragment is never mistaken for a verdict.
FRAGMENT_HUES = ("#056dff", "#9233ff", "#007e9e", "#ad6500", "#558000", "#eb0046", "#008554")

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


#: Plate geometry, derived rather than typed, so the margins are equal by construction.
#: The first version had a left gap of 32 against a right gap of 48 and frame insets of
#: 12/16/10/16; nobody notices that in a picture, which is why it is arithmetic now.
PAD = 18            # outer margin, all four sides
WELL_R = 14         # well radius
PITCH = 38          # centre-to-centre
FRAME_INSET = 13    # frame edge to the outermost well edge, all four sides
AXIS_GAP = 17       # well edge to its row letter or column number
TITLE_H = 26        # title band, when there is a title
CAPTION_H = 24      # caption band


def plate_svg(selected: Iterable[str], *, present: Iterable[str] | None = None,
              title: str | None = None, note: str | None = None,
              accent: str = PASS, legend: bool = True) -> str:
    """A 96-well map with `selected` wells filled.

    `present` is the rest of the plate's occupied wells -- the kit's other modules -- drawn
    faintly, so the reader sees the plate they are holding rather than a sparse scatter of
    dots. Wells in `selected` that are not valid 96-well addresses are listed under the
    plate instead of being dropped silently: a missing well is a missing pipette step.
    """
    chosen, bad = {}, []
    for well in selected:
        parsed = _parse_well(well)
        (bad.append(str(well)) if parsed is None
         else chosen.setdefault(parsed, str(well).strip().upper()))
    occupied = {p for p in (_parse_well(w) for w in (present or ())) if p}

    # Left to right: pad, row letters, frame inset, wells, frame inset, pad.
    axis_w = AXIS_GAP + 9
    frame_x = PAD + axis_w
    first_cx = frame_x + FRAME_INSET + WELL_R
    last_cx = first_cx + PITCH * 11
    frame_w = (last_cx + WELL_R + FRAME_INSET) - frame_x
    width = frame_x + frame_w + PAD

    # Top to bottom: pad, title, column numbers, frame inset, wells, frame inset, caption.
    head = PAD + (TITLE_H if title else 0)
    frame_y = head + axis_w
    first_cy = frame_y + FRAME_INSET + WELL_R
    last_cy = first_cy + PITCH * 7
    frame_h = (last_cy + WELL_R + FRAME_INSET) - frame_y
    legend_h = 26 if legend else 0
    height = frame_y + frame_h + CAPTION_H + legend_h + (18 if bad else 0)

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
           f'width="{width}" height="{height}" role="img" '
           f'style="max-width:100%;font-family:ui-sans-serif,system-ui,sans-serif">']
    if title:
        out.append(f'<text x="{PAD}" y="{PAD + 15}" font-size="14" font-weight="600" '
                   f'fill="currentColor">{_esc(title)}</text>')

    out.append(f'<rect x="{frame_x}" y="{frame_y}" width="{frame_w}" height="{frame_h}" '
               f'rx="10" fill="none" stroke="{MUTED}" stroke-opacity="{MUTED_OPACITY}"/>')

    for c in COLUMNS:
        out.append(f'<text x="{first_cx + PITCH * (c - 1)}" y="{frame_y - 7}" '
                   f'font-size="10" text-anchor="middle" fill="currentColor" '
                   f'opacity=".6">{c}</text>')
    for r, row in enumerate(ROWS):
        out.append(f'<text x="{frame_x - 9}" y="{first_cy + PITCH * r + 4}" font-size="10" '
                   f'text-anchor="middle" fill="currentColor" opacity=".6">{row}</text>')

    for r, row in enumerate(ROWS):
        for c in COLUMNS:
            cx, cy = first_cx + PITCH * (c - 1), first_cy + PITCH * r
            key = (row, c)
            if key in chosen:
                out.append(f'<circle cx="{cx}" cy="{cy}" r="{WELL_R}" fill="{accent}"/>')
                # Inside the well, because a printed protocol has no hover.
                out.append(f'<text x="{cx}" y="{cy + 4}" font-size="9.5" font-weight="700" '
                           f'text-anchor="middle" fill="#fff">{row}{c}</text>')
            elif key in occupied:
                out.append(f'<circle cx="{cx}" cy="{cy}" r="{WELL_R}" fill="{MUTED}" '
                           f'fill-opacity="{FILL_OPACITY}" stroke="{MUTED}" '
                           f'stroke-opacity="{MUTED_OPACITY}"/>')
            else:
                out.append(f'<circle cx="{cx}" cy="{cy}" r="{WELL_R}" fill="none" '
                           f'stroke="{MUTED}" stroke-opacity="{MUTED_OPACITY}" '
                           f'stroke-dasharray="2 3"/>')

    y = frame_y + frame_h + 17
    out.append(f'<text x="{PAD}" y="{y}" font-size="11" fill="currentColor" opacity=".75">'
               f'{_esc(note or f"{len(chosen)} well(s) to pick")}</text>')
    if legend:
        y += legend_h
        swatches = [(accent, "to pick", None),
                    (MUTED, "other kit module", MUTED),
                    ("none", "empty", MUTED)]
        x = PAD + 6
        for fill, label, stroke in swatches:
            dash = ' stroke-dasharray="2 3"' if label == "empty" else ""
            opacity = (f' fill-opacity="{FILL_OPACITY}"' if fill == MUTED else "")
            out.append(f'<circle cx="{x}" cy="{y - 4}" r="6" fill="{fill}"{opacity}'
                       + (f' stroke="{stroke}" stroke-opacity="{MUTED_OPACITY}"{dash}'
                          if stroke else "") + '/>')
            out.append(f'<text x="{x + 11}" y="{y}" font-size="10.5" fill="currentColor" '
                       f'opacity=".7">{_esc(label)}</text>')
            x += 20 + 7 * len(label)
    if bad:
        out.append(f'<text x="{PAD}" y="{height - 6}" font-size="11" fill="{FAIL}">'
                   f'not a 96-well address, not drawn: {_esc(", ".join(bad[:6]))}</text>')
    out.append("</svg>")
    return "".join(out)


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


#: Assembly geometry, same discipline: every band derived from one pad.
TRACK = 760
BAR_H = 34
BAND = 20           # vertical rhythm between bands


def assembly_svg(oligos, *, title: str | None = None, orderable: bool | None = None,
                 verdict: str | None = None) -> str:
    """The fragments end to end, with each junction's overhang drawn between them.

    `oligos` is the design's oligo table. `orderable` and `verdict` describe the **final
    artifact** check, not the coding-span QC: a design can pass the latter and fail the
    former, and drawing one where the reader expects the other is a mistake this package has
    already made once.
    """
    rows = list(oligos.to_dict("records")) if hasattr(oligos, "to_dict") else list(oligos)
    if not rows:
        return ('<div style="font:13px ui-sans-serif,system-ui;opacity:.7">'
                'no fragments to draw</div>')

    spans = [(int(r.get("cds_start", 0) or 0), int(r.get("cds_end", 0) or 0)) for r in rows]
    total = max((e for _, e in spans), default=0) or 1

    # Bands: pad, title, verdict, overhang labels, bar, lengths, ids, caption, pad.
    y = PAD
    title_y = y + 15 if title else None
    if title:
        y += TITLE_H
    verdict_y = y + 12 if verdict else None
    if verdict:
        y += BAND
    oh_y = y + 11
    bar_y = y + BAND
    len_y = bar_y + BAR_H + 15
    id_y = len_y + 14
    caption_y = id_y + BAND
    height = caption_y + PAD - 6
    width = TRACK + PAD * 2

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
           f'width="{width}" height="{height}" role="img" '
           f'style="max-width:100%;font-family:ui-sans-serif,system-ui,sans-serif">']
    if title:
        out.append(f'<text x="{PAD}" y="{title_y}" font-size="14" font-weight="600" '
                   f'fill="currentColor">{_esc(title)}</text>')
    if verdict:
        tone = FAIL if orderable is False else PASS if orderable else "currentColor"
        out.append(f'<text x="{PAD}" y="{verdict_y}" font-size="11.5" font-weight="700" '
                   f'fill="{tone}">{_esc(verdict)}</text>')

    for i, (row, (start_nt, end_nt)) in enumerate(zip(rows, spans)):
        x = PAD + TRACK * start_nt / total
        w = max(2.0, TRACK * (end_nt - start_nt) / total)
        hue = FRAGMENT_HUES[i % len(FRAGMENT_HUES)]
        out.append(f'<rect x="{x:.1f}" y="{bar_y}" width="{w:.1f}" height="{BAR_H}" '
                   f'rx="4" fill="{hue}" opacity=".85"/>')

        ident = str(row.get("fragment_id") or row.get("order_fragment_id")
                    or f"F{i + 1}").split("_")[-1]
        # A narrow fragment gets its label underneath rather than losing it. A 19-mer
        # splits into seven, and the middle ones are too thin to write inside.
        if w > 11 * len(ident) + 10:
            out.append(f'<text x="{x + w / 2:.1f}" y="{bar_y + 22}" font-size="11.5" '
                       f'font-weight="700" text-anchor="middle" fill="#fff">'
                       f'{_esc(ident)}</text>')
        else:
            out.append(f'<text x="{x + w / 2:.1f}" y="{id_y}" font-size="10" '
                       f'font-weight="700" text-anchor="middle" fill="{hue}">'
                       f'{_esc(ident)}</text>')

        length = row.get("oligo_length") or (end_nt - start_nt)
        out.append(f'<text x="{x + w / 2:.1f}" y="{len_y}" font-size="10" '
                   f'text-anchor="middle" fill="currentColor" opacity=".65">'
                   f'{length} nt</text>')

        # The junction as both fragments see it on the coding strand. The physical 3' end
        # of one fragment is the reverse complement of the next fragment's 5' end -- F1
        # ends GAAT where F2 begins ATTC -- so drawing the physical end labels one junction
        # with two different words. This is the shared identity, and the one
        # `junction_overhangs` reports.
        exit_oh = row.get("oh3_coding_site_5to3") or row.get("three_prime_end_overhang")
        if exit_oh and i < len(rows) - 1:
            jx = x + w
            out.append(f'<path d="M {jx - 7:.1f} {bar_y + 5} '
                       f'q 7 {BAR_H / 2 - 5} 0 {BAR_H - 10} l 14 0 '
                       f'q -7 -{BAR_H / 2 - 5} 0 -{BAR_H - 10} z" '
                       f'fill="currentColor" opacity=".28"/>')
            out.append(f'<text x="{jx:.1f}" y="{oh_y}" font-size="10.5" '
                       f'font-family="ui-monospace,monospace" font-weight="700" '
                       f'text-anchor="middle" fill="currentColor">{_esc(exit_oh)}</text>')

    ends = (rows[0].get("oh5_coding_site_5to3") or rows[0].get("five_prime_end_overhang"),
            rows[-1].get("oh3_coding_site_5to3") or rows[-1].get("three_prime_end_overhang"))
    if ends[0]:
        out.append(f'<text x="{PAD}" y="{oh_y}" font-size="10.5" '
                   f'font-family="ui-monospace,monospace" text-anchor="start" '
                   f'fill="currentColor" opacity=".75">{_esc(ends[0])}</text>')
    if ends[1]:
        out.append(f'<text x="{PAD + TRACK}" y="{oh_y}" font-size="10.5" '
                   f'font-family="ui-monospace,monospace" text-anchor="end" '
                   f'fill="currentColor" opacity=".75">{_esc(ends[1])}</text>')

    enzyme = rows[0].get("wrap_enzyme") or "BsaI"
    out.append(f'<text x="{PAD}" y="{caption_y}" font-size="11" fill="currentColor" '
               f'opacity=".65">{len(rows)} fragments, {total} nt coding sequence. '
               f'Four-base {_esc(enzyme)} junctions are labelled on the coding strand, '
               f'not as the physical single-stranded end. Lengths are the ordered oligos, '
               f'wrapper included.</text>')
    out.append("</svg>")
    return "".join(out)
