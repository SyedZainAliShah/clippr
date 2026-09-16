"""V5 — digest and ligate CLIPPR's exported oligos from first principles, and check the product.

This deliberately shares nothing with the code it checks. It does not import
`clippr.assembly`, never calls `split_cds`, `build_oligos` or `reassemble`, and does not read
the `payload_5to3` column -- reading the payload would be reading the answer. It takes the
`oligo_sequence_5to3` strings a user would actually order, finds the Type IIS sites in them,
computes the cuts, recovers the sticky ends, joins compatible ends, and compares the product
with the coding sequence the design claims.

**Enzyme definition, frozen here rather than imported.** BsaI recognises GGTCTC and cuts the
top strand 1 nt downstream of the site and the bottom strand 5 nt downstream, leaving a
four-base 5' overhang: GGTCTC(1/5). The same site read on the other strand appears as
GAGACC. This is stated explicitly because taking it from `Bio.Restriction` -- which the
production code also uses -- would make a shared dependency look like independent
corroboration.

**Scope.** This simulates double-stranded fragments joining by complementary four-base 5'
overhangs. It does not model single-stranded pool synthesis, ligase kinetics, or any yield.
It establishes that the exported sequences encode the intended product under the stated
model, not that an assembly would succeed at the bench.

    python validation/independent/v5_digest_ligate.py [--limit N]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from journal import load_designs, provenance                       # noqa: E402

#: recognition site -> (bases from site end to the top-strand cut, overhang length)
ENZYMES = {"BsaI": ("GGTCTC", 1, 4), "BbsI": ("GAAGAC", 2, 4), "SapI": ("GCTCTTC", 1, 3)}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

#: Destination overhangs as *coding sites*, 5' then 3'. Declared in
#: `assembly_spec.DESTINATION_OVERHANGS`; written out here rather than imported so the
#: expected product is a second statement of the contract, not a restatement of the code.
#:
#: This is what makes the check an equality rather than a containment. `cds in product` was
#: the earlier predicate, and it accepted a design whose entry overhang had been changed from
#: CTCA to AAAA, and one with a GGTCTC inserted just inside the entry overhang: the CDS was
#: still in there, so both passed. An assembly that puts the right CDS in the wrong vector
#: context is not a correct assembly.
DESTINATIONS = {
    "level_minus1": ("ACAT", "TTGT"),
    "level0": ("CTCA", "CGAG"),
    "level1": ("GGAG", "CGCT"),
}
DEFAULT_DESTINATION = "level0"

#: No released fragment may carry another recognition site: it would be cut again, and the
#: ordered oligo would not yield the fragment the design assumes.
FORBIDDEN = ("BsaI", "BbsI", "SapI")


def rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


class DigestError(Exception):
    """The oligo does not present the two sites a Type IIS release requires."""


def digest(oligo: str, enzyme: str = "BsaI") -> tuple[str, str, str]:
    """Release the internal fragment, returning (5' overhang, core, 3' overhang).

    The oligo must present exactly one forward site and exactly one reverse-orientation site,
    with the forward one first. Every occurrence is counted rather than taking the outermost
    pair: an enzyme cuts at every site it finds, so an oligo with a third site does not
    release the fragment the design assumes, and the released core is checked for further
    sites for the same reason.
    """
    site, spacer, overhang = ENZYMES[enzyme]
    # Every occurrence, not just the outermost pair. `find`/`rfind` silently ignored extra
    # sites, so an oligo with a spurious GGTCTC inside it digested as though it were clean --
    # the enzyme would cut there too, and the released fragment would not be the intended one.
    forwards = _occurrences(oligo, site)
    reverses = _occurrences(oligo, rc(site))
    if len(forwards) > 1 or len(reverses) > 1:
        raise DigestError(f"{enzyme} cuts more than twice: {len(forwards)} forward and "
                          f"{len(reverses)} reverse sites at {forwards} / {reverses}; the "
                          f"released fragment is not the intended one")
    forward = forwards[0] if forwards else -1
    reverse = reverses[0] if reverses else -1
    if forward < 0 or reverse < 0:
        raise DigestError(f"{enzyme} site missing: forward={forward}, reverse={reverse}")
    if forward >= reverse:
        raise DigestError(f"{enzyme} sites are not in releasing orientation "
                          f"(forward at {forward}, reverse at {reverse})")

    top_cut = forward + len(site) + spacer
    # The reverse-strand site cuts the same distance from its own 3' end, measured leftwards.
    bottom_cut = reverse - spacer
    if bottom_cut - top_cut < overhang * 2:
        raise DigestError("released fragment is shorter than its own overhangs")
    fragment = oligo[top_cut:bottom_cut]
    for other in FORBIDDEN:
        pattern = ENZYMES[other][0]
        for strand, needle in (("forward", pattern), ("reverse", rc(pattern))):
            at = fragment.find(needle)
            if at >= 0:
                raise DigestError(f"released fragment carries a {other} site on the {strand} "
                                  f"strand at {at}; it would be cut again")
    return fragment[:overhang], fragment[overhang:-overhang], fragment[-overhang:]


def _occurrences(text: str, needle: str) -> list[int]:
    """Every start index, including overlapping ones."""
    found, at = [], text.find(needle)
    while at >= 0:
        found.append(at)
        at = text.find(needle, at + 1)
    return found


def ligate(pieces: list[tuple[str, str, str]], start_overhang: str) -> tuple[str, list[int]]:
    """Join pieces by complementary overhangs, each used exactly once.

    Walks the chain rather than trusting input order, so a design whose oligos happened to
    be listed in the wrong order would still assemble -- and one whose ends do not form a
    single chain fails here rather than silently producing a shorter product.
    """
    by_start: dict[str, list[int]] = {}
    for i, (five, _core, _three) in enumerate(pieces):
        by_start.setdefault(five, []).append(i)

    order, used, current = [], set(), start_overhang
    product = current
    while True:
        candidates = [i for i in by_start.get(current, []) if i not in used]
        if not candidates:
            break
        if len(candidates) > 1:
            raise DigestError(f"overhang {current} is presented by {len(candidates)} "
                              f"fragments; the assembly is ambiguous")
        index = candidates[0]
        used.add(index)
        order.append(index)
        _five, core, three = pieces[index]
        product += core + three
        current = three
    if len(used) != len(pieces):
        raise DigestError(f"only {len(used)} of {len(pieces)} fragments joined into one "
                          f"chain; ends do not match up")
    return product, order


def check(row: dict, destination: str = DEFAULT_DESTINATION) -> dict:
    """One design: digest its exported oligos, ligate, compare with the expected construct.

    The comparison is an equality against `five + cds + three`, where the flanking
    overhangs come from the declared destination rather than from the observed result.
    """
    oligos = [o["oligo_sequence_5to3"] for o in row["oligos"]]
    pieces = [digest(o) for o in oligos]
    starts = {p[0] for p in pieces}
    ends = {p[2] for p in pieces}
    # The chain begins at the overhang no fragment produces as its 3' end.
    entry = sorted(starts - ends)
    if len(entry) != 1:
        raise DigestError(f"expected exactly one entry overhang, found {entry}")
    product, order = ligate(pieces, entry[0])

    cds = row["cds"]
    five, three = DESTINATIONS[destination]
    expected = five + cds + three
    contains = cds in product
    return {
        "target": row["target"],
        "architecture": row["architecture"],
        "fragments": len(pieces),
        "assembly_order": order,
        "in_listed_order": order == sorted(order),
        "product_nt": len(product),
        "cds_nt": len(cds),
        "expected_nt": len(expected),
        "destination": destination,
        "cds_recovered": contains,
        # The gate. Containment was the old predicate and it passed two designs that were
        # demonstrably wrong; equality against the independently declared destination is what
        # actually states "this order yields this construct".
        "product_exact": product == expected,
        "entry_overhang_expected": five,
        "exit_overhang_expected": three,
        "leading_context": product[:product.find(cds)] if contains else None,
        "trailing_context": product[product.find(cds) + len(cds):] if contains else None,
        "junction_overhangs_seen": [p[2] for p in pieces[:-1]],
    }


def corrupt_controls(row: dict) -> list[tuple[str, str]]:
    """Deliberate damage that must be detected, with the reason each should fail."""
    base = [o["oligo_sequence_5to3"] for o in row["oligos"]]
    out = []

    one_base_lost = list(base)
    one_base_lost[1] = one_base_lost[1][:20] + one_base_lost[1][21:]
    out.append(("one base lost from fragment 1", one_base_lost))

    mismatched = list(base)
    site, spacer, oh = ENZYMES["BsaI"]
    cut = mismatched[2].find(site) + len(site) + spacer
    mismatched[2] = (mismatched[2][:cut] + "TTTT" + mismatched[2][cut + oh:])
    out.append(("junction overhang mutated on fragment 2", mismatched))

    flipped = list(base)
    f = flipped[0].find(site)
    flipped[0] = flipped[0][:f] + rc(site) + flipped[0][f + len(site):]
    out.append(("recognition site reversed on fragment 0", flipped))

    extra = list(base)
    extra[1] = extra[1][:60] + site + extra[1][60:]
    out.append(("extra recognition site inside fragment 1", extra))

    truncated = base[:-1]
    out.append(("terminal fragment missing", truncated))

    # The two corruptions the orchestrator found passing. Both leave the CDS intact inside
    # the product, so the old containment predicate accepted them; both are wrong assemblies.
    wrong_entry = list(base)
    entry_cut = wrong_entry[0].find(site) + len(site) + spacer
    wrong_entry[0] = (wrong_entry[0][:entry_cut] + "AAAA"
                      + wrong_entry[0][entry_cut + oh:])
    out.append(("entry overhang changed from the destination's", wrong_entry))

    site_inside = list(base)
    inside_cut = site_inside[0].find(site) + len(site) + spacer + oh
    site_inside[0] = (site_inside[0][:inside_cut] + site + site_inside[0][inside_cut:])
    out.append(("extra recognition site just inside the entry overhang", site_inside))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--corpus", default=str(ROOT.parent / "Users" / "syedz"))
    ap.add_argument("--jsonl", default=str(ROOT / "work" / "independent"
                                        / "corpus200_pinned.jsonl"))
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    rows = load_designs(args.jsonl)
    if args.limit:
        rows = rows[:args.limit]
    print(f"{len(rows)} exported designs; enzyme model BsaI GGTCTC(1/5), 4-base 5' overhang")

    ok, failures, by_arch, contexts = 0, [], {}, {}
    for row in rows:
        try:
            result = check(row)
        except (DigestError, KeyError) as exc:
            failures.append((row["target"], f"{type(exc).__name__}: {exc}"))
            continue
        if not result["product_exact"]:
            why = ("assembled product does not contain the CDS" if not result["cds_recovered"]
                   else f"product is not {result['destination']} + CDS + exit exactly "
                        f"(got {result['leading_context']!r} ... "
                        f"{result['trailing_context']!r}, expected "
                        f"{result['entry_overhang_expected']!r} ... "
                        f"{result['exit_overhang_expected']!r})")
            failures.append((row["target"], why))
            continue
        ok += 1
        by_arch[result["architecture"]] = by_arch.get(result["architecture"], 0) + 1
        contexts.setdefault((result["leading_context"], result["trailing_context"]), 0)
        contexts[(result["leading_context"], result["trailing_context"])] += 1

    print(f"\nproduct equals destination + CDS + exit exactly: {ok}/{len(rows)}   "
          f"by architecture {by_arch}")
    print("terminal context around the CDS (leading, trailing):")
    for (lead, trail), n in sorted(contexts.items(), key=lambda kv: -kv[1]):
        print(f"   {lead!r} ... {trail!r}   x{n}")
    for target, why in failures[:10]:
        print(f"   FAILED {target}: {why}")

    print("\ncorrupt controls on the first design (each must fail):")
    detected = 0
    for name, oligos in corrupt_controls(rows[0]):
        row = dict(rows[0], oligos=[{"oligo_sequence_5to3": o} for o in oligos])
        try:
            # The same acceptance predicate the real designs are gated on. Judging controls
            # by a weaker test than the corpus is how two wrong assemblies passed.
            result = check(row)
            verdict = ("NOT DETECTED" if result["product_exact"]
                       else "detected: product is not destination + CDS + exit exactly")
            detected += not result["product_exact"]
        except (DigestError, KeyError) as exc:
            verdict = f"detected: {type(exc).__name__}: {str(exc)[:56]}"
            detected += 1
        print(f"   {name:44s} {verdict}")

    out = ROOT / "work" / "independent" / "v5_digest_ligate.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(
        {"provenance": provenance(__file__, args.jsonl),
         "designs": len(rows), "reconstructed": ok, "by_architecture": by_arch,
         "failures": failures, "corrupt_controls_detected": detected,
         "corrupt_controls_total": len(corrupt_controls(rows[0])),
         "model": "BsaI GGTCTC(1/5), 4-base 5' overhang; double-stranded fragments only",
         "independence": ("reads only oligo_sequence_5to3 from exported records; imports no "
                          "clippr assembly code and does not read payload_5to3")},
        indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    # Compared against the number of controls there actually are, not a literal. Two
    # controls were added after this line was written and `detected == 5` silently became
    # unsatisfiable: the script reported every control detected in its JSON and exited 1
    # anyway, so the exit code said "failing" while the evidence said "clean".
    total_controls = len(corrupt_controls(rows[0]))
    return 0 if ok == len(rows) and detected == total_controls else 1


if __name__ == "__main__":
    raise SystemExit(main())
