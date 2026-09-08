"""Split a coding sequence into orderable fragments and wrap them for Golden Gate.

An ordered fragment is three pieces:

    [pad][enzyme site][spacer]  payload  [spacer][enzyme site, reverse-complemented][pad]

The payload already carries its own overhangs -- the 4 nucleotides a cut leaves are part
of the coding sequence on both sides of the junction (see `arelf`), so consecutive
payloads overlap by exactly that much and the arms do not add them.

**The 3' arm must be reverse-complemented.** Both Type IIS sites have to face *inward*, so
the enzyme cuts the payload out and leaves the arms behind. Writing the site in its
forward spelling at the 3' end points it the wrong way: the enzyme then cuts outward, into
the pad, and the fragment is destroyed. The failure is invisible in a sequence file and
shows up as an assembly that simply does not work, so `wrap_fragment` builds the 3' arm by
reverse-complementing a whole forward arm rather than by writing the site out, and
`test_assembly.py` asserts the forward site never appears at the 3' end.

Enzyme geometry is read from Biopython rather than tabulated: the spacer is
`fst5 - len(site)` and the overhang length is `abs(ovhg)`, so BsaI (1 nt spacer, 4 nt
overhang), BbsI (2, 4), BsmBI (1, 4) and SapI (1, 3) all follow without special cases.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, asdict

from . import constants as C

_COMPLEMENT = str.maketrans("ACGT", "TGCA")

#: Filler outside each enzyme site. A Type IIS enzyme binds poorly at a blunt end, so it
#: needs a few bases to sit on. These two 7-mers are arbitrary -- they were recovered from
#: the deposited designs so our oligos reproduce them byte for byte, and either may be
#: replaced with any sequence that introduces no new restriction site.
PAD_5: str = "ACAGCCA"
PAD_3: str = "TATCGGC"

DEFAULT_ENZYME: str = "BsaI"


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]


def enzyme_geometry(enzyme: str = DEFAULT_ENZYME) -> tuple[str, int, int]:
    """(recognition site, spacer length, overhang length) for a Type IIS enzyme."""
    from Bio import Restriction

    try:
        e = getattr(Restriction, enzyme)
    except AttributeError as exc:
        raise ValueError(f"unknown restriction enzyme {enzyme!r}") from exc
    site = str(e.site).upper()
    spacer = e.fst5 - len(site)
    overhang = abs(e.ovhg)
    if spacer < 0 or overhang <= 0:
        raise ValueError(f"{enzyme} does not cut outside its site; not a Type IIS enzyme")
    return site, spacer, overhang


@dataclass(frozen=True)
class Fragment:
    """One orderable piece of the coding sequence."""

    fragment_id: str
    assembly_order: int
    aa_start: int
    aa_end: int
    aa_length: int
    cds_start: int
    cds_end: int
    payload: str
    oh5_coding_site: str
    oh3_coding_site: str

    def as_dict(self) -> dict:
        return asdict(self)


def split_cds(
    cds: str,
    cuts: Sequence[int],
    overhangs: Sequence[str] | None = None,
    destination: tuple[str, str] | None = None,
    enzyme: str = DEFAULT_ENZYME,
) -> list[Fragment]:
    """Cut `cds` at the given residue positions into overlapping payloads.

    `cuts` are residue indices, so fragment boundaries land on codon boundaries. Passing
    `overhangs` cross-checks the junction sequences the caller believes it chose against
    what the coding sequence actually reads: a mismatch means codon optimisation moved a
    locked site, which would otherwise surface only as a failed assembly at the bench.

    `destination` is the vector's (5' coding site, 3' coding site) pair; it defaults to
    level 0. Those two are not part of the coding sequence -- they are contributed by the
    backbone -- so they are appended to the outer ends only.
    """
    cds = cds.upper().replace("U", "T")
    _, _, oh_len = enzyme_geometry(enzyme)
    dest5, dest3 = destination or C.DESTINATION_OVERHANGS["level0"]

    cuts = [int(c) for c in cuts]
    if cuts != sorted(set(cuts)):
        raise ValueError("cuts must be strictly ascending with no repeats")
    if len(cds) % 3:
        raise ValueError(f"coding sequence of {len(cds)} nt is not a whole number of codons")
    for c in cuts:
        if 3 * c - oh_len < 0 or 3 * c >= len(cds):
            raise ValueError(f"cut at residue {c} leaves no room for a {oh_len}-nt overhang")

    if overhangs is not None:
        if len(overhangs) != len(cuts):
            raise ValueError(f"{len(overhangs)} overhangs for {len(cuts)} cuts")
        for c, want in zip(cuts, overhangs):
            got = cds[3 * c - oh_len:3 * c]
            if got != str(want).upper().replace("U", "T"):
                raise ValueError(
                    f"cut at residue {c}: coding sequence reads {got}, but the chosen "
                    f"overhang is {want}. The locked site was not preserved."
                )

    starts = [0] + [3 * c for c in cuts]
    ends = [3 * c for c in cuts] + [len(cds)]

    fragments = []
    last = len(starts) - 1
    for k, (s, e) in enumerate(zip(starts, ends)):
        oh5 = dest5 if k == 0 else cds[s - oh_len:s]
        oh3 = dest3 if k == last else cds[e - oh_len:e]
        payload = oh5 + cds[s:e] + (dest3 if k == last else "")
        fragments.append(Fragment(
            fragment_id=f"F{k + 1}",
            assembly_order=k + 1,
            aa_start=s // 3,
            aa_end=e // 3,
            aa_length=(e - s) // 3,
            cds_start=s,
            cds_end=e,
            payload=payload,
            oh5_coding_site=oh5,
            oh3_coding_site=oh3,
        ))
    return fragments


def wrap_fragment(
    payload: str,
    enzyme: str = DEFAULT_ENZYME,
    five_pad: str = PAD_5,
    three_pad: str = PAD_3,
) -> str:
    """Add inward-facing Type IIS sites to a payload, giving the sequence to order.

    Both arms are built as `pad + site + spacer` and the 3' one is reverse-complemented
    whole. Constructing it that way rather than writing the site out at the 3' end is what
    makes the sites face inward; see the module docstring.
    """
    site, spacer, _ = enzyme_geometry(enzyme)
    five_arm = f"{five_pad}{site}{'A' * spacer}"
    three_arm = reverse_complement(f"{three_pad}{site}{'A' * spacer}")
    return f"{five_arm}{payload.upper()}{three_arm}"


def build_oligos(
    cds: str,
    cuts: Sequence[int],
    overhangs: Sequence[str] | None = None,
    destination: tuple[str, str] | None = None,
    enzyme: str = DEFAULT_ENZYME,
    prefix: str = "",
):
    """The ordered fragments for one assembly, as a DataFrame ready to export."""
    import pandas as pd

    fragments = split_cds(cds, cuts, overhangs, destination, enzyme)
    rows = []
    for f in fragments:
        oligo = wrap_fragment(f.payload, enzyme)
        rows.append({
            "order_fragment_id": f"{prefix}_{f.fragment_id}" if prefix else f.fragment_id,
            "fragment_id": f.fragment_id,
            "assembly_order": f.assembly_order,
            "aa_start_0based": f.aa_start,
            "aa_end_0based": f.aa_end,
            "aa_length": f.aa_length,
            "cds_start": f.cds_start,
            "cds_end": f.cds_end,
            "payload_5to3": f.payload,
            "oligo_sequence_5to3": oligo,
            "oligo_length": len(oligo),
            "oh5_coding_site_5to3": f.oh5_coding_site,
            "oh3_coding_site_5to3": f.oh3_coding_site,
            "five_prime_end_overhang": f.oh5_coding_site,
            "three_prime_end_overhang": reverse_complement(f.oh3_coding_site),
            "wrap_enzyme": enzyme,
        })
    return pd.DataFrame(rows)


def reassemble(fragments: Sequence[Fragment], enzyme: str = DEFAULT_ENZYME) -> str:
    """Join payloads back into the coding sequence, as the ligation would.

    Consecutive payloads share their junction overhang, so each is appended with that
    overlap removed. Round-tripping this against the input is the cheapest possible check
    that a split did not lose or duplicate bases.
    """
    _, _, oh_len = enzyme_geometry(enzyme)
    if not fragments:
        raise ValueError("no fragments to reassemble")
    first, last = fragments[0], fragments[-1]
    out = first.payload[len(first.oh5_coding_site):]
    for f in fragments[1:]:
        out += f.payload[oh_len:]
    return out[:-len(last.oh3_coding_site)]
