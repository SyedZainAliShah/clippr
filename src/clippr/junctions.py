"""Which internal junction overhangs a design may actually use, derived not assumed.

A junction overhang is not free real estate. It sits inside the coding sequence, shared by the
two modules it joins, so changing it changes codons in both — and only the four-base sequences
reachable by *synonymous* substitution are available at all. Everything here derives that set
from the architecture and the genetic code rather than proposing overhangs and hoping.

**Published kit rules versus CLIPPR extensions**, kept apart deliberately:

  *published*   the deposited kit's own junction overhangs, one per block transition, and the
                external destination ends. These are facts about the deposited vectors.
  *CLIPPR*      any *other* synonymous four-mer at the same position. The kit does not
                sanction these; this package derives them, and a design using one is a CLIPPR
                extension that no published protocol has validated.

**An assignment is global, not per-position.** The kit gives every B→C join the same overhang,
so a junction role is the unit of choice: changing B→C changes it everywhere B meets C, and
every affected module is rebuilt consistently. Choosing ends independently per position would
produce a set that cannot assemble.

**External destination ends stay fixed** unless a caller supplies a different destination.
Internal junctions are redesignable; the ends that meet the backbone are not, because they are
a property of the vector rather than of this design.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

from .assembly_spec import DESTINATION_OVERHANGS
from .overhangs import palindromic, reverse_complement

INTERFACE = 4


@dataclass(frozen=True)
class JunctionRole:
    """One block transition, and the overhang the deposited kit gives it."""

    left_block: str
    right_block: str
    deposited: str

    @property
    def name(self) -> str:
        return f"{self.left_block}->{self.right_block}"


@dataclass(frozen=True)
class JunctionSite:
    """Where a junction sits in one product, and the codons it occupies."""

    role: str
    start: int
    codon_start: int
    codon_end: int
    offset: int
    amino_acids: str
    deposited: str


def _synonyms(genetic_code: int = 1) -> dict[str, list[str]]:
    from Bio.Data import CodonTable

    table = CodonTable.unambiguous_dna_by_id[genetic_code].forward_table
    out: dict[str, list[str]] = {}
    for codon, aa in table.items():
        out.setdefault(aa, []).append(codon)
    return {aa: sorted(codons) for aa, codons in out.items()}


def site_for(product: str, frame: int, start: int, deposited: str,
             role: str, genetic_code: int = 1) -> JunctionSite:
    """The codon context an overhang occupies, in the product's own reading frame."""
    from Bio.Data import CodonTable

    forward = CodonTable.unambiguous_dna_by_id[genetic_code].forward_table
    codon_start = frame + ((start - frame) // 3) * 3
    codon_end = frame + ((start + INTERFACE - 1 - frame) // 3) * 3 + 3
    codons = [product[i:i + 3] for i in range(codon_start, codon_end, 3)]
    return JunctionSite(
        role=role, start=start, codon_start=codon_start, codon_end=codon_end,
        offset=start - codon_start,
        amino_acids="".join(forward.get(c, "*") for c in codons),
        deposited=deposited)


def permitted_overhangs(site: JunctionSite, genetic_code: int = 1,
                        forbidden: tuple[str, ...] = ()) -> dict[str, str]:
    """Every four-mer reachable at this site without changing the protein.

    Returns `{overhang: codon_string}` so a caller can rebuild the affected bases rather than
    re-deriving them. The deposited overhang is always present: it is reachable by definition,
    and an incumbent that dropped out of its own candidate set would be a bug.

    `forbidden` removes any four-mer that would create a recognition site *within the
    overhang itself*; sites spanning the wider context are caught later, on the assembled
    product, because they depend on neighbours this function cannot see.
    """
    synonyms = _synonyms(genetic_code)
    if "*" in site.amino_acids:
        return {site.deposited: ""}

    options = [synonyms.get(aa, []) for aa in site.amino_acids]
    if not all(options):
        return {site.deposited: ""}

    found: dict[str, str] = {}
    for combination in itertools.product(*options):
        codons = "".join(combination)
        overhang = codons[site.offset:site.offset + INTERFACE]
        if len(overhang) != INTERFACE:
            continue
        if any(pattern in overhang or reverse_complement(pattern) in overhang
               for pattern in forbidden):
            continue
        found.setdefault(overhang, codons)
    found.setdefault(site.deposited, found.get(site.deposited, ""))
    return found


def valid_assignment(overhangs, destination: str = "level0") -> tuple[bool, str | None]:
    """Whether a set of internal overhangs can assemble alongside the destination ends.

    Three ways a set fails, all from the declared assembly model rather than from taste:

      * a repeated overhang -- two junctions that anneal to each other, so the order of the
        product is not determined
      * a palindrome -- anneals to itself
      * an overhang whose reverse complement is also in the set -- the two ends are
        complementary, so they ligate to each other

    The destination's own ends are included in the check: an internal junction that clashes
    with the backbone is as broken as one that clashes with another junction.
    """
    internal = [str(o).upper() for o in overhangs]
    ends = list(DESTINATION_OVERHANGS[destination])
    full = internal + ends

    seen: set[str] = set()
    for overhang in full:
        if overhang in seen:
            return False, f"{overhang} appears twice"
        seen.add(overhang)
    for overhang in full:
        if palindromic(overhang):
            return False, f"{overhang} is palindromic and anneals to itself"
    for overhang in full:
        partner = reverse_complement(overhang)
        if partner != overhang and partner in seen:
            return False, f"{overhang} and {partner} are complementary"
    return True, None


def roles_in(compiled: dict, inventory) -> list[JunctionRole]:
    """The distinct block transitions a compiled product uses, with their kit overhangs."""
    roles: dict[str, JunctionRole] = {}
    for junction in compiled["junctions"]:
        left = inventory.modules[junction["left"]].block
        right = inventory.modules[junction["right"]].block
        role = JunctionRole(left, right, junction["overhang"])
        roles.setdefault(role.name, role)
    return [roles[name] for name in sorted(roles)]


def sites_in(compiled: dict, inventory, genetic_code: int = 1) -> list[JunctionSite]:
    """Every junction position in a product, with the codon context it occupies."""
    frame = compiled["coding_interval"][0]
    product = compiled["product"]
    out = []
    for junction in compiled["junctions"]:
        left = inventory.modules[junction["left"]].block
        right = inventory.modules[junction["right"]].block
        out.append(site_for(product, frame, junction["start"], junction["overhang"],
                            f"{left}->{right}", genetic_code))
    return out
