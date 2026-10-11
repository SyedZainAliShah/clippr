"""Where to cut a PPR coding sequence for Golden Gate assembly.

A synthetic PPR is tandem repeats of one 31-residue template, so its coding sequence is
highly repetitive and must be ordered as several fragments and re-assembled. Two facts
set the geometry, both established by measurement rather than assumption
(see `validation/compare_cuts.py`):

1. Cuts are codon-aligned. A cut at residue `c` splits the CDS at nucleotide `3c`, so
   every fragment translates on its own.
2. The 4-nt overhang a cut produces is `cds[3c-4 : 3c]` -- the third nucleotide of codon
   `c-2` followed by the whole of codon `c-1`. Nothing downstream of the cut contributes.

Consequence: the overhang is a function of two adjacent residues and their codon choice.
Positions where those two residues are codon-degenerate can be tuned for ligation
fidelity without changing the protein; positions where they are not are fixed. This
module finds the tunable positions from the genetic code alone.

The ARELF motif recurs once per repeat and is invariant across every design, which makes
it a stable coordinate system for describing cut positions -- `offsets_from_motif()`
reports each tunable cut as a residue offset from the preceding ARELF. It is an anchor,
not a restriction: cuts are permitted anywhere the degeneracy criterion is met.
"""
from __future__ import annotations

import heapq
from collections.abc import Sequence

from . import constants as C


def find_motifs(protein: str) -> list[int]:
    """Zero-based indices of every ARELF occurrence in `protein`."""
    motif = C.ARELF_MOTIF
    return [i for i in range(len(protein) - len(motif) + 1)
            if protein[i:i + len(motif)] == motif]


def cut_context(cds: str, cut_aa: int, offset: int = 0) -> str:
    """The 4-nt overhang a codon-aligned cut at residue `cut_aa` would produce.

    `offset` shifts the cut by that many codons, so `cut_context(cds, a, k)` is the
    overhang for a cut at residue `a + k`. Pass an ARELF start as `cut_aa` to read
    offsets off the motif.
    """
    c = cut_aa + offset
    start = 3 * c - 4
    if start < 0 or 3 * c > len(cds):
        raise ValueError(
            f"cut at residue {c} is outside the coding sequence "
            f"(needs 3*{c}-4 >= 0 and 3*{c} <= {len(cds)})"
        )
    return cds[start:3 * c].upper()


def achievable_overhangs(protein: str, cut_aa: int) -> list[str]:
    """Every 4-nt overhang reachable at this cut by synonymous substitution.

    The overhang is (third base of codon `cut_aa-2`) + (codon `cut_aa-1`), so the
    reachable set is the product of the wobble bases of the first residue and the full
    codons of the second.
    """
    if cut_aa < 2 or cut_aa > len(protein):
        raise ValueError(f"cut at residue {cut_aa} has no two preceding residues")
    a, b = protein[cut_aa - 2], protein[cut_aa - 1]
    try:
        cod_a, cod_b = C.SYNONYMOUS_CODONS[a], C.SYNONYMOUS_CODONS[b]
    except KeyError as e:
        raise ValueError(f"no codons for residue {e.args[0]!r} at cut {cut_aa}") from e
    wobble = sorted({c[2] for c in cod_a})
    return sorted({w + c for w in wobble for c in cod_b})


def safe_overhangs(
    protein: str,
    cut_aa: int,
    enzymes: Sequence[str] | str = C.DEFAULT_ENZYME_PROFILE,
) -> list[str]:
    """Overhangs for which at least one *local* synonymous realization avoids the active
    restriction sites. Final feasibility is checked on the completed coding sequence.

    A **feasibility filter, not an optimiser, and local rather than global**: it asks
    whether an assembly junction is realizable in its immediate codon context, and leaves
    the choice among realizable sequences to the codon optimiser. Handing the optimiser a
    combination with no synonymous realization makes it fail outright rather than search.

    Locking an overhang freezes the wobble base of codon `cut_aa-2` and the whole of codon
    `cut_aa-1`, which can make a site unavoidable. The case measured on this corpus:
    overhang ``AGAC`` fixes codon `cut_aa-1` to ``GAC`` (Asp) and forces codon `cut_aa-2`
    to end in ``A`` -- and when that residue is arginine, both options (``AGA``, ``CGA``)
    give ``...GAGAC``. If the *following* residue's codons all begin with C or G -- true
    of Ala, Gly, Pro, Gln, Glu, Asp, Val and His -- the result is ``GAGACC`` or
    ``GAGACG``: **BsaI or BsmBI on the reverse strand**, with no synonymous escape.

    Note the site is on the reverse strand. Both orientations are searched for exactly
    that reason; a forward-only scan would miss this entire class.

    **Scope.** The window is one codon either side of the locked region. That is large
    enough to identify the junction-associated site classes considered here; it is not a
    complete restriction-site analysis. It does not model interactions further out, or
    sites formed between two different junctions -- both need the whole sequence, which
    does not exist yet at this stage. This is a *necessary* condition for feasibility, not
    a sufficient one: the codon optimiser still enforces the real constraint across the
    finished sequence, and `qc.synthesis_qc` checks it independently afterwards.

    What absence of a recognition site means is also bounded: it is the modelled site that
    is absent, which is the right hard constraint when domestication is required. It is
    not a claim that no enzyme could ever cut, which would have to account for star
    activity, methylation and reaction conditions.

    `achievable_overhangs` deliberately stays unfiltered, because validation needs to ask
    what is reachable in principle. This is what a design should choose from.
    """
    from Bio import Restriction

    sites = []
    for name in C.enzymes_for(enzymes):
        site = str(getattr(Restriction, name).site).upper()
        sites.append(site)
        rc = site.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        if rc != site:
            sites.append(rc)

    before = C.SYNONYMOUS_CODONS.get(protein[cut_aa - 2], ())
    after = (C.SYNONYMOUS_CODONS.get(protein[cut_aa], ())
             if cut_aa < len(protein) else ("",))

    out = []
    for oh in achievable_overhangs(protein, cut_aa):
        # codon cut_aa-2 must end in oh[0]; codon cut_aa-1 is fixed to oh[1:]
        lefts = [c for c in before if c[2] == oh[0]]
        if not lefts:
            # No codon for this residue can supply that wobble base, so the overhang has
            # no realization at all. Substituting an empty left context here would turn
            # "impossible" into "unconstrained" and let it pass -- the exact inversion
            # this filter exists to prevent. Unreachable from achievable_overhangs, which
            # derives oh[0] from these same codons, but the contract is stated here.
            continue
        if any(not any(s in left + oh[1:] + right for s in sites)
               for left in lefts for right in after):
            out.append(oh)
    return out


def realization_count(
    protein: str,
    cut_aa: int,
    overhang: str,
    enzymes: Sequence[str] | str = C.DEFAULT_ENZYME_PROFILE,
) -> int:
    """How many local synonymous realizations of this overhang avoid every active site.

    Diagnostic only -- it is deliberately **not** used to rank overhangs. Two candidates
    can have near-identical predicted fidelity while one leaves 64 codon realizations and
    the other leaves 2, and that is worth reporting to whoever reads the design. Feeding
    it into selection would mix constraint satisfaction with preference ranking, so the
    number is recorded and left for a human to weigh.
    """
    from Bio import Restriction

    sites = []
    for name in C.enzymes_for(enzymes):
        site = str(getattr(Restriction, name).site).upper()
        sites.append(site)
        rc = site.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        if rc != site:
            sites.append(rc)

    oh = overhang.upper().replace("U", "T")
    before = C.SYNONYMOUS_CODONS.get(protein[cut_aa - 2], ())
    after = (C.SYNONYMOUS_CODONS.get(protein[cut_aa], ())
             if cut_aa < len(protein) else ("",))
    lefts = [c for c in before if c[2] == oh[0]]
    return sum(1 for left in lefts for right in after
               if not any(s in left + oh[1:] + right for s in sites))


def tunable_cuts(protein: str, min_options: int = 2) -> list[int]:
    """Residue positions where the overhang can be varied at least `min_options` ways.

    A cut with exactly one achievable overhang is a liability: it cannot be moved off a
    sequence that clashes with the rest of the reaction, so the default excludes those
    and nothing else. Raising the threshold discards usable positions -- the corpus
    designs cut at KM, which offers only two options, and those designs assemble -- so
    prefer to let `overhangs.best_set` weigh a narrow position against the reaction it
    sits in rather than filtering it out here.
    """
    return [c for c in range(2, len(protein) + 1)
            if len(achievable_overhangs(protein, c)) >= min_options]


def offsets_from_motif(protein: str, cuts: Sequence[int] | None = None) -> list[int]:
    """Distinct residue offsets of `cuts` from the ARELF motif that precedes each.

    Describes cut positions in repeat-relative coordinates, which is how they are
    reported in the literature and how they stay comparable across architectures.
    """
    cuts = tunable_cuts(protein) if cuts is None else cuts
    motifs = find_motifs(protein)
    out = set()
    for c in cuts:
        prior = [m for m in motifs if m <= c]
        if prior:
            out.add(c - max(prior))
    return sorted(out)


def candidate_cuts(
    protein: str,
    n_fragments: int,
    min_fragment_aa: int = 40,
    max_fragment_aa: int = 120,
    min_options: int = 2,
    max_results: int = 20000,
) -> list[list[int]]:
    """Every way to split `protein` into `n_fragments` tunable, length-legal fragments.

    Fragment bounds are in residues; the defaults keep every fragment inside the length
    an oligo-pool vendor will synthesise as one piece (120 aa = 360 nt plus overhangs
    and flanks). Returns cut positions only -- the fragments are the spans between them.

    Raises if more than `max_results` splits qualify. Returning the first `max_results`
    would look like a sample but is not one: combinations are generated in ascending cut
    order, so a truncated list is biased towards splits that cut early. Narrow the
    fragment bounds instead.
    """
    if n_fragments < 1:
        raise ValueError("n_fragments must be at least 1")
    n = len(protein)
    if n_fragments == 1:
        return [[]] if min_fragment_aa <= n <= max_fragment_aa else []
    if not n_fragments * min_fragment_aa <= n <= n_fragments * max_fragment_aa:
        raise ValueError(
            f"{n} residues cannot split into {n_fragments} fragments of "
            f"{min_fragment_aa}-{max_fragment_aa} aa"
        )

    allowed = set(tunable_cuts(protein, min_options))
    out: list[list[int]] = []

    def extend(prefix: list[int], last: int, remaining: int) -> None:
        """Grow a split left to right, never generating one that cannot be completed.

        Filtering after `itertools.combinations` is not viable: a 612-residue 19S protein
        split seven ways would iterate C(610, 6) -- about 5e13 -- before rejecting almost
        all of them. Pruning at each step keeps the walk inside the legal region.
        """
        if len(out) > max_results:
            raise ValueError(
                f"more than {max_results:,} splits satisfy {min_fragment_aa}-"
                f"{max_fragment_aa} aa per fragment; narrow the fragment bounds "
                f"(mean fragment here is {n / n_fragments:.0f} aa)")
        if remaining == 0:
            if min_fragment_aa <= n - last <= max_fragment_aa:
                out.append(list(prefix))
            return
        lo = last + min_fragment_aa
        # leave room for every fragment still to come, plus the final one
        hi = min(last + max_fragment_aa, n - remaining * min_fragment_aa)
        for cut in range(lo, hi + 1):
            if cut in allowed:
                prefix.append(cut)
                extend(prefix, cut, remaining - 1)
                prefix.pop()

    extend([], 0, n_fragments - 1)
    return out


def balanced_cuts(
    protein: str,
    n_fragments: int,
    min_fragment_aa: int = 40,
    max_fragment_aa: int = 120,
    window: int = 4,
    limit: int = 200,
    min_options: int = 2,
) -> list[list[int]]:
    """Legal splits nearest to even spacing, most balanced first.

    Enumerating every split is infeasible for the larger architectures and unnecessary:
    cost and synthesis risk both track the *longest* fragment, so the useful candidates
    all sit near even spacing. This searches outward from the ideal positions by
    increasing total displacement and stops at `limit`, which bounds the work regardless
    of protein length.
    """
    n = len(protein)
    if n_fragments < 2:
        return [[]]
    allowed = set(tunable_cuts(protein, min_options))
    ideal = [round((i + 1) * n / n_fragments) for i in range(n_fragments - 1)]

    per: list[list[int]] = []
    for target in ideal:
        opts = sorted((p for p in range(target - window, target + window + 1)
                       if p in allowed and 2 <= p < n),
                      key=lambda p: (abs(p - target), p))
        if not opts:
            return []
        per.append(opts)

    def cost(idx: tuple[int, ...]) -> int:
        return sum(abs(per[i][j] - ideal[i]) for i, j in enumerate(idx))

    start = tuple(0 for _ in per)
    heap = [(cost(start), start)]
    seen = {start}
    out: list[list[int]] = []
    while heap and len(out) < limit:
        _, idx = heapq.heappop(heap)
        cuts = [per[i][j] for i, j in enumerate(idx)]
        if cuts == sorted(set(cuts)):
            bounds = [0] + cuts + [n]
            if all(min_fragment_aa <= b - a <= max_fragment_aa
                   for a, b in zip(bounds, bounds[1:])):
                out.append(cuts)
        for i in range(len(per)):
            if idx[i] + 1 < len(per[i]):
                nxt = idx[:i] + (idx[i] + 1,) + idx[i + 1:]
                if nxt not in seen:
                    seen.add(nxt)
                    heapq.heappush(heap, (cost(nxt), nxt))
    return out
