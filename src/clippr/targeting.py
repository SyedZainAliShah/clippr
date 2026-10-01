"""Choose *which* window of a supplied sequence a PPR should read.

Until now the caller named the target outright: nine, fourteen or nineteen bases, picked by
hand. That is the entire input, and it costs two things a wet-lab reviewer asked about.

**Freedom given away.** A 5'UTR is longer than the site a PPR reads, so most of the sequence a
caller supplies represents a choice they made silently. Measured on a 200 nt stretch of the
chloroplast `rbcL` coding span, with the whole `rbcL` feature declared intended, **99 of 192**
possible 9-mers carry no other exact occurrence in the annotation, and **187 of 187** 14-mers
carry none. A window chosen for that property is usually available; choosing by hand throws the
chance away.

That 99 is 98 if the declaration is narrowed to the searched interval rather than the whole
feature, and the single window between the two counts is `UUCCACCUG`, which occurs twice inside
`rbcL` and nowhere else. Both counts are right for what they declare, which is the reason intent
is declared here rather than assumed -- and why `TargetWindow.multiply_intended` names a locus
that absorbed more than one occurrence instead of quietly absorbing it. Both are pinned in
`tests/test_targeting.py`.

**Occurrence evidence that cannot be read.** Given a bare k-mer, a host scan cannot separate the
site the caller meant from a coincidence elsewhere, because nobody told it which site was meant.
This module makes that declarable, and the declaration is recorded with the result.

Two things make the accounting harder than it looks, and both have already produced a wrong
number in this project:

  - **A gene name is not a location.** The Chlamydomonas chloroplast annotation holds 109
    features under 100 distinct names; `psbA` appears twice, byte-identical, at two genomic
    locations. An occurrence keyed on the name merges the two silently. Every count here is
    therefore keyed on feature *identity*.
  - **A window can recur inside the caller's own sequence.** Excluding a whole interval hides
    that. So the window's own position is the only occurrence treated as itself, and a second
    copy anywhere -- including inside the supplied sequence -- stays visible.

**What this does not do.** It does not predict binding: an exact match is a sequence coincidence
until a bench result says otherwise, and nothing here measures affinity. It does not rank on GC
or on overhang feasibility -- window GC is RNA, while the band `synthesis_profile` governs is the
coding DNA of a different molecule, and overhang feasibility is a property of the design that
follows, not a proxy available here. It produces a plain sequence, so the existing entry points
consume it unchanged.
"""
from __future__ import annotations

import hashlib
from collections.abc import Sequence
from dataclasses import dataclass, field

from .biology import ARCHITECTURES
# `_normalise` and `_context` are reused rather than reimplemented on purpose: if this module
# decided for itself what counts as a match, it could disagree with `offtarget` about whether a
# sequence occurs, and two modules disagreeing on that would be invisible in both.
from .offtarget import Transcript, _context, _normalise

#: Where a supplied sequence came from. This is never inferred: see `SearchSequence.origin`.
ORIGINS = ("designed", "native", "unknown")

#: How an offset in this module is to be read. Recorded on every result, because a
#: source-relative offset is **not** a genomic coordinate -- `psbA`'s feature spans five joined
#: segments, so `feature.start + offset` is wrong across a join and wrong again on the minus
#: strand. No conversion to genomic coordinates is offered here, deliberately.
COORDINATE_CONVENTION = "0-based offset from the start of the supplied sequence, read 5'->3'"


def _rna(seq: str) -> str:
    """Upper-case RNA. Downstream entry points take RNA, so windows are produced in it."""
    return str(seq).strip().upper().replace("T", "U")


@dataclass(frozen=True)
class Occurrence:
    """One exact occurrence of a window, keyed by what it was found *in*.

    `source` is a stable identity and `label` is for display. They are separate because a gene
    name is not an identity -- nine names in the shipped annotation cover two features each.
    """

    source: str
    label: str
    kind: str
    offset: int
    context: str = ""

    def __str__(self) -> str:
        return f"{self.label} ({self.kind}) at {self.offset}"


@dataclass(frozen=True)
class IntendedLocus:
    """A place the caller declares the target is *supposed* to occur.

    Declared as a locus rather than as a fixed occurrence because the intended offset moves with
    the window: sliding a 9-mer along `rbcL` moves its own position in `rbcL` with it, so a
    declaration pinned to one offset would be wrong for every window but one.

    `start`/`end` bound the declaration within that source. Absorbing a whole source is allowed
    (the default) but it is the blunt option: two copies of a gene are two loci, declared
    separately, so a result records which of them the caller meant.
    """

    source: str
    start: int = 0
    end: int | None = None

    def covers(self, occurrence: Occurrence) -> bool:
        if occurrence.source != self.source:
            return False
        if occurrence.offset < self.start:
            return False
        return self.end is None or occurrence.offset < self.end


@dataclass(frozen=True)
class SearchSequence:
    """The stretch a window is chosen from, with enough provenance to audit the choice.

    `origin` is required and has no default. Where the sequence came from decides how its
    occurrences should be read: a designed UTR's intended site sits in a plasmid and is absent
    from a wild-type reference, while a native transcript's intended site is *in* the reference
    and any scan will find it. Defaulting this would make a result's meaning depend on an
    assumption nobody recorded, so "unknown" must be chosen explicitly too.
    """

    identifier: str
    sequence: str
    origin: str

    def __post_init__(self) -> None:
        if self.origin not in ORIGINS:
            raise ValueError(
                f"origin must be one of {ORIGINS}, got {self.origin!r}. It is not inferred: "
                f"a designed sequence's intended site is absent from a wild-type reference "
                f"and a native one's is present, which changes how every occurrence reads.")
        _normalise(self.sequence)       # refuses anything that is not plain nucleotides

    @property
    def rna(self) -> str:
        return _rna(self.sequence)

    @property
    def content_sha256(self) -> str:
        """Hash of the normalised sequence, so a result can name what it was computed from."""
        return hashlib.sha256(self.rna.encode()).hexdigest()


@dataclass(frozen=True)
class TargetWindow:
    """One candidate window: what a PPR would read, where it came from, and what else carries it.

    `occurrences` is every exact occurrence found, never a filtered list. Classification is
    reported through `self_occurrence`, `intended` and `other_occurrences`, so a caller can
    always recount from the raw evidence.
    """

    sequence: str
    offset: int
    source_id: str
    source_sha256: str
    occurrences: tuple[Occurrence, ...] = ()
    intended_loci: tuple[IntendedLocus, ...] = ()
    coordinate_convention: str = COORDINATE_CONVENTION

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def intended_declared(self) -> bool:
        return bool(self.intended_loci)

    @property
    def intended(self) -> tuple[Occurrence, ...]:
        """The occurrences a declared locus covers -- all of them, not one per locus."""
        return tuple(o for o in self.occurrences
                     if any(loc.covers(o) for loc in self.intended_loci))

    @property
    def multiply_intended(self) -> tuple[str, ...]:
        """Loci that absorbed more than one occurrence of this window.

        A locus covering two copies is the case Codex's review warned about: excluding an
        interval can hide a second occurrence inside it. Naming the locus keeps it visible
        instead, so a caller can narrow the declaration if the repeat was not expected.
        """
        out = []
        for loc in self.intended_loci:
            n = sum(1 for o in self.occurrences if loc.covers(o))
            if n > 1:
                out.append(f"{loc.source} covers {n} occurrences of this window")
        return tuple(out)

    @property
    def self_occurrence(self) -> Occurrence | None:
        """This window at its own position. The one occurrence that is the window itself."""
        for o in self.occurrences:
            if o.source == f"search:{self.source_id}" and o.offset == self.offset:
                return o
        return None

    @property
    def other_occurrences(self) -> tuple[Occurrence, ...]:
        """Everything the window also matches, apart from itself and declared intent.

        A repeat inside the supplied sequence counts here. It is a real finding: the caller's
        own construct carrying the site twice is not the same experiment as carrying it once.
        """
        skip = {(o.source, o.offset) for o in self.intended}
        me = self.self_occurrence
        if me is not None:
            skip.add((me.source, me.offset))
        return tuple(o for o in self.occurrences if (o.source, o.offset) not in skip)


@dataclass(frozen=True)
class Selection:
    """What was chosen, what else qualified, and what did not -- with the reason in each case."""

    chosen: TargetWindow | None
    alternatives: tuple[TargetWindow, ...] = ()
    rejected: tuple[tuple[TargetWindow, str], ...] = ()
    reference: str = ""
    lengths: tuple[int, ...] = ()
    intended_declared: bool = False
    tie_count: int = 0
    reason: str | None = None
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def found(self) -> bool:
        return self.chosen is not None


def feature_identity(index: int, transcript: Transcript) -> str:
    """A stable id for an annotated feature. The index is what the name cannot supply."""
    return f"feature:{index}:{transcript.name}"


def occurrences_in_reference(target: str, transcripts: Sequence[Transcript]
                             ) -> tuple[Occurrence, ...]:
    """Every exact sense-strand occurrence in the annotation, each keyed by feature identity.

    Deliberately per-feature rather than `offtarget.find_in_transcripts` over the whole list:
    that function reports a gene *name*, and two features sharing a name are then
    indistinguishable. A test pins that this returns the same total count, so the two agree on
    how many and differ only in what they can tell apart.
    """
    probe = _normalise(target)
    out: list[Occurrence] = []
    for index, tr in enumerate(transcripts):
        i = tr.sequence.find(probe)
        while i != -1:
            out.append(Occurrence(source=feature_identity(index, tr), label=tr.name,
                                  kind=tr.kind, offset=i,
                                  context=_context(tr.sequence, i, len(probe))))
            i = tr.sequence.find(probe, i + 1)
    return tuple(out)


def occurrences_in_source(target: str, source: SearchSequence) -> tuple[Occurrence, ...]:
    """Every occurrence inside the supplied sequence, including repeats."""
    probe, hay = _normalise(target), _normalise(source.sequence)
    out: list[Occurrence] = []
    i = hay.find(probe)
    while i != -1:
        out.append(Occurrence(source=f"search:{source.identifier}", label=source.identifier,
                              kind="search_sequence", offset=i,
                              context=_context(hay, i, len(probe))))
        i = hay.find(probe, i + 1)
    return tuple(out)


def enumerate_windows(source: SearchSequence, length: int) -> tuple[tuple[int, str], ...]:
    """Every (offset, window) of `length` in the supplied sequence, left to right."""
    if length not in set(ARCHITECTURES.values()):
        raise ValueError(
            f"window length must be one of {sorted(ARCHITECTURES.values())}, got {length}. "
            f"A PPR needs one repeat per base, so a length with no architecture cannot be "
            f"built -- see `biology.ARCHITECTURES`.")
    rna = source.rna
    if len(rna) < length:
        return ()
    return tuple((i, rna[i:i + length]) for i in range(len(rna) - length + 1))


def screen_window(offset: int, window: str, source: SearchSequence,
                  transcripts: Sequence[Transcript] | None = None,
                  intended_loci: Sequence[IntendedLocus] = ()) -> TargetWindow:
    """Build one window's full occurrence record. `transcripts=None` means source-only."""
    found = occurrences_in_source(window, source)
    if transcripts is not None:
        found = found + occurrences_in_reference(window, transcripts)
    return TargetWindow(sequence=_rna(window), offset=offset, source_id=source.identifier,
                        source_sha256=source.content_sha256, occurrences=found,
                        intended_loci=tuple(intended_loci))


def loci_for_gene(name: str, transcripts: Sequence[Transcript]) -> tuple[IntendedLocus, ...]:
    """Every annotated feature carrying `name`, as separately declared loci.

    One locus per feature, never one per name. `psbA` returns **two**, because the annotation
    holds two byte-identical copies at different genomic locations -- a caller who means both
    says so by passing both, and a caller who means one passes one. Returning a single merged
    locus would make that distinction unexpressible, which is how the figure behind this module
    was wrong the first time.
    """
    return tuple(IntendedLocus(source=feature_identity(i, tr))
                 for i, tr in enumerate(transcripts) if tr.name == name)


def select_window(source: SearchSequence, lengths: Sequence[int],
                  transcripts: Sequence[Transcript] | None = None,
                  intended_loci: Sequence[IntendedLocus] = (),
                  max_other_occurrences: int = 0,
                  reference: str | None = None,
                  keep_alternatives: int = 10) -> Selection:
    """Choose the window with the least other-occurrence evidence against it.

    `lengths` is required and explicit. Asking for several is allowed and is recorded, because
    "a 9-mer is not clean here, a 14-mer is" is an answer worth getting -- but it is never
    assumed: a caller who wants one architecture passes one length.

    **Ranking, and what each term is.** Candidates are ordered by (1) fewest other occurrences,
    which is the only term carrying evidence; then (2) shortest window, because fewer bases means
    fewer modules to order; then (3) lowest offset, which is arbitrary and exists only so the
    result is deterministic. Terms 2 and 3 are **not** biological preferences and must not be
    reported as though a shorter or earlier window binds better. `tie_count` says how many
    candidates shared the winning evidence, so a caller can see when the choice came down to
    the arbitrary term.

    Returns a `Selection` with `chosen=None` and a stated reason when nothing qualifies. That is
    an answer, not an error: a repetitive stretch may genuinely contain no usable window.
    """
    if not lengths:
        raise ValueError("lengths is required: pass the architecture length(s) to consider, "
                         f"any of {sorted(ARCHITECTURES.values())}")
    wanted = tuple(dict.fromkeys(int(n) for n in lengths))

    loci = tuple(intended_loci)
    candidates: list[TargetWindow] = []
    for length in wanted:
        for offset, window in enumerate_windows(source, length):
            candidates.append(screen_window(offset, window, source, transcripts, loci))

    # Never name a reference that was not screened. `transcripts=None` screens the supplied
    # sequence and nothing else, so it reports exactly that -- an earlier version of this line
    # returned the chloroplast accession in that case, which would have put a reference a result
    # never consulted into its own provenance.
    if transcripts is None:
        ref = "none -- supplied sequence only"
    else:
        ref = reference if reference is not None else "supplied annotation"
    notes: list[str] = []
    if transcripts is None:
        notes.append("no reference annotation was screened -- occurrences are from the "
                     "supplied sequence only")
    if not loci:
        notes.append("no intended locus was declared, which is not the same as none existing; "
                     "for a native source every occurrence of the target in its own gene will "
                     "therefore count against it")

    if not candidates:
        return Selection(chosen=None, reference=ref, lengths=wanted,
                         intended_declared=bool(loci), notes=tuple(notes),
                         reason=f"the supplied sequence is {len(source.rna)} nt, shorter than "
                                f"every requested window length {wanted}")

    ok, rejected = [], []
    for w in candidates:
        n = len(w.other_occurrences)
        if n <= max_other_occurrences:
            ok.append(w)
        else:
            where = ", ".join(str(o) for o in w.other_occurrences[:3])
            rejected.append((w, f"{n} other exact occurrence(s): {where}"
                                f"{' ...' if n > 3 else ''}"))

    if not ok:
        return Selection(chosen=None, alternatives=(), rejected=tuple(rejected), reference=ref,
                         lengths=wanted, intended_declared=bool(loci), notes=tuple(notes),
                         reason=f"no window of {wanted} in {source.identifier} has "
                                f"{max_other_occurrences} or fewer other exact occurrences; "
                                f"{len(candidates)} were screened")

    ok.sort(key=lambda w: (len(w.other_occurrences), w.length, w.offset))
    best = (len(ok[0].other_occurrences), ok[0].length)
    ties = sum(1 for w in ok if (len(w.other_occurrences), w.length) == best)
    if ties > 1:
        notes.append(f"{ties} candidates tied on evidence and length; the lowest offset was "
                     f"taken, which is arbitrary")

    notes.extend(ok[0].multiply_intended)
    return Selection(chosen=ok[0], alternatives=tuple(ok[1:1 + keep_alternatives]),
                     rejected=tuple(rejected), reference=ref, lengths=wanted,
                     intended_declared=bool(loci), tie_count=ties, notes=tuple(notes))


def window_report(selection: Selection) -> str:
    """A readable account of the choice, carrying its own caveats."""
    lines: list[str] = []
    if not selection.found:
        lines.append(f"NO WINDOW SELECTED -- {selection.reason}")
    else:
        w = selection.chosen
        lines.append(f"selected {w.sequence}  ({w.length} nt, offset {w.offset} in "
                     f"{w.source_id})")
        lines.append(f"  other exact occurrences: {len(w.other_occurrences)}")
        for o in w.other_occurrences[:5]:
            lines.append(f"    {o}")
        if w.intended_declared:
            lines.append(f"  declared intended: {len(w.intended)}")
        lines.append(f"  {len(selection.alternatives)} alternative(s) scored as well or worse; "
                     f"{len(selection.rejected)} rejected")
    lines.append(f"  reference screened: {selection.reference}")
    lines.append(f"  lengths considered: {selection.lengths}")
    lines.append(f"  coordinates: {COORDINATE_CONVENTION}")
    for n in selection.notes:
        lines.append(f"  note: {n}")
    lines.append("An exact match is a sequence occurrence, not demonstrated binding. Nothing "
                 "here predicts affinity.")
    return "\n".join(lines)
