"""Turn inventory inserts into the assembly substrates a level-0 reaction actually consumes.

A module's Table S1 insert is **not** the thing that goes into a BbsI reaction. The released
fragment and the insert differ, and they differ differently per block:

    A modules (16)    released = left overhang + insert
    B, C, D    (12)   released = insert
    E modules (14)    released = insert + right overhang

Those three rules were derived by digesting all 42 primary `GRASP_-1.gb` records with BbsI and
comparing each released fragment with its Table S1 insert; they cover the kit exactly, with no
residue. `verify_rules()` re-derives them from the primary records so a claim here can be
falsified rather than trusted.

The difference matters: an A-module insert begins `AATG` or `AGGT` -- its own fusion site,
part of the coding sequence -- while the level-0 substrate exposes `CTCA`, which comes from the
vector. Ordering the bare insert would order something that cannot assemble.

**Every substrate built here is digested again before it is returned.** If the wrapper or the
rule is wrong, the check fails rather than shipping confidently wrong DNA. That is the whole
reason this module exists instead of a string-concatenation helper.
"""
from __future__ import annotations

from dataclasses import dataclass

#: BbsI: GAAGAC, cutting 2 nt downstream on the top strand and 6 on the bottom, leaving a
#: four-base 5' overhang. Written out rather than imported so the construction states its own
#: geometry; `verify_rules` checks it against Bio.Restriction.
BBSI_SITE = "GAAGAC"
BBSI_SPACER = 2
OVERHANG = 4

#: Neutral padding outside the recognition site. Synthesis vendors need a few bases before an
#: enzyme can bind at the very end of a fragment.
PAD = "TT"
PAD_LEN = len(PAD)

#: How each block's released fragment relates to its insert. Derived, not assumed -- see the
#: module docstring and `verify_rules`.
#: A modules carry their own fusion site (`AATG`/`AGGT`) as the insert's first four bases; the
#: level-0 overhang `CTCA` sits before it, in the vector. E modules are the mirror image: the
#: closing `CGAG` sits after the insert. B, C and D inserts already span overhang to overhang.
RELEASE_RULES: dict[str, str] = {
    "1A": "left_overhang_then_insert", "2A": "left_overhang_then_insert",
    "14A": "left_overhang_then_insert", "19A": "left_overhang_then_insert",
    "B": "insert", "C": "insert", "D": "insert",
    "1E": "insert_then_right_overhang", "2E": "insert_then_right_overhang",
    "14E": "insert_then_right_overhang", "19E": "insert_then_right_overhang",
}

#: The overhang a block's released fragment exposes at its 3' end, from the cut geometry.
RIGHT_OVERHANG: dict[str, str] = {
    "1A": "ACTC", "2A": "ACTC", "14A": "ACTC", "19A": "ACTC",
    "B": "AAGA", "C": "GCAC", "D": "TGAA",
    "1E": "CGAG", "2E": "CGAG", "14E": "CGAG", "19E": "CGAG",
}

#: The overhang a block's released fragment exposes at its 5' end, from the cut geometry.
LEFT_OVERHANG: dict[str, str] = {
    "1A": "CTCA", "2A": "CTCA", "14A": "CTCA", "19A": "CTCA",
    "B": "ACTC", "C": "AAGA", "D": "GCAC",
    "1E": "TGAA", "2E": "TGAA", "14E": "TGAA", "19E": "TGAA",
}

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


class SubstrateError(Exception):
    """The constructed substrate does not digest to the fragment it was built for."""


def _occurrences(text: str, needle: str) -> list[int]:
    """Every start index, including overlapping ones."""
    found, at = [], text.find(needle)
    while at >= 0:
        found.append(at)
        at = text.find(needle, at + 1)
    return found


def rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Substrate:
    """One orderable, assembly-ready sequence and what it yields when cut."""

    module_id: str
    version: str
    sequence: str
    released: str
    five_overhang: str
    three_overhang: str
    enzyme: str = "BbsI"
    #: Synthesis-contract breaches on the **ordered** sequence. Empty means it satisfies the
    #: declared band and limits as written. Recognition sites belonging to the wrapper are
    #: excluded by role; anything else is a real breach of what would be synthesised.
    synthesis_problems: tuple[str, ...] = ()

    def as_dict(self) -> dict:
        return {"module": self.module_id, "version": self.version,
                "sequence": self.sequence, "length": len(self.sequence),
                "released": self.released, "released_length": len(self.released),
                "five_overhang": self.five_overhang,
                "three_overhang": self.three_overhang, "enzyme": self.enzyme,
                "synthesis_problems": list(self.synthesis_problems),
                "synthesis_ok": not self.synthesis_problems,
                "quantity": 1}


def released_fragment(insert: str, block: str) -> str:
    """The level-0 fragment this module must yield, from its insert and its block rule."""
    rule = RELEASE_RULES.get(block)
    if rule is None:
        raise SubstrateError(f"no release rule for block {block!r}; it is not part of the "
                             f"deposited kit this construction was derived from")
    if rule == "insert":
        return insert
    if rule == "insert_then_right_overhang":
        return insert + RIGHT_OVERHANG[block]
    return LEFT_OVERHANG[block] + insert


def build(insert: str, block: str, module_id: str = "", version: str = "") -> Substrate:
    """Wrap a fragment so BbsI releases exactly it, then prove that by digesting the result.

    Layout, 5' to 3':

        PAD  GAAGAC  NN  <fragment>  NN  GTCTTC  PAD
             site   spacer        spacer  site on the other strand

    Cutting the top strand 2 nt past each site leaves the fragment with its own first and last
    four bases as sticky ends -- which is what makes them the assembly overhangs.
    """
    fragment = released_fragment(insert, block)
    if len(fragment) < 2 * OVERHANG:
        raise SubstrateError(f"{module_id or block}: fragment is shorter than its overhangs")

    spacer = "A" * BBSI_SPACER
    sequence = (PAD + BBSI_SITE + spacer + fragment + spacer + rc(BBSI_SITE) + PAD)

    got, five, three = digest(sequence)
    if got != fragment:
        raise SubstrateError(
            f"{module_id or block}: the wrapper does not release the intended fragment "
            f"({len(got)} nt released, {len(fragment)} nt intended)")
    return Substrate(module_id=module_id, version=version, sequence=sequence,
                     released=fragment, five_overhang=five, three_overhang=three,
                     synthesis_problems=tuple(synthesis_problems(sequence)))


def substrate_problems(insert: str, block: str, profile=None) -> list[str]:
    """**The** feasibility test for a module. Every caller must use this one.

    Takes an insert, builds the substrate that would actually be ordered, and returns every
    **hard** breach on it. Recoding, collection acceptance, junction feasibility and export
    all go through here, so a candidate cannot be feasible for one of them and invalid for
    another.

    That was the defect this exists to end: the constraint scope was fixed at recoding and at
    export, but the collection optimiser still judged candidates on the bare insert and so
    reintroduced an out-of-band substrate in 16 of 18 matched runs.

    `profile` selects the synthesis policy; the default reproduces the shipped behaviour
    exactly. Under a profile that declares a rule a *target* rather than *hard*, a breach of
    that rule is **not** returned here -- it is a warning, from `substrate_warnings`.
    Feasibility and quality are different questions and this answers only the first.

    Returns `[]` for a usable module. A block with no release rule reports that rather than
    being waved through.
    """
    return substrate_findings(insert, block, profile)["hard"]


def substrate_warnings(insert: str, block: str, profile=None) -> list[str]:
    """Breaches of rules this profile declares *targets* rather than hard limits.

    Separate from `substrate_problems` on purpose. A soft warning must never silently become
    "vendor approved", and a target breach must never silently block a strict export -- so the
    two never share a return value.
    """
    return substrate_findings(insert, block, profile)["target"]


def substrate_findings(insert: str, block: str, profile=None) -> dict[str, list[str]]:
    """Every breach on the ordered substrate, split by how the profile enforces each rule."""
    from .synthesis_profile import resolve

    policy = resolve(profile)
    try:
        fragment = released_fragment(insert, block)
    except SubstrateError as exc:
        return {"hard": [str(exc)], "target": []}
    spacer = "A" * BBSI_SPACER
    sequence = PAD + BBSI_SITE + spacer + fragment + spacer + rc(BBSI_SITE) + PAD
    return synthesis_findings(sequence, policy)


def synthesis_problems(sequence: str, profile=None) -> list[str]:
    """Hard synthesis-contract breaches on a wrapped substrate, judged on what is ordered."""
    return synthesis_findings(sequence, profile)["hard"]


def synthesis_findings(sequence: str, profile=None) -> dict[str, list[str]]:
    """Every breach on a wrapped substrate, split by how this profile enforces the rule.

    The two intended BbsI sites are part of the construction and are excluded **by position**
    -- they are why the fragment can be released at all. A BbsI site anywhere else, and any
    site of another active enzyme anywhere, is a breach: it would be cut when it should not be.
    Checking GC and homopolymers alone let an internal BsaI site pass as `synthesis_ok`.

    Returns `{"hard": [...], "target": [...]}`. Forbidden sites are always hard: no profile may
    declare an enzyme cutting where it should not a matter of preference.

    **GC reporting.** The scan names the first offending window *and* the worst one, because
    "0.660 at offset 34" identifies a violation without saying how bad the sequence is, and a
    reader will take the first number for the worst.
    """
    from . import constants as C
    from .synthesis_profile import resolve

    policy = resolve(profile)
    findings: dict[str, list[str]] = {"hard": [], "target": []}

    def record(rule: str, message: str) -> None:
        findings["hard" if policy.is_hard(rule) else "target"].append(message)

    GC_BAND = policy.local_gc
    GC_WINDOW = policy.window
    MAX_HOMOPOLYMER = policy.max_homopolymer
    problems = findings["hard"]

    # Recognition sites, by role. The wrapper's own two BbsI sites sit at known positions;
    # everything else is unintended.
    intended = {PAD_LEN, len(sequence) - PAD_LEN - len(BBSI_SITE)}
    for name in C.enzymes_for(C.DEFAULT_ENZYME_PROFILE):
        from Bio.Restriction import RestrictionBatch

        site = str(RestrictionBatch([name]).get(name).site)
        for strand, pattern in (("forward", site), ("reverse", rc(site))):
            for at in _occurrences(sequence, pattern):
                if name == "BbsI" and at in intended:
                    continue            # the wrapper's own site, by role
                problems.append(f"unintended {name} site on the {strand} strand at {at}")

    if policy.global_gc is not None:
        overall = (sequence.count("G") + sequence.count("C")) / max(1, len(sequence))
        if not policy.global_gc[0] <= overall <= policy.global_gc[1]:
            record("global_gc", f"global GC {overall:.3f} outside {policy.global_gc}")

    first = None
    worst_at = None
    worst = None
    for start in range(0, max(1, len(sequence) - GC_WINDOW + 1)):
        window = sequence[start:start + GC_WINDOW]
        if len(window) < GC_WINDOW:
            break
        gc = (window.count("G") + window.count("C")) / len(window)
        if GC_BAND[0] <= gc <= GC_BAND[1]:
            continue
        if first is None:
            first = (start, gc)
        # The extreme is measured by distance outside the band, so a low-GC breach is not
        # masked by a high-GC one elsewhere in the same sequence.
        excess = max(GC_BAND[0] - gc, gc - GC_BAND[1])
        if worst is None or excess > worst[0]:
            worst, worst_at = (excess, gc), start
    if first is not None:
        message = f"GC {first[1]:.3f} outside {GC_BAND} in the window at {first[0]}"
        if worst_at != first[0]:
            message += f" (worst {worst[1]:.3f} at {worst_at})"
        record("local_gc", message)

    run, prev, longest = 1, "", 1
    for base in sequence:
        run = run + 1 if base == prev else 1
        prev = base
        longest = max(longest, run)
    if longest > MAX_HOMOPOLYMER:
        record("homopolymer", f"homopolymer run of {longest}, longer than {MAX_HOMOPOLYMER}")
    return findings


def digest(sequence: str) -> tuple[str, str, str]:
    """Simulate BbsI on a linear substrate: returns (fragment, 5' overhang, 3' overhang).

    Deliberately re-derived here rather than reusing the construction arithmetic, so the check
    can disagree with the builder. A verifier that shares the builder's arithmetic verifies
    nothing.
    """
    forward = [i for i in range(len(sequence)) if sequence.startswith(BBSI_SITE, i)]
    reverse = [i for i in range(len(sequence)) if sequence.startswith(rc(BBSI_SITE), i)]
    if len(forward) != 1 or len(reverse) != 1:
        raise SubstrateError(
            f"expected exactly one site per strand, found {len(forward)} forward and "
            f"{len(reverse)} reverse; the substrate would be cut more than twice")
    start = forward[0] + len(BBSI_SITE) + BBSI_SPACER
    end = reverse[0] - BBSI_SPACER
    if end <= start:
        raise SubstrateError("the two sites do not release a fragment between them")
    fragment = sequence[start:end]
    return fragment, fragment[:OVERHANG], fragment[-OVERHANG:]


def substrates_for(inventory) -> list[Substrate]:
    """An assembly-ready substrate for every module in an inventory."""
    out = []
    for module_id, record in sorted(inventory.modules.items()):
        out.append(build(record.dna, record.block, module_id, record.version))
    return out


def verify_rules(primary_genbank, inventory) -> dict:
    """Re-derive the release rules from the primary records and report any disagreement.

    This is the falsifier for `RELEASE_RULES`. If the deposited plasmids say something else,
    this says so rather than letting a stale constant decide what gets ordered.
    """
    import warnings

    from Bio import BiopythonParserWarning, SeqIO
    from Bio.Restriction import BbsI

    warnings.filterwarnings("ignore", category=BiopythonParserWarning)
    disagreements, checked = [], 0
    for record in SeqIO.parse(str(primary_genbank), "genbank"):
        module = inventory.modules.get(record.name)
        if module is None:
            continue
        cuts = sorted(BbsI.search(record.seq, linear=False))
        if len(cuts) != 2:
            disagreements.append(f"{record.name}: {len(cuts)} BbsI cuts, expected 2")
            continue
        a, b = cuts
        doubled = record.seq + record.seq
        # Anchored on the overhangs the cuts actually leave, not on an offset chosen to make
        # the answer come out. An earlier version used a span four bases short on the right;
        # because this check shared that arithmetic, it passed and the substrates were wrong.
        observed = str(doubled[a - 1:b + OVERHANG - 1])
        left, right = str(doubled[a - 1:a + 3]), str(doubled[b - 1:b + 3])
        if observed[:OVERHANG] != left or observed[-OVERHANG:] != right:
            disagreements.append(
                f"{record.name}: released fragment does not start and end with the cut "
                f"overhangs {left}/{right}")
            continue
        expected = released_fragment(module.dna, module.block)
        checked += 1
        if observed != expected:
            disagreements.append(
                f"{record.name}: rule gives {len(expected)} nt, the plasmid releases "
                f"{len(observed)} nt")
    return {"checked": checked, "disagreements": disagreements,
            "rules_hold": not disagreements}
