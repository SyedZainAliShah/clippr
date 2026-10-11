"""Build a PPR from the deposited GRASP module kit instead of ordering new DNA.

**Why a second realisation route exists.** Everything else in this package designs DNA to be
synthesised. But the GRASP authors deposited a 42-plasmid kit, and a lab that holds it can
assemble many PPRs from parts it already owns. For such a lab, "order 906 nt of new DNA" is
the wrong answer to a question that has a cheaper one. So a target gets two routes, judged by
the same audit:

    target ─┬─ de novo synthesis   maximum sequence freedom, costs DNA
            └─ GRASP module kit    no synthesis, constrained to the deposited parts

**Provenance, stated plainly.** The *idea* of compiling a target into an ordered part list is
the reference implementation's, and we found it by reading that implementation. The *data* here
is not: module identity, overhangs and plate positions are derived from the paper's published
Supplementary Table S1 by `tools/derive_parts.py`, and `validation/compare_parts.py` checks our
selections against the compositions published in Table S2. No reference source is incorporated.
This should be described as a GRASP-compatible route informed by the published implementation,
never as independently invented.

**The inventory turns out to be exactly sized for its job**, which is derivable from the
overhangs alone and worth stating because it explains every design decision below. The modules
chain through a graph with a single branch point:

    AATG ┐                                     ┌─ CTTC ── 2A ─┐
    AGGT ┴─ 1A ─ ACTC ─ B ─ AAGA ─ C ─ GCAC ─ D ─ TGAA ─┼─ GTGA ── 14A ─┼─ ACTC (loop)
                                                        ├─ CACG ── 19A ─┘
                                                        └─ TTCG ── 2E (end)

One loop of `B C D` plus a linker pair contributes five modules, and *n* bases need *n+1*
modules, so 9, 14 and 19 bases need 10, 15 and 20 modules -- two, three and four sub-assemblies.
Since Golden Gate needs unique overhangs within one reaction, each sub-assembly may use `B C D`
only once, and each internal join needs its own linker pair. Nineteen bases need three internal
linkers, and the kit contains exactly three: `1E/2A`, `14E/14A`, `19E/19A`. The kit is built for
19S and nothing longer.

**How a module encodes a base.** A module named `L<last>5<fifth>` carries the 5th residue for
its own position and the last residue for the position before it, so the base at position *k* is
read by the pair (fifth of module *k*, last of module *k+1*). That was derived from Table S2 and
reproduces 28 of its 31 published variants; the two that differ are inconsistencies within the
published table itself, not in this mapping, and are named in `validation/compare_parts.py`.

**What this route does not do.** It selects the PPR modules. It does not model the acceptor
plasmids, the non-PPR elements of the transcriptional unit, or the assembly protocol beyond
naming the hierarchy. Describe the output as a *GRASP-compatible PPR module realisation*, never
as a complete construct.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

PARTS_JSON = Path(__file__).resolve().parent / "parts.json"

#: Blocks that carry a base-reading pair, in the order they chain.
LOOP_BLOCKS = ("B", "C", "D")

#: Linker pairs, in the order the kit intends them to be used. Each internal join between
#: sub-assemblies consumes one, and there are exactly three, which caps the kit at 19S.
LINKERS = (("1E", "2A"), ("14E", "14A"), ("19E", "19A"))

#: Bases per sub-assembly: one `B C D` loop plus its linker.
BASES_PER_BLOCK_RUN = 5

#: The shortest architecture, 9S, is two sub-assemblies joined by one linker. A single run
#: (1A B C D 2E) is arithmetically well-formed but is not a PPR the kit builds, so the lower
#: bound has to be stated -- divisibility alone let 4-base targets through.
MIN_RUNS = 2

#: Every target length the kit can realise, derived rather than listed: N linkers allow N+1
#: runs, and a run of `BASES_PER_BLOCK_RUN` modules reads one base fewer than it has modules.
BUILDABLE_LENGTHS = tuple(runs * BASES_PER_BLOCK_RUN - 1
                          for runs in range(MIN_RUNS, len(LINKERS) + 2))


@dataclass(frozen=True)
class Module:
    """One deposited plasmid from the kit."""

    plasmid_id: str
    block: str
    last: str | None
    fifth: str | None
    five_overhang: str
    three_overhang: str
    plate: str
    vector: str
    insert_length: int

    def __str__(self) -> str:
        pair = f"{self.fifth or '-'}/{self.last or '-'}"
        return (f"{self.plasmid_id:24s} plate {self.plate:>3}  {self.five_overhang}.."
                f"{self.three_overhang}  5th/last {pair}")


@lru_cache(maxsize=1)
def inventory() -> tuple[Module, ...]:
    """The 42 deposited modules, derived from the paper's Table S1."""
    data = json.loads(PARTS_JSON.read_text(encoding="utf-8"))
    return tuple(Module(plasmid_id=m["plasmid_id"], block=m["block"], last=m["last"],
                        fifth=m["fifth"], five_overhang=m["five_overhang"],
                        three_overhang=m["three_overhang"], plate=m["plate"],
                        vector=m["vector"], insert_length=m["insert_length"])
                 for m in data["modules"])


@lru_cache(maxsize=1)
def code_to_base() -> dict[str, str]:
    return json.loads(PARTS_JSON.read_text(encoding="utf-8"))["code_to_base"]


def base_to_code(base: str) -> str:
    """The (5th, last) residue pair that reads this base."""
    base = base.upper().replace("T", "U")
    for code, b in code_to_base().items():
        if b == base:
            return code
    raise ValueError(f"no PPR code reads {base!r}")


def block_chain(n_bases: int) -> list[str]:
    """The ordered blocks realising an n-base target, derived from the overhang graph.

    Raises for any length the kit cannot build, rather than returning a partial answer.
    """
    runs = (n_bases + 1) // BASES_PER_BLOCK_RUN
    if runs > len(LINKERS) + 1:
        raise ValueError(
            f"a {n_bases}-base target needs {runs - 1} internal linkers but the kit "
            f"contains only {len(LINKERS)} ({', '.join(a for a, _ in LINKERS)}), so the "
            f"deposited inventory cannot build it")
    if n_bases not in BUILDABLE_LENGTHS:
        listed = ", ".join(str(n) for n in BUILDABLE_LENGTHS[:-1])
        raise ValueError(
            f"the kit builds targets of {listed} or {BUILDABLE_LENGTHS[-1]} bases; "
            f"{n_bases} is not one of them")

    chain = ["1A"]
    for i in range(runs):
        chain += list(LOOP_BLOCKS)
        if i < runs - 1:
            chain += list(LINKERS[i])
    chain.append("2E")
    return chain


#: The two 5' fusion sites the deposited start modules offer. `AATG` carries the initiating
#: methionine itself; `AGGT` continues an N-terminal fusion from the backbone. Both exist for
#: every 1A variant, and the choice belongs to the design context.
FUSION_SITES = ("AATG", "AGGT")
DEFAULT_FUSION_SITE = "AATG"


def _find(block: str, fifth: str | None, last: str | None,
          fusion_site: str | None = None) -> Module | None:
    """The deposited module for this block and specificity pair.

    `fusion_site` disambiguates the start block, whose `AATG` and `AGGT` variants are
    identical in block, specificity residues and length and differ only in their 5' overhang.
    Without it the first listed variant always won, which is why the two `AGGT` modules were
    unreachable and went unexercised by all 200 corpus targets.
    """
    for m in inventory():
        if m.block != block:
            continue
        if fifth is not None and m.fifth != fifth:
            continue
        if last is not None and m.last != last:
            continue
        if fusion_site is not None and m.five_overhang in FUSION_SITES                 and m.five_overhang != fusion_site:
            continue
        return m
    return None


@dataclass(frozen=True)
class PartsPlan:
    """A target realised as deposited modules, or the reason it cannot be."""

    target: str
    modules: tuple[Module, ...] = ()
    sub_assemblies: tuple[tuple[Module, ...], ...] = ()
    reason: str | None = None

    @property
    def available(self) -> bool:
        """A route is available only when it actually produced modules.

        Defining this as `reason is None` alone let a default-constructed plan report as
        available and then crash `report` on `modules[0]`.
        """
        return self.reason is None and bool(self.modules)

    @property
    def plates(self) -> tuple[str, ...]:
        return tuple(m.plate for m in self.modules)


def select(target_rna: str, fusion_site: str = DEFAULT_FUSION_SITE) -> PartsPlan:
    """Compile a target RNA into deposited modules, or say why the kit cannot build it.

    Returns a plan rather than raising, because "this target needs synthesis" is an answer the
    caller should be able to act on, not an error.
    """
    target = str(target_rna).strip().upper().replace("T", "U")
    if not target or set(target) - set("ACGU"):
        return PartsPlan(target=target_rna, reason="not a plain RNA sequence")

    try:
        chain = block_chain(len(target))
    except ValueError as exc:
        return PartsPlan(target=target, reason=str(exc))

    codes = [base_to_code(b) for b in target]         # each is (5th, last)
    chosen: list[Module] = []
    for i, block in enumerate(chain):
        # module i carries the 5th residue for base i, and the last residue for base i-1
        fifth = codes[i][0] if i < len(codes) else None
        last = codes[i - 1][1] if i > 0 else None
        module = _find(block, fifth, last, fusion_site)
        if module is None:
            return PartsPlan(
                target=target,
                reason=(f"no deposited module for block {block} with 5th residue "
                        f"{fifth or '-'} and last residue {last or '-'}"))
        chosen.append(module)

    # split at each linker boundary: a sub-assembly may use B/C/D only once, because Golden
    # Gate needs unique overhangs within a single reaction
    subs: list[tuple[Module, ...]] = []
    current: list[Module] = []
    for module in chosen:
        current.append(module)
        if module.block in {a for a, _ in LINKERS}:
            subs.append(tuple(current))
            current = []
    if current:
        subs.append(tuple(current))

    return PartsPlan(target=target, modules=tuple(chosen), sub_assemblies=tuple(subs))


def report(plan: PartsPlan) -> str:
    """A readable pick list, or a clear statement that this route is unavailable."""
    if not plan.available:
        why = plan.reason or "no modules were selected"
        return (f"GRASP kit route unavailable for {plan.target}: {why}\n"
                f"Use the de novo synthesis route instead.")

    lines = [f"GRASP kit route for {plan.target} — {len(plan.modules)} modules "
             f"in {len(plan.sub_assemblies)} sub-assemblies",
             "", "no new PPR DNA needs synthesising; every module below is a deposited "
             "plasmid", ""]
    for i, sub in enumerate(plan.sub_assemblies, 1):
        lines.append(f"  sub-assembly {i}  ({sub[0].five_overhang} -> "
                     f"{sub[-1].three_overhang})")
        for m in sub:
            lines.append(f"    {m}")
        lines.append("")
    lines += [
        f"plate positions: {', '.join(plan.plates)}",
        f"vector: {plan.modules[0].vector} (all modules)",
        "",
        "This selects the PPR modules only. The acceptor plasmids, the non-PPR elements of "
        "the transcriptional unit and the assembly protocol are not modelled here — this is a "
        "GRASP-compatible PPR module realisation, not a complete construct.",
        "Module identity and overhangs derive from Dennis et al. 2025 Supplementary Table S1.",
    ]
    return "\n".join(lines)
