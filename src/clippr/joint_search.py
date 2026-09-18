"""Joint junction/codon search, with an observed Pareto set rather than a single answer.

Searches over **global junction assignments**: a choice of overhang for each junction class,
applied everywhere that class occurs. Choosing ends independently per position would produce a
set that cannot assemble, so the class is the unit of choice.

A *junction class* is a codon context plus the overhang the deposited kit puts there. The kit
gives every `1A->B`, `2A->B`, `14A->B` and `19A->B` join the same `ACTC` in the same `EL`
context at the same offset, so they are one physical interface and move together.

**Three declared objectives**, directions fixed before any search runs:

  ``fidelity``       predicted ligation fidelity -- **maximise**
  ``adaptation``     mean codon adaptation over the inventory -- **maximise**
  ``synthesis_fitness``  composite synthesis quality of the ordered substrates -- **maximise**
                   (replaced ``repeat_burden``, which was identically zero across this kit)

For a multi-reaction compilation, **every reaction's fidelity is reported and the minimum is
the search summary**. Not a product across reactions: a product reads as a predicted overall
yield, and it is not one -- see `docs/objectives.md` §5.

**Observed, not optimal.** The returned front is nondominated among the candidates actually
evaluated. Calling it the Pareto front would claim the feasible space was enumerated, and a
bounded beam search does not enumerate it.

**Candidates are evaluated as completed sequences**, not as overhang scores. A junction set
with excellent predicted fidelity whose codons cannot satisfy the synthesis constraints is
infeasible, and scoring it on fidelity alone would rank it highly on the strength of a design
that does not exist.
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field

from .assembly_spec import BLOCK_JOIN_OVERHANGS
from .inventories import INTERFACE, Inventory, compile_target, derive
from .junctions import permitted_overhangs, sites_in, valid_assignment
from .synthesis_profile import resolve as _resolve_profile

#: Objective directions, declared before searching. `docs/objectives.md` §6.
#: The three objectives, and their directions, declared before any search runs.
#:
#: `synthesis_fitness` replaced `repeat_burden` on 2026-09-17. The old axis was identically
#: zero over both the deposited and recoded inventories -- and not because of its k: it is
#: zero at k = 20, 16, 12 and 10, registering only at k = 8. A search reporting three
#: objectives was ordering by two. `repeat_burden` is still computed and reported as a
#: diagnostic, because a repetitive candidate is a real problem; it simply cannot rank the
#: candidates this kit produces.
DIRECTIONS = {"fidelity": "max", "adaptation": "max", "synthesis_fitness": "max"}

#: How far below the best observed feasible fidelity a candidate may sit and still be
#: considered for the default recommendation. A CLIPPR engineering choice, stated so the
#: recommendation is a rule rather than a preference applied after seeing results.
FIDELITY_TOLERANCE = 0.01


@dataclass(frozen=True)
class JunctionClass:
    """One physical interface, the roles that share it, and what may replace it."""

    stage: str
    key: str
    amino_acids: str
    offset: int
    deposited: str
    roles: tuple[str, ...]
    options: tuple[tuple[str, str], ...]        # (overhang, codons)

    @property
    def codon_span(self) -> int:
        return 3 * len(self.amino_acids)


@dataclass
class Candidate:
    """One assignment, its objective vector, and whether it is usable at all."""

    assignment: dict[str, str]
    feasible: bool
    reason: str | None = None
    objectives: dict | None = None
    per_reaction_fidelity: dict[str, list[dict]] = field(default_factory=dict)
    unscorable_reactions: list[str] = field(default_factory=list)
    inventory_version: str | None = None

    def as_dict(self) -> dict:
        return {"assignment": dict(sorted(self.assignment.items())),
                "feasible": self.feasible, "reason": self.reason,
                "objectives": self.objectives,
                "per_reaction_fidelity": self.per_reaction_fidelity,
                "unscorable_reactions": self.unscorable_reactions,
                "inventory_version": self.inventory_version}


def junction_classes(inv: Inventory, targets, genetic_code: int = 1,
                     forbidden: tuple[str, ...] = ("GGTCTC", "GAAGAC", "GCTCTTC"),
                     level0_only: bool = True) -> list[JunctionClass]:
    """Every distinct interface these targets use, and the overhangs each one permits.

    Grouped by codon context and deposited overhang rather than by role name, because roles
    sharing a context share a physical interface and must move together.
    """
    grouped: dict[str, dict] = {}
    for target in targets:
        compiled = compile_target(inv, target)
        if not compiled.get("available"):
            continue
        for site in sites_in(compiled, inv, genetic_code):
            key = f"{site.amino_acids}@{site.offset}:{site.deposited}"
            entry = grouped.setdefault(key, {"site": site, "roles": set()})
            entry["roles"].add(site.role)

    classes = []
    for key in sorted(grouped):
        site = grouped[key]["site"]
        options = permitted_overhangs(site, genetic_code, forbidden)
        classes.append(JunctionClass(
            stage=("level1" if site.deposited in BLOCK_JOIN_OVERHANGS else "level0"),
            key=key, amino_acids=site.amino_acids, offset=site.offset,
            deposited=site.deposited, roles=tuple(sorted(grouped[key]["roles"])),
            options=tuple(sorted(options.items()))))

    if level0_only:
        # A block-join class cannot change any reaction this package can score: its ends
        # belong to a level-1 reaction whose context is not supplied. Leaving them in the
        # search spends budget on changes whose effect is unmeasurable.
        classes = [c for c in classes if c.stage == "level0"]
    return classes


def apply_assignment(inv: Inventory, classes, assignment: dict[str, str],
                     label: str) -> Inventory:
    """Rebuild every module an assignment touches, on both sides of each junction.

    The codon span straddles the junction: the left module carries `offset + 4` of its bases
    and the right module carries the remainder. Rewriting only one side would leave the two
    modules unable to anneal, which is the failure this function exists to avoid.
    """
    replacements: dict[str, str] = {}
    by_key = {c.key: c for c in classes}

    # An assignment naming something this inventory does not have is a caller error, and it
    # was reaching here as a bare KeyError from a dict lookup -- a refusal, but one whose
    # message named only the missing string and not what it was missing from. The most likely
    # source is a front computed against a different inventory or a different target set, so
    # the error says which and lists what is actually available.
    unknown = sorted(set(assignment) - set(by_key))
    if unknown:
        raise ValueError(
            f"assignment names junction class(es) this inventory does not have: "
            f"{', '.join(unknown)}. Available: {', '.join(sorted(by_key)) or '(none)'}. "
            f"A front computed over different targets or a different inventory will do this.")

    for key, overhang in assignment.items():
        cls = by_key[key]
        if overhang == cls.deposited:
            continue
        options = dict(cls.options)
        if overhang not in options:
            raise ValueError(
                f"junction class {key} cannot take overhang {overhang!r}; its options are "
                f"{', '.join(sorted(options)) or '(none)'}. An overhang absent from the "
                f"options was never scored, so applying it would publish an unmeasured design.")
        codons = options[overhang]
        if not codons:
            continue
        left_tail = codons[:cls.offset + INTERFACE]
        right_head = codons[cls.offset:]

        for module_id, record in inv.modules.items():
            dna = replacements.get(module_id, record.dna)
            if record.three_interface == cls.deposited:
                dna = dna[:len(dna) - len(left_tail)] + left_tail
            if record.five_interface == cls.deposited:
                dna = right_head + dna[len(right_head):]
            if dna != replacements.get(module_id, record.dna):
                replacements[module_id] = dna

    return derive(inv, label, replacements, interface_assignment=label)


def evaluate(inv: Inventory, classes, assignment: dict[str, str], targets, table, *,
             k: int = 20, matrix: str = "BsaI-HFv2",
             destination: str = "level0", profile=None) -> Candidate:
    """Score one assignment on completed sequences, or say why it is infeasible."""
    from Bio.Seq import Seq

    from .objectives import codon_adaptation, duplicated_kmers
    from .overhangs import set_fidelity
    from .substrates import substrate_problems

    overhangs = sorted(set(assignment.values()))
    ok, why = valid_assignment(overhangs, destination)
    if not ok:
        return Candidate(assignment, False, f"invalid junction set: {why}")

    rebuilt = apply_assignment(inv, classes, assignment, label="candidate")

    # Proteins must survive. A synonymous overhang keeps the protein by construction, but
    # construction is an argument and this is a check.
    #
    # Constraint violations are judged as *introduced*, not inherited. One deposited module
    # already carries TTTTT and fails CLIPPR's homopolymer rule; treating that as infeasible
    # would make the deposited inventory unusable as a baseline, and the incumbent has to be
    # evaluable or there is nothing to compare against.
    for module_id, record in rebuilt.modules.items():
        original = inv.modules[module_id]
        s, e = original.coding_interval
        if len(record.dna) != len(original.dna):
            return Candidate(assignment, False, f"{module_id}: length changed")
        if str(Seq(record.dna[s:e]).translate()) != str(Seq(original.dna[s:e]).translate()):
            return Candidate(assignment, False, f"{module_id}: protein changed")
        # Only the *specific* problems the original already had are excused. The earlier
        # form -- `if problems and not module_constraint_problems(original.dna)` -- exempted
        # every problem in a module as soon as the original had any problem at all, so a
        # module that arrived with a homopolymer could acquire a brand-new BsaI site and
        # still pass.
        # Judged on the substrate that would be ordered, through the one shared validator.
        # Comparing insert-level problems accepted a candidate whose delivered substrate had
        # a 0.660 GC window.
        inherited = set(substrate_problems(original.dna, original.block, profile))
        introduced = [p for p in substrate_problems(record.dna, record.block, profile)
                      if p not in inherited]
        if introduced:
            return Candidate(assignment, False, f"{module_id}: {introduced[0]}")

    # Each reaction is scored with the enzyme that actually cuts it. A level-0 reaction is a
    # BbsI reaction -- all 42 deposited plasmids carry two BbsI sites and no BsaI site -- so
    # scoring it against BsaI-HFv2 describes a reaction that does not happen.
    #
    # A reaction is scored only when its full participant list is known. Level 0 qualifies.
    # Level 1 does not: its block release geometry is established, but the deposited BsaI
    # reaction co-assembles parts this package does not compile, and ligation fidelity is a
    # property of the whole competing overhang set. The PPR-only figure is carried as a named
    # diagnostic beside the reaction, never as its fidelity -- adding just the *known* omitted
    # ends to the 19S set moves it from 0.996046 to 0.742067, below level 0.
    #
    # Every scored reaction carries the overhang set it was scored on, so a number here can
    # always be traced back to the strands it came from rather than taken on trust.
    per_reaction: dict[str, list[dict]] = {}
    unscorable: list[str] = []
    for target in targets:
        compiled = compile_target(rebuilt, target)
        if not compiled.get("available"):
            return Candidate(assignment, False, f"{target}: no longer compiles")
        rows = []
        for reaction in compiled["stages"]["reactions"]:
            if not reaction["scorable"]:
                unscorable.append(f"{target} reaction {reaction['reaction']} "
                                  f"({reaction['stage']})")
                row = {"reaction": reaction["reaction"],
                       "stage": reaction["stage"], "enzyme": reaction["enzyme"],
                       "fidelity": None,
                       "geometry_established": reaction.get("geometry_established", False),
                       "participants_established": reaction.get("participants_established",
                                                                False),
                       "note": reaction["note"]}
                subset = reaction.get("block_subset_overhangs")
                if subset:
                    row["block_subset"] = {
                        "overhangs": list(subset),
                        "matrix": reaction["matrix"],
                        "fidelity": round(set_fidelity(subset, reaction["matrix"]), 6),
                        "is": ("a diagnostic over the PPR block junctions only, not this "
                               "reaction's fidelity: the competing set is incomplete")}
                rows.append(row)
                continue
            ends = reaction["reaction_overhangs"]
            valid, reason = valid_assignment(reaction["internal_overhangs"],
                                             reaction["stage"])
            if not valid:
                return Candidate(assignment, False,
                                 f"{target} reaction {reaction['reaction']}: {reason}")
            rows.append({"reaction": reaction["reaction"], "stage": reaction["stage"],
                         "enzyme": reaction["enzyme"],
                         "fidelity": round(set_fidelity(ends, reaction["matrix"]), 6),
                         "matrix": reaction["matrix"],
                         "overhangs": list(ends),
                         "ends_from": reaction["ends_from"],
                         "note": reaction["note"]})
        per_reaction[target] = rows

    cais = []
    for record in rebuilt.modules.values():
        s, e = record.coding_interval
        cais.append(codon_adaptation(record.dna[s:e], table)["cai"])

    # Repetition is measured on the sequences that would actually be **ordered**. For a
    # reusable inventory that is the distinct module inserts, not the assembled products: a
    # product's repeats are mostly the unavoidable consequence of reusing the same module,
    # and charging them to the module would penalise reuse for working.
    #
    # The product figure is kept as a separate named diagnostic rather than dropped.
    distinct = {record.dna for record in rebuilt.modules.values()}
    burden = sum(duplicated_kmers(dna, k=k) for dna in sorted(distinct))
    product_burden = sum(
        duplicated_kmers(compile_target(rebuilt, t)["product"], k=k) for t in targets)

    scored = [row["fidelity"] for rows in per_reaction.values() for row in rows
              if row["fidelity"] is not None]
    if not scored:
        return Candidate(assignment, False,
                         "no reaction in this compilation can be scored with the context "
                         "supplied")
    # The third axis, measured on the ordered substrates rather than the inserts. Its spread
    # is what makes it an axis: `repeat_burden` was identically zero across this kit, so
    # dominance was being decided by two objectives while three were reported.
    from .synthesis_fitness import inventory_fitness

    fitness = inventory_fitness(rebuilt, k=k, profile=profile)
    objectives = {
        # The minimum across *scorable* reactions, labelled. Never a product: a product reads
        # as a predicted overall yield and is not one.
        "fidelity": round(min(scored), 6),
        "adaptation": round(sum(cais) / len(cais), 6),
        "synthesis_fitness": fitness["mean"],
    }
    # Diagnostics: reported, never part of the dominance test. `repeat_burden` is retained
    # because it is the term that fires on a genuinely repetitive candidate.
    objectives["repeat_burden"] = burden
    objectives["product_repeat_burden_diagnostic"] = product_burden
    objectives["synthesis_fitness_spread"] = fitness["spread"]
    if fitness["degenerate"]:
        # Said out loud rather than discovered later. An axis with no spread is not an axis.
        objectives["synthesis_fitness_degenerate"] = True
    return Candidate(assignment, True, None, objectives, per_reaction,
                     sorted(set(unscorable)), rebuilt.version)


def search(inv: Inventory, targets, table, *, beam_width: int = 8, profile=None,
           max_evaluations: int = 20, wall_seconds: float = 120.0,
           k: int = 20, matrix: str = "BsaI-HFv2", destination: str = "level0",
           genetic_code: int = 1) -> dict:
    """Bounded beam search over junction assignments.

    Starts from the deposited assignment as incumbent and expands by changing one junction
    class at a time. This search genuinely uses no RNG -- candidates are enumerated from the
    permitted-overhang sets in sorted order -- so the same inputs and budget give the same
    result without a seed. That is not true of `library_search`, whose proposals are drawn
    from a seeded generator.
    """
    from .objectives import nondominated

    classes = junction_classes(inv, targets, genetic_code)
    incumbent = {c.key: c.deposited for c in classes}
    started = time.perf_counter()

    evaluated: dict[tuple, Candidate] = {}

    def score(assignment: dict[str, str]) -> Candidate:
        key = tuple(sorted(assignment.items()))
        if key not in evaluated:
            evaluated[key] = evaluate(inv, classes, assignment, targets, table,
                                      k=k, matrix=matrix, destination=destination,
                                      profile=profile)
        return evaluated[key]

    base = score(incumbent)
    beam = [incumbent]
    # A budget that stops the search before it expands anything at all is still an exhausted
    # budget. Setting the flag only inside the loop reported `complete` for a search that had
    # evaluated nothing but the incumbent -- which reads as "no improvement exists" when the
    # truth is "nothing was tried".
    budget_exhausted = len(evaluated) >= max_evaluations

    while len(evaluated) < max_evaluations:
        if time.perf_counter() - started >= wall_seconds:
            budget_exhausted = True
            break
        expansions = []
        for assignment in beam:
            for cls in classes:
                for overhang, _codons in cls.options:
                    if overhang == assignment[cls.key]:
                        continue
                    nxt = dict(assignment)
                    nxt[cls.key] = overhang
                    expansions.append(nxt)
        if not expansions:
            break

        fresh = []
        for assignment in expansions:
            if len(evaluated) >= max_evaluations:
                budget_exhausted = True
                break
            if time.perf_counter() - started >= wall_seconds:
                budget_exhausted = True
                break
            candidate = score(assignment)
            if candidate.feasible:
                fresh.append(assignment)
        if not fresh:
            break
        # Rank the beam by a declared scalar; the Pareto set is computed separately over
        # every evaluated candidate, so beam ranking never decides what is reported.
        fresh.sort(key=lambda a: (-score(a).objectives["fidelity"],
                                  -score(a).objectives["adaptation"],
                                  -score(a).objectives["synthesis_fitness"]))
        beam = fresh[:beam_width]

    candidates = list(evaluated.values())
    feasible = [c for c in candidates if c.feasible]
    front_indices = nondominated([c.objectives for c in feasible], DIRECTIONS)
    front = [feasible[i] for i in front_indices]

    fidelities = {round(c.objectives["fidelity"], 6) for c in feasible}
    recommended = recommend(front)
    return {
        # What this front is *about*. A front is a set of assignments plus objective values
        # measured against one specific inventory under one specific context. Applied to a
        # different inventory the assignment may still rebuild -- the overhang names exist in
        # both -- while every objective value becomes false. Selection verifies this block.
        "binding": {
            "inventory_version": inv.version,
            "inventory_label": inv.label,
            "interface_assignment": inv.interface_assignment,
            "targets": sorted(targets),
            "codon_table_sha256": _table_identity(table),
            "genetic_code": genetic_code,
            "matrix": matrix,
            "destination": destination,
            "k": k,
            "fidelity_tolerance": FIDELITY_TOLERANCE,
            # Recorded so selection can reconstruct the policy this front was measured
            # under. Without it, selection silently rechecked every candidate against
            # the default and 7 of 16 broad-policy points came back infeasible.
            "synthesis_profile": _resolve_profile(profile).as_dict(),
            "source_sha256": _source_identity(),
        },
        "classes": [{"key": c.key, "stage": c.stage, "amino_acids": c.amino_acids,
                     "offset": c.offset, "deposited": c.deposited, "roles": list(c.roles),
                     "permitted": len(c.options)} for c in classes],
        "search_scope": ("level-0 junction classes only; block joins belong to a level-1 "
                         "reaction whose block-plasmid context is not supplied"),
        "targets": list(targets),
        "incumbent": base.as_dict(),
        "evaluated": len(candidates),
        "feasible": len(feasible),
        "infeasible": len(candidates) - len(feasible),
        "observed_front": [c.as_dict() for c in front],
        "recommended": recommended.as_dict() if recommended else None,
        "recommendation_policy": (
            "feasible only; within {tol} of the best observed feasible fidelity; then "
            "highest adaptation; then highest synthesis fitness; then lexicographic "
            "assignment"
        ).format(tol=FIDELITY_TOLERANCE),
        "fidelity_degenerate": len(fidelities) <= 1,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "budget_seconds": wall_seconds, "max_evaluations": max_evaluations,
        "budget_exhausted": budget_exhausted,
        "completion": "budget_exhausted" if budget_exhausted else "complete",
        "scope": ("observed nondominated set among evaluated candidates; not the globally "
                  "optimal front, because the feasible space was not enumerated"),
    }


def _table_identity(table) -> str:
    """Content hash of the codon table as the search actually received it."""
    import hashlib

    return hashlib.sha256(
        json.dumps(table, sort_keys=True, default=str).encode()).hexdigest()


def _source_identity() -> str:
    from .manifest import source_fingerprint

    return source_fingerprint()["sha256"]


def binding_mismatches(binding: dict, inv: Inventory, table, *, targets=None,
                       genetic_code: int = 1, matrix: str = "BsaI-HFv2",
                       destination: str = "level0", k: int = 20) -> list[str]:
    """Every way the supplied context differs from the one the front was measured under.

    Empty means the front's objective values describe this inventory. Anything else means
    they do not, and publishing them alongside these sequences would report one design's
    numbers for another's DNA.
    """
    if not binding:
        return ["the front carries no binding block, so it cannot be matched to an inventory"]

    checks = [
        ("inventory_version", binding.get("inventory_version"), inv.version),
        ("interface_assignment", binding.get("interface_assignment"),
         inv.interface_assignment),
        ("codon_table_sha256", binding.get("codon_table_sha256"), _table_identity(table)),
        ("genetic_code", binding.get("genetic_code"), genetic_code),
        ("matrix", binding.get("matrix"), matrix),
        ("destination", binding.get("destination"), destination),
        ("k", binding.get("k"), k),
        ("source_sha256", binding.get("source_sha256"), _source_identity()),
    ]
    out = [f"{name}: front has {was!r}, this context is {now!r}"
           for name, was, now in checks if was is not None and was != now]
    if targets is not None and binding.get("targets") not in (None, sorted(targets)):
        out.append(f"targets: front has {binding['targets']}, this context is "
                   f"{sorted(targets)}")
    return out


def recommend(front: list[Candidate]) -> Candidate | None:
    """The single default, chosen by a rule declared before the search ran."""
    feasible = [c for c in front if c.feasible and c.objectives]
    if not feasible:
        return None
    best = max(c.objectives["fidelity"] for c in feasible)
    near = [c for c in feasible
            if c.objectives["fidelity"] >= best - FIDELITY_TOLERANCE]
    near.sort(key=lambda c: (-c.objectives["adaptation"], -c.objectives["synthesis_fitness"],
                             tuple(sorted(c.assignment.items()))))
    return near[0]
