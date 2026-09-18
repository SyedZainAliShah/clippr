"""Optimise an inventory as a collection, rather than one module at a time.

Takes an existing inventory and a context; returns a new compatible inventory and a comparison.
Callable without redesigning interfaces and without rerunning target selection.

**Why a collection objective is not the sum of per-module ones.** Recoding each module in
isolation maximises each module's own codon adaptation and is blind to what the modules share
with each other. A collection objective can trade a little adaptation in one module for less
shared sequence across the set. Whether that trade is worth making is a measurement, reported
here, not an assumption built into the code.

**Three modes, all through one interface**, so they can be compared on equal terms:

  ``greedy``     accept only improvements; refine each module against the rest
  ``anneal``     the same neighbourhood, explored with a Metropolis acceptance rule
  ``random``     a control: same budget, proposals accepted without regard to the objective

The control exists because an optimiser that beats nothing has not been shown to work. If
``anneal`` cannot beat ``random`` on the same budget, that is the finding.

**All three modes are reproducible for a fixed seed and settings, and none is seed-free
deterministic.** Even ``greedy`` draws its candidate replacements from a seeded RNG -- it is
deterministic in its *acceptance* rule, not in the proposals it gets to accept. A run reported
without its seed, processing order and budget is not reproducible.

**Removing before evaluating.** When a module's replacement is scored, the module's own current
contribution is taken out of the collection first. Otherwise every candidate is measured
against a collection that still contains the thing it is replacing, and the incumbent's own
sharing counts against its replacement.

**Hard constraints stay outside the objective** -- `docs/objectives.md` §1. A proposal that
violates one is rejected outright, never traded against a better score.
"""
from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass, field

from .inventories import INTERFACE, Inventory, derive

MODES = ("greedy", "anneal", "random")

#: Normalised terms and their weights for the scalar objective. Both terms are mapped to
#: "higher is better" in [0, 1]-ish ranges before weighting, so the weights mean what they
#: look like. These are CLIPPR engineering defaults chosen from our own measurements, not
#: values taken from any other implementation.
DEFAULT_WEIGHTS = {"adaptation": 1.0, "sharing": 1.0}


@dataclass
class Proposal:
    """One candidate replacement for one module, and why it was or was not usable."""

    module_id: str
    dna: str
    feasible: bool
    reason: str | None = None
    score: float | None = None


@dataclass
class SearchResult:
    """Everything a caller needs to judge the search, including that it did nothing."""

    inventory: Inventory
    mode: str
    incumbent_objectives: dict
    final_objectives: dict
    accepted: list[str] = field(default_factory=list)
    proposals: int = 0
    evaluations: int = 0
    infeasible: int = 0
    elapsed_seconds: float = 0.0
    budget_seconds: float = 0.0
    budget_exhausted: bool = False
    seed: int = 0

    @property
    def improved(self) -> bool:
        return self.final_objectives["scalar"] > self.incumbent_objectives["scalar"] + 1e-12

    def as_dict(self) -> dict:
        return {"mode": self.mode, "seed": self.seed,
                "inventory_version": self.inventory.version,
                "incumbent": self.incumbent_objectives, "final": self.final_objectives,
                "improved": self.improved,
                "accepted_modules": self.accepted, "accepted": len(self.accepted),
                "proposals": self.proposals, "evaluations": self.evaluations,
                "infeasible": self.infeasible,
                "elapsed_seconds": round(self.elapsed_seconds, 3),
                "budget_seconds": self.budget_seconds,
                "budget_exhausted": self.budget_exhausted,
                "completion": "budget_exhausted" if self.budget_exhausted else "complete"}


def collection_objectives(records: dict[str, str], frames: dict[str, int], table: dict, *,
                          k: int = 20, weights: dict | None = None) -> dict:
    """The collection's objective vector and the scalar the search minimises against.

    `adaptation` is the mean CAI over the modules' own coding spans. `sharing` is the k-mer
    collection penalty, normalised by the number of distinct record pairs so the term does not
    grow simply because the inventory is larger.
    """
    from .objectives import codon_adaptation, collection_penalty

    weights = {**DEFAULT_WEIGHTS, **(weights or {})}
    cais = []
    for module_id, dna in records.items():
        frame = frames[module_id]
        n = (len(dna) - frame) // 3
        cais.append(codon_adaptation(dna[frame:frame + 3 * n], table)["cai"])
    adaptation = sum(cais) / len(cais) if cais else 0.0

    penalty = collection_penalty(records, k=k)["penalty"]
    pairs = max(1, len(records) * (len(records) - 1) // 2)
    per_pair = penalty / pairs
    # Mapped to "higher is better" so both terms point the same way before weighting.
    sharing_term = 1.0 / (1.0 + per_pair)

    scalar = (weights["adaptation"] * adaptation + weights["sharing"] * sharing_term) \
        / (weights["adaptation"] + weights["sharing"])
    return {"adaptation": round(adaptation, 6), "sharing_penalty": penalty,
            "sharing_per_pair": round(per_pair, 4), "sharing_term": round(sharing_term, 6),
            "scalar": round(scalar, 6), "k": k, "modules": len(records)}


def _propose(module_id: str, record, table, rng: random.Random,
             genetic_code: int, profile=None) -> Proposal:
    """One synonymous replacement for one module, already checked against the contract."""
    from Bio.Seq import Seq

    from .codons import optimize_cds
    from .recoding import locked_interface_sites
    from .substrates import substrate_problems

    start, end = record.coding_interval
    protein = str(Seq(record.dna[start:end]).translate())
    if "*" in protein:
        return Proposal(module_id, record.dna, False, "stop codon in frame")

    # The policy reaches the solver, not only the validator below. Accepting a profile here
    # and then letting `optimize_cds` fall back to its own defaults is how a requested policy
    # became a strict search with a permissive check.
    from .synthesis_profile import resolve as _resolve

    bounds = _resolve(profile).solver_bounds()
    try:
        result = optimize_cds(protein, locked_sites=locked_interface_sites(
            record.dna, record.frame), codon_table=table, genetic_code=genetic_code,
            seed=rng.randrange(1 << 30), unique_kmer_size=None, **bounds)
    except Exception as exc:                          # noqa: BLE001 - recorded, not hidden
        return Proposal(module_id, record.dna, False, f"{type(exc).__name__}: {exc}")
    if not result.get("constraints_ok"):
        return Proposal(module_id, record.dna, False, "optimiser reported constraints unmet")

    candidate = record.dna[:start] + result["cds"] + record.dna[end:]
    # Re-verified against the module contract, not trusted from the optimiser's verdict.
    if len(candidate) != len(record.dna):
        return Proposal(module_id, candidate, False, "length changed")
    if candidate[:INTERFACE] != record.five_interface \
            or candidate[-INTERFACE:] != record.three_interface:
        return Proposal(module_id, candidate, False, "interface changed")
    if str(Seq(candidate[start:end]).translate()) != protein:
        return Proposal(module_id, candidate, False, "protein changed")
    # The substrate that would actually be ordered, not the bare insert. Judging candidates
    # on the insert let this optimiser reintroduce an out-of-band substrate in 16 of 18
    # matched runs, on an inventory the recoder had already cleaned.
    problems = substrate_problems(candidate, record.block, profile)
    if problems:
        return Proposal(module_id, candidate, False, problems[0])
    return Proposal(module_id, candidate, True)


def optimise_library(inv: Inventory, codon_table: dict, *, mode: str = "greedy",
                     profile=None,
                     label: str | None = None, genetic_code: int = 1, seed: int = 42,
                     k: int = 20, weights: dict | None = None,
                     wall_seconds: float = 600.0, max_proposals: int = 400,
                     initial_temperature: float = 0.02,
                     final_temperature: float = 0.0005,
                     reverse_order: bool = False) -> SearchResult:
    """Refine an inventory as a collection. Returns the incumbent when nothing improves.

    Both a time budget and a proposal budget apply; the first reached stops new work. The
    result says which, because "no improvement found" and "ran out of budget" are different
    claims and only one of them is about the search space.
    """
    if mode not in MODES:
        raise ValueError(f"mode {mode!r} is not one of {MODES}")

    from .codons import complete_table

    table = complete_table(codon_table, genetic_code)
    frames = {m: r.frame for m, r in inv.modules.items()}
    records = {m: r.dna for m, r in inv.modules.items()}
    incumbent = collection_objectives(records, frames, table, k=k, weights=weights)

    rng = random.Random(seed)
    order = sorted(inv.modules, reverse=reverse_order)
    started = time.perf_counter()
    accepted: list[str] = []
    proposals = evaluations = infeasible = 0
    current = dict(records)
    current_score = incumbent["scalar"]
    best, best_score = dict(records), current_score
    budget_exhausted = False
    observed_deltas: list[float] = []
    temperature_scale = initial_temperature

    # Every mode sweeps the inventory until the same budget is spent. Capping greedy at one
    # pass while annealing swept four was not a comparison of methods -- it was a comparison
    # of budgets, and the plan requires them matched. Greedy accepts only improvements, so
    # extra sweeps can only help it; giving them is the fair test.
    sweeps = max(1, math.ceil(max_proposals / max(1, len(order))))
    visits = [m for _ in range(sweeps) for m in order]

    for step, module_id in enumerate(visits):
        if time.perf_counter() - started >= wall_seconds or proposals >= max_proposals:
            budget_exhausted = step < len(visits) - 1
            break

        proposal = _propose(module_id, inv.modules[module_id], table, rng,
                            genetic_code, profile)
        proposals += 1
        if not proposal.feasible:
            infeasible += 1
            continue

        # Remove this module's own contribution before scoring its replacement: otherwise
        # the candidate is measured against a collection that still contains the incumbent,
        # and the incumbent's sharing counts against the thing replacing it.
        trial = dict(current)
        trial[module_id] = proposal.dna
        scored = collection_objectives(trial, frames, table, k=k, weights=weights)
        evaluations += 1
        delta = scored["scalar"] - current_score
        if delta:
            observed_deltas.append(abs(delta))

        if mode == "greedy":
            take = delta > 0
        elif mode == "random":
            take = rng.random() < 0.5
        else:
            if delta > 0:
                take = True
            else:
                # The schedule is calibrated to the objective's own scale. A fixed absolute
                # temperature accepted ~90% of worsening moves here, because this objective's
                # deltas are ~1e-3 while the temperature was 2e-2: exp(-1e-3/2e-2) ~ 0.95.
                # That is a random walk wearing an annealing schedule. Anchoring T to the
                # median observed |delta| makes the acceptance probability mean something
                # whatever the objective's units happen to be.
                if len(observed_deltas) >= 5:
                    temperature_scale = statistics.median(observed_deltas)
                fraction = step / max(1, len(visits) - 1)
                temperature = temperature_scale * (
                    (final_temperature / initial_temperature) ** fraction)
                take = temperature > 0 and rng.random() < math.exp(delta / temperature)

        if take:
            current, current_score = trial, scored["scalar"]
            accepted.append(module_id)
            # Best-so-far is tracked independently of the accepted state, so an annealing
            # run that wanders uphill and back still returns the best thing it ever saw.
            if current_score > best_score:
                best, best_score = dict(current), current_score

    elapsed = time.perf_counter() - started
    final = collection_objectives(best, frames, table, k=k, weights=weights)
    replacements = {m: dna for m, dna in best.items() if dna != records[m]}
    return SearchResult(
        inventory=derive(inv, label or f"{inv.label}+{mode}", replacements),
        mode=mode, incumbent_objectives=incumbent, final_objectives=final,
        accepted=accepted, proposals=proposals, evaluations=evaluations,
        infeasible=infeasible, elapsed_seconds=elapsed, budget_seconds=wall_seconds,
        budget_exhausted=budget_exhausted, seed=seed)
