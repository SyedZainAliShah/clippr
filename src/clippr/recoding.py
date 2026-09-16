"""Recode an inventory for a host context, with every interface frozen.

Takes an inventory and a design context, returns a new inventory version. The protein, the
length and both four-base interfaces of every module are preserved, so a recoded module drops
into an assembly built from unrecoded neighbours.

**Why freezing the interfaces is enough.** A module's insert is not a standalone coding
sequence: codons run across the joins. But the leading and trailing partial codons always fall
inside the four-base interfaces, so freezing those also freezes every boundary-spanning base.
A neighbour's codons cannot be disturbed by recoding this module.

**An unchanged module is a delivered module.** If no candidate improves on the incumbent, the
incumbent is kept and marked `unchanged`. Coverage is counted over every module the inventory
holds, not over the ones that happened to move — otherwise a weak result would look like a
strong one on a smaller denominator.

Budgets are both time and attempts, and the first reached stops new work. A module skipped for
budget is reported as not attempted, never silently folded into "unchanged".
"""
from __future__ import annotations

import time
from dataclasses import dataclass

from . import constants as C
from .inventories import INTERFACE, Inventory, derive

#: The synthesis contract, applied to the **whole module** rather than to the sub-span the
#: optimiser was handed. `optimize_cds` constrains the coding sequence it is given; a module's
#: frozen flanking bases sit outside that span, so a GC window straddling them can exceed the
#: band while the optimiser correctly reports its own span clean. The ordered fragment is the
#: whole module, so the whole module is what has to satisfy the contract.
GC_BAND = (0.35, 0.65)
GC_WINDOW = 50
MAX_HOMOPOLYMER = 4
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass
class RecodingReport:
    """What happened to every module, and what it cost."""

    inventory: Inventory
    improved: list[str]
    unchanged: list[str]
    infeasible: list[str]
    not_attempted: list[str]
    objectives: dict[str, dict]
    elapsed_seconds: float
    budget_seconds: float
    budget_exhausted: bool

    def as_dict(self) -> dict:
        return {"inventory_version": self.inventory.version,
                "inventory_label": self.inventory.label,
                "modules": len(self.inventory),
                "improved": self.improved, "unchanged": self.unchanged,
                "infeasible": self.infeasible, "not_attempted": self.not_attempted,
                "objectives": self.objectives,
                "elapsed_seconds": round(self.elapsed_seconds, 3),
                "budget_seconds": self.budget_seconds,
                "budget_exhausted": self.budget_exhausted,
                "coverage": (f"{len(self.improved) + len(self.unchanged)}"
                             f"/{len(self.inventory)} realised")}


def module_constraint_problems(dna: str) -> list[str]:
    """Constraint breaches measured over the complete module sequence.

    Applied to the candidate the recoding would actually deliver, not to the sub-span the
    optimiser optimised. Two of the 42 modules first recoded here passed the optimiser's own
    check at 0.640 and reached 0.660 over the full module, because a window straddling the
    frozen flanking bases includes bases the optimiser never saw.
    """
    from Bio.Restriction import RestrictionBatch

    problems = []
    for name in C.enzymes_for(C.DEFAULT_ENZYME_PROFILE):
        site = str(RestrictionBatch([name]).get(name).site)
        for strand, pattern in (("forward", site),
                                ("reverse", site.translate(_COMPLEMENT)[::-1])):
            at = dna.find(pattern)
            if at >= 0:
                problems.append(f"{name} site on the {strand} strand at {at}")

    for i in range(0, max(1, len(dna) - GC_WINDOW + 1)):
        window = dna[i:i + GC_WINDOW]
        if len(window) < GC_WINDOW:
            break
        gc = (window.count("G") + window.count("C")) / len(window)
        if not GC_BAND[0] <= gc <= GC_BAND[1]:
            problems.append(f"GC {gc:.3f} outside {GC_BAND} in the window at {i}")
            break

    run, prev = 1, ""
    for base in dna:
        run = run + 1 if base == prev else 1
        prev = base
        if run > MAX_HOMOPOLYMER:
            problems.append(f"homopolymer run longer than {MAX_HOMOPOLYMER}")
            break
    return problems


def locked_interface_sites(dna: str, frame: int) -> dict[int, str]:
    """The interface bases that lie inside this module's coding span, as optimiser locks.

    Only the part of each interface inside the coding span is locked. Locking all four bases
    would address positions past the end of the sequence actually being optimised, which the
    optimiser cannot honour and which silently shifts every constraint after it.
    """
    n = (len(dna) - frame) // 3
    start, end = frame, frame + 3 * n
    locked: dict[int, str] = {}
    head_end = min(INTERFACE, end)
    if head_end > start:
        locked[0] = dna[start:head_end]
    tail_from = max(len(dna) - INTERFACE, start)
    if end > tail_from:
        locked[tail_from - start] = dna[tail_from:end]
    return locked


def recode_inventory(inv: Inventory, codon_table: dict, *, label: str,
                     genetic_code: int = 1, seeds: int = 4,
                     wall_seconds: float = 600.0,
                     per_module_seconds: float = 10.0) -> RecodingReport:
    """Recode every module in `inv`, keeping proteins, lengths and interfaces.

    `seeds` candidates are tried per module and the best feasible one by codon adaptation is
    kept, provided it beats the incumbent. Every candidate is re-verified here — interfaces,
    length and protein — rather than trusted from the optimiser's own verdict.
    """
    from Bio.Seq import Seq

    from .codons import complete_table, optimize_cds
    from .objectives import codon_adaptation
    from .substrates import substrate_problems

    table = complete_table(codon_table, genetic_code)
    started = time.perf_counter()
    replacements: dict[str, str] = {}
    improved, unchanged, infeasible, not_attempted = [], [], [], []
    objectives: dict[str, dict] = {}

    for module_id in sorted(inv.modules):
        record = inv.modules[module_id]
        if time.perf_counter() - started >= wall_seconds:
            not_attempted.append(module_id)
            continue

        start, end = record.coding_interval
        original_coding = record.dna[start:end]
        protein = str(Seq(original_coding).translate())
        if "*" in protein:
            infeasible.append(module_id)
            objectives[module_id] = {"reason": "module carries a stop codon in its frame"}
            continue

        before = codon_adaptation(original_coding, table)
        locked = locked_interface_sites(record.dna, record.frame)

        best, best_score = None, before["cai"]
        for seed in range(42, 42 + seeds):
            if time.perf_counter() - started >= wall_seconds:
                break
            attempt_started = time.perf_counter()
            try:
                result = optimize_cds(protein, locked_sites=locked, codon_table=table,
                                      genetic_code=genetic_code, seed=seed,
                                      enzymes=C.DEFAULT_ENZYME_PROFILE,
                                      unique_kmer_size=None)
            except Exception:                        # noqa: BLE001 - recorded, not hidden
                continue
            if time.perf_counter() - attempt_started > per_module_seconds:
                continue
            if not result.get("constraints_ok"):
                continue

            candidate = record.dna[:start] + result["cds"] + record.dna[end:]
            # Re-verified here rather than trusted: the optimiser reporting success is not
            # the same as the candidate satisfying this module's contract.
            if len(candidate) != len(record.dna):
                continue
            if candidate[:INTERFACE] != record.five_interface:
                continue
            if candidate[-INTERFACE:] != record.three_interface:
                continue
            if str(Seq(candidate[start:end]).translate()) != protein:
                continue
            # Checked on the sequence that is actually synthesised -- the released
            # fragment -- not on the bare insert. An E module's fragment is its insert plus
            # the GC-rich `CGAG` overhang, which pushed a window to 0.660 against a 0.65
            # band on an insert that passed at 0.640. The same class of mistake as
            # constraining only the optimiser's sub-span.
            if substrate_problems(candidate, record.block):
                continue

            score = codon_adaptation(candidate[start:end], table)["cai"]
            if score > best_score:
                best, best_score = candidate, score

        after = before if best is None else codon_adaptation(best[start:end], table)
        objectives[module_id] = {"cai_before": round(before["cai"], 4),
                                 "cai_after": round(after["cai"], 4),
                                 "codons_scored": after["codons"],
                                 "coding_interval": [start, end]}
        if best is None:
            unchanged.append(module_id)
        else:
            replacements[module_id] = best
            improved.append(module_id)

    elapsed = time.perf_counter() - started
    return RecodingReport(
        inventory=derive(inv, label, replacements),
        improved=improved, unchanged=unchanged, infeasible=infeasible,
        not_attempted=not_attempted, objectives=objectives,
        elapsed_seconds=elapsed, budget_seconds=wall_seconds,
        budget_exhausted=bool(not_attempted))
