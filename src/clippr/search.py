"""Explore the design space, then justify one answer — instead of taking the first that works.

**The weakness this closes.** `design_oneshot` ranks assembly plans by predicted ligation
fidelity, keeps the first whose coding sequence satisfies every constraint, and discards the
rest unexamined. The result is defensible but it cannot answer the obvious question: *is this
design good relative to what was available?* The alternatives were computed and thrown away.

**Why ranking on fidelity alone is not enough.** Fidelity is scored *before* the coding
sequence exists, so a plan that looks best can turn out to force poor codon usage, while a
slightly worse-scoring plan yields a much better sequence:

    plan   fidelity   resulting sequence
    A        0.92     poor
    B        0.91     excellent

Taking A because 0.92 > 0.91 is a decision made on incomplete information. This module
evaluates several plans all the way through codon optimisation and QC, so the comparison is
between finished designs rather than between predictions about them.

**A declared priority hierarchy, not arbitrary weights.** A weighted sum such as
`0.4·fidelity + 0.3·codon + 0.3·synthesis` hides its assumptions in three numbers nobody can
defend. This module uses a lexicographic order instead, so every assumption is visible and
arguable:

    L1  hard feasibility     every constraint satisfied and QC not FAIL. Non-negotiable.
    L2  fidelity band        within `fidelity_tolerance` of the best feasible fidelity.
    L3  sequence quality     prefer QC PASS over WARNING, then higher optimiser score.
    L4  tie-break            higher fidelity.

L2 is an epsilon-constraint rather than an absolute floor because predicted fidelity is capped
by the destination pair, so what matters is distance from the best *reachable* value, not from
1.000.

**One answer, with the alternatives shown.** A Pareto front hands a wet lab a multi-objective
arbitration problem it did not ask for and is not equipped to settle. `certificate()` instead
reports the selected design, how far it sits from the best available value on each objective,
whether anything evaluated dominates it, and what the nearest alternatives would cost. The
decision rule is explicit and the alternatives are visible; the choice is still made.

**This does not replace `design_oneshot`.** That function stays exactly as it is, including its
byte-level agreement with the reference corpus on 200 designs. `explore()` reuses its own plan
enumeration via `design.plan_candidates`, so the candidates examined here are precisely the ones
it would have discarded. The two entry points can legitimately disagree: that disagreement is
the point, and `certificate()` says so when it happens.

**Measured: the trade-off surface is degenerate here, and that is a finding.** Running this
over all three architectures with `budget=8`, predicted fidelity was **0.8297 for every single
candidate** and QC was PASS for every single candidate. Only the optimiser score varied:

    architecture   candidates   fidelity spread   QC        optimiser score spread
    9S                      8   none (0.8297)     all PASS  -150.5 to -139.7
    14S                     8   none (0.8297)     all PASS  -247.7 to -227.1
    19S                     8   none (0.8297)     all PASS  -333.2 to -322.0

There are two reasons and both are structural. Predicted fidelity is capped by the destination
overhang pair, and every junction set this pipeline proposes reaches that cap, so the L2 band
never binds and `fidelity_tolerance` has no effect in this configuration. QC passes universally
because `arelf.safe_overhangs` has already excluded the overhangs that would have caused
trouble.

So **the hierarchy is currently single-objective in practice**. Do not describe this module as
a multi-objective optimiser: on this pipeline, at these destination overhangs, there is nothing
to trade off. What it genuinely delivers is narrower and still worth having:

  * it selects the best of N *finished* designs on the one axis that varies, which changed the
    answer for all three architectures tested — a 1-4% improvement in the optimiser objective
    over the plan `design_oneshot` happened to reach first;
  * it can answer "is this design good relative to the alternatives?" with evidence rather
    than assertion;
  * it establishes the degeneracy above, which is a real property of the design space and was
    previously assumed rather than measured.

The hierarchy is kept in full because it is the correct structure the moment any of those
objectives starts discriminating — a different destination pair, a stricter enzyme profile, a
harder protein — and because a decision rule that only works when the answer is easy is not a
decision rule.

**Cost.** One full codon optimisation per candidate, so runtime is roughly `budget` times a
single design: measured 4.5 s for 9S, 10.5 s for 14S and 48.3 s for 19S at `budget=8`. The
default budget of 6 keeps a one-off design interactive; raise it for a library run.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from . import constants as C
from .assembly import build_oligos
from .codons import optimize_cds
from .qc import synthesis_qc

#: QC verdicts in preference order. FAIL is excluded at L1, never traded away at L3.
_QC_RANK = {"PASS": 0, "WARNING": 1, "FAIL": 2}

#: Default width of the L2 fidelity band, in absolute predicted-fidelity units.
FIDELITY_TOLERANCE = 0.02

#: Candidates to evaluate in full. Each costs one codon optimisation.
DEFAULT_BUDGET = 6


@dataclass(frozen=True)
class Candidate:
    """One assembly plan carried all the way through to a finished, scored design."""

    cuts: tuple[int, ...]
    overhangs: tuple[str, ...]
    fidelity: float
    cds: str
    constraints_ok: bool
    objectives_score: float
    qc: dict = field(repr=False)

    @property
    def feasible(self) -> bool:
        """L1. A FAIL verdict is a hard exclusion, not a cost to be traded off."""
        return bool(self.constraints_ok) and self.qc.get("status") != "FAIL"

    @property
    def qc_rank(self) -> int:
        return _QC_RANK.get(self.qc.get("status", "FAIL"), 2)

    def dominates(self, other: "Candidate") -> bool:
        """Better or equal on every objective and strictly better on at least one.

        Objectives are maximise-fidelity, maximise-optimiser-score, minimise-QC-rank.
        """
        at_least = (self.fidelity >= other.fidelity
                    and self.objectives_score >= other.objectives_score
                    and self.qc_rank <= other.qc_rank)
        strictly = (self.fidelity > other.fidelity
                    or self.objectives_score > other.objectives_score
                    or self.qc_rank < other.qc_rank)
        return at_least and strictly


def explore(protein: str, codon_table, *, n_fragments: int, destination, matrix: str,
            genetic_code: int = 1, seed: int = 42,
            enzyme_profile=C.DEFAULT_ENZYME_PROFILE,
            budget: int = DEFAULT_BUDGET, on_progress=None) -> list[Candidate]:
    """Carry the top `budget` assembly plans through codon optimisation and QC.

    Evaluates finished designs rather than predictions about them. Plans come from
    `design.plan_candidates`, the same enumeration `design_oneshot` uses.
    """
    from .design import plan_candidates

    # Ask for as many plans as the budget can evaluate. Left at its default the enumeration
    # stops after three ceiling-scoring plans, which is all `design_oneshot` needs and far
    # too few to call a search.
    plans, _ceiling = plan_candidates(protein, n_fragments, destination, matrix,
                                      enzyme_profile, want=budget)
    out: list[Candidate] = []
    for i, (cuts, overhangs, fidelity) in enumerate(plans[:budget], 1):
        locked = {3 * c - 4: o for c, o in zip(cuts, overhangs)}
        opt = optimize_cds(protein, locked_sites=locked, codon_table=codon_table,
                           genetic_code=genetic_code, seed=seed, enzymes=enzyme_profile)
        out.append(Candidate(
            cuts=tuple(cuts), overhangs=tuple(overhangs), fidelity=float(fidelity),
            cds=opt["cds"], constraints_ok=bool(opt["constraints_ok"]),
            objectives_score=float(opt.get("objectives_score", 0.0)),
            qc=synthesis_qc(opt["cds"]),
        ))
        if on_progress:
            on_progress(i, min(budget, len(plans)), out[-1])
    return out


def select(candidates: list[Candidate], *,
           fidelity_tolerance: float = FIDELITY_TOLERANCE
           ) -> tuple[Candidate, list[Candidate]]:
    """Apply the hierarchy. Returns (selected, ranked alternatives).

    Raises if nothing is feasible, because there is no defensible answer in that case and
    silently returning the least-bad design would hide a failure.
    """
    if not candidates:
        raise ValueError("no candidates to select from")
    feasible = [c for c in candidates if c.feasible]
    if not feasible:
        raise ValueError(
            f"none of the {len(candidates)} evaluated designs satisfied every constraint "
            f"with a QC verdict better than FAIL")

    best_fid = max(c.fidelity for c in feasible)
    band = [c for c in feasible if c.fidelity >= best_fid - fidelity_tolerance]
    band.sort(key=lambda c: (c.qc_rank, -c.objectives_score, -c.fidelity))
    return band[0], band[1:] + [c for c in feasible if c not in band]


def certificate(selected: Candidate, candidates: list[Candidate], *,
                fidelity_tolerance: float = FIDELITY_TOLERANCE, limit: int = 3) -> str:
    """Say what was considered, why this won, and what the alternatives would cost.

    Deliberately not a Pareto plot: it names the decision rule, the selected design's
    distance from the best available value on each objective, and the price of the nearest
    alternatives.
    """
    feasible = [c for c in candidates if c.feasible]
    lines = [
        f"Evaluated {len(candidates)} complete designs, {len(feasible)} feasible.",
        "",
        "Decision rule, applied in order:",
        "  1. every constraint satisfied and QC not FAIL",
        f"  2. predicted fidelity within {fidelity_tolerance:.3f} of the best feasible value",
        "  3. prefer QC PASS, then the higher codon/uniqueness optimiser score",
        "  4. tie-break on higher fidelity",
        "",
    ]
    if not feasible:
        return "\n".join(lines + ["No feasible design; nothing selected."])

    best_fid = max(c.fidelity for c in feasible)
    best_obj = max(c.objectives_score for c in feasible)
    lines += [
        "Selected:",
        f"  cuts {list(selected.cuts)}  overhangs {list(selected.overhangs)}",
        f"  fidelity          {selected.fidelity:.4f}"
        + (f"  ({best_fid - selected.fidelity:+.4f} vs the best feasible)"
           if selected.fidelity < best_fid else "   (best available)"),
        f"  optimiser score   {selected.objectives_score:.1f}"
        + (f"  ({selected.objectives_score - best_obj:+.1f} vs the best feasible)"
           if selected.objectives_score < best_obj else "   (best available)"),
        f"  QC                {selected.qc.get('status')}",
    ]

    dominators = [c for c in feasible if c.dominates(selected)]
    lines.append(f"  dominated by      {len(dominators)} of the evaluated designs"
                 + ("" if dominators else "  (not dominated)"))

    others = [c for c in feasible if c is not selected][:limit]
    if others:
        lines += ["", f"Nearest alternatives ({len(others)} of {len(feasible) - 1}):"]
        for c in others:
            trade = []
            d_fid = c.fidelity - selected.fidelity
            d_obj = c.objectives_score - selected.objectives_score
            if abs(d_fid) > 1e-9:
                trade.append(f"fidelity {d_fid:+.4f}")
            if abs(d_obj) > 1e-9:
                trade.append(f"optimiser score {d_obj:+.1f}")
            if c.qc.get("status") != selected.qc.get("status"):
                trade.append(f"QC {c.qc.get('status')}")
            lines.append(f"  cuts {list(c.cuts)}: "
                         + (", ".join(trade) if trade else "identical on every objective"))
    else:
        lines += ["", "No feasible alternative was found within the budget."]

    lines += [
        "",
        "Predicted fidelity comes from published ligation matrices, not measured assembly, "
        "and the optimiser score is comparable across candidates for this protein only. "
        "Neither is an expression or synthesis-success claim.",
    ]
    return "\n".join(lines)


def design_searched(target_rna: str, *, budget: int = DEFAULT_BUDGET,
                    fidelity_tolerance: float = FIDELITY_TOLERANCE,
                    on_progress=None, **kwargs) -> dict:
    """`design_oneshot`, but choosing among fully evaluated alternatives.

    Returns the same shape as `design_oneshot` plus `certificate`, `alternatives` and
    `n_evaluated`, and a `differs_from_oneshot` flag saying whether exploring changed the
    answer. Accepts the same keyword arguments as `design_oneshot`.
    """
    from .design import design_oneshot

    base = design_oneshot(target_rna, **kwargs)

    # Re-optimise under the settings the base design actually resolved, not the ones the
    # caller passed: `organism` may have selected the table, and `extra_blacklist` is merged
    # into the enzyme profile. Using the requested arguments here would silently explore a
    # different constraint set from the one being compared against.
    cands = explore(
        base["protein"], base["codon_table"],
        n_fragments=base["n_fragments"], destination=kwargs.get("destination"),
        matrix=kwargs.get("matrix", "BsaI-HFv2"),
        genetic_code=base["genetic_code"], seed=kwargs.get("seed", 42),
        enzyme_profile=base["enzyme_profile_effective"],
        budget=budget, on_progress=on_progress,
    )
    chosen, alts = select(cands, fidelity_tolerance=fidelity_tolerance)

    out = dict(base)
    if chosen.cds != base["cds"]:
        oligos = build_oligos(chosen.cds, list(chosen.cuts), list(chosen.overhangs),
                              kwargs.get("destination"),
                              kwargs.get("enzyme", "BsaI"))
        out.update({"cds": chosen.cds, "cuts": list(chosen.cuts),
                    "junction_overhangs": list(chosen.overhangs),
                    "fidelity": chosen.fidelity, "qc": chosen.qc, "oligos": oligos})
    out.update({
        "certificate": certificate(chosen, cands,
                                   fidelity_tolerance=fidelity_tolerance),
        "alternatives": alts,
        "n_evaluated": len(cands),
        "differs_from_oneshot": chosen.cds != base["cds"],
    })
    return out
