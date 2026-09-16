"""A synthesis-quality objective that does not degenerate, and why the old one did.

**The problem it replaces.** CLIPPR's third Pareto axis was `repeat_burden` -- duplicated
k-mers within a module. Over the deposited and recoded inventories it is **identically zero**,
and not because of its parameter: measured at k = 20, 16, 12 and 10 it is zero for all 42
modules, and only at k = 8 does anything register. A repeated 10-mer inside a ~100 nt module is
genuinely rare. The search was therefore two-objective while reporting three, and the third
contributed no ordering information at all.

Changing `k` would have looked like a fix and been none. The measure is wrong, not its setting.

**What replaces it.** A composite over quantities that actually vary across our corpus,
measured on the **ordered substrate** rather than the insert:

    GC centrality        how near each window sits to the middle of the declared band.
                         Across our 42 substrates GC spans 0.417-0.574, so this discriminates.
    homopolymer headroom how far the longest run sits below the cap. Spans 3-4 today.
    internal repetition  duplicated k-mers. Degenerate now, retained because it is the thing
                         that fires on a genuinely bad candidate, which is when it matters.
    collection sharing   k-mers shared with the rest of the inventory, when one is supplied.

**Why a continuous preference and not a feasibility test.** A constraint check returns the same
answer for every feasible candidate, so it cannot rank them. Each term here is a distance, so
ranking still moves when everything comfortably passes. That idea is taken from the GRASP
reference's `prefer_ideal`; the formulation below is CLIPPR's own.

**The score is never reported alone.** `synthesis_fitness` returns its components beside the
total, because an aggregate whose parts cannot be inspected is a number nobody can check -- and
the weights are a declared engineering choice, not a measurement.
"""
from __future__ import annotations

from dataclasses import dataclass, field

#: Declared weights. They are a choice, not a finding, and they are reported with every score
#: so that a reader can reweigh rather than take these on trust. They sum to 1.0 so the total
#: lands in [0, 1] and is comparable across inventories of different sizes.
WEIGHTS: dict[str, float] = {
    "gc_centrality": 0.45,
    "homopolymer_headroom": 0.25,
    "internal_repetition": 0.15,
    "collection_sharing": 0.15,
}

#: k for both repetition terms. 20 for continuity with the QC this project already reports;
#: the choice is recorded with every score because two k values are not comparable.
DEFAULT_K = 20


@dataclass(frozen=True)
class Fitness:
    """A composite synthesis score in [0, 1], with the parts that produced it."""

    total: float
    components: dict[str, float]
    weights: dict[str, float]
    k: int
    scope: str = "ordered substrate"
    detail: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {"total": self.total, "components": self.components,
                "weights": self.weights, "k": self.k, "scope": self.scope,
                "detail": self.detail}


def _windows(sequence: str, width: int) -> list[str]:
    if len(sequence) < width:
        return [sequence] if sequence else []
    return [sequence[i:i + width] for i in range(len(sequence) - width + 1)]


def _gc(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / max(1, len(sequence))


def gc_centrality(sequence: str, band: tuple[float, float], width: int) -> float:
    """1.0 when every window sits at the band's midpoint, 0.0 at its edge or beyond.

    Measured against the **declared band's** centre rather than a fixed 0.5: a policy that
    permits 0.15-0.85 is not asking for the same sequence as one permitting 0.35-0.65, and
    scoring both toward 0.5 would impose a preference no profile declared.

    The worst window sets the score, not the mean. One bad window is what a vendor's model
    reacts to, and averaging would let a long comfortable sequence hide it.
    """
    lo, hi = band
    centre = (lo + hi) / 2
    half = (hi - lo) / 2
    if half <= 0:
        return 0.0
    worst = 0.0
    for window in _windows(sequence, width) or [sequence]:
        distance = abs(_gc(window) - centre) / half        # 0 at centre, 1 at the edge
        worst = max(worst, min(1.0, distance))
    return round(1.0 - worst, 6)


def homopolymer_headroom(sequence: str, cap: int) -> float:
    """1.0 for no run longer than 1, falling to 0.0 at the cap and beyond.

    A ramp rather than a threshold, so a candidate with a run of 2 is preferred to one with a
    run of 4 even though a feasibility test calls both acceptable.
    """
    if cap < 1:
        return 0.0
    longest, run, previous = 1, 1, ""
    for base in sequence:
        run = run + 1 if base == previous else 1
        previous = base
        longest = max(longest, run)
    return round(max(0.0, min(1.0, (cap - longest) / max(1, cap - 1))), 6)


def internal_repetition(sequence: str, k: int) -> float:
    """1.0 when every k-mer is unique, falling as duplicates accumulate.

    Degenerate on today's corpus, which is exactly why it is one term of four rather than the
    whole objective. It is retained because a candidate that *is* repetitive is a real problem,
    and an objective that cannot see it would have to be rebuilt later.
    """
    words = [sequence[i:i + k] for i in range(max(0, len(sequence) - k + 1))]
    if not words:
        return 1.0
    duplicates = len(words) - len(set(words))
    return round(1.0 - min(1.0, duplicates / len(words)), 6)


def collection_sharing(sequence: str, others: list[str], k: int) -> float:
    """1.0 when this sequence shares no k-mer with the rest of the collection.

    Normalised by this sequence's own k-mer count, so a long module is not penalised simply
    for having more windows to share.
    """
    words = {sequence[i:i + k] for i in range(max(0, len(sequence) - k + 1))}
    if not words or not others:
        return 1.0
    elsewhere: set[str] = set()
    for other in others:
        elsewhere.update(other[i:i + k] for i in range(max(0, len(other) - k + 1)))
    shared = len(words & elsewhere)
    return round(1.0 - min(1.0, shared / len(words)), 6)


def synthesis_fitness(sequence: str, *, profile=None, others: list[str] | None = None,
                      k: int = DEFAULT_K, weights: dict[str, float] | None = None) -> Fitness:
    """Composite synthesis quality for one ordered sequence. Higher is better.

    `sequence` must be the **ordered substrate**, not the bare insert: the rule this project
    states in `docs/objectives.md` §0 is that a contract is evaluated on what leaves the
    building, and a quality score is no exception.
    """
    from .synthesis_profile import resolve

    policy = resolve(profile)
    used = dict(weights or WEIGHTS)
    total_weight = sum(used.values())
    if total_weight <= 0:
        raise ValueError("weights must sum to something positive")

    components = {
        "gc_centrality": gc_centrality(sequence, policy.local_gc, policy.window),
        "homopolymer_headroom": homopolymer_headroom(sequence, policy.max_homopolymer),
        "internal_repetition": internal_repetition(sequence, k),
        "collection_sharing": collection_sharing(sequence, others or [], k),
    }
    total = sum(used.get(name, 0.0) * value for name, value in components.items())
    return Fitness(total=round(total / total_weight, 6), components=components,
                   weights=used, k=k,
                   detail={"profile": policy.name, "band": list(policy.local_gc),
                           "window": policy.window, "max_homopolymer": policy.max_homopolymer,
                           "collection_size": len(others or [])})


def inventory_fitness(inventory, *, profile=None, k: int = DEFAULT_K) -> dict:
    """Mean composite fitness over an inventory's ordered substrates, with its spread.

    The spread is reported because a mean alone cannot say whether an objective discriminates,
    and discriminating is the entire reason this replaced `repeat_burden`. A standard deviation
    of zero here means the axis has degenerated again and should be reported, not hidden.
    """
    import statistics

    from .substrates import build

    scores, per_module = [], {}
    sequences = {module_id: build(record.dna, record.block).sequence
                 for module_id, record in inventory.modules.items()}
    for module_id, sequence in sorted(sequences.items()):
        others = [s for m, s in sequences.items() if m != module_id]
        got = synthesis_fitness(sequence, profile=profile, others=others, k=k)
        per_module[module_id] = got.as_dict()
        scores.append(got.total)

    spread = statistics.pstdev(scores) if len(scores) > 1 else 0.0
    return {
        "mean": round(statistics.mean(scores), 6) if scores else 0.0,
        "min": round(min(scores), 6) if scores else 0.0,
        "max": round(max(scores), 6) if scores else 0.0,
        "spread": round(spread, 6),
        "degenerate": spread < 1e-9,
        "modules": len(scores),
        "k": k,
        "per_module": per_module,
    }
