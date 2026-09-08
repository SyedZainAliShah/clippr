"""Design a set of mutually orthogonal PPR targets for a regulator library.

For a CLIPPR library, PPR_i must bind 5'UTR_i and not 5'UTR_j. Cross-talk destroys
independent control. Because the 5'UTRs are designed rather than given by biology, the
targets can be chosen to be maximally distinguishable in the first place - a freedom
GRASP does not have, since it targets whichever transcript nature supplies.

This is max-min distance selection over sequence space: choose N sequences maximising the
smallest pairwise distance between any two. Two metrics are offered:

  uniform   plain Hamming distance. Safe default, no assumptions.
  weighted  Hamming weighted by per-position cognate effects measured in the GRASP
            perturbation assay. Opt-in - see the caveat on POSITION_WEIGHTS below.

Nothing here predicts binding. It maximises sequence separation, which is a necessary but
not sufficient condition for orthogonality.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field

BASES = "ACGU"

#: Per-position cognate effects measured in the GRASP perturbation assay and re-analysed
#: by this project. Keys are that assay's position indices, 5' to 3' over the contacted
#: bases; values are the estimated effect on editing proportion.
#:
#: CAVEAT. These were measured on a 10-motif protein against a 15-nt RNA in an editing
#: assay, with three protein variants per position and wide intervals. Using them to
#: weight a designed target of different length assumes the positional pattern transfers,
#: which we have not established. The global heterogeneity test was p = 0.0066 and the
#: refinement did not improve held-out cognate ranking. Treat weighting as a hypothesis,
#: not a calibrated model - which is why `uniform` is the default.
POSITION_WEIGHTS: dict[int, float] = {
    1: 0.0725, 2: 0.1082, 3: 0.0641, 4: 0.0482, 5: 0.0218,
    6: 0.0184, 7: 0.0494, 8: 0.0120, 9: 0.0232, 10: 0.0096,
}


@dataclass
class TargetSet:
    targets: list[str]
    metric: str
    min_distance: float
    mean_distance: float
    worst_pair: tuple[int, int]
    distances: list[list[float]] = field(repr=False, default_factory=list)

    def report(self) -> str:
        lines = [
            f"{len(self.targets)} targets, length {len(self.targets[0])}, metric '{self.metric}'",
            f"  minimum pairwise distance : {self.min_distance:.4f}",
            f"  mean pairwise distance    : {self.mean_distance:.4f}",
            f"  closest pair              : {self.worst_pair[0]} / {self.worst_pair[1]}",
            "",
        ]
        for i, t in enumerate(self.targets):
            lines.append(f"  [{i:2}] {t}")
        return "\n".join(lines)


def _weights(length: int, metric: str) -> list[float]:
    if metric == "uniform":
        return [1.0] * length
    if metric == "weighted":
        w = [POSITION_WEIGHTS.get(i + 1, 0.0) for i in range(length)]
        total = sum(w) or 1.0
        return [x * length / total for x in w]      # normalise so scales are comparable
    raise ValueError(f"unknown metric {metric!r}; use 'uniform' or 'weighted'")


def distance(a: str, b: str, w: list[float]) -> float:
    """Weighted Hamming distance between two equal-length targets."""
    if len(a) != len(b):
        raise ValueError("targets must be the same length")
    return sum(wi for x, y, wi in zip(a, b, w) if x != y)


def _matrix(targets: list[str], w: list[float]) -> list[list[float]]:
    n = len(targets)
    m = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = distance(targets[i], targets[j], w)
            m[i][j] = m[j][i] = d
    return m


def _score(targets: list[str], w: list[float]) -> tuple[float, tuple[int, int]]:
    """Minimum pairwise distance and the pair achieving it."""
    best, pair = float("inf"), (0, 1)
    for i in range(len(targets)):
        for j in range(i + 1, len(targets)):
            d = distance(targets[i], targets[j], w)
            if d < best:
                best, pair = d, (i, j)
    return best, pair


def design_orthogonal_set(
    n: int,
    length: int = 9,
    metric: str = "uniform",
    seed: int = 0,
    iterations: int = 20000,
    forbid: list[str] | None = None,
) -> TargetSet:
    """Choose `n` targets of `length` nt maximising the minimum pairwise distance.

    Greedy seed followed by local search: repeatedly replace one member of the closest
    pair with a random alternative, keeping the change only if the minimum improves.

    `forbid` optionally excludes sequences (e.g. motifs already used elsewhere in the
    construct, or sequences with unwanted structure).
    """
    if n < 2:
        raise ValueError("need at least 2 targets")
    if n > 4 ** length:
        raise ValueError(f"cannot fit {n} distinct targets of length {length}")
    rng = random.Random(seed)
    w = _weights(length, metric)
    banned = set(forbid or [])

    def sample() -> str:
        while True:
            s = "".join(rng.choice(BASES) for _ in range(length))
            if s not in banned:
                return s

    # greedy seed: each new target is the best of a batch of candidates
    targets = [sample()]
    while len(targets) < n:
        best_cand, best_d = None, -1.0
        for _ in range(200):
            c = sample()
            if c in targets:
                continue
            d = min(distance(c, t, w) for t in targets)
            if d > best_d:
                best_cand, best_d = c, d
        targets.append(best_cand)

    # local search on the binding constraint
    cur, pair = _score(targets, w)
    for _ in range(iterations):
        i = pair[rng.randint(0, 1)]
        old = targets[i]
        cand = sample()
        if cand in targets:
            continue
        targets[i] = cand
        new, new_pair = _score(targets, w)
        if new > cur:
            cur, pair = new, new_pair
        else:
            targets[i] = old

    m = _matrix(targets, w)
    off = [m[i][j] for i in range(n) for j in range(i + 1, n)]
    lo, worst = _score(targets, w)
    return TargetSet(
        targets=targets, metric=metric, min_distance=lo,
        mean_distance=sum(off) / len(off), worst_pair=worst, distances=m,
    )


def crosstalk_report(targets: list[str], metric: str = "uniform") -> str:
    """Pairwise distance matrix for an existing target set."""
    w = _weights(len(targets[0]), metric)
    m = _matrix(targets, w)
    head = "      " + "".join(f"{i:>7}" for i in range(len(targets)))
    rows = [head]
    for i, row in enumerate(m):
        rows.append(f"  [{i:2}]" + "".join(
            "      ." if i == j else f"{row[j]:7.2f}" for j in range(len(targets))))
    return "\n".join(rows)
