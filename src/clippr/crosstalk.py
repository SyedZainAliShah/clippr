"""Two-tier cross-talk assessment for a library of targets.

**Tier A — sequence separation. Hard, and the only criterion that gates anything.**
Weighted Hamming distance between targets. It makes no biological claim: it says two
sequences differ in *k* positions, which is geometry, and needs no model to be true. A
library that passes here is separated in sequence space whatever anyone later learns about
PPR binding.

**Tier B — predicted affinity. An annotation, and never a gate.**
Given a PPR specificity scoring table, score the PPR designed for target A against target
B. High relative score means the pair deserves a look; it does not mean they cross-react,
and a low score does not mean they are safe.

**Why tier B cannot be promoted to a criterion.** The available scoring table (Yan et al.,
as distributed with PPRmatcher) was derived from experiments on **P-type** PPR motifs. The
GRASP scaffold this package designs is **S-type**. The four code pairs GRASP uses do score
their cognate bases positively in that table, which is reassuring, but the applicability of
a P-type model to an S-type scaffold has not been established. Treating its output as a
guarantee would import exactly the kind of unvalidated assumption the rest of this package
refuses to make.

**The scoring table is not distributed with CLIPPR.** PPRmatcher carries no licence, so
redistributing its table would be a rights problem regardless of how useful it is. Obtain
it yourself and pass the path; without one, this module reports tier A alone and says so.

The useful question is not "is the model right" but **"does adding it change which library
you would build"** -- `compare_tiers` answers that directly.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .ppr import rna_to_code

#: Where a scoring table can be obtained. Not vendored; see the module docstring.
SCORE_TABLE_SOURCE = "https://github.com/ian-small/PPRmatcher (file: Yan.tsv)"


def load_ppr_scores(path: str | Path) -> dict[str, dict[str, float]]:
    """Read a PPR specificity table: code pair -> base -> score.

    Expects the PPRmatcher `Yan.tsv` layout -- a comment line, a header naming the bases,
    then one row per 5th/last residue pair. Any table in that shape works; nothing here is
    specific to that file beyond its format.
    """
    rows = [ln for ln in Path(path).read_text(encoding="utf-8").splitlines() if ln.strip()]
    header = next((i for i, ln in enumerate(rows) if ln.lower().startswith("5th")), None)
    if header is None:
        raise ValueError(f"{path}: no header row starting '5th/last' found")
    bases = [b.strip().upper().replace("U", "T") for b in rows[header].split("\t")[1:]]
    if set(bases) != set("ACGT"):
        raise ValueError(f"{path}: header names {bases}, expected the four bases")

    table: dict[str, dict[str, float]] = {}
    for line in rows[header + 1:]:
        parts = line.split("\t")
        if len(parts) != len(bases) + 1:
            continue
        code = parts[0].strip().upper()
        try:
            table[code] = {b: float(v) for b, v in zip(bases, parts[1:])}
        except ValueError:
            continue
    if not table:
        raise ValueError(f"{path}: no scoring rows parsed")
    return table


def predicted_affinity(designed_for: str, against: str,
                       scores: dict[str, dict[str, float]]) -> float:
    """Score the PPR designed for one target against another, relative to its own.

    1.0 means the PPR scores the other target exactly as well as its own; near 0 means far
    worse. A *relative* number, because the raw sum has no meaningful scale.

    This is a model-based annotation, not a binding prediction.
    """
    codes = rna_to_code(designed_for)
    other = str(against).strip().upper().replace("U", "T")
    own = str(designed_for).strip().upper().replace("U", "T")
    if len(other) != len(codes):
        raise ValueError(
            f"cannot score a {len(codes)}-repeat PPR against a {len(other)}-base target")

    missing = [c for c in codes if c not in scores]
    if missing:
        raise ValueError(f"scoring table has no entry for code pair(s) {sorted(set(missing))}")

    def total(seq: str) -> float:
        return sum(scores[c][b] for c, b in zip(codes, seq))

    best = total(own)
    if best == 0:
        return 0.0
    return total(other) / best


@dataclass(frozen=True)
class PairRisk:
    """One ordered pair of targets, with both tiers reported separately."""

    a: str
    b: str
    hamming: int
    affinity_a_to_b: float | None = None
    affinity_b_to_a: float | None = None

    @property
    def worst_affinity(self) -> float | None:
        vals = [v for v in (self.affinity_a_to_b, self.affinity_b_to_a) if v is not None]
        return max(vals) if vals else None

    @property
    def flag(self) -> str:
        """A word for the reader. Driven by tier A; tier B can only add 'review'."""
        if self.hamming <= 2:
            return "TOO CLOSE"
        worst = self.worst_affinity
        if worst is not None and worst >= 0.75:
            return "review"
        return "ok"


def assess(targets, scores=None) -> list[PairRisk]:
    """Every same-length pair, closest first, with tier B attached when available.

    Tier A here is plain, unweighted Hamming: a count of differing positions. The
    position weights in `orthogonal` belong to target *selection*, where a float score
    ranks candidate sets; a risk report wants the countable thing a reader can verify
    by eye.
    """
    targets = [str(t).strip().upper() for t in targets if str(t).strip()]
    out: list[PairRisk] = []
    for i, a in enumerate(targets):
        for b in targets[i + 1:]:
            if len(a) != len(b):
                continue                       # Hamming is undefined across lengths
            pair = PairRisk(a, b, sum(x != y for x, y in zip(a, b)))
            if scores:
                try:
                    pair = PairRisk(a, b, pair.hamming,
                                    predicted_affinity(a, b, scores),
                                    predicted_affinity(b, a, scores))
                except ValueError:
                    pass                        # unscoreable pair keeps tier A only
            out.append(pair)
    return sorted(out, key=lambda p: (p.hamming, -(p.worst_affinity or 0)))


def compare_tiers(targets, scores, min_hamming: int = 5,
                  affinity_review: float = 0.75) -> dict:
    """Would the predicted-affinity tier point you at pairs Hamming does not?

    The honest way to use a model whose applicability is unproven: not "is it right", but
    "does it change the answer". If the two tiers agree, the annotation costs nothing and
    adds confidence. If they disagree, the disagreeing pairs are a model-sensitive region
    that deserves a human look -- which is exactly what should be reported, rather than
    silently trusting or silently ignoring the model.

    `tiers_disagree` never means the gate moved. Tier A alone accepts and rejects; a
    disagreement is a reading recommendation for a person, nothing more.
    """
    pairs = assess(targets, scores)
    hamming_fail = {(p.a, p.b) for p in pairs if p.hamming < min_hamming}
    affinity_flag = {(p.a, p.b) for p in pairs
                     if p.worst_affinity is not None and p.worst_affinity >= affinity_review}
    return {
        "n_pairs": len(pairs),
        "hamming_flagged": sorted(hamming_fail),
        "affinity_flagged": sorted(affinity_flag),
        "agree": sorted(hamming_fail & affinity_flag),
        "affinity_only": sorted(affinity_flag - hamming_fail),
        "hamming_only": sorted(hamming_fail - affinity_flag),
        "tiers_disagree": bool(affinity_flag - hamming_fail),
    }


def report(targets, scores=None, min_hamming: int = 5, limit: int = 12) -> str:
    """A readable two-tier summary, with the caveat attached to the numbers."""
    pairs = assess(targets, scores)
    if not pairs:
        return "fewer than two targets of any one length; nothing to compare"

    has_b = any(p.worst_affinity is not None for p in pairs)
    head = f"{'target A':<22}{'target B':<22}{'hamming':>8}"
    head += f"{'A->B':>8}{'B->A':>8}  flag" if has_b else "  flag"
    lines = [head, "-" * (len(head) + 4)]
    for p in pairs[:limit]:
        row = f"{p.a:<22}{p.b:<22}{p.hamming:>8}"
        if has_b:
            row += (f"{p.affinity_a_to_b:>8.2f}{p.affinity_b_to_a:>8.2f}"
                    if p.affinity_a_to_b is not None else f"{'-':>8}{'-':>8}")
        lines.append(row + f"  {p.flag}")
    if len(pairs) > limit:
        lines.append(f"... {len(pairs) - limit} further pairs, all better separated")

    lines.append("")
    lines.append(f"Tier A, Hamming distance, is the hard criterion: pairs below "
                 f"{min_hamming} of {len(pairs[0].a)} positions are too close.")
    if has_b:
        lines.append("Tier B, predicted affinity, is a model-based annotation and gates "
                     "nothing. The scoring table is derived from P-type PPR experiments "
                     "while this scaffold is S-type, so its applicability here is "
                     "unestablished. 'review' means look, not fail.")
    else:
        lines.append(f"Tier B not shown: no scoring table supplied. Obtain one from "
                     f"{SCORE_TABLE_SOURCE} and pass it to add predicted-affinity "
                     f"annotation.")
    return "\n".join(lines)
