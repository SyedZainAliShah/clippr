"""The measurable objectives, implemented exactly as `docs/objectives.md` defines them.

Each function returns the value *and* what it was computed over, because two objective values
taken over different intervals or different *k* are not comparable and a bare float hides
that. Reporting the scope alongside the number is what makes a later comparison checkable.

Nothing here decides feasibility. Hard constraints live apart from objectives, and a candidate
that violates one is infeasible no matter how it scores -- see `docs/objectives.md` §1.
"""
from __future__ import annotations

import math
from collections import Counter

_COMPLEMENT = str.maketrans("ACGT", "TGCA")

#: One codon each, so their relative weight is 1 by construction and no optimiser can move
#: them. Including them would dilute every score with a constant.
_SINGLE_CODON_AA = ("M", "W")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def relative_adaptiveness(codon_table: dict) -> dict[str, float]:
    """w(c) = f(c) / max f over synonyms, keyed by DNA codon."""
    weights: dict[str, float] = {}
    for _aa, codons in codon_table.items():
        if not codons:
            continue
        top = max(codons.values())
        if top <= 0:
            continue
        for codon, value in codons.items():
            weights[codon.replace("U", "T")] = value / top
    return weights


def codon_adaptation(dna: str, codon_table: dict, genetic_code: int = 1) -> dict:
    """Geometric mean of relative synonymous weights over `dna`.

    The awkward cases are decided here rather than left to chance, because independent
    implementations differ on them and a silent difference looks like a disagreement about
    the sequence:

      * methionine and tryptophan are excluded -- one codon each, weight 1 by construction
      * stop codons are excluded
      * a codon absent from the table is excluded and counted in `missing_weights`, never
        treated as weight 1
      * a zero weight gives CAI 0 rather than `-inf`, and the offending codons are reported;
        one zero must not make two otherwise different sequences incomparable
    """
    from Bio.Data import CodonTable

    forward = CodonTable.unambiguous_dna_by_id[genetic_code].forward_table
    weights = relative_adaptiveness(codon_table)

    used, missing, zeros = [], 0, []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa = forward.get(codon)
        if aa is None or aa in _SINGLE_CODON_AA:
            continue
        if codon not in weights:
            missing += 1
            continue
        if weights[codon] <= 0:
            zeros.append(codon)
            continue
        used.append(weights[codon])

    if zeros:
        cai = 0.0
    elif used:
        cai = math.exp(sum(math.log(w) for w in used) / len(used))
    else:
        cai = 0.0
    return {"cai": cai, "codons": len(used), "nt": len(dna),
            "missing_weights": missing, "zero_weight_codons": sorted(set(zeros))}


def duplicated_kmers(seq: str, k: int = 20) -> int:
    """Occurrences beyond the first, summed over every distinct k-mer."""
    if len(seq) < k:
        return 0
    windows = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    return len(windows) - len(set(windows))


def distinct_repeated_words(seq: str, k: int = 20) -> int:
    """How many distinct k-mers occur more than once. Not the same as `duplicated_kmers`."""
    if len(seq) < k:
        return 0
    counts = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    return sum(1 for n in counts.values() if n > 1)


def collection_penalty(records: dict[str, str], k: int = 20,
                       include_reverse_complement: bool = False) -> dict:
    """Shared k-mer mass between distinct records: sum over pairs of sum_w min(c_i, c_j).

    Measured over **distinct physical records only**. Reusing one selected record across
    several target products is what a reusable inventory is for; counting that reuse as
    sharing would penalise the feature for working as designed.

    This is an engineering objective for steering a search. It is **not** a recombination
    probability and carries no biological interpretation.
    """
    ids = sorted(records)
    counts: dict[str, Counter] = {}
    for module_id in ids:
        seq = records[module_id]
        words = [seq[i:i + k] for i in range(max(0, len(seq) - k + 1))]
        if include_reverse_complement:
            rc = reverse_complement(seq)
            words += [rc[i:i + k] for i in range(max(0, len(rc) - k + 1))]
        counts[module_id] = Counter(words)

    penalty = 0
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            left, right = counts[ids[a]], counts[ids[b]]
            smaller, larger = (left, right) if len(left) <= len(right) else (right, left)
            penalty += sum(min(n, larger[w]) for w, n in smaller.items() if w in larger)
    return {"penalty": penalty, "k": k, "records": len(ids),
            "include_reverse_complement": include_reverse_complement}


def worst_shared_tract(records: dict[str, str]) -> dict:
    """Longest exact substring common to any two distinct records, and which pair.

    Reported independently of `collection_penalty`. Lowering the penalty does not necessarily
    lower this: the penalty aggregates many short coincidences, this is a single extreme, and
    a search can improve one while leaving the other untouched.
    """
    from .homology import longest_shared

    ids = sorted(records)
    worst, pair = 0, None
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            length = longest_shared(records[ids[a]], records[ids[b]])[0]
            if length > worst:
                worst, pair = length, (ids[a], ids[b])
    return {"nt": worst, "pair": pair, "records": len(ids)}


def dominates(a: dict, b: dict, directions: dict[str, str],
              tolerance: float = 1e-9) -> bool:
    """True when `a` dominates `b`: no worse everywhere, strictly better somewhere.

    `directions` maps each objective to "max" or "min". The tolerance keeps a floating-point
    wobble from manufacturing dominance -- differences within it count as equal.
    """
    strictly_better = False
    for name, direction in directions.items():
        x, y = a[name], b[name]
        if direction == "max":
            if x < y - tolerance:
                return False
            if x > y + tolerance:
                strictly_better = True
        elif direction == "min":
            if x > y + tolerance:
                return False
            if x < y - tolerance:
                strictly_better = True
        else:
            raise ValueError(f"objective {name!r} has direction {direction!r}, "
                             f"expected 'max' or 'min'")
    return strictly_better


def nondominated(candidates: list[dict], directions: dict[str, str],
                 tolerance: float = 1e-9) -> list[int]:
    """Indices of the observed nondominated set.

    Observed, not globally optimal: nondominated among the candidates actually evaluated. It
    is the Pareto front only if the feasible space was exhaustively enumerated, which for a
    real inventory it is not.
    """
    keep = []
    for i, candidate in enumerate(candidates):
        if not any(dominates(other, candidate, directions, tolerance)
                   for j, other in enumerate(candidates) if i != j):
            keep.append(i)
    return keep
