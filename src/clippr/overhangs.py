"""Golden Gate overhang fidelity, from the published mis-ligation matrices.

Source data: Pryor et al. (2020) PLoS ONE 15(9):e0238592, supplementary tables S1
(BsaI-HFv2) and S4 (BbsI-HF). Each is a 256x256 table of observed ligation counts;
columns are the 256 four-base overhangs in alphabetical order and rows are their reverse
complements.

Fidelity follows the published definition (Pryor et al., Materials and Methods): p(O) is
the fraction of O's ligations that pair it with its Watson-Crick partner rather than with
any other strand present, where N_total counts O's ligations "to any overhangs in the set
and its WC pair". Both strands of every junction are in the tube, so that competitor pool
is the set *together with* its reverse complements. The set fidelity is the product of
p(O), reported as the geometric mean of the two directional products so that it does not
depend on which strand you call "top".

Note this is a predicted surrogate, not a measured assembly efficiency.

**The aggregation convention, and how far it matters (settled 2026-09-14).**
Pryor's definition fixes the *denominator* -- competitors are the set together with its WC
pairs -- but not how the per-overhang probabilities are combined into one number. Two
combinations are equally faithful to it and equally orientation-invariant:

    A  geometric mean of the two per-strand directional products        <- implemented here
    B  one product of per-junction ratios, numerator and denominator
       each summed over both strands of the junction

They are not interchangeable in principle. On random valid four-member sets (seed 0) they
differ by at most 0.0038 (BsaI-HFv2) and 0.0034 (BbsI-HF), and they can *reverse a ranking*:
A scores ['GCAC','GTGA','TTGA','CATA'] above ['GGTT','AAAC','AGGA','ATCC'] (0.8341118 vs
0.8340968) while B scores it below (0.8315602 vs 0.8339764). A scalar gap that small is
therefore not an argument that the choice is safe, and an earlier draft of this module made
exactly that argument.

What licenses A is a direct measurement on the pools a design actually chooses from: running
`_plan_fragments` under both formulas across all 200 corpus targets, spanning 9S, 14S and 19S,
the plans `design_oneshot` consumes are **identical in 200 of 200 cases** -- same cuts, same
overhangs, same order. The reversals exist among arbitrary four-member sets; they do not reach
the reachable candidate space. So the defensible claim is "this choice changes no design CLIPPR
produces on this corpus", never "the choice is immaterial". If the candidate space changes --
a new destination, a different enzyme profile, another architecture -- that measurement is the
one to repeat.
"""
from __future__ import annotations

import hashlib
import itertools
import math
from collections.abc import Sequence
from functools import lru_cache
from pathlib import Path

import numpy as np

#: The matrices ship *inside* the package, not beside it. Resolving them relative to the
#: repository root works when running from a checkout and fails after `pip install`, where
#: there is no repository -- the wheel would install the code without the data it needs.
DATA = Path(__file__).resolve().parent / "data" / "matrices"
MATRICES = {"BsaI-HFv2": "BsaI-HFv2.xlsx", "BbsI-HF": "BbsI-HF.xlsx"}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]


@lru_cache(maxsize=8)
def load_matrix(name: str = "BsaI-HFv2") -> tuple[np.ndarray, list[str]]:
    """Parse a Pryor supplementary table to a 256x256 count array plus its labels.

    `counts[i, j]` is the observed number of ligations between top-strand overhang
    `labels[j]` and bottom-strand overhang `reverse_complement(labels[i])`; the table's
    row labels are the reverse complements of its column labels, in the same order,
    which this function checks rather than assumes.

    Caches a .npz beside the source so repeat loads do not re-parse the workbook.
    """
    if name not in MATRICES:
        raise ValueError(f"unknown matrix {name!r}; have {sorted(MATRICES)}")
    # The tag versions the cache layout; bump it whenever the stored arrays change, so a
    # stale .npz from an older layout is regenerated rather than misread.
    cache = DATA / f"{name}.v2.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        return z["counts"], [str(x) for x in z["labels"]]

    import pandas as pd
    src = DATA / MATRICES[name]
    if not src.exists():
        raise FileNotFoundError(
            f"{src} missing. Download from the Pryor et al. 2020 supplement:\n"
            f"  https://journals.plos.org/plosone/article/file?"
            f"id=10.1371/journal.pone.0238592.s001&type=supplementary  (BsaI-HFv2)\n"
            f"  ...s004... (BbsI-HF)"
        )
    raw = pd.read_excel(src, header=None)
    labels = [str(x).strip().upper() for x in raw.iloc[0, 1:]]
    rows = [str(x).strip().upper() for x in raw.iloc[1:, 0]]
    counts = raw.iloc[1:, 1:].to_numpy(dtype=np.int64)
    if counts.shape != (256, 256):
        raise ValueError(f"{name}: expected a 256x256 table, got {counts.shape}")
    if rows != [reverse_complement(c) for c in labels]:
        raise ValueError(f"{name}: row labels are not the reverse complements of the columns")
    try:
        np.savez_compressed(cache, counts=counts, labels=np.array(labels))
    except OSError:
        # An installed package often sits in a read-only site-packages. The cache only
        # saves about a second of parsing, so losing it is not worth failing a design over.
        pass
    return counts, labels


@lru_cache(maxsize=8)
def _index(matrix: str) -> tuple[np.ndarray, dict[str, int]]:
    counts, labels = load_matrix(matrix)
    return counts, {o: i for i, o in enumerate(labels)}


def _directional(overhangs: list[str], counts, idx) -> float:
    """Product over the set of p(O) = N_correct / N_total, each O read as top strand.

    A junction puts *both* of its strands in the tube, so the strands competing for O are
    the set members and their reverse complements together -- Pryor's "any overhangs in
    the set and its WC pair". Summing over the reverse complements alone made the score
    depend on which strand each junction happened to be named by: reverse-complementing
    one member of ["CTCA", "CTCG"] moved it from 0.830 to 1.000.

    Row i of the table is labelled `reverse_complement(labels[i])`, so bottom strand
    `reverse_complement(P)` is row `idx[P]`, bottom strand `P` is row
    `idx[reverse_complement(P)]`, and O's correct pairing sits on the diagonal.
    """
    product = 1.0
    for o in overhangs:
        correct = counts[idx[o], idx[o]]
        total = sum(counts[idx[p], idx[o]] + counts[idx[reverse_complement(p)], idx[o]]
                    for p in overhangs)
        if total == 0:
            return 0.0
        product *= correct / total
    return product


def fidelity_components(overhangs, matrix: str = "BsaI-HFv2") -> tuple[float, float]:
    """The forward and reverse directional products, before pooling.

    The count table is exactly symmetric, but these two still differ: p(O) and
    p(reverse_complement(O)) share a numerator and not a denominator, because the two
    strands of a junction have different affinities for the rest of the pool. Exposed so
    a result can be cross-checked against tools that report a single direction, such as
    the NEBridge Ligase Fidelity Viewer.
    """
    ohs = _clean(overhangs)
    counts, idx = _index(matrix)
    return (_directional(ohs, counts, idx),
            _directional([reverse_complement(o) for o in ohs], counts, idx))


def _clean(overhangs) -> list[str]:
    ohs = [str(o).strip().upper().replace("U", "T") for o in overhangs]
    if not ohs:
        raise ValueError("empty overhang set")
    for o in ohs:
        if len(o) != 4 or set(o) - set("ACGT"):
            raise ValueError(f"not a four-base DNA overhang: {o!r}")
    return ohs


def set_fidelity(overhangs, matrix: str = "BsaI-HFv2") -> float:
    """Predicted fidelity of a set of overhangs, in [0, 1]. Higher is better.

    Orientation-invariant: the geometric mean of the forward and reverse directional
    products, so reverse-complementing any member -- or the whole set -- leaves the
    score unchanged. Which strand a design calls "top" is arbitrary, and a fidelity
    score that moved when that label flipped would not be measuring the reaction.
    """
    fwd, rev = fidelity_components(overhangs, matrix)
    return math.sqrt(fwd * rev)


#: Bump whenever the *formula* in `_directional`, `fidelity_components` or `set_fidelity`
#: changes. It is part of the planning cache key, so a stale entry cannot outlive a scorer
#: change and resurrect the overhang selection the old formula preferred. Version 2 is the
#: pooled-competitor denominator adopted 2026-09-11; version 1 summed competitors over the
#: reverse complements of the set alone and chose different junctions on 5 of 5 targets.
SCORER_VERSION = 2


@lru_cache(maxsize=8)
def matrix_fingerprint(matrix: str = "BsaI-HFv2") -> str:
    """Digest of the ligation table's contents, for cache keys that must notice a data swap.

    Hashes the counts and their labels rather than the file, so a re-exported workbook with
    identical numbers is correctly treated as the same input.
    """
    counts, labels = load_matrix(matrix)
    h = hashlib.sha256(np.ascontiguousarray(counts).tobytes())
    h.update("".join(labels).encode())
    return h.hexdigest()[:16]


#: Results of `best_set`, keyed by everything that can change them. Per-process and
#: deliberately unbounded: the reachable key space is small -- 200 corpus proteins collapse
#: to a handful of distinct candidate pools -- and an eviction policy would be complexity
#: bought for no measured need.
_BEST_SET_CACHE: dict[tuple, tuple[tuple[str, ...], float]] = {}
_CACHE_STATS = {"hits": 0, "misses": 0}


def cache_stats() -> dict:
    """Hits, misses and size for the `best_set` cache, for reporting reuse."""
    return {**_CACHE_STATS, "entries": len(_BEST_SET_CACHE)}


def clear_cache() -> None:
    """Empty the `best_set` cache. Used to measure cold cost and to prove equivalence."""
    _BEST_SET_CACHE.clear()
    _CACHE_STATS.update(hits=0, misses=0)


def palindromic(overhang: str) -> bool:
    """A palindromic overhang ligates to itself and cannot be used."""
    o = overhang.upper().replace("U", "T")
    return o == reverse_complement(o)


def valid_set(overhangs) -> bool:
    """No palindromes, no repeats, no member that is another's reverse complement."""
    ohs = [str(o).upper().replace("U", "T") for o in overhangs]
    if len(set(ohs)) != len(ohs):
        return False
    if any(palindromic(o) for o in ohs):
        return False
    seen = set()
    for o in ohs:
        if reverse_complement(o) in seen:
            return False
        seen.add(o)
    return True


def reaction_overhangs(junctions: Sequence[str], level: str = "level0") -> list[str]:
    """The physical single-strand overhangs competing in one assembly reaction.

    `DESTINATION_OVERHANGS` stores each destination as a pair of *coding sites*, and the
    3' coding site is the reverse complement of the overhang the enzyme actually leaves.
    This is the one place that knows that, so callers get the strands that are really in
    the tube instead of rederiving the flip at each site.

    Retracted 2026-09-11: this helper was previously documented as guarding a "false
    all-clear", on the grounds that scoring level 0's stored ("CTCA", "CGAG") verbatim
    gave 1.000 while the present strands CTCA/CTCG gave 0.794. That gap was a defect in
    `set_fidelity`, not a hazard in the convention -- CGAG is the reverse complement of
    CTCG, so the two spellings name the same two junctions and an orientation-invariant
    scorer must and now does return 0.794 for both. The level 0 measurement itself is
    unaffected.
    """
    from . import constants as C

    try:
        five, three_coding = C.DESTINATION_OVERHANGS[level]
    except KeyError as e:
        raise ValueError(
            f"unknown level {level!r}; have {sorted(C.DESTINATION_OVERHANGS)}"
        ) from e
    return [five, *junctions, reverse_complement(three_coding)]


def enumerate_candidates(cut_contexts: Sequence[str]) -> list[list[str]]:
    """Per junction, the usable overhangs among those its cut can reach.

    `cut_contexts[i]` is a comma-separated list of the overhangs achievable at junction
    i, as produced by `arelf.achievable_overhangs`. Palindromic overhangs are dropped:
    a palindrome is its own reverse complement, so it ligates to itself and no choice of
    partners can rescue it.
    """
    out = []
    for i, ctx in enumerate(cut_contexts):
        options = [o.strip().upper().replace("U", "T")
                   for o in (ctx.split(",") if isinstance(ctx, str) else ctx) if o.strip()]
        usable = [o for o in dict.fromkeys(options) if not palindromic(o)]
        if not usable:
            raise ValueError(f"junction {i}: every achievable overhang is palindromic")
        out.append(usable)
    return out


def best_set(
    candidates: Sequence[Sequence[str]],
    fixed: Sequence[str] | None = None,
    matrix: str = "BsaI-HFv2",
    max_combinations: int = 500_000,
) -> tuple[list[str], float]:
    """Choose one overhang per junction, maximising fidelity of the whole reaction.

    `candidates[i]` is the list of overhangs achievable at junction i; `fixed` holds
    destination overhangs that are not ours to choose but still compete in the reaction.

    Exhaustive. Our measurements on GRASP architectures show the reachable space is small
    and the optimum is near 1.0, so enumeration is both simpler than a stochastic search
    and provably optimal. Raises if the space exceeds `max_combinations` rather than
    silently degrading to a heuristic.
    """
    fixed = list(fixed or [])
    # Exhaustive enumeration over the candidate product is the whole cost here: measured at
    # 19S, three splits took 12.0 s against 0.005 s to build the pools they search. The key
    # covers every input that can change the answer -- the ordered pools, the fixed
    # destination overhangs, the table's contents and the scoring formula's version.
    key = (tuple(tuple(c) for c in candidates), tuple(fixed), matrix, max_combinations,
           SCORER_VERSION, matrix_fingerprint(matrix))
    cached = _BEST_SET_CACHE.get(key)
    if cached is not None:
        _CACHE_STATS["hits"] += 1
        # A fresh list every time: the caller owns what it gets back and may sort or mutate
        # it, and a shared list would let one caller corrupt every later hit.
        return list(cached[0]), cached[1]
    _CACHE_STATS["misses"] += 1

    total = 1
    for c in candidates:
        total *= max(len(c), 1)
    if total > max_combinations:
        raise ValueError(
            f"{total:,} combinations exceeds max_combinations={max_combinations:,}; "
            "narrow the candidate lists rather than switching to a heuristic"
        )
    best, best_score = None, -1.0
    for combo in itertools.product(*candidates):
        full = list(combo) + fixed
        if not valid_set(full):
            continue
        score = set_fidelity(full, matrix)
        if score > best_score:
            best, best_score = list(combo), score
            if best_score >= 1.0:
                # 1.0 is the maximum the metric can take, so nothing later can beat it.
                # Exhaustive search is still exhaustive: this returns a proven optimum,
                # it only stops enumerating combinations that cannot improve on it.
                break
    if best is None:
        raise ValueError("no valid overhang combination; all candidates conflict")
    _BEST_SET_CACHE[key] = (tuple(best), best_score)
    return best, best_score
