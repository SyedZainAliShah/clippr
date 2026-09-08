"""Golden Gate overhang fidelity, from the published mis-ligation matrices.

Source data: Pryor et al. (2020) PLoS ONE 15(9):e0238592, supplementary tables S1
(BsaI-HFv2) and S4 (BbsI-HF). Each is a 256x256 table of observed ligation counts;
columns are the 256 four-base overhangs in alphabetical order and rows are their reverse
complements.

Fidelity follows the published definition: for a set of overhangs, p(O) is the fraction
of O's ligations that pair it with its Watson-Crick partner rather than with any other
member of the set, and the set fidelity is the product of p(O) over the set. The score
reported here is the orientation-invariant geometric mean of the two directional
products, so it does not depend on which strand you call "top".

Note this is a predicted surrogate, not a measured assembly efficiency.
"""
from __future__ import annotations

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

    Every top strand in the reaction can anneal to any bottom strand present, so the
    competitors for O are the reverse complements of all set members, O's own included.

    Row i of the table is labelled `reverse_complement(labels[i])`, so the row carrying
    bottom strand `reverse_complement(P)` is row `idx[P]`, and the correct pairing of O
    with its own complement sits on the diagonal.
    """
    product = 1.0
    for o in overhangs:
        correct = counts[idx[o], idx[o]]
        total = sum(counts[idx[p], idx[o]] for p in overhangs)
        if total == 0:
            return 0.0
        product *= correct / total
    return product


def fidelity_components(overhangs, matrix: str = "BsaI-HFv2") -> tuple[float, float]:
    """The forward and reverse directional products, before pooling.

    The assay's count table is not exactly symmetric, so these are two noisy estimates
    of the same quantity. Exposed so a result can be cross-checked against tools that
    report a single direction, such as the NEBridge Ligase Fidelity Viewer.
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

    Use this rather than assembling the set by hand. `DESTINATION_OVERHANGS` stores each
    destination as a pair of *coding sites*, and the 3' coding site is the reverse
    complement of the overhang the enzyme actually leaves. Scoring the stored pair
    verbatim compares two sequences that never meet in the tube: level 0's ("CTCA",
    "CGAG") scores a clean 1.000 that way, while the strands really present -- CTCA and
    CTCG -- score 0.794. The false reading is the more reassuring one, which is what
    makes it worth removing from the caller's hands.
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
    return best, best_score
