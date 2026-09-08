"""Synthesis QC for a coding sequence, with thresholds calibrated to actually separate.

The baseline pipeline flagged **WARNING on 200 of 200** designs with zero failures. A
flag that always fires carries no information: it cannot rank designs, cannot gate an
order, and trains whoever reads it to ignore it. This module exists to fix that, so the
thresholds below are chosen against measured distributions and every one is documented
with where it sits.

**Repeats are the axis that discriminates.** Measured over 400 sequences -- the 200
reference CDSs and the 200 our own optimiser produces for the same proteins:

    metric              reference (min/med/max)      ours (min/med/max)
    repeat_fraction     0.04 / 0.44 / 0.76           0.00 / 0.00 / 0.02
    longest_repeat      20   / 38   / 62             17   / 19   / 20
    longest_homopolymer  4   /  4   /  5              4   /  4   /  4
    gc_pct              39.1 / 46.0 / 48.5           58.8 / 59.8 / 60.6

Under the defaults here the baseline population splits **14 PASS / 86 WARNING /
100 FAIL** -- a real distribution rather than a constant.

**Two honesty notes.**

*Our own sequences all pass, and that is partly circular.* `codons.optimize_cds` defaults
to making every 20-mer unique, which caps `longest_repeat` at 19-20 by construction, so it
is aimed at exactly what this module measures. The margins are not tight -- WARNING sits
at 25 nt and FAIL at 50 -- and uniquifying is only an objective, so it guarantees nothing.
But 200/200 PASS should be read as the two components agreeing, not as independent
confirmation.

*GC is deliberately not calibrated from this corpus.* The two populations sit at 39-48%
and 59-61% because one targets E. coli and the other Chlamydomonas. Neither is defective;
that is the organism. Cut-points drawn from the combined distribution would flag every
correct Chlamydomonas sequence. The GC and homopolymer gates here are synthesis-feasibility
limits instead, and they fire on neither population -- which is the right behaviour for a
gate whose job is to catch sequences no vendor will make, not to rank ordinary ones.
"""
from __future__ import annotations

from collections.abc import Mapping

from . import constants as _C

#: Enzyme sites that must not appear inside a fragment, by the default profile.
#: `constants.ENZYME_PROFILES` records why each is excluded -- assembly chemistry for
#: BsaI/BbsI, the iGEM RFC[1000] standard for SapI, downstream MoClo headroom for BsmBI.
BLACKLIST_ENZYMES: tuple[str, ...] = _C.enzymes_for()

#: Default cut-points. `repeat_*` discriminate; the rest are feasibility floors.
#:
#: Percentiles are over the combined 400-sequence distribution described in the module
#: docstring. Where a gate is marked "inert here", it did not fire on either population --
#: it guards a real synthesis failure this corpus happens not to produce.
DEFAULT_THRESHOLDS: dict[str, float] = {
    # discriminating -- calibrated
    "longest_repeat_warn": 25,      # p63; direct repeats past ~25 nt trouble most vendors
    "longest_repeat_fail": 50,      # p91
    "repeat_fraction_warn": 0.10,   # p54
    "repeat_fraction_fail": 0.50,   # p75
    # feasibility floors -- inert on this corpus, by design
    "gc_pct_warn_low": 35.0,
    "gc_pct_warn_high": 65.0,
    "gc_pct_fail_low": 25.0,
    "gc_pct_fail_high": 75.0,
    "gc_window_warn_low": 25.0,
    "gc_window_warn_high": 75.0,
    "gc_window_fail_low": 15.0,
    "gc_window_fail_high": 85.0,
    "homopolymer_warn": 6,
    "homopolymer_fail": 9,
}

REPEAT_K: int = 20
GC_WINDOW: int = 50
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _homopolymers(cds: str, minimum: int = 5) -> tuple[int, int]:
    """(longest run of one base, number of runs at least `minimum` long)."""
    longest = run = 1
    count = 0
    for i, b in enumerate(cds):
        if i and b == cds[i - 1]:
            run += 1
        else:
            if run >= minimum:
                count += 1
            run = 1
        longest = max(longest, run)
    if run >= minimum:
        count += 1
    return longest, count


def _repeats(cds: str, k: int = REPEAT_K) -> tuple[float, int, int]:
    """(fraction of positions inside a repeated k-mer, distinct repeated k-mers, longest).

    `repeat_fraction` is the discriminating measure: it counts how much of the sequence is
    covered by duplicated k-mers rather than only the single worst stretch, so a sequence
    repetitive throughout scores worse than one with a single long duplication.

    **It is a sequence-complexity score, not a synthesis-failure probability.** k-mer
    uniqueness is a recognised sequence-design constraint -- DNA Chisel ships
    `UniquifyAllKmers` for it -- but nothing here is calibrated against vendor outcomes,
    and no design in this project has been synthesised. Read it as "how repetitive is
    this", never as "how likely is this to fail".
    """
    n = len(cds)
    if n < k:
        return 0.0, 0, 0
    positions: dict[str, list[int]] = {}
    for i in range(n - k + 1):
        positions.setdefault(cds[i:i + k], []).append(i)

    covered: set[int] = set()
    duplicated = 0
    for pos in positions.values():
        if len(pos) > 1:
            duplicated += 1
            for p in pos:
                covered.update(range(p, p + k))

    longest = 0
    for size in range(k // 2, n):
        seen: set[str] = set()
        hit = False
        for i in range(n - size + 1):
            frag = cds[i:i + size]
            if frag in seen:
                hit = True
                break
            seen.add(frag)
        if hit:
            longest = size
        else:
            break
    return len(covered) / n, duplicated, longest


def _blacklist(cds: str, enzymes=BLACKLIST_ENZYMES) -> list[str]:
    """Enzyme sites present on either strand."""
    from Bio import Restriction

    hits = []
    for name in enzymes:
        site = str(getattr(Restriction, name).site).upper()
        rc = site.translate(_COMPLEMENT)[::-1]
        n = cds.count(site) + (cds.count(rc) if rc != site else 0)
        if n:
            hits.append(f"{name}x{n}")
    return hits


def synthesis_qc(cds: str, thresholds: Mapping[str, float] | None = None) -> dict:
    """Assess whether a coding sequence is likely to synthesise, and why not.

    Returns a verdict of PASS, WARNING or FAIL alongside every measurement behind it, so
    a caller can rank designs on the numbers rather than only on the label. `thresholds`
    overrides any subset of `DEFAULT_THRESHOLDS`.
    """
    cds = str(cds).upper().replace("U", "T")
    if not cds:
        raise ValueError("empty coding sequence")
    if set(cds) - set("ACGT"):
        raise ValueError(f"non-ACGT bases in sequence: {sorted(set(cds) - set('ACGT'))}")

    t = {**DEFAULT_THRESHOLDS, **(thresholds or {})}
    n = len(cds)

    gc_pct = 100 * sum(1 for b in cds if b in "GC") / n
    windows = [100 * sum(1 for b in cds[i:i + GC_WINDOW] if b in "GC") / GC_WINDOW
               for i in range(max(1, n - GC_WINDOW + 1))] or [gc_pct]
    longest_homopolymer, homopolymer_count = _homopolymers(cds)
    repeat_fraction, duplicated, longest_repeat = _repeats(cds)
    blacklist_hits = _blacklist(cds)

    warnings: list[str] = []
    failures: list[str] = []

    if blacklist_hits:
        failures.append(f"blacklisted enzyme site: {', '.join(blacklist_hits)}")

    if longest_repeat >= t["longest_repeat_fail"]:
        failures.append(f"longest exact repeat {longest_repeat} nt "
                        f">= {t['longest_repeat_fail']:.0f}")
    elif longest_repeat >= t["longest_repeat_warn"]:
        warnings.append(f"longest exact repeat {longest_repeat} nt "
                        f">= {t['longest_repeat_warn']:.0f}")

    if repeat_fraction >= t["repeat_fraction_fail"]:
        failures.append(f"{repeat_fraction:.0%} of the sequence lies in a repeated "
                        f"{REPEAT_K}-mer")
    elif repeat_fraction >= t["repeat_fraction_warn"]:
        warnings.append(f"{repeat_fraction:.0%} of the sequence lies in a repeated "
                        f"{REPEAT_K}-mer")

    if not t["gc_pct_fail_low"] <= gc_pct <= t["gc_pct_fail_high"]:
        failures.append(f"GC {gc_pct:.1f}% outside "
                        f"{t['gc_pct_fail_low']:.0f}-{t['gc_pct_fail_high']:.0f}%")
    elif not t["gc_pct_warn_low"] <= gc_pct <= t["gc_pct_warn_high"]:
        warnings.append(f"GC {gc_pct:.1f}% outside "
                        f"{t['gc_pct_warn_low']:.0f}-{t['gc_pct_warn_high']:.0f}%")

    lo, hi = min(windows), max(windows)
    if lo < t["gc_window_fail_low"] or hi > t["gc_window_fail_high"]:
        failures.append(f"local GC {lo:.0f}-{hi:.0f}% in a {GC_WINDOW}-nt window")
    elif lo < t["gc_window_warn_low"] or hi > t["gc_window_warn_high"]:
        warnings.append(f"local GC {lo:.0f}-{hi:.0f}% in a {GC_WINDOW}-nt window")

    if longest_homopolymer >= t["homopolymer_fail"]:
        failures.append(f"homopolymer run of {longest_homopolymer}")
    elif longest_homopolymer >= t["homopolymer_warn"]:
        warnings.append(f"homopolymer run of {longest_homopolymer}")

    return {
        "status": "FAIL" if failures else ("WARNING" if warnings else "PASS"),
        "length_nt": n,
        "gc_pct": round(gc_pct, 2),
        "gc_window_min": round(lo, 2),
        "gc_window_max": round(hi, 2),
        "longest_homopolymer": longest_homopolymer,
        "homopolymer_count": homopolymer_count,
        "repeated_kmer_fraction": round(repeat_fraction, 4),
        "repeated_kmers": duplicated,
        "longest_repeat": longest_repeat,
        "blacklist_hits": blacklist_hits,
        "warnings": warnings,
        "failures": failures,
    }
