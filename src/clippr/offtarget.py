"""Does a designed target also occur somewhere it should not?

A PPR does not know which RNA you meant. If the sequence you designed it against also
appears in an endogenous chloroplast transcript, the protein will bind there too, and the
regulator stops being specific to your construct.

`orthogonal.py` asks whether your targets differ from *each other*. This asks the question
that module cannot: whether a target collides with **the host**.

**The result is uncomfortable and it is the point of this module.** The Chlamydomonas
chloroplast genome is 203,828 bp and 34.5% GC. A 9-nt target is expected to occur roughly
one to several times in it by chance alone, and AT-rich targets far more often than that,
because the genome itself is AT-rich. Length is what buys specificity: a 14-nt target is
expected essentially never, and a 19-nt target never. Measure before assuming 9S is safe
in vivo -- see `architecture_advice()`.

Nothing here predicts binding affinity. It reports sequence occurrence, which is a
necessary condition for an off-target interaction, not a sufficient one.
"""
from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path

CACHE = Path(__file__).resolve().parents[2] / "data" / "genomes"

#: Chlamydomonas reinhardtii chloroplast, RefSeq complete genome.
CHLOROPLAST = "NC_005353.1"

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Hit:
    """One place a target sequence occurs in the genome."""

    strand: str          # "+" or "-"
    position: int        # 0-based, on the plus strand
    mismatches: int
    context: str         # the genomic sequence with flanks, match in upper case

    def __str__(self) -> str:
        kind = "exact" if not self.mismatches else f"{self.mismatches} mismatch"
        return f"{self.strand}{self.position:>7}  {kind:12s} {self.context}"


def load_genome(accession: str = CHLOROPLAST) -> str:
    """Fetch a genome from NCBI, caching it so a design run works offline afterwards."""
    CACHE.mkdir(parents=True, exist_ok=True)
    cached = CACHE / f"{accession}.txt"
    if cached.exists():
        return cached.read_text(encoding="utf-8").strip()

    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={accession}&rettype=fasta&retmode=text")
    req = urllib.request.Request(url, headers={"User-Agent": "clippr"})
    text = urllib.request.urlopen(req, timeout=60).read().decode()
    seq = "".join(text.splitlines()[1:]).upper()
    if not seq or set(seq) - set("ACGTN"):
        raise ValueError(f"{accession}: fetched sequence is not plain nucleotides")
    cached.write_text(seq, encoding="utf-8")
    return seq


def _normalise(target: str) -> str:
    t = str(target).strip().upper().replace("U", "T")
    if not t or set(t) - set("ACGT"):
        raise ValueError(f"not a plain nucleotide sequence: {target!r}")
    return t


def find_hits(target: str, genome: str, max_mismatches: int = 0,
              flank: int = 6) -> list[Hit]:
    """Every occurrence of `target` in either strand, allowing up to `max_mismatches`.

    Both strands are searched because a chloroplast transcript may come from either, and
    a PPR binds the RNA that is actually transcribed.
    """
    t = _normalise(target)
    n, g = len(t), genome
    hits: list[Hit] = []

    def _ctx(i: int) -> str:
        lo, hi = max(0, i - flank), min(len(g), i + n + flank)
        return g[lo:i].lower() + g[i:i + n] + g[i + n:hi].lower()

    for strand, probe in (("+", t), ("-", reverse_complement(t))):
        if max_mismatches == 0:
            # str.find runs in C; the explicit position loop below is ~100x slower and
            # exact matching is the common case, so it gets its own path.
            i = g.find(probe)
            while i != -1:
                hits.append(Hit(strand, i, 0, _ctx(i)))
                i = g.find(probe, i + 1)
            continue

        for i in range(len(g) - n + 1):
            mm = 0
            for a, b in zip(g[i:i + n], probe):
                if a != b:
                    mm += 1
                    if mm > max_mismatches:
                        break
            else:
                hits.append(Hit(strand, i, mm, _ctx(i)))
    return hits


def expected_by_chance(target: str, genome: str) -> float:
    """How often a sequence like this would occur by chance, given the genome's own bases.

    Uses the genome's mononucleotide composition rather than assuming equal bases. That
    matters here: the chloroplast genome is 34.5% GC, so an AT-rich target occurs far more
    often than a uniform model predicts, and a uniform model would understate the risk for
    exactly the targets most likely to be chosen.
    """
    t = _normalise(target)
    total = len(genome)
    freq = {b: genome.count(b) / total for b in "ACGT"}
    p = 1.0
    for base in t:
        p *= freq.get(base, 0.0)
    positions = (total - len(t) + 1) * 2      # both strands
    return p * positions


def scan(target: str, genome: str | None = None, max_mismatches: int = 1) -> dict:
    """Assess one target against the host genome.

    Returns the exact hits, the near hits, what chance alone would predict, and a verdict
    that is deliberately conservative: any exact occurrence is a concern worth reading,
    not a pass/fail the caller can ignore.
    """
    g = load_genome() if genome is None else genome
    t = _normalise(target)
    exact = find_hits(t, g, 0)
    near = [h for h in find_hits(t, g, max_mismatches) if h.mismatches]
    expected = expected_by_chance(t, g)

    if exact:
        verdict = "OCCURS IN HOST"
    elif near:
        verdict = "NEAR MATCH IN HOST"
    else:
        verdict = "not found in host"

    return {
        "target": target,
        "length": len(t),
        "verdict": verdict,
        "exact_hits": exact,
        "near_hits": near,
        "n_exact": len(exact),
        "n_near": len(near),
        "expected_by_chance": expected,
        "genome_length": len(g),
    }


def architecture_advice(genome: str | None = None) -> dict[int, float]:
    """Expected chance occurrences of a target of each architecture length.

    The number that decides whether a 9-nt target can be specific in this host at all.
    Computed against the genome's real base composition, for an average target; an AT-rich
    target will be worse and a GC-rich one better.
    """
    g = load_genome() if genome is None else genome
    total = len(g)
    freq = {b: g.count(b) / total for b in "ACGT"}
    mean_p = sum(f * f for f in freq.values())      # per-position match probability
    positions = total * 2
    return {n: mean_p ** n * positions for n in (9, 14, 19)}


def report(results: list[dict]) -> str:
    """A readable summary for a set of scanned targets."""
    lines = [f"{'target':<22}{'len':>4}{'exact':>7}{'near':>6}{'expected':>10}  verdict",
             "-" * 78]
    for r in results:
        lines.append(f"{r['target']:<22}{r['length']:>4}{r['n_exact']:>7}"
                     f"{r['n_near']:>6}{r['expected_by_chance']:>10.2f}  {r['verdict']}")
    flagged = [r for r in results if r["n_exact"]]
    if flagged:
        lines.append("")
        lines.append(f"{len(flagged)} target(s) occur in the host genome:")
        for r in flagged:
            for h in r["exact_hits"][:3]:
                lines.append(f"  {r['target']}  {h}")
    return "\n".join(lines)
