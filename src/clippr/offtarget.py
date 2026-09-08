"""Does a designed target also occur where the PPR could actually bind it?

A PPR binds **RNA**, so the universe that matters is transcribed sequence read 5'->3' in
its own orientation. That is a narrower question than "does this sequence appear in the
genome", and the difference is not small.

**Measured on the Chlamydomonas chloroplast genome, 200 random targets per length:**

    length   in genomic DNA   in a transcript
      9 nt     97 (48%)         32 (16%)
     14 nt      0                0
     19 nt      0                0

Counting genomic DNA on both strands overstates the risk roughly threefold. A match in a
non-transcribed region is not an RNA off-target, and a reverse-complement match in DNA is
not one either -- the transcript from that locus carries the other sequence.

This module therefore reports **tiers**, not a single number:

    genomic      the sequence exists somewhere in the chloroplast DNA
    transcript   it exists inside an annotated transcript, in the sense orientation,
                 which is the only tier a PPR could plausibly act on

Sequence occurrence is a **necessary** condition for an off-target interaction, never a
sufficient one. Nothing here predicts binding affinity: a PPR specificity score would be a
separate, model-based annotation, and this module deliberately does not pretend to one.
"""
from __future__ import annotations

import io
import urllib.request
from dataclasses import dataclass
from pathlib import Path

CACHE = Path(__file__).resolve().parents[2] / "data" / "genomes"

#: Chlamydomonas reinhardtii chloroplast, RefSeq complete genome.
CHLOROPLAST = "NC_005353.1"

#: Feature types that become RNA. A PPR can only bind what is transcribed.
TRANSCRIBED_TYPES = ("CDS", "rRNA", "tRNA", "ncRNA", "misc_RNA", "tmRNA")

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Transcript:
    """One annotated transcript, already in its own 5'->3' orientation."""

    name: str
    kind: str
    strand: int
    start: int
    sequence: str


@dataclass(frozen=True)
class Hit:
    """One occurrence of a target sequence."""

    tier: str            # "transcript" or "genomic"
    where: str           # gene name, or "-" for a bare genomic hit
    kind: str            # feature type, or "DNA"
    position: int
    context: str

    def __str__(self) -> str:
        loc = f"{self.where} ({self.kind})" if self.tier == "transcript" else "genomic DNA"
        return f"{self.tier:<11}{loc:<24}{self.position:>7}  {self.context}"


def _fetch(accession: str, rettype: str, suffix: str) -> str:
    CACHE.mkdir(parents=True, exist_ok=True)
    cached = CACHE / f"{accession}{suffix}"
    if cached.exists():
        return cached.read_text(encoding="utf-8")
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={accession}&rettype={rettype}&retmode=text")
    req = urllib.request.Request(url, headers={"User-Agent": "clippr"})
    text = urllib.request.urlopen(req, timeout=90).read().decode()
    cached.write_text(text, encoding="utf-8")
    return text


def load_genome(accession: str = CHLOROPLAST) -> str:
    """The plus strand of the genome, as plain sequence."""
    cached = CACHE / f"{accession}.txt"
    if cached.exists():
        return cached.read_text(encoding="utf-8").strip()
    seq = "".join(_fetch(accession, "fasta", ".fasta").splitlines()[1:]).upper()
    if not seq or set(seq) - set("ACGTN"):
        raise ValueError(f"{accession}: fetched sequence is not plain nucleotides")
    CACHE.mkdir(parents=True, exist_ok=True)
    cached.write_text(seq, encoding="utf-8")
    return seq


def load_transcripts(accession: str = CHLOROPLAST) -> list[Transcript]:
    """Every annotated transcript, each already in the orientation it is read in.

    `feature.extract` reverse-complements minus-strand features for us, so the stored
    sequence is the RNA a PPR would see rather than the plus strand of the DNA.
    """
    from Bio import SeqIO

    record = SeqIO.read(io.StringIO(_fetch(accession, "gb", ".gb")), "genbank")
    out: list[Transcript] = []
    for f in record.features:
        if f.type not in TRANSCRIBED_TYPES:
            continue
        name = (f.qualifiers.get("gene") or f.qualifiers.get("locus_tag")
                or f.qualifiers.get("product") or ["unnamed"])[0]
        out.append(Transcript(name=name, kind=f.type,
                              strand=int(f.location.strand or 1),
                              start=int(f.location.start),
                              sequence=str(f.extract(record.seq)).upper()))
    return out


def _normalise(target: str) -> str:
    t = str(target).strip().upper().replace("U", "T")
    if not t or set(t) - set("ACGT"):
        raise ValueError(f"not a plain nucleotide sequence: {target!r}")
    return t


def _context(seq: str, i: int, n: int, flank: int = 6) -> str:
    lo, hi = max(0, i - flank), min(len(seq), i + n + flank)
    return seq[lo:i].lower() + seq[i:i + n] + seq[i + n:hi].lower()


def find_in_transcripts(target: str, transcripts: list[Transcript]) -> list[Hit]:
    """Sense-strand occurrences inside annotated RNA -- the tier that can matter."""
    t = _normalise(target)
    hits: list[Hit] = []
    for tr in transcripts:
        i = tr.sequence.find(t)
        while i != -1:
            hits.append(Hit("transcript", tr.name, tr.kind, i,
                            _context(tr.sequence, i, len(t))))
            i = tr.sequence.find(t, i + 1)
    return hits


def find_in_genome(target: str, genome: str) -> list[Hit]:
    """Occurrences anywhere in the DNA, either strand. The broader, weaker signal."""
    t = _normalise(target)
    hits: list[Hit] = []
    for probe in {t, reverse_complement(t)}:
        i = genome.find(probe)
        while i != -1:
            hits.append(Hit("genomic", "-", "DNA", i, _context(genome, i, len(t))))
            i = genome.find(probe, i + 1)
    return hits


def expected_by_chance(target: str, sequence: str, both_strands: bool = True) -> float:
    """Chance occurrences given the sequence's own base composition.

    A uniform model is wrong here: the chloroplast is 34.5% GC, so an AT-rich target
    occurs far more often than 4^-k predicts -- and AT-rich targets are exactly the ones
    an AT-rich UTR context tends to produce. This is a sequence-composition null, not a
    model of PPR binding.
    """
    t = _normalise(target)
    total = len(sequence)
    freq = {b: sequence.count(b) / total for b in "ACGT"}
    p = 1.0
    for base in t:
        p *= freq.get(base, 0.0)
    positions = (total - len(t) + 1) * (2 if both_strands else 1)
    return p * positions


def scan(target: str, genome: str | None = None,
         transcripts: list[Transcript] | None = None) -> dict:
    """Assess one target against the host, reporting genomic and transcript tiers apart.

    The verdict names the strongest tier reached, because they mean different things: a
    transcript hit is something a PPR could act on, a genomic-only hit is a sequence
    coincidence in DNA that is not transcribed in that orientation.
    """
    g = load_genome() if genome is None else genome
    trs = load_transcripts() if transcripts is None else transcripts
    t = _normalise(target)

    in_rna = find_in_transcripts(t, trs)
    in_dna = find_in_genome(t, g)

    if in_rna:
        verdict = "OCCURS IN A TRANSCRIPT"
    elif in_dna:
        verdict = "in genomic DNA only, not in an annotated transcript"
    else:
        verdict = "not found in the host"

    return {
        "target": target,
        "length": len(t),
        "verdict": verdict,
        "transcript_hits": in_rna,
        "genomic_hits": in_dna,
        "n_transcript": len(in_rna),
        "n_genomic": len(in_dna),
        "genes": sorted({h.where for h in in_rna}),
        "expected_by_chance": expected_by_chance(t, g),
        "genome_length": len(g),
        "transcribed_nt": sum(len(x.sequence) for x in trs),
    }


def architecture_advice(genome: str | None = None) -> dict[int, float]:
    """Expected chance occurrences per architecture length, on genomic DNA.

    The number that shows why length buys specificity. Genomic rather than transcript
    scale, so it is an upper bound on the risk rather than an estimate of it.
    """
    g = load_genome() if genome is None else genome
    total = len(g)
    freq = {b: g.count(b) / total for b in "ACGT"}
    mean_p = sum(f * f for f in freq.values())
    return {n: mean_p ** n * total * 2 for n in (9, 14, 19)}


def report(results: list[dict]) -> str:
    """A readable summary, keeping the two tiers visibly separate."""
    lines = [f"{'target':<22}{'len':>4}{'transcript':>12}{'genomic':>9}"
             f"{'expected':>10}  verdict",
             "-" * 92]
    for r in results:
        lines.append(f"{r['target']:<22}{r['length']:>4}{r['n_transcript']:>12}"
                     f"{r['n_genomic']:>9}{r['expected_by_chance']:>10.2f}  {r['verdict']}")
    flagged = [r for r in results if r["n_transcript"]]
    if flagged:
        lines.append("")
        lines.append("occurrences inside annotated transcripts — the tier a PPR could act on:")
        for r in flagged:
            for h in r["transcript_hits"][:3]:
                lines.append(f"  {r['target']}  {h}")
    lines.append("")
    lines.append("Sequence occurrence is necessary for an off-target interaction, not "
                 "sufficient. No binding affinity is predicted here.")
    return "\n".join(lines)
