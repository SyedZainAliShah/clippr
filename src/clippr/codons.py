"""Codon optimisation for a synthetic PPR coding sequence.

Two things make this harder than a normal codon-optimisation job.

**The organelle is not a detail.** Chlamydomonas *nuclear* Leu usage is CUG 0.73 and
GC-rich; its *chloroplast* usage is UUA 0.33 and AT-rich. Optimising a chloroplast
construct against a nuclear table produces materially wrong DNA, and the genetic code
differs too (NCBI table 1 versus 11). Both the table and the code are therefore
parameters with no silent default: `optimize_cds` will not guess.

**The protein is tandem repeats of one 31-residue template**, so the naive
"use the single best codon for every residue" objective emits the same nucleotides once
per repeat. That is a synthesis-grade repeat problem, not a cosmetic one -- every one of
the 200 reference designs came back flagged. `optimize_cds` therefore defaults to making
every 20-mer in the sequence unique, measured over four targets spanning 9S and 19S:

    setting             dup 20-mers   longest exact repeat    CAI     GC
    best codon only             425                  95 nt   0.843  0.610
    unique 15-mers                0                  14 nt   0.761  0.595
    unique 20-mers                0                  19 nt   0.793  0.600
    reference (E. coli)         129                  44 nt   0.411  0.447

Uniquifying costs about six percent of codon adaptation and removes nearly all repeats.
The "longest repeat" figure is a ceiling rather than an independent result -- unique
20-mers caps it at 19 by construction -- so the load-bearing comparison is 95 nt against
44 nt for an E. coli-optimised baseline. k=20 is the default because it dominates k=15: same duplicate
count, less adaptation given up.

Across the full 200-design corpus (`validation/compare_codons.py`), **197 of 200 come out
with no repeated 20-mer at all; 3 retain at least one.** Uniquifying is an *objective*,
not a constraint, so it is best-effort and cannot be promised -- treat the three as
designs to inspect, not as failures. All 200 satisfy every hard constraint.

The baseline row scores that DNA against the *Chlamydomonas* table. Its low CAI is
not a defect in its optimiser; it targeted E. coli, which is a different host from this
lab's. That mismatch is the reason this module exists.

**CAI is not an expression claim.** It measures how closely codon usage matches a
reference table, and a higher value does not establish that a sequence is better
expressed. The defensible statement is "the baseline sequence is optimised for a
different host", not "our sequence expresses better" -- that would need measurement.
"""
from __future__ import annotations

import json
import random
from collections.abc import Mapping, Sequence
from pathlib import Path

from . import constants as C

CACHE = Path(__file__).resolve().parents[2] / "data" / "codon_tables"

#: Backwards-compatible alias for the default profile. Prefer naming a profile --
#: `constants.ENZYME_PROFILES` records *why* each site is excluded, which is not the same
#: claim for all four: BsaI/BbsI are assembly chemistry, SapI is an iGEM RFC[1000]
#: standard requirement, BsmBI is a preference about keeping later MoClo levels open.
BLACKLIST_ENZYMES: tuple[str, ...] = C.enzymes_for()

#: Chlamydomonas reinhardtii, NCBI taxonomy id. The wet lab's organism.
CHLAMYDOMONAS_TAXID: int = 3055


def _table_name(genetic_code: int) -> str:
    """NCBI table id -> the name Biopython and DNA Chisel index tables by.

    Both take a name, not a number: `reverse_translate(table=1)` raises `KeyError: '1'`.
    Resolving through Biopython keeps the public API numeric, which is how NCBI, the
    literature and the rest of this package refer to genetic codes.
    """
    from Bio.Data import CodonTable

    try:
        names = CodonTable.unambiguous_dna_by_id[genetic_code].names
    except KeyError as e:
        raise ValueError(f"unknown NCBI genetic code table {genetic_code}") from e
    return next(n for n in names if n)


def _as_dna(table: Mapping[str, Mapping[str, float]]) -> dict[str, dict[str, float]]:
    """Kazusa tables are spelled in RNA; DNA Chisel wants DNA codons."""
    return {aa: {c.replace("U", "T"): float(f) for c, f in codons.items()}
            for aa, codons in table.items()}


def complete_table(table: Mapping[str, Mapping[str, float]], genetic_code: int = 1
                   ) -> dict[str, dict[str, float]]:
    """Fill in codons an amino acid never used, at frequency zero.

    A table derived from a small gene set, or supplied by hand, can legitimately omit a
    rare codon -- a CDS FASTA that never uses AGT yields a serine entry without it. DNA
    Chisel indexes its table by whatever codon it is currently looking at, so such a gap
    surfaces far from its cause as a bare `KeyError: 'AGT'` inside the solver. Filling
    the gaps with zero keeps the meaning (never use this codon) and removes the crash.
    """
    from Bio.Data import CodonTable

    forward = CodonTable.unambiguous_dna_by_id[genetic_code].forward_table
    by_aa: dict[str, list[str]] = {}
    for codon, aa in forward.items():
        by_aa.setdefault(aa, []).append(codon)

    out = {}
    for aa, codons in _as_dna(table).items():
        full = {c: 0.0 for c in by_aa.get(aa, ())}
        full.update(codons)
        out[aa] = full
    return out


def table_from_kazusa(taxid: int = CHLAMYDOMONAS_TAXID) -> dict[str, dict[str, float]]:
    """Download a codon usage table from Kazusa, caching it under `data/codon_tables/`.

    Defaults to Chlamydomonas reinhardtii nuclear. The cache means a design run is
    reproducible offline and does not depend on Kazusa being up.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    cached = CACHE / f"kazusa_{taxid}.json"
    if cached.exists():
        return json.loads(cached.read_text(encoding="utf-8"))

    import python_codon_tables as pct

    table = _as_dna(pct.download_codons_table(taxid=taxid))
    if not table:
        raise ValueError(f"Kazusa returned no codon table for taxid {taxid}")
    cached.write_text(json.dumps(table, indent=1, sort_keys=True), encoding="utf-8")
    return table


def table_from_cds_fasta(path: str | Path, genetic_code: int = 1
                         ) -> dict[str, dict[str, float]]:
    """Build a codon usage table by counting codons across a CDS FASTA.

    Organism-agnostic and offline, so it covers the cases Kazusa does not: a specific
    chloroplast genome, a strain, or an in-house gene set. Frequencies are normalised
    within each amino acid, which is what a codon-optimiser consumes.
    """
    from Bio import SeqIO
    from Bio.Data import CodonTable

    forward = CodonTable.unambiguous_dna_by_id[genetic_code].forward_table
    counts: dict[str, dict[str, int]] = {}
    n_records = n_codons = 0
    for record in SeqIO.parse(str(path), "fasta"):
        seq = str(record.seq).upper().replace("U", "T")
        n_records += 1
        for i in range(0, len(seq) - len(seq) % 3, 3):
            codon = seq[i:i + 3]
            aa = forward.get(codon)
            if aa is None:          # stop codon or an ambiguity code
                continue
            counts.setdefault(aa, {}).setdefault(codon, 0)
            counts[aa][codon] += 1
            n_codons += 1
    if not n_codons:
        raise ValueError(f"{path}: no usable codons found in {n_records} records")
    return {aa: {c: n / sum(cod.values()) for c, n in cod.items()}
            for aa, cod in counts.items()}


def table_from_csv(path: str | Path) -> dict[str, dict[str, float]]:
    """Read a user-supplied table with columns `codon,frequency`.

    Frequencies are re-normalised within each amino acid, so a table given as raw counts,
    per-thousand values or fractions all load the same way.
    """
    import csv

    from Bio.Data import CodonTable

    forward = CodonTable.unambiguous_dna_by_id[1].forward_table
    counts: dict[str, dict[str, float]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        missing = {"codon", "frequency"} - {(f or "").lower() for f in reader.fieldnames or []}
        if missing:
            raise ValueError(f"{path}: missing column(s) {sorted(missing)}")
        for row in reader:
            row = {(k or "").lower(): v for k, v in row.items()}
            codon = row["codon"].strip().upper().replace("U", "T")
            aa = forward.get(codon)
            if aa is None:
                continue
            counts.setdefault(aa, {})[codon] = float(row["frequency"])
    if not counts:
        raise ValueError(f"{path}: no usable codon rows")
    return {aa: {c: v / total for c, v in cod.items()}
            for aa, cod in counts.items() if (total := sum(cod.values()))}


def optimize_cds(
    protein: str,
    locked_sites: Mapping[int, str] | None = None,
    codon_table: Mapping[str, Mapping[str, float]] | None = None,
    genetic_code: int = 1,
    seed: int = 42,
    enzymes: Sequence[str] | str = C.DEFAULT_ENZYME_PROFILE,
    gc_bounds: tuple[float, float] = (0.35, 0.65),
    gc_window: int = 50,
    max_homopolymer: int = 4,
    unique_kmer_size: int | None = 20,
) -> dict:
    """Design a coding sequence for `protein` under the assembly's constraints.

    `locked_sites` maps a nucleotide start offset to a sequence that must appear there
    verbatim -- the chosen Golden Gate overhangs, which are fixed once
    `overhangs.best_set` has picked them. `codon_table` is required: see the module
    docstring for why there is no default. `unique_kmer_size` defaults to 20 for the
    repeat reasons in the module docstring; pass `None` to optimise codons alone.

    Returns `{'cds', 'constraints_ok', 'summary', 'objectives_score'}`. A design that
    cannot satisfy every constraint is returned with `constraints_ok=False` and the
    solver's own report in `summary`, rather than raising -- a failure is a result to
    record across a sweep, not a crash.
    """
    import dnachisel as dc

    if codon_table is None:
        raise ValueError(
            "codon_table is required: a nuclear table applied to a chloroplast construct "
            "produces materially wrong DNA. Use table_from_kazusa(), "
            "table_from_cds_fasta() or table_from_csv()."
        )
    if not protein:
        raise ValueError("empty protein")

    codon_table = complete_table(codon_table, genetic_code)
    uncovered = sorted(set(protein) - set(codon_table))
    if uncovered:
        raise ValueError(
            f"codon table has no entry for {', '.join(uncovered)}, which the protein "
            f"needs. A table counted from a CDS set only covers the amino acids that "
            f"set actually used."
        )

    random.seed(seed)
    try:
        import numpy as np
        np.random.seed(seed)
    except ImportError:      # pragma: no cover - numpy is a hard dependency
        pass

    table_name = _table_name(genetic_code)
    sequence = dc.reverse_translate(protein, table=table_name)
    constraints = [
        dc.EnforceTranslation(genetic_table=table_name, translation=protein),
        dc.EnforceGCContent(mini=gc_bounds[0], maxi=gc_bounds[1], window=gc_window),
    ]
    constraints += [dc.AvoidPattern(dc.EnzymeSitePattern(e))
                    for e in C.enzymes_for(enzymes)]
    constraints += [dc.AvoidPattern(dc.HomopolymerPattern(b, max_homopolymer + 1))
                    for b in "ACGT"]

    for start, locked in (locked_sites or {}).items():
        locked = str(locked).upper().replace("U", "T")
        if start < 0 or start + len(locked) > len(sequence):
            raise ValueError(
                f"locked site at {start} (+{len(locked)}) falls outside a "
                f"{len(sequence)}-nt coding sequence"
            )
        constraints.append(dc.EnforceSequence(sequence=locked,
                                              location=(start, start + len(locked))))

    objectives = [dc.CodonOptimize(codon_usage_table=codon_table,
                                   method="use_best_codon")]
    if unique_kmer_size:
        objectives.append(dc.UniquifyAllKmers(k=unique_kmer_size, boost=2.0))

    problem = dc.DnaOptimizationProblem(sequence=sequence, constraints=constraints,
                                        objectives=objectives, logger=None)
    ok = True
    try:
        problem.resolve_constraints()
        problem.optimize()
    except Exception as exc:                      # solver gave up on this protein
        return {"cds": problem.sequence, "constraints_ok": False,
                "summary": f"{type(exc).__name__}: {exc}", "objectives_score": None}

    cds = problem.sequence
    try:
        ok = problem.all_constraints_pass()
        summary = problem.constraints_text_summary()
    except Exception as exc:                      # pragma: no cover
        ok, summary = False, f"{type(exc).__name__}: {exc}"

    return {
        "cds": cds,
        "constraints_ok": bool(ok),
        "summary": summary,
        "objectives_score": problem.objectives_evaluations().scores_sum(),
    }
