"""Write a finished design out: order CSV, FASTA, annotated GenBank, cost estimate.

Nothing here decides anything -- by the time a design reaches this module every sequence
is fixed. These functions only serialise, so the tests that matter are round-trips: a
GenBank file we write must re-parse to the same sequence and the same feature spans.

The GenBank record is the artefact a wet lab actually opens, so it is annotated rather
than bare: the CDS, the start codon, the N-terminal helix, every PPR repeat labelled with
the RNA base it reads, both specificity residues of each repeat, and every fragment with
its junction overhang. All of that is derived from `ppr`, `arelf` and `assembly`; none of
it is transcribed from another tool's output.
"""
from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from . import constants as C

#: IDT oPools list price for the 50 pmol scale, in EUR.
#:
#: Measured constant at 109.00 across all 200 corpus designs -- 4-oligo and 7-oligo pools
#: priced identically, and total_bases from 1038 to 2064 made no difference. So this
#: number does not vary with anything a designer controls, and comparing designs on cost
#: is meaningless at this scale. It is a list price, not a quote.
OPOOL_LIST_PRICE_EUR: float = 109.00

#: Optional 5' phosphorylation, charged per oligo. Derived from the corpus:
#: 6.52 EUR / 4 oligos and 11.41 / 7 both give 1.63.
PHOSPHORYLATION_EUR_PER_OLIGO: float = 1.63


def write_oligo_csv(oligos, path: str | Path) -> Path:
    """Write the order table. Column order is preserved as given."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    oligos.to_csv(path, index=False)
    return path


def write_fasta(oligos, path: str | Path) -> Path:
    """Write the oligos as FASTA, one record per orderable fragment.

    The description carries what someone reading the file needs in order to check an
    order without opening the CSV beside it: wrap enzyme, both overhangs, and length.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for _, r in oligos.iterrows():
        seq = r["oligo_sequence_5to3"]
        lines.append(f">{r['order_fragment_id']}|{r['wrap_enzyme']}|"
                     f"{r['oh5_coding_site_5to3']}..{r['oh3_coding_site_5to3']}|{len(seq)}bp")
        lines += [seq[i:i + 60] for i in range(0, len(seq), 60)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_gene_fasta(cds: str, target_rna: str, path: str | Path,
                     architecture: str | None = None) -> Path:
    """Write the assembled coding sequence as a single FASTA record."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arch = architecture or f"{len(cds) // 3}aa"
    header = f">CLIPPR_{target_rna}|binder_CDS|{arch}|{len(cds) // 3}aa|{len(cds)}nt"
    body = [cds[i:i + 60] for i in range(0, len(cds), 60)]
    path.write_text("\n".join([header] + body) + "\n", encoding="utf-8")
    return path


def _repeat_features(cds: str, target_rna: str):
    """One feature per PPR repeat, plus its two specificity residues."""
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    from .ppr import describe

    d = describe(target_rna)
    feats = []
    n_term = len(C.N_TERMINAL)
    repeat_aa = len(C.REPEAT_TEMPLATE.format(fifth="X", last="Y"))

    feats.append(SeqFeature(
        SimpleLocation(0, 3), type="misc_feature",
        qualifiers={"label": ["start codon"],
                    "note": ["Initiating methionine; not a destination overhang."]}))
    feats.append(SeqFeature(
        SimpleLocation(3, 3 * n_term), type="misc_feature",
        qualifiers={"label": ["N-terminal solvating helix"],
                    "note": [f"{n_term} aa; precedes repeat 1."]}))

    for i, (code, base) in enumerate(zip(d["code_pairs"], d["target_rna"]), start=1):
        start = n_term + (i - 1) * repeat_aa
        feats.append(SeqFeature(
            SimpleLocation(3 * start, 3 * (start + repeat_aa)), type="misc_feature",
            qualifiers={"label": [f"PPR{i} {base} ({code})"],
                        "note": [f"{repeat_aa}-aa PPR repeat {i} reading {base}; "
                                 f"5th/last residues {code[0]}/{code[1]}."]}))
        fifth = start + 1
        last = start + repeat_aa - 4
        feats.append(SeqFeature(
            SimpleLocation(3 * fifth, 3 * fifth + 3), type="misc_feature",
            qualifiers={"label": [f"PPR{i} 5th {code[0]}"],
                        "note": [f"Specificity residue reading {base}."]}))
        feats.append(SeqFeature(
            SimpleLocation(3 * last, 3 * last + 3), type="misc_feature",
            qualifiers={"label": [f"PPR{i} last {code[1]}"],
                        "note": [f"Specificity residue reading {base}."]}))
    return feats


def write_genbank(
    cds: str,
    fragments: Sequence,
    target_rna: str,
    path: str | Path,
    genetic_code: int = 1,
    enzyme: str = "BsaI",
    destination: tuple[str, str] | None = None,
) -> Path:
    """Write the coding sequence as an annotated GenBank record.

    `fragments` are `assembly.Fragment` objects; each becomes a feature carrying its
    junction overhang, so the assembly plan is legible in a sequence viewer rather than
    only in a spreadsheet.
    """
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    dest5, dest3 = destination or C.DESTINATION_OVERHANGS["level0"]

    record = SeqRecord(
        Seq(cds),
        id=f"CLIPPR_{target_rna}"[:16],
        name=f"CLIPPR_{target_rna}"[:16],
        description=f"CLIPPR PPR binder CDS targeting {target_rna}.",
        annotations={"molecule_type": "DNA", "topology": "linear",
                     "data_file_division": "SYN",
                     "comment": (
                         f"CLIPPR PPR binder ORF for {target_rna}. Wrap enzyme {enzyme}. "
                         f"Destination overhangs {dest5}/{dest3} sit outside this record; "
                         f"they are contributed by the backbone, not the CDS.")},
    )

    record.features.append(SeqFeature(
        SimpleLocation(0, len(cds)), type="CDS",
        qualifiers={"label": [f"CLIPPR binder CDS ({target_rna})"],
                    "codon_start": ["1"], "transl_table": [str(genetic_code)],
                    "translation": [str(Seq(cds).translate(table=genetic_code))]}))
    record.features += _repeat_features(cds, target_rna)

    for f in fragments:
        record.features.append(SeqFeature(
            SimpleLocation(f.cds_start, f.cds_end), type="misc_feature",
            qualifiers={"label": [f"fragment {f.fragment_id}"],
                        "note": [f"Order {f.assembly_order}; residues {f.aa_start}-"
                                 f"{f.aa_end}; 3' junction overhang "
                                 f"{f.oh3_coding_site}."]}))

    from Bio import SeqIO
    with open(path, "w", encoding="utf-8") as fh:
        SeqIO.write(record, fh, "genbank")
    return path


def opool_quote(oligos, phosphorylate_5prime: bool = False) -> dict:
    """Estimate the cost of ordering these fragments as an IDT oPool.

    **A list price, not a quote.** The pool price is flat: measured at 109.00 EUR for
    every one of the 200 corpus designs, whether the pool held 4 oligos or 7 and whether
    it totalled 1,038 or 2,064 bases. Cost therefore does not discriminate between
    designs at this scale, and should not be used to choose one.
    """
    lengths = [len(s) for s in oligos["oligo_sequence_5to3"]]
    if not lengths:
        raise ValueError("no oligos to quote")
    phospho = round(PHOSPHORYLATION_EUR_PER_OLIGO * len(lengths), 2)
    return {
        "vendor": "IDT",
        "product": "oPools DNA",
        "currency": "EUR",
        "n_oligos": len(lengths),
        "total_bases": sum(lengths),
        "min_oligo_nt": min(lengths),
        "max_oligo_nt": max(lengths),
        "dna_eur": OPOOL_LIST_PRICE_EUR,
        "phospho_eur": phospho,
        "total_eur": round(OPOOL_LIST_PRICE_EUR + (phospho if phosphorylate_5prime else 0), 2),
        "phosphorylate_5prime": phosphorylate_5prime,
        "list_price_not_a_quote": True,
    }
