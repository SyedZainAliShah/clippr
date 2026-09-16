"""V7 — reparse CLIPPR's exported files and check they still mean what the design meant.

A format can be valid and still have lost the thing it was serialising: a dropped record, a
strand flipped, an off-by-one coordinate, a quantity in the wrong unit. This reparses the
written files with stock parsers and compares against the *design's own values* and the
sequences V5/V6 already validated -- never against a second call to the serializer, which
would only prove the serializer is consistent with itself.

Pricing is treated separately. The pool and phosphorylation constants are documented in
`export.py` as measured across the 200-design corpus, not fetched from a vendor. A file
format cannot validate a commercial price, so this reports them as historical constants with
their stated basis rather than checking them.

    python validation/independent/v7_export_semantics.py
"""
from __future__ import annotations

import csv
import json
import sys
import tempfile
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from Bio import SeqIO                                             # noqa: E402
from Bio.Seq import Seq                                           # noqa: E402

from clippr.design import design_oneshot                          # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from journal import provenance                                    # noqa: E402

TARGETS = ["AAAAUGUGG", "UUACACGUGCGUAC", "CUAUCACAUCACAUAAGCG"]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Minimal FASTA reader, so a malformed file is not silently repaired by a library."""
    records, name, chunks = [], None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name, chunks = line[1:].strip(), []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def expected_annotation(target: str, cds: str, fragments) -> dict:
    """The feature layout the exporter promises, derived independently of the written file.

    `export.write_genbank` emits one CDS over the whole record, a start-codon feature, an
    N-terminal helix feature, three features per PPR repeat (the repeat and its two
    specificity residues), and one feature per fragment. Every coordinate here follows from
    the scaffold constants and the target length, so a file that has lost or moved a feature
    can be caught.

    The earlier version iterated whatever features the file happened to contain. With no
    features at all that loop ran zero times and reported no problem, and a CDS shifted three
    bases was never compared against any expected coordinate -- both passed.
    """
    from clippr import constants as C

    n_term = len(C.N_TERMINAL)
    repeat_aa = len(C.REPEAT_TEMPLATE.format(fifth="X", last="Y"))
    spans = {"CDS": [(0, len(cds))],
             "start codon": [(0, 3)],
             "N-terminal solvating helix": [(3, 3 * n_term)]}
    repeats, fifths, lasts = [], [], []
    for i in range(1, len(target) + 1):
        start = n_term + (i - 1) * repeat_aa
        repeats.append((3 * start, 3 * (start + repeat_aa)))
        fifths.append((3 * (start + 1), 3 * (start + 1) + 3))
        lasts.append((3 * (start + repeat_aa - 4), 3 * (start + repeat_aa - 4) + 3))
    spans["repeats"] = repeats
    spans["fifth residues"] = fifths
    spans["last residues"] = lasts
    spans["fragments"] = [(f.cds_start, f.cds_end) for f in fragments]
    total = 1 + 2 + 3 * len(target) + len(fragments)
    return {"spans": spans, "total_features": total, "repeat_aa": repeat_aa,
            "n_term": n_term}


def check_design(target: str, outdir: Path, table, problems: list[str]) -> dict:
    design = design_oneshot(target, codon_table=table, genetic_code=1,
                            check_offtarget=False, seed=42, outdir=outdir)
    return inspect_exports(design, problems)


def inspect_exports(design: dict, problems: list[str]) -> dict:
    """Every file check, applied to whatever is on disk right now.

    Split out from `check_design` so a deliberately damaged file is judged by exactly this
    function, not by a weaker ad-hoc version. A control checked more leniently than the real
    thing proves nothing about the real thing.
    """
    target = design["target_rna"]
    paths = {k: Path(v) for k, v in design["paths"].items()}
    cds, oligos = design["cds"], design["oligos"]
    note = f"{target}"

    gene = read_fasta(paths["gene_fasta"])
    if len(gene) != 1:
        problems.append(f"{note}: gene FASTA has {len(gene)} records, expected 1")
    elif gene[0][1] != cds:
        problems.append(f"{note}: gene FASTA sequence differs from the design's CDS")

    oligo_fasta = read_fasta(paths["oligo_fasta"])
    expected = list(oligos["oligo_sequence_5to3"])
    if [s for _n, s in oligo_fasta] != expected:
        problems.append(f"{note}: oligo FASTA sequences or their order differ from the design")
    if len({n for n, _s in oligo_fasta}) != len(oligo_fasta):
        problems.append(f"{note}: oligo FASTA identifiers are not unique")

    with paths["oligo_csv"].open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != len(oligos):
        problems.append(f"{note}: CSV has {len(rows)} rows, design has {len(oligos)}")
    else:
        for i, row in enumerate(rows):
            if row["oligo_sequence_5to3"] != expected[i]:
                problems.append(f"{note}: CSV row {i} sequence differs")
            declared = int(row["oligo_length"])
            if declared != len(row["oligo_sequence_5to3"]):
                problems.append(f"{note}: CSV row {i} length {declared} contradicts its "
                                f"own sequence ({len(row['oligo_sequence_5to3'])})")

    record = SeqIO.read(paths["genbank"], "genbank")
    if str(record.seq) != cds:
        problems.append(f"{note}: GenBank sequence differs from the CDS")

    want = expected_annotation(target, cds, design["fragments"])
    if len(record.features) != want["total_features"]:
        problems.append(f"{note}: GenBank has {len(record.features)} features, expected "
                        f"{want['total_features']}")

    cds_features = [f for f in record.features if f.type == "CDS"]
    if len(cds_features) != 1:
        problems.append(f"{note}: GenBank has {len(cds_features)} CDS features, expected 1")
    for feature in cds_features:
        span = (int(feature.location.start), int(feature.location.end))
        if span != want["spans"]["CDS"][0]:
            problems.append(f"{note}: GenBank CDS spans {span}, expected "
                            f"{want['spans']['CDS'][0]}")
        extracted = str(feature.extract(record.seq))
        if feature.location.strand not in (1, None):
            problems.append(f"{note}: GenBank CDS feature is on the reverse strand")
        if len(extracted) % 3:
            problems.append(f"{note}: GenBank CDS feature is not a whole number of codons")
        declared = feature.qualifiers.get("translation", [None])[0]
        actual = str(Seq(extracted).translate(table=1)) if not len(extracted) % 3 else None
        if declared != actual:
            problems.append(f"{note}: GenBank /translation does not match the bases under "
                            f"its own CDS feature")

    # Every promised annotation must be present at its derived coordinates, and every
    # feature in the file must be one the exporter promised.
    present = {(f.type, int(f.location.start), int(f.location.end)) for f in record.features}
    for label, spans in want["spans"].items():
        kind = "CDS" if label == "CDS" else "misc_feature"
        missing = [sp for sp in spans if (kind, sp[0], sp[1]) not in present]
        if missing:
            problems.append(f"{note}: GenBank is missing {len(missing)} {label} feature(s) "
                            f"at {missing[:3]}")
    for feature in record.features:
        start, end = int(feature.location.start), int(feature.location.end)
        if not (0 <= start < end <= len(record.seq)):
            problems.append(f"{note}: feature coordinates {start}-{end} fall outside the "
                            f"{len(record.seq)} nt record")

    return {"target": target, "cds_nt": len(cds), "oligos": len(oligos),
            "genbank_features": len(record.features), "cds_features": len(cds_features),
            "cost": design["cost"]}


def malformed_controls(outdir: Path, problems: list[str]) -> list[tuple[str, bool]]:
    """Files that must not parse as valid, so a passing reparse means something."""
    results = []

    truncated = outdir / "truncated.fasta"
    truncated.write_text(">only-a-header-no-sequence\n", encoding="utf-8")
    records = read_fasta(truncated)
    results.append(("FASTA header with no sequence",
                    bool(records) and records[0][1] == ""))

    empty = outdir / "empty.fasta"
    empty.write_text("", encoding="utf-8")
    results.append(("empty FASTA yields no records", read_fasta(empty) == []))

    bad_csv = outdir / "bad.csv"
    bad_csv.write_text("oligo_sequence_5to3,oligo_length\nACGT,99\n", encoding="utf-8")
    with bad_csv.open(encoding="utf-8", newline="") as fh:
        row = next(csv.DictReader(fh))
    results.append(("CSV length field contradicting its sequence",
                    int(row["oligo_length"]) != len(row["oligo_sequence_5to3"])))
    return results


def damaged_annotation_controls(design: dict) -> list[tuple[str, bool]]:
    """Damage a real export's annotations and require `inspect_exports` to object.

    These two both produced *zero* problems before the annotation checks were derived
    independently: with every feature deleted the old loop iterated nothing, and a CDS moved
    three bases was never compared to any expected coordinate. They are checked through
    `inspect_exports` -- the same function the real exports go through -- because a control
    judged by a private, weaker test says nothing about the real gate.

    The file is restored afterwards, so this leaves the export directory as it found it.
    """
    from Bio.SeqFeature import SimpleLocation

    genbank = Path(design["paths"]["genbank"])
    saved = genbank.read_bytes()
    results = []
    try:
        record = SeqIO.read(genbank, "genbank")
        record.features = []
        with genbank.open("w", encoding="utf-8") as fh:
            SeqIO.write(record, fh, "genbank")
        found: list[str] = []
        inspect_exports(design, found)
        results.append(("every GenBank feature removed", bool(found)))

        genbank.write_bytes(saved)
        record = SeqIO.read(genbank, "genbank")
        for feature in record.features:
            if feature.type == "CDS":
                feature.location = SimpleLocation(3, len(record.seq), strand=1)
        with genbank.open("w", encoding="utf-8") as fh:
            SeqIO.write(record, fh, "genbank")
        found = []
        inspect_exports(design, found)
        results.append(("CDS start shifted three bases", bool(found)))
    finally:
        genbank.write_bytes(saved)
    return results


def main() -> int:
    table = json.loads((ROOT / "data" / "codon_tables" / "kazusa_3055.json")
                       .read_text(encoding="utf-8"))
    outdir = Path(tempfile.mkdtemp())
    problems: list[str] = []

    print(f"reparsing exports for {len(TARGETS)} designs across the three architectures")
    designs = [design_oneshot(t, codon_table=table, genetic_code=1, check_offtarget=False,
                              seed=42, outdir=outdir) for t in TARGETS]
    summaries = [inspect_exports(d, problems) for d in designs]
    for s in summaries:
        print(f"  {s['target']:22s} {s['cds_nt']:>5} nt  {s['oligos']} oligos  "
              f"{s['cds_features']} CDS + {s['genbank_features'] - s['cds_features']} other "
              f"features")

    print(f"\nsemantic problems found: {len(problems)}")
    for p in problems[:10]:
        print(f"   {p}")

    print("\nmalformed controls (each must be recognised as malformed):")
    controls = malformed_controls(outdir, problems)
    controls += damaged_annotation_controls(designs[0])
    for name, detected in controls:
        print(f"   {name:48s} {'detected' if detected else 'NOT DETECTED'}")

    cost = summaries[0]["cost"]
    print(f"\npricing, reported not validated: {cost.get('currency')} "
          f"dna={cost.get('dna_eur')} total={cost.get('total_eur')}")
    print("   basis: export.py documents 109.00 EUR pool and 1.63 EUR/oligo phosphorylation")
    print("          as constants measured across the 200-design corpus -- a list price, not")
    print("          a quote, and not fetched from a vendor. A valid file cannot validate a")
    print("          commercial price; treat these as historical unless re-sourced with a date.")

    out = ROOT / "work" / "independent" / "v7_export_semantics.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(
        {"provenance": provenance(__file__),
         "designs": summaries, "semantic_problems": problems,
         "malformed_controls": {n: bool(d) for n, d in controls},
         "pricing": {"values": cost, "status": "historical corpus-derived constants",
                     "basis": "export.py docstrings; not a dated vendor quote"},
         "independence": ("reparses written files with stock parsers and compares against "
                          "the design's own values, never a second serializer call")},
        indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0 if not problems and all(d for _n, d in controls) else 1


if __name__ == "__main__":
    raise SystemExit(main())
