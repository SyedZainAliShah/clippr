"""Reconstructing what the deposited kit actually assembles, from sequences you supply.

`parts.py` answers which modules realise a target. It deliberately carries no insert DNA --
see `NOTICE.md` -- so it can say "use these ten plasmids" but not "here is the sequence you
would get". This module is the one place where the kit's sequences enter, and they enter
from a file the user supplies.

**Why this exists.** The kit-versus-synthesis homology comparison has now produced two wrong
numbers, and both failed the same way: the sequence being compared was never an artefact
anyone could inspect. The first summed every module two members shared and called it a
contiguous tract. The second compared module identity at matching slot indices, which cannot
see two members using the same modules at different positions -- the actual worst case. A
number computed from an exported sequence can be checked; a number computed from bookkeeping
cannot.

**Scope, stated rather than assumed.** What is reconstructed here is the joined *insert*
sequence: the modules, spliced at their four-base overlaps. It is not an expression
construct. The acceptor backbone is not part of the kit records and is not inferred; a
vector-level claim needs vector sequence supplied separately, and until it is, this result
describes inserts only.
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path

#: Table S1's overhang columns use a prime (U+02B9), not an apostrophe. Both are accepted,
#: because a file round-tripped through a spreadsheet often loses the distinction.
_COLUMN_ALIASES = {
    "plasmid id": "plasmid_id",
    "insert sequence": "insert",
    "5ʹ overhang": "five_overhang",
    "5' overhang": "five_overhang",
    "3ʹ overhang": "three_overhang",
    "3' overhang": "three_overhang",
}
_REQUIRED = ("plasmid_id", "insert")


class InsertsUnavailable(Exception):
    """Sequence reconstruction was asked for without usable insert records.

    Distinct from a malformed file: the module selection remains valid and usable, only the
    sequence-level result is unavailable. Callers should report the plan and say so, rather
    than presenting an absent reconstruction as a failed design.
    """


class InsertRecordsInvalid(Exception):
    """The supplied records cannot be trusted, with the offending rows named."""


@dataclass(frozen=True)
class InsertSource:
    """Where the insert sequences came from, so a product can be tied back to a file."""

    path: str
    sha256: str
    n_records: int

    def as_dict(self) -> dict:
        return {"path": self.path, "sha256": self.sha256, "n_records": self.n_records}


def _normalise_columns(columns) -> dict[str, str]:
    found = {}
    for col in columns:
        key = _COLUMN_ALIASES.get(str(col).strip().lower())
        if key and key not in found:
            found[key] = col
    return found


def load_inserts(path: str | Path, *, check_lengths: bool = True
                 ) -> tuple[dict[str, str], InsertSource]:
    """Read insert sequences from a supplied Table S1, validating before returning them.

    Returns the sequences by plasmid id plus the source record. Every rejection names the
    rows responsible: a file that is silently wrong is worse than one that refuses to load,
    because the reconstruction it feeds would look entirely plausible.

    `check_lengths` cross-checks each sequence against the insert length already recorded in
    `parts.json`. That is what catches a file that parses cleanly but is not this kit.
    """
    import pandas as pd

    src = Path(path)
    if not src.is_file():
        raise InsertsUnavailable(
            f"no insert records at {src}. Supply the paper's Supplementary Table S1; the "
            f"package does not ship it (see NOTICE.md). Module selection still works "
            f"without it -- only sequence reconstruction is unavailable.")

    raw = src.read_bytes()
    table = pd.read_excel(src) if src.suffix.lower() in (".xlsx", ".xls") \
        else pd.read_csv(src)
    columns = _normalise_columns(table.columns)
    missing = [c for c in _REQUIRED if c not in columns]
    if missing:
        raise InsertRecordsInvalid(
            f"{src} is missing required column(s) {missing}; found {list(table.columns)}")

    sequences: dict[str, str] = {}
    duplicates, bad_alphabet, empty = [], [], []
    for _, row in table.iterrows():
        pid = str(row[columns["plasmid_id"]]).strip()
        if not pid or pid.lower() == "nan":
            continue
        seq = str(row[columns["insert"]]).upper()
        cleaned = re.sub(r"\s", "", seq)
        if not cleaned or cleaned == "NAN":
            empty.append(pid)
            continue
        if re.search(r"[^ACGT]", cleaned):
            bad_alphabet.append(f"{pid}:{sorted(set(re.findall(r'[^ACGT]', cleaned)))}")
            continue
        if pid in sequences and sequences[pid] != cleaned:
            duplicates.append(pid)
            continue
        sequences[pid] = cleaned

    problems = []
    if duplicates:
        problems.append(f"conflicting duplicate ids {sorted(set(duplicates))}")
    if bad_alphabet:
        problems.append(f"non-ACGT sequence in {sorted(set(bad_alphabet))}")
    if empty:
        problems.append(f"empty insert sequence for {sorted(set(empty))}")
    if problems:
        raise InsertRecordsInvalid(f"{src}: " + "; ".join(problems))
    if not sequences:
        raise InsertsUnavailable(f"{src} parsed but contained no usable insert records")

    if check_lengths:
        _check_against_index(src, sequences)

    return sequences, InsertSource(str(src), hashlib.sha256(raw).hexdigest(), len(sequences))


def _check_against_index(src: Path, sequences: dict[str, str]) -> None:
    """Confirm the supplied file describes the kit `parts.json` already indexes.

    `parts.json` records every module's insert *length* without its sequence, which makes it
    an independent check on the file: a table that parses but disagrees on lengths is a
    different dataset, and reconstructing from it would produce a confident wrong answer.
    """
    from .parts import inventory

    known = {m.plasmid_id: m.insert_length for m in inventory()}
    shared = sorted(set(known) & set(sequences))
    if not shared:
        raise InsertRecordsInvalid(
            f"{src}: none of its {len(sequences)} plasmid ids appear in the kit index; "
            f"this does not look like the GRASP module table")
    wrong = [f"{pid} has {len(sequences[pid])} nt, index says {known[pid]}"
             for pid in shared if len(sequences[pid]) != known[pid]]
    if wrong:
        raise InsertRecordsInvalid(
            f"{src}: insert length disagrees with the kit index for {len(wrong)} module(s): "
            + "; ".join(wrong[:5]) + (" ..." if len(wrong) > 5 else ""))


class ScopeMismatch(Exception):
    """Two products were compared over a region on which they do not agree."""


class OverlapMismatch(Exception):
    """Two modules that should join do not share the four bases they declare."""


@dataclass(frozen=True)
class Junction:
    """One four-base join, and where it sits in the product it belongs to."""

    left: str
    right: str
    overhang: str
    start: int

    def as_dict(self) -> dict:
        return {"left": self.left, "right": self.right, "overhang": self.overhang,
                "start": self.start, "end": self.start + len(self.overhang)}


@dataclass(frozen=True)
class StageProduct:
    """What one assembly reaction yields, with the identity of every join inside it."""

    name: str
    sequence: str
    modules: tuple[str, ...]
    junctions: tuple[Junction, ...]

    def as_dict(self) -> dict:
        return {"name": self.name, "length": len(self.sequence),
                "modules": list(self.modules),
                "junctions": [j.as_dict() for j in self.junctions],
                "sequence": self.sequence}


@dataclass(frozen=True)
class Reconstruction:
    """Every stage the kit route passes through, and the scope of what that covers."""

    target: str
    stages: tuple[StageProduct, ...]
    final: StageProduct
    source: InsertSource
    scope: str

    def as_dict(self) -> dict:
        return {"target": self.target, "scope": self.scope,
                "source": self.source.as_dict(),
                "stages": [s.as_dict() for s in self.stages],
                "final": self.final.as_dict()}


#: What a reconstruction covers. Stated on every result rather than left to the reader,
#: because an insert-only product presented without it reads as an expression construct.
INSERT_SCOPE = ("joined module inserts only; no acceptor backbone supplied, so this is not "
                "an expression construct and carries no vector context")


def _join(modules, inserts: dict[str, str], *, where: str) -> tuple[str, list[Junction]]:
    """Splice consecutive inserts at their shared four bases, checking every join.

    The overlap is verified three ways at once -- the left insert's last four bases, the
    right insert's first four, and the overhang both modules declare in the kit index. A
    join that passes all three is not a coincidence of the sequence file.
    """
    ids = [m.plasmid_id for m in modules]
    absent = [i for i in ids if i not in inserts]
    if absent:
        raise InsertsUnavailable(
            f"{where}: no insert sequence supplied for {absent}; reconstruction is "
            f"unavailable for this target, though its module selection stands")

    sequence = inserts[ids[0]]
    junctions: list[Junction] = []
    for left, right in zip(modules, modules[1:]):
        a, b = inserts[left.plasmid_id], inserts[right.plasmid_id]
        observed = a[-4:]
        if observed != b[:4]:
            raise OverlapMismatch(
                f"{where}: {left.plasmid_id} ends {observed} but {right.plasmid_id} "
                f"starts {b[:4]}; they cannot be joined")
        if not (observed == left.three_overhang == right.five_overhang):
            raise OverlapMismatch(
                f"{where}: {left.plasmid_id}/{right.plasmid_id} overlap {observed} "
                f"disagrees with the declared overhangs "
                f"{left.three_overhang}/{right.five_overhang}")
        junctions.append(Junction(left.plasmid_id, right.plasmid_id, observed,
                                  len(sequence) - 4))
        sequence += b[4:]
    return sequence, junctions


def reconstruct(plan, inserts: dict[str, str], source: InsertSource) -> Reconstruction:
    """Build every stage product for a `parts.select` plan from supplied insert sequences.

    Stages follow the plan's own `sub_assemblies`, which exist because Golden Gate needs
    unique overhangs within one reaction -- they are not a presentation choice. The final
    product joins those sub-assembly products at their linker overhangs.

    Raises rather than returning a partial result: a reconstruction that skipped a join it
    could not verify would be exactly the unverifiable artefact this module exists to avoid.
    """
    if not plan.available:
        raise InsertsUnavailable(
            f"no kit route for {plan.target}: {plan.reason or 'no modules were selected'}")

    stages = []
    for i, sub in enumerate(plan.sub_assemblies, 1):
        seq, junctions = _join(sub, inserts, where=f"sub-assembly {i}")
        stages.append(StageProduct(f"sub-assembly {i}", seq,
                                   tuple(m.plasmid_id for m in sub), tuple(junctions)))

    seq, junctions = _join(plan.modules, inserts, where="final assembly")
    final = StageProduct("final", seq, tuple(m.plasmid_id for m in plan.modules),
                         tuple(junctions))
    return Reconstruction(plan.target, tuple(stages), final, source, INSERT_SCOPE)


#: Coding frame declared by the first module's 5' fusion site. Both deposited variants are
#: four-base MoClo CDS sites whose first base is the last base of the codon upstream of the
#: insert: `AATG` carries the initiator itself, `AGGT` continues an N-terminal fusion. The
#: frame is therefore read from the part, never searched for -- picking whichever frame had
#: no stop codons would be a heuristic dressed as a result.
_FUSION_SITE_OFFSET = {"AATG": 1, "AGGT": 1}


@dataclass(frozen=True)
class Translation:
    """The protein a reconstructed product encodes, and what the product does not contain."""

    frame: int
    protein: str
    stop_codons: int
    upstream_bases_needed: int
    note: str

    def as_dict(self) -> dict:
        return {"frame": self.frame, "length_aa": len(self.protein),
                "stop_codons": self.stop_codons,
                "upstream_bases_needed": self.upstream_bases_needed,
                "note": self.note, "protein": self.protein}


def translate(rec: Reconstruction, modules) -> Translation:
    """Translate a reconstructed product in the frame its first fusion site declares.

    `modules` is the plan's module tuple, whose first entry names the fusion site. An
    unrecognised site raises rather than guessing: a wrong frame produces a plausible-looking
    protein, which is the failure mode this refuses to risk.

    The leading bases before the frame belong to a codon the acceptor backbone completes, so
    they are reported as `upstream_bases_needed` rather than silently trimmed. The protein
    returned is therefore what the insert encodes *from its own first whole codon*, which is
    one residue short of the assembled construct's N-terminus.
    """
    from Bio.Seq import Seq

    site = modules[0].five_overhang
    if site not in _FUSION_SITE_OFFSET:
        raise OverlapMismatch(
            f"first module {modules[0].plasmid_id} declares 5' site {site}, which is not a "
            f"known CDS fusion site {sorted(_FUSION_SITE_OFFSET)}; the reading frame cannot "
            f"be derived from the part and will not be guessed")

    frame = _FUSION_SITE_OFFSET[site]
    coding = rec.final.sequence[frame:]
    remainder = len(coding) % 3
    if remainder:
        raise OverlapMismatch(
            f"{rec.target}: {len(coding)} coding bases after the {site} fusion site is not a "
            f"whole number of codons ({remainder} left over); the product is malformed")

    protein = str(Seq(coding).translate())
    return Translation(
        frame=frame,
        protein=protein,
        stop_codons=protein.count("*"),
        upstream_bases_needed=3 - frame,
        note=(f"frame declared by the {site} fusion site; the {frame} base(s) before it "
              f"complete a codon whose remaining {3 - frame} base(s) come from the acceptor "
              f"backbone, which is not supplied here"),
    )


def common_coding_region(rec: Reconstruction, translation: Translation,
                         cds: str) -> tuple[str, str, dict]:
    """The stretch the two realisation routes can honestly be compared over.

    A kit product and a synthesis CDS are not the same object. The kit insert opens with a
    fusion site whose first base belongs to a codon the backbone completes, and the synthesis
    CDS carries two N-terminal residues the kit realisation does not. Comparing them whole --
    901 nt against 906 -- measures that offset rather than any shared DNA.

    Both routes encode the same protein from the kit's second codon onward, so that is the
    declared region: the kit product from its second whole codon, and the CDS from its
    fourth. Returns both slices and the coordinates they were taken at, so a reader can cut
    the same region themselves.
    """
    from Bio.Seq import Seq

    kit_start = translation.frame + 3
    cds_start = 9
    kit_region, cds_region = rec.final.sequence[kit_start:], cds[cds_start:]

    # The basis sentence is a claim, so check it rather than print it. Handed a CDS for a
    # different target this returned a plausible tract length under the assertion "both
    # encode the same protein" -- numbers with an unjustified sentence attached, which is
    # the failure this module was written to stop.
    # Truncating each region to whole codons before translating hid a real disagreement:
    # a CDS with one extra base compared 897 nt against 898 and still reported a verified
    # shared protein, because the odd base was silently dropped.
    ragged = [n for n, r in (("kit", kit_region), ("cds", cds_region)) if len(r) % 3]
    if ragged or len(kit_region) != len(cds_region):
        raise ScopeMismatch(
            f"{rec.target}: regions are not comparable -- kit {len(kit_region)} nt, cds "
            f"{len(cds_region)} nt"
            + (f"; not a whole number of codons: {ragged}" if ragged else "")
            + ". Supply a CDS whose coding region matches, or compare explicitly as "
              "unrelated sequences.")

    kit_protein = str(Seq(kit_region).translate())
    cds_protein = str(Seq(cds_region).translate())
    if kit_protein != cds_protein:
        raise ScopeMismatch(
            f"{rec.target}: the kit product and the supplied CDS do not encode the same "
            f"protein over the declared region (kit {kit_start}:, cds {cds_start}:; "
            f"{len(kit_protein)} aa vs {len(cds_protein)} aa). They are not comparable "
            f"here -- check the CDS belongs to this target.")

    return kit_region, cds_region, {
        "kit_start": kit_start,
        "cds_start": cds_start,
        "common_protein_aa": len(kit_protein),
        "basis": ("kit product from its second whole codon; CDS from its fourth residue. "
                  "Verified to encode the same protein over this region"),
    }


def compare_to_synthesis(rec: Reconstruction, translation: Translation, cds: str) -> dict:
    """Longest exact DNA shared between the kit and synthesis routes for one target.

    Uses `homology.longest_shared`, the same implementation the library-level comparison
    uses, so the two numbers mean the same thing. Returns the matching sequence and where it
    sits in each route, because a length with no coordinates cannot be checked.
    """
    from .homology import longest_shared

    kit_region, cds_region, coords = common_coding_region(rec, translation, cds)
    length, start_a, start_b = longest_shared(kit_region, cds_region)[:3]
    return {
        "target": rec.target,
        "region": coords,
        "kit_region_nt": len(kit_region),
        "cds_region_nt": len(cds_region),
        "longest_shared_nt": length,
        "kit_offset": start_a,
        "cds_offset": start_b,
        "match": kit_region[start_a:start_a + length] if length else "",
    }


def write(rec: Reconstruction, path: str | Path, *,
          translation: Translation | None = None) -> Path:
    """Write the reconstruction as JSON: every stage, junction, coordinate and its source.

    JSON rather than a sequence file because the point is that an independent checker can
    re-derive the product from the raw records and compare, which needs the joins and their
    coordinates, not just the final string.
    """
    import json

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = rec.as_dict()
    if translation is not None:
        payload["translation"] = translation.as_dict()
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out


def _module_spans(stage: StageProduct) -> list[tuple[str, int, int]]:
    """Where each module sits in its stage product, as half-open [start, end) coordinates.

    Consecutive modules overlap by the four bases of their shared overhang, so the spans
    genuinely overlap. That is the assembly, not an off-by-one: the junction bases belong to
    both neighbours, which is why they anneal.
    """
    if not stage.junctions:
        return [(stage.modules[0], 0, len(stage.sequence))]
    spans = []
    starts = [0] + [j.start for j in stage.junctions]
    ends = [j.start + len(j.overhang) for j in stage.junctions] + [len(stage.sequence)]
    for module, start, end in zip(stage.modules, starts, ends):
        spans.append((module, start, end))
    return spans


def write_fasta(rec: Reconstruction, path: str | Path, *, stages: bool = True) -> Path:
    """Write the reconstructed products as FASTA.

    Every description repeats the scope, because a FASTA file travels without its JSON and a
    bare insert-only sequence read as an expression construct is the misuse this package is
    built to prevent.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    records = list(rec.stages) + [rec.final] if stages else [rec.final]
    lines = []
    for stage in records:
        name = stage.name.replace(" ", "_")
        lines.append(f">{rec.target}|{name} length={len(stage.sequence)} "
                     f"modules={'+'.join(stage.modules)} source={rec.source.sha256[:16]} "
                     f"scope={rec.scope}")
        lines.extend(stage.sequence[i:i + 60] for i in range(0, len(stage.sequence), 60))
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def write_genbank(rec: Reconstruction, path: str | Path, *,
                  translation: Translation | None = None) -> Path:
    """Write the final product as annotated GenBank: one feature per module and per junction.

    The CDS feature is emitted only when a `Translation` is supplied, and it carries that
    translation's own frame and protein rather than letting the writer re-derive them. A
    GenBank file that computed its own frame could disagree with the JSON beside it, and the
    reader would have no way to tell which was the design.
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio import SeqIO

    stage = rec.final
    record = SeqRecord(Seq(stage.sequence), id=rec.target[:16], name=rec.target[:16],
                       description=f"CLIPPR kit reconstruction of {rec.target}")
    record.annotations = {
        "molecule_type": "DNA", "topology": "linear", "data_file_division": "SYN",
        "comment": (f"{rec.scope}. Reconstructed from supplied insert records "
                    f"{rec.source.sha256[:16]} ({rec.source.path}). Module spans overlap by "
                    f"the four bases of each shared overhang."),
    }
    for module, start, end in _module_spans(stage):
        record.features.append(SeqFeature(
            SimpleLocation(start, end), type="misc_feature",
            qualifiers={"label": [module], "note": ["kit module insert"]}))
    for junction in stage.junctions:
        record.features.append(SeqFeature(
            SimpleLocation(junction.start, junction.start + len(junction.overhang)),
            type="misc_binding",
            qualifiers={"label": [f"{junction.left}/{junction.right}"],
                        "note": [f"Golden Gate overhang {junction.overhang}"]}))
    if translation is not None:
        end = translation.frame + 3 * len(translation.protein)
        record.features.append(SeqFeature(
            SimpleLocation(translation.frame, end), type="CDS",
            qualifiers={"label": [f"{rec.target} PPR"], "codon_start": [1],
                        "translation": [translation.protein],
                        "note": [translation.note]}))

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as fh:
        SeqIO.write(record, fh, "genbank")
    return out
