"""Deposited-kit reconstruction: what it builds, and what it refuses to build.

The malformed-input tests build their own small tables and always run. The tests that need
real kit sequences skip when Supplementary Table S1 has not been supplied, because the
package deliberately does not ship it -- see `NOTICE.md`.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from clippr.parts import inventory, select
from clippr.products import (InsertRecordsInvalid, InsertsUnavailable, OverlapMismatch,
                             compare_to_synthesis, load_inserts, reconstruct, translate,
                             write)

TABLE_S1 = Path(__file__).resolve().parents[1] / "data" / "grasp_supp" / "Table S1.xlsx"
SIX = ["AAAAUGUGG", "GCUAAAGAC", "UUACACGUG", "CGUACGUAC", "AUCGAUCGA", "GGCCAAUUG"]
#: Lengths of the joined insert product, by architecture. Not CDS lengths -- see
#: `test_scope_is_stated_and_not_a_cds`.
EXPECTED_NT = {9: 901, 14: 1366, 19: 1831}


@pytest.fixture(scope="module")
def kit():
    if not TABLE_S1.is_file():
        pytest.skip(f"supply {TABLE_S1.name} to run kit reconstruction tests")
    return load_inserts(TABLE_S1)


def _csv(tmp_path, rows, name="table.csv"):
    path = tmp_path / name
    lines = ["Plasmid ID,Insert Sequence"] + [f"{a},{b}" for a, b in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class TestRejections:
    def test_absent_file_is_unavailable_not_a_crash(self, tmp_path):
        with pytest.raises(InsertsUnavailable, match="Module selection still works"):
            load_inserts(tmp_path / "nothing.xlsx")

    def test_missing_column_is_named(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("Plasmid ID,Something Else\na,b\n", encoding="utf-8")
        with pytest.raises(InsertRecordsInvalid, match="missing required column"):
            load_inserts(path)

    def test_non_acgt_sequence_is_named(self, tmp_path):
        real = inventory()[0]
        path = _csv(tmp_path, [(real.plasmid_id, "ACGTXNACGT")])
        with pytest.raises(InsertRecordsInvalid, match="non-ACGT"):
            load_inserts(path)

    def test_conflicting_duplicate_ids_are_named(self, tmp_path):
        real = inventory()[0]
        path = _csv(tmp_path, [(real.plasmid_id, "ACGT"), (real.plasmid_id, "TGCA")])
        with pytest.raises(InsertRecordsInvalid, match="conflicting duplicate"):
            load_inserts(path)

    def test_a_table_that_is_not_this_kit_is_refused(self, tmp_path):
        path = _csv(tmp_path, [("not-a-grasp-plasmid", "ACGTACGT")])
        with pytest.raises(InsertRecordsInvalid, match="does not look like"):
            load_inserts(path)

    def test_length_disagreement_with_the_index_is_refused(self, tmp_path):
        """The check that catches a file which parses cleanly but is the wrong dataset."""
        real = inventory()[0]
        path = _csv(tmp_path, [(real.plasmid_id, "ACGT" * 10)])
        with pytest.raises(InsertRecordsInvalid, match="insert length disagrees"):
            load_inserts(path)

    def test_a_bad_overlap_refuses_to_join(self, kit):
        sequences, source = kit
        plan = select("AAAAUGUGG")
        broken = dict(sequences)
        second = plan.modules[1].plasmid_id
        broken[second] = "TTTT" + broken[second][4:]
        with pytest.raises(OverlapMismatch, match="cannot be joined"):
            reconstruct(plan, broken, source)

    def test_a_missing_module_sequence_is_unavailable(self, kit):
        sequences, source = kit
        plan = select("AAAAUGUGG")
        without = {k: v for k, v in sequences.items() if k != plan.modules[2].plasmid_id}
        with pytest.raises(InsertsUnavailable, match="module selection stands"):
            reconstruct(plan, without, source)

    def test_a_target_with_no_kit_route_is_refused(self, kit):
        sequences, source = kit
        with pytest.raises(InsertsUnavailable, match="no kit route"):
            reconstruct(select("ACGUACGUACGU"), sequences, source)


class TestReconstruction:
    @pytest.mark.parametrize("target", SIX)
    def test_the_six_targets_reconstruct(self, kit, target):
        sequences, source = kit
        plan = select(target)
        rec = reconstruct(plan, sequences, source)
        assert len(rec.final.sequence) == EXPECTED_NT[len(target)]
        assert rec.final.modules == tuple(m.plasmid_id for m in plan.modules)
        assert len(rec.final.junctions) == len(plan.modules) - 1
        assert len(rec.stages) == len(plan.sub_assemblies)

    def test_every_junction_is_where_it_says_it_is(self, kit):
        """Coordinates must index the product, or an exported result cannot be checked."""
        sequences, source = kit
        rec = reconstruct(select("UUACACGUG"), sequences, source)
        for j in rec.final.junctions:
            assert rec.final.sequence[j.start:j.start + 4] == j.overhang

    def test_independent_rejoin_from_raw_records_agrees(self, kit):
        """Re-derive the product without the helper under test, from the exported recipe."""
        sequences, source = kit
        plan = select("AUCGAUCGA")
        rec = reconstruct(plan, sequences, source)

        rebuilt = sequences[rec.final.modules[0]]
        for pid in rec.final.modules[1:]:
            nxt = sequences[pid]
            assert rebuilt[-4:] == nxt[:4]
            rebuilt += nxt[4:]
        assert rebuilt == rec.final.sequence


class TestTranslationAndScope:
    @pytest.mark.parametrize("target,aa", [("AAAAUGUGG", 300), ("UUACACGUGCGUAC", 455)])
    def test_frame_comes_from_the_fusion_site_and_has_no_stops(self, kit, target, aa):
        sequences, source = kit
        plan = select(target)
        tr = translate(reconstruct(plan, sequences, source), plan.modules)
        assert tr.frame == 1 and tr.stop_codons == 0 and len(tr.protein) == aa

    def test_the_product_is_an_insert_not_a_cds(self, kit):
        """It opens mid-codon: two bases of the first codon are the backbone's, not ours."""
        sequences, source = kit
        plan = select("AAAAUGUGG")
        tr = translate(reconstruct(plan, sequences, source), plan.modules)
        assert tr.upstream_bases_needed == 2
        assert "acceptor backbone" in tr.note

    def test_scope_is_stated_on_every_reconstruction(self, kit):
        sequences, source = kit
        rec = reconstruct(select("AAAAUGUGG"), sequences, source)
        assert "not an expression construct" in rec.scope

    def test_an_unknown_fusion_site_is_not_guessed(self, kit):
        """A wrong frame yields a plausible protein, so it must refuse rather than search."""
        from dataclasses import replace
        sequences, source = kit
        plan = select("AAAAUGUGG")
        rec = reconstruct(plan, sequences, source)
        odd = (replace(plan.modules[0], five_overhang="TTTT"),) + plan.modules[1:]
        with pytest.raises(OverlapMismatch, match="will not be guessed"):
            translate(rec, odd)


class TestComparison:
    def test_the_compared_regions_are_the_same_length(self, kit):
        """Both routes encode the same protein over the declared region, so both slices
        must cover it; unequal lengths would mean the region was mis-declared."""
        from test_design import TOY_TABLE
        from clippr.design import design_oneshot

        sequences, source = kit
        plan = select("AAAAUGUGG")
        rec = reconstruct(plan, sequences, source)
        tr = translate(rec, plan.modules)
        cds = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE,
                             check_offtarget=False, seed=42)["cds"]
        out = compare_to_synthesis(rec, tr, cds)
        assert out["kit_region_nt"] == out["cds_region_nt"]
        assert out["match"] == "" or len(out["match"]) == out["longest_shared_nt"]


class TestComparisonRefusesMismatches:
    def test_a_different_targets_cds_is_rejected(self, kit):
        """It returned a plausible tract under 'both encode the same protein'."""
        from clippr.design import design_oneshot
        from clippr.products import ScopeMismatch
        from test_design import TOY_TABLE

        sequences, source = kit
        plan = select("CGUACGUAC")
        rec = reconstruct(plan, sequences, source)
        tr = translate(rec, plan.modules)
        wrong = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE,
                               check_offtarget=False, seed=42)["cds"]
        with pytest.raises(ScopeMismatch, match="do not encode the same protein"):
            compare_to_synthesis(rec, tr, wrong)

    def test_a_ragged_region_is_rejected(self, kit):
        """One extra base made 897 nt compare against 898 and still report agreement,
        because both were truncated to whole codons before translating."""
        from clippr.design import design_oneshot
        from clippr.products import ScopeMismatch
        from test_design import TOY_TABLE

        sequences, source = kit
        plan = select("CGUACGUAC")
        rec = reconstruct(plan, sequences, source)
        tr = translate(rec, plan.modules)
        cds = design_oneshot("CGUACGUAC", codon_table=TOY_TABLE,
                             check_offtarget=False, seed=42)["cds"]
        with pytest.raises(ScopeMismatch, match="not comparable"):
            compare_to_synthesis(rec, tr, cds + "A")

    def test_the_matching_target_reports_its_verified_protein(self, kit):
        from clippr.design import design_oneshot
        from test_design import TOY_TABLE

        sequences, source = kit
        plan = select("CGUACGUAC")
        rec = reconstruct(plan, sequences, source)
        tr = translate(rec, plan.modules)
        right = design_oneshot("CGUACGUAC", codon_table=TOY_TABLE,
                               check_offtarget=False, seed=42)["cds"]
        out = compare_to_synthesis(rec, tr, right)
        assert out["region"]["common_protein_aa"] > 0
        assert "Verified" in out["region"]["basis"]


class TestExport:
    def test_export_carries_provenance_and_coordinates(self, kit, tmp_path):
        sequences, source = kit
        plan = select("AAAAUGUGG")
        rec = reconstruct(plan, sequences, source)
        tr = translate(rec, plan.modules)
        payload = json.loads(write(rec, tmp_path / "kit.json",
                                   translation=tr).read_text(encoding="utf-8"))
        assert payload["source"]["sha256"] == source.sha256
        assert payload["translation"]["frame"] == 1
        assert payload["final"]["junctions"][0]["overhang"]
        assert "not an expression construct" in payload["scope"]

    def test_fasta_repeats_the_scope_on_every_record(self, kit, tmp_path):
        """A FASTA travels without its JSON, so each header must carry the caveat itself."""
        from clippr.products import write_fasta

        sequences, source = kit
        plan = select("AAAAUGUGG")
        rec = reconstruct(plan, sequences, source)
        text = write_fasta(rec, tmp_path / "kit.fasta").read_text(encoding="utf-8")
        headers = [l for l in text.splitlines() if l.startswith(">")]
        assert len(headers) == len(rec.stages) + 1
        assert all("not an expression construct" in h for h in headers)
        assert all(source.sha256[:16] in h for h in headers)

    def test_fasta_final_record_is_the_product(self, kit, tmp_path):
        from clippr.products import write_fasta

        sequences, source = kit
        plan = select("AAAAUGUGG")
        rec = reconstruct(plan, sequences, source)
        text = write_fasta(rec, tmp_path / "kit.fasta").read_text(encoding="utf-8")
        blocks = text.split(">")[1:]
        final = "".join(blocks[-1].splitlines()[1:])
        assert final == rec.final.sequence

    def test_genbank_cds_translation_matches_its_own_span(self, kit, tmp_path):
        """The declared translation must be derivable from the file, not just asserted in it.

        A GenBank whose /translation disagrees with the bases under its CDS feature is worse
        than one with no CDS at all: it looks authoritative and is wrong.
        """
        from Bio import SeqIO
        from Bio.Seq import Seq

        from clippr.products import write_genbank

        sequences, source = kit
        plan = select("AAAAUGUGG")
        rec = reconstruct(plan, sequences, source)
        tr = translate(rec, plan.modules)
        path = write_genbank(rec, tmp_path / "kit.gb", translation=tr)
        record = SeqIO.read(path, "genbank")

        assert str(record.seq) == rec.final.sequence
        cds = [f for f in record.features if f.type == "CDS"]
        assert len(cds) == 1
        assert str(Seq(cds[0].extract(record.seq)).translate()) == \
            cds[0].qualifiers["translation"][0] == tr.protein

    def test_genbank_module_spans_overlap_by_the_overhang(self, kit, tmp_path):
        """Consecutive modules share their four junction bases. That is the assembly."""
        from Bio import SeqIO

        from clippr.products import write_genbank

        sequences, source = kit
        plan = select("AAAAUGUGG")
        rec = reconstruct(plan, sequences, source)
        path = write_genbank(rec, tmp_path / "kit.gb")
        record = SeqIO.read(path, "genbank")

        modules = [f for f in record.features if f.type == "misc_feature"]
        assert len(modules) == len(rec.final.modules)
        for left, right in zip(modules, modules[1:]):
            assert int(left.location.end) - int(right.location.start) == 4
        junctions = [f for f in record.features if f.type == "misc_binding"]
        assert len(junctions) == len(rec.final.junctions)
        for feature, junction in zip(junctions, rec.final.junctions):
            assert str(feature.extract(record.seq)) == junction.overhang

    def test_genbank_without_a_translation_omits_the_cds(self, kit, tmp_path):
        """No frame was supplied, so none is invented."""
        from Bio import SeqIO

        from clippr.products import write_genbank

        sequences, source = kit
        plan = select("AAAAUGUGG")
        rec = reconstruct(plan, sequences, source)
        record = SeqIO.read(write_genbank(rec, tmp_path / "kit.gb"), "genbank")
        assert not [f for f in record.features if f.type == "CDS"]
        assert "not an expression construct" in record.annotations["comment"]
