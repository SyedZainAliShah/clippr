"""Unit tests for serialisation: CSV, FASTA, GenBank and the cost estimate."""
from __future__ import annotations

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from clippr.assembly import build_oligos, split_cds
from clippr.export import (
    OPOOL_LIST_PRICE_EUR,
    PHOSPHORYLATION_EUR_PER_OLIGO,
    opool_quote,
    write_fasta,
    write_gene_fasta,
    write_genbank,
    write_oligo_csv,
)
from clippr.ppr import describe

TARGET = "AAAAUGUGG"


@pytest.fixture(scope="module")
def design():
    protein = describe(TARGET)["aa_sequence"]
    # a plain back-translation is enough; these tests are about serialising, not design
    from clippr import constants as C
    cds = "".join(C.SYNONYMOUS_CODONS[a][0] for a in protein)
    cuts = [88, 171, 238]
    frags = split_cds(cds, cuts)
    return cds, frags, build_oligos(cds, cuts, prefix=TARGET)


class TestOligoCsv:
    def test_round_trips(self, design, tmp_path):
        import pandas as pd
        _, _, oligos = design
        p = write_oligo_csv(oligos, tmp_path / "o.csv")
        assert list(pd.read_csv(p)["oligo_sequence_5to3"]) == list(
            oligos["oligo_sequence_5to3"])

    def test_creates_missing_directories(self, design, tmp_path):
        _, _, oligos = design
        p = write_oligo_csv(oligos, tmp_path / "deep" / "er" / "o.csv")
        assert p.exists()


class TestFasta:
    def test_one_record_per_oligo(self, design, tmp_path):
        _, _, oligos = design
        recs = list(SeqIO.parse(write_fasta(oligos, tmp_path / "o.fa"), "fasta"))
        assert len(recs) == len(oligos)

    def test_sequences_match_the_table(self, design, tmp_path):
        _, _, oligos = design
        recs = list(SeqIO.parse(write_fasta(oligos, tmp_path / "o.fa"), "fasta"))
        assert [str(r.seq) for r in recs] == list(oligos["oligo_sequence_5to3"])

    def test_header_carries_enzyme_and_overhangs(self, design, tmp_path):
        _, _, oligos = design
        first = next(SeqIO.parse(write_fasta(oligos, tmp_path / "o.fa"), "fasta"))
        assert "BsaI" in first.id and ".." in first.id

    def test_wrapped_lines_do_not_corrupt_the_sequence(self, design, tmp_path):
        """60-column wrapping must not drop or duplicate bases."""
        _, _, oligos = design
        recs = list(SeqIO.parse(write_fasta(oligos, tmp_path / "o.fa"), "fasta"))
        assert len(str(recs[0].seq)) == oligos.iloc[0]["oligo_length"]

    def test_gene_fasta_round_trips(self, design, tmp_path):
        cds, _, _ = design
        rec = next(SeqIO.parse(write_gene_fasta(cds, TARGET, tmp_path / "g.fa"), "fasta"))
        assert str(rec.seq) == cds


@pytest.fixture(scope="module")
def record(design, tmp_path_factory):
    cds, frags, _ = design
    p = write_genbank(cds, frags, TARGET, tmp_path_factory.mktemp("gb") / "d.gb")
    return SeqIO.read(p, "genbank")


class TestGenBank:
    def test_sequence_round_trips(self, design, record):
        cds, _, _ = design
        assert str(record.seq) == cds

    def test_single_cds_spanning_everything(self, design, record):
        cds, _, _ = design
        feats = [f for f in record.features if f.type == "CDS"]
        assert len(feats) == 1
        assert int(feats[0].location.start) == 0
        assert int(feats[0].location.end) == len(cds)

    def test_translation_qualifier_is_correct(self, design, record):
        cds, _, _ = design
        feats = [f for f in record.features if f.type == "CDS"]
        assert feats[0].qualifiers["translation"][0] == str(Seq(cds).translate())

    def test_every_feature_lies_inside_the_record(self, design, record):
        cds, _, _ = design
        for f in record.features:
            assert 0 <= int(f.location.start) <= int(f.location.end) <= len(cds)

    def test_one_feature_per_repeat(self, record):
        d = describe(TARGET)
        labels = [f.qualifiers.get("label", [""])[0] for f in record.features]
        for i, base in enumerate(d["target_rna"], start=1):
            assert any(l.startswith(f"PPR{i} {base}") for l in labels)

    def test_specificity_residues_match_the_code(self, design, record):
        """The annotation must agree with the protein, not merely be well-formed."""
        cds, _, _ = design
        protein = str(Seq(cds).translate())
        d = describe(TARGET)
        by_label = {f.qualifiers.get("label", [""])[0]: f for f in record.features}
        for i, code in enumerate(d["code_pairs"], start=1):
            for which, expect in (("5th", code[0]), ("last", code[1])):
                f = by_label[f"PPR{i} {which} {expect}"]
                assert protein[int(f.location.start) // 3] == expect

    def test_one_feature_per_fragment(self, design, record):
        _, frags, _ = design
        labels = [f.qualifiers.get("label", [""])[0] for f in record.features]
        for fr in frags:
            assert f"fragment {fr.fragment_id}" in labels

    def test_destination_noted_as_outside_the_record(self, record):
        assert "outside this record" in record.annotations["comment"]


class TestQuote:
    def test_counts_and_lengths(self, design):
        _, _, oligos = design
        q = opool_quote(oligos)
        assert q["n_oligos"] == len(oligos)
        assert q["total_bases"] == sum(oligos["oligo_length"])
        assert q["min_oligo_nt"] == min(oligos["oligo_length"])
        assert q["max_oligo_nt"] == max(oligos["oligo_length"])

    def test_pool_price_is_flat(self, design):
        """Measured constant across 200 designs; cost cannot discriminate between them."""
        _, _, oligos = design
        assert opool_quote(oligos)["dna_eur"] == OPOOL_LIST_PRICE_EUR
        assert opool_quote(oligos.head(2))["dna_eur"] == OPOOL_LIST_PRICE_EUR

    def test_phosphorylation_is_per_oligo(self, design):
        _, _, oligos = design
        q = opool_quote(oligos)
        assert q["phospho_eur"] == pytest.approx(
            PHOSPHORYLATION_EUR_PER_OLIGO * len(oligos), abs=0.01)

    def test_phosphorylation_excluded_unless_requested(self, design):
        _, _, oligos = design
        assert opool_quote(oligos)["total_eur"] == OPOOL_LIST_PRICE_EUR
        withp = opool_quote(oligos, phosphorylate_5prime=True)
        assert withp["total_eur"] > OPOOL_LIST_PRICE_EUR

    def test_flagged_as_not_a_quote(self, design):
        _, _, oligos = design
        assert opool_quote(oligos)["list_price_not_a_quote"] is True

    def test_empty_rejected(self, design):
        _, _, oligos = design
        with pytest.raises(ValueError, match="no oligos"):
            opool_quote(oligos.head(0))
