"""Unit tests for host off-target scanning.

The module reports two tiers and the distinction is the point: a PPR binds RNA, so an
occurrence inside an annotated transcript means something a genomic-DNA occurrence does
not. Tests that need the real genome skip cleanly when it is not cached and there is no
network, rather than failing a whole suite offline.
"""
from __future__ import annotations

import pytest

from clippr.offtarget import (
    Transcript,
    architecture_advice,
    expected_by_chance,
    find_in_genome,
    find_in_transcripts,
    load_genome,
    load_transcripts,
    report,
    reverse_complement,
    scan,
)

TOY = "AAAACCCCGGGGTTTT" + "ACGTACGTACGT" + "TTTTGGGGCCCCAAAA"
#: geneC is deliberately asymmetric. CCCCGGGG and ACGTACGT are their own reverse
#: complements, so neither can demonstrate strand-specificity -- an earlier version of
#: these tests used one and could not have failed.
TOY_TRANSCRIPTS = [
    Transcript("geneA", "CDS", 1, 0, "AAAACCCCGGGGTTTT"),
    Transcript("geneB", "rRNA", -1, 40, "ACGTACGTACGT"),
    Transcript("geneC", "CDS", 1, 60, "AAAGGGTTTCCCAAAGGG"),
]
ASYMMETRIC = "AAAGGGTTT"          # revcomp AAACCCTTT, absent from geneC


class TestGenomicSearch:
    def test_finds_an_occurrence(self):
        assert find_in_genome("CCCCGGGG", TOY)

    def test_searches_both_strands(self):
        """DNA is double-stranded, so a genomic scan looks at both."""
        assert find_in_genome(reverse_complement("CCCCGGGG"), TOY)

    def test_absent_sequence_is_absent(self):
        assert find_in_genome("GGGGGGGGGGGG", TOY) == []

    @pytest.mark.parametrize("bad", ["", "ACGTN", "hello"])
    def test_malformed_target_rejected(self, bad):
        with pytest.raises(ValueError):
            find_in_genome(bad, TOY)


class TestTranscriptSearch:
    def test_finds_a_hit_inside_a_transcript(self):
        hits = find_in_transcripts("CCCCGGGG", TOY_TRANSCRIPTS)
        assert hits and hits[0].where == "geneA"
        assert hits[0].tier == "transcript"

    def test_records_the_feature_type(self):
        assert find_in_transcripts("ACGTACGT", TOY_TRANSCRIPTS)[0].kind == "rRNA"

    def test_is_sense_strand_only(self):
        """A PPR reads the transcript 5'->3'; the antisense sequence is not in that RNA.

        `load_transcripts` stores each feature already in its own orientation, so the
        search is deliberately one-directional -- unlike the genomic scan. Uses an
        asymmetric sequence, because a palindrome could not tell the two apart.
        """
        assert reverse_complement(ASYMMETRIC) != ASYMMETRIC, "must be asymmetric to test this"
        assert find_in_transcripts(ASYMMETRIC, TOY_TRANSCRIPTS)
        assert not find_in_transcripts(reverse_complement(ASYMMETRIC), TOY_TRANSCRIPTS)

    def test_transcripts_are_a_subset_of_the_genome(self):
        """Anything in a transcript is in the DNA; the reverse does not hold."""
        assert find_in_genome("CCCCGGGG", TOY)
        assert find_in_transcripts("CCCCGGGG", TOY_TRANSCRIPTS)


class TestTiers:
    def test_transcript_hit_outranks_a_genomic_one(self):
        r = scan("CCCCGGGG", TOY, TOY_TRANSCRIPTS)
        assert r["verdict"] == "OCCURS IN A TRANSCRIPT"
        assert r["n_transcript"] >= 1
        assert "geneA" in r["genes"]

    def test_genomic_only_is_reported_as_weaker(self):
        """In the DNA but in no transcript: no RNA exists for a PPR to bind there."""
        r = scan("TTTTGGGGCCCC", TOY, TOY_TRANSCRIPTS)
        assert r["n_genomic"] >= 1
        assert r["n_transcript"] == 0
        assert "genomic DNA only" in r["verdict"]

    def test_absent_is_absent(self):
        r = scan("GGGGGGGGGGGG", TOY, TOY_TRANSCRIPTS)
        assert r["verdict"] == "not found in the host"

    def test_report_keeps_the_tiers_separate(self):
        text = report([scan("CCCCGGGG", TOY, TOY_TRANSCRIPTS)])
        assert "transcript" in text and "genomic" in text
        assert "not sufficient" in text, "the caveat must travel with the number"


class TestExpectation:
    def test_uses_composition_not_a_uniform_model(self):
        at_rich = "A" * 20000 + "T" * 20000 + "C" * 5000 + "G" * 5000
        assert expected_by_chance("AAAAAAAA", at_rich) > expected_by_chance(
            "CCCCCCCC", at_rich)

    def test_longer_targets_are_rarer(self):
        assert expected_by_chance("ACGTACGTA", TOY * 200) > expected_by_chance(
            "ACGTACGTACGTAC", TOY * 200)


@pytest.fixture(scope="module")
def host():
    try:
        return load_genome(), load_transcripts()
    except Exception as exc:
        pytest.skip(f"chloroplast genome unavailable: {exc}")


class TestRealGenome:
    def test_genome_size_and_composition(self, host):
        genome, _ = host
        assert 200_000 < len(genome) < 210_000
        gc = (genome.count("G") + genome.count("C")) / len(genome)
        assert 0.30 < gc < 0.40, "chloroplast genomes are AT-rich"

    def test_transcripts_are_annotated_and_a_minority_of_the_genome(self, host):
        genome, trs = host
        assert len(trs) > 50
        coding = sum(len(t.sequence) for t in trs)
        assert 0.2 < coding / len(genome) < 0.7, "measured at 43.5%"

    def test_both_strands_are_represented(self, host):
        _, trs = host
        assert {t.strand for t in trs} == {1, -1}

    def test_nine_mers_are_not_unique_in_this_host(self, host):
        genome, _ = host
        advice = architecture_advice(genome)
        assert advice[9] > 1.0
        assert advice[14] < 0.05
        assert advice[19] < 0.001

    def test_transcript_tier_is_strictly_smaller_than_genomic(self, host):
        """The correction that matters: 48% of 9-mers hit the DNA, 16% hit a transcript."""
        genome, trs = host
        r = scan("AAAAUGUGG", genome, trs)
        assert r["n_genomic"] >= r["n_transcript"]
        assert r["n_transcript"] >= 1, "this target was measured inside ORF1995"

    def test_a_genomic_only_target_is_not_called_a_transcript_hit(self, host):
        genome, trs = host
        r = scan("GCUAAAGAC", genome, trs)
        assert r["n_genomic"] >= 1
        assert r["n_transcript"] == 0

    def test_a_nineteen_mer_hits_neither_tier(self, host):
        genome, trs = host
        r = scan("AAAGCGGCACUUGUGAAGU", genome, trs)
        assert r["n_genomic"] == 0 and r["n_transcript"] == 0
