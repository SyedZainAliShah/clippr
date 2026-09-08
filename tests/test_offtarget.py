"""Unit tests for host off-target scanning.

The genome is fetched and cached on first use, so tests that need it skip cleanly when
there is no cached copy and no network rather than failing a whole suite offline.
"""
from __future__ import annotations

import pytest

from clippr.offtarget import (
    Hit,
    architecture_advice,
    expected_by_chance,
    find_hits,
    load_genome,
    report,
    reverse_complement,
    scan,
)

#: A tiny synthetic genome, so the search logic is tested without a network call.
TOY = "AAAACCCCGGGGTTTT" + "ACGTACGTACGT" + "TTTTGGGGCCCCAAAA"


class TestSearch:
    def test_finds_an_exact_occurrence(self):
        hits = find_hits("CCCCGGGG", TOY)
        assert any(h.strand == "+" and h.mismatches == 0 for h in hits)

    def test_finds_the_reverse_strand(self):
        """A chloroplast transcript may come from either strand."""
        probe = reverse_complement("CCCCGGGG")
        hits = find_hits(probe, TOY)
        assert any(h.strand == "-" for h in hits)

    def test_absent_sequence_is_absent(self):
        assert find_hits("GGGGGGGGGGGG", TOY) == []

    def test_mismatch_tolerance_widens_the_search(self):
        exact = find_hits("ACGTACGT", TOY, 0)
        loose = find_hits("ACGTACGA", TOY, 1)
        assert len(loose) > len([h for h in find_hits("ACGTACGA", TOY, 0)])
        assert exact

    def test_mismatch_count_is_reported(self):
        hits = find_hits("ACGTACGA", TOY, 1)
        assert all(h.mismatches <= 1 for h in hits)
        assert any(h.mismatches == 1 for h in hits)

    def test_context_shows_the_match_in_upper_case(self):
        h = find_hits("CCCCGGGG", TOY)[0]
        assert "CCCCGGGG" in h.context
        assert any(c.islower() for c in h.context), "flanks should be lower case"

    def test_fast_and_slow_paths_agree(self):
        """Exact matching takes a str.find path; it must find exactly what the loop does."""
        fast = {(h.strand, h.position) for h in find_hits("ACGT", TOY, 0)}
        slow = {(h.strand, h.position) for h in find_hits("ACGT", TOY, 1)
                if h.mismatches == 0}
        assert fast == slow

    def test_rna_spelling_accepted(self):
        assert find_hits("CCCCGGGG", TOY) == find_hits("CCCCGGGG".replace("T", "U"), TOY)

    @pytest.mark.parametrize("bad", ["", "ACGTN", "hello"])
    def test_malformed_target_rejected(self, bad):
        with pytest.raises(ValueError):
            find_hits(bad, TOY)


class TestExpectation:
    def test_uses_the_genome_composition_not_a_uniform_model(self):
        """An AT-rich genome makes AT-rich targets commoner; a uniform model would hide that."""
        at_rich = "A" * 20000 + "T" * 20000 + "C" * 5000 + "G" * 5000
        assert expected_by_chance("AAAAAAAA", at_rich) > expected_by_chance("CCCCCCCC", at_rich)

    def test_longer_targets_are_rarer(self):
        assert expected_by_chance("ACGTACGTA", TOY * 200) > expected_by_chance(
            "ACGTACGTACGTAC", TOY * 200)

    def test_never_negative(self):
        assert expected_by_chance("ACGTACGTA", TOY * 100) >= 0


class TestScan:
    def test_verdict_reflects_the_hits(self):
        r = scan("CCCCGGGG", TOY, max_mismatches=0)
        assert r["verdict"] == "OCCURS IN HOST"
        assert r["n_exact"] >= 1

    def test_clean_target_says_so(self):
        r = scan("GGGGGGGGGGGG", TOY, max_mismatches=0)
        assert r["verdict"] == "not found in host"
        assert r["n_exact"] == 0

    def test_near_matches_are_separated_from_exact(self):
        r = scan("ACGTACGA", TOY, max_mismatches=1)
        assert all(h.mismatches for h in r["near_hits"])

    def test_reports_every_documented_key(self):
        r = scan("ACGT", TOY)
        for key in ("target", "length", "verdict", "exact_hits", "near_hits",
                    "n_exact", "n_near", "expected_by_chance", "genome_length"):
            assert key in r

    def test_report_renders(self):
        text = report([scan("CCCCGGGG", TOY, 0), scan("GGGGGGGGGGGG", TOY, 0)])
        assert "OCCURS IN HOST" in text and "verdict" in text


@pytest.fixture(scope="module")
def genome():
    try:
        return load_genome()
    except Exception as exc:                       # no cache and no network
        pytest.skip(f"chloroplast genome unavailable: {exc}")


class TestRealGenome:
    def test_genome_is_the_expected_size(self, genome):
        assert 200_000 < len(genome) < 210_000
        assert not set(genome) - set("ACGTN")

    def test_chloroplast_is_at_rich(self, genome):
        gc = (genome.count("G") + genome.count("C")) / len(genome)
        assert 0.30 < gc < 0.40, "chloroplast genomes are AT-rich; this one measured 34.5%"

    def test_nine_mers_are_not_unique_in_this_host(self, genome):
        """The finding this module exists to surface.

        A 9-nt target is expected to occur several times by chance in a 204 kb genome, so
        9S architectures cannot be assumed specific in vivo without checking.
        """
        advice = architecture_advice(genome)
        assert advice[9] > 1.0, "a 9-mer should be expected more than once"
        assert advice[14] < 0.05
        assert advice[19] < 0.001

    def test_length_is_what_buys_specificity(self, genome):
        advice = architecture_advice(genome)
        assert advice[9] > advice[14] > advice[19]

    def test_a_known_corpus_target_occurs_in_the_host(self, genome):
        """AAAAUGUGG was measured at 5 exact occurrences; a regression here matters."""
        r = scan("AAAAUGUGG", genome, max_mismatches=0)
        assert r["n_exact"] >= 1
        assert r["verdict"] == "OCCURS IN HOST"

    def test_a_nineteen_mer_does_not(self, genome):
        r = scan("AAAGCGGCACUUGUGAAGU", genome, max_mismatches=0)
        assert r["n_exact"] == 0
