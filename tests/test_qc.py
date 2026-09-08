"""Unit tests for synthesis QC.

The point of this module is that its verdict *varies*, so the tests are built around
making each state happen deliberately rather than around one golden sequence.
"""
from __future__ import annotations

import pytest

from clippr.qc import DEFAULT_THRESHOLDS, REPEAT_K, synthesis_qc

import random


def _clean_sequence(length: int = 900, gc: str = "ACGT") -> str:
    """A random sequence that is genuinely clean.

    A 900-nt random draw carries a blacklisted enzyme site often enough that picking one
    seed and hoping is not good enough -- the first attempt here contained a BbsI site.
    Reject and redraw until the sequence is actually free of them.
    """
    for seed in range(200):
        rng = random.Random(seed)
        seq = "".join(rng.choice(gc) for _ in range(length))
        if not synthesis_qc(seq)["blacklist_hits"]:
            return seq
    raise AssertionError("could not draw a blacklist-free sequence")


#: A clean sequence: mixed composition, no long run, no repeated 20-mer, no enzyme site.
CLEAN = _clean_sequence()

#: A 30-nt unit tiled 30 times -- the failure mode a tandem-repeat protein produces.
REPETITIVE = "ACGTTGCAAGGCTTACGGATCCATGCATG" * 31


def test_clean_sequence_passes():
    assert synthesis_qc(CLEAN)["status"] == "PASS"


def test_repetitive_sequence_fails():
    q = synthesis_qc(REPETITIVE)
    assert q["status"] == "FAIL"
    assert q["repeated_kmer_fraction"] > 0.5


class TestAllStatesReachable:
    """A flag that cannot reach every state is the defect this module was written to fix."""

    def test_pass(self):
        assert synthesis_qc(CLEAN)["status"] == "PASS"

    def test_warning(self):
        # one duplicated stretch of ~30 nt, but most of the sequence unique
        seq = CLEAN[:400] + CLEAN[:30] + CLEAN[400:]
        q = synthesis_qc(seq)
        assert q["status"] == "WARNING", q

    def test_fail(self):
        assert synthesis_qc(REPETITIVE)["status"] == "FAIL"


class TestRepeats:
    def test_no_repeat_in_random_sequence(self):
        q = synthesis_qc(CLEAN)
        assert q["repeated_kmer_fraction"] == 0.0
        assert q["repeated_kmers"] == 0

    def test_duplicated_kmer_is_counted(self):
        seq = CLEAN[:300] + CLEAN[:REPEAT_K] + CLEAN[300:]
        assert synthesis_qc(seq)["repeated_kmers"] >= 1

    def test_longest_repeat_is_measured(self):
        chunk = CLEAN[:40]
        q = synthesis_qc(CLEAN[:300] + chunk + CLEAN[300:])
        assert q["longest_repeat"] >= 40

    def test_repeat_fraction_is_bounded(self):
        for seq in (CLEAN, REPETITIVE):
            assert 0.0 <= synthesis_qc(seq)["repeated_kmer_fraction"] <= 1.0

    def test_fraction_captures_spread_not_just_the_worst_stretch(self):
        """Two sequences can share a longest repeat yet differ in how much is repetitive."""
        one_long = CLEAN[:200] + CLEAN[:60] + CLEAN[200:]
        spread = synthesis_qc(REPETITIVE)["repeated_kmer_fraction"]
        assert spread > synthesis_qc(one_long)["repeated_kmer_fraction"]


class TestBlacklist:
    def test_enzyme_site_is_a_hard_failure(self):
        seq = CLEAN[:300] + "GGTCTC" + CLEAN[300:]
        q = synthesis_qc(seq)
        assert q["status"] == "FAIL"
        assert any("BsaI" in h for h in q["blacklist_hits"])

    def test_reverse_strand_site_is_caught(self):
        seq = CLEAN[:300] + "GAGACC" + CLEAN[300:]
        assert synthesis_qc(seq)["blacklist_hits"]

    def test_clean_sequence_has_no_hits(self):
        assert synthesis_qc(CLEAN)["blacklist_hits"] == []


class TestFeasibilityGates:
    def test_extreme_gc_fails(self):
        assert synthesis_qc("GC" * 450)["status"] == "FAIL"

    def test_low_gc_fails(self):
        assert synthesis_qc("AT" * 450)["status"] == "FAIL"

    @staticmethod
    def _with_run(n: int) -> str:
        """Insert a run of exactly `n` A's, flanked so it cannot merge with neighbours.

        Splicing into a random sequence is not enough: the first attempt here landed
        beside an existing AA and produced a 9-run where 6 was intended.
        """
        return CLEAN[:300] + "C" + "A" * n + "C" + CLEAN[300:]

    def test_long_homopolymer_is_flagged(self):
        q = synthesis_qc(self._with_run(10))
        assert q["longest_homopolymer"] == 10
        assert q["status"] == "FAIL"

    def test_moderate_homopolymer_warns_not_fails(self):
        q = synthesis_qc(self._with_run(6))
        assert q["longest_homopolymer"] == 6
        assert q["status"] == "WARNING"
        assert "homopolymer" in " ".join(q["warnings"])

    def test_chlamydomonas_gc_is_not_penalised(self):
        """~60% GC is normal for this host and must not be treated as a defect.

        Our optimiser measures 58.8-60.6% on the corpus, and `codons.optimize_cds`
        constrains GC to 65%, so the warn line sits just above what the host produces.
        """
        seq = _clean_sequence(gc="GCGCGCATAT")        # 60% G/C by construction
        q = synthesis_qc(seq)
        assert 55 <= q["gc_pct"] <= 65
        assert not any("GC" in w for w in q["warnings"] + q["failures"])


class TestThresholds:
    def test_defaults_are_documented_keys(self):
        for k in ("longest_repeat_warn", "longest_repeat_fail",
                  "repeat_fraction_warn", "repeat_fraction_fail"):
            assert k in DEFAULT_THRESHOLDS

    def test_overrides_take_effect(self):
        strict = synthesis_qc(CLEAN, {"repeat_fraction_warn": -1})
        assert strict["status"] != "PASS"

    def test_partial_override_keeps_the_rest(self):
        q = synthesis_qc(REPETITIVE, {"longest_repeat_fail": 10_000})
        assert q["status"] == "FAIL"      # still fails on repeat_fraction

    def test_loosening_everything_passes_a_bad_sequence(self):
        q = synthesis_qc(REPETITIVE, {"longest_repeat_warn": 10_000,
                                      "longest_repeat_fail": 10_000,
                                      "repeat_fraction_warn": 2,
                                      "repeat_fraction_fail": 2})
        assert q["status"] == "PASS"


class TestInput:
    def test_rna_spelling_accepted(self):
        assert synthesis_qc("ACGU" * 100)["length_nt"] == 400

    def test_empty_rejected(self):
        with pytest.raises(ValueError, match="empty"):
            synthesis_qc("")

    def test_ambiguous_bases_rejected(self):
        with pytest.raises(ValueError, match="non-ACGT"):
            synthesis_qc("ACGTN" * 20)

    def test_every_documented_key_is_returned(self):
        q = synthesis_qc(CLEAN)
        for key in ("status", "longest_homopolymer", "homopolymer_count",
                    "repeated_kmer_fraction", "gc_pct", "gc_window_min", "gc_window_max",
                    "blacklist_hits", "warnings", "failures"):
            assert key in q

    def test_reasons_accompany_every_non_pass(self):
        for seq in (REPETITIVE, "GC" * 450):
            q = synthesis_qc(seq)
            assert q["warnings"] or q["failures"]
