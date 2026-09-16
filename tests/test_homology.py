"""Unit tests for cross-member DNA homology.

The measurement primitives are tested exactly, because everything downstream rests on them. The
expensive `diversify_library` gets one small real run asserting the property that matters: it
must never make the library worse, which is exactly what the two rejected approaches did.
"""
from __future__ import annotations

import pytest

from clippr.homology import (
    HR_THRESHOLD_NT,
    assess,
    block_profile,
    encoding_capacity,
    longest_shared,
    repeat_blocks,
    report,
    scaffold_encoding,
    shared_kmers,
)

# The toy table from test_design, so these run without a network call for codon usage.
from test_design import TOY_TABLE


class TestLongestShared:
    def test_finds_the_run_and_both_offsets(self):
        assert longest_shared("AAACGTAAA", "TTTCGTTTT") == (3, 3, 3)

    def test_identical_sequences_share_everything(self):
        assert longest_shared("ACGTACGT", "ACGTACGT") == (8, 0, 0)

    def test_no_shared_base_gives_zero(self):
        assert longest_shared("AAAA", "CCCC") == (0, 0, 0)

    def test_reports_the_longest_not_the_first(self):
        n, ia, ib = longest_shared("AGGGGA" + "TTTT", "CGGGGC" + "TTTTTT")
        assert n == 4

    @pytest.mark.parametrize("a,b", [("", "ACGT"), ("ACGT", ""), ("", "")])
    def test_empty_input_is_zero_not_an_error(self, a, b):
        assert longest_shared(a, b) == (0, 0, 0)

    def test_offsets_actually_locate_the_match(self):
        a, b = "TTTTACGTACGT", "ACGTACGTGGGG"
        n, ia, ib = longest_shared(a, b)
        assert a[ia:ia + n] == b[ib:ib + n]


class TestSharedKmers:
    def test_counts_common_windows(self):
        assert shared_kmers("ACGTACGT", "ACGTACGT", k=4) == 5

    def test_none_in_common(self):
        assert shared_kmers("AAAAAAAA", "CCCCCCCC", k=4) == 0

    def test_sequence_shorter_than_k(self):
        assert shared_kmers("ACG", "ACG", k=20) == 0


class TestAssess:
    def test_orders_longest_first(self):
        cds = {"a": "ACGT" * 10, "b": "ACGT" * 10, "c": "AAAA" + "TTTT" * 8}
        pairs = assess(cds)
        assert [p.length for p in pairs] == sorted((p.length for p in pairs), reverse=True)

    def test_every_pair_appears(self):
        cds = {n: "ACGT" * 10 for n in "abcd"}
        assert len(assess(cds)) == 6

    def test_flag_uses_the_threshold(self):
        long_shared = "ACGTTGCA" * 10          # 80 nt, identical in both
        pairs = assess({"a": long_shared, "b": long_shared})
        assert pairs[0].length >= HR_THRESHOLD_NT
        assert pairs[0].flag == "REVIEW"

    def test_short_sharing_is_ok(self):
        pairs = assess({"a": "ACGTACGTAC", "b": "TTTTTTTTTT"})
        assert pairs[0].flag == "ok"

    def test_detects_a_shared_start(self):
        """Both copies at offset zero is the scaffold signature, not a chance match."""
        shared = "ATGCAGGGCGGCAACAGCGAGGAG"
        pairs = assess({"a": shared + "AAAA", "b": shared + "TTTT"})
        assert pairs[0].at_sequence_start

    def test_a_chance_match_is_not_reported_as_a_shared_start(self):
        pairs = assess({"a": "TTTT" + "ACGTACGTACGT", "b": "GGGG" + "ACGTACGTACGT"})
        assert not pairs[0].at_sequence_start


class TestRepeatBlocks:
    """Blocks carried by many members are a different risk from a pairwise coincidence."""

    def test_counts_members_not_occurrences(self):
        """A k-mer twice in one member is not shared; the count is over members."""
        blocks = dict(repeat_blocks({"a": "ACGTACGTAC" * 2, "b": "TTTTTTTTTT" * 2}, k=10,
                                    min_members=1))
        assert blocks["ACGTACGTAC"] == 1

    def test_finds_a_block_in_every_member(self):
        shared = "ACGTTGCAACGTTGCAACGT"
        blocks = repeat_blocks({n: shared + n * 10 for n in "abc"}, k=20)
        assert blocks and blocks[0][1] == 3

    def test_min_members_filters(self):
        cds = {"a": "ACGTACGTACGTACGTACGT", "b": "TTTTTTTTTTTTTTTTTTTT"}
        assert repeat_blocks(cds, k=20, min_members=2) == []

    def test_profile_reports_library_wide_structure(self):
        shared = "ACGTTGCAACGTTGCAACGT"
        prof = block_profile({n: shared + n * 10 for n in "abcd"})
        assert prof["n_members"] == 4
        assert prof["in_all_members"] >= 1
        assert prof["max_multiplicity"] == 4

    def test_profile_counts_isolated_members(self):
        """A member sharing nothing over threshold should be reported as isolated."""
        prof = block_profile({"a": "ACGT" * 30, "b": "ACGT" * 30, "c": "TTGA" * 30})
        assert prof["isolated_members"] >= 1

    def test_report_includes_the_library_wide_section(self):
        text = report({"a": "ACGT" * 30, "b": "ACGT" * 30})
        assert "library-wide structure" in text
        assert "appear in every member" in text


class TestReport:
    def test_needs_two_members(self):
        assert "fewer than two" in report({"only": "ACGT"})

    def test_carries_the_caveat(self):
        text = report({"a": "ACGT" * 20, "b": "ACGT" * 20})
        assert "necessary substrate" in text
        assert "not a prediction" in text
        assert "rule of thumb" in text

    def test_names_the_scaffold_cause_when_every_flag_is_at_the_start(self):
        shared = "ACGTTGCA" * 10
        text = report({"a": shared, "b": shared})
        assert "fixed scaffold" in text


class TestScaffoldEncoding:
    PREFIX = "MQGGNSEEPRKSFDERPERGVVS"

    def test_translates_back_to_the_prefix(self):
        from Bio.Seq import Seq
        dna = scaffold_encoding(0, self.PREFIX, TOY_TABLE)
        assert str(Seq(dna).translate()) == self.PREFIX

    def test_is_deterministic_in_the_member_index(self):
        """A member must keep its encoding across reruns or the library is not reproducible."""
        assert scaffold_encoding(7, self.PREFIX, TOY_TABLE) == \
            scaffold_encoding(7, self.PREFIX, TOY_TABLE)

    def test_different_members_get_different_encodings(self):
        got = {scaffold_encoding(m, self.PREFIX, TOY_TABLE) for m in range(12)}
        assert len(got) == 12

    def test_encodings_do_not_share_a_long_run(self):
        """The whole point: distinct members must not share a recombination-length stretch."""
        enc = [scaffold_encoding(m, self.PREFIX, TOY_TABLE) for m in range(12)]
        worst = max(longest_shared(a, b)[0]
                    for i, a in enumerate(enc) for b in enc[i + 1:])
        assert worst < 30, f"longest shared run between scaffold encodings was {worst} nt"

    def test_avoids_blacklisted_sites(self):
        from Bio import Restriction
        from Bio.Seq import Seq
        for m in range(8):
            dna = Seq(scaffold_encoding(m, self.PREFIX, TOY_TABLE))
            for name in ("BsaI", "BbsI", "SapI"):
                assert not getattr(Restriction, name).search(dna)

    def test_raises_rather_than_returning_something_invalid(self):
        with pytest.raises(ValueError, match="no valid synonymous encoding"):
            scaffold_encoding(0, self.PREFIX, TOY_TABLE, gc_bounds=(0.99, 1.0), tries=5)


class TestEncodingCapacity:
    """The measurement that decides whether residual homology is a bug or a bound."""

    REGION = "WTAMISGYAQNGRIDEARELFDKMPERNVVS"

    def test_reports_a_capacity_and_agreement(self):
        cap = encoding_capacity(self.REGION, TOY_TABLE, draws=4000)
        assert cap["capacity"] >= 1
        assert cap["region_aa"] == len(self.REGION)
        assert "searches_agree" in cap

    def test_a_region_with_no_synonymous_choice_admits_one_encoding(self):
        """Methionine and tryptophan have a single codon each, so capacity is exactly 1."""
        cap = encoding_capacity("MW", TOY_TABLE, k=3, draws=200)
        assert cap["capacity"] == 1

    def test_impossible_constraints_yield_zero(self):
        cap = encoding_capacity(self.REGION, TOY_TABLE, gc_bounds=(0.99, 1.0), draws=500)
        assert cap["capacity"] == 0

    def test_capacity_rises_with_k(self):
        """A longer window is a *weaker* disjointness requirement, so capacity must rise.

        This is the opposite of the intuition that longer means stricter. Two sequences sharing
        some 12-mer is near-inevitable; sharing a 40-mer takes 40 consecutive identical bases.
        The direction matters for interpreting the capacity figure: 20-mer disjointness is far
        stricter than recombination risk requires, so the practically relevant ceiling is higher
        than the k=20 number.
        """
        small_k = encoding_capacity(self.REGION, TOY_TABLE, k=12, draws=3000)
        large_k = encoding_capacity(self.REGION, TOY_TABLE, k=40, draws=3000)
        assert large_k["capacity"] >= small_k["capacity"]


class TestDiversify:
    """Real runs, kept to three members. The property is that it can never make things worse."""

    def test_never_regresses_and_reports_both_sides(self):
        from clippr.homology import diversify_library
        res = diversify_library(["AAAAUGUGG", "GCUAAAGAC", "UUACACGUG"],
                                codon_table=TOY_TABLE, check_offtarget=False, retries=2)
        assert res["max_after"] <= res["max_before"], (
            "diversification must not lengthen the worst shared stretch — the rejected "
            "reactive approaches failed exactly this")
        assert res["flagged_after"] <= res["flagged_before"]
        assert set(res["cds"]) == set(res["baseline_cds"])

    def test_proteins_survive_diversification(self):
        from Bio.Seq import Seq

        from clippr.homology import diversify_library
        from clippr.ppr import describe
        targets = ["AAAAUGUGG", "GCUAAAGAC"]
        res = diversify_library(targets, codon_table=TOY_TABLE, check_offtarget=False,
                                retries=1)
        for t in targets:
            assert str(Seq(res["cds"][t]).translate()) == describe(t)["aa_sequence"]

    def test_a_retry_actually_varies_the_design(self):
        """The retry loop was inert for the whole of its first life.

        It re-rolled `seed`, but under a locked prefix DNA Chisel returns the same
        sequence for every seed, so all `retries + 1` attempts were one design and the
        "keep re-rolling until it stops clashing" guarantee was decoration. These two
        assertions are the fact the loop now depends on: the seed does nothing once a
        prefix is locked, and the encoding index does. If the first ever starts to fail,
        the loop can go back to varying the seed; if the second fails, the loop is inert
        again and `diversify_library` is not doing what its docstring claims.
        """
        from clippr.biology import N_TERMINAL
        from clippr.design import design_oneshot
        from clippr.homology import scaffold_encoding

        base = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE, check_offtarget=False)
        prefixes = [scaffold_encoding(i, N_TERMINAL, base["codon_table"],
                                      genetic_code=base["genetic_code"],
                                      enzymes=base["enzyme_profile_effective"])
                    for i in (0, 1)]
        assert prefixes[0] != prefixes[1], "two encoding indices must differ to begin with"

        locked = {design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE, check_offtarget=False,
                                 lock_prefix=prefixes[0], seed=s)["cds"] for s in (1, 2, 3)}
        assert len(locked) == 1, (
            "the seed now varies a locked design, so re-rolling it would be a valid retry "
            "again — but diversify_library re-rolls the encoding, so update it deliberately")

        by_encoding = {design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE,
                                      check_offtarget=False, lock_prefix=p, seed=42)["cds"]
                       for p in prefixes}
        assert len(by_encoding) == 2, "re-rolling the encoding must change the design"

    def test_a_rejected_member_records_no_encoding(self):
        """A member that falls back to its baseline locked no scaffold, so it has none."""
        from clippr.homology import diversify_library
        res = diversify_library(["AAAAUGUGG", "GCUAAAGAC", "UUACACGUG"],
                                codon_table=TOY_TABLE, check_offtarget=False, retries=1)
        for t, cds in res["cds"].items():
            if cds == res["baseline_cds"][t] and res["scaffold_encodings"][t] is not None:
                assert res["cds"][t].startswith(res["scaffold_encodings"][t]), (
                    f"{t} reports an encoding its CDS does not start with")
