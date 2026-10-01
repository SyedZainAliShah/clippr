"""Window selection, held to the acceptance cases the review asked for.

The ones that exist because a real defect hid there:

  - `TestFeatureIdentity` -- a gene name is not a location. The figure behind this module was
    wrong because occurrences were keyed on the name, and the annotation carries two `psbA`.
  - `TestRepeatsInsideAnInterval` -- excluding an interval can hide a second occurrence in it.
  - `TestAgreesWithOfftarget` -- the identity-preserving scan must not disagree with the already
    tested one about *whether* a sequence occurs, only about what it can tell apart.
"""
from __future__ import annotations

import pytest

from clippr.offtarget import Transcript, find_in_transcripts
from clippr.targeting import (COORDINATE_CONVENTION, IntendedLocus, Occurrence, SearchSequence,
                              enumerate_windows, feature_identity, loci_for_gene,
                              occurrences_in_reference, occurrences_in_source, window_report,
                              screen_window, select_window)

#: Two features sharing a name and sequence, the way the chloroplast annotation carries `psbA`.
TWIN_A = Transcript(name="dup", kind="CDS", strand=1, start=100,
                    sequence="ACGTACGTAACCGGTTACGTACGTAA")
TWIN_B = Transcript(name="dup", kind="CDS", strand=-1, start=9000,
                    sequence="ACGTACGTAACCGGTTACGTACGTAA")
SOLO = Transcript(name="solo", kind="CDS", strand=1, start=500,
                  sequence="GGGCCCGGGCCCAAATTTGGGCCCAAA")

#: A 60 nt sequence whose every 9-, 14- and 19-mer occurs exactly once in it. Needed because a
#: repetitive sequence has NO acceptable window -- every candidate recurs inside the source --
#: which is correct behaviour and useless as a fixture for testing the selected case. Generated
#: with a fixed seed and verified distinct at all three lengths.
UNIQUE = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG"


def seq_of(s: str, origin: str = "designed") -> SearchSequence:
    return SearchSequence(identifier="probe", sequence=s, origin=origin)


class TestSearchSequence:
    def test_origin_is_required_and_validated(self):
        with pytest.raises(ValueError, match="origin must be one of"):
            SearchSequence(identifier="x", sequence="ACGTACGTA", origin="synthetic")

    def test_origin_is_never_inferred(self):
        """A designed source's intended site is absent from a wild-type reference and a native
        one's is present, so guessing would change what every occurrence means."""
        for origin in ("designed", "native", "unknown"):
            assert SearchSequence("x", "ACGTACGTA", origin).origin == origin

    def test_refuses_non_nucleotides(self):
        with pytest.raises(ValueError):
            SearchSequence("x", "ACGTXYZ", "designed")

    def test_hash_identifies_the_content_not_the_name(self):
        a = SearchSequence("one", "ACGTACGTA", "designed")
        b = SearchSequence("two", "ACGUACGUA", "designed")      # same sequence, RNA spelling
        assert a.content_sha256 == b.content_sha256

    def test_rna_is_what_downstream_consumes(self):
        assert seq_of("ACGTACGTA").rna == "ACGUACGUA"


class TestWindowEnumeration:
    def test_every_offset_appears_once(self):
        got = enumerate_windows(seq_of("A" * 20), 9)
        assert len(got) == 12
        assert [o for o, _ in got] == list(range(12))

    def test_a_length_with_no_architecture_is_refused(self):
        """A PPR needs one repeat per base, so a 10-mer is not buildable."""
        with pytest.raises(ValueError, match="window length must be one of"):
            enumerate_windows(seq_of("A" * 40), 10)

    def test_a_sequence_shorter_than_the_window_yields_nothing(self):
        assert enumerate_windows(seq_of("ACGTA"), 9) == ()


class TestFeatureIdentity:
    def test_two_features_sharing_a_name_get_different_identities(self):
        a, b = feature_identity(0, TWIN_A), feature_identity(1, TWIN_B)
        assert a != b
        assert "dup" in a and "dup" in b

    def test_occurrences_distinguish_the_two_copies(self):
        """The defect this module exists to prevent: a name-keyed count merges two places."""
        occ = occurrences_in_reference("ACGTACGTA", [TWIN_A, TWIN_B])
        assert len(occ) == 4                      # twice in each copy
        assert len({o.source for o in occ}) == 2  # but two distinct sources
        assert len({o.label for o in occ}) == 1   # and only one name, which is the trap

    def test_loci_for_gene_returns_one_locus_per_feature(self):
        loci = loci_for_gene("dup", [TWIN_A, TWIN_B, SOLO])
        assert len(loci) == 2
        assert len(loci_for_gene("solo", [TWIN_A, TWIN_B, SOLO])) == 1

    def test_declaring_one_copy_does_not_absolve_the_other(self):
        trs = [TWIN_A, TWIN_B]
        one = loci_for_gene("dup", trs)[:1]
        w = screen_window(0, "ACGTACGTA", seq_of("ACGTACGTA"), trs, one)
        assert any(o.source == feature_identity(1, TWIN_B) for o in w.other_occurrences), \
            "the undeclared copy must still count against the window"

    def test_declaring_both_copies_clears_them(self):
        trs = [TWIN_A, TWIN_B]
        w = screen_window(0, "ACGTACGTA", seq_of("ACGTACGTA"), trs, loci_for_gene("dup", trs))
        assert w.other_occurrences == ()


class TestRepeatsInsideAnInterval:
    def test_a_repeat_in_the_supplied_sequence_is_visible(self):
        """Codex's control: the same 9-mer at offsets 0 and 9."""
        src = seq_of("ACGTACGTAACGTACGTA")
        w = screen_window(0, "ACGTACGTA", src)
        assert len(w.occurrences) == 2
        assert w.self_occurrence.offset == 0
        assert [o.offset for o in w.other_occurrences] == [9]

    def test_a_locus_absorbing_two_occurrences_says_so(self):
        """Excluding a whole interval can hide a second copy. Naming it keeps it visible."""
        trs = [TWIN_A]
        w = screen_window(0, "ACGTACGTA", seq_of("ACGTACGTA"), trs, loci_for_gene("dup", trs))
        assert w.other_occurrences == ()
        assert w.multiply_intended, "a locus covering two occurrences must be reported"
        assert "2 occurrences" in w.multiply_intended[0]

    def test_a_bounded_locus_does_not_absorb_beyond_its_end(self):
        trs = [TWIN_A]
        narrow = [IntendedLocus(source=feature_identity(0, TWIN_A), start=0, end=5)]
        w = screen_window(0, "ACGTACGTA", seq_of("ACGTACGTA"), trs, narrow)
        assert len(w.intended) == 1
        assert len(w.other_occurrences) == 1


class TestAgreesWithOfftarget:
    def test_same_total_count_as_the_tested_scanner(self):
        trs = [TWIN_A, TWIN_B, SOLO]
        for probe in ("ACGTACGTA", "GGGCCCGGG", "TTTTTTTTT"):
            mine = occurrences_in_reference(probe, trs)
            theirs = find_in_transcripts(probe, trs)
            assert len(mine) == len(theirs), f"{probe}: disagreement about whether it occurs"

    def test_offsets_match_too(self):
        trs = [TWIN_A, TWIN_B]
        mine = sorted(o.offset for o in occurrences_in_reference("ACGTACGTA", trs))
        theirs = sorted(h.position for h in find_in_transcripts("ACGTACGTA", trs))
        assert mine == theirs


class TestSelection:
    def test_picks_a_window_with_no_other_occurrence(self):
        sel = select_window(seq_of(UNIQUE), lengths=(9,),
                            transcripts=[SOLO], intended_loci=loci_for_gene("solo", [SOLO]))
        assert sel.found
        assert sel.chosen.other_occurrences == ()
        assert sel.chosen.length == 9

    def test_no_acceptable_window_is_an_answer_not_an_error(self):
        trs = [TWIN_A, TWIN_B]
        sel = select_window(seq_of(TWIN_A.sequence, "native"), lengths=(9,), transcripts=trs)
        assert not sel.found
        assert sel.chosen is None
        assert "or fewer other exact occurrences" in sel.reason
        assert sel.rejected, "the rejected candidates and their reasons must survive"

    def test_rejection_carries_its_reason(self):
        trs = [TWIN_A, TWIN_B]
        sel = select_window(seq_of(TWIN_A.sequence, "native"), lengths=(9,), transcripts=trs)
        _window, why = sel.rejected[0]
        assert "other exact occurrence" in why

    def test_a_sequence_too_short_says_so(self):
        sel = select_window(seq_of("ACGTA"), lengths=(9,))
        assert not sel.found
        assert "shorter than" in sel.reason

    def test_lengths_are_required(self):
        with pytest.raises(ValueError, match="lengths is required"):
            select_window(seq_of("A" * 40), lengths=())

    def test_several_lengths_are_allowed_and_recorded(self):
        sel = select_window(seq_of(UNIQUE), lengths=(9, 14))
        assert sel.lengths == (9, 14)

    def test_shorter_wins_a_tie_on_evidence_and_is_not_called_biological(self):
        """Fewer bases means fewer modules to order. It is a cost argument, not a binding one."""
        sel = select_window(seq_of(UNIQUE), lengths=(14, 9))
        assert sel.found and sel.chosen.length == 9

    def test_ties_are_counted_so_an_arbitrary_choice_is_visible(self):
        sel = select_window(seq_of(UNIQUE), lengths=(9,))
        assert sel.tie_count >= 1
        if sel.tie_count > 1:
            assert any("arbitrary" in n for n in sel.notes)

    def test_selection_is_deterministic(self):
        src = seq_of(UNIQUE)
        a = select_window(src, lengths=(9, 14), transcripts=[SOLO])
        b = select_window(src, lengths=(9, 14), transcripts=[SOLO])
        assert (a.chosen is None) == (b.chosen is None)
        if a.chosen:
            assert (a.chosen.sequence, a.chosen.offset) == (b.chosen.sequence, b.chosen.offset)


class TestProvenanceAndHonesty:
    def test_an_unscreened_reference_is_never_named(self):
        """A result must not carry a reference it did not consult."""
        sel = select_window(seq_of("A" * 40), lengths=(9,), transcripts=None)
        assert "none" in sel.reference
        assert any("no reference annotation was screened" in n for n in sel.notes)

    def test_absent_intent_is_recorded_as_undeclared_not_as_absent(self):
        sel = select_window(seq_of("A" * 40), lengths=(9,))
        assert sel.intended_declared is False
        assert any("not the same as none existing" in n for n in sel.notes)

    def test_the_window_carries_the_source_hash_and_convention(self):
        src = seq_of(UNIQUE)
        sel = select_window(src, lengths=(9,))
        assert sel.chosen.source_sha256 == src.content_sha256
        assert sel.chosen.coordinate_convention == COORDINATE_CONVENTION

    def test_raw_occurrences_are_retained_so_a_caller_can_recount(self):
        trs = [TWIN_A, TWIN_B]
        w = screen_window(0, "ACGTACGTA", seq_of("ACGTACGTA"), trs, loci_for_gene("dup", trs))
        assert len(w.occurrences) == 5           # 1 in the source, 2 in each twin
        assert w.other_occurrences == ()         # classified, not discarded

    def test_report_refuses_to_claim_binding(self):
        sel = select_window(seq_of("A" * 40), lengths=(9,))
        text = window_report(sel)
        assert "not demonstrated binding" in text
        assert "predicts affinity" in text

    def test_report_states_the_no_window_case_plainly(self):
        text = window_report(select_window(seq_of("ACGTA"), lengths=(9,)))
        assert "NO WINDOW SELECTED" in text

    def test_no_genomic_coordinate_is_offered(self):
        """`start + offset` is wrong across a join and on the minus strand, so it is absent."""
        occ = occurrences_in_source("ACGTACGTA", seq_of("ACGTACGTA"))
        assert not any(hasattr(o, attr) for o in occ
                       for attr in ("genomic_start", "genomic_position", "absolute_offset"))


class TestFeedsTheExistingRoutes:
    def test_the_chosen_sequence_is_a_plain_rna_target(self):
        """The whole point of the shape: existing entry points take it unchanged."""
        from clippr.ppr import describe

        sel = select_window(seq_of(UNIQUE), lengths=(9,))
        assert sel.found
        assert set(sel.chosen.sequence) <= set("ACGU")
        assert describe(sel.chosen.sequence)["architecture"] == "9S"

    def test_the_chosen_sequence_reaches_the_kit_route(self):
        from clippr.parts import select as select_parts

        sel = select_window(seq_of(UNIQUE), lengths=(9,))
        plan = select_parts(sel.chosen.sequence)
        assert plan.target == sel.chosen.sequence

    def test_an_intended_occurrence_declaration_does_not_change_the_sequence(self):
        """Declaring intent changes the evidence, never the DNA that gets ordered."""
        src = seq_of(UNIQUE)
        bare = select_window(src, lengths=(9,), transcripts=[SOLO],
                             max_other_occurrences=99)
        declared = select_window(src, lengths=(9,), transcripts=[SOLO],
                                 intended_loci=loci_for_gene("solo", [SOLO]),
                                 max_other_occurrences=99)
        assert bare.chosen.sequence == declared.chosen.sequence


class TestOccurrenceClassification:
    def test_self_occurrence_is_the_window_itself(self):
        w = screen_window(3, "ACGTACGTA", seq_of("GGGACGTACGTA"))
        assert w.self_occurrence is not None
        assert w.self_occurrence.offset == 3
        assert w.self_occurrence not in w.other_occurrences

    def test_an_occurrence_outside_every_locus_stays_other(self):
        loci = [IntendedLocus(source="feature:0:nothing")]
        w = screen_window(0, "ACGTACGTA", seq_of("ACGTACGTAACGTACGTA"), None, loci)
        assert len(w.other_occurrences) == 1

    def test_occurrence_str_names_what_and_where(self):
        o = Occurrence(source="feature:7:psbA", label="psbA", kind="CDS", offset=42)
        assert "psbA" in str(o) and "42" in str(o)


@pytest.fixture(scope="module")
def host():
    from clippr.offtarget import load_transcripts
    try:
        return load_transcripts()
    except Exception as exc:                            # pragma: no cover - needs the cache
        pytest.skip(f"chloroplast annotation unavailable: {exc}")


class TestAgainstTheRealAnnotation:
    """The figures `PROPOSAL_TO_CODEX.md` published, pinned as behaviour.

    They were measured with a throwaway script. A number nothing recomputes goes stale -- this
    project shipped a test badge that was wrong by 644 for exactly that reason -- so the numbers
    that justified this module are regressions now.
    """

    def _src(self, host, gene):
        feat = next(t for t in host if t.name == gene)
        return SearchSequence(f"{gene}_200", feat.sequence[:200], origin="native")

    def test_psba_is_duplicated_in_the_annotation(self, host):
        assert len(loci_for_gene("psbA", host)) == 2

    def test_declaring_both_psba_copies_leaves_128_of_192(self, host):
        sel = select_window(self._src(host, "psbA"), lengths=(9,), transcripts=host,
                            intended_loci=loci_for_gene("psbA", host))
        assert len(sel.rejected) == 64
        assert sel.found

    def test_declaring_one_psba_copy_leaves_nothing(self, host):
        """The inverted repeat means every window has a twin. 0 of 192, and that is the answer."""
        sel = select_window(self._src(host, "psbA"), lengths=(9,), transcripts=host,
                            intended_loci=loci_for_gene("psbA", host)[:1])
        assert not sel.found
        assert len(sel.rejected) == 192

    def test_rbcl_leaves_99_of_192(self, host):
        """99, not the 98 the proposal published. The difference is one window, below."""
        sel = select_window(self._src(host, "rbcL"), lengths=(9,), transcripts=host,
                            intended_loci=loci_for_gene("rbcL", host))
        assert len(sel.rejected) == 93

    def test_the_one_window_that_explains_99_versus_98(self, host):
        """`UUCCACCUG` occurs twice inside rbcL and nowhere else.

        The proposal counted the second copy as an off-target, because its predicate intended
        only the exact window position. Declaring the whole feature intended counts it as part
        of the gene being targeted. Both are defensible; the point is that the divergence is
        reported rather than silently decided, so the experimenter picks.
        """
        src = self._src(host, "rbcL")
        window = screen_window(142, "UUCCACCUG", src, host, loci_for_gene("rbcL", host))
        assert window.other_occurrences == ()
        assert window.multiply_intended, "the second copy inside the gene must be surfaced"
        assert "2 occurrences" in window.multiply_intended[0]

    def test_a_narrow_locus_recovers_the_proposals_stricter_count(self, host):
        """Bounding the locus to the searched interval makes the second copy count again."""
        src = self._src(host, "rbcL")
        feature = f"feature:{[t.name for t in host].index('rbcL')}:rbcL"
        narrow = [IntendedLocus(source=feature, start=0, end=192)]
        sel = select_window(src, lengths=(9,), transcripts=host, intended_loci=narrow)
        assert len(sel.rejected) == 94, "one more rejection than the whole-feature declaration"
