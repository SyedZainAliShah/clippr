"""Controls for the intron advisory.

**The advisory never changes a design.** It reports where an intron *could* go and whether
putting one there would still satisfy this project's own synthesis rules. Every number it
reports is a prediction from someone else's experiments, in someone else's construct, in an
organism that is only sometimes ours -- so the thing that matters most here is that it says
so, and that it refuses the three hosts where an intron would destroy the gene rather than
merely fail to help.

Background for the thresholds: `docs/intron_feasibility.md`.
"""
from __future__ import annotations

import pytest

from ppropt import introns
from ppropt.synthesis_profile import STRICT_LEGACY

#: The canonical 9S Chlamydomonas CDS this repository produces for AAAAUGUGG at seed 42.
#: Pinned as a literal so these tests cost milliseconds rather than a design run.
from test_design import TOY_TABLE  # noqa: F401  (kept for parity with the other suites)


class TestTheIntronItself:
    def test_rbcs2i1_is_the_published_sequence(self):
        """145 nt, canonical GT..AG boundaries, from X04472.1:1307..1451."""
        assert len(introns.RBCS2I1) == 145
        assert introns.RBCS2I1.startswith("GT")
        assert introns.RBCS2I1.endswith("AG")

    def test_rbcs2i1_satisfies_this_projects_own_synthesis_rules(self):
        """The finding that makes any of this feasible: no exemption is needed.

        If this ever fails, the advisory is recommending a sequence our own policy would
        refuse, and it must stop recommending it.
        """
        problems = introns.sequence_problems(introns.RBCS2I1, STRICT_LEGACY)
        assert problems == [], problems

    def test_unspliced_it_would_destroy_the_reading_frame(self):
        """Why the three spliceosome-free hosts are refused rather than warned."""
        assert len(introns.RBCS2I1) % 3 == 1
        assert introns.frames_with_stops(introns.RBCS2I1) == [0, 1, 2]


class TestHostApplicability:
    @pytest.mark.parametrize("organism", sorted(introns.SPLICEOSOMAL_HOSTS))
    def test_a_nuclear_eukaryote_is_applicable(self, organism):
        assert introns.applicable(organism) is True

    @pytest.mark.parametrize("organism", ["e_coli", "c_reinhardtii_chloroplast",
                                          "n_tabacum_chloroplast"])
    def test_a_host_without_a_spliceosome_is_refused(self, organism):
        """Not a warning. An unspliced intron frameshifts and hits a stop in every frame."""
        assert introns.applicable(organism) is False
        why = introns.why_not(organism)
        assert why and "spliceosome" in why.lower()

    def test_evidence_is_scoped_to_the_host_it_was_measured_in(self):
        """Yeast and Arabidopsis have a spliceosome; rbcS2i1 was never shown to work there.

        Applicable and supported-by-evidence are different claims, and conflating them is
        how a Chlamydomonas result silently becomes a recommendation for Arabidopsis.
        """
        assert introns.has_evidence("c_reinhardtii_nuclear") is True
        assert introns.has_evidence("s_cerevisiae") is False
        assert introns.has_evidence("a_thaliana_nuclear") is False


class TestInsertionSites:
    def test_sites_are_between_two_guanines(self):
        """Baier's NG/GN rule. G|G spliced at 132% of control, C|C at 3%."""
        cds = "AAGGTTGGCC"
        assert introns.insertion_sites(cds) == [3, 7]

    def test_no_site_is_reported_at_a_sequence_boundary(self):
        """An intron before the first base or after the last is not an intron.

        Position 0 and position len(cds) are excluded. Everything between them is interior,
        including position 1 of "GG", which sits between the two bases.
        """
        assert introns.insertion_sites("") == []
        assert introns.insertion_sites("G") == []
        for cds in ("GG", "GGGG", "ATGGGCGGC"):
            sites = introns.insertion_sites(cds)
            assert all(0 < p < len(cds) for p in sites), sites

    def test_a_gc_rich_coding_sequence_offers_many_sites(self):
        """The reason this is feasible at all: Chlamydomonas codon choice makes G|G common."""
        cds = "ATG" + "GGC" * 100
        sites = introns.insertion_sites(cds)
        assert len(sites) > 100


class TestTheAdvisory:
    @pytest.fixture(scope="class")
    @classmethod
    def cds(cls):
        """A moderate-GC coding sequence, so junctions stay inside the band.

        Deliberately not a real PPROPT sequence: a real Chlamydomonas CDS runs at 0.606 GC
        and the intron at 0.579, and almost every junction between them breaches the ceiling.
        That case is pinned separately in TestTheRealChlamydomonasCase, which is the finding
        rather than the happy path.
        """
        return "ATG" + "ATGGAACTTAAGGGTACAATTGGAAATCTT" * 20

    def test_it_never_returns_a_modified_sequence(self, cds):
        """An advisory that edits the design is not an advisory."""
        out = introns.advise(cds, "c_reinhardtii_nuclear")
        assert "cds" not in out
        for value in out.values():
            assert value != cds

    def test_every_suggested_position_is_a_real_site(self, cds):
        out = introns.advise(cds, "c_reinhardtii_nuclear")
        sites = set(introns.insertion_sites(cds))
        assert out["suggested"]
        assert all(p in sites for p in out["suggested"])

    def test_suggested_positions_respect_the_exon_length_rule(self, cds):
        """Exons of 500 bp or shorter, against a native average of 240 nt."""
        out = introns.advise(cds, "c_reinhardtii_nuclear")
        assert all(length <= introns.MAX_EXON_NT for length in out["exon_lengths"])

    def test_a_refused_host_gets_no_suggestions_and_a_reason(self, cds):
        out = introns.advise(cds, "e_coli")
        assert out["applicable"] is False
        assert out["suggested"] == []
        assert "spliceosome" in out["why_not"].lower()

    def test_every_advisory_carries_the_disclaimer(self, cds):
        """Including the ones that succeed. The design itself has no intron in it."""
        for organism in ("c_reinhardtii_nuclear", "e_coli", "s_cerevisiae"):
            out = introns.advise(cds, organism)
            assert out["disclaimer"]
            assert "does not take introns into consideration" in out["disclaimer"]

    def test_a_site_that_would_break_the_gc_band_is_not_suggested(self):
        """The straddling-window problem, again.

        A 50-nt window spanning the exon/intron boundary sees GC-rich coding on one side.
        That is the same failure that refused every de novo design before the ordered
        fragment was constrained, and it must not be reintroduced by an advisory that
        recommends a position our own policy would then reject.
        """
        gc_rich = "ATG" + "GCGCCGCCGGCGCCGCCGGCGCCGCCGGCG" * 8
        out = introns.advise(gc_rich, "c_reinhardtii_nuclear")
        for position in out["suggested"]:
            assert introns.junction_problems(
                gc_rich, position, introns.RBCS2I1, STRICT_LEGACY) == []

    def test_a_breach_the_sequence_already_had_is_not_blamed_on_the_intron(self):
        """The distinction junction_problems exists to draw.

        This fixture sits at 0.33 GC, below the 0.35 floor, so windows near any junction
        breach whatever you insert. Blaming that on the intron would report "nowhere to put
        an intron" when the real finding is "this sequence is already out of band".
        """
        already_bad = "ATG" + "ATGGAACTTAAGGGTACAATTGGAAATCTT" * 20
        assert introns.sequence_problems(already_bad, STRICT_LEGACY), (
            "fixture must already breach the band for this test to mean anything")
        out = introns.advise(already_bad, "c_reinhardtii_nuclear")
        assert out["suggested"], "a pre-existing breach vetoed every candidate position"

    def test_the_evidence_field_names_the_host_and_the_effect(self, cds):
        out = introns.advise(cds, "c_reinhardtii_nuclear")
        assert "Chlamydomonas" in out["evidence"]
        assert out["intron"]["length_nt"] == 145


class TestTheRealChlamydomonasCase:
    """The finding this advisory exists to surface, pinned against a real design.

    The coding sequence satisfies the GC band. The intron satisfies the GC band. Almost
    every junction between them does not, because a 50-nt window straddling the seam sees
    GC-rich coding on one side and the intron's GC-rich 5' end on the other.

    This is the same straddling-window failure that refused every de novo design until the
    constraint was moved onto the ordered fragment, and it has the same remedy: choose the
    codons with the intron in place rather than after it. Until that is built, the honest
    answer is that no placement satisfies both published rules at once -- which is what the
    advisory now says, instead of returning an empty list and leaving the reader to guess.
    """

    #: The CDS this repository produces for AAAAUGUGG, c_reinhardtii_nuclear, seed 42.
    #: Pinned as a literal so this costs milliseconds rather than a design run.
    CDS = (
        "ATGCAGGGCGGCAACAGCGAGGAGCCCCGCAAAAGCTTCGATGAACGTCCTGAGCGCGGGGTAGTGAGCT"
        "GGACCGCCATGATCAGCGGTTATGCCCAGAACGGCCGGATCGACGAGGCAAGGGAGCTGTTCGACAAGAT"
        "GCCCGAGAGGAACGTGGTGAGCTGGACCGCGATGATCAGCGGTTACGCACAGAACGGCCGCATTGACGAG"
        "GCTCGTGAGCTATTCGACAAGATGCCCGAGCGCAACGTGGTGAGCTGGACGGCTATGATTAGCGGCTACG"
        "CACAGAATGGACGCATCGACGAGGCCCGCGAGCTGTTCGATAAGATGCCCGAGCGGAACGTGGTGAGCTG"
        "GACTGCCATGATCAGCGGCTATGCTCAGAACGGCCGCATAGACGAGGCCCGCGAACTGTTCGACAAGATG"
        "CCGGAGCGTAACGTGGTGAGCTGGAACGCAATG"
    )

    def test_the_coding_sequence_itself_is_in_band(self):
        """Otherwise the next test would be measuring the wrong thing."""
        assert introns.sequence_problems(self.CDS, STRICT_LEGACY) == []

    def test_embedding_the_intron_in_an_exon_oligo_almost_never_works(self):
        """Route A. The straddling window is synthesised as one piece, and it breaches.

        If this ever starts passing broadly, the finding behind
        `docs/intron_feasibility.md` has changed and the document needs revisiting.
        """
        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        assert out["sites_available"] > 30, "G|G sites should be common in GC-rich codons"
        a = out["route_a_embedded"]
        assert a["sites_clean"] < out["sites_available"] / 10, (
            f"route A was expected to be nearly always refused; "
            f"{a['sites_clean']} of {out['sites_available']} passed")

    def test_ordering_the_intron_as_its_own_part_does_work(self):
        """Route B, and the reason the answer is yes rather than no.

        The straddling window never exists at synthesis time, because the two sides are
        ordered as separate oligos. What must pass instead is the intron's own wrapped
        oligo, and it does -- with BbsI.
        """
        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        b = out["route_b_own_part"]
        assert b["enzyme"] == "BbsI"
        assert b["sites_clean"] > 10 * out["route_a_embedded"]["sites_clean"], (
            "route B should open an order of magnitude more positions than route A")
        assert out["suggested"], "a viable placement should be found"

    def test_it_explains_itself_rather_than_returning_an_empty_list(self):
        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        if not out["suggested"]:
            assert out.get("blocked_by"), (
                "no suggestions and no reason is the one unacceptable answer")
            assert "codons" in out["blocked_by"]


class TestTheEnzymeChoiceIsLoadBearing:
    """One base of recognition site decides whether an intron can be ordered at all.

    rbcS2i1's own worst 50-nt window is 0.640, against a 0.65 ceiling. Wrapping it shifts
    the window frame until some frame includes the enzyme site: BsaI's GGTCTC carries four
    G or C, BbsI's GAAGAC carries three. That is 0.660 against 0.640.

    An intron is not coding, so unlike a fragment it cannot be recoded to compensate, and
    the arbitrary pads cannot absorb it either. The enzyme is the only lever, which is why
    PART_ENZYME is a declared constant rather than an incidental default.

    It is also what the published systems do: MoInClo flanks its intron modules with BbsI,
    and GRASP assembles level 0 with BpiI, the same recognition site.
    """

    def test_the_intron_alone_sits_just_under_the_ceiling(self):
        from ppropt.synthesis_profile import STRICT_LEGACY as P
        worst = max(
            (introns.RBCS2I1[i:i + P.window].count("G")
             + introns.RBCS2I1[i:i + P.window].count("C")) / P.window
            for i in range(len(introns.RBCS2I1) - P.window + 1))
        assert worst <= P.local_gc[1], "the intron must pass on its own"
        assert P.local_gc[1] - worst < 0.02, (
            "the margin is supposed to be tiny; if it grew, the enzyme choice may no "
            "longer be load-bearing and this reasoning needs rechecking")

    def test_bsai_refuses_the_intron_and_bbsi_accepts_it(self):
        assert introns.part_problems(introns.RBCS2I1, "BsaI"), (
            "BsaI was expected to push the intron over the ceiling")
        assert introns.part_problems(introns.RBCS2I1, "BbsI") == [], (
            "BbsI was expected to keep it inside; if this fails the viable route is gone")

    def test_the_declared_part_enzyme_is_the_one_that_works(self):
        assert introns.part_problems(introns.RBCS2I1, introns.PART_ENZYME) == []


class TestPlacementIsChosenForAssembly:
    """Where the introns go decides whether the thing assembles at all.

    The overhang at an intron junction is the four bases before it. Every candidate sits
    between two Gs, so the obvious walk -- take the first clean site past each spacing
    target -- keeps picking the same four-mer. Identical overhangs in one Golden Gate
    reaction let fragments swap places.

    Measured on the three architectures with fixed codons: the naive walk gives 0.18 on 9S,
    0.03 on 14S and 0.003 on 19S, against a baseline near 0.83. Choosing the set for
    fidelity recovers 0.80 on 9S, 0.68 on 14S and 0.48 on 19S. So selection is not a
    refinement here; without it the advice is actively harmful.
    """

    CDS = TestTheRealChlamydomonasCase.CDS

    def test_no_two_junctions_share_an_overhang(self):
        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        ohs = [self.CDS[p - 4:p] for p in out["suggested"]]
        assert len(set(ohs)) == len(ohs), f"repeated overhang in {ohs}: fragments can swap"

    def test_the_chosen_set_beats_the_first_fit_walk(self):
        """The whole point of _select. If these converge, selection has stopped working."""
        from ppropt.assembly import reverse_complement
        from ppropt.assembly_spec import DESTINATION_OVERHANGS
        from ppropt.overhangs import set_fidelity

        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        clean = out["route_b_own_part"]["clean_positions"]
        assert out["suggested"], "expected a placement to compare against"

        naive, target = [], introns.NATIVE_EXON_NT
        while target < len(self.CDS):
            nxt = next((p for p in clean if p >= target
                        and p > (naive[-1] if naive else 0)
                        and p - (naive[-1] if naive else 0) <= introns.MAX_EXON_NT), None)
            if nxt is None:
                break
            naive.append(nxt)
            target = nxt + introns.NATIVE_EXON_NT

        dest = DESTINATION_OVERHANGS["level0"]
        fixed = [dest[0], reverse_complement(dest[1])]
        # Both sides must be scored on the enzyme the reaction uses. Pinning the naive walk
        # to a literal matrix made this comparison meaningless the moment selection moved
        # to BbsI-HF: it reported a BsaI-HFv2 walk beating a BbsI-HF selection.
        naive_score = set_fidelity(fixed + [self.CDS[p - 4:p] for p in naive],
                                   introns.PART_MATRIX)
        assert out["predicted_fidelity"] > naive_score, (
            f"selection {out['predicted_fidelity']} did not beat first-fit {naive_score:.4f}")

    def test_a_poor_placement_is_declared_rather_than_handed_over(self):
        """A number this far below baseline must not arrive without a sentence."""
        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        if out["predicted_fidelity"] is not None and out["predicted_fidelity"] < 0.75:
            assert out.get("fidelity_warning"), (
                "a placement under 0.75 was offered with no warning attached")


class TestTheIntronCountIsSearched:
    """How many introns to place is a result, not an input derived from spacing.

    Each extra intron adds a junction to the same reaction and fidelity falls steeply. The
    first version fixed the count at `round(len(cds) / 240) - 1`, which on 14S meant five
    introns scoring 0.517 where three score 0.687, and on 19S seven scoring 0.307. It was
    picking the worst end of the curve and reporting it as the recommendation.
    """

    CDS = TestTheRealChlamydomonasCase.CDS

    def test_a_curve_is_reported_not_a_single_count(self):
        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        curve = out["fidelity_by_intron_count"]
        assert curve, "the trade-off across counts should be reported"
        assert len(out["suggested"]) in curve

    def test_the_chosen_count_is_the_most_that_clears_the_floor(self):
        """More introns is more expression in the literature, so take the most affordable."""
        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        curve = out["fidelity_by_intron_count"]
        affordable = [k for k, f in curve.items() if f >= introns.FIDELITY_FLOOR]
        if affordable:
            assert len(out["suggested"]) == max(affordable), (
                f"chose {len(out['suggested'])} introns from curve {curve} with floor "
                f"{introns.FIDELITY_FLOOR}")

    def test_fidelity_does_not_improve_with_more_introns(self):
        """The shape the search exists to exploit. If it inverted, the reasoning is wrong."""
        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        curve = out["fidelity_by_intron_count"]
        counts = sorted(curve)
        assert curve[counts[0]] >= curve[counts[-1]], (
            f"more introns scored better than fewer, which contradicts the measurement: "
            f"{curve}")


class TestTheStrongerIntronIsRefused:
    """Why the advisory does not offer LHCBM1i2, which the literature prefers.

    Baier 2020 screened 33 native introns and LHCBM1i2 was the strongest, about 8-fold
    against rbcS2i1's 5.5. It cannot be used here: a 50-nt window at offset 14 reaches 0.680
    against this project's 0.65 ceiling, and the breach is *interior* to the intron rather
    than at a wrapper seam. No wrapping enzyme or pad moves an interior window, and an intron
    is not coding so it cannot be recoded. The only lever would be the band itself.
    """

    def test_it_is_the_published_sequence(self):
        assert len(introns.LHCBM1I2) == 253
        assert introns.LHCBM1I2.startswith("GT") and introns.LHCBM1I2.endswith("AG")

    def test_it_breaches_the_band_on_its_own(self):
        problems = introns.sequence_problems(introns.LHCBM1I2, STRICT_LEGACY)
        assert problems, "if this now passes, the stronger intron has become available"
        assert "0.68" in problems[0]

    def test_no_wrapping_enzyme_rescues_an_interior_breach(self):
        """The contrast with rbcS2i1, whose breach was at the seam and BbsI fixed."""
        for enzyme in ("BsaI", "BbsI", "SapI"):
            assert introns.part_problems(introns.LHCBM1I2, enzyme), (
                f"{enzyme} appeared to rescue an interior window, which should be impossible")
        assert introns.part_problems(introns.RBCS2I1, "BbsI") == [], (
            "rbcS2i1 with BbsI is the case that does work; it anchors the contrast")


class TestAnUnknownHostIsNotGivenABiologicalVerdict:
    """A mistyped organism must not collect a confident claim about its spliceosome."""

    def test_an_unknown_host_is_reported_as_unknown(self):
        why = introns.why_not("c_reinhardti_nuclear")        # one letter short
        assert why and "not one of the hosts" in why
        assert "has no nuclear spliceosome" not in why

    def test_a_known_host_without_a_spliceosome_still_says_so(self):
        why = introns.why_not("e_coli")
        assert why and "spliceosome" in why.lower()
        assert "not one of the hosts" not in why


class TestTheFloorIsAKnifeEdgeAndSaysSo:
    """Four of five seeds land at 0.7507 against a 0.75 floor. That is not a verdict.

    A threshold clearing by seven ten-thousandths is deciding more than its precision
    justifies: a floor of 0.751 would have dropped an intron from four of those five
    placements. The result therefore carries the margin, and a placement scraping over the
    line is flagged rather than presented as a clean pass.
    """

    CDS = TestTheRealChlamydomonasCase.CDS

    def test_the_margin_to_the_floor_is_reported(self):
        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        if out["predicted_fidelity"] is not None:
            assert out["margin_to_floor"] is not None
            assert abs(out["margin_to_floor"]
                       - (out["predicted_fidelity"] - introns.FIDELITY_FLOOR)) < 1e-6

    def test_a_placement_scraping_the_line_is_flagged(self):
        """Synthesised rather than hoped for: drive the floor up to just under the result."""
        out = introns.advise(self.CDS, "c_reinhardtii_nuclear")
        if out["predicted_fidelity"] is None:
            pytest.skip("no placement to probe")
        tight = round(out["predicted_fidelity"] - 0.005, 4)
        original = introns.FIDELITY_FLOOR
        try:
            introns.FIDELITY_FLOOR = tight
            probed = introns.advise(self.CDS, "c_reinhardtii_nuclear")
            assert probed["margin_to_floor"] < 0.01
            assert probed.get("fidelity_warning"), (
                "a placement within 0.01 of the floor was offered with no caveat")
        finally:
            introns.FIDELITY_FLOOR = original
