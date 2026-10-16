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
