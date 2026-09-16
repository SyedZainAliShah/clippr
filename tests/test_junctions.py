"""Junction candidate derivation and assembly-model validity."""
from __future__ import annotations

from pathlib import Path

import pytest

from clippr import inventories as inv
from clippr.junctions import (JunctionSite, permitted_overhangs, site_for, sites_in,
                              valid_assignment)

TABLE_S1 = Path(__file__).resolve().parents[1] / "data" / "grasp_supp" / "Table S1.xlsx"


@pytest.fixture(scope="module")
def deposited():
    if not TABLE_S1.is_file():
        pytest.skip(f"supply {TABLE_S1.name} to run junction tests")
    return inv.load_deposited(TABLE_S1)


class TestValidAssignment:
    def test_the_deposited_ends_are_acceptable(self):
        assert valid_assignment(["ACTC", "AAGA", "GCAC", "TGAA", "CTTC"])[0]

    def test_a_repeated_overhang_is_refused(self):
        ok, why = valid_assignment(["ACTC", "ACTC"])
        assert not ok and "twice" in why

    def test_a_palindrome_is_refused(self):
        """AATT anneals to itself, so the product's order is not determined."""
        ok, why = valid_assignment(["AATT"])
        assert not ok and "palindromic" in why

    def test_a_complementary_pair_is_refused(self):
        ok, why = valid_assignment(["AAGA", "TCTT"])
        assert not ok and "complementary" in why

    def test_a_clash_with_the_destination_is_refused(self):
        """The backbone's own ends are in the reaction too."""
        ok, why = valid_assignment(["CTCA"])
        assert not ok and "twice" in why

    def test_a_different_destination_changes_what_clashes(self):
        assert valid_assignment(["CTCA"], destination="level1")[0]
        assert not valid_assignment(["GGAG"], destination="level1")[0]


class TestPermittedOverhangs:
    def _site(self, amino_acids="EL", offset=2, deposited="ACTC"):
        return JunctionSite(role="t", start=0, codon_start=0, codon_end=3 * len(amino_acids),
                            offset=offset, amino_acids=amino_acids, deposited=deposited)

    def test_the_incumbent_is_always_among_its_own_options(self):
        """An incumbent missing from its own candidate set would be a bug, not a result."""
        options = permitted_overhangs(self._site())
        assert "ACTC" in options

    def test_every_option_preserves_the_protein(self):
        from Bio.Seq import Seq

        site = self._site("AR", offset=2, deposited="AAGA")
        for overhang, codons in permitted_overhangs(site).items():
            if not codons:
                continue
            assert str(Seq(codons).translate()) == "AR"
            assert codons[site.offset:site.offset + 4] == overhang

    def test_forbidden_patterns_are_excluded(self):
        site = self._site("AR", offset=0, deposited="GCAC")
        without = permitted_overhangs(site)
        with_filter = permitted_overhangs(site, forbidden=("GC",))
        assert len(with_filter) < len(without)
        assert not any("GC" in o for o in with_filter if o != site.deposited)

    def test_a_stop_codon_context_offers_only_the_incumbent(self):
        site = self._site("*L", offset=2, deposited="ACTC")
        assert permitted_overhangs(site) == {"ACTC": ""}


@pytest.mark.usefixtures("deposited")
class TestRealSites:
    def test_every_site_reports_its_codon_context(self, deposited):
        compiled = inv.compile_target(deposited, "AAAAUGUGG")
        for site in sites_in(compiled, deposited):
            assert len(site.amino_acids) in (2, 3)
            assert 0 <= site.offset < 3
            assert site.codon_end - site.codon_start == 3 * len(site.amino_acids)

    def test_the_overhang_sits_where_the_site_says_it_does(self, deposited):
        compiled = inv.compile_target(deposited, "AAAAUGUGG")
        product = compiled["product"]
        for site in sites_in(compiled, deposited):
            assert product[site.start:site.start + 4] == site.deposited
            span = product[site.codon_start:site.codon_end]
            assert span[site.offset:site.offset + 4] == site.deposited

    def test_each_role_has_exactly_one_codon_context(self, deposited):
        """A global assignment is only coherent if a role means one physical interface."""
        seen = {}
        for target in ("AAAAUGUGG", "UUACACGUGCGUAC", "CUAUCACAUCACAUAAGCG"):
            compiled = inv.compile_target(deposited, target)
            for site in sites_in(compiled, deposited):
                signature = (site.amino_acids, site.offset, site.deposited)
                seen.setdefault(site.role, set()).add(signature)
        assert all(len(v) == 1 for v in seen.values())

    def test_real_sites_offer_real_alternatives(self, deposited):
        """If every site permitted only its incumbent, the search would be pointless."""
        compiled = inv.compile_target(deposited, "AAAAUGUGG")
        counts = [len(permitted_overhangs(s)) for s in sites_in(compiled, deposited)]
        assert min(counts) > 1

    def test_site_for_is_consistent_with_a_hand_worked_case(self):
        # frame 1, so codons are [1,4) ATG, [4,7) AAA, [7,10) CCC.
        # A junction at 4 starts exactly on a codon boundary: offset 0, spanning AAA and CCC.
        product = "XATGAAACCC"
        site = site_for(product, frame=1, start=4, deposited=product[4:8], role="t")
        assert (site.codon_start, site.codon_end) == (4, 10)
        assert site.offset == 0
        assert site.amino_acids == "KP"

    def test_a_junction_off_the_codon_boundary_reports_its_offset(self):
        # A junction at 5 starts one base into the AAA codon.
        product = "XATGAAACCC"
        site = site_for(product, frame=1, start=5, deposited=product[5:9], role="t")
        assert (site.codon_start, site.codon_end) == (4, 10)
        assert site.offset == 1
        assert product[site.codon_start:site.codon_end][1:5] == site.deposited
