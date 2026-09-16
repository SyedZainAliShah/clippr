"""Assembly substrates: the release rules, the wrapper, and the digest that proves it.

Expected values are derived from the deposited plasmids or worked by hand, never from a run of
the builder. `verify_rules` is the falsifier for the rules themselves and re-derives them from
the primary records with arithmetic the builder does not share.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from clippr import inventories as inv
from clippr import substrates as sub

ROOT = Path(__file__).resolve().parents[1]
TABLE_S1 = ROOT / "data" / "grasp_supp" / "Table S1.xlsx"


#: Primary GenBank records, read for verification and never redistributed (`NOTICE.md`).
#: Resolved rather than hardcoded: an absolute developer path pins a repository to one
#: machine, and these tests skip cleanly where the records are absent.
def _primary_records() -> Path:
    override = os.environ.get("CLIPPR_GRASP_PRIMARY")
    if override:
        return Path(override)
    return ROOT.parent / "work" / "provenance" / "grasp_primary"


PRIMARY = _primary_records() / "GRASP_-1.gb"

#: The level-0 overhang chain, from the primary cut geometry.
LEVEL0 = {"CTCA", "ACTC", "AAGA", "GCAC", "TGAA", "CGAG"}


@pytest.fixture(scope="module")
def deposited():
    if not TABLE_S1.is_file():
        pytest.skip(f"supply {TABLE_S1.name} to run substrate tests")
    return inv.load_deposited(TABLE_S1)


class TestDigest:
    def test_a_wrapped_fragment_comes_back_intact(self):
        fragment = "CTCA" + "ATGAAACCC" * 3 + "ACTC"
        built = sub.build(fragment, "B")
        assert built.released == fragment

    def test_the_exposed_ends_are_the_fragment_ends(self):
        fragment = "CTCA" + "ATGAAACCC" * 3 + "ACTC"
        built = sub.build(fragment, "B")
        assert built.five_overhang == "CTCA"
        assert built.three_overhang == "ACTC"

    def test_the_substrate_is_longer_than_what_it_releases(self):
        """It has to be: the wrapper is cut off."""
        fragment = "CTCA" + "ATGAAACCC" * 3 + "ACTC"
        built = sub.build(fragment, "B")
        assert len(built.sequence) > len(built.released)

    def test_a_second_site_is_refused(self):
        """An extra site means the enzyme cuts somewhere the design did not intend."""
        fragment = "CTCA" + "GAAGAC" + "ATGAAACCC" + "ACTC"
        with pytest.raises(sub.SubstrateError, match="more than twice"):
            sub.build(fragment, "B")

    def test_an_unknown_block_is_refused_rather_than_guessed(self):
        with pytest.raises(sub.SubstrateError, match="no release rule"):
            sub.build("ACGT" * 10, "Z9")


class TestReleaseRules:
    def test_an_a_module_gains_its_left_overhang(self):
        """Its insert starts with the fusion site; CTCA comes from the vector."""
        insert = "AATG" + "CCC" * 10 + "ACTC"
        assert sub.released_fragment(insert, "1A") == "CTCA" + insert

    def test_a_bcd_module_already_spans_overhang_to_overhang(self):
        insert = "ACTC" + "CCC" * 10 + "AAGA"
        assert sub.released_fragment(insert, "B") == insert

    def test_an_e_module_gains_its_right_overhang(self):
        insert = "TGAA" + "CCC" * 10 + "CTTC"
        assert sub.released_fragment(insert, "1E") == insert + "CGAG"


@pytest.mark.usefixtures("deposited")
class TestAgainstThePrimaryRecords:
    def test_the_rules_hold_for_every_deposited_module(self, deposited):
        """The falsifier. If the plasmids disagree, this says so."""
        if not PRIMARY.is_file():
            pytest.skip("primary GenBank records not available")
        got = sub.verify_rules(PRIMARY, deposited)
        assert got["checked"] == 42
        assert got["rules_hold"], got["disagreements"][:3]

    def test_every_module_builds_a_substrate(self, deposited):
        built = sub.substrates_for(deposited)
        assert len(built) == len(deposited)

    def test_every_substrate_exposes_level_zero_ends(self, deposited):
        for s in sub.substrates_for(deposited):
            assert s.five_overhang in LEVEL0
            assert s.three_overhang in LEVEL0

    def test_the_exposed_ends_form_the_level_zero_chain(self, deposited):
        """CTCA -> ACTC -> AAGA -> GCAC -> TGAA -> CGAG, and nothing else."""
        pairs = {(s.five_overhang, s.three_overhang)
                 for s in sub.substrates_for(deposited)}
        assert pairs == {("CTCA", "ACTC"), ("ACTC", "AAGA"), ("AAGA", "GCAC"),
                         ("GCAC", "TGAA"), ("TGAA", "CGAG")}

    def test_no_substrate_carries_a_stray_recognition_site(self, deposited):
        from Bio.Restriction import BbsI, BsaI
        from Bio.Seq import Seq

        for s in sub.substrates_for(deposited):
            assert len(BbsI.search(Seq(s.sequence), linear=True)) == 2
            assert len(BsaI.search(Seq(s.sequence), linear=True)) == 0

    def test_a_substrate_is_not_its_insert(self, deposited):
        """The defect this module exists to fix: ordering bare inserts."""
        for s in sub.substrates_for(deposited):
            assert s.sequence != deposited.modules[s.module_id].dna


@pytest.mark.usefixtures("deposited")
class TestSynthesisContractOnWhatIsOrdered:
    """The band applies to the ordered sequence, not to a sub-span of it.

    `pPR-1_19E_LN5N` passed at 0.640 as a bare insert and reached 0.660 as the released
    fragment, because an E module's fragment is its insert plus the GC-rich `CGAG` overhang.
    Checking the insert alone let that reach an order reporting success.
    """

    def test_a_substrate_reports_its_own_contract_breaches(self):
        # an insert engineered to sit just under the band, pushed over by the CGAG overhang
        insert = "TGAA" + "GCCGCC" * 8 + "CTTC"
        built = sub.build(insert, "1E")
        assert built.synthesis_problems
        assert any("GC" in p for p in built.synthesis_problems)

    def test_a_clean_substrate_reports_none(self):
        insert = "ACTC" + "ATGCAC" * 12 + "AAGA"
        assert sub.build(insert, "B").synthesis_problems == ()

    def test_the_intended_enzyme_sites_are_not_counted_against_it(self):
        """Both BbsI sites are the construction; counting them would fail every substrate."""
        insert = "ACTC" + "ATGCAC" * 12 + "AAGA"
        built = sub.build(insert, "B")
        assert "GAAGAC" in built.sequence
        assert not any("BbsI" in p for p in built.synthesis_problems)

    def test_a_long_homopolymer_in_the_fragment_is_caught(self):
        insert = "ACTC" + "ATGCAC" * 6 + "A" * 8 + "ATGCAC" * 6 + "AAGA"
        assert any("homopolymer" in p for p in sub.build(insert, "B").synthesis_problems)

    def test_the_shipped_inventory_has_no_contract_breaches(self, deposited):
        """Recoding now targets the released fragment, so this must come out empty."""
        from clippr import inventories as _inv

        path = Path(__file__).resolve().parents[1] / "work" / "phaseb" / "inventory_recoded.json"
        if not path.is_file():
            pytest.skip("no recoded inventory on disk")
        offenders = {s.module_id: s.synthesis_problems
                     for s in sub.substrates_for(_inv.load(path)) if s.synthesis_problems}
        assert offenders == {}
