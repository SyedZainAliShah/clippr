"""Level-1 reactions: the geometry, and the guard against inventing one.

This package reported a level-1 reaction unscorable until 2026-09-16, on the stated grounds
that "its overhangs are set by the block plasmid and its destination, which this package does
not have". Both clauses were wrong. The block plasmids were in the provenance directory, and
they are not needed: the released block fragment is the joined module inserts, so a level-1
reaction's ends are the compiled product's own first and last interfaces.

Expected values come from the module table and from the deposited block plasmids, never from a
run of `compile_target`. `validation/experiments/level1_geometry.py` is the full derivation and
its two falsifiers; the tests against the primary records here are the same falsifiers, so a
change to the release rules that breaks the geometry fails a test rather than a document.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from clippr import assembly_spec as spec
from clippr import inventories as inv

ROOT = Path(__file__).resolve().parents[1]


#: Primary GenBank records, read for verification and never redistributed (`NOTICE.md`).
#: Resolved rather than hardcoded: an absolute developer path pins a repository to one
#: machine, and these tests skip cleanly where the records are absent.
def _primary_records() -> Path:
    override = os.environ.get("CLIPPR_GRASP_PRIMARY")
    if override:
        return Path(override)
    return ROOT.parent / "work" / "provenance" / "grasp_primary"

TABLE_S1 = ROOT / "data" / "grasp_supp" / "Table S1.xlsx"
PRIMARY = _primary_records()
BLOCKS = PRIMARY / "9SDYW_level0.gb"
PRODUCTS = PRIMARY / "9SrpoaDYW_variants_level1.gb"

#: Worked from the module table: 1A enters, each E module hands off, 2E exits.
EXPECTED = {
    "AAAAUGUGG": ["AATG", "CTTC", "TTCG"],
    "UUACACGUGCGUAC": ["AATG", "CTTC", "GTGA", "TTCG"],
    "CUAUCACAUCACAUAAGCG": ["AATG", "CTTC", "GTGA", "CACG", "TTCG"],
}


@pytest.fixture(scope="module")
def table():
    import json

    from clippr.codons import complete_table
    return complete_table(json.loads(
        (ROOT / "data" / "codon_tables" / "kazusa_3055.json").read_text("utf-8")), 1)


@pytest.fixture(scope="module")
def deposited():
    if not TABLE_S1.is_file():
        pytest.skip(f"supply {TABLE_S1.name} to run level-1 tests")
    return inv.load_deposited(TABLE_S1)


def level1_of(compiled) -> dict:
    found = [r for r in compiled["stages"]["reactions"] if r["stage"] == "level1"]
    assert len(found) == 1, found
    return found[0]


@pytest.mark.usefixtures("deposited")
class TestTheGeometryIsEstablishedButTheReactionIsNot:
    """Two different statuses, kept apart on purpose.

    The block *release geometry* is established: the deposited block plasmids and Table S1
    agree on every interface. The reaction's *participant list* is not, because the deposited
    BsaI reaction co-assembles parts this package does not compile.

    Ligation fidelity is a property of the whole competing overhang set, so a PPR-only figure
    is not this reaction's fidelity. That is not a formality: adding only the known omitted
    ends to the 19S set moves the prediction from 0.996046 to 0.742067 -- below level 0.
    """

    @pytest.mark.parametrize("target", sorted(EXPECTED))
    def test_the_block_subset_is_the_products_own_ends(self, deposited, target):
        reaction = level1_of(inv.compile_target(deposited, target))
        assert reaction["block_subset_overhangs"] == EXPECTED[target]
        assert reaction["ends_from"] == "product"
        assert reaction["geometry_established"]

    def test_the_reaction_is_not_reported_scorable(self, deposited):
        """The participants are unknown, so there is no complete-reaction fidelity."""
        for target in EXPECTED:
            reaction = level1_of(inv.compile_target(deposited, target))
            assert reaction["scorable"] is False
            assert reaction["participants_established"] is False
            assert reaction["reaction_overhangs"] is None

    def test_the_subset_is_never_presented_as_the_reactions_fidelity(self, deposited, table):
        """The regression for the defect this correction fixes.

        A subset score sitting in `fidelity` or `reaction_overhangs` would read as the
        reaction's own number while describing a tube that does not exist.
        """
        from clippr import joint_search as js

        targets = ["CUAUCACAUCACAUAAGCG"]
        classes = js.junction_classes(deposited, targets)
        got = js.evaluate(deposited, classes, {c.key: c.deposited for c in classes},
                          targets, table)
        level1 = [r for rows in got.per_reaction_fidelity.values() for r in rows
                  if r["stage"] == "level1"]
        assert level1
        for reaction in level1:
            assert reaction["fidelity"] is None
            assert reaction["block_subset"]["overhangs"]
            assert "diagnostic" in reaction["block_subset"]["is"]
        assert got.unscorable_reactions, "level 1 must be reported unscorable"

    def test_the_known_omitted_ends_change_the_number(self):
        """Why a subset score cannot stand in for the reaction, in one assertion.

        Values are the sensitivity calculation, not a reconstructed reaction: the true
        participant list is still unknown and may move it further.
        """
        from clippr.assembly_spec import LEVEL1_KNOWN_COASSEMBLED_ENDS
        from clippr.overhangs import set_fidelity

        subset = EXPECTED["CUAUCACAUCACAUAAGCG"]
        with_known = subset + list(LEVEL1_KNOWN_COASSEMBLED_ENDS)
        assert set_fidelity(subset, "BsaI-HFv2") == pytest.approx(0.996046, abs=1e-6)
        assert set_fidelity(with_known, "BsaI-HFv2") == pytest.approx(0.742067, abs=1e-6)
        # and the consequence that retires the "level 0 is the bottleneck" headline
        assert set_fidelity(with_known, "BsaI-HFv2") < 0.751683

    def test_the_entry_follows_the_chosen_fusion_variant(self, deposited):
        """`1A` exists as an `AATG` and an `AGGT` variant, so a constant here would be wrong."""
        for fusion in ("AATG", "AGGT"):
            reaction = level1_of(inv.compile_target(deposited, "AAAAUGUGG",
                                                    fusion_site=fusion))
            assert reaction["block_subset_overhangs"][0] == fusion

    def test_it_is_a_bsai_reaction(self, deposited):
        reaction = level1_of(inv.compile_target(deposited, "AAAAUGUGG"))
        assert reaction["enzyme"] == "BsaI"
        assert reaction["matrix"] == "BsaI-HFv2"

    def test_a_level_zero_reaction_is_scored_and_complete(self, deposited):
        """Level 0's participants *are* known, so it carries a real fidelity."""
        compiled = inv.compile_target(deposited, "AAAAUGUGG")
        level0 = [r for r in compiled["stages"]["reactions"] if r["stage"] == "level0"]
        assert level0
        for reaction in level0:
            assert reaction["scorable"] is True
            assert reaction["participants_established"] is True
            assert reaction["ends_from"] == "destination vector"
            assert reaction["reaction_overhangs"][0] == "CTCA"


class TestAgainstTheDepositedBlockPlasmids:
    """The falsifiers. If the plasmids disagree with the module table, these say so."""

    @pytest.fixture(scope="class")
    @classmethod
    def released(cls):
        if not BLOCKS.is_file():
            pytest.skip("primary block-plasmid records not available")
        import sys

        sys.path.insert(0, str(ROOT / "validation" / "experiments"))
        from level1_geometry import bsai_release, is_block
        from Bio import SeqIO

        out = {}
        for record in SeqIO.parse(BLOCKS, "genbank"):
            if is_block(record):
                out[record.name] = bsai_release(str(record.seq))
        return out

    def test_every_block_plasmid_releases_one_fragment(self, released):
        assert len(released) == 28

    def test_the_released_ends_are_the_module_interfaces(self, released):
        classes = {(five, three) for _, five, three in released.values()}
        assert classes == {("AGGT", "CTTC"), ("CTTC", "TTCG")}

    def test_the_joined_blocks_are_the_compiled_product_length(self, released):
        """499 + 406 - 4 = 901, and 901 is what this package compiles for 9S."""
        left = next(f for f, five, three in released.values() if three == "CTTC")
        right = next(f for f, five, three in released.values() if five == "CTTC")
        assert len(left) + len(right) - 4 == 901

    def test_the_joined_blocks_occur_in_a_deposited_level_one_product(self, released):
        """The falsifier that matters: the derivation must reproduce a real construct."""
        if not PRODUCTS.is_file():
            pytest.skip("primary level-1 product records not available")
        from Bio import SeqIO

        left = next(f for f, five, three in released.values() if three == "CTTC")
        right = next(f for f, five, three in released.values() if five == "CTTC")
        joined = left + right[4:]
        products = [str(r.seq).upper() for r in SeqIO.parse(PRODUCTS, "genbank")]
        assert any(joined in (p + p) for p in products)


class TestWhatIsStillNotCovered:
    def test_the_coassembled_parts_are_declared_unmodelled(self):
        """A level-1 figure is the PPR block junctions only, and the source must say so."""
        assert spec.LEVEL1_COASSEMBLED_PARTS_UNMODELLED is True

    def test_the_structured_status_agrees_with_the_prose(self):
        """The defect this file was rewritten for: a caveat in prose and `True` in a flag.

        `context_available=True` said the missing context was resolved while the note beside
        it said the opposite. A reader believes the field.
        """
        assert spec.LEVEL1_CONTEXT_AVAILABLE is False
        assert spec.LEVEL1_PARTICIPANTS_ESTABLISHED is False
        assert spec.ASSEMBLY_STAGES["level1"]["context_available"] is False
        assert spec.ASSEMBLY_STAGES["level1"]["geometry_established"] is True

    def test_the_two_level_ones_are_not_confused(self):
        """The acceptor vector's own ends and a block reaction's ends are different sets."""
        assert spec.DESTINATION_OVERHANGS["level1"] == ("GGAG", "CGCT")
        assert spec.ASSEMBLY_STAGES["level1"]["destination"] is None
        assert spec.ASSEMBLY_STAGES["level1"]["flanks_from_product"] is True
