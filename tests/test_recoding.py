"""Fixed-interface recoding: what it preserves, what it refuses, and how it reports coverage."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from clippr import inventories as inv
from clippr.recoding import (GC_BAND, locked_interface_sites, module_constraint_problems,
                             recode_inventory)

TABLE_S1 = Path(__file__).resolve().parents[1] / "data" / "grasp_supp" / "Table S1.xlsx"
CODON_TABLE = Path(__file__).resolve().parents[1] / "data" / "codon_tables" / "kazusa_3055.json"


@pytest.fixture(scope="module")
def deposited():
    if not TABLE_S1.is_file():
        pytest.skip(f"supply {TABLE_S1.name} to run recoding tests")
    return inv.load_deposited(TABLE_S1)


@pytest.fixture(scope="module")
def table():
    return json.loads(CODON_TABLE.read_text(encoding="utf-8"))


class TestLockedSites:
    def test_only_the_interface_inside_the_coding_span_is_locked(self):
        """Locking all four bases would address positions past the optimised sequence."""
        dna = "AATG" + "CCC" * 10 + "ACTC"
        locked = locked_interface_sites(dna, frame=1)
        start = 1
        n = (len(dna) - start) // 3
        end = start + 3 * n
        for offset, text in locked.items():
            assert 0 <= offset
            assert offset + len(text) <= end - start

    def test_the_locked_text_matches_the_sequence_it_came_from(self):
        dna = "AATG" + "CCC" * 10 + "ACTC"
        for offset, text in locked_interface_sites(dna, frame=1).items():
            assert dna[1 + offset:1 + offset + len(text)] == text


class TestWholeModuleConstraints:
    """The optimiser constrains the span it is given; the ordered fragment is the module."""

    def test_a_clean_sequence_has_no_problems(self):
        assert module_constraint_problems("ATGC" * 30) == []

    def test_a_forbidden_site_is_found_on_either_strand(self):
        forward = module_constraint_problems("ATGC" * 10 + "GGTCTC" + "ATGC" * 10)
        reverse = module_constraint_problems("ATGC" * 10 + "GAGACC" + "ATGC" * 10)
        assert any("BsaI" in p for p in forward)
        assert any("BsaI" in p for p in reverse)

    def test_gc_outside_the_band_is_found(self):
        problems = module_constraint_problems("GCGC" * 30)
        assert any("GC" in p and str(GC_BAND[1]) in p for p in problems)

    def test_a_long_homopolymer_is_found(self):
        assert any("homopolymer" in p
                   for p in module_constraint_problems("ATGCATGC" + "A" * 6 + "ATGCATGC"))


@pytest.mark.usefixtures("deposited", "table")
class TestRecoding:
    def test_every_module_gets_a_realisation(self, deposited, table):
        report = recode_inventory(deposited, table, label="t", seeds=2, wall_seconds=600)
        realised = len(report.improved) + len(report.unchanged)
        assert realised == len(deposited), "coverage is over every module, not the moved ones"
        assert len(report.inventory) == len(deposited)

    def test_proteins_lengths_and_interfaces_are_preserved(self, deposited, table):
        from Bio.Seq import Seq

        report = recode_inventory(deposited, table, label="t", seeds=2, wall_seconds=600)
        for module_id, before in deposited.modules.items():
            after = report.inventory.modules[module_id]
            assert len(after.dna) == len(before.dna)
            assert after.dna[:4] == before.dna[:4]
            assert after.dna[-4:] == before.dna[-4:]
            s, e = before.coding_interval
            assert (str(Seq(after.dna[s:e]).translate())
                    == str(Seq(before.dna[s:e]).translate()))

    def test_no_delivered_module_violates_the_whole_module_contract(self, deposited, table):
        """Two modules first passed the optimiser's own check and failed this one."""
        report = recode_inventory(deposited, table, label="t", seeds=2, wall_seconds=600)
        offenders = {m: module_constraint_problems(r.dna)
                     for m, r in report.inventory.modules.items()
                     if module_constraint_problems(r.dna)}
        assert offenders == {}

    def test_a_zero_budget_attempts_nothing_and_says_so(self, deposited, table):
        report = recode_inventory(deposited, table, label="t", seeds=1, wall_seconds=0.0)
        assert report.not_attempted and report.budget_exhausted
        assert not report.improved
        # The inventory is still complete: unattempted modules keep their incumbents.
        assert len(report.inventory) == len(deposited)

    def test_the_result_is_a_new_version_of_the_same_identities(self, deposited, table):
        report = recode_inventory(deposited, table, label="t", seeds=2, wall_seconds=600)
        assert set(report.inventory.modules) == set(deposited.modules)
        assert report.inventory.version != deposited.version

    def test_the_recoded_inventory_still_compiles(self, deposited, table):
        report = recode_inventory(deposited, table, label="t", seeds=2, wall_seconds=600)
        got = inv.compile_target(report.inventory, "AAAAUGUGG")
        assert got["available"] and got["product_nt"] == 901
