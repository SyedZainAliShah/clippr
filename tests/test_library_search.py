"""The collection optimiser: what it guarantees, and what it honestly reports failing to do."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from clippr import inventories as inv
from clippr.library_search import MODES, collection_objectives, optimise_library
from clippr.recoding import module_constraint_problems

TABLE_S1 = Path(__file__).resolve().parents[1] / "data" / "grasp_supp" / "Table S1.xlsx"
CODON_TABLE = Path(__file__).resolve().parents[1] / "data" / "codon_tables" / "kazusa_3055.json"


@pytest.fixture(scope="module")
def deposited():
    if not TABLE_S1.is_file():
        pytest.skip(f"supply {TABLE_S1.name} to run library search tests")
    return inv.load_deposited(TABLE_S1)


@pytest.fixture(scope="module")
def table():
    return json.loads(CODON_TABLE.read_text(encoding="utf-8"))


class TestCollectionObjectives:
    def test_sharing_is_normalised_by_pair_count(self):
        """Otherwise the term grows with inventory size rather than with actual sharing."""
        from clippr.codons import complete_table
        t = complete_table({"K": {"AAA": 1.0, "AAG": 0.5}}, 1)
        two = collection_objectives({"a": "AAAAAAAAAAAA", "b": "AAAAAAAAAAAA"},
                                    {"a": 0, "b": 0}, t, k=4)
        four = collection_objectives({k: "AAAAAAAAAAAA" for k in "abcd"},
                                     {k: 0 for k in "abcd"}, t, k=4)
        assert two["sharing_per_pair"] == pytest.approx(four["sharing_per_pair"])

    def test_the_scalar_reports_the_k_it_used(self):
        from clippr.codons import complete_table
        t = complete_table({"K": {"AAA": 1.0}}, 1)
        got = collection_objectives({"a": "AAAAAA"}, {"a": 0}, t, k=7)
        assert got["k"] == 7


@pytest.mark.usefixtures("deposited", "table")
class TestSearchContract:
    def test_an_unknown_mode_is_refused(self, deposited, table):
        with pytest.raises(ValueError, match="not one of"):
            optimise_library(deposited, table, mode="telepathy")

    def test_a_zero_budget_returns_the_incumbent_and_says_so(self, deposited, table):
        got = optimise_library(deposited, table, mode="greedy", max_proposals=0,
                               wall_seconds=0.0)
        assert got.inventory.version == deposited.version
        assert not got.improved
        assert got.as_dict()["completion"] == "budget_exhausted"

    def test_budget_exhausted_is_not_a_claim_about_the_search_space(self, deposited, table):
        """'Ran out of budget' and 'no improvement exists' are different statements."""
        got = optimise_library(deposited, table, mode="greedy", max_proposals=2,
                               wall_seconds=600)
        assert got.budget_exhausted
        assert got.proposals <= 2

    def test_cache_free_evaluations_are_counted_separately_from_proposals(self, deposited,
                                                                         table):
        got = optimise_library(deposited, table, mode="greedy", max_proposals=6,
                               wall_seconds=600)
        assert got.evaluations <= got.proposals
        assert got.evaluations + got.infeasible <= got.proposals

    @pytest.mark.parametrize("mode", MODES)
    def test_the_search_never_introduces_a_violation(self, deposited, table, mode):
        """The guarantee is about what the search produces, not what it inherited.

        One deposited module (`pPR-1_D_LD5T`) contains `TTTTT` and so already fails CLIPPR's
        max-homopolymer-4 rule. A module the search leaves untouched keeps whatever the input
        had; what must never happen is the search *creating* a violation.
        """
        inherited = {m for m, r in deposited.modules.items()
                     if module_constraint_problems(r.dna)}
        got = optimise_library(deposited, table, mode=mode, seed=42, max_proposals=20,
                               wall_seconds=600)
        introduced = {m: module_constraint_problems(r.dna)
                      for m, r in got.inventory.modules.items()
                      if module_constraint_problems(r.dna) and m not in inherited}
        assert introduced == {}

    def test_the_deposited_inventory_is_reported_against_the_synthesis_profile(self,
                                                                              deposited):
        """A user ordering the deposited kit should learn this before, not after.

        CLIPPR's homopolymer limit of 4 is an engineering choice, not a published rule, so
        this is a statement about fit to our synthesis profile -- not a defect in the kit.
        """
        offenders = {m: module_constraint_problems(r.dna)
                     for m, r in deposited.modules.items()
                     if module_constraint_problems(r.dna)}
        assert set(offenders) == {"pPR-1_D_LD5T"}
        assert "homopolymer" in offenders["pPR-1_D_LD5T"][0]

    @pytest.mark.parametrize("mode", MODES)
    def test_proteins_and_interfaces_survive_every_mode(self, deposited, table, mode):
        from Bio.Seq import Seq

        got = optimise_library(deposited, table, mode=mode, seed=42, max_proposals=20,
                               wall_seconds=600)
        for module_id, before in deposited.modules.items():
            after = got.inventory.modules[module_id]
            assert len(after.dna) == len(before.dna)
            assert after.dna[:4] == before.dna[:4] and after.dna[-4:] == before.dna[-4:]
            s, e = before.coding_interval
            assert (str(Seq(after.dna[s:e]).translate())
                    == str(Seq(before.dna[s:e]).translate()))

    def test_the_same_seed_reproduces(self, deposited, table):
        a = optimise_library(deposited, table, mode="anneal", seed=7, max_proposals=20,
                             wall_seconds=600)
        b = optimise_library(deposited, table, mode="anneal", seed=7, max_proposals=20,
                             wall_seconds=600)
        assert a.inventory.version == b.inventory.version

    def test_a_different_seed_explores_differently(self, deposited, table):
        a = optimise_library(deposited, table, mode="anneal", seed=7, max_proposals=20,
                             wall_seconds=600)
        b = optimise_library(deposited, table, mode="anneal", seed=8, max_proposals=20,
                             wall_seconds=600)
        assert a.inventory.version != b.inventory.version

    def test_the_result_never_scores_below_the_incumbent(self, deposited, table):
        """Best-so-far is kept apart from the accepted state, so wandering cannot lose ground."""
        for mode in MODES:
            got = optimise_library(deposited, table, mode=mode, seed=42, max_proposals=20,
                                   wall_seconds=600)
            assert (got.final_objectives["scalar"]
                    >= got.incumbent_objectives["scalar"] - 1e-12)

    def test_the_output_still_compiles(self, deposited, table):
        got = optimise_library(deposited, table, mode="greedy", seed=42, max_proposals=20,
                               wall_seconds=600)
        compiled = inv.compile_target(got.inventory, "AAAAUGUGG")
        assert compiled["available"] and compiled["product_nt"] == 901
