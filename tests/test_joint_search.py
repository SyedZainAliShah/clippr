"""The joint junction/codon search: its contract, its honesty, and its front."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from clippr import inventories as inv, joint_search as js

TABLE_S1 = Path(__file__).resolve().parents[1] / "data" / "grasp_supp" / "Table S1.xlsx"
CODON_TABLE = Path(__file__).resolve().parents[1] / "data" / "codon_tables" / "kazusa_3055.json"
TARGETS = ["AAAAUGUGG", "UUACACGUGCGUAC"]


@pytest.fixture(scope="module")
def deposited():
    if not TABLE_S1.is_file():
        pytest.skip(f"supply {TABLE_S1.name} to run joint search tests")
    return inv.load_deposited(TABLE_S1)


@pytest.fixture(scope="module")
def table():
    from clippr.codons import complete_table
    return complete_table(json.loads(CODON_TABLE.read_text(encoding="utf-8")), 1)


@pytest.fixture(scope="module")
def classes(deposited):
    return js.junction_classes(deposited, TARGETS)


class TestDeclaredContract:
    def test_objective_directions_are_declared_not_inferred(self):
        assert js.DIRECTIONS == {"fidelity": "max", "adaptation": "max",
                                 "repeat_burden": "min"}


@pytest.mark.usefixtures("deposited", "table", "classes")
class TestClasses:
    def test_roles_sharing_a_context_are_one_class(self, classes):
        """1A->B, 2A->B, 14A->B and 19A->B are one physical interface, so they move together."""
        by_key = {c.key: c for c in classes}
        assert len(by_key["EL@2:ACTC"].roles) > 1

    def test_every_class_offers_its_incumbent(self, classes):
        for cls in classes:
            assert cls.deposited in dict(cls.options)


@pytest.mark.usefixtures("deposited", "table", "classes")
class TestApplyAssignment:
    def test_the_incumbent_assignment_changes_nothing(self, deposited, classes):
        incumbent = {c.key: c.deposited for c in classes}
        rebuilt = js.apply_assignment(deposited, classes, incumbent, "same")
        assert all(rebuilt.modules[m].dna == deposited.modules[m].dna
                   for m in deposited.modules)

    def test_changing_a_class_rewrites_both_sides_of_the_junction(self, deposited, classes):
        """Rewriting one side would leave the two modules unable to anneal."""
        cls = next(c for c in classes if c.key == "AR@2:AAGA")
        other = next(o for o, _c in cls.options if o != cls.deposited)
        assignment = {c.key: c.deposited for c in classes}
        assignment[cls.key] = other
        rebuilt = js.apply_assignment(deposited, classes, assignment, "changed")

        left = [m for m, r in deposited.modules.items() if r.three_interface == cls.deposited]
        right = [m for m, r in deposited.modules.items() if r.five_interface == cls.deposited]
        assert left and right
        assert all(rebuilt.modules[m].three_interface == other for m in left)
        assert all(rebuilt.modules[m].five_interface == other for m in right)

    def test_a_changed_assignment_still_compiles(self, deposited, classes):
        cls = next(c for c in classes if c.key == "AR@2:AAGA")
        other = next(o for o, _c in cls.options if o != cls.deposited)
        assignment = {c.key: c.deposited for c in classes}
        assignment[cls.key] = other
        rebuilt = js.apply_assignment(deposited, classes, assignment, "changed")
        compiled = inv.compile_target(rebuilt, "AAAAUGUGG")
        assert compiled["available"] and compiled["product_nt"] == 901


@pytest.mark.usefixtures("deposited", "table", "classes")
class TestEvaluate:
    def test_the_incumbent_is_feasible(self, deposited, classes, table):
        incumbent = {c.key: c.deposited for c in classes}
        got = js.evaluate(deposited, classes, incumbent, TARGETS, table)
        assert got.feasible, got.reason
        assert set(js.DIRECTIONS) <= set(got.objectives)

    def test_every_reaction_reports_its_stage_and_enzyme(self, deposited, classes, table):
        incumbent = {c.key: c.deposited for c in classes}
        got = js.evaluate(deposited, classes, incumbent, TARGETS, table)
        assert set(got.per_reaction_fidelity) == set(TARGETS)
        for rows in got.per_reaction_fidelity.values():
            assert all(r["stage"] in ("level0", "level1") for r in rows)
            assert all(r["enzyme"] in ("BbsI", "BsaI") for r in rows)

    def test_level_zero_is_scored_with_bbsi(self, deposited, classes, table):
        """All 42 deposited plasmids carry two BbsI sites and no BsaI site."""
        incumbent = {c.key: c.deposited for c in classes}
        got = js.evaluate(deposited, classes, incumbent, TARGETS, table)
        level0 = [r for rows in got.per_reaction_fidelity.values() for r in rows
                  if r["stage"] == "level0"]
        assert level0
        assert all(r["enzyme"] == "BbsI" and r["matrix"] == "BbsI-HF" for r in level0)

    def test_level_one_is_not_given_a_fidelity(self, deposited, classes, table):
        """Corrected 2026-09-16, twice, in opposite directions.

        It first asserted level 1 was unscorable because "its ends come from a block plasmid
        this package does not hold" -- wrong; the block geometry is established. It then
        asserted a fidelity -- also wrong, because ligation fidelity is a property of the whole
        competing overhang set and the deposited BsaI reaction co-assembles parts this package
        does not compile.

        What holds: the geometry is known, the participants are not, and a PPR-only figure is
        a diagnostic in its own field.
        """
        incumbent = {c.key: c.deposited for c in classes}
        got = js.evaluate(deposited, classes, incumbent, TARGETS, table)
        level1 = [r for rows in got.per_reaction_fidelity.values() for r in rows
                  if r["stage"] == "level1"]
        assert level1
        assert all(r["enzyme"] == "BsaI" for r in level1)
        assert all(r["fidelity"] is None for r in level1)
        assert all(r["geometry_established"] and not r["participants_established"]
                   for r in level1)
        assert all(r["block_subset"]["matrix"] == "BsaI-HFv2" for r in level1)
        assert got.unscorable_reactions

    def test_the_level_one_score_is_never_invented(self, deposited, classes, table):
        """A number must name the set it came from, whether it is a score or a diagnostic."""
        incumbent = {c.key: c.deposited for c in classes}
        got = js.evaluate(deposited, classes, incumbent, TARGETS, table)
        for rows in got.per_reaction_fidelity.values():
            for reaction in rows:
                if reaction["fidelity"] is not None:
                    assert reaction.get("overhangs"), reaction
                if reaction.get("block_subset"):
                    assert reaction["block_subset"]["overhangs"]

    def test_the_summary_is_the_minimum_over_scorable_reactions(self, deposited, classes,
                                                                table):
        """A product across reactions reads as a predicted yield, and is not one."""
        incumbent = {c.key: c.deposited for c in classes}
        got = js.evaluate(deposited, classes, incumbent, TARGETS, table)
        scored = [r["fidelity"] for rows in got.per_reaction_fidelity.values()
                  for r in rows if r["fidelity"] is not None]
        assert got.objectives["fidelity"] == pytest.approx(min(scored))

    def test_repeat_burden_is_measured_on_what_is_ordered(self, deposited, classes, table):
        """The objective was documented as order sequences and computed over products."""
        from clippr.objectives import duplicated_kmers

        incumbent = {c.key: c.deposited for c in classes}
        got = js.evaluate(deposited, classes, incumbent, TARGETS, table)
        distinct = {r.dna for r in deposited.modules.values()}
        assert got.objectives["repeat_burden"] == sum(
            duplicated_kmers(d, k=20) for d in distinct)
        assert "product_repeat_burden_diagnostic" in got.objectives

    def test_an_invalid_junction_set_is_infeasible_with_a_reason(self, deposited, classes,
                                                                 table):
        assignment = {c.key: c.deposited for c in classes}
        first = classes[0]
        assignment[first.key] = "AATT"          # palindromic
        got = js.evaluate(deposited, classes, assignment, TARGETS, table)
        assert not got.feasible and "palindromic" in got.reason

    def test_only_level_zero_classes_are_searched(self, deposited):
        """A block-join class cannot change any reaction this package can score."""
        level0 = js.junction_classes(deposited, TARGETS)
        everything = js.junction_classes(deposited, TARGETS, level0_only=False)
        assert all(c.stage == "level0" for c in level0)
        assert any(c.stage == "level1" for c in everything)
        assert len(level0) < len(everything)

    def test_evaluate_rejects_a_new_failure_on_an_already_failing_module(
            self, deposited, classes, table):
        """The guard must be exercised through `evaluate`, not through a private re-check.

        Asserting that my own list of introduced messages is non-empty tests my list. This
        drives a protein-preserving candidate through the production function and requires it
        to be refused for the *new* violation, while the inherited one is still excused.
        """
        from unittest.mock import patch

        from Bio.Data import CodonTable
        from Bio.Seq import Seq

        from clippr.recoding import module_constraint_problems

        offender = next(m for m, r in deposited.modules.items()
                        if module_constraint_problems(r.dna))
        record = deposited.modules[offender]
        start, end = record.coding_interval

        synonyms = {}
        for codon, aa in CodonTable.unambiguous_dna_by_id[1].forward_table.items():
            synonyms.setdefault(aa, []).append(codon)

        # Drive GC up synonymously, away from the inherited homopolymer, until a GC window
        # breaches the band. Protein and both interfaces are untouched.
        run = record.dna.index("TTTTT")
        dna = record.dna
        for i in range(start, end, 3):
            if i < 4 or i + 3 > len(dna) - 4 or (i < run + 5 and i + 3 > run):
                continue
            aa = str(Seq(dna[i:i + 3]).translate())
            richest = max(synonyms[aa], key=lambda c: c.count("G") + c.count("C"))
            dna = dna[:i] + richest + dna[i + 3:]

        assert (str(Seq(dna[start:end]).translate())
                == str(Seq(record.dna[start:end]).translate())), "probe changed the protein"
        inherited = set(module_constraint_problems(record.dna))
        introduced = [p for p in module_constraint_problems(dna) if p not in inherited]
        assert introduced, "the probe failed to introduce a new problem"

        mutated = inv.derive(deposited, "probe", {offender: dna})
        incumbent = {c.key: c.deposited for c in classes}
        with patch.object(js, "apply_assignment", return_value=mutated):
            got = js.evaluate(deposited, classes, incumbent, TARGETS, table)

        assert not got.feasible, "evaluate accepted a newly introduced violation"
        assert offender in got.reason


@pytest.mark.usefixtures("deposited", "table")
class TestSearch:
    def test_it_is_deterministic(self, deposited, table):
        a = js.search(deposited, TARGETS, table, max_evaluations=12, wall_seconds=300)
        b = js.search(deposited, TARGETS, table, max_evaluations=12, wall_seconds=300)
        assert a["observed_front"] == b["observed_front"]

    def test_the_front_holds_only_feasible_candidates(self, deposited, table):
        got = js.search(deposited, TARGETS, table, max_evaluations=20, wall_seconds=300)
        assert all(c["feasible"] for c in got["observed_front"])

    def test_it_never_calls_the_front_globally_optimal(self, deposited, table):
        got = js.search(deposited, TARGETS, table, max_evaluations=12, wall_seconds=300)
        assert "not the globally optimal front" in got["scope"]

    def test_a_one_evaluation_budget_reports_exhaustion(self, deposited, table):
        """'Ran out of budget' must never read as 'no improvement exists'."""
        got = js.search(deposited, TARGETS, table, max_evaluations=1, wall_seconds=300)
        assert got["budget_exhausted"] and got["completion"] == "budget_exhausted"

    def test_the_recommendation_respects_the_fidelity_tolerance(self, deposited, table):
        got = js.search(deposited, TARGETS, table, max_evaluations=30, wall_seconds=300)
        if not got["recommended"]:
            pytest.skip("no feasible candidate in this budget")
        best = max(c["objectives"]["fidelity"] for c in got["observed_front"])
        assert (got["recommended"]["objectives"]["fidelity"]
                >= best - js.FIDELITY_TOLERANCE)

    def test_degeneracy_is_reported_rather_than_dressed_as_a_trade_off(self, deposited,
                                                                      table):
        got = js.search(deposited, TARGETS, table, max_evaluations=20, wall_seconds=300)
        fidelities = {c["objectives"]["fidelity"]
                      for c in got["observed_front"]}
        assert got["fidelity_degenerate"] == (len(fidelities) <= 1)

    def test_cache_hits_are_not_counted_as_evaluations(self, deposited, table):
        got = js.search(deposited, TARGETS, table, max_evaluations=15, wall_seconds=300)
        assert got["evaluated"] == got["feasible"] + got["infeasible"]
        assert got["evaluated"] <= 15
