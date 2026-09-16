"""The six connected tasks and the result package they write."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from clippr import workflow as w
from clippr.ordering import HISTORICAL_OPOOL

TABLE_S1 = Path(__file__).resolve().parents[1] / "data" / "grasp_supp" / "Table S1.xlsx"
CODON_TABLE = Path(__file__).resolve().parents[1] / "data" / "codon_tables" / "kazusa_3055.json"
TARGETS = ["AAAAUGUGG", "UUACACGUGCGUAC"]


@pytest.fixture(scope="module")
def table():
    return json.loads(CODON_TABLE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def deposited_path():
    if not TABLE_S1.is_file():
        pytest.skip(f"supply {TABLE_S1.name} to run workflow tests on the real kit")
    return TABLE_S1


class TestPackage:
    def test_a_run_where_everything_failed_still_has_provenance(self, tmp_path):
        """A directory of failures with no record of what produced them is the bad case."""
        failed = w.WorkflowResult(task="t", ok=False, summary="all failed",
                                  failures=[{"target": "X", "reason": "no route"}])
        package = json.loads(
            w.write_package([failed], tmp_path).read_text(encoding="utf-8"))
        assert package["ok"] is False
        assert package["software"]["dependencies"]
        assert package["failures"] == [{"task": "t", "target": "X", "reason": "no route"}]

    def test_artefact_paths_are_relative_to_the_package(self, tmp_path):
        """Absolute paths pin a package to one machine and carry a home directory with it."""
        inside = tmp_path / "sub" / "thing.json"
        inside.parent.mkdir(parents=True)
        inside.write_text("{}", encoding="utf-8")
        result = w.WorkflowResult(task="t", ok=True, summary="",
                                  artefacts={"thing": str(inside)})
        package = json.loads(
            w.write_package([result], tmp_path).read_text(encoding="utf-8"))
        assert package["tasks"][0]["artefacts"]["thing"] == "sub/thing.json"

    def test_an_artefact_outside_the_package_is_copied_in(self, tmp_path):
        """Marking it and moving on left the package naming a file it did not contain."""
        outside = tmp_path / "elsewhere" / "thing.json"
        outside.parent.mkdir(parents=True)
        outside.write_text("{}", encoding="utf-8")
        result = w.WorkflowResult(task="t", ok=True, summary="",
                                  artefacts={"thing": str(outside)})
        root = tmp_path / "pkg"
        package = json.loads(w.write_package([result], root).read_text(encoding="utf-8"))
        rel = package["tasks"][0]["artefacts"]["thing"]
        assert not rel.startswith("[")
        assert (root / rel).is_file(), "the package names a file it does not contain"

    def test_a_missing_artefact_prevents_a_complete_success(self, tmp_path):
        result = w.WorkflowResult(task="t", ok=True, summary="",
                                  artefacts={"gone": str(tmp_path / "nope.json")})
        package = json.loads(
            w.write_package([result], tmp_path / "pkg").read_text(encoding="utf-8"))
        assert package["ok"] is False
        assert package["missing_artefacts"]

    def test_the_package_states_its_scope(self, tmp_path):
        package = json.loads(w.write_package([], tmp_path).read_text(encoding="utf-8"))
        assert "nothing here is biologically validated" in package["scope"]
        assert "no order was placed" in package["scope"]


class TestOrderTask:
    def test_it_runs_without_any_supplied_input(self, tmp_path):
        got = w.plan_order(["ACGT" * 30] * 6, HISTORICAL_OPOOL, tmp_path)
        assert got.ok and got.data["pool_count"] == 1

    def test_cost_is_unavailable_without_a_dated_price(self, tmp_path):
        got = w.plan_order(["ACGT" * 30] * 6, HISTORICAL_OPOOL, tmp_path)
        assert got.data["cost"]["available"] is False
        assert "eligible" in got.summary

    def test_it_never_implies_vendor_approval(self, tmp_path):
        got = w.plan_order(["ACGT" * 30] * 6, HISTORICAL_OPOOL, tmp_path)
        assert "not obtained" in got.data["vendor_approval"]

    def test_an_ineligible_order_reports_the_constraint(self, tmp_path):
        got = w.plan_order(["ACGT" * 30], HISTORICAL_OPOOL, tmp_path)
        assert not got.ok
        assert any("cannot fill even one pool" in f["reason"] for f in got.failures)

    def test_a_multi_pool_order_succeeds_and_costs_the_sum_of_its_pools(self, tmp_path):
        """The per-pool maximum was applied to the whole order, refusing every split."""
        from decimal import Decimal

        from clippr.ordering import PriceTier, ProductProfile

        profile = ProductProfile(
            name="t", vendor="v", region="", currency="EUR",
            min_oligo_nt=20, max_oligo_nt=100, min_oligos=2, max_oligos=4,
            price_tiers=(PriceTier(2, 4, Decimal("10")),),
            priced_on="2026-09-15", price_status="current")
        got = w.plan_order(["ACGT" * 5] * 9, profile, tmp_path)
        assert got.ok and got.data["pool_count"] == 3
        assert got.data["cost"]["total"] == "30.00"
        assert all(r["eligible"] for r in got.data["per_pool_eligibility"])


@pytest.mark.usefixtures("deposited_path", "table")
class TestConnectedTasks:
    def test_load_and_compile_reports_the_inventory_it_used(self, deposited_path, tmp_path):
        got = w.load_and_compile(deposited_path, TARGETS, tmp_path)
        assert got.ok
        assert got.data["inventory_version"] and got.data["modules"] == 42
        assert Path(got.artefacts["compiled_targets"]).is_file()

    def test_it_surfaces_synthesis_profile_notes(self, deposited_path, tmp_path):
        """The deposited kit has one module outside our profile; say so before ordering."""
        got = w.load_and_compile(deposited_path, TARGETS, tmp_path)
        assert "pPR-1_D_LD5T" in got.data["synthesis_profile_notes"]

    def test_an_unbuildable_target_is_a_failure_not_a_crash(self, deposited_path, tmp_path):
        got = w.load_and_compile(deposited_path, ["ACGU"], tmp_path)
        assert not got.ok and got.failures

    def test_recode_then_compile_chains_through_saved_artefacts(self, deposited_path,
                                                                table, tmp_path):
        recoded = w.recode_for_host(deposited_path, table, tmp_path / "r",
                                    label="t", seeds=1, wall_seconds=600)
        assert recoded.ok
        compiled = w.load_and_compile(recoded.artefacts["inventory"], TARGETS,
                                      tmp_path / "c")
        assert compiled.ok
        assert compiled.data["inventory_version"] == recoded.data["to_version"]

    def test_explore_warns_that_choosing_changes_the_interface_version(self, deposited_path,
                                                                       table, tmp_path):
        recoded = w.recode_for_host(deposited_path, table, tmp_path / "r",
                                    label="t", seeds=1, wall_seconds=600)
        got = w.explore_interfaces(recoded.artefacts["inventory"], TARGETS, table,
                                   tmp_path / "e", max_evaluations=8, wall_seconds=300)
        assert "must be recompiled" in got.data["selecting_an_alternative"]
        assert "not the globally optimal front" in got.data["scope"]

    def test_optimise_reports_returning_the_incumbent(self, deposited_path, table,
                                                      tmp_path):
        recoded = w.recode_for_host(deposited_path, table, tmp_path / "r",
                                    label="t", seeds=1, wall_seconds=600)
        got = w.optimise_collection(recoded.artefacts["inventory"], table, tmp_path / "o",
                                    mode="greedy", max_proposals=0, wall_seconds=0.0)
        assert "incumbent" in got.summary


@pytest.mark.usefixtures("deposited_path", "table")
class TestFrontBinding:
    """A front's objectives describe one inventory. Applied to another they are false.

    The concrete case: a front from the recoded inventory, selected against the deposited one,
    reported adaptation 0.672749 for sequences actually scoring 0.219935.
    """

    def _front(self, deposited_path, table, tmp_path):
        recoded = w.recode_for_host(deposited_path, table, tmp_path / "r",
                                    label="t", seeds=1, wall_seconds=600)
        front = w.explore_interfaces(recoded.artefacts["inventory"], TARGETS, table,
                                     tmp_path / "f", max_evaluations=10, wall_seconds=300)
        return recoded, front

    def test_a_front_records_what_it_was_measured_against(self, deposited_path, table,
                                                          tmp_path):
        recoded, front = self._front(deposited_path, table, tmp_path)
        data = json.loads(Path(front.artefacts["front"]).read_text(encoding="utf-8"))
        binding = data["binding"]
        assert binding["inventory_version"] == recoded.data["to_version"]
        for field in ("codon_table_sha256", "genetic_code", "matrix", "k", "source_sha256"):
            assert binding[field] is not None

    def test_a_mismatched_inventory_is_refused(self, deposited_path, table, tmp_path):
        from clippr import inventories as inv

        _recoded, front = self._front(deposited_path, table, tmp_path)
        other = inv.load_deposited(deposited_path)
        inv.save(other, tmp_path / "deposited.json")
        got = w.select_interface(tmp_path / "deposited.json", front.artefacts["front"],
                                 tmp_path / "sel", codon_table=table)
        assert not got.ok
        assert any("inventory_version" in f["reason"] for f in got.failures)

    def test_the_matching_inventory_is_accepted_and_recomputed(self, deposited_path, table,
                                                               tmp_path):
        recoded, front = self._front(deposited_path, table, tmp_path)
        got = w.select_interface(recoded.artefacts["inventory"], front.artefacts["front"],
                                 tmp_path / "sel", codon_table=table)
        assert got.ok
        assert got.data["objectives_recomputed"] is True
        assert got.data["objective_drift"] == {}

    def test_objectives_are_recomputed_not_copied(self, deposited_path, table, tmp_path):
        """A value carried across from the front is a value nobody checked against what shipped."""
        recoded, front = self._front(deposited_path, table, tmp_path)
        got = w.select_interface(recoded.artefacts["inventory"], front.artefacts["front"],
                                 tmp_path / "sel", codon_table=table)
        record = json.loads(Path(got.artefacts["selection"]).read_text(encoding="utf-8"))
        assert "objectives_as_recorded_in_the_front" in record
        assert record["objectives"] is not None


class TestPackageInstances:
    def test_two_runs_of_one_task_keep_separate_artefacts(self, tmp_path):
        """Keying on the task name alone let the second copy overwrite the first."""
        a, b = tmp_path / "a.json", tmp_path / "b.json"
        a.write_text('{"which": "first"}', encoding="utf-8")
        b.write_text('{"which": "second"}', encoding="utf-8")
        results = [w.WorkflowResult(task="same", ok=True, summary="", artefacts={"f": str(a)}),
                   w.WorkflowResult(task="same", ok=True, summary="", artefacts={"f": str(b)})]
        root = tmp_path / "pkg"
        package = json.loads(w.write_package(results, root).read_text(encoding="utf-8"))

        paths = [t["artefacts"]["f"] for t in package["tasks"]]
        assert paths[0] != paths[1], "one instance overwrote the other"
        assert json.loads((root / paths[0]).read_text())["which"] == "first"
        assert json.loads((root / paths[1]).read_text())["which"] == "second"
        assert package["artefact_collisions"] == []

    def test_artefacts_carry_content_hashes(self, tmp_path):
        a = tmp_path / "a.json"
        a.write_text("{}", encoding="utf-8")
        result = w.WorkflowResult(task="t", ok=True, summary="", artefacts={"f": str(a)})
        package = json.loads(
            w.write_package([result], tmp_path / "pkg").read_text(encoding="utf-8"))
        assert package["tasks"][0]["artefact_sha256"]["f"]


@pytest.mark.usefixtures("deposited_path")
class TestAssemblyReadyOrdering:
    def test_order_items_are_assembly_ready_by_default(self, deposited_path, tmp_path):
        """A bare insert is not an assembly substrate."""
        from clippr import inventories as inv

        inv.save(inv.load_deposited(deposited_path), tmp_path / "inv.json")
        got = w.order_items_for(tmp_path / "inv.json", tmp_path / "items")
        assert got.data["form"] == "assembly_ready"
        assert got.data["enzyme"] == "BbsI"

    def test_every_substrate_exposes_a_level_zero_end_pair(self, deposited_path, tmp_path):
        from clippr import inventories as inv

        inv.save(inv.load_deposited(deposited_path), tmp_path / "inv.json")
        got = w.order_items_for(tmp_path / "inv.json", tmp_path / "items")
        level0 = {"CTCA", "ACTC", "AAGA", "GCAC", "TGAA", "CGAG"}
        for five, three in got.data["exposed_ends"]:
            assert five in level0 and three in level0

    def test_substrates_differ_from_their_inserts(self, deposited_path, tmp_path):
        from clippr import inventories as inv

        inv.save(inv.load_deposited(deposited_path), tmp_path / "inv.json")
        items = json.loads(Path(
            w.order_items_for(tmp_path / "inv.json", tmp_path / "items")
            .artefacts["order_items"]).read_text(encoding="utf-8"))
        assert all(i["sequence"] != i["released"] for i in items)

    def test_raw_inserts_remain_available_and_are_labelled(self, deposited_path, tmp_path):
        from clippr import inventories as inv

        inv.save(inv.load_deposited(deposited_path), tmp_path / "inv.json")
        got = w.order_items_for(tmp_path / "inv.json", tmp_path / "raw",
                                form="raw_inserts")
        items = json.loads(Path(got.artefacts["order_items"]).read_text(encoding="utf-8"))
        assert all(i["form"] == "raw_insert" for i in items)
        assert all("not an assembly substrate" in i["note"] for i in items)

    def test_an_unknown_form_is_refused(self, deposited_path, tmp_path):
        from clippr import inventories as inv

        inv.save(inv.load_deposited(deposited_path), tmp_path / "inv.json")
        with pytest.raises(ValueError, match="assembly_ready"):
            w.order_items_for(tmp_path / "inv.json", tmp_path / "x", form="guesswork")


@pytest.mark.usefixtures("deposited_path", "table")
class TestSearchContextIsPropagated:
    """A front records the scoring knobs it used; selection must honour them.

    Validating a `k=10` front with the default `k=20` refused a perfectly valid search, and
    `allow_mismatch=True` suppresses the check rather than honouring the context.
    """

    def test_a_non_default_k_round_trips_without_allow_mismatch(self, deposited_path,
                                                                table, tmp_path):
        recoded = w.recode_for_host(deposited_path, table, tmp_path / "r",
                                    label="t", seeds=1, wall_seconds=600)
        front = w.explore_interfaces(recoded.artefacts["inventory"], ["AAAAUGUGG"], table,
                                     tmp_path / "f", max_evaluations=8, wall_seconds=300,
                                     k=10)
        binding = json.loads(
            Path(front.artefacts["front"]).read_text(encoding="utf-8"))["binding"]
        assert binding["k"] == 10

        got = w.select_interface(recoded.artefacts["inventory"], front.artefacts["front"],
                                 tmp_path / "sel", codon_table=table)
        assert got.ok, got.failures
        assert got.data["objective_drift"] == {}


@pytest.mark.usefixtures("deposited_path")
class TestOrderRejectsContractBreaches:
    def test_a_failing_substrate_makes_the_order_step_fail(self, deposited_path, tmp_path,
                                                           monkeypatch):
        """Exporting a sequence outside the band while reporting success is the defect."""
        from clippr import inventories as inv
        from clippr import substrates as sub

        inv.save(inv.load_deposited(deposited_path), tmp_path / "inv.json")
        real = sub.build

        def failing(insert, block, module_id="", version=""):
            built = real(insert, block, module_id, version)
            return sub.Substrate(**{**built.__dict__,
                                    "synthesis_problems": ("GC 0.660 outside band",)})

        monkeypatch.setattr(sub, "build", failing)
        got = w.order_items_for(tmp_path / "inv.json", tmp_path / "items")
        assert not got.ok
        assert any(f["stage"] == "synthesis" for f in got.failures)
