"""Versioned inventories: loading, validation, compilation, versioning and reload.

The tests needing real kit sequences skip when Supplementary Table S1 has not been supplied,
because the package deliberately does not ship it -- see `NOTICE.md`. The structural tests
build their own small inventories and always run.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from clippr import inventories as inv
from clippr.parts import FUSION_SITES, inventory as kit_index

TABLE_S1 = Path(__file__).resolve().parents[1] / "data" / "grasp_supp" / "Table S1.xlsx"


@pytest.fixture(scope="module")
def deposited():
    if not TABLE_S1.is_file():
        pytest.skip(f"supply {TABLE_S1.name} to run inventory tests")
    return inv.load_deposited(TABLE_S1)


def _record(module_id="m1", dna="AATGCCCGGGACTC", frame=1, **kw):
    return inv.ModuleRecord(
        module_id=module_id, version=kw.get("version", "test"), dna=dna, frame=frame,
        block=kw.get("block", "1A"), fifth=kw.get("fifth"), last=kw.get("last"),
        five_interface=dna[:4], three_interface=dna[-4:],
        source_sha256=kw.get("source_sha256", "0" * 64))


class TestModuleRecord:
    def test_the_coding_interval_is_whole_codons_from_the_frame(self):
        r = _record(dna="A" * 10, frame=1)
        assert r.coding_interval == (1, 10)      # 9 bases = 3 codons

    def test_a_frame_that_leaves_a_partial_codon_truncates_it(self):
        r = _record(dna="A" * 11, frame=1)
        assert r.coding_interval == (1, 10)

    def test_identity_is_the_sequence_not_the_label(self):
        assert _record(dna="ACGTACGTACGT").sha256 == _record(dna="ACGTACGTACGT").sha256
        assert _record(dna="ACGTACGTACGT").sha256 != _record(dna="ACGTACGTACGA").sha256


class TestVersioning:
    def test_the_version_follows_the_contents(self):
        a = inv.Inventory(label="x", interface_assignment="d", context={},
                          modules={"m1": _record()})
        b = inv.Inventory(label="different-label", interface_assignment="d", context={},
                          modules={"m1": _record()})
        assert a.version == b.version, "the label is not part of the content identity"

    def test_a_changed_sequence_changes_the_version(self):
        a = inv.Inventory(label="x", interface_assignment="d", context={},
                          modules={"m1": _record(dna="AATGCCCGGGACTC")})
        b = inv.Inventory(label="x", interface_assignment="d", context={},
                          modules={"m1": _record(dna="AATGCCCGGTACTC")})
        assert a.version != b.version

    def test_a_changed_interface_assignment_changes_the_version(self):
        a = inv.Inventory(label="x", interface_assignment="deposited", context={},
                          modules={"m1": _record()})
        b = inv.Inventory(label="x", interface_assignment="redesigned", context={},
                          modules={"m1": _record()})
        assert a.version != b.version


class TestDerive:
    def test_an_untouched_module_is_carried_over_as_unchanged(self):
        base = inv.Inventory(label="base", interface_assignment="d", context={},
                             modules={"m1": _record()})
        got = inv.derive(base, "next", {})
        assert got.modules["m1"].unchanged is True
        assert got.modules["m1"].dna == base.modules["m1"].dna

    def test_a_replaced_module_records_what_it_came_from(self):
        base = inv.Inventory(label="base", interface_assignment="d", context={},
                             modules={"m1": _record(dna="AATGCCCGGGACTC")})
        got = inv.derive(base, "next", {"m1": "AATGCCCGGTACTC"})
        assert got.modules["m1"].unchanged is False
        assert got.modules["m1"].derived_from.startswith("test:")

    def test_replacing_with_the_identical_sequence_counts_as_unchanged(self):
        base = inv.Inventory(label="base", interface_assignment="d", context={},
                             modules={"m1": _record(dna="AATGCCCGGGACTC")})
        got = inv.derive(base, "next", {"m1": "AATGCCCGGGACTC"})
        assert got.modules["m1"].unchanged is True


@pytest.mark.usefixtures("deposited")
class TestDepositedInventory:
    def test_every_indexed_module_loads(self, deposited):
        assert len(deposited) == len(kit_index()) == 42

    def test_it_validates_clean(self, deposited):
        assert inv.validate(deposited) == []

    def test_every_module_has_a_derived_frame(self, deposited):
        assert all(r.frame in (0, 1, 2) for r in deposited.modules.values())

    def test_every_module_has_a_witness_context(self, deposited):
        """40 of 42 is not 42: the two AGGT start variants need their own witness."""
        witnesses = inv.witness_targets()
        assert set(witnesses) >= {m.plasmid_id for m in kit_index()}

    def test_both_fusion_site_variants_are_reachable(self, deposited):
        first = {site: inv.compile_target(deposited, "AAAAUGUGG",
                                          fusion_site=site)["modules"][0][0]
                 for site in FUSION_SITES}
        assert first["AATG"] != first["AGGT"]
        assert first["AATG"].endswith("AATG") and first["AGGT"].endswith("AGGT")


@pytest.mark.usefixtures("deposited")
class TestCompile:
    def test_a_compiled_product_joins_at_the_shared_overhangs(self, deposited):
        got = inv.compile_target(deposited, "AAAAUGUGG")
        total = sum(len(deposited.modules[m].dna) for m, _v in got["modules"])
        assert got["product_nt"] == total - 4 * (len(got["modules"]) - 1)

    def test_the_scope_says_it_is_not_an_expression_construct(self, deposited):
        assert "not an expression construct" in inv.compile_target(
            deposited, "AAAAUGUGG")["boundary_scope"]

    def test_an_unbuildable_target_is_reported_not_raised(self, deposited):
        got = inv.compile_target(deposited, "ACGU")
        assert got["available"] is False and got["reason"]

    def test_a_missing_module_names_itself(self, deposited):
        broken = inv.derive(deposited, "broken", {})
        needed = [m for m, _v in inv.compile_target(deposited, "AAAAUGUGG")["modules"]]
        del broken.modules[needed[0]]
        with pytest.raises(inv.IncompatibleInventory, match=needed[0]):
            inv.compile_target(broken, "AAAAUGUGG")

    def test_mixed_interface_assignments_are_refused(self, deposited):
        """Modules whose boundaries do not meet cannot ligate, so they must not compile."""
        mixed = inv.derive(deposited, "mixed", {})
        needed = [m for m, _v in inv.compile_target(deposited, "AAAAUGUGG")["modules"]]
        victim = mixed.modules[needed[1]]
        mixed.modules[needed[1]] = inv.ModuleRecord(
            **{**victim.__dict__, "dna": "TTTT" + victim.dna[4:], "five_interface": "TTTT"})
        with pytest.raises(inv.IncompatibleInventory, match="interface assignment"):
            inv.compile_target(mixed, "AAAAUGUGG")


@pytest.mark.usefixtures("deposited")
class TestSaveAndLoad:
    def test_a_saved_inventory_reloads_identically(self, deposited, tmp_path):
        path = inv.save(deposited, tmp_path / "inv.json")
        back = inv.load(path)
        assert back.version == deposited.version
        assert all(back.modules[m].dna == deposited.modules[m].dna
                   for m in deposited.modules)

    def test_a_reloaded_inventory_compiles_the_same_product(self, deposited, tmp_path):
        back = inv.load(inv.save(deposited, tmp_path / "inv.json"))
        assert (inv.compile_target(back, "AAAAUGUGG")["product"]
                == inv.compile_target(deposited, "AAAAUGUGG")["product"])

    def test_an_edited_file_is_refused(self, deposited, tmp_path):
        """The recorded version must agree with what the records hash to."""
        path = inv.save(deposited, tmp_path / "inv.json")
        data = json.loads(path.read_text(encoding="utf-8"))
        victim = sorted(data["modules"])[0]
        dna = data["modules"][victim]["dna"]
        data["modules"][victim]["dna"] = ("T" if dna[0] != "T" else "G") + dna[1:]
        path.write_text(json.dumps(data), encoding="utf-8")
        with pytest.raises(inv.InventoryInvalid, match="edited since it was written"):
            inv.load(path)

    def test_an_unreadable_schema_is_refused(self, deposited, tmp_path):
        path = tmp_path / "inv.json"
        path.write_text(json.dumps({"schema_version": 999}), encoding="utf-8")
        with pytest.raises(inv.InventoryInvalid, match="schema"):
            inv.load(path)
