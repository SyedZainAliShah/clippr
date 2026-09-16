"""What the run manifest must record, and must not quietly omit."""
from __future__ import annotations

import json

import pytest

from clippr import design_oneshot, manifest
from test_design import TOY_TABLE


@pytest.fixture(scope="module")
def design():
    return design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE, check_offtarget=False, seed=42)


class TestCodonTableIdentity:
    def test_effective_hash_is_the_table_the_optimiser_received(self, design):
        """Definitional: `effective_sha256` must be `complete_table`'s output, not the input.

        For an already-complete DNA-spelled table the two are equal, so this cannot be
        checked by asserting they differ -- an earlier version of this test did exactly
        that and failed.
        """
        import hashlib
        from clippr.codons import complete_table
        rec = manifest.codon_table_record(design["codon_table"], 1)
        want = hashlib.sha256(json.dumps(complete_table(design["codon_table"], 1),
                                         sort_keys=True,
                                         separators=(",", ":")).encode()).hexdigest()
        assert rec["effective_sha256"] == want

    def test_completion_is_visible_when_the_table_omits_a_codon(self, design):
        """The case where the optimiser's table really does differ from the supplied one."""
        partial = {aa: dict(codons) for aa, codons in design["codon_table"].items()}
        aa = next(a for a, codons in partial.items() if len(codons) > 1)
        del partial[aa][sorted(partial[aa])[0]]
        rec = manifest.codon_table_record(partial, 1)
        assert rec["supplied"]["canonical_sha256"] != rec["effective_sha256"]

    def test_a_changed_table_changes_the_recorded_hash(self, design):
        """A hash that does not move when the input moves is recording nothing."""
        before = manifest.build(design)["inputs"]["codon_table"]["effective_sha256"]
        altered = {aa: dict(codons) for aa, codons in design["codon_table"].items()}
        first = sorted(altered)[0]
        codon = sorted(altered[first])[0]
        altered[first][codon] = altered[first][codon] / 2 + 0.01
        moved = dict(design, codon_table=altered)
        after = manifest.build(moved)["inputs"]["codon_table"]["effective_sha256"]
        assert before != after

    def test_source_file_is_hashed_when_given(self, design, tmp_path):
        src = tmp_path / "table.json"
        src.write_text(json.dumps(design["codon_table"]), encoding="utf-8")
        rec = manifest.build(design, codon_table_source=src)["inputs"]["codon_table"]
        assert rec["supplied"]["path"] == str(src)
        assert rec["supplied"]["sha256"]

    def test_taxon_and_compartment_are_recorded(self, design):
        rec = manifest.build(design, taxid=3055,
                             compartment="nuclear")["inputs"]["codon_table"]
        assert (rec["taxid"], rec["compartment"], rec["genetic_code"]) == (3055, "nuclear", 1)


class TestScreening:
    def test_a_skipped_screen_is_not_a_clean_screen(self, design):
        """`offtarget=None` alone cannot distinguish 'disabled' from 'found nothing'."""
        rec = manifest.build(design)["screening"]
        assert rec["assessed"] is False
        assert any("not assessed" in s for s in rec["status"])
        assert "reference_genome" not in rec, "named a reference that was never read"


class TestScreeningStates:
    """An attempted-and-failed screen is not an assessed one."""

    def _unavailable(self, design):
        off = {"verdict": "not checked", "error": "OSError: reference unavailable",
               "n_transcript": None, "n_genomic": None, "genes": []}
        return dict(design, offtarget=off)

    def test_an_unavailable_screen_is_not_assessed(self, design):
        """It carries a result dict like a successful scan; treating that as evidence of
        assessment printed `assessed: true` beside the status 'not assessed'."""
        rec = manifest.build(self._unavailable(design))["screening"]
        assert rec["assessed"] is False
        assert rec["designs_unavailable"] == 1 and rec["designs_screened"] == 0
        assert rec["unavailable_errors"]
        assert "reference_genome" not in rec

    def test_a_mixed_batch_counts_each_state(self, design):
        rec = manifest.build([design, self._unavailable(design)])["screening"]
        assert (rec["designs_disabled"], rec["designs_unavailable"],
                rec["designs_screened"]) == (1, 1, 0)
        assert rec["assessed"] is False


class TestRunIdentity:
    def test_identity_notices_an_edit_to_the_source(self, design, monkeypatch):
        """HEAD plus dirty *filenames* cannot tell two edits of one file apart."""
        before = manifest.run_identity(design["codon_table"], 1)["source_sha256"]
        real = manifest.source_fingerprint
        monkeypatch.setattr(manifest, "source_fingerprint",
                            lambda: {"files": 0, "sha256": "changed"})
        after = manifest.run_identity(design["codon_table"], 1)["source_sha256"]
        monkeypatch.setattr(manifest, "source_fingerprint", real)
        assert before != after and len(before) == 64

    def test_identity_carries_the_table_and_scorer(self, design):
        ident = manifest.run_identity(design["codon_table"], 1)
        assert ident["matrix_fingerprint"] and ident["scorer_version"]
        assert ident["codon_table_effective_sha256"]


class TestCoverage:
    def test_completion_and_failures_are_carried(self, design):
        m = manifest.build(design, completion="partial", failures=["BADTARGET"])
        assert m["completion"] == "partial" and m["failures"] == ["BADTARGET"]

    def test_stage_seconds_sum_across_designs(self, design):
        m = manifest.build([design, design])
        assert m["configuration"]["n_designs"] == 2
        assert m["stage_seconds"]["optimization"] == pytest.approx(
            2 * design["timings"]["optimization"])

    def test_empty_run_is_refused(self):
        with pytest.raises(ValueError, match="at least one design"):
            manifest.build([])


class TestWrite:
    def test_written_manifest_round_trips(self, design, tmp_path):
        path = manifest.write(manifest.build(design), tmp_path / "pkg" / "manifest.json")
        assert json.loads(path.read_text(encoding="utf-8"))["schema_version"] == 1


class TestBatchHomogeneity:
    """A manifest states the codon table, genetic code and enzyme profile once.

    If the designs it describes disagree on any of them, that single statement is false for
    every design after the first, and no downstream reader could detect it.
    """

    def _pair(self, **override):
        from test_design import TOY_TABLE
        from clippr.design import design_oneshot

        a = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE, check_offtarget=False, seed=42)
        b = dict(a)
        b["target_rna"] = "GCUAAAGAC"
        b.update(override)
        return [a, b]

    def test_a_homogeneous_batch_is_accepted(self):
        from clippr import manifest

        record = manifest.build(self._pair())
        assert record["configuration"]["n_designs"] == 2

    def test_disagreeing_genetic_code_is_refused_and_named(self):
        from clippr import manifest

        with pytest.raises(ValueError, match="genetic_code"):
            manifest.build(self._pair(genetic_code=11))

    def test_disagreeing_enzyme_profile_is_refused_and_named(self):
        from clippr import manifest

        with pytest.raises(ValueError, match="enzyme_profile_effective"):
            manifest.build(self._pair(enzyme_profile_effective=("BsaI",)))

    def test_key_order_alone_is_not_a_disagreement(self):
        """A dict written in a different order is the same table, not a different one."""
        from clippr import manifest

        pair = self._pair()
        table = pair[0]["codon_table"]
        pair[1]["codon_table"] = {k: table[k] for k in reversed(list(table))}
        assert manifest.build(pair)["configuration"]["n_designs"] == 2


class TestSourceFingerprintScope:
    """The fingerprint must cover the shipped data, not only the code.

    `parts.json` and `scaffold.json` determine a design as much as any module does, so a
    fingerprint blind to them would report the same identity for two runs built on different
    scaffolds.
    """

    def test_shipped_data_is_in_scope(self):
        from clippr.manifest import source_fingerprint

        fingerprint = source_fingerprint()
        assert fingerprint["by_suffix"][".json"] >= 2
        assert fingerprint["files"] == sum(fingerprint["by_suffix"].values())

    def test_a_changed_data_file_changes_the_digest(self, tmp_path):
        """The falsifier: if this passes with the file absent, the scope claim is empty."""
        from pathlib import Path

        from clippr import manifest
        from clippr.manifest import source_fingerprint

        package = Path(manifest.__file__).resolve().parent
        before = source_fingerprint()["sha256"]
        probe = package / "_fingerprint_probe.json"
        try:
            probe.write_text('{"probe": 1}', encoding="utf-8")
            during = source_fingerprint()["sha256"]
        finally:
            probe.unlink(missing_ok=True)
        assert during != before
        assert source_fingerprint()["sha256"] == before
