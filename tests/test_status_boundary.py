"""Controls for the one readiness boundary: what may be written as a normal order file.

Every test here is a negative case that passed before the boundary was unified. The
project had two predicates: `library._ready` fell back to `status != "fail"`, while
`design_oneshot` chose the diagnostic export route only on `status == "fail"` and then
computed readiness from a stricter conjunction. Anything in between, `not_assessed`
above all, was refused by one and accepted by the other.

These use explicit verdict fixtures. They do not claim the default example naturally
returns these statuses; they claim that when it does, nothing orderable is written.
"""
from __future__ import annotations

import copy
import json

import pytest

from clippr import design as design_mod
from clippr.design import design_oneshot, synthesis_readiness
from clippr.library import LibraryResult, design_library, write_library

from test_design import TOY_TABLE

TARGETS = ["AAAAUGUGG", "GCUAAAGAC", "UUACACGUG"]


class _Absent:
    def __repr__(self):
        return "<absent>"


_ABSENT = _Absent()


@pytest.fixture(scope="module")
def real_designs():
    """Real designs, so every field a consumer reads is genuinely present."""
    lib = design_library(TARGETS, codon_table=TOY_TABLE, check_offtarget=False)
    assert not lib.failed
    return lib.designs


def _with_verdict(design, **overrides):
    """A real design with only its verdict fields replaced."""
    out = copy.deepcopy(design)
    for key, value in overrides.items():
        if value is _ABSENT:
            out.pop(key, None)
        elif key == "status":
            out.setdefault("final_artifacts", {})["status"] = value
        else:
            out[key] = value
    return out


def _accepted(design):
    """A fixture the boundary must accept, so the controls are not all negative."""
    out = _with_verdict(design, status="pass", constraints_ok=True,
                        synthesis_ready=_ABSENT)
    out["final_artifacts"]["hard"] = []
    return out


class TestReadinessPredicate:
    """The predicate itself, independent of any export."""

    def test_pass_with_constraints_is_ready(self, real_designs):
        ready, why = synthesis_readiness(_accepted(real_designs["AAAAUGUGG"]))
        assert ready and why == []

    @pytest.mark.parametrize("status", ["not_assessed", "fail", "unknown",
                                       "PASS", "", None])
    def test_only_known_accepted_states_are_ready(self, real_designs, status):
        """`status != "fail"` admitted every one of these except `fail`."""
        d = _with_verdict(real_designs["AAAAUGUGG"], status=status,
                          constraints_ok=True, synthesis_ready=_ABSENT)
        d["final_artifacts"]["hard"] = []
        ready, why = synthesis_readiness(d)
        assert not ready, f"{status!r} must not be orderable"
        assert why, "a refusal must say why"

    def test_absent_assessment_is_not_ready(self, real_designs):
        d = _with_verdict(real_designs["AAAAUGUGG"], final_artifacts=_ABSENT,
                          synthesis_ready=_ABSENT, constraints_ok=True)
        ready, why = synthesis_readiness(d)
        assert not ready
        assert any("no final-artifact assessment" in w for w in why)

    def test_unsatisfied_optimiser_is_not_ready(self, real_designs):
        """The fallback path ignored this; the strict path did not."""
        d = _with_verdict(real_designs["AAAAUGUGG"], status="pass",
                          constraints_ok=False, synthesis_ready=_ABSENT)
        d["final_artifacts"]["hard"] = []
        ready, why = synthesis_readiness(d)
        assert not ready
        assert any("optimiser" in w for w in why)

    def test_absent_constraint_result_is_not_ready(self, real_designs):
        d = _with_verdict(real_designs["AAAAUGUGG"], status="pass",
                          constraints_ok=_ABSENT, synthesis_ready=_ABSENT)
        d["final_artifacts"]["hard"] = []
        ready, why = synthesis_readiness(d)
        assert not ready

    def test_a_stale_ready_flag_cannot_override_the_evidence(self, real_designs):
        """A flag that disagrees with the record is itself a reason to refuse.

        The old `_ready` returned the flag whenever it was present, so an edited or
        stale `synthesis_ready=True` was authoritative over a hard breach.
        """
        d = _with_verdict(real_designs["AAAAUGUGG"], status="fail",
                          constraints_ok=True, synthesis_ready=True)
        d["final_artifacts"]["hard"] = ["GC 0.660 outside (0.35, 0.65)"]
        ready, why = synthesis_readiness(d)
        assert not ready
        assert any("disagrees" in w for w in why)

    def test_a_hard_breach_with_a_pass_status_is_contradictory(self, real_designs):
        d = _with_verdict(real_designs["AAAAUGUGG"], status="pass",
                          constraints_ok=True, synthesis_ready=_ABSENT)
        d["final_artifacts"]["hard"] = ["GC 0.660 outside (0.35, 0.65)"]
        ready, why = synthesis_readiness(d)
        assert not ready


class TestLibraryExportBoundary:
    """What `write_library` is allowed to put on disk."""

    @staticmethod
    def _lib(designs):
        out = LibraryResult()
        out.designs = designs
        return out

    def test_not_assessed_writes_no_order_sheet(self, real_designs, tmp_path):
        """The finding's first control: readiness removed, status `not_assessed`."""
        d = _with_verdict(real_designs["AAAAUGUGG"], status="not_assessed",
                          synthesis_ready=_ABSENT, constraints_ok=True)
        lib = self._lib({"AAAAUGUGG": d})

        assert lib.accepted == {}
        assert set(lib.rejected) == {"AAAAUGUGG"}

        paths = write_library(lib, tmp_path)
        assert "order_sheet" not in paths
        assert "opool" not in paths
        assert not (tmp_path / "library_order_sheet.csv").exists()
        assert (tmp_path / "library_REJECTED_oligos.csv").exists()

    def test_unknown_status_writes_no_order_sheet(self, real_designs, tmp_path):
        d = _with_verdict(real_designs["AAAAUGUGG"], status="probably fine",
                          synthesis_ready=_ABSENT, constraints_ok=True)
        paths = write_library(self._lib({"AAAAUGUGG": d}), tmp_path)
        assert "order_sheet" not in paths

    def test_mixed_acceptance_exports_only_the_accepted(self, real_designs, tmp_path):
        import pandas as pd

        good = _accepted(real_designs["AAAAUGUGG"])
        bad = _with_verdict(real_designs["GCUAAAGAC"], status="not_assessed",
                            constraints_ok=True, synthesis_ready=_ABSENT)
        lib = self._lib({"AAAAUGUGG": good, "GCUAAAGAC": bad})

        assert set(lib.accepted) == {"AAAAUGUGG"}
        paths = write_library(lib, tmp_path)

        sheet = pd.read_csv(paths["order_sheet"])
        assert set(sheet["target_rna"]) == {"AAAAUGUGG"}
        rejected = pd.read_csv(paths["rejected_oligos"])
        assert set(rejected["target_rna"]) == {"GCUAAAGAC"}

    def test_the_manifest_states_rejection_coverage(self, real_designs, tmp_path):
        """`completion` counted crashes only, so an all-rejected run said "complete"."""
        d = _with_verdict(real_designs["AAAAUGUGG"], status="not_assessed",
                          synthesis_ready=_ABSENT, constraints_ok=True)
        write_library(self._lib({"AAAAUGUGG": d}), tmp_path)

        record = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
        assert record["completion"] != "complete"
        assert record["synthesis_coverage"] == {"requested": 1, "accepted": 0,
                                                "rejected": 1, "not_assessed": 1,
                                                "failed": 0}
        assert record["rejections"], "a rejected design must be named in the manifest"

    def test_an_accepted_library_manifest_is_complete(self, real_designs, tmp_path):
        """The positive side: the coverage block must not report a false rejection."""
        write_library(self._lib({"AAAAUGUGG": _accepted(real_designs["AAAAUGUGG"])}),
                      tmp_path)
        record = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
        assert record["completion"] == "complete"
        assert record["synthesis_coverage"]["accepted"] == 1
        assert record["rejections"] == []

    def test_the_summary_states_the_rejection(self, real_designs):
        d = _with_verdict(real_designs["AAAAUGUGG"], status="not_assessed",
                          synthesis_ready=_ABSENT, constraints_ok=True)
        text = self._lib({"AAAAUGUGG": d}).summary()
        assert "0 of 1" in text and "not orderable" in text

    def test_the_qc_row_is_scoped_and_carries_the_final_verdict(self, real_designs):
        """`qc` said PASS on a rejected design with no final column beside it."""
        d = _with_verdict(real_designs["AAAAUGUGG"], status="not_assessed",
                          synthesis_ready=_ABSENT, constraints_ok=True)
        table = self._lib({"AAAAUGUGG": d}).qc_table()

        assert "cds_qc" in table.columns, "the CDS-scoped label must say its scope"
        assert "final_artifacts" in table.columns
        assert "synthesis_ready" in table.columns
        assert table["final_artifacts"].iloc[0] == "not_assessed"
        assert bool(table["synthesis_ready"].iloc[0]) is False


class TestDirectDesignExportBoundary:
    """`design_oneshot` chose its export route on `status == "fail"` alone."""

    @staticmethod
    def _unassessed(monkeypatch):
        real = design_mod.final_artifact_report

        def patched(*a, **k):
            out = real(*a, **k)
            out["status"] = "not_assessed"
            out["hard"] = []
            return out

        monkeypatch.setattr(design_mod, "final_artifact_report", patched)

    def test_not_assessed_writes_no_normal_order_files(self, tmp_path, monkeypatch):
        self._unassessed(monkeypatch)
        out = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE,
                             check_offtarget=False, outdir=str(tmp_path))

        assert out["synthesis_ready"] is False
        assert "oligo_csv" not in out["paths"]
        assert "oligo_fasta" not in out["paths"]
        assert "diagnostic_oligo_csv" in out["paths"]
        assert all("REJECTED" in p.name for p in tmp_path.glob("*oligos.csv"))

    def test_the_refusal_reasons_are_recorded(self, tmp_path, monkeypatch):
        self._unassessed(monkeypatch)
        out = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE,
                             check_offtarget=False)
        assert out["not_ready_because"], "a refusal must carry its reasons"
