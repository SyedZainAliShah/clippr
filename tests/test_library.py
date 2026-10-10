"""Unit tests for library-scale design, and for the audit export."""
from __future__ import annotations

import pytest

from clippr.library import design_library, write_library

# The toy table from test_design, so these run without a network call for codon usage.
from test_design import TOY_TABLE

TARGETS = ["AAAAUGUGG", "GCUAAAGAC", "UUACACGUG"]


@pytest.fixture(scope="module")
def lib():
    return design_library(TARGETS, codon_table=TOY_TABLE, check_offtarget=False)


class TestLibrary:
    def test_designs_every_target(self, lib):
        assert set(lib.targets) == set(TARGETS)
        assert not lib.failed

    def test_one_combined_order_sheet(self, lib):
        """Counted over accepted designs: a rejected one is not an order row.

        `oligos()` used to concatenate every design regardless of verdict, which is how a
        rejected design reached a normal order sheet. The full set is still reachable for
        diagnostics through `include_rejected=True`.
        """
        every = lib.oligos(include_rejected=True)
        assert len(every) == sum(len(r["oligos"]) for r in lib.designs.values())
        assert "target_rna" in every.columns, "rows must say which design they belong to"
        assert len(lib.oligos()) == sum(len(r["oligos"]) for r in lib.accepted.values())

    def test_qc_table_has_one_row_per_design(self, lib):
        assert len(lib.qc_table()) == len(TARGETS)

    def test_a_failure_does_not_lose_the_rest(self):
        """One bad target in fifty should not cost you the other forty-nine."""
        out = design_library(TARGETS[:2] + ["NOTAVALIDTARGET"],
                             codon_table=TOY_TABLE, check_offtarget=False)
        assert len(out.designs) == 2
        assert "NOTAVALIDTARGET" in out.failed

    def test_progress_callback_fires(self):
        seen = []
        design_library(TARGETS[:2], codon_table=TOY_TABLE, check_offtarget=False,
                       on_progress=lambda i, n, t: seen.append((i, n, t)))
        assert [s[0] for s in seen] == [1, 2]

    def test_summary_reports_the_fragment_count(self, lib):
        assert "12 across all designs" in lib.summary()


class TestTotals:
    def test_counts_every_fragment(self, lib):
        """Totals describe the whole library, including designs that were rejected."""
        assert lib.totals["n_oligos"] == len(lib.oligos(include_rejected=True))

    def test_counts_every_design(self, lib):
        assert lib.totals["n_designs"] == len(lib.designs)


class TestCrosstalk:
    def test_reports_pairwise_separation(self, lib):
        assert "targets of length 9" in lib.crosstalk()

    def test_closest_pairs_are_ordered(self, lib):
        pairs = lib.closest_pairs()
        assert [p[2] for p in pairs] == sorted(p[2] for p in pairs)

    def test_only_compares_equal_lengths(self):
        """Hamming distance is undefined between targets of different length."""
        mixed = design_library(["AAAAUGUGG", "GCUAAAGACUUGCA"],
                               codon_table=TOY_TABLE, check_offtarget=False)
        assert mixed.closest_pairs() == []

    def test_two_tier_report_when_scores_are_supplied(self, lib, tmp_path):
        from clippr.crosstalk import load_ppr_scores

        p = tmp_path / "scores.tsv"
        p.write_text("For motif types: P\n5th/last\tA\tC\tG\tU\n"
                     "TN\t0.70\t-0.51\t0.53\t-0.22\nNN\t-0.44\t0.77\t-0.07\t0.62\n"
                     "TD\t0.01\t-1.0\t0.84\t-0.02\nND\t-0.53\t0.66\t0.55\t0.83\n",
                     encoding="utf-8")
        text = lib.crosstalk(scores=load_ppr_scores(p))
        assert "hamming" in text and "gates nothing" in text

    def test_a_weighted_metric_with_scores_is_refused_not_ignored(self, lib):
        """The two-tier report counts positions unweighted.

        Accepting both and honouring one is the silent-lie failure mode: the caller would
        believe a weighting had been applied that never was.
        """
        with pytest.raises(ValueError, match="would be ignored"):
            lib.crosstalk(metric="weighted", scores={"TN": {"A": 1.0}})


def _verdict(design: dict, *, ready: bool, assessed: bool = True) -> dict:
    """A copy of a real design carrying a stated final-artifact verdict.

    Only the decision is constructed. Everything the QC table, GenBank export and order
    sheet read is a genuine design, so these tests cannot pass against a fixture that has
    drifted from what `design_oneshot` returns.
    """
    out = dict(design)
    if assessed:
        out["synthesis_ready"] = ready
        out["final_artifacts"] = {
            "status": "pass" if ready else "fail",
            "hard": [] if ready else ["F1: GC 0.660 outside (0.35, 0.65) at 1"],
            "advisory": [],
        }
    else:
        out.pop("synthesis_ready", None)
        out.pop("final_artifacts", None)
    return out


def _library(**designs):
    from clippr.library import LibraryResult
    return LibraryResult(designs=dict(designs))


def _accepted_library(lib):
    """The first real design, marked accepted."""
    target, design = next(iter(lib.designs.items()))
    return _library(**{target: _verdict(design, ready=True)})


class TestWrite:
    def test_a_rejected_library_writes_no_order_sheet(self, lib, tmp_path):
        """Was `test_writes_the_order_sheet_and_qc`.

        Every design in `lib` is rejected by the final-artifact contract, so the normal
        order files must be absent and the rejected ones present. The analysis artefacts
        are still written, because a rejection is worth inspecting.
        """
        assert lib.accepted == {}, "if these start passing, this test has to change again"
        paths = write_library(lib, tmp_path)
        for key in ("qc_table", "summary", "genbank_dir",
                    "rejected_oligos", "rejected_reasons"):
            assert key in paths
        assert "order_sheet" not in paths and "opool" not in paths

    def test_an_accepted_library_writes_the_order_sheet_and_qc(self, lib, tmp_path):
        paths = write_library(_accepted_library(lib), tmp_path)
        for key in ("order_sheet", "opool", "qc_table", "summary", "genbank_dir"):
            assert key in paths
        assert "rejected_oligos" not in paths

    def test_opool_has_only_what_a_vendor_needs(self, lib, tmp_path):
        import pandas as pd
        lib = _accepted_library(lib)
        paths = write_library(lib, tmp_path)
        pool = pd.read_csv(paths["opool"])
        assert list(pool.columns) == ["Pool name", "Sequence"]
        assert len(pool) == len(lib.oligos())

    def test_one_genbank_per_design(self, lib, tmp_path):
        from pathlib import Path
        paths = write_library(lib, tmp_path)
        assert len(list(Path(paths["genbank_dir"]).glob("*.gb"))) == len(TARGETS)

    def test_the_package_carries_its_own_manifest(self, lib, tmp_path):
        """A directory of CSVs with no record of what produced them is not reproducible."""
        import json
        from pathlib import Path

        paths = write_library(lib, tmp_path)
        record = json.loads(Path(paths["manifest"]).read_text(encoding="utf-8"))
        # These targets are rejected by the final-artifact contract under the strict
        # policy, so "complete" would be false. `completion` used to count crashes only
        # and said "complete" on a library with nothing orderable in it.
        assert record["completion"] == "partial"
        assert record["synthesis_coverage"] == {"requested": len(TARGETS), "accepted": 0,
                                                "rejected": len(TARGETS),
                                                "not_assessed": 0, "failed": 0}
        assert {r["target"] for r in record["rejections"]} == set(TARGETS)
        assert record["configuration"]["n_designs"] == len(TARGETS)
        assert record["software"]["dependencies"]
        assert {d["target"] for d in record["designs"]} == set(TARGETS)

    def test_a_run_that_designed_nothing_still_records_what_ran(self, tmp_path):
        """Every target failed, so there is no codon table to report -- but a package with
        failures and no provenance statement is exactly what the manifest exists to prevent.
        """
        import json
        from pathlib import Path

        lib = design_library(["NOTATARGET"], codon_table=TOY_TABLE, check_offtarget=False)
        assert not lib.designs and lib.failed
        paths = write_library(lib, tmp_path)
        record = json.loads(Path(paths["manifest"]).read_text(encoding="utf-8"))
        assert record["completion"] == "failed"
        assert [f["target"] for f in record["failures"]] == ["NOTATARGET"]
        assert record["software"]["dependencies"]


class TestAuditExport:
    def test_audit_exports_to_csv(self, lib, tmp_path):
        audit = lib.designs["AAAAUGUGG"]["audit"]
        import pandas as pd
        df = pd.read_csv(audit.write_csv(tmp_path / "audit.csv"))
        assert len(df) == len(audit.overhang_decisions)

    def test_every_row_is_self_contained(self, lib):
        """A lab-notebook entry must not depend on remembering the active profile."""
        df = lib.designs["AAAAUGUGG"]["audit"].to_dataframe()
        for col in ("target_rna", "enzyme_profile", "enzymes_excluded", "organism",
                    "overhang", "status", "reason"):
            assert col in df.columns
        assert df["enzyme_profile"].nunique() == 1

    def test_statuses_are_the_documented_three(self, lib):
        df = lib.designs["AAAAUGUGG"]["audit"].to_dataframe()
        assert set(df["status"]) <= {"selected", "considered", "rejected"}


class TestRejectedDesignsNeverBecomeOrders:
    """A design that fails the final-artifact contract must not reach an order sheet.

    The public synthesis route used to bypass the gate entirely: `design_library` kept every
    design that did not raise, `oligos()` concatenated all of them, `write_library` exported
    them as a normal order sheet, and the workflow reported `ok=True`.
    """

    def test_an_all_rejected_library_exports_no_order_sheet(self, lib, tmp_path):
        target, design = next(iter(lib.designs.items()))
        rejected = _library(**{target: _verdict(design, ready=False)})
        assert rejected.accepted == {}
        paths = write_library(rejected, tmp_path)
        assert "order_sheet" not in paths and "opool" not in paths
        assert "rejected_oligos" in paths and "rejected_reasons" in paths
        assert not (tmp_path / "library_order_sheet.csv").exists()

    def test_a_mixed_library_orders_only_the_accepted_design(self, lib, tmp_path):
        (good, gd), (bad, bd) = list(lib.designs.items())[:2]
        mixed = _library(**{good: _verdict(gd, ready=True),
                            bad: _verdict(bd, ready=False)})
        write_library(mixed, tmp_path)
        sheet = (tmp_path / "library_order_sheet.csv").read_text(encoding="utf-8")
        assert good in sheet
        assert bad not in sheet, "a rejected design reached the normal order sheet"
        assert bad in (tmp_path / "library_REJECTED_oligos.csv").read_text(encoding="utf-8")

    def test_a_design_with_no_assessment_is_not_treated_as_ready(self, lib):
        target, design = next(iter(lib.designs.items()))
        unknown = _library(**{target: _verdict(design, ready=True, assessed=False)})
        assert unknown.accepted == {}, "absent information is not readiness"
        assert "no final-artifact assessment" in unknown.rejection_reasons()[target][0]

    def test_oligos_excludes_rejected_designs_by_default(self, lib):
        (good, gd), (bad, bd) = list(lib.designs.items())[:2]
        mixed = _library(**{good: _verdict(gd, ready=True),
                            bad: _verdict(bd, ready=False)})
        assert set(mixed.oligos()["target_rna"]) == {good}
        both = mixed.oligos(include_rejected=True)
        assert set(both["target_rna"]) == {good, bad}
        assert set(both["synthesis_ready"]) == {True, False}

    def test_the_public_workflow_does_not_call_a_rejected_library_complete(
            self, lib, tmp_path, monkeypatch):
        """Through `design_for_synthesis`, not just the direct result's keys."""
        from clippr import workflow as w
        target, design = next(iter(lib.designs.items()))
        rejected = _library(**{target: _verdict(design, ready=False)})
        monkeypatch.setattr("clippr.library.design_library", lambda *a, **k: rejected)
        result = w.design_for_synthesis([target], TOY_TABLE, tmp_path / "out")
        assert result.ok is False
        assert result.failures, "a rejected design must appear as a failure"
        assert "order_sheet" not in result.artefacts
        assert not (tmp_path / "out" / "library_order_sheet.csv").exists()
        assert result.data["rejected"] == [target]
        assert result.data["accepted"] == []
