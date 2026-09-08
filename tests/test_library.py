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
        oligos = lib.oligos()
        assert len(oligos) == sum(len(r["oligos"]) for r in lib.designs.values())
        assert "target_rna" in oligos.columns, "rows must say which design they belong to"

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

    def test_summary_mentions_the_saving(self, lib):
        assert "saved" in lib.summary()


class TestPooling:
    def test_pooling_is_cheaper_than_separate_pools(self, lib):
        """The reason library mode exists: the price is per pool, not per design."""
        c = lib.cost
        assert c["separate_pools_eur"] > c["pooled_eur"]
        assert c["saving_eur"] == c["separate_pools_eur"] - c["pooled_eur"]

    def test_pooled_price_does_not_grow_with_the_library(self, lib):
        small = design_library(TARGETS[:1], codon_table=TOY_TABLE, check_offtarget=False)
        assert small.cost["pooled_eur"] == lib.cost["pooled_eur"]

    def test_flagged_as_not_a_quote(self, lib):
        assert lib.cost["list_price_not_a_quote"] is True

    def test_counts_every_fragment(self, lib):
        assert lib.cost["n_oligos"] == len(lib.oligos())


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


class TestWrite:
    def test_writes_the_order_sheet_and_qc(self, lib, tmp_path):
        paths = write_library(lib, tmp_path)
        for key in ("order_sheet", "opool", "qc_table", "summary", "genbank_dir"):
            assert key in paths

    def test_opool_has_only_what_a_vendor_needs(self, lib, tmp_path):
        import pandas as pd
        paths = write_library(lib, tmp_path)
        pool = pd.read_csv(paths["opool"])
        assert list(pool.columns) == ["Pool name", "Sequence"]
        assert len(pool) == len(lib.oligos())

    def test_one_genbank_per_design(self, lib, tmp_path):
        from pathlib import Path
        paths = write_library(lib, tmp_path)
        assert len(list(Path(paths["genbank_dir"]).glob("*.gb"))) == len(TARGETS)


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
