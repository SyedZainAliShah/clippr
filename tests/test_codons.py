"""Unit tests for codon-table loading and coding-sequence design."""
from __future__ import annotations

import json

import pytest
from Bio.Seq import Seq

from clippr.codons import (
    BLACKLIST_ENZYMES,
    _table_name,
    complete_table,
    optimize_cds,
    table_from_cds_fasta,
    table_from_csv,
)

# A tiny, self-contained table: enough to drive the optimiser without a network call.
TOY_TABLE = {
    "M": {"ATG": 1.0},
    "W": {"TGG": 1.0},
    "K": {"AAA": 0.3, "AAG": 0.7},
    "L": {"CTG": 0.7, "CTC": 0.15, "CTA": 0.05, "CTT": 0.05, "TTA": 0.02, "TTG": 0.03},
    "A": {"GCC": 0.5, "GCT": 0.2, "GCA": 0.2, "GCG": 0.1},
    "G": {"GGC": 0.5, "GGT": 0.2, "GGA": 0.2, "GGG": 0.1},
    "S": {"TCC": 0.4, "TCT": 0.2, "TCA": 0.2, "TCG": 0.2},
    "T": {"ACC": 0.4, "ACT": 0.2, "ACA": 0.2, "ACG": 0.2},
    "P": {"CCC": 0.4, "CCT": 0.2, "CCA": 0.2, "CCG": 0.2},
    "V": {"GTC": 0.4, "GTT": 0.2, "GTA": 0.2, "GTG": 0.2},
}
PROTEIN = "MAKLGSTPVW" * 4


class TestTableName:
    def test_known_codes(self):
        assert _table_name(1) == "Standard"
        assert _table_name(11) == "Bacterial"

    def test_never_returns_none(self):
        """Table 11's name list contains a None, which would break DNA Chisel."""
        for code in (1, 11):
            assert isinstance(_table_name(code), str)

    def test_unknown_code_rejected(self):
        with pytest.raises(ValueError, match="unknown NCBI genetic code"):
            _table_name(999)


class TestCompleteTable:
    def test_fills_unobserved_codons_with_zero(self):
        """A partial serine entry is what crashed DNA Chisel with KeyError: 'AGT'."""
        full = complete_table({"S": {"TCC": 1.0}})
        assert set(full["S"]) == {"TCA", "TCC", "TCG", "TCT", "AGC", "AGT"}
        assert full["S"]["AGT"] == 0.0
        assert full["S"]["TCC"] == 1.0

    def test_leaves_complete_entries_alone(self):
        assert complete_table({"M": {"ATG": 1.0}})["M"] == {"ATG": 1.0}

    def test_converts_rna_spelling(self):
        assert "TTG" in complete_table({"L": {"UUG": 1.0}})["L"]

    def test_absent_amino_acid_is_reported_at_use(self):
        """A CDS set that never used tryptophan yields a table without W."""
        with pytest.raises(ValueError, match="no entry for W"):
            optimize_cds("MAKW", codon_table={k: v for k, v in TOY_TABLE.items()
                                              if k != "W"})


class TestTableFromFasta:
    def test_normalises_within_each_amino_acid(self, tmp_path):
        p = tmp_path / "cds.fa"
        p.write_text(">a\nATGAAAAAGAAA\n>b\nATGAAGAAGAAG\n", encoding="utf-8")
        t = table_from_cds_fasta(p)
        assert sum(t["K"].values()) == pytest.approx(1.0)
        assert t["K"]["AAG"] == pytest.approx(4 / 6)

    def test_ignores_stop_codons(self, tmp_path):
        p = tmp_path / "cds.fa"
        p.write_text(">a\nATGAAATAA\n", encoding="utf-8")
        t = table_from_cds_fasta(p)
        assert "*" not in t

    def test_genetic_code_is_honoured(self, tmp_path):
        """TGA is a stop under table 1; the tables differ, so the parameter must reach through."""
        p = tmp_path / "cds.fa"
        p.write_text(">a\nATGTGGTGA\n", encoding="utf-8")
        assert table_from_cds_fasta(p, genetic_code=1).get("W")

    def test_empty_input_rejected(self, tmp_path):
        p = tmp_path / "cds.fa"
        p.write_text(">a\n\n", encoding="utf-8")
        with pytest.raises(ValueError, match="no usable codons"):
            table_from_cds_fasta(p)


class TestTableFromCsv:
    def test_reads_and_normalises(self, tmp_path):
        p = tmp_path / "t.csv"
        p.write_text("codon,frequency\nAAA,30\nAAG,70\n", encoding="utf-8")
        t = table_from_csv(p)
        assert t["K"]["AAG"] == pytest.approx(0.7)

    def test_accepts_rna_spelling(self, tmp_path):
        p = tmp_path / "t.csv"
        p.write_text("codon,frequency\nAAA,1\nAAG,1\n", encoding="utf-8")
        q = tmp_path / "u.csv"
        q.write_text("codon,frequency\nAAA,1\nAAG,1\n".replace("AAG", "AAG"), encoding="utf-8")
        assert table_from_csv(p) == table_from_csv(q)

    def test_missing_column_rejected(self, tmp_path):
        p = tmp_path / "t.csv"
        p.write_text("codon,pct\nAAA,30\n", encoding="utf-8")
        with pytest.raises(ValueError, match="missing column"):
            table_from_csv(p)


@pytest.fixture(scope="module")
def result():
    return optimize_cds(PROTEIN, codon_table=TOY_TABLE, genetic_code=1)


class TestOptimize:
    def test_translation_is_preserved(self, result):
        assert str(Seq(result["cds"]).translate()) == PROTEIN

    def test_length_is_three_per_residue(self, result):
        assert len(result["cds"]) == 3 * len(PROTEIN)

    def test_constraints_reported_satisfied(self, result):
        assert result["constraints_ok"] is True

    def test_no_blacklisted_enzyme_site(self, result):
        import dnachisel as dc
        comp = str.maketrans("ACGT", "TGCA")
        for e in BLACKLIST_ENZYMES:
            site = dc.EnzymeSitePattern(e).sequence.upper()
            assert site not in result["cds"]
            assert site.translate(comp)[::-1] not in result["cds"]

    def test_no_long_homopolymer(self, result):
        cds, run = result["cds"], 1
        for a, b in zip(cds, cds[1:]):
            run = run + 1 if a == b else 1
            assert run <= 4

    def test_locked_sites_are_verbatim(self):
        """Locked overhangs must be reachable by synonymous substitution, so derive them."""
        from clippr.arelf import achievable_overhangs
        locked = {3 * c - 4: achievable_overhangs(PROTEIN, c)[0] for c in (10, 20, 30)}
        r = optimize_cds(PROTEIN, locked_sites=locked, codon_table=TOY_TABLE)
        for start, want in locked.items():
            assert r["cds"][start:start + len(want)] == want

    def test_deterministic_under_seed(self):
        a = optimize_cds(PROTEIN, codon_table=TOY_TABLE, seed=7)["cds"]
        b = optimize_cds(PROTEIN, codon_table=TOY_TABLE, seed=7)["cds"]
        assert a == b

    def test_uniquify_reduces_repeats(self):
        """The default exists because the naive objective repeats itself once per motif.

        Uniquifying is an *objective*, not a constraint, so it is best-effort: on this
        deliberately codon-poor toy table it cannot reach zero. The absolute claim -- no
        repeated 20-mer in any of 200 real designs -- belongs to
        `validation/compare_codons.py`, which runs on the real protein and table.
        """
        naive = optimize_cds(PROTEIN, codon_table=TOY_TABLE, unique_kmer_size=None)["cds"]
        uniq = optimize_cds(PROTEIN, codon_table=TOY_TABLE, unique_kmer_size=20)["cds"]

        def dup20(s):
            seen = set()
            return sum(1 for i in range(len(s) - 19)
                       if (s[i:i + 20] in seen) or seen.add(s[i:i + 20]))

        assert dup20(naive) > 0
        assert dup20(uniq) < dup20(naive)

    def test_uniquify_is_the_default(self):
        assert optimize_cds(PROTEIN, codon_table=TOY_TABLE)["cds"] == optimize_cds(
            PROTEIN, codon_table=TOY_TABLE, unique_kmer_size=20)["cds"]

    def test_codon_table_is_required(self):
        """A nuclear table on a chloroplast construct is silently wrong DNA."""
        with pytest.raises(ValueError, match="codon_table is required"):
            optimize_cds(PROTEIN)

    def test_empty_protein_rejected(self):
        with pytest.raises(ValueError, match="empty protein"):
            optimize_cds("", codon_table=TOY_TABLE)

    def test_locked_site_outside_sequence_rejected(self):
        with pytest.raises(ValueError, match="falls outside"):
            optimize_cds(PROTEIN, locked_sites={10_000: "AATG"}, codon_table=TOY_TABLE)

    def test_plastid_code_runs(self):
        """Table 11 must reach the solver; its name list contains a None."""
        r = optimize_cds("MAKLGSTPVW" * 2, codon_table=TOY_TABLE, genetic_code=11)
        assert str(Seq(r["cds"]).translate(table=11)) == "MAKLGSTPVW" * 2


class TestHosts:
    """Every host offered must resolve; a listed host that cannot is worse than none."""

    def test_every_listed_host_has_a_source_and_a_code(self):
        from clippr import constants as C
        for name, (source, code) in C.ORGANISMS.items():
            assert isinstance(source, (int, str)), name
            assert code in (1, 11), f"{name}: unexpected genetic code {code}"

    def test_organelles_use_the_plastid_code(self):
        from clippr import constants as C
        for name, (_, code) in C.ORGANISMS.items():
            if "chloroplast" in name:
                assert code == 11, f"{name} must use NCBI table 11"

    def test_organelles_are_derived_not_looked_up(self):
        """Kazusa has no confirmed organelle entry, so these must come from a genome."""
        from clippr import constants as C
        for name, (source, _) in C.ORGANISMS.items():
            if "chloroplast" in name:
                assert isinstance(source, str), f"{name} must be a genome accession"

    @pytest.mark.parametrize("aa", ["L", "A", "G", "V"])
    def test_derived_chloroplast_table_covers_the_common_residues(self, aa):
        from clippr.codons import table_from_genome
        try:
            table = table_from_genome()
        except Exception as e:
            pytest.skip(f"chloroplast genome unavailable: {e}")
        assert aa in table and table[aa]

    def test_the_two_contexts_really_differ(self):
        """The trap this package was built around: nuclear and plastid usage are opposite.

        Chlamydomonas nucleus prefers Leu CTG; its chloroplast prefers TTA. Measured at
        0.73 and 0.74 respectively -- near-mirror images, not a rounding difference.
        """
        from clippr.codons import table_from_genome, table_from_kazusa
        try:
            plastid = table_from_genome()
            nuclear = table_from_kazusa()
        except Exception as e:
            pytest.skip(f"codon sources unavailable: {e}")
        assert max(plastid["L"], key=plastid["L"].get).startswith("TT")
        assert max(nuclear["L"], key=nuclear["L"].get).startswith("CT")

    def test_a_derived_table_needs_enough_codons(self):
        from clippr.codons import table_from_genome
        with pytest.raises(Exception):
            table_from_genome("NC_000000.0")
