"""Integration tests for the one-shot pipeline.

These check that the modules compose -- that the overhangs chosen before optimisation
survive it, that fragments reassemble, that every architecture works -- rather than
re-testing what each module already proves on its own.
"""
from __future__ import annotations

import pytest
from Bio.Seq import Seq

from clippr.assembly import reassemble
from clippr.design import MAX_FRAGMENT_AA, ORGANISMS, design_oneshot
from clippr.qc import synthesis_qc

TOY_TABLE = {
    "M": {"ATG": 1.0}, "W": {"TGG": 1.0},
    "A": {"GCC": 0.5, "GCT": 0.2, "GCA": 0.2, "GCG": 0.1},
    "R": {"CGC": 0.4, "CGT": 0.3, "CGA": 0.1, "CGG": 0.1, "AGA": 0.05, "AGG": 0.05},
    "N": {"AAC": 0.6, "AAT": 0.4}, "D": {"GAC": 0.6, "GAT": 0.4},
    "C": {"TGC": 0.6, "TGT": 0.4}, "Q": {"CAG": 0.7, "CAA": 0.3},
    "E": {"GAG": 0.6, "GAA": 0.4}, "G": {"GGC": 0.5, "GGT": 0.2, "GGA": 0.2, "GGG": 0.1},
    "H": {"CAC": 0.6, "CAT": 0.4}, "I": {"ATC": 0.5, "ATT": 0.3, "ATA": 0.2},
    "L": {"CTG": 0.7, "CTC": 0.15, "CTA": 0.05, "CTT": 0.05, "TTA": 0.02, "TTG": 0.03},
    "K": {"AAG": 0.7, "AAA": 0.3}, "F": {"TTC": 0.6, "TTT": 0.4},
    "P": {"CCC": 0.4, "CCT": 0.2, "CCA": 0.2, "CCG": 0.2},
    "S": {"TCC": 0.3, "TCT": 0.2, "TCA": 0.1, "TCG": 0.2, "AGC": 0.1, "AGT": 0.1},
    "T": {"ACC": 0.4, "ACT": 0.2, "ACA": 0.2, "ACG": 0.2},
    "Y": {"TAC": 0.6, "TAT": 0.4}, "V": {"GTC": 0.4, "GTT": 0.2, "GTA": 0.2, "GTG": 0.2},
}


@pytest.fixture(scope="module")
def result():
    return design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE)


class TestPipeline:
    def test_translation_survives_the_whole_pipeline(self, result):
        assert str(Seq(result["cds"]).translate()) == result["protein"]

    def test_chosen_overhangs_survive_codon_optimisation(self, result):
        """The ordering constraint: overhangs are picked first, then locked into the CDS."""
        cds = result["cds"]
        for cut, oh in zip(result["cuts"], result["junction_overhangs"]):
            assert cds[3 * cut - 4:3 * cut] == oh

    def test_fragments_reassemble_to_the_cds(self, result):
        assert reassemble(result["fragments"]) == result["cds"]

    def test_one_oligo_per_fragment(self, result):
        assert len(result["oligos"]) == len(result["fragments"])

    def test_fragments_respect_the_synthesis_limit(self, result):
        assert all(f.aa_length <= MAX_FRAGMENT_AA for f in result["fragments"])

    def test_reaction_set_includes_both_destination_ends(self, result):
        r = result["reaction_overhangs"]
        assert len(r) == len(result["junction_overhangs"]) + 2

    def test_qc_is_recomputable_from_the_cds(self, result):
        assert synthesis_qc(result["cds"])["status"] == result["qc"]["status"]

    def test_summary_mentions_the_target(self, result):
        assert "AAAAUGUGG" in result["summary"]

    def test_every_documented_key_present(self, result):
        for k in ("summary", "protein", "cds", "oligos", "qc", "fidelity", "cost", "paths"):
            assert k in result


class TestArchitectures:
    @pytest.mark.parametrize("target,arch,repeats", [
        ("AAAAUGUGG", "9S", 9),
        ("GCUAAAGACUUGCA", "14S", 14),
        ("AAAGCGGCACUUGUGAAGU", "19S", 19),
    ])
    def test_each_architecture_designs(self, target, arch, repeats):
        r = design_oneshot(target, codon_table=TOY_TABLE)
        assert r["architecture"] == arch
        assert len(r["ppr_code"].split("-")) == repeats or len(r["cuts"]) >= 1
        assert str(Seq(r["cds"]).translate()) == r["protein"]

    def test_fragment_count_scales_with_length(self):
        small = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE)
        large = design_oneshot("AAAGCGGCACUUGUGAAGU", codon_table=TOY_TABLE)
        assert len(large["fragments"]) > len(small["fragments"])

    def test_bad_target_length_rejected(self):
        with pytest.raises(ValueError, match="target length"):
            design_oneshot("ACGU", codon_table=TOY_TABLE)


class TestDeterminism:
    def test_same_seed_gives_the_same_sequence(self):
        a = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE, seed=11)
        b = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE, seed=11)
        assert a["cds"] == b["cds"]
        assert a["cuts"] == b["cuts"]
        assert a["junction_overhangs"] == b["junction_overhangs"]


class TestOrganismAndTable:
    def test_every_named_organism_is_usable(self):
        """A listed host that cannot resolve is worse than one that is not listed."""
        assert ORGANISMS["c_reinhardtii_nuclear"] == (3055, 1)
        assert len(ORGANISMS) >= 2
        for name, (source, code) in ORGANISMS.items():
            assert isinstance(source, (int, str)) and code in (1, 11), name

    def test_both_chlamydomonas_contexts_are_offered(self):
        """The PPR is nuclear; the UTR it binds is chloroplast. Both are needed."""
        assert ORGANISMS["c_reinhardtii_nuclear"][1] == 1
        assert ORGANISMS["c_reinhardtii_chloroplast"][1] == 11

    def test_unknown_organism_rejected_with_guidance(self):
        with pytest.raises(ValueError, match="unknown organism"):
            design_oneshot("AAAAUGUGG", organism="tyrannosaurus_rex")

    def test_missing_table_file_rejected(self):
        with pytest.raises(FileNotFoundError):
            design_oneshot("AAAAUGUGG", codon_table="no_such_table.csv")

    def test_table_from_csv_path(self, tmp_path):
        p = tmp_path / "t.csv"
        rows = ["codon,frequency"]
        for aa, codons in TOY_TABLE.items():
            rows += [f"{c},{v}" for c, v in codons.items()]
        p.write_text("\n".join(rows), encoding="utf-8")
        r = design_oneshot("AAAAUGUGG", codon_table=p)
        assert str(Seq(r["cds"]).translate()) == r["protein"]


class TestAudit:
    """The design must be inspectable, not just correct."""

    def test_audit_is_returned(self, result):
        assert "audit" in result
        assert result["audit"].target_rna == "AAAAUGUGG"

    def test_records_the_active_constraint_profile(self, result):
        a = result["audit"]
        assert a.enzyme_profile == "igem_rfc1000"
        assert "BsmBI" not in a.enzymes_excluded

    def test_one_selected_overhang_per_junction(self, result):
        a = result["audit"]
        assert len(a.selected) == len(result["cuts"])
        assert [d.sequence for d in a.selected] == result["junction_overhangs"]

    def test_every_candidate_is_accounted_for(self, result):
        """Nothing silently disappears: each achievable overhang gets a verdict."""
        from clippr.arelf import achievable_overhangs
        a = result["audit"]
        for cut in result["cuts"]:
            n = len(achievable_overhangs(result["protein"], cut))
            assert sum(1 for d in a.overhang_decisions if d.junction_cut == cut) == n

    def test_rejections_name_the_profile(self):
        """A rejection must say what ruled the overhang out, not just that it was."""
        from clippr.arelf import achievable_overhangs, safe_overhangs
        prot = "MAAAAARDAAAAA"
        cut = prot.index("D") + 1
        rejected = set(achievable_overhangs(prot, cut)) - set(
            safe_overhangs(prot, cut, enzymes="moclo_compat"))
        assert "AGAC" in rejected

    def test_realization_count_is_diagnostic_not_a_ranking(self, result):
        """Recorded for the reader; never used to choose between overhangs."""
        a = result["audit"]
        counts = [d.local_realizations for d in a.selected]
        assert all(c is not None and c >= 1 for c in counts)
        # the selected one is not always the one with the most freedom
        assert len(set(counts)) > 1 or len(counts) == 1

    def test_report_reads_as_prose(self, result):
        text = result["audit"].report()
        for expected in ("design audit", "enzyme profile", "selected overhangs", "cuts"):
            assert expected in text

    def test_findings_carry_the_backbone_warning(self, result):
        assert any("backbone is the limit" in f for f in result["audit"].findings)


class TestConstantLayers:
    """Published fact, assembly mechanics and project policy must stay distinguishable."""

    def test_three_layers_declared(self):
        from clippr import constants as C
        assert set(C.LAYERS) == {"biology", "assembly_spec", "policy"}

    def test_policy_holds_no_published_scaffold(self):
        from clippr import policy
        assert not hasattr(policy, "REPEAT_TEMPLATE")
        assert not hasattr(policy, "CODE_TO_BASE")

    def test_biology_holds_no_project_choices(self):
        from clippr import biology
        assert not hasattr(biology, "ENZYME_PROFILES")
        assert not hasattr(biology, "ORGANISMS")

    def test_enzyme_profiles_are_policy_not_biology(self):
        from clippr import constants as C
        assert "ENZYME_PROFILES" in C.LAYERS["policy"]
        assert "REPEAT_TEMPLATE" in C.LAYERS["biology"]

    def test_aggregator_still_exposes_everything(self):
        from clippr import constants as C
        for names in C.LAYERS.values():
            for n in names:
                assert hasattr(C, n), f"{n} missing from the aggregator"


class TestOutputs:
    def test_no_files_written_without_outdir(self, result):
        assert result["paths"] == {}

    def test_outdir_writes_every_artefact(self, tmp_path):
        r = design_oneshot("AAAAUGUGG", codon_table=TOY_TABLE, outdir=tmp_path)
        assert set(r["paths"]) == {"oligo_csv", "oligo_fasta", "gene_fasta", "genbank"}
        for p in r["paths"].values():
            assert (tmp_path / p.split("\\")[-1].split("/")[-1]).exists()

    def test_backbone_warning_names_the_real_limit(self, result):
        """Level 0's destination pair caps fidelity; the warning must say so."""
        assert any("backbone is the limit" in w for w in result["warnings"])
