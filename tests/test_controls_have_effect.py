"""Every parameter offered must actually change the outcome.

A control that changes nothing is worse than no control: it implies a decision was made
when none was. This project has already shipped one -- `enzyme_profile` sat in the notebook
form for a while without being passed through, so switching it did nothing at all.

Each test varies one parameter and asserts the result differs in the way that parameter
should cause. Where a parameter can legitimately be a no-op, the test is built so it cannot
be: the EagI case below picks a site that *is* present, because blacklisting a site the
sequence never had proves nothing.
"""
from __future__ import annotations

import pytest
from Bio import Restriction

from clippr import design_oneshot
from test_design import TOY_TABLE

BASE = dict(target_rna="AAAAUGUGG", seed=42, check_offtarget=False,
            codon_table=TOY_TABLE)


def rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.fixture(scope="module")
def ref():
    return design_oneshot(**BASE)


class TestControlsHaveEffect:
    def test_enzyme_changes_the_ordered_oligo(self, ref):
        alt = design_oneshot(**{**BASE, "enzyme": "BbsI"})
        a = ref["oligos"].iloc[0]["oligo_sequence_5to3"]
        b = alt["oligos"].iloc[0]["oligo_sequence_5to3"]
        assert a != b
        assert str(Restriction.BsaI.site) in a
        assert str(Restriction.BbsI.site) in b

    def test_matrix_changes_predicted_fidelity(self, ref):
        alt = design_oneshot(**{**BASE, "matrix": "BbsI-HF"})
        assert abs(alt["fidelity"] - ref["fidelity"]) > 1e-9

    def test_enzyme_profile_changes_what_is_excluded(self, ref):
        alt = design_oneshot(**{**BASE, "enzyme_profile": "assembly"})
        assert alt["audit"].enzymes_excluded != ref["audit"].enzymes_excluded

    def test_extra_blacklist_removes_a_site_that_was_there(self, ref):
        """Uses a site the default sequence actually contains.

        Blacklisting an absent site is a legitimate no-op and proves nothing, which is
        how an earlier version of this test wrongly reported the parameter as ignored.
        """
        site = str(Restriction.EagI.site)
        before = ref["cds"].count(site) + ref["cds"].count(rc(site))
        assert before > 0, "pick a site the unconstrained sequence really contains"

        alt = design_oneshot(**{**BASE, "extra_blacklist": "EagI"})
        after = alt["cds"].count(site) + alt["cds"].count(rc(site))
        assert after == 0
        assert "EagI" in alt["audit"].enzymes_excluded

    def test_destination_changes_the_reaction(self, ref):
        alt = design_oneshot(**{**BASE, "destination": ("GGAG", "CGCT")})
        assert alt["reaction_overhangs"][0] == "GGAG"
        assert alt["reaction_overhangs"] != ref["reaction_overhangs"]

    def test_n_fragments_is_obeyed(self, ref):
        alt = design_oneshot(**{**BASE, "n_fragments": 5})
        assert len(alt["oligos"]) == 5 != len(ref["oligos"])

    def test_seed_reproduces_exactly(self, ref):
        assert design_oneshot(**BASE)["cds"] == ref["cds"]

    def test_check_offtarget_toggles_the_report(self, ref):
        assert ref["offtarget"] is None
        on = design_oneshot(**{**BASE, "check_offtarget": True})
        assert on["offtarget"] is not None

    def test_codon_table_changes_the_sequence(self, ref):
        sparse = {aa: {next(iter(c)): 1.0} for aa, c in TOY_TABLE.items()}
        alt = design_oneshot(**{**BASE, "codon_table": sparse})
        assert alt["cds"] != ref["cds"]

    def test_genetic_code_reaches_the_audit(self, ref):
        assert ref["audit"].genetic_code == 1
        alt = design_oneshot(**{**BASE, "genetic_code": 11})
        assert alt["audit"].genetic_code == 11


class TestNotebookControlsAreWired:
    """Every #@param in the notebook must be passed to something, not merely displayed."""

    def test_every_form_control_is_used(self):
        import json
        from pathlib import Path

        nb_path = Path(__file__).resolve().parents[1] / "notebooks" / "CLIPPR_designer.ipynb"
        if not nb_path.exists():
            pytest.skip("notebook not present")
        nb = json.loads(nb_path.read_text(encoding="utf-8"))
        source = "".join("".join(c["source"]) for c in nb["cells"]
                         if c["cell_type"] == "code")

        declared = {line.split("=")[0].strip()
                    for line in source.splitlines() if "#@param" in line}
        assert declared, "no form controls found"

        for name in declared:
            # a control is wired if it is read somewhere other than its own declaration
            uses = [l for l in source.splitlines()
                    if name in l and "#@param" not in l]
            assert uses, f"{name} is declared in a form but never used"
