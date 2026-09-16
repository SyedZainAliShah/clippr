"""Unit tests for the GRASP kit realisation route.

The load-bearing properties are physical rather than cosmetic: the chosen modules must chain
by their overhangs, no sub-assembly may repeat an overhang (Golden Gate would misligate), and
decoding the chosen modules must give back the target that asked for them. A pick list that
looks plausible and does not assemble is the failure mode worth guarding.
"""
from __future__ import annotations

import pytest

from clippr.parts import (
    BUILDABLE_LENGTHS,
    LINKERS,
    LOOP_BLOCKS,
    PartsPlan,
    base_to_code,
    block_chain,
    code_to_base,
    inventory,
    report,
    select,
)

NINE, FOURTEEN, NINETEEN = "AAAAUGUGG", "GCUAAAGACUUGCA", "AAAGCGGCACUUGUGAAGU"


class TestInventory:
    def test_holds_the_deposited_kit(self):
        assert len(inventory()) == 42

    def test_every_module_is_in_the_published_vector(self):
        assert {m.vector for m in inventory()} == {"pAGM1311"}

    def test_specificity_residues_are_only_the_published_options(self):
        """A 5th residue of D would be a plasmid the kit does not contain."""
        assert {m.fifth for m in inventory() if m.fifth} == {"N", "T"}
        assert {m.last for m in inventory() if m.last} == {"D", "N"}

    def test_plate_positions_are_unique(self):
        plates = [m.plate for m in inventory()]
        assert len(plates) == len(set(plates))

    def test_the_code_matches_the_published_mapping(self):
        assert code_to_base() == {"TN": "A", "NN": "C", "TD": "G", "ND": "U"}

    def test_every_base_has_a_code(self):
        assert {base_to_code(b) for b in "ACGU"} == {"TN", "NN", "TD", "ND"}

    def test_an_unreadable_base_is_rejected(self):
        with pytest.raises(ValueError, match="no PPR code reads"):
            base_to_code("X")


class TestBlockChain:
    @pytest.mark.parametrize("n,modules,subs", [(9, 10, 2), (14, 15, 3), (19, 20, 4)])
    def test_the_three_architectures_resolve(self, n, modules, subs):
        chain = block_chain(n)
        assert len(chain) == modules
        assert chain[0] == "1A" and chain[-1] == "2E"

    def test_a_length_the_kit_cannot_build_is_refused(self):
        with pytest.raises(ValueError, match="9, 14 or 19"):
            block_chain(12)

    def test_the_buildable_set_is_exactly_the_three_architectures(self):
        """Derived from the linker count, so it cannot drift from the inventory."""
        assert BUILDABLE_LENGTHS == (9, 14, 19)

    @pytest.mark.parametrize("n", [n for n in range(1, 40) if n not in (9, 14, 19)])
    def test_no_other_length_is_accepted(self, n):
        """Divisibility alone let n=4 through: 5 modules is a whole run but not a PPR.

        Parametrising every length up to 39 rather than a hand-picked few is the difference
        between catching that and missing it, which is what happened.
        """
        with pytest.raises(ValueError):
            block_chain(n)

    def test_beyond_the_kit_says_which_resource_ran_out(self):
        """24 bases would need a fourth internal linker; the kit contains three."""
        with pytest.raises(ValueError, match="internal linkers"):
            block_chain(24)

    def test_each_internal_join_uses_a_distinct_linker(self):
        chain = block_chain(19)
        used = [b for b in chain if b in {a for a, _ in LINKERS}]
        assert len(used) == len(set(used)) == 3

    def test_loop_blocks_repeat_once_per_sub_assembly(self):
        chain = block_chain(19)
        for block in LOOP_BLOCKS:
            assert chain.count(block) == 4


class TestSelection:
    @pytest.mark.parametrize("target", [NINE, FOURTEEN, NINETEEN])
    def test_every_architecture_is_realisable(self, target):
        plan = select(target)
        assert plan.available, plan.reason
        assert len(plan.modules) == len(target) + 1

    @pytest.mark.parametrize("target", [NINE, FOURTEEN, NINETEEN])
    def test_overhangs_chain_end_to_end(self, target):
        """Each module's 3' overhang must be the next one's 5', or it will not assemble."""
        mods = select(target).modules
        for a, b in zip(mods, mods[1:]):
            assert a.three_overhang == b.five_overhang, (
                f"{a.plasmid_id} ends {a.three_overhang} but {b.plasmid_id} starts "
                f"{b.five_overhang}")

    @pytest.mark.parametrize("target", [NINE, FOURTEEN, NINETEEN])
    def test_no_sub_assembly_repeats_an_overhang(self, target):
        """Golden Gate needs unique overhangs in one reaction; a repeat misligates.

        This is why the kit splits a long PPR into sub-assemblies at all, so it is the
        property most worth pinning.
        """
        for i, sub in enumerate(select(target).sub_assemblies, 1):
            overhangs = [sub[0].five_overhang] + [m.three_overhang for m in sub]
            assert len(overhangs) == len(set(overhangs)), \
                f"sub-assembly {i} of {target} repeats an overhang: {overhangs}"

    @pytest.mark.parametrize("target", [NINE, FOURTEEN, NINETEEN])
    def test_the_modules_decode_back_to_the_target(self, target):
        """The round trip. If this holds, the pick list really encodes what was asked for."""
        mods = select(target).modules
        decoded = "".join(code_to_base()[f"{a.fifth}{b.last}"]
                          for a, b in zip(mods, mods[1:]))
        assert decoded == target

    def test_sub_assembly_count_follows_the_architecture(self):
        assert len(select(NINE).sub_assemblies) == 2
        assert len(select(FOURTEEN).sub_assemblies) == 3
        assert len(select(NINETEEN).sub_assemblies) == 4

    def test_every_module_is_from_the_inventory(self):
        ids = {m.plasmid_id for m in inventory()}
        assert all(m.plasmid_id in ids for m in select(NINE).modules)


class TestUnavailable:
    """A target the kit cannot build is an answer, not an error."""

    def test_a_bad_length_returns_a_reason_rather_than_raising(self):
        plan = select("ACGUACGUACGU")
        assert not plan.available
        assert "9, 14 or 19" in plan.reason

    def test_a_target_beyond_the_kit_is_refused_with_its_reason(self):
        plan = select("ACGU" * 6)          # 24 bases
        assert not plan.available
        assert "linker" in plan.reason

    def test_a_non_rna_string_is_refused(self):
        assert not select("hello world").available

    def test_a_four_base_target_is_refused_rather_than_silently_built(self):
        """The regression: 4+1 is a whole run, so the old guard produced a confident plan.

        A wet lab following that output would have pulled five plasmids for an architecture
        the kit cannot assemble.
        """
        plan = select("ACGU")
        assert not plan.available
        assert plan.reason and "is not one of them" in plan.reason

    def test_an_empty_target_is_refused(self):
        assert not select("").available

    def test_the_report_says_to_use_the_other_route(self):
        text = report(select("ACGUACGUACGU"))
        assert "unavailable" in text and "de novo" in text


class TestReport:
    def test_lists_plate_positions_for_the_bench(self):
        text = report(select(NINE))
        for plate in select(NINE).plates:
            assert plate in text

    def test_states_the_scope_limit(self):
        """The route selects PPR modules, not a complete construct, and must say so."""
        text = report(select(NINE))
        assert "not a complete construct" in text
        assert "acceptor" in text

    def test_credits_the_published_source(self):
        assert "Table S1" in report(select(NINE))

    def test_says_no_synthesis_is_needed(self):
        assert "no new PPR DNA needs synthesising" in report(select(NINE))


class TestPlanShape:
    def test_an_unavailable_plan_carries_no_modules(self):
        plan = PartsPlan(target="ACGU", reason="too short")
        assert not plan.available and plan.modules == () and plan.plates == ()

    def test_a_plan_with_no_modules_is_not_available_even_without_a_reason(self):
        """`reason is None` alone made this report as available and then crash report()."""
        assert not PartsPlan(target="ACGU").available

    def test_report_on_such_a_plan_explains_rather_than_crashing(self):
        text = report(PartsPlan(target="ACGU"))
        assert "unavailable" in text
        assert "no modules were selected" in text
