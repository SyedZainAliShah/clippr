"""The workflow boundaries a synthesis policy has to survive.

An independent review found the profile reached every *signature* and few of the *consumers*.
Each test below reproduces one of its findings and fails on the code as it stood:

    F1  `order_items_for(profile=...)` accepted the argument and exported under the default,
        so a broad request rejected compliant modules and a tighter request reported success
        over modules that breach the band it was handed.
    F2  a front recorded no policy, selection rechecked under the default, and a candidate the
        verifier rejected still returned `ok=True` with `objectives=None`.
    F3  `solver_bounds()` returned identical arguments for `hard` and `target`, so a declared
        target was a hard constraint and nothing was ever soft.
    F4  the collection solver received a profile and none of its bounds.
    F5  `frozen=True` left `enforcement` mutable and `narrowed()` shared it, so editing a
        derived policy rewrote its parent.
    F6  completeness was `bool(extra) and bool(evidence)`, so one overhang plus prose saying
        the list was incomplete returned `participants_established=True`.

The lesson they share, and the reason these live at the boundaries rather than in the unit
suites: a parameter in a signature is not evidence that the consumer uses it. Every one of
these defects sat under a green test suite.
"""
from __future__ import annotations

import json
import random
from pathlib import Path
from unittest.mock import patch

import pytest

from clippr import assembly_spec, codons
from clippr import inventories as inv
from clippr import objectives
from clippr import joint_search as js
from clippr import library_search as ls
from clippr import substrates as sub
from clippr import synthesis_profile as sp
from clippr import workflow as w

ROOT = Path(__file__).resolve().parents[1]
CODON_TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"
RECODED = ROOT / "work" / "phaseb" / "inventory_recoded.json"
TARGETS = ["AAAAUGUGG", "UUACACGUGCGUAC", "CUAUCACAUCACAUAAGCG"]

needs_inventory = pytest.mark.skipif(
    not RECODED.is_file(), reason="a recoded inventory is needed for the boundary tests")

#: Deliberately not one of the three registered profiles. A policy assembled by a caller has
#: to survive the same round trips as a shipped one, and only a custom one proves that.
TIGHTER = sp.STRICT_LEGACY.narrowed(name="tighter-than-shipped", local_gc=(0.40, 0.60))


@pytest.fixture(scope="module")
def raw_table():
    return json.loads(CODON_TABLE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def library():
    return inv.load(RECODED)


# --------------------------------------------------------------------------- F3: target != hard

class TestATargetSteersTheSolverRatherThanGatingIt:
    def test_the_two_enforcements_send_different_arguments(self):
        advisory = sp.STRICT_LEGACY.narrowed(
            name="advisory", enforcement={"local_gc": sp.TARGET, "homopolymer": sp.TARGET,
                                          "forbidden_sites": sp.HARD})
        assert sp.STRICT_LEGACY.solver_bounds() != advisory.solver_bounds()
        assert advisory.solver_bounds()["soft_rules"] == ("homopolymer", "local_gc")
        assert sp.STRICT_LEGACY.solver_bounds()["soft_rules"] == ()

    def test_an_unsatisfiable_target_still_returns_a_sequence(self):
        """The review's counterexample. WW can only be TGGTGG, so a cap of 1 is impossible.

        Declared as a target the solver must deliver the sequence and carry the breach.
        Anything else means target is a synonym for hard.
        """
        advisory = sp.BROAD_EXPERIMENTAL.narrowed(
            name="advisory-cap-one", max_homopolymer=1,
            enforcement={"homopolymer": sp.TARGET, "local_gc": sp.TARGET,
                         "forbidden_sites": sp.HARD})
        got = codons.optimize_cds("WW", codon_table={"W": {"TGG": 1.0}},
                                  unique_kmer_size=None, **advisory.solver_bounds())
        assert got["cds"] == "TGGTGG"
        assert got["constraints_ok"] is True
        assert "homopolymer" in got["soft_rules"]

    def test_a_strict_target_reaches_the_same_optimum_as_a_wide_limit(self, raw_table):
        """The finding the four-arm benchmark produced, pinned on one module.

        Enforcing the strict band as a *target* reaches the same codon adaptation as widening
        the band to 0.15-0.85, because the solver is free to trade the band against adaptation
        instead of being walled off. Over the full kit that is worth +0.201091 mean CAI on 40
        of 42 modules. If this ever stops holding, the soft path has collapsed back to hard and
        `docs/synthesis_policy.md` section B is wrong again.
        """
        strict_soft = sp.STRICT_LEGACY.narrowed(
            name="strict-as-target",
            enforcement={"local_gc": sp.TARGET, "homopolymer": sp.TARGET,
                         "forbidden_sites": sp.HARD})
        protein = "MAKEFVRLLTQDGHIPWNKACYSEM" * 2
        table = codons.complete_table(raw_table, 1)
        soft = codons.optimize_cds(protein, codon_table=table, seed=7,
                                   unique_kmer_size=None, **strict_soft.solver_bounds())
        wide = codons.optimize_cds(protein, codon_table=table, seed=7, unique_kmer_size=None,
                                   **sp.BROAD_EXPERIMENTAL.solver_bounds())
        assert soft["constraints_ok"] and wide["constraints_ok"]
        cai = (lambda cds: objectives.codon_adaptation(cds, table)["cai"])
        assert cai(soft["cds"]) == pytest.approx(cai(wide["cds"]), abs=1e-9)

        # And the strict band as a *hard* limit does not reach it -- the control, without
        # which the assertion above would also pass if every arm simply scored the same.
        hard = codons.optimize_cds(protein, codon_table=table, seed=7, unique_kmer_size=None,
                                   **sp.STRICT_LEGACY.solver_bounds())
        assert cai(hard["cds"]) < cai(wide["cds"])

    def test_the_same_cap_declared_hard_does_refuse(self):
        """The control. Without it the test above would prove nothing about enforcement."""
        strict = sp.BROAD_EXPERIMENTAL.narrowed(name="hard-cap-one", max_homopolymer=1)
        got = codons.optimize_cds("WW", codon_table={"W": {"TGG": 1.0}},
                                  unique_kmer_size=None, **strict.solver_bounds())
        assert got["constraints_ok"] is False

    def test_an_enzyme_site_can_never_be_declared_a_preference(self):
        """The solver always refused to soften a site; the profile used to advertise it anyway.

        That mismatch is worse than either behaviour alone: `solver_bounds()` reported
        `forbidden_sites` among the soft rules and the recorded binding carried it, so an
        artefact could claim sites were advisory while the solver had treated them as
        absolute. Constructing such a policy is a caller error and now raises.
        """
        with pytest.raises(ValueError, match="not a preference"):
            sp.STRICT_LEGACY.narrowed(
                name="tries-to-soften-sites",
                enforcement={"forbidden_sites": sp.TARGET, "local_gc": sp.TARGET})


# ------------------------------------------------------------------- F5: a policy is immutable

class TestAPolicyCannotBeEditedInPlace:
    def test_the_enforcement_map_refuses_mutation(self):
        with pytest.raises(TypeError):
            sp.STRICT_LEGACY.enforcement["local_gc"] = sp.TARGET   # type: ignore[index]

    def test_a_derived_policy_does_not_share_its_parents_map(self):
        parent = sp.STRICT_LEGACY.narrowed(name="parent")
        child = parent.narrowed(name="child", enforcement={"local_gc": sp.TARGET})
        assert parent.enforcement is not child.enforcement
        assert parent.how("local_gc") == sp.HARD
        assert child.how("local_gc") == sp.TARGET

    def test_every_shipped_policy_survives_its_recorded_form(self):
        for profile in sp.PROFILES.values():
            assert sp.from_dict(profile.as_dict()).version == profile.version

    def test_a_custom_policy_survives_its_recorded_form(self):
        """The reload that matters: selection rebuilds from a record, not from the registry."""
        back = sp.from_dict(TIGHTER.as_dict())
        assert back.version == TIGHTER.version
        assert back.local_gc == (0.40, 0.60)
        assert back.solver_bounds() == TIGHTER.solver_bounds()


# --------------------------------------------------------------- F6: completeness is a contract

class TestCompletenessIsAContractNotACitation:
    def test_prose_and_one_overhang_do_not_establish_completeness(self):
        """The review's counterexample, verbatim."""
        got = assembly_spec.level1_participants(
            ["AATG", "CTTC", "TTCG"], roles={"linker": ["TGTG"]},
            evidence="Only the linker start is supplied; DYW and bridging parts are missing.")
        assert got["participants_established"] is False
        assert set(got["roles_missing"]) >= {"final_cassette", "editing_domain", "bridge"}

    def test_evidence_alone_is_not_enough(self):
        got = assembly_spec.level1_participants(["AATG", "CTTC", "TTCG"],
                                                evidence="Pryor 2020, Table 2")
        assert got["participants_established"] is False
        assert "no participant supplied" in got["why_incomplete"]

    #: A genuine chain: every internal end shared by two parts, AATG and CGCT terminal. It has
    #: to be genuine, or the citation is never the deciding factor -- an earlier version of the
    #: test below used a role map whose ends did not join, so completeness failed on geometry
    #: and the test stayed green when the citation requirement was removed entirely.
    CHAIN = {"linker": ["CTTC", "GTGA"], "editing_domain": ["GTGA", "CACG"],
             "bridge": ["CACG", "TTCG"], "final_cassette": ["TTCG", "CGCT"]}
    BLOCKS = ["AATG", "CTTC"]

    def test_an_otherwise_complete_set_is_refused_for_the_citation_alone(self):
        got = assembly_spec.level1_participants(self.BLOCKS, roles=self.CHAIN)
        assert got["participants_established"] is False
        assert got["why_incomplete"] == (
            "no evidence supplied naming where these participants came from")

    def test_the_same_set_with_a_citation_is_complete(self):
        """The contract has to be satisfiable, or refusing everything would pass every test."""
        got = assembly_spec.level1_participants(self.BLOCKS, roles=self.CHAIN,
                                                evidence="Pryor 2020, Table 2")
        assert got["participants_established"] is True
        assert got["why_incomplete"] is None
        assert got["roles_missing"] == []


# ---------------------------------------------------------- F1: export honours what it is given

@needs_inventory
class TestExportHonoursTheRequestedPolicy:
    def test_two_policies_over_one_inventory_give_two_verdicts(self, tmp_path):
        """The pairing is the test, and it has to be a pair.

        Either call on its own passes even when the profile is ignored entirely: the shipped
        inventory happens to satisfy the default, so a permissive export succeeds whether or
        not the policy arrived. Only the disagreement between two calls over the same
        sequences shows that the argument reached the substrate builder.

        This is also why the tighter band exists as a fixture rather than the broad profile
        being used alone -- measured on both saved inventories, strict and broad give
        identical export verdicts (0 failures each), so neither distinguishes the other.
        """
        broad = w.order_items_for(RECODED, tmp_path / "broad", profile=sp.BROAD_EXPERIMENTAL)
        tight = w.order_items_for(RECODED, tmp_path / "tight", profile=TIGHTER)
        assert broad.ok is True and broad.failures == []
        assert tight.ok is False and tight.failures

    def test_a_tighter_policy_fails_on_its_own_violations(self, tmp_path, library):
        """Previously ok=True over every module that breaches the requested band."""
        expected = {module_id for module_id, record in library.modules.items()
                    if sub.substrate_problems(record.dna, record.block, TIGHTER)}
        assert expected, "this fixture no longer violates the tighter band"
        got = w.order_items_for(RECODED, tmp_path / "tight", profile=TIGHTER)
        assert got.ok is False
        assert {f["module"] for f in got.failures} == expected

    def test_the_same_module_is_judged_differently_under_each(self, library):
        """The mechanism underneath: a substrate carries the policy it was built with."""
        breaching = next(module_id for module_id, record in sorted(library.modules.items())
                         if sub.substrate_problems(record.dna, record.block, TIGHTER))
        record = library.modules[breaching]
        assert sub.build(record.dna, record.block, profile=TIGHTER).synthesis_problems
        assert not sub.build(record.dna, record.block,
                             profile=sp.BROAD_EXPERIMENTAL).synthesis_problems


# ------------------------------------------------------ F4: the collection solver gets the bounds

@needs_inventory
class TestTheCollectionSolverReceivesTheRequestedBounds:
    def test_a_proposal_is_solved_under_the_profile_not_the_default(self, raw_table, library):
        module_id = sorted(library.modules)[0]
        seen: list = []
        original = codons.optimize_cds

        def observe(*args, **kwargs):
            seen.append(kwargs.get("gc_bounds", "NOT PASSED"))
            return original(*args, **kwargs)

        with patch.object(codons, "optimize_cds", side_effect=observe):
            ls._propose(module_id, library.modules[module_id],
                        codons.complete_table(raw_table, 1), random.Random(42), 1,
                        sp.BROAD_EXPERIMENTAL)
        assert seen, "the solver was never called"
        assert seen[0] == sp.BROAD_EXPERIMENTAL.local_gc, seen


# ----------------------------------------------- F2: the front carries its policy into selection

@needs_inventory
class TestSelectionKeepsThePolicyAndVerifiesBeforePublishing:
    @pytest.fixture(scope="class")
    @classmethod
    def front(cls, tmp_path_factory, raw_table):
        """One front measured under a deliberately non-default policy, shared by the class."""
        out = tmp_path_factory.mktemp("front")
        return w.explore_interfaces(RECODED, TARGETS, raw_table, out,
                                    profile=sp.BROAD_EXPERIMENTAL,
                                    max_evaluations=60, wall_seconds=900)

    def test_a_front_records_the_policy_it_was_measured_under(self, front):
        payload = json.loads(Path(front.artefacts["front"]).read_text(encoding="utf-8"))
        recorded = payload["binding"].get("synthesis_profile")
        assert recorded, "a front that records no policy cannot be rechecked under it"
        assert recorded["name"] == sp.BROAD_EXPERIMENTAL.name
        assert recorded["version"] == sp.BROAD_EXPERIMENTAL.version

    def test_every_point_on_the_front_selects_under_that_policy(self, front, tmp_path,
                                                               raw_table):
        """The review asked for all of them, not a sample.

        A front measured as feasible under one policy must stay selectable under it, or the
        front is describing a run nobody can reproduce.
        """
        payload = json.loads(Path(front.artefacts["front"]).read_text(encoding="utf-8"))
        assert payload["observed_front"], "an empty front proves nothing"
        for index in range(len(payload["observed_front"])):
            got = w.select_interface(RECODED, front.artefacts["front"],
                                     tmp_path / f"sel_{index}", index=index,
                                     codon_table=raw_table)
            assert got.ok is True, f"point {index}: {got.failures}"
            assert got.data["synthesis_profile"] == sp.BROAD_EXPERIMENTAL.name
            assert got.data["objectives"], f"point {index} published no objectives"

    def test_a_candidate_that_fails_verification_publishes_nothing(self, front, tmp_path,
                                                                   raw_table, library):
        """Previously ok=True, objectives=None, and the inventory written regardless."""
        payload = json.loads(Path(front.artefacts["front"]).read_text(encoding="utf-8"))
        table = codons.complete_table(raw_table, 1)
        classes = js.junction_classes(library, payload["targets"])

        # Declare the strict policy on a front measured under the broad one, then pick a point
        # that is genuinely infeasible under what selection will now use.
        index = next((i for i, candidate in enumerate(payload["observed_front"])
                      if not js.evaluate(library, classes, candidate["assignment"],
                                         payload["targets"], table,
                                         profile=sp.STRICT_LEGACY).feasible), None)
        if index is None:
            pytest.skip("this front has no strict-infeasible member to exercise the guard")

        payload["binding"]["synthesis_profile"] = sp.STRICT_LEGACY.as_dict()
        mixed = tmp_path / "front_declaring_strict.json"
        mixed.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        outdir = tmp_path / "refused"
        got = w.select_interface(RECODED, mixed, outdir, index=index, codon_table=raw_table)
        assert got.ok is False
        assert got.artefacts == {}, "nothing may be returned after a failed verification"
        assert got.failures[0]["stage"] == "verification"
        assert not (outdir / "inventory_selected.json").exists(), \
            "a rejected inventory must not be left on disk for the next step to find"


# ----------------------------------------------- the level-1 gap, and what filling it produces

class TestLevel1ReadinessIsReachableAndScoresCorrectly:
    """`level1_participants` was tested, library-reachable, and called by nothing.

    A gap a user can see but no workflow can fill is half a fix, so `workflow.level1_readiness`
    exposes it: refuse with the missing roles, or score the reaction when the list is complete.
    """

    CHAIN = {"linker": ["CTTC", "GTGA"], "editing_domain": ["GTGA", "CACG"],
             "bridge": ["CACG", "TTCG"], "final_cassette": ["TTCG", "CGCT"]}
    BLOCKS = ["AATG", "CTTC"]

    def test_an_incomplete_reaction_is_refused_and_names_what_is_missing(self, tmp_path):
        got = w.level1_readiness(self.BLOCKS, tmp_path / "partial",
                                 roles={"linker": ["TGTG"]}, evidence="partial")
        assert got.ok is False
        assert set(got.data["roles_missing"]) == {"final_cassette", "editing_domain", "bridge"}
        assert got.data["scored"] is None, "an incomplete reaction must never be scored"
        assert got.failures[0]["stage"] == "participants"

    def test_a_complete_reaction_is_scored(self, tmp_path):
        got = w.level1_readiness(self.BLOCKS, tmp_path / "complete", roles=self.CHAIN,
                                 evidence="Pryor 2020, Table 2")
        assert got.ok is True
        assert got.data["scored"]["fidelity"] > 0.99

    def test_fidelity_is_scored_over_distinct_junctions_not_contributed_ends(self, tmp_path):
        """The error this task shipped for about ten minutes, and why it is worth a test.

        Every internal overhang is contributed twice, once by the part on each side. Scoring
        the list as supplied makes each junction compete against a perfect copy of itself, and
        the result is not a pessimistic number but a meaningless one: 0.0039 against 0.9940 on
        this very chain. The failure mode is quiet, because a fidelity of 0.0039 looks like a
        badly chosen overhang set rather than a counting error.
        """
        from clippr.overhangs import set_fidelity

        got = w.level1_readiness(self.BLOCKS, tmp_path / "distinct", roles=self.CHAIN,
                                 evidence="Pryor 2020, Table 2")
        scored = got.data["scored"]
        assert scored["junctions"] == sorted(set(scored["junctions"])), "duplicates crept back"
        assert scored["contributed_ends"] > len(scored["junctions"])
        assert scored["fidelity"] == pytest.approx(set_fidelity(scored["junctions"]))

        every_end = ["AATG", "CTTC", "CTTC", "GTGA", "GTGA", "CACG",
                     "CACG", "TTCG", "TTCG", "CGCT"]
        assert set_fidelity(every_end) < 0.01, "the wrong calculation is no longer wrong"

    def test_the_artefact_records_the_refusal_too(self, tmp_path):
        """A refusal is an outcome worth keeping: it says what to go and find."""
        outdir = tmp_path / "refused"
        got = w.level1_readiness(self.BLOCKS, outdir, roles={"linker": ["TGTG"]},
                                 evidence="partial")
        payload = json.loads(Path(got.artefacts["readiness"]).read_text(encoding="utf-8"))
        assert payload["participants_established"] is False
        assert payload["scored"] is None
        assert payload["why_incomplete"]
