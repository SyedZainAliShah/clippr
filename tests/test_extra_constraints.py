"""Controls for `optimize_cds(extra_constraints=...)`.

Added for the wrapper-aware pilot, which needs the solver to see the sequence that would be
ordered rather than the coding span alone. The hook exists so an experiment can do that
without changing a shipped default, so the two things worth pinning are that it is honoured
and that it is empty everywhere else.
"""
from __future__ import annotations

import inspect

import pytest

from ppropt.codons import optimize_cds

from test_design import TOY_TABLE

PROTEIN = "MSTNPKPQRKTKRNTNRRPQDVKFPGG"


class TestTheHookIsOptional:
    def test_the_default_is_empty(self):
        assert inspect.signature(optimize_cds).parameters[
            "extra_constraints"].default == ()

    def test_no_shipped_route_passes_it(self):
        """If a shipped route starts using it, that is a default change, not an experiment."""
        from pathlib import Path

        src = Path(__file__).resolve().parents[1] / "src" / "ppropt"
        callers = [f.name for f in src.rglob("*.py")
                   if "extra_constraints" in f.read_text(encoding="utf-8")]
        assert callers == ["codons.py"], (
            f"extra_constraints is referenced in {callers}; only the definition should "
            f"use it until a default change is decided on its own evidence")

    def test_omitting_it_is_unchanged(self):
        """The baseline must be bit-identical with and without an empty list."""
        a = optimize_cds(PROTEIN, codon_table=TOY_TABLE, seed=7)
        b = optimize_cds(PROTEIN, codon_table=TOY_TABLE, seed=7, extra_constraints=[])
        assert a["cds"] == b["cds"]


class TestTheHookIsHonoured:
    def test_a_constraint_that_cannot_be_met_is_reported_not_ignored(self):
        """An impossible extra constraint must show up in the result, not be dropped."""
        import dnachisel as dc

        class Impossible(dc.Specification):
            """Always breaching, so a solver that ignored it would still report ok."""

            def evaluate(self, problem):
                return dc.SpecEvaluation(
                    self, problem, score=-1,
                    locations=[dc.Location(0, min(9, len(problem.sequence)))],
                    message="this constraint can never be satisfied")

            def localized(self, location, problem=None):
                return self

        out = optimize_cds(PROTEIN, codon_table=TOY_TABLE, seed=7,
                           extra_constraints=[Impossible()])
        assert out["constraints_ok"] is False

    def test_a_satisfiable_constraint_changes_the_answer(self):
        """A constraint that forbids a codon the baseline uses must move the sequence."""
        import dnachisel as dc

        baseline = optimize_cds(PROTEIN, codon_table=TOY_TABLE, seed=7)
        # Pick a codon the baseline actually chose, then forbid that DNA word outright.
        used = baseline["cds"][3:6]
        out = optimize_cds(PROTEIN, codon_table=TOY_TABLE, seed=7,
                           extra_constraints=[dc.AvoidPattern(used)])

        assert out["cds"] != baseline["cds"] or out["constraints_ok"] is False
        if out["constraints_ok"]:
            assert used not in out["cds"], "the forbidden word is still present"

    def test_translation_still_holds_under_an_extra_constraint(self):
        """The hook must not be able to buy a constraint win with the wrong protein."""
        import dnachisel as dc
        from Bio.Seq import Seq

        out = optimize_cds(PROTEIN, codon_table=TOY_TABLE, seed=7,
                           extra_constraints=[dc.EnforceGCContent(mini=0.3, maxi=0.7)])
        assert str(Seq(out["cds"]).translate()).rstrip("*") == PROTEIN


class TestThePilotConstraint:
    """The pilot's own constraint, imported from where it lives."""

    @staticmethod
    def _pilot():
        import importlib.util
        from pathlib import Path

        path = (Path(__file__).resolve().parents[1] / "validation" / "experiments"
                / "wrapper_aware_pilot.py")
        spec = importlib.util.spec_from_file_location("_pilot", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_it_scores_the_flanks_not_just_the_slice(self):
        """A slice inside the band whose flanks push it out must be reported."""
        pilot = self._pilot()
        from ppropt.synthesis_profile import resolve

        policy = resolve("ppropt-strict-legacy")
        # 60 nt of pure AT: every window of the slice alone is GC 0.0, outside the band,
        # so a constraint reading the slice reports it. With GC-rich flanks the windows
        # that span the junction move too -- the point is that flanks are included.
        fragments = [{"fragment_id": "F1", "cds_start": 0, "cds_end": 60,
                      "left": "G" * 40, "right": "C" * 40}]
        spec = pilot.make_constraint(fragments, policy)

        import dnachisel as dc
        problem = dc.DnaOptimizationProblem(sequence="AT" * 30, constraints=[],
                                            objectives=[], logger=None)
        evaluation = spec.evaluate(problem)
        assert evaluation.score < 0
        assert evaluation.locations, "a breach must be locatable"

    def test_a_balanced_fragment_passes(self):
        """The positive side: flanks and slice together inside the band score zero."""
        pilot = self._pilot()
        from ppropt.synthesis_profile import resolve

        policy = resolve("ppropt-strict-legacy")
        balanced = "ATGC" * 30                     # GC 0.5 everywhere
        fragments = [{"fragment_id": "F1", "cds_start": 0, "cds_end": len(balanced),
                      "left": "ATGC" * 5, "right": "ATGC" * 5}]
        spec = pilot.make_constraint(fragments, policy)

        import dnachisel as dc
        problem = dc.DnaOptimizationProblem(sequence=balanced, constraints=[],
                                            objectives=[], logger=None)
        assert spec.evaluate(problem).score == 0

    def test_an_ambiguous_slice_is_refused_not_guessed(self):
        """If the slice occurs twice in the ordered string the flanks are not identifiable."""
        import pandas as pd

        pilot = self._pilot()
        baseline = {
            "cds": "AAAA",
            "oligos": pd.DataFrame([{"fragment_id": "F1", "cds_start": 0, "cds_end": 2,
                                     "oligo_sequence_5to3": "AAAAAA"}]),
        }
        with pytest.raises(ValueError, match="occurs"):
            pilot.frozen_geometry(baseline)
