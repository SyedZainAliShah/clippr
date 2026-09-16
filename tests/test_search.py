"""Unit tests for design-space exploration and the decision certificate.

The properties worth pinning are about the *decision rule*, not about which design wins: that
hard feasibility can never be traded away, that the hierarchy is applied in the declared order,
and that the certificate reports the alternatives rather than asserting the answer is good.

Selection is tested on constructed candidates rather than real designs, because the measured
degeneracy documented in `search.py` means real candidates differ on only one axis and so
cannot exercise the hierarchy at all.
"""
from __future__ import annotations

import pytest

from clippr.search import (
    Candidate,
    certificate,
    explore,
    select,
)

# The toy table from test_design, so these run without a network call for codon usage.
from test_design import TOY_TABLE


def make(fid: float, score: float, status: str = "PASS", ok: bool = True,
         cuts=(76, 151)) -> Candidate:
    return Candidate(cuts=tuple(cuts), overhangs=("CGAC", "AATG"), fidelity=fid,
                     cds="ATG" * 10, constraints_ok=ok, objectives_score=score,
                     qc={"status": status})


class TestFeasibility:
    def test_qc_fail_is_never_feasible(self):
        assert not make(0.99, 0.0, status="FAIL").feasible

    def test_unsatisfied_constraints_are_never_feasible(self):
        assert not make(0.99, 0.0, ok=False).feasible

    def test_warning_is_feasible(self):
        """WARNING is a design to look at, not a design to exclude."""
        assert make(0.8, 0.0, status="WARNING").feasible

    def test_selection_refuses_when_nothing_is_feasible(self):
        """Returning the least-bad infeasible design would hide a failure."""
        with pytest.raises(ValueError, match="satisfied every constraint"):
            select([make(0.99, 0.0, status="FAIL"), make(0.99, 1.0, ok=False)])

    def test_a_perfect_infeasible_design_never_wins(self):
        """L1 is absolute: the best scores on every other axis cannot buy feasibility."""
        best_but_broken = make(1.0, 999.0, status="FAIL")
        modest = make(0.5, -100.0)
        chosen, _ = select([best_but_broken, modest], fidelity_tolerance=1.0)
        assert chosen is modest


class TestHierarchy:
    def test_fidelity_band_excludes_distant_candidates(self):
        """L2: a much better sequence does not justify falling out of the fidelity band."""
        high_fid = make(0.90, -200.0, cuts=(1, 2))
        low_fid_great_seq = make(0.50, 0.0, cuts=(3, 4))
        chosen, _ = select([high_fid, low_fid_great_seq], fidelity_tolerance=0.02)
        assert chosen is high_fid

    def test_inside_the_band_sequence_quality_decides(self):
        """L3: this is the case the module exists for — fidelity ties, sequence differs."""
        a = make(0.900, -200.0, cuts=(1, 2))
        b = make(0.895, -100.0, cuts=(3, 4))
        chosen, _ = select([a, b], fidelity_tolerance=0.02)
        assert chosen is b, "a 0.005 fidelity concession should buy a much better sequence"

    def test_pass_beats_warning_before_score(self):
        """L3 puts the QC verdict above the optimiser score, not alongside it."""
        warn_but_better = make(0.9, 0.0, status="WARNING", cuts=(1, 2))
        passing = make(0.9, -500.0, status="PASS", cuts=(3, 4))
        chosen, _ = select([warn_but_better, passing])
        assert chosen is passing

    def test_fidelity_breaks_exact_ties(self):
        """L4: everything else equal, take the higher fidelity."""
        lo = make(0.900, -100.0, cuts=(1, 2))
        hi = make(0.910, -100.0, cuts=(3, 4))
        chosen, _ = select([lo, hi], fidelity_tolerance=0.02)
        assert chosen is hi

    def test_widening_the_band_can_change_the_winner(self):
        """The tolerance is a real control, not decoration."""
        high_fid = make(0.90, -200.0, cuts=(1, 2))
        low_fid_great = make(0.50, 0.0, cuts=(3, 4))
        strict, _ = select([high_fid, low_fid_great], fidelity_tolerance=0.02)
        loose, _ = select([high_fid, low_fid_great], fidelity_tolerance=0.50)
        assert strict is high_fid and loose is low_fid_great

    def test_every_feasible_candidate_is_returned_somewhere(self):
        cands = [make(0.9, -1.0, cuts=(1, 2)), make(0.5, -2.0, cuts=(3, 4)),
                 make(0.9, -3.0, cuts=(5, 6))]
        chosen, alts = select(cands, fidelity_tolerance=0.02)
        assert len([chosen, *alts]) == 3


class TestDominance:
    def test_strictly_better_dominates(self):
        assert make(0.9, 0.0).dominates(make(0.8, -1.0))

    def test_equal_does_not_dominate(self):
        assert not make(0.9, 0.0).dominates(make(0.9, 0.0))

    def test_a_mixed_trade_does_not_dominate(self):
        """Better fidelity but worse sequence is a trade-off, not dominance."""
        assert not make(0.9, -10.0).dominates(make(0.8, 0.0))

    def test_worse_qc_blocks_dominance(self):
        assert not make(0.9, 0.0, status="WARNING").dominates(make(0.9, -1.0, status="PASS"))


class TestCertificate:
    def test_states_the_decision_rule(self):
        cands = [make(0.9, -1.0, cuts=(1, 2)), make(0.9, -2.0, cuts=(3, 4))]
        chosen, _ = select(cands)
        text = certificate(chosen, cands)
        assert "Decision rule" in text
        for level in ("QC not FAIL", "within", "prefer QC PASS", "tie-break"):
            assert level in text

    def test_reports_distance_from_the_best_available(self):
        """The question the module exists to answer must appear in its output."""
        cands = [make(0.9, -1.0, cuts=(1, 2)), make(0.95, -50.0, cuts=(3, 4))]
        chosen, _ = select(cands, fidelity_tolerance=0.10)
        text = certificate(chosen, cands)
        assert "vs the best feasible" in text or "best available" in text

    def test_names_the_alternatives_and_their_price(self):
        cands = [make(0.9, -1.0, cuts=(1, 2)), make(0.9, -2.0, cuts=(3, 4))]
        chosen, _ = select(cands)
        text = certificate(chosen, cands)
        assert "Nearest alternatives" in text
        assert "optimiser score" in text

    def test_says_so_when_there_are_no_alternatives(self):
        cands = [make(0.9, -1.0)]
        chosen, _ = select(cands)
        assert "No feasible alternative" in certificate(chosen, cands)

    def test_reports_dominance(self):
        cands = [make(0.9, -1.0, cuts=(1, 2)), make(0.9, -2.0, cuts=(3, 4))]
        chosen, _ = select(cands)
        assert "dominated by" in certificate(chosen, cands)

    def test_carries_the_caveat(self):
        cands = [make(0.9, -1.0)]
        chosen, _ = select(cands)
        text = certificate(chosen, cands)
        assert "not measured assembly" in text
        assert "this protein only" in text

    def test_infeasible_set_yields_an_honest_certificate(self):
        cands = [make(0.9, 0.0, status="FAIL")]
        assert "No feasible design" in certificate(cands[0], cands)


@pytest.fixture(scope="module")
def candidates():
    """One real exploration, shared across the tests below — each costs a codon run."""
    from clippr.ppr import describe
    protein = describe("AAAAUGUGG")["aa_sequence"]
    return explore(protein, TOY_TABLE, n_fragments=4, destination=None,
                   matrix="BsaI-HFv2", budget=3)


class TestExplore:
    """One real run. Slow enough to keep to a single small case."""

    def test_evaluates_the_requested_number(self, candidates):
        assert len(candidates) == 3

    def test_every_candidate_is_a_finished_design(self, candidates):
        """The point of the module: compare finished designs, not predictions about them."""
        for c in candidates:
            assert c.cds and len(c.cds) % 3 == 0
            assert c.qc.get("status") in ("PASS", "WARNING", "FAIL")

    def test_candidates_are_distinct_plans(self, candidates):
        assert len({c.cuts for c in candidates}) == len(candidates)

    def test_selection_works_on_real_candidates(self, candidates):
        chosen, alts = select(candidates)
        assert chosen in candidates
        assert len(alts) == len(candidates) - 1
