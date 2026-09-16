"""The third Pareto axis, and the measurement that condemned the one it replaced.

`repeat_burden` was identically zero over both the deposited and the recoded inventory, so a
search reporting three objectives was ordering by two. The first test below pins that fact,
because the tempting fix — lowering `k` — does not work either, and someone will try it.

The replacement is a composite over quantities that do vary on this corpus. `test_it_is_not
_degenerate_on_the_real_corpus` is the one that matters: an axis with no spread is not an axis,
and if a future inventory flattens it, that test fails rather than the search quietly returning
to two objectives.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from clippr import synthesis_fitness as sf
from clippr import synthesis_profile as sp

ROOT = Path(__file__).resolve().parents[1]
TABLE_S1 = ROOT / "data" / "grasp_supp" / "Table S1.xlsx"
RECODED = ROOT / "work" / "phaseb" / "inventory_recoded.json"


class TestWhyTheOldAxisWasReplaced:
    @pytest.mark.skipif(not TABLE_S1.is_file(), reason="needs the deposited kit")
    def test_repeat_burden_is_zero_at_every_plausible_k(self):
        """The fix that looks obvious and is not.

        The reference uses k=16 where we used 20, which invites the conclusion that the
        parameter was the problem. It is zero at 20, 16, 12 and 10 alike.
        """
        from clippr import inventories as inv
        from clippr.objectives import duplicated_kmers

        deposited = inv.load_deposited(TABLE_S1)
        for k in (20, 16, 12, 10):
            total = sum(duplicated_kmers(m.dna, k) for m in deposited.modules.values())
            assert total == 0, f"k={k} is no longer degenerate; revisit the axis rationale"


class TestTheComponents:
    def test_gc_centrality_is_one_at_the_band_centre(self):
        """Centre of (0.35, 0.65) is 0.50, so an even sequence scores 1.0."""
        # Every 50-nt window must be exactly 0.50 GC. "ATGC" repeated does not manage that:
        # 50 is not a multiple of 4, so its windows oscillate between 0.48 and 0.52.
        even = "AG" * 60
        assert sf.gc_centrality(even, (0.35, 0.65), 50) == pytest.approx(1.0, abs=1e-9)

    def test_gc_centrality_is_zero_at_the_band_edge(self):
        assert sf.gc_centrality("GC" * 60, (0.35, 0.65), 50) == pytest.approx(0.0, abs=1e-9)
        assert sf.gc_centrality("AT" * 60, (0.35, 0.65), 50) == pytest.approx(0.0, abs=1e-9)

    def test_gc_centrality_follows_the_declared_band_not_a_fixed_half(self):
        """A profile permitting 0.15-0.85 is not asking for the same sequence as 0.35-0.65.

        Scoring both toward 0.5 would impose a preference no profile declared.
        """
        # 0.75 GC throughout: outside the narrow band entirely, inside the wide one.
        sequence = "AGGC" * 30
        narrow = sf.gc_centrality(sequence, (0.35, 0.65), 50)
        wide = sf.gc_centrality(sequence, (0.15, 0.85), 50)
        assert wide > narrow, "a wider band must be more forgiving of the same sequence"

    def test_the_worst_window_sets_the_score_not_the_mean(self):
        """One bad window is what a vendor's model reacts to; averaging would hide it."""
        mostly_even = "ATGC" * 25 + "GCGCGCGCGCGCGC" * 4
        worst = sf.gc_centrality(mostly_even, (0.35, 0.65), 50)
        assert worst < sf.gc_centrality("ATGC" * 40, (0.35, 0.65), 50)

    def test_homopolymer_headroom_is_a_ramp_not_a_threshold(self):
        """A run of 2 must beat a run of 4 even though both are feasible at cap 4."""
        short = sf.homopolymer_headroom("ATGCATGCAAGCAT", cap=4)
        longer = sf.homopolymer_headroom("ATGCATGCAAAAGC", cap=4)
        assert short > longer
        assert 0.0 <= longer <= 1.0

    def test_homopolymer_headroom_bottoms_out_at_the_cap(self):
        assert sf.homopolymer_headroom("A" * 8 + "TGCA", cap=4) == pytest.approx(0.0)

    def test_internal_repetition_fires_on_a_repetitive_sequence(self):
        """Degenerate on today's corpus, which is why it is one term of four and not the axis."""
        unique = "".join("ATGC"[(i * 7 + i // 4) % 4] for i in range(200))
        repetitive = "ATGCATGCATGCATGCATGCATGC" * 8
        assert sf.internal_repetition(repetitive, 20) < sf.internal_repetition(unique, 20)
        assert sf.internal_repetition("ATGCATGCATGCATGCATGCAT", 20) == pytest.approx(1.0)

    def test_collection_sharing_falls_when_a_neighbour_shares_kmers(self):
        sequence = "ATGCATGCATGCATGCATGCATGCATGCATGC"
        alone = sf.collection_sharing(sequence, [], 20)
        beside_itself = sf.collection_sharing(sequence, [sequence], 20)
        assert alone == pytest.approx(1.0)
        assert beside_itself < alone

    def test_an_empty_collection_is_not_a_penalty(self):
        """One module on its own shares nothing; that is not a defect."""
        assert sf.collection_sharing("ATGC" * 20, [], 20) == pytest.approx(1.0)


class TestTheCompositeScore:
    def test_the_total_is_bounded_and_reports_its_parts(self):
        got = sf.synthesis_fitness("ATGC" * 40)
        assert 0.0 <= got.total <= 1.0
        assert set(got.components) == set(sf.WEIGHTS)
        assert got.weights == sf.WEIGHTS
        assert got.k == sf.DEFAULT_K

    def test_the_weights_are_declared_with_every_score(self):
        """They are an engineering choice, not a measurement, so they travel with the number."""
        got = sf.synthesis_fitness("ATGC" * 40).as_dict()
        assert got["weights"] and abs(sum(got["weights"].values()) - 1.0) < 1e-9
        assert got["scope"] == "ordered substrate"

    def test_the_profile_changes_the_score(self):
        """The band is part of the measurement, so two policies must not agree by accident."""
        sequence = "AGGC" * 30
        strict = sf.synthesis_fitness(sequence, profile=sp.STRICT_LEGACY).total
        broad = sf.synthesis_fitness(sequence, profile=sp.BROAD_EXPERIMENTAL).total
        assert strict != broad

    def test_zero_weights_are_refused_rather_than_dividing_by_zero(self):
        with pytest.raises(ValueError, match="positive"):
            sf.synthesis_fitness("ATGC" * 40, weights={"gc_centrality": 0.0})


@pytest.mark.skipif(not TABLE_S1.is_file(), reason="needs the deposited kit")
class TestOnTheRealCorpus:
    @staticmethod
    def kit():
        from clippr import inventories as inv
        return inv.load_deposited(TABLE_S1)

    def test_it_is_not_degenerate_on_the_real_corpus(self):
        """**The** reason this axis exists. An axis with no spread cannot rank anything."""
        got = sf.inventory_fitness(self.kit())
        assert got["modules"] == 42
        assert not got["degenerate"]
        assert got["spread"] > 0.01, got
        assert got["max"] > got["min"] + 0.1, "the axis barely separates the corpus"

    def test_degeneracy_is_reported_rather_than_hidden(self):
        """If a future inventory flattens it, the flag says so instead of it passing silently."""
        got = sf.inventory_fitness(self.kit())
        assert "degenerate" in got and "spread" in got

    def test_every_module_gets_a_score_with_its_parts(self):
        got = sf.inventory_fitness(self.kit())
        assert len(got["per_module"]) == 42
        for payload in got["per_module"].values():
            assert set(payload["components"]) == set(sf.WEIGHTS)
            assert 0.0 <= payload["total"] <= 1.0


class TestItIsWiredIntoTheSearch:
    def test_the_declared_directions_use_it(self):
        from clippr import joint_search as js

        assert js.DIRECTIONS == {"fidelity": "max", "adaptation": "max",
                                 "synthesis_fitness": "max"}
        assert "repeat_burden" not in js.DIRECTIONS, "the degenerate axis must not rank"

    @pytest.mark.skipif(not RECODED.is_file(), reason="needs the recoded inventory")
    def test_the_search_reports_it_and_keeps_the_old_axis_as_a_diagnostic(self):
        import json

        from clippr import inventories as inv
        from clippr import joint_search as js
        from clippr.codons import complete_table

        table = complete_table(json.loads(
            (ROOT / "data" / "codon_tables" / "kazusa_3055.json").read_text("utf-8")), 1)
        library = inv.load(RECODED)
        targets = ["AAAAUGUGG"]
        classes = js.junction_classes(library, targets)
        got = js.evaluate(library, classes, {c.key: c.deposited for c in classes},
                          targets, table)
        assert got.feasible
        assert 0.0 <= got.objectives["synthesis_fitness"] <= 1.0
        # retained, reported, and not ranking
        assert "repeat_burden" in got.objectives
        assert "synthesis_fitness_spread" in got.objectives
