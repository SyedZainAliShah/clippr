"""One validator, every acceptance path.

Three probes, each of which passed before `substrates.substrate_problems` became the single
feasibility test, and each of which must now be refused:

  * greedy seed 42 on the recoded inventory, which returned a module whose ordered substrate
    had a 0.660 GC window (16 of 18 matched runs did)
  * a protein-preserving candidate injected straight into the joint search, bypassing the
    recoder -- the check must not depend on the recoder being right
  * an internal BsaI site, which `synthesis_problems` did not look for at all

They live together rather than in the three modules they exercise because they are one defect:
a constraint checked on a sub-span of what is delivered is a check of something else.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
from Bio.Data import CodonTable
from Bio.Restriction import BsaI
from Bio.Seq import Seq

from clippr import inventories as inv
from clippr import joint_search as js
from clippr import substrates as sub
from clippr.library_search import MODES, optimise_library

ROOT = Path(__file__).resolve().parents[1]
TABLE_S1 = ROOT / "data" / "grasp_supp" / "Table S1.xlsx"
CODON_TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"
RECODED = ROOT / "work" / "phaseb" / "inventory_recoded.json"

#: The module the 0.660 window was found in. An E module: its released fragment is the insert
#: plus the GC-rich `CGAG` overhang, so it sits closest to the ceiling of all 42.
PROBE_MODULE = "pPR-1_19E_LN5N"


def worst_gc(sequence: str) -> float:
    """The highest GC fraction over any sliding window, computed here rather than imported.

    A helper borrowed from the code under test would only prove the check agrees with itself.
    """
    from clippr.recoding import GC_WINDOW

    windows = [sequence[i:i + GC_WINDOW]
               for i in range(max(1, len(sequence) - GC_WINDOW + 1))]
    return max((w.count("G") + w.count("C")) / len(w)
               for w in windows if len(w) == GC_WINDOW)


@pytest.fixture(scope="module")
def deposited():
    if not TABLE_S1.is_file():
        pytest.skip(f"supply {TABLE_S1.name} to run substrate-scope tests")
    return inv.load_deposited(TABLE_S1)


@pytest.fixture(scope="module")
def table():
    return json.loads(CODON_TABLE.read_text(encoding="utf-8"))


def gc_maximised(dna: str, coding: tuple[int, int]) -> str:
    """The same protein, every codon swapped for its most GC-rich synonym.

    Deterministic and independent of the recoder, so the probe cannot drift when the recoder
    changes. On `pPR-1_19E_LN5N` it lands at GC 0.660 -- the observed breach, not a
    hypothetical one.
    """
    forward = CodonTable.unambiguous_dna_by_id[1].forward_table
    synonyms: dict[str, list[str]] = {}
    for codon, aa in forward.items():
        synonyms.setdefault(aa, []).append(codon)
    richest = {aa: max(codons, key=lambda c: (sum(b in "GC" for b in c), c))
               for aa, codons in synonyms.items()}

    start, end = coding
    cds = dna[start:end]
    swapped = "".join(richest.get(str(Seq(cds[i:i + 3]).translate()), cds[i:i + 3])
                      for i in range(0, len(cds), 3))
    return dna[:start] + swapped + dna[end:]


@pytest.fixture(scope="module")
def injected(deposited):
    """A candidate that keeps the protein and the length, and breaches on the substrate."""
    record = deposited.modules[PROBE_MODULE]
    return gc_maximised(record.dna, record.coding_interval)


class TestTheInjectedCandidate:
    def test_it_really_does_preserve_the_protein(self, deposited, injected):
        """Otherwise it would be rejected for the wrong reason and prove nothing."""
        record = deposited.modules[PROBE_MODULE]
        start, end = record.coding_interval
        assert len(injected) == len(record.dna)
        assert (str(Seq(injected[start:end]).translate())
                == str(Seq(record.dna[start:end]).translate()))

    def test_the_bare_insert_is_where_the_breach_hides(self, deposited, injected):
        """The insert alone is under the ceiling; the released fragment is not.

        This is the whole defect in two assertions. Checking the left-hand one was the bug.
        """
        record = deposited.modules[PROBE_MODULE]
        fragment = sub.released_fragment(injected, record.block)
        assert len(fragment) > len(injected)
        assert worst_gc(fragment) > worst_gc(injected)

    def test_the_shared_validator_refuses_it(self, deposited, injected):
        problems = sub.substrate_problems(injected, deposited.modules[PROBE_MODULE].block)
        assert problems
        assert any("GC" in p for p in problems)

    def test_the_joint_search_refuses_it_without_consulting_the_recoder(
            self, deposited, table, injected, monkeypatch):
        """Injected past `apply_assignment`, so a correct recoder cannot mask the check."""
        from clippr.codons import complete_table

        targets = ["AAAAUGUGG"]
        classes = js.junction_classes(deposited, targets)
        assignment = {c.key: c.deposited for c in classes}
        rigged = inv.derive(deposited, "injected", {PROBE_MODULE: injected})
        monkeypatch.setattr(js, "apply_assignment", lambda *a, **k: rigged)

        got = js.evaluate(deposited, classes, assignment, targets,
                          complete_table(table, 1))
        assert not got.feasible
        assert PROBE_MODULE in got.reason and "GC" in got.reason

    def test_the_clean_deposited_module_is_still_accepted(self, deposited):
        """The check must reject the probe, not everything that reaches it."""
        record = deposited.modules[PROBE_MODULE]
        assert sub.substrate_problems(record.dna, record.block) == []


class TestAnInternalRecognitionSite:
    """`synthesis_problems` scored GC and homopolymers and looked for no enzyme sites at all.

    A BsaI site inside a level-0 insert is not hypothetical: the block that carries it is cut
    by BsaI at level 1, so the internal site would cut the assembled block apart.
    """

    def test_an_internal_bsai_site_is_caught(self):
        insert = "ACTC" + "ATGCGT" * 3 + "GGTCTC" + "ATGCGT" * 5 + "AAGA"
        assert BsaI.search(Seq(insert), linear=True)          # the probe is what it claims
        problems = sub.substrate_problems(insert, "B")
        assert any("BsaI" in p for p in problems), problems

    def test_the_reverse_strand_counts_too(self):
        """A site on the bottom strand cuts exactly as well as one on the top."""
        insert = "ACTC" + "ATGCGT" * 3 + "GAGACC" + "ATGCGT" * 5 + "AAGA"
        problems = sub.substrate_problems(insert, "B")
        assert any("BsaI" in p and "reverse" in p for p in problems), problems

    def test_the_wrappers_own_bbsi_sites_are_not_counted(self):
        """Both are the construction. Counting them would fail all 42 substrates."""
        insert = "ACTC" + "ATGCAC" * 12 + "AAGA"
        assert sub.substrate_problems(insert, "B") == []


def offenders_in(inventory) -> dict[str, tuple[str, ...]]:
    return {s.module_id: s.synthesis_problems
            for s in sub.substrates_for(inventory) if s.synthesis_problems}


@pytest.mark.usefixtures("deposited", "table")
class TestTheOptimiserJudgesWhatItShips:
    """Judged as *introduced*, not absolute.

    One deposited module already carries a run of five Ts and so fails CLIPPR's homopolymer
    rule on arrival. Demanding a clean output would make the deposited inventory unusable as
    a baseline; what the optimiser owes is that it adds nothing.
    """

    def test_greedy_seed_42_introduces_no_unorderable_module(self, deposited, table):
        """The named reproducer. `greedy/s42/fwd` returned a 0.660 GC substrate."""
        got = optimise_library(deposited, table, mode="greedy", seed=42,
                               max_proposals=42, wall_seconds=600)
        inherited = offenders_in(deposited)
        assert {m: p for m, p in offenders_in(got.inventory).items()
                if p != inherited.get(m)} == {}

    def test_no_mode_introduces_one(self, deposited, table):
        inherited = offenders_in(deposited)
        for mode in MODES:
            got = optimise_library(deposited, table, mode=mode, seed=42,
                                   max_proposals=20, wall_seconds=600)
            new = {m: p for m, p in offenders_in(got.inventory).items()
                   if p != inherited.get(m)}
            assert new == {}, f"{mode}: {new}"

    def test_the_recoded_inventory_on_disk_is_orderable(self):
        """The artefact the rest of the workflow consumes, not a fresh run of the search."""
        if not RECODED.is_file():
            pytest.skip("no recoded inventory on disk")
        assert offenders_in(inv.load(RECODED)) == {}
