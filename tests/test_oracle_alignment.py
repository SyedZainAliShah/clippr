"""Comparing codon adaptation across systems requires the same coding region.

The first version of the matched oracle benchmark scored `cai(ours)` against `cai(theirs)`
directly. Only **4 of 42** modules have identical translated spans; in the other 38 the
reference carries additional boundary amino acids, so the two numbers described different
regions. The reported mean and win count were 0.634781 and 29; over the shared span they are
0.658672 and 22.

Neither representation is wrong. Comparing them unaligned was, and it is the same error as
enforcing a constraint on a sub-span of what ships -- moved from the constraint into the
comparison, which is why it survived four passes that were all looking at constraints.

The numbers below are asserted exactly so that a regression to unaligned scoring fails here
rather than being published.
"""
from __future__ import annotations

import csv
import json
import math
import statistics
import sys
from pathlib import Path

import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation" / "experiments"))

from matched_oracle_benchmark import align_coding_spans        # noqa: E402

REFERENCE_CSV = ROOT / "work" / "m3" / "reference" / "output" / "optimized_library.csv"
INVENTORY = ROOT / "work" / "phaseb" / "inventory_recoded.json"
CODON_TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"


class TestTheAlignmentItself:
    def test_an_identical_span_aligns_at_zero(self):
        cds = "ATGAAACCCGGG"
        got = align_coding_spans(cds, cds)
        assert got["aligned"] and got["identical_span"]
        assert got["reference_nt_bounds"] == [0, len(cds)]
        assert got["reference_trimmed_aa"] == 0

    def test_a_contained_span_reports_its_bounds(self):
        """The real shape: the reference carries extra boundary residues on both sides."""
        inner = "AAACCCGGG"
        outer = "ATG" + inner + "TAT"
        got = align_coding_spans(inner, outer)
        assert got["aligned"] and not got["identical_span"]
        assert got["reference_nt_bounds"] == [3, 3 + len(inner)]
        assert got["reference_span"] == inner
        assert got["reference_trimmed_aa"] == 2

    def test_an_absent_protein_is_refused_not_guessed(self):
        got = align_coding_spans("ATGAAACCC", "ATGGGGTTTCAT")
        assert not got["aligned"]
        assert "absent" in got["reason"]

    def test_an_ambiguous_protein_is_refused_rather_than_resolved(self):
        """Two valid offsets means no defensible choice, so it must refuse."""
        inner = "AAACCC"
        got = align_coding_spans(inner, inner + "GGG" + inner)
        assert not got["aligned"]
        assert "2 times" in got["reason"]

    def test_alignment_is_on_protein_not_nucleotides(self):
        """Synonymous codons must still align; that is the whole point of the comparison."""
        ours = "AAACCC"                      # Lys Pro
        theirs = "ATG" + "AAGCCG" + "TAT"    # Met Lys Pro Tyr, both synonymous
        got = align_coding_spans(ours, theirs)
        assert got["aligned"]
        assert got["reference_span"] == "AAGCCG"


@pytest.mark.skipif(not REFERENCE_CSV.is_file() or not INVENTORY.is_file(),
                    reason="the saved reference output and recoded inventory are needed")
class TestAgainstTheSavedCorpus:
    @staticmethod
    def rows():
        from clippr import inventories as inv

        lib = inv.load(INVENTORY)
        raw = json.loads(CODON_TABLE.read_text(encoding="utf-8"))
        weights = {c.replace("U", "T"): v / max(d.values())
                   for aa, d in raw.items() for c, v in d.items() if max(d.values()) > 0}
        code = CodonTable.unambiguous_dna_by_id[1].forward_table

        def cai(dna):
            used = [weights[c] for i in range(0, len(dna) - 2, 3) for c in [dna[i:i + 3]]
                    if c in weights and code.get(c) not in (None, "M", "W")]
            return (math.exp(sum(math.log(w) for w in used) / len(used))
                    if used and min(used) > 0 else 0.0)

        out = []
        for row in csv.DictReader(REFERENCE_CSV.open(encoding="utf-8")):
            name = row["optimized_part_id"].removesuffix("_v1")
            record = lib.modules["pPR-1_" + name]
            ours = record.dna[slice(*record.coding_interval)]
            got = align_coding_spans(ours, row["optimized_cds"])
            out.append({"module": name, "alignment": got, "ours": cai(ours),
                        "reference_aligned": cai(got["reference_span"]) if got["aligned"] else None,
                        "reference_whole": cai(row["optimized_cds"])})
        return out

    def test_most_spans_are_not_identical(self):
        """The fact that makes unaligned comparison wrong, asserted rather than assumed."""
        rows = self.rows()
        identical = sum(r["alignment"]["identical_span"] for r in rows)
        assert len(rows) == 42
        assert identical == 4, f"{identical}/42 identical; the corpus changed"

    def test_every_module_aligns_unambiguously(self):
        assert all(r["alignment"]["aligned"] for r in self.rows())

    def test_the_aligned_comparison_is_the_reported_one(self):
        rows = self.rows()
        ours = statistics.mean(r["ours"] for r in rows)
        reference = statistics.mean(r["reference_aligned"] for r in rows)
        ahead = sum(r["ours"] > r["reference_aligned"] for r in rows)
        assert ours == pytest.approx(0.621004, abs=1e-6)
        assert reference == pytest.approx(0.658672, abs=1e-6)
        assert ahead == 22

    def test_the_unaligned_comparison_would_report_something_else(self):
        """The regression guard: scoring whole CDSs gives a different, wrong answer."""
        rows = self.rows()
        reference = statistics.mean(r["reference_whole"] for r in rows)
        ahead = sum(r["ours"] > r["reference_whole"] for r in rows)
        assert reference == pytest.approx(0.634781, abs=1e-6)
        assert ahead == 29
        assert ahead != 22, "the two comparisons must not coincide, or this proves nothing"
