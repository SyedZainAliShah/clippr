"""Unit tests for fragment splitting and Golden Gate wrapping.

The reverse-complement tests come first deliberately: a 3' site written in its forward
spelling is the one failure here that is invisible in a sequence file and only shows up as
an assembly that does not work.
"""
from __future__ import annotations

import pytest

from clippr.assembly import (
    DEFAULT_ENZYME,
    PAD_3,
    PAD_5,
    Fragment,
    build_oligos,
    enzyme_geometry,
    reassemble,
    reverse_complement,
    split_cds,
    wrap_fragment,
)

CDS = "".join(["ATG", "AAA", "CCC", "GGG", "TTT", "ACA", "GCA", "TCA",
               "GAT", "CGT", "TAC", "AAG", "CTG", "GTC", "TTA", "CCA"] * 4)
CUTS = [16, 32, 48]


class TestReverseComplementTrap:
    """The named trap: both Type IIS sites must face inward."""

    def test_forward_site_appears_exactly_once(self):
        site, _, _ = enzyme_geometry("BsaI")
        oligo = wrap_fragment("ATGAAACCCGGG")
        assert oligo.count(site) == 1, "the 3' site is in its forward spelling"

    def test_three_prime_arm_carries_the_reverse_complement(self):
        site, _, _ = enzyme_geometry("BsaI")
        oligo = wrap_fragment("ATGAAACCCGGG")
        assert reverse_complement(site) in oligo
        assert oligo.index(site) < oligo.index(reverse_complement(site))

    def test_sites_point_inward_at_the_payload(self):
        """Read from each end, both sites must run towards the payload, not away."""
        site, spacer, _ = enzyme_geometry("BsaI")
        payload = "ATGAAACCCGGG"
        oligo = wrap_fragment(payload)
        assert oligo.startswith(PAD_5 + site)
        assert oligo.endswith(reverse_complement(PAD_3 + site + "A" * spacer))

    def test_a_forward_three_prime_arm_would_be_caught(self):
        """Sanity-check the test itself: the wrong construction must fail the assertion."""
        site, spacer, _ = enzyme_geometry("BsaI")
        wrong = f"{PAD_5}{site}{'A'*spacer}ATGAAACCCGGG{'A'*spacer}{site}{PAD_3}"
        assert wrong.count(site) == 2

    def test_every_enzyme_wraps_inward(self):
        for enzyme in ("BsaI", "BbsI", "BsmBI", "SapI"):
            site, _, _ = enzyme_geometry(enzyme)
            oligo = wrap_fragment("ATGAAACCCGGG", enzyme=enzyme)
            assert oligo.count(site) == 1, enzyme
            assert reverse_complement(site) in oligo, enzyme


class TestEnzymeGeometry:
    @pytest.mark.parametrize("enzyme,site,spacer,overhang", [
        ("BsaI", "GGTCTC", 1, 4),
        ("BbsI", "GAAGAC", 2, 4),
        ("BsmBI", "CGTCTC", 1, 4),
        ("SapI", "GCTCTTC", 1, 3),
    ])
    def test_published_geometry(self, enzyme, site, spacer, overhang):
        assert enzyme_geometry(enzyme) == (site, spacer, overhang)

    def test_unknown_enzyme_rejected(self):
        with pytest.raises(ValueError, match="unknown restriction enzyme"):
            enzyme_geometry("NotAnEnzyme")

    def test_blunt_cutter_rejected(self):
        """EcoRV cuts inside its own site, so it cannot wrap a fragment."""
        with pytest.raises(ValueError, match="not a Type IIS"):
            enzyme_geometry("EcoRV")


class TestSplit:
    def test_fragment_count_is_cuts_plus_one(self):
        assert len(split_cds(CDS, CUTS)) == len(CUTS) + 1

    def test_spans_tile_the_sequence(self):
        frags = split_cds(CDS, CUTS)
        assert frags[0].cds_start == 0
        assert frags[-1].cds_end == len(CDS)
        for a, b in zip(frags, frags[1:]):
            assert a.cds_end == b.cds_start

    def test_junction_overhang_is_shared(self):
        """Consecutive payloads must overlap by exactly the overhang."""
        frags = split_cds(CDS, CUTS)
        for a, b in zip(frags, frags[1:]):
            assert a.payload[-4:] == b.payload[:4] == a.oh3_coding_site

    def test_outer_ends_carry_the_destination(self):
        frags = split_cds(CDS, CUTS, destination=("CTCA", "CGAG"))
        assert frags[0].payload.startswith("CTCA")
        assert frags[-1].payload.endswith("CGAG")

    def test_three_prime_end_overhang_is_flipped(self):
        """The reported end overhang is the strand the enzyme leaves, not the coding site."""
        df = build_oligos(CDS, CUTS, destination=("CTCA", "CGAG"))
        assert df.iloc[-1]["oh3_coding_site_5to3"] == "CGAG"
        assert df.iloc[-1]["three_prime_end_overhang"] == "CTCG"

    def test_overhangs_are_cross_checked(self):
        good = [CDS[3 * c - 4:3 * c] for c in CUTS]
        split_cds(CDS, CUTS, overhangs=good)
        with pytest.raises(ValueError, match="locked site was not preserved"):
            split_cds(CDS, CUTS, overhangs=["AAAA"] + good[1:])

    def test_overhang_count_must_match(self):
        with pytest.raises(ValueError, match="overhangs for"):
            split_cds(CDS, CUTS, overhangs=["AAAA"])

    def test_unsorted_cuts_rejected(self):
        with pytest.raises(ValueError, match="strictly ascending"):
            split_cds(CDS, [32, 16])

    def test_repeated_cut_rejected(self):
        with pytest.raises(ValueError, match="strictly ascending"):
            split_cds(CDS, [16, 16])

    def test_non_codon_sequence_rejected(self):
        with pytest.raises(ValueError, match="whole number of codons"):
            split_cds(CDS + "AT", CUTS)

    @pytest.mark.parametrize("cut", [0, 1000])
    def test_out_of_range_cut_rejected(self, cut):
        with pytest.raises(ValueError, match="no room"):
            split_cds(CDS, [cut])


class TestRoundTrip:
    def test_reassembly_reproduces_the_cds(self):
        assert reassemble(split_cds(CDS, CUTS)) == CDS

    @pytest.mark.parametrize("cuts", [[16], [16, 32], [16, 32, 48], [8, 24, 40, 56]])
    def test_round_trip_for_various_splits(self, cuts):
        assert reassemble(split_cds(CDS, cuts)) == CDS

    def test_single_fragment_round_trips(self):
        assert reassemble(split_cds(CDS, [])) == CDS

    def test_empty_rejected(self):
        with pytest.raises(ValueError, match="no fragments"):
            reassemble([])


class TestOligos:
    def test_oligo_is_payload_plus_both_arms(self):
        df = build_oligos(CDS, CUTS)
        site, spacer, _ = enzyme_geometry(DEFAULT_ENZYME)
        arm = len(PAD_5) + len(site) + spacer
        for _, r in df.iterrows():
            assert r["oligo_length"] == len(r["payload_5to3"]) + 2 * arm
            assert r["payload_5to3"] in r["oligo_sequence_5to3"]

    def test_assembly_order_is_sequential(self):
        df = build_oligos(CDS, CUTS)
        assert list(df["assembly_order"]) == [1, 2, 3, 4]

    def test_prefix_names_the_order_ids(self):
        df = build_oligos(CDS, CUTS, prefix="AAAAUGUGG")
        assert df.iloc[0]["order_fragment_id"] == "AAAAUGUGG_F1"

    def test_aa_lengths_sum_to_the_protein(self):
        df = build_oligos(CDS, CUTS)
        assert df["aa_length"].sum() == len(CDS) // 3

    def test_required_columns_present(self):
        df = build_oligos(CDS, CUTS)
        for col in ("fragment_id", "assembly_order", "oligo_sequence_5to3",
                    "oligo_length", "oh5_coding_site_5to3", "oh3_coding_site_5to3",
                    "aa_length"):
            assert col in df.columns
