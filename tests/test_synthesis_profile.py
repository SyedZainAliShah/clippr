"""One versioned synthesis policy, and the two failures it exists to prevent.

**Failure one: the same threshold in two places.** `recoding.GC_BAND` was the validator's band
while `codons.optimize_cds` kept its own `gc_bounds` default that `recode_inventory` never
passed. Changing one moved which candidates were *accepted* without moving which were
*searched for*, so the obvious repair -- "widening the band is one line" -- produced 0.657195
where a genuinely broadened policy reaches 0.819921. The regression for that is
`test_a_profile_reaches_the_solver_not_just_the_validator`, and it asserts both numbers.

**Failure two: an advisory rule enforced as a rejection.** A vendor design *target* refused as
a hard limit cost about 0.16 CAI, and a soft warning silently read as vendor approval would be
the same error inverted. Hard and target findings therefore never share a return value.

The strict profile is the default and must reproduce today's shipped results exactly; several
tests below pin that, because a drifting default would silently rewrite published figures.
"""
from __future__ import annotations

import json
import statistics
from pathlib import Path

import pytest

from clippr import substrates as sub
from clippr import synthesis_profile as sp

ROOT = Path(__file__).resolve().parents[1]
TABLE_S1 = ROOT / "data" / "grasp_supp" / "Table S1.xlsx"
CODON_TABLE = ROOT / "data" / "codon_tables" / "kazusa_3055.json"

#: A B-module insert whose substrate is clean under every shipped profile.
CLEAN = "ACTC" + "ATGCAC" * 12 + "AAGA"
#: The same shape, carrying a run of five -- a breach of the homopolymer rule only.
RUN_OF_FIVE = "ACTC" + "ATGCAC" * 6 + "A" * 5 + "ATGCAC" * 6 + "AAGA"


class TestTheProfileObject:
    def test_the_default_is_the_strict_legacy_policy(self):
        """Changing this line changes every published number, so it is pinned."""
        assert sp.DEFAULT_PROFILE is sp.STRICT_LEGACY
        assert sp.STRICT_LEGACY.local_gc == (0.35, 0.65)
        assert sp.STRICT_LEGACY.window == 50
        assert sp.STRICT_LEGACY.max_homopolymer == 4
        assert sp.STRICT_LEGACY.global_gc is None
        assert all(sp.STRICT_LEGACY.is_hard(r)
                   for r in ("local_gc", "homopolymer", "forbidden_sites"))

    def test_the_strict_profile_does_not_claim_a_source_it_does_not_have(self):
        """Its band matches Twist's published guidance. That is not where it came from.

        Back-dating a citation would be the same class of error as every other in this
        project's history, so the field says so in words.
        """
        assert "unrecorded" in sp.STRICT_LEGACY.source.lower()
        assert "twist" not in sp.STRICT_LEGACY.source.lower()

    def test_a_sourced_profile_carries_its_source_product_and_date(self):
        twist = sp.TWIST_CODON_GUIDANCE
        assert twist.source.startswith("https://")
        assert twist.read_on == "2026-09-16"
        assert "codon-optimisation" in twist.product or "codon" in twist.product.lower()
        assert twist.unresolved, "a profile with no unresolved caveats is claiming too much"

    def test_the_version_changes_when_a_verdict_changing_field_does(self):
        assert sp.STRICT_LEGACY.version != sp.BROAD_EXPERIMENTAL.version
        assert sp.STRICT_LEGACY.version == sp.STRICT_LEGACY.version
        widened = sp.STRICT_LEGACY.narrowed(local_gc=(0.30, 0.70))
        assert widened.version != sp.STRICT_LEGACY.version

    def test_a_profile_cannot_be_mutated_after_a_result_is_bound_to_it(self):
        with pytest.raises(Exception):
            sp.STRICT_LEGACY.local_gc = (0.1, 0.9)          # type: ignore[misc]

    def test_an_unlisted_rule_is_hard_because_silence_must_not_relax(self):
        quiet = sp.SynthesisProfile(name="quiet", source="test", product="test",
                                    read_on=None, local_gc=(0.35, 0.65), window=50,
                                    max_homopolymer=4)
        assert quiet.is_hard("local_gc") and quiet.is_hard("anything_at_all")

    def test_a_nonsense_enforcement_is_refused(self):
        with pytest.raises(ValueError, match="enforcement"):
            sp.SynthesisProfile(name="bad", source="t", product="t", read_on=None,
                                local_gc=(0.35, 0.65), window=50, max_homopolymer=4,
                                enforcement={"local_gc": "maybe"})

    def test_an_inverted_band_is_refused(self):
        with pytest.raises(ValueError, match="band"):
            sp.SynthesisProfile(name="bad", source="t", product="t", read_on=None,
                                local_gc=(0.8, 0.2), window=50, max_homopolymer=4)

    def test_an_unknown_name_is_refused_rather_than_defaulted(self):
        with pytest.raises(ValueError, match="unknown synthesis profile"):
            sp.resolve("whatever-sounds-plausible")
        assert sp.resolve(None) is sp.DEFAULT_PROFILE
        assert sp.resolve("broad-experimental") is sp.BROAD_EXPERIMENTAL


class TestHardVersusTarget:
    def test_the_same_breach_is_hard_or_advisory_by_profile(self):
        """One sequence, one rule, two policies, two different answers. That is the point."""
        strict = sub.substrate_findings(RUN_OF_FIVE, "B", sp.STRICT_LEGACY)
        advisory = sub.substrate_findings(RUN_OF_FIVE, "B", sp.TWIST_CODON_GUIDANCE)
        assert strict["hard"] and not strict["target"]
        assert advisory["target"] and not advisory["hard"]
        assert "homopolymer" in strict["hard"][0]
        assert "homopolymer" in advisory["target"][0]

    def test_a_target_breach_never_appears_as_a_hard_problem(self):
        """A soft warning silently blocking an export is the failure in one direction."""
        assert sub.substrate_problems(RUN_OF_FIVE, "B", sp.TWIST_CODON_GUIDANCE) == []
        assert sub.substrate_warnings(RUN_OF_FIVE, "B", sp.TWIST_CODON_GUIDANCE)

    def test_a_hard_breach_never_appears_only_as_a_warning(self):
        """And silently becoming 'vendor approved' is the failure in the other."""
        assert sub.substrate_problems(RUN_OF_FIVE, "B", sp.STRICT_LEGACY)
        assert sub.substrate_warnings(RUN_OF_FIVE, "B", sp.STRICT_LEGACY) == []

    def test_no_profile_can_make_a_forbidden_site_advisory(self):
        """An enzyme cutting where it should not is never a matter of preference."""
        lenient = sp.STRICT_LEGACY.narrowed(
            name="tries-to-soften-sites",
            enforcement={"local_gc": sp.TARGET, "homopolymer": sp.TARGET,
                         "forbidden_sites": sp.TARGET})
        internal_bsai = "ACTC" + "ATGCGT" * 3 + "GGTCTC" + "ATGCGT" * 5 + "AAGA"
        problems = sub.substrate_problems(internal_bsai, "B", lenient)
        assert any("BsaI" in p for p in problems), problems

    def test_a_clean_substrate_is_clean_under_every_shipped_profile(self):
        for profile in sp.PROFILES.values():
            found = sub.substrate_findings(CLEAN, "B", profile)
            assert found == {"hard": [], "target": []}, (profile.name, found)


class TestGcReporting:
    def test_the_global_band_is_only_checked_when_a_profile_declares_one(self):
        """The strict profile has no global band; Twist's guidance pairs one with the window."""
        gc_rich = "ACTC" + "GCCGCC" * 12 + "AAGA"
        strict = sub.substrate_findings(gc_rich, "B", sp.STRICT_LEGACY)
        twist = sub.substrate_findings(gc_rich, "B", sp.TWIST_CODON_GUIDANCE)
        assert not any("global GC" in m for m in strict["hard"] + strict["target"])
        assert any("global GC" in m for m in twist["target"])

    def test_the_worst_window_is_reported_and_not_just_the_first(self):
        """A reader takes the first number for the worst; it usually is not.

        The audit flagged this: every reference sequence was reported at 0.660 because the
        scan stopped there, which names a violation without saying how bad the sequence is.
        """
        mild_then_severe = ("ACTC" + "ATGCAC" * 9 + "GCGCGC" * 12 + "ATGCAC" * 2 + "AAGA")
        found = sub.substrate_findings(mild_then_severe, "B", sp.STRICT_LEGACY)
        message = next(m for m in found["hard"] if m.startswith("GC "))
        assert "worst" in message, message


@pytest.mark.skipif(not TABLE_S1.is_file() or not CODON_TABLE.is_file(),
                    reason="the deposited kit and codon table are needed")
class TestTheProfileReachesEveryStage:
    @staticmethod
    def recoded(profile):
        from clippr import inventories as inv
        from clippr.codons import complete_table
        from clippr.objectives import codon_adaptation
        from clippr.recoding import recode_inventory

        raw = json.loads(CODON_TABLE.read_text(encoding="utf-8"))
        table = complete_table(raw, 1)
        report = recode_inventory(inv.load_deposited(TABLE_S1), raw,
                                  label="profile-regression", seeds=4, profile=profile)
        cais = [codon_adaptation(m.dna[slice(*m.coding_interval)], table)["cai"]
                for m in report.inventory.modules.values()]
        return report, statistics.mean(cais)

    def test_a_profile_reaches_the_solver_not_just_the_validator(self):
        """**The** regression for the "one line" defect.

        Patching only the validator's constant gives 0.657195 because the solver keeps
        searching the narrow space. A profile that reaches both reaches 0.819921. Asserting
        the broad figure is what proves the solver actually received the policy -- the
        validator alone cannot produce it.
        """
        _, strict_mean = self.recoded(sp.STRICT_LEGACY)
        _, broad_mean = self.recoded(sp.BROAD_EXPERIMENTAL)
        assert strict_mean == pytest.approx(0.618830, abs=1e-5)
        assert broad_mean == pytest.approx(0.819921, abs=1e-5)
        assert broad_mean > 0.70, "the solver did not receive the broadened policy"

    def test_the_default_reproduces_the_shipped_result(self):
        """No argument at all must equal the strict profile, or published numbers moved."""
        _, explicit = self.recoded(sp.STRICT_LEGACY)
        _, implicit = self.recoded(None)
        assert implicit == pytest.approx(explicit, abs=1e-12)
        assert implicit == pytest.approx(0.618830, abs=1e-5)

    def test_the_report_records_the_policy_that_produced_it(self):
        """A CAI figure without its policy is not reproducible; the two differ by 0.20."""
        report, _ = self.recoded(sp.BROAD_EXPERIMENTAL)
        recorded = report.as_dict()["profile"]
        assert recorded["name"] == "broad-experimental"
        assert recorded["version"] == sp.BROAD_EXPERIMENTAL.version

    def test_every_delivered_module_is_feasible_under_its_own_profile(self):
        """Whatever policy produced an inventory, its output must satisfy that policy."""
        for profile in (sp.STRICT_LEGACY, sp.BROAD_EXPERIMENTAL):
            report, _ = self.recoded(profile)
            offenders = {m: sub.substrate_problems(r.dna, r.block, profile)
                         for m, r in report.inventory.modules.items()
                         if sub.substrate_problems(r.dna, r.block, profile)}
            assert offenders == {}, (profile.name, offenders)
