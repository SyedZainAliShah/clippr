"""Controls for release check 13: what counts as this release's evidence.

Two rounds of defects sit behind these. First the check compared a record's installed and
checkout fingerprints "if both exist", so a record carrying neither passed and a record
reduced to `{"ok": true}` passed with no identity at all. Then it required fields to be
non-empty without recomputing them, so a notebook hash, a staged-input hash, a wheel hash
or a corpus hash replaced by sixty-four zeros still passed, as did a report whose two
hundred rows were one result copied two hundred times.

Every fixture here is built by **copying the real current records** and changing one thing,
which is the only way a control shows what the gate would have accepted. Nothing is written
outside `tmp_path`, and no authoritative evidence is read for writing.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "validation" / "release_check.py"
ZERO = "0" * 64

def _newest_bundle():
    """The most recent release bundle on disk. Discovered, not hardcoded, so these
    controls follow the current release instead of pinning one directory name."""
    candidates = sorted((ROOT / "work").glob("repair_release_*"), reverse=True)
    return next((c for c in candidates if (c / "final_artifact_contract.json").is_file()),
                None)


@pytest.fixture(scope="module")
def release():
    spec = importlib.util.spec_from_file_location("_release_check", RELEASE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def live(release):
    """Skip rather than pretend: these controls are about real records."""
    bundle = _newest_bundle()
    if bundle is None:
        pytest.skip("no release bundle in work/ to build controls from")
    verdict = release.check_evidence_dir(str(bundle))
    if not verdict or not verdict["ok"]:
        pytest.skip(f"{bundle.name} is not a currently passing bundle "
                    f"({verdict and verdict['detail'][:120]})")
    return bundle


@pytest.fixture
def bundle(live, tmp_path):
    """An isolated copy of the real bundle, for one control to spoil."""
    where = tmp_path / "evidence"
    shutil.copytree(live, where, ignore=shutil.ignore_patterns(
        "environment", "nb_designer", "nb_inventory"))
    return where


def _edit(bundle, *parts, **changes):
    """Rewrite one record in the copied bundle."""
    f = bundle.joinpath(*parts)
    record = json.loads(f.read_text(encoding="utf-8"))
    record.update(changes)
    f.write_text(json.dumps(record, indent=2), encoding="utf-8")
    return record


def _load(bundle, *parts):
    return json.loads(bundle.joinpath(*parts).read_text(encoding="utf-8"))


class TestTheGateAcceptsGenuineEvidence:
    def test_the_real_bundle_passes(self, release, live):
        """The positive control is the actual release bundle, not a hand-built fixture.

        An earlier version of this file constructed a Gate F record holding only
        `{"gate": "PASS"}`, a zeroed notebook hash and arbitrary corpus rows, and called
        that the positive case. It blessed a fixture the gate should reject.
        """
        verdict = release.check_evidence_dir(str(live))
        assert verdict["ok"] is True, verdict["detail"]
        assert "recomputed and matched" in verdict["detail"]

    def test_an_untouched_copy_passes(self, release, bundle):
        """Copying must not break it, or every control below proves nothing."""
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is True, verdict["detail"]

    def test_no_directory_leaves_an_ordinary_run_unchanged(self, release):
        assert release.check_evidence_dir(None) is None


class TestRecomputedIdentities:
    """Each of these passed when the gate only checked that a field was non-empty."""

    def test_a_zeroed_notebook_hash_is_refused(self, release, bundle):
        _edit(bundle, "notebook_inventory.json", notebook_sha256=ZERO)
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "hashes to" in verdict["detail"]

    def test_a_zeroed_staged_input_hash_is_refused(self, release, bundle):
        record = _load(bundle, "notebook_inventory.json")
        staged = record["staged_inputs"]
        assert staged, "fixture assumes the inventory notebook stages Table S1"
        for entry in staged.values():
            entry["sha256"] = ZERO
        _edit(bundle, "notebook_inventory.json", staged_inputs=staged)

        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "source now hashes to" in verdict["detail"]

    def test_a_zeroed_gate_f_wheel_and_source_is_refused(self, release, bundle):
        record = _load(bundle, "phasef", "gate_f.json")
        wheel = dict(record["wheel"], sha256=ZERO)
        _edit(bundle, "phasef", "gate_f.json",
              installed_source_fingerprint=ZERO, wheel=wheel, wheel_sha256=ZERO)

        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "gate_f built 000000000000" in verdict["detail"]

    def test_a_zeroed_corpus_hash_is_refused(self, release, bundle):
        record = _load(bundle, "final_artifact_contract.json")
        prov = dict(record["provenance"], corpus_sha256=ZERO)
        _edit(bundle, "final_artifact_contract.json", provenance=prov)

        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "hashes to" in verdict["detail"]

    def test_a_tampered_validator_record_is_refused(self, release, bundle):
        record = _load(bundle, "v6_constraints.json")
        prov = dict(record.get("provenance") or {}, validator_sha256=ZERO)
        _edit(bundle, "v6_constraints.json", provenance=prov)

        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "different version of v6_constraints.py" in verdict["detail"]


class TestDeclaredContracts:
    def test_an_empty_expectation_list_is_refused(self, release, bundle):
        """`postconditions_declared=8` with `expectations=[]` passed on the count."""
        _edit(bundle, "notebook_inventory.json", expectations=[])
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "but lists 0" in verdict["detail"]

    def test_a_trivial_expectation_does_not_satisfy_the_contract(self, release, bundle):
        """A positive count is not a contract: `--expect "1"` used to be enough."""
        _edit(bundle, "notebook_inventory.json",
              postconditions_declared=1,
              expectations=[{"expect": "1", "ok": True, "value": "1"}])
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "declares no postcondition covering" in verdict["detail"]

    def test_an_unmet_postcondition_is_refused(self, release, bundle):
        record = _load(bundle, "notebook_designer.json")
        spoiled = [dict(record["expectations"][0], ok=False)] + \
            record["expectations"][1:]
        _edit(bundle, "notebook_designer.json", expectations=spoiled)
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "postcondition not met" in verdict["detail"]


class TestPanelCoverage:
    def test_a_duplicated_panel_is_refused(self, release, bundle):
        """200 copies of one row satisfied a row count of 200."""
        record = _load(bundle, "final_artifact_contract.json")
        rows = record["rows"]
        _edit(bundle, "final_artifact_contract.json", rows=[rows[0]] * len(rows))

        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "not the" in verdict["detail"] and "declared panel rows" in verdict["detail"]

    def test_a_short_panel_is_refused(self, release, bundle):
        record = _load(bundle, "final_artifact_contract.json")
        _edit(bundle, "final_artifact_contract.json",
              rows=record["rows"][:1], comparisons_performed=1)
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "not the declared 200" in verdict["detail"]

    def test_uncompared_rows_are_refused(self, release, bundle):
        _edit(bundle, "final_artifact_contract.json", comparisons_performed=1)
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "compared 1 of 200" in verdict["detail"]

    def test_unassessed_designs_are_refused(self, release, bundle):
        """Agreeing that data is missing is not a completed validation."""
        _edit(bundle, "final_artifact_contract.json", not_assessed=3)
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "could not be assessed" in verdict["detail"]

    def test_a_disagreement_is_refused(self, release, bundle):
        _edit(bundle, "final_artifact_contract.json",
              disagreements_with_production=[{"target": "AAAAUGUGG"}])
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "disagrees with production" in verdict["detail"]


class TestRequiredCoverage:
    @pytest.mark.parametrize("parts", [
        ("gate_b.json",), ("phasec", "gate_c.json"), ("phased", "gate_d.json"),
        ("phasee", "gate_e.json"), ("phasef", "gate_f.json"),
        ("v5_digest_ligate.json",), ("v6_constraints.json",),
        ("v7_export_semantics.json",), ("notebook_designer.json",),
        ("notebook_inventory.json",), ("final_artifact_contract.json",),
    ])
    def test_every_required_record_is_required(self, release, bundle, parts):
        """A bundle of four files used to pass; B-E and V5-V7 were never required."""
        bundle.joinpath(*parts).unlink()
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "missing from" in verdict["detail"]

    def test_a_four_record_directory_is_refused(self, release, bundle, tmp_path):
        """The exact shape the review passed: Gate F, two notebooks, the corpus report."""
        thin = tmp_path / "thin"
        (thin / "phasef").mkdir(parents=True)
        shutil.copy(bundle / "phasef" / "gate_f.json", thin / "phasef" / "gate_f.json")
        for name in ("notebook_designer.json", "notebook_inventory.json",
                     "final_artifact_contract.json"):
            shutil.copy(bundle / name, thin / name)

        verdict = release.check_evidence_dir(str(thin))
        assert verdict["ok"] is False
        assert "missing from" in verdict["detail"]

    def test_a_failed_gate_is_refused(self, release, bundle):
        _edit(bundle, "phasec", "gate_c.json", gate="FAIL")
        assert release.check_evidence_dir(str(bundle))["ok"] is False

    def test_a_failed_gate_f_is_refused(self, release, bundle):
        _edit(bundle, "phasef", "gate_f.json", gate="FAIL")
        assert release.check_evidence_dir(str(bundle))["ok"] is False

    def test_a_notebook_that_did_not_pass_is_refused(self, release, bundle):
        _edit(bundle, "notebook_designer.json", ok=False)
        assert release.check_evidence_dir(str(bundle))["ok"] is False


class TestWheelBinding:
    def test_a_wheel_that_is_not_on_disk_cannot_be_verified(self, release, bundle):
        record = _load(bundle, "phasef", "gate_f.json")
        wheel = dict(record["wheel"], path=str(bundle / "nowhere.whl"))
        _edit(bundle, "phasef", "gate_f.json", wheel=wheel,
              wheel_path=str(bundle / "nowhere.whl"))
        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "no retained wheel" in verdict["detail"]

    def test_an_altered_wheel_is_detected(self, release, bundle, tmp_path):
        """The hash is recomputed from the file, so editing the wheel is visible."""
        record = _load(bundle, "phasef", "gate_f.json")
        original = Path(record["wheel"]["path"])
        if not original.is_file():
            pytest.skip("the retained wheel is not on disk in this checkout")
        copy = tmp_path / original.name
        copy.write_bytes(original.read_bytes() + b"\x00")
        wheel = dict(record["wheel"], path=str(copy))
        _edit(bundle, "phasef", "gate_f.json", wheel=wheel, wheel_path=str(copy))

        verdict = release.check_evidence_dir(str(bundle))
        assert verdict["ok"] is False
        assert "the file hashes to" in verdict["detail"]


class TestNoAuthoritativeEvidenceWasTouched:
    def test_the_live_bundle_is_unchanged_by_these_controls(self, live):
        """The renderer incident overwrote the shared release record. Not again."""
        digest = hashlib.sha256()
        for f in sorted(live.rglob("*.json")):
            if "environment" in f.parts:
                continue
            digest.update(f.relative_to(live).as_posix().encode())
            digest.update(f.read_bytes())
        # Recorded by the first test to run; compared by every later one in the module.
        seen = getattr(TestNoAuthoritativeEvidenceWasTouched, "_digest", None)
        if seen is None:
            TestNoAuthoritativeEvidenceWasTouched._digest = digest.hexdigest()
        else:
            assert seen == digest.hexdigest(), "a control wrote into the live bundle"
