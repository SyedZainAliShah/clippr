"""What a run was: its inputs, its software, and what it did or did not check.

A design is only reproducible if the things it depended on are recorded, and the ones that
bite are not the obvious ones. The codon table is the clearest case: it can arrive supplied,
downloaded, cached or bundled, and what reaches DNA Chisel is `complete_table(table,
genetic_code)`, which rewrites U to T and fills omitted codons at frequency zero. For a
table that is already complete and spelled in DNA -- the Kazusa tables are -- that is the
identity, so the two hashes agree. The point of recording both is that you cannot tell
which case you are in without them.

Screening is the other. `design_oneshot(check_offtarget=False)` returns `offtarget=None`,
which is indistinguishable from a clean screen unless the manifest says the check did not
run. A reference genome is named only when one was actually read.

Nothing here is inferred from a design's outputs. If an input was not supplied it is
recorded as unknown rather than guessed.
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

#: Bump when the recorded structure changes in a way a reader must notice.
SCHEMA_VERSION = 1

_DEPENDENCIES = ("numpy", "pandas", "biopython", "dnachisel", "openpyxl",
                 "python_codon_tables")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_table(table) -> str:
    """Hash a codon table by canonical JSON, so key order cannot change the digest."""
    return _sha256_bytes(json.dumps(table, sort_keys=True, separators=(",", ":")).encode())


def _software() -> dict:
    versions = {}
    for name in _DEPENDENCIES:
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            versions[name] = None
    try:
        clippr_version = version("clippr")
    except PackageNotFoundError:            # running from a checkout without an install
        clippr_version = None
    return {
        "clippr": clippr_version,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "dependencies": versions,
        "source": _source_revision(),
        "source_fingerprint": source_fingerprint(),
    }


#: Data files shipped inside the package that determine a design as much as the code
#: does. `parts.json` is the kit index every module selection reads; `scaffold.json`
#: holds the protein constants every CDS is built from. A fingerprint covering only
#: `.py` files is unchanged by an edit to either, so two runs built on different
#: scaffolds would claim the same source identity -- the precise failure this record
#: exists to make impossible.
_SOURCE_SUFFIXES = (".py", ".json")


def source_fingerprint() -> dict:
    """Content digest of every source and data file in the installed package.

    `_source_revision` records git HEAD plus the *names* of dirty files, which cannot
    distinguish two different edits of the same dirty file -- and on an uncommitted
    checkout that is the normal case, not the exception. Hashing the contents identifies
    the source that actually ran, committed or not, installed or checked out.

    Both code and shipped data are in scope: see `_SOURCE_SUFFIXES`. The per-suffix
    counts are reported so a reader can see what was covered rather than assuming.
    """
    root = Path(__file__).resolve().parent
    files = sorted(f for f in root.rglob("*")
                   if f.is_file() and f.suffix in _SOURCE_SUFFIXES
                   and "__pycache__" not in f.parts)
    digest = hashlib.sha256()
    for f in files:
        digest.update(f.relative_to(root).as_posix().encode())
        digest.update(b"")
        digest.update(f.read_bytes())
    counts = {suffix: sum(1 for f in files if f.suffix == suffix)
              for suffix in _SOURCE_SUFFIXES}
    return {"files": len(files), "by_suffix": counts, "sha256": digest.hexdigest()}


def run_identity(codon_table, genetic_code: int, *, matrix: str = "BsaI-HFv2",
                 **extra) -> dict:
    """Everything that would invalidate a resumed journal, in one record.

    Source contents, dependency versions, the ligation table's contents, the scoring
    formula's version and the optimiser's effective codon table. Pass the experiment's own
    configuration hash and settings through `extra`.
    """
    from .overhangs import SCORER_VERSION, matrix_fingerprint

    return {
        "source_sha256": source_fingerprint()["sha256"],
        "dependencies": _software()["dependencies"],
        "matrix": matrix,
        "matrix_fingerprint": matrix_fingerprint(matrix),
        "scorer_version": SCORER_VERSION,
        "codon_table_effective_sha256": codon_table_record(
            codon_table, genetic_code)["effective_sha256"],
        "genetic_code": genetic_code,
        **extra,
    }


def _source_revision() -> dict:
    """Git revision and dirty state, when the code is a checkout and git is available.

    An installed wheel has no repository, so this is `{"available": false}` there rather
    than a fabricated identity. Never raises: provenance that crashes a design run is worse
    than provenance that admits it is missing.
    """
    root = Path(__file__).resolve().parents[2]
    if not (root / "pyproject.toml").is_file():
        return {"available": False, "reason": "not a source checkout"}
    try:
        run = lambda *a: subprocess.run(  # noqa: E731 - one-line local, not an API
            ["git", "-C", str(root), *a], capture_output=True, text=True, timeout=10)
        head = run("rev-parse", "HEAD")
        if head.returncode != 0:
            return {"available": False, "reason": "git rev-parse failed"}
        dirty = run("status", "--porcelain")
        modified = [ln[3:] for ln in dirty.stdout.splitlines() if ln.strip()]
        return {"available": True, "revision": head.stdout.strip(),
                "dirty": bool(modified), "modified_files": sorted(modified)}
    except (OSError, subprocess.SubprocessError) as exc:
        return {"available": False, "reason": f"{type(exc).__name__}: {exc}"}


def codon_table_record(table, genetic_code: int, *, source: str | Path | None = None,
                       taxid: int | None = None, compartment: str | None = None) -> dict:
    """Identify the codon table a run used, as supplied and as the optimiser saw it.

    `source` is the file the table came from, if any; its bytes are hashed so a run can be
    tied to an artefact on disk. `effective_sha256` covers `complete_table(...)`, which is
    what actually reached DNA Chisel.
    """
    from .codons import complete_table

    # Two hashes of the same input, because they answer different questions. `sha256` ties
    # the run to a file on disk; `canonical_sha256` is serialisation-independent and is the
    # only one comparable with `effective_sha256`. Without the second, a reader sees a file
    # hash differing from the effective hash and concludes the table was transformed, when
    # the difference is only JSON formatting.
    raw: dict = {"path": None, "sha256": None,
                 "canonical_sha256": _sha256_table(table)}
    if source is not None:
        path = Path(source)
        raw["path"] = str(path)
        if path.is_file():
            raw["sha256"] = _sha256_bytes(path.read_bytes())

    return {
        "taxid": taxid,
        "compartment": compartment,
        "genetic_code": genetic_code,
        "supplied": raw,
        "effective_sha256": _sha256_table(complete_table(table, genetic_code)),
        "effective_note": ("sha256 of complete_table(table, genetic_code), the optimiser's "
                           "input. Equal to supplied.canonical_sha256 when the table was "
                           "already complete and spelled in DNA, which is the usual case"),
    }


def _screening(designs: list[dict]) -> dict:
    """Whether host screening ran, and against what.

    Three states, not two. A screen that was *attempted and failed* -- no network, no cached
    genome -- carries a result dictionary like a successful one, so treating any dictionary
    as evidence of assessment reported `assessed: true` beside the status "not assessed",
    which is self-contradictory output about the one thing this record exists to state.
    Only a scan that actually completed counts, and only its reference is named.
    """
    def state(d):
        off = d.get("offtarget")
        if off is None:
            return "disabled"
        return "unavailable" if off.get("verdict") == "not checked" else "completed"

    by_state = {s: [d for d in designs if state(d) == s]
                for s in ("completed", "unavailable", "disabled")}
    completed = by_state["completed"]
    references = sorted({d["offtarget"].get("reference") for d in completed
                         if d["offtarget"].get("reference")})
    record = {
        "assessed": bool(completed),
        "designs_screened": len(completed),
        "designs_unavailable": len(by_state["unavailable"]),
        "designs_disabled": len(by_state["disabled"]),
        "status": sorted({d["audit"].offtarget_status for d in designs if d.get("audit")}),
    }
    errors = sorted({d["offtarget"].get("error") for d in by_state["unavailable"]
                     if d["offtarget"].get("error")})
    if errors:
        record["unavailable_errors"] = errors
    if references:
        record["reference_genome"] = [_genome_record(a) for a in references]
    return record


def _genome_record(accession: str) -> dict:
    from .offtarget import CACHE
    cached = CACHE / f"{accession}.gb"
    return {"accession": accession,
            "sha256": _sha256_bytes(cached.read_bytes()) if cached.is_file() else None}


#: Settings the manifest reports once for the whole batch. Every design must agree on them.
_BATCH_SETTINGS = ("genetic_code", "codon_table", "enzyme_profile_effective")


def _require_homogeneous(designs: list[dict]) -> None:
    """Refuse to describe a heterogeneous batch with one set of settings.

    The manifest states the codon table, genetic code and enzyme profile once, reading them
    off the first design. If a later design used different ones, that record is not merely
    incomplete -- it is false about every design after the first, and nothing downstream
    could detect it. Naming the disagreement is the only honest option; silently reporting
    the first design's settings is not.
    """
    for setting in _BATCH_SETTINGS:
        seen = {_canonical(d.get(setting)) for d in designs}
        if len(seen) > 1:
            offenders = sorted({str(d.get("target_rna", "?")) for d in designs
                                if _canonical(d.get(setting))
                                != _canonical(designs[0].get(setting))})
            raise ValueError(
                f"designs in this batch disagree on {setting!r}, which the manifest reports "
                f"once for the whole run; {len(seen)} distinct values, first differing at "
                f"{offenders[:3]}. Build one manifest per homogeneous batch.")


def _canonical(value) -> str:
    """A serialisation-independent key, so dict ordering never reads as a disagreement."""
    try:
        return json.dumps(value, sort_keys=True, default=str)
    except TypeError:
        return repr(value)


def build(designs, *, codon_table_source: str | Path | None = None,
          taxid: int | None = None, compartment: str | None = None,
          route: str = "synthesis", completion: str = "complete",
          failures: list | None = None) -> dict:
    """Assemble the run manifest for one or more `design_oneshot` results.

    `completion` is the caller's honest statement about coverage -- "complete", "partial",
    "timed out". A run that lost targets must say so here and list them in `failures`,
    because a manifest that only describes successes silently shrinks the denominator.
    """
    designs = [designs] if isinstance(designs, dict) else list(designs)
    if not designs:
        raise ValueError("a manifest needs at least one design")

    first = designs[0]
    _require_homogeneous(designs)
    stage_seconds: dict[str, float] = {}
    for d in designs:
        for stage, seconds in d.get("timings", {}).items():
            stage_seconds[stage] = stage_seconds.get(stage, 0.0) + seconds

    return {
        "schema_version": SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "route": route,
        "completion": completion,
        "software": _software(),
        "inputs": {
            "codon_table": codon_table_record(
                first["codon_table"], first["genetic_code"],
                source=codon_table_source, taxid=taxid, compartment=compartment),
        },
        "configuration": {
            "enzyme_profile": list(first["enzyme_profile_effective"]),
            "n_designs": len(designs),
        },
        "screening": _screening(designs),
        "designs": [{
            "target": d["target_rna"],
            "architecture": d["architecture"],
            "cds_sha256": _sha256_bytes(d["cds"].encode()),
            "qc": d["qc"]["status"],
            "fidelity": d["fidelity"],
            "timings": d.get("timings", {}),
        } for d in designs],
        "stage_seconds": stage_seconds,
        "failures": list(failures or []),
    }


def write(manifest: dict, path: str | Path) -> Path:
    """Write a manifest as sorted, indented JSON so two runs diff cleanly."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out
