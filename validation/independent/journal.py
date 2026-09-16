"""Read design rows out of an experiment journal, ignoring its bookkeeping lines.

A journal written by `clippr.experiment` is not a plain list of designs. It opens with a run
identity record, and appends a checkpoint after each row plus one at the end of every run, so
a 200-design corpus is a 403-line file. A checker that reads every line blindly sees those
records as malformed designs -- which is exactly what happened to V5 and V6 the first time
they were pointed at a journal rather than at a flat corpus.

The reserved keys are written out here rather than imported from `clippr.experiment`, so the
independent checkers stay independent of production code. They are a file format, not a
measurement, but restating them costs one line and keeps the claim simple.
"""
from __future__ import annotations

import json
from pathlib import Path

#: Bookkeeping records. A line carrying any of these is not a design.
RESERVED = ("__run_identity__", "__failure__", "__run__")


def load_designs(path: str | Path) -> list[dict]:
    """Every design row in the journal, in file order, with bookkeeping lines dropped."""
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        if any(key in record for key in RESERVED):
            continue
        rows.append(record)
    return rows


def identity(path: str | Path) -> dict | None:
    """The journal's run identity record, if it has one."""
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.strip():
            record = json.loads(line)
            return record.get("__run_identity__")
    return None


def provenance(validator: str, jsonl: str | Path | None = None) -> dict:
    """What a validator result is evidence *about*, recorded inside the result itself.

    A release gate that reads `{"clean": 200, "designs": 200}` learns nothing about which
    corpus those 200 rows came from, which version of the validator produced the verdict, or
    which package source built the designs. A stale result then satisfies a new release
    directory, which is the failure this record prevents.

    `validator` is the checker's own file path; its bytes are hashed, so editing the checker
    invalidates results produced by the previous version rather than silently inheriting
    their verdict.
    """
    import sys

    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root / "src"))
    from clippr.manifest import source_fingerprint

    record = {
        "validator": Path(validator).name,
        "validator_sha256": _sha256(Path(validator)),
        "package_source_sha256": source_fingerprint()["sha256"],
    }
    if jsonl is not None:
        record["corpus_path"] = str(jsonl)
        record["corpus_sha256"] = _sha256(Path(jsonl))
        record["corpus_identity"] = identity(jsonl)
    return record


def _sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest()


def coverage(path: str | Path) -> dict:
    """Designs, journalled failures and completion, counted separately.

    `load_designs` returns only successful rows, so "every row passes" says nothing about
    whether every requested target produced a row. A run that lost fifty targets and passed
    the rest looks identical to a complete one through that function alone. Coverage is
    therefore a separate gate, and a validator that reports a pass rate should report this
    beside it.
    """
    designs, failures, runs = 0, [], 0
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        if "__failure__" in record:
            failures.append(record["__failure__"])
        elif "__run__" in record:
            runs += 1
        elif "__run_identity__" in record:
            continue
        else:
            designs += 1
    return {"designs": designs, "failures": failures,
            "failure_count": len(failures), "checkpoints": runs,
            "complete": not failures}
