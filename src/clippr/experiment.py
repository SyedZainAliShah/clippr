"""Long runs that survive being interrupted, and stop when told to.

A 200-target corpus takes hours. Three things go wrong at that length, and all three have
happened here: the machine is interrupted and the work is lost; a resumed run silently mixes
rows from two different configurations; and a run that timed out is reported as if it had
finished.

`run` addresses exactly those. Rows are appended and flushed as they complete, so an
interrupted run keeps what it had. The identity of the run -- source, configuration, inputs --
is written as the first line and checked on resume, so rows produced under a different
configuration are refused rather than merged. A run stopped by its budget returns
`status="timed out"` with the rows it finished, never `"complete"`.

**The budget is wall time, not CPU time.** Stage timings elsewhere in this package are also
wall, and a user waiting on a design waits in wall time. It is checked before starting each
item, so a single item that overruns the budget still finishes rather than being killed
mid-write; the budget bounds how long the run keeps starting new work, not the longest
possible run.

**Serial by design.** `codons.optimize_cds` seeds the global `random` and `numpy.random`
state, so two designs running concurrently in one process would interleave their draws and
neither would be reproducible from its seed. Parallelism here needs separate processes and a
demonstration that seeded outputs are unchanged; it is not a free speedup.
"""
from __future__ import annotations

import json
import time
import uuid
from collections.abc import Callable, Iterable
from contextlib import contextmanager
from pathlib import Path

#: Marks the first line of a record file. A row can never collide with it: rows are written
#: by the caller's worker, and this key names the run rather than any one item.
IDENTITY_KEY = "__run_identity__"
#: Records that are not results. Failures are journalled so a resumed run still knows what
#: went wrong before it, and each run appends its elapsed time so a budget spans resumes
#: rather than restarting its clock every time the process does.
FAILURE_KEY = "__failure__"
RUN_KEY = "__run__"
#: One record per line; written explicitly so an escape cannot be lost in transit.
NEWLINE = chr(10)


@contextmanager
def stage(timings: dict[str, float], name: str):
    """Accumulate wall time for one pipeline stage into `timings`.

    Wall, not CPU, for the same reason the run budget is: these stages call out to DNA
    Chisel and the filesystem, and the number a user waits through is the wall clock.
    Accumulates rather than assigns, so a stage entered once per retry reports its total.
    """
    start = time.perf_counter()
    try:
        yield
    finally:
        timings[name] = timings.get(name, 0.0) + time.perf_counter() - start


class IncompatibleResume(Exception):
    """Existing rows were produced under a different identity than the one requested."""


def _repair_tail(path: Path) -> str | None:
    """Move an unparseable final line out of the journal before anything is appended.

    The reader already ignored a truncated tail, but `run` then opened the same file in
    append mode and the next record fused onto the broken bytes -- producing a line like
    `{"key": "3", "va{"key": "3", "v": 6}` that no later resume could parse. Ignoring bad
    bytes in memory is not the same as not writing after them.

    The quarantined text is kept beside the journal rather than deleted; it is the only
    evidence of what the interrupted run was part-way through writing.
    """
    if not path.is_file() or path.stat().st_size == 0:
        return None
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    if not lines:
        return None
    try:
        json.loads(lines[-1])
    except json.JSONDecodeError:
        pass
    else:
        # The last record parses, but a file that does not end in a newline still fuses the
        # next append onto it. Parsing the tail is not the same as being able to write after
        # it -- the first version of this repair checked only the former.
        if not text.endswith(NEWLINE):
            path.write_text(text + NEWLINE, encoding="utf-8")
        return None
    # A fixed `.partial` name meant the second interruption destroyed the first one's
    # evidence -- and a journal interrupted twice is exactly the case where that evidence
    # matters most. Number the quarantines instead, never reusing one.
    stem = path.with_suffix(path.suffix + ".partial")
    quarantine, n = stem, 1
    while quarantine.exists():
        quarantine = stem.with_name(f"{stem.name}.{n}")
        n += 1
    quarantine.write_text(lines[-1], encoding="utf-8")
    path.write_text(NEWLINE.join(lines[:-1]) + NEWLINE, encoding="utf-8")
    return str(quarantine)


def _read_existing(path: Path, identity: dict) -> tuple[list[dict], list[dict], float]:
    """Completed rows, journalled failures and seconds already spent, under this identity.

    A trailing line that is not valid JSON is a run interrupted mid-write. It is dropped
    with a note rather than taken as corruption of the whole journal, because the records
    before it are complete and discarding them would lose hours of valid work.
    """
    if not path.is_file() or path.stat().st_size == 0:
        return [], [], 0.0
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    header = json.loads(lines[0])
    if IDENTITY_KEY not in header:
        raise IncompatibleResume(
            f"{path} has no identity header; it was not written by this runner")
    found = header[IDENTITY_KEY]
    if found != identity:
        differing = sorted({k for k in set(found) | set(identity)
                            if found.get(k) != identity.get(k)})
        raise IncompatibleResume(
            f"{path} was produced under a different identity; fields differing: "
            f"{', '.join(differing) or '(ordering only)'}. Delete it to start over rather "
            f"than mixing runs.")

    rows, failures, runs = [], [], {}
    for i, line in enumerate(lines[1:], start=2):
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            if i == len(lines):                      # interrupted trailing write
                break
            raise
        if FAILURE_KEY in record:
            failures.append(record[FAILURE_KEY])
        elif RUN_KEY in record:
            rid = record[RUN_KEY].get("run_id", "legacy")
            runs[rid] = max(runs.get(rid, 0.0),
                            float(record[RUN_KEY].get("wall_seconds", 0.0)))
        else:
            rows.append(record)
    # Each run checkpoints its own elapsed time repeatedly; take the last value per run and
    # sum across runs, so a killed run still contributes the time it had already logged.
    return rows, failures, sum(runs.values())


def run(items: Iterable, worker: Callable[[object], dict], *, path: str | Path,
        identity: dict, key: Callable[[object], str] = str,
        wall_seconds: float | None = None,
        on_progress: Callable[[int, int, str], None] | None = None) -> dict:
    """Apply `worker` to each item, recording rows as they complete.

    Returns the rows plus an honest account of coverage: how many were requested, how many
    were produced, which failed and why, and whether the budget stopped it. `requested` is
    the full input count, so a failure cannot shrink the denominator.

    `identity` should name everything that would invalidate a resume -- source revision,
    configuration, input hashes. `manifest.build(...)` produces a suitable record; pass a
    subset of it rather than inventing a second vocabulary.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Two items with the same key are one unit of work. Left alone the second would run
    # again and append a duplicate row, because the completed set was only consulted, never
    # updated, inside the loop.
    items, seen_keys = list(items), set()
    deduplicated = []
    for item in items:
        k = key(item)
        if k in seen_keys:
            continue
        seen_keys.add(k)
        deduplicated.append(item)
    duplicates_dropped = len(items) - len(deduplicated)
    items = deduplicated

    quarantined = _repair_tail(path)
    done, failures, already_spent = _read_existing(path, identity)
    resumed = len(done)
    have = {r["key"] for r in done}
    if not path.is_file() or path.stat().st_size == 0:
        with path.open("w", encoding="utf-8") as fh:
            fh.write(json.dumps({IDENTITY_KEY: identity}, sort_keys=True) + "\n")

    started = time.perf_counter()
    status = "complete"
    prior_failures = len(failures)
    run_id = uuid.uuid4().hex

    def spent() -> float:
        """Seconds across every run of this journal, not just this process."""
        return already_spent + (time.perf_counter() - started)

    clean_exit = False
    persisted = 0.0
    with path.open("a", encoding="utf-8") as fh:
        def checkpoint(final: bool = False) -> None:
            """Persist elapsed time. Written after every item and again on any exit.

            Checkpointing only between items loses the time spent inside the item that was
            interrupted: a run killed 30 s into its second item recorded 4 s, and the resume
            then started fresh work believing almost its whole budget remained. The `finally`
            below covers every catchable exit, KeyboardInterrupt included.

            `clean_exit` separates a total that is exact from one that is a **lower bound**.
            An abruptly terminated process (SIGKILL, power loss) cannot write anything, so its
            last checkpoint under-reports by up to one item's duration; a reader must treat
            `wall_seconds_all_runs` as a floor unless `clean_exit` is true.
            """
            nonlocal persisted
            persisted = time.perf_counter() - started
            fh.write(json.dumps({RUN_KEY: {
                "run_id": run_id,
                "wall_seconds": persisted,
                "clean_exit": bool(final and clean_exit),
                "basis": "exact" if (final and clean_exit) else "lower bound",
            }}) + NEWLINE)
            fh.flush()

        try:
            for i, item in enumerate(items, 1):
                k = key(item)
                if k in have:
                    continue
                if wall_seconds is not None and spent() >= wall_seconds:
                    status = "timed out"
                    break
                try:
                    row = {"key": k, **worker(item)}
                except Exception as exc:              # noqa: BLE001 - recorded, not hidden
                    failure = {"key": k, "error": f"{type(exc).__name__}: {exc}"}
                    failures.append(failure)
                    fh.write(json.dumps({FAILURE_KEY: failure}) + NEWLINE)
                    checkpoint()
                    continue
                done.append(row)
                fh.write(json.dumps(row, default=str) + NEWLINE)
                checkpoint()
                if on_progress:
                    on_progress(i, len(items), k)
            clean_exit = True
        finally:
            checkpoint(final=True)

    # The total reported is the total persisted, not a fresh reading taken after the file
    # closed. They differed by the cost of closing, so every resume inherited a slightly
    # smaller budget than the previous run had actually reported spending -- a small leak,
    # but one that compounds across resumes and makes the returned and recorded totals
    # disagree for no reason a reader could discover.
    elapsed_total = already_spent + persisted
    return {
        "rows": done,
        # `requested` counts this call's items; `completed` counts every row in the file,
        # including any resumed from an earlier call. Reporting only those two invites
        # reading "requested 24, completed 72" as a coverage claim, so the split is explicit.
        "requested": len(items),
        "completed": len(done),
        "resumed": resumed,
        "completed_this_run": len(done) - resumed,
        "failures": failures,
        "status": status,
        "new_failures": len(failures) - prior_failures,
        "duplicates_dropped": duplicates_dropped,
        "wall_seconds": elapsed_total - already_spent,
        "wall_seconds_all_runs": elapsed_total,
        "budget_wall_seconds": wall_seconds,
        # The budget is checked between items, so an item already running when it expires
        # still finishes. Saying "complete" without saying "and over budget" would let a run
        # that blew its cap read as one that respected it.
        "budget_exceeded": bool(wall_seconds is not None and elapsed_total > wall_seconds),
        "path": str(path),
        "quarantined_tail": quarantined,
    }
