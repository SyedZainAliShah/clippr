"""The three ways a long run goes wrong: lost work, mixed runs, timeouts called success."""
from __future__ import annotations

import json
from unittest.mock import patch

import pytest

from clippr.experiment import IDENTITY_KEY, IncompatibleResume, run

IDENTITY = {"revision": "abc123", "config": "nine-mer"}


def double(n):
    return {"value": n * 2}


class TestResume:
    def test_completed_rows_are_not_recomputed(self, tmp_path):
        path = tmp_path / "rows.jsonl"
        seen = []

        def counting(n):
            seen.append(n)
            return double(n)

        run([1, 2, 3], counting, path=path, identity=IDENTITY)
        assert seen == [1, 2, 3]
        second = run([1, 2, 3, 4], counting, path=path, identity=IDENTITY)
        assert seen == [1, 2, 3, 4], "a resumed run repeated work it had already done"
        assert second["completed"] == 4

    def test_rows_survive_an_interrupted_run(self, tmp_path):
        path = tmp_path / "rows.jsonl"

        def fails_on_third(n):
            if n == 3:
                raise KeyboardInterrupt("simulated interruption")
            return double(n)

        with pytest.raises(KeyboardInterrupt):
            run([1, 2, 3], fails_on_third, path=path, identity=IDENTITY)
        resumed = run([1, 2, 3], double, path=path, identity=IDENTITY)
        assert resumed["completed"] == 3
        assert [r["value"] for r in resumed["rows"]] == [2, 4, 6]

    def test_resumed_rows_are_counted_apart_from_new_ones(self, tmp_path):
        """`requested` is this call's items; `completed` is every row in the file."""
        path = tmp_path / "rows.jsonl"
        run([1, 2], double, path=path, identity=IDENTITY)
        second = run([3], double, path=path, identity=IDENTITY)
        assert second["requested"] == 1
        assert second["completed"] == 3
        assert second["resumed"] == 2 and second["completed_this_run"] == 1

    def test_a_different_identity_is_refused(self, tmp_path):
        path = tmp_path / "rows.jsonl"
        run([1, 2], double, path=path, identity=IDENTITY)
        with pytest.raises(IncompatibleResume, match="config"):
            run([1, 2], double, path=path, identity={**IDENTITY, "config": "fourteen-mer"})

    def test_a_foreign_file_is_refused(self, tmp_path):
        path = tmp_path / "rows.jsonl"
        path.write_text('{"key": "1", "value": 2}\n', encoding="utf-8")
        with pytest.raises(IncompatibleResume, match="no identity header"):
            run([1], double, path=path, identity=IDENTITY)

    def test_identity_is_recorded_first(self, tmp_path):
        path = tmp_path / "rows.jsonl"
        run([1], double, path=path, identity=IDENTITY)
        first = json.loads(path.read_text(encoding="utf-8").splitlines()[0])
        assert first[IDENTITY_KEY] == IDENTITY


class TestDurability:
    def test_failures_survive_a_resume(self, tmp_path):
        """A failure held only in memory is invisible to the next run."""
        path = tmp_path / "rows.jsonl"

        def fails_on_two(n):
            if n == 2:
                raise ValueError("no design for 2")
            return double(n)

        run([1, 2, 3], fails_on_two, path=path, identity=IDENTITY)
        second = run([4], double, path=path, identity=IDENTITY)
        assert [f["key"] for f in second["failures"]] == ["2"]
        assert second["new_failures"] == 0

    def test_a_duplicate_key_cannot_run_twice(self, tmp_path):
        seen = []

        def counting(n):
            seen.append(n)
            return double(n)

        out = run([1, 1, 2], counting, path=tmp_path / "rows.jsonl", identity=IDENTITY)
        assert seen == [1, 2] and out["duplicates_dropped"] == 1
        assert out["completed"] == 2

    def test_an_interrupted_trailing_write_survives_two_resumes(self, tmp_path):
        """Checking only the first resume's return value missed the durable defect.

        The reader ignored the broken tail, but the next record was appended straight onto
        it -- fusing them into a line no later resume could parse. The file has to be
        inspected, and a second resume attempted, for this to fail when it is broken.
        """
        import json as J
        path = tmp_path / "rows.jsonl"
        run([1, 2], double, path=path, identity=IDENTITY)
        with path.open("a", encoding="utf-8") as fh:
            fh.write('{"key": "3", "val')                      # truncated mid-write

        first = run([3], double, path=path, identity=IDENTITY)
        assert first["quarantined_tail"], "the broken tail was not moved aside"
        second = run([4], double, path=path, identity=IDENTITY)

        assert second["completed"] == 4
        assert sorted(r["key"] for r in second["rows"]) == ["1", "2", "3", "4"]
        for line in path.read_text(encoding="utf-8").splitlines():
            J.loads(line)                                      # every line must parse

    def test_interrupted_time_is_not_lost(self, tmp_path):
        """Elapsed time written only on a clean exit gave a killed run a fresh budget."""
        import time
        path = tmp_path / "rows.jsonl"

        def dies_partway(n):
            time.sleep(0.05)
            if n == 3:
                raise KeyboardInterrupt("simulated kill before the closing record")
            return double(n)

        with pytest.raises(KeyboardInterrupt):
            run(range(10), dies_partway, path=path, identity=IDENTITY, wall_seconds=5)
        resumed = run(range(10), lambda n: (time.sleep(0.05), double(n))[1],
                      path=path, identity=IDENTITY, wall_seconds=5)
        assert resumed["wall_seconds_all_runs"] > resumed["wall_seconds"], (
            "the interrupted run's time was reset rather than carried")


class TestBudget:
    def test_a_timed_out_run_does_not_report_completion(self, tmp_path):
        import time

        def slow(n):
            time.sleep(0.05)
            return double(n)

        out = run(range(100), slow, path=tmp_path / "rows.jsonl", identity=IDENTITY,
                  wall_seconds=0.15)
        assert out["status"] == "timed out"
        assert out["completed"] < out["requested"]
        assert out["rows"], "a timed-out run should keep the work it finished"

    def test_the_budget_spans_resumes(self, tmp_path):
        """Each resume restarting its clock made any cap unreachable by running twice."""
        import time
        path = tmp_path / "rows.jsonl"

        def slow(n):
            time.sleep(0.05)
            return double(n)

        first = run(range(4), slow, path=path, identity=IDENTITY, wall_seconds=0.12)
        second = run(range(20), slow, path=path, identity=IDENTITY, wall_seconds=0.12)
        assert second["wall_seconds_all_runs"] > first["wall_seconds"]
        assert second["status"] == "timed out"

    def test_an_overrun_is_reported_even_when_every_item_finishes(self, tmp_path):
        """The check is between items, so 'complete' must not imply 'within budget'."""
        import time
        out = run([1], lambda n: (time.sleep(0.2), double(n))[1],
                  path=tmp_path / "rows.jsonl", identity=IDENTITY, wall_seconds=0.01)
        assert out["status"] == "complete" and out["budget_exceeded"] is True

    def test_an_unbudgeted_run_completes(self, tmp_path):
        out = run([1, 2], double, path=tmp_path / "rows.jsonl", identity=IDENTITY)
        assert out["status"] == "complete" and out["budget_wall_seconds"] is None


class TestCoverage:
    def test_a_failure_does_not_shrink_the_denominator(self, tmp_path):
        def fails_on_two(n):
            if n == 2:
                raise ValueError("no design for 2")
            return double(n)

        out = run([1, 2, 3], fails_on_two, path=tmp_path / "rows.jsonl", identity=IDENTITY)
        assert out["requested"] == 3 and out["completed"] == 2
        assert out["failures"] == [{"key": "2", "error": "ValueError: no design for 2"}]
        assert out["status"] == "complete", "failures are reported, not treated as a timeout"


class TestBudgetSurvivesInterruption:
    """Time spent inside an interrupted item must not vanish from the budget.

    Checkpointing only *between* items lost whatever the interrupted item had consumed. A run
    killed 30 s into its second item recorded 4 s, and the resume then began fresh work
    believing nearly its whole budget remained.
    """

    def test_a_keyboard_interrupt_still_records_its_time(self, tmp_path):
        path = tmp_path / "clock.jsonl"
        clock = [0.0]

        def worker(n):
            clock[0] += 4.0 if n == 1 else 30.0
            if n == 2:
                raise KeyboardInterrupt
            return {"v": n}

        with patch("clippr.experiment.time.perf_counter", side_effect=lambda: clock[0]):
            with pytest.raises(KeyboardInterrupt):
                run([1, 2], worker, path=path, identity=IDENTITY, wall_seconds=5)
            spent_before_resume = clock[0]

            def quick(n):
                clock[0] += 1.0
                return {"v": n}

            resumed = run([3], quick, path=path, identity=IDENTITY, wall_seconds=5)

        assert resumed["wall_seconds_all_runs"] >= spent_before_resume
        assert resumed["completed_this_run"] == 0, (
            "the budget was already exhausted, so no new item should have started")

    def test_an_interrupted_total_is_labelled_a_lower_bound(self, tmp_path):
        path = tmp_path / "labelled.jsonl"

        def worker(n):
            if n == 2:
                raise KeyboardInterrupt
            return {"v": n}

        with pytest.raises(KeyboardInterrupt):
            run([1, 2], worker, path=path, identity=IDENTITY)

        records = [json.loads(line) for line in path.read_text().splitlines()]
        checkpoints = [r["__run__"] for r in records if "__run__" in r]
        assert checkpoints, "no checkpoint was written"
        assert checkpoints[-1]["clean_exit"] is False
        assert checkpoints[-1]["basis"] == "lower bound"

    def test_a_clean_run_is_labelled_exact(self, tmp_path):
        path = tmp_path / "exact.jsonl"
        run([1, 2], lambda n: {"v": n}, path=path, identity=IDENTITY)
        records = [json.loads(line) for line in path.read_text().splitlines()]
        final = [r["__run__"] for r in records if "__run__" in r][-1]
        assert final["clean_exit"] is True
        assert final["basis"] == "exact"


class TestQuarantinePreservation:
    """A journal interrupted twice must keep both partial writes.

    The quarantine filename was fixed, so the second interruption overwrote the first one's
    evidence -- and a journal interrupted twice is precisely when that evidence matters.
    """

    def test_a_second_interruption_does_not_destroy_the_first_quarantine(self, tmp_path):
        path = tmp_path / "twice.jsonl"
        run([1], lambda n: {"v": n}, path=path, identity=IDENTITY)

        path.write_text(path.read_text() + '{"key": "2", "v', encoding="utf-8")
        first = run([2], lambda n: {"v": n}, path=path, identity=IDENTITY)

        path.write_text(path.read_text() + '{"key": "3", "part', encoding="utf-8")
        second = run([3], lambda n: {"v": n}, path=path, identity=IDENTITY)

        assert first["quarantined_tail"] and second["quarantined_tail"]
        assert first["quarantined_tail"] != second["quarantined_tail"]
        from pathlib import Path as _P
        assert _P(first["quarantined_tail"]).read_text() == '{"key": "2", "v'
        assert _P(second["quarantined_tail"]).read_text() == '{"key": "3", "part'
