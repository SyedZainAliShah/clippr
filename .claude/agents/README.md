# The four roles, and the one thing that is not a model

*Set up 2026-10-16, after a regression reached `main` because the agent that wrote the change
was also the agent that decided it was fine.*

## The problem this solves

Every defect this project has shipped was produced by a competent agent who was sure the work
was correct. Not carelessness — the ordinary failure of being author, reviewer and judge at
the same time. Concretely, in one session:

- a change was pushed after running the notebook tests but not the suite, and six tests were
  broken on `main` for three commits;
- three of thirteen `file:line` citations in a document pointed at the wrong line;
- a shell heredoc mangled an escape and silently corrupted a generator, twice;
- a run was described as in progress when the process list showed it had never started.

None of these would have been caught by an instruction to be careful. Prose asking for
diligence is read by the same mind that already believed it was being diligent.

## The structure

| Role | Can write? | Decides whether the work is correct? |
|---|---|---|
| `ppropt-orchestrator` | no | **no** — relays the gate |
| `ppropt-researcher` | no | no — verifies claims only |
| `ppropt-worker` | yes | no |
| `ppropt-tester` | **no** | no — reports the gate |
| `validation/gate.py` | — | **yes** |

The only thing that decides is a script. That is the whole idea.

Each agent runs on `opus` at `effort: max`, in its own context, and is told what it may not
do as firmly as what it may. The tester cannot edit files, so it cannot fix a failing test to
make a run green. The worker cannot run the gate, so it cannot grade itself. The orchestrator
cannot write production code, so it has nothing to defend.

## The gate

`validation/gate.py` runs seven checks and exits non-zero if any fails:

| check | catches |
|---|---|
| suite | anything the 1231 tests cover, counted from JUnit XML |
| python parses | a source file corrupted by shell escaping |
| notebooks current | a generator edited without rebuilding its notebook |
| notebook line endings | CRLF, which `.gitattributes` turns into a CI failure |
| citations resolve | a `file:line` past the end of the file, or a missing file |
| banned phrases | wording `CLAUDE.md` forbids |
| working tree | reports what is uncommitted (advisory) |

```bash
cd C:/CodePractice/Drylab/clippr
PYTHONUTF8=1 .venv/Scripts/python.exe validation/gate.py
```

`--quick` skips the suite and reports `INCOMPLETE`, which **exits 1**. A gate that has not
run the suite has not gated anything, so it cannot report success.

It reads counts from JUnit XML rather than pytest's stdout because the summary line does not
reliably survive redirection in this environment — the specific quirk that once let an
unfinished run look like a passing one.

## What this does and does not guarantee

**It guarantees** that the verdict is a number, not a judgement. There is no wording to
soften and no reviewer to persuade. An agent cannot argue six failures down to zero, because
the count comes from a file pytest wrote.

**It does not guarantee** that `gate.py` is immutable. Any agent with write access to the
repository could edit it, and saying otherwise would be exactly the kind of overclaim the
gate exists to prevent. What exists instead is tamper-evidence:

- `gate.py --self` prints its SHA-256. The orchestrator records it at the start of a task and
  compares at the end; a change is reported **above** the result, not alongside it.
- `git diff validation/gate.py` shows any edit.
- The worker's definition forbids touching it, and the tester's definition requires flagging
  a changed hash as a finding in its own right.

That is honest tamper-evidence rather than tamper-proofing. The remaining protection is that
a human sees the hash.

**It does not check that a citation says what the prose claims.** The script can only confirm
the line exists. Whether it supports the sentence built on it is `ppropt-researcher`'s job,
and the reason that role quotes every line rather than asserting it checked.

## Using it

Ask for the orchestrator by name and give it the task:

> Use ppropt-orchestrator to add the intron advisory output to the designer notebook.

It will record the gate hash, delegate research if the task rests on a claim, brief the
worker, run the tester, loop on failure up to three times, and report the gate's verdict in
the gate's own words.

For a small, self-contained change the full chain is overkill — run the gate yourself. The
chain earns its cost when a change touches an invariant, rests on a published number, or
alters anything a reader will quote.

## The baseline

As of `a8f541b`: **1231 tests, 1197 passed, 34 skipped, 0 failures, 0 errors.**

The 34 skips are expected: thirty-three release-bundle controls declining to validate against
a bundle built from different source, and one target that passes the final-artifact contract.
If that number moves, a test has stopped protecting something, and that is worth a sentence
rather than a shrug.
