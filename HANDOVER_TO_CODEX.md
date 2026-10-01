# Handover to Codex — 2026-09-17, revised 2026-09-28 and 2026-10-01

Everything done since the GC-policy review, what it settled, and what is genuinely still open.
Written to be read by someone who has the repository and none of the conversation.

## What changed on 2026-10-01

Two commits, both driven by the wet lab reviewer's feedback and by a decision that costing
should not exist at all. Neither changes a design this tool produces.

**`44446fe` — cost estimation removed.** We cannot obtain quotes, so every number the pricing
code produced was a list price of unrecorded age, and most of the surrounding code existed to
say so. Checked before cutting: cost fed **no objective** — `objectives`, `joint_search`,
`library_search`, `synthesis_fitness` and `qc` carry no reference to it. Gone: `PriceTier`,
`estimate_cost`, `opool_quote`, `_total_cost`, `LibraryResult.cost`, and the 109.00 EUR /
1.63 EUR constants. Kept: `check_eligibility`, `plan_pools`, `write_order_files`, none of
which ever needed a price.

Two consequences worth a reviewer's attention:

  - `LibraryResult.cost` is now `LibraryResult.totals` (designs, oligos, bases). The summary
    line loses the pooled-versus-separate comparison, which was a property of the price model
    rather than of any design.
  - `phasee_gate` measured the pooling heuristic's distance from optimal **in euros**. Under a
    flat per-pool tier that was a restatement of pool count, so it now measures pool count
    directly — the quantity the heuristic actually controls. Gap is zero at every n from 2 to
    10, as before.

**`b33a383` — the notebook split in two and cut down.** The reviewer could not make sense of
the inventory results, and the cause was structural: one document carried two complete
workflows sharing nothing but a Setup cell, so a reader reached the inventory route holding
the wrong mental model. It is now `notebooks/CLIPPR_inventory.ipynb`, with its own header,
host form and setup; both builders share `tools/notebook_kit.py` so a fix to the install cell
cannot land in one and not the other.

The designer notebook went from 40 cells to 22, and 3,256 words of prose to 1,573. "Explore
the alternatives" is **cut**: it carried a slider the reviewer did not understand, and its own
markdown recorded that predicted fidelity came out identical for every candidate in this
configuration. The design audit moved to the end, as machinery rather than an answer. The
off-target cell prints the verdict first instead of a background-rate table the reviewer asked
the purpose of.

**The tests badge read 283; the measured figure is 927.** It is now a named constant beside
the command that produced it. That badge is the same failure the thread below names: a number
in a document that nothing recomputed, wrong by 644 and green the whole time.

**Later the same day: the merge, and a proposal that failed review.** `synthesis-policy-profile`
was fast-forwarded onto `main` (`3100ff5..b6d0f68`) and the clean-install release check re-run,
taking it to 10/10 with nothing needing a human. `PROPOSAL_TO_CODEX.md` then scoped open issues
14 and 15 — and its first revision carried six defects, every one of them a number or a
categorical claim that had never been checked against the thing it described. Codex found all
six. Revision 2 corrects them in place with the retractions listed at the top, and the two that
are facts about this repository rather than about the proposal are recorded in item 14 below. The
headline figure moved from 67% to **51%**. This document already records several numbers that
shipped unchecked under a green suite; the only new part here is that a proposal is covered by no
suite at all, so nothing but a reviewer could have caught it.

---

## What changed in this pass, in one place

Your second review returned six findings and the verdict *do not merge*. All six reproduced,
all six are repaired, and each is now pinned by a regression in `tests/test_policy_boundaries.py`
that fails when the defect is reintroduced — verified by `validation/falsify_boundaries.py`,
which reverts each repair in a scratch worktree and checks the suite goes red. **10 of 10
reverted defects are caught.**

Fixing them turned up five more defects that no review had found, all of the same kind — a
number or a claim that looked plausible and had never been checked against the thing it
described:

| found while | defect |
|---|---|
| writing the regressions | a profile could **advertise** enzyme sites as soft while every consumer treated them as absolute; the recorded binding carried that claim |
| writing the regressions | `select_interface` wrote the selected inventory to disk **before** verifying it, so a rejected candidate left a file for the next step to find |
| falsifying the regressions | three of my own regressions did not discriminate what they claimed to — the falsifier caught them, not me |
| wiring `level1_readiness` | fidelity scored over **every contributed end** rather than the distinct junctions, returning 0.0039 where the correct figure is 0.9940 |
| adding a notebook policy selector | the notebook imported a module that exists only on this branch while its Setup cell installs from `main` — **every Colab session broke**, and the in-tree notebook test could not see it. Reverted; see §7 item 9 |

**Two published conclusions are withdrawn and replaced by measurements**, both because the
experiment behind them did not test what it said:

1. *"Soft enforcement has no independent value."* The soft arms were never soft. Rerun against
   the repaired solver, soft enforcement at the **strict** thresholds reaches mean CAI
   **0.819921** — identical module for module to widening the band — against 0.618830 hard.
   The retracted claim was the opposite of the truth. §3.
2. *The runtime table.* Its script measured no chain despite saying it did, and reported a sum
   of medians no run ever took. Rewritten; end to end is **4.70 s** default, **6.26 s** broad. §5.

And one comparison ran for the first time and **we lose it**: cross-module sharing, 0.964147
ours against 0.911273 theirs. §7, item 7.

**The thread running through all of it**, and the reason the regressions are written at the
workflow boundaries rather than as unit tests: *a parameter in a signature is not evidence that
the consumer uses it; a number in a document is not evidence that anything computed it; and a
passing test is not evidence that it tested what users get.* Every defect above sat under a
green suite, and the last one sat under a green suite that was structurally incapable of
seeing it.

**Revision note, 2026-09-28.** Two claims in the first version of this document were wrong and
are corrected here rather than silently edited: the State table said the notebook had a policy
selector and 22/22 cells executing, and open-issue 9 was marked Closed. Both described work
reverted the following day. The notebook figures below are re-measured, not recalled.

**State**

| | |
|---|---|
| branch | `synthesis-policy-profile` at `b6d0f68`, **level with `origin/main`**. The merge landed as a server-side fast-forward, `3100ff5..b6d0f68`; open issue 1 is closed |
| commits | `74e4242` · `a7181fe` · `2f936bc` · `35898a1` · `fa98fc6` · `3b737a9` · `44446fe` cost removed · `b33a383` notebook split · `b6d0f68` handover and notebook step labels — all pushed, and all now on `main` |
| tests | **927 passed of 927**, 0 failed. Was 949; the 22 that went were the ones asserting the costing feature that was removed |
| release check | **10 / 10**, 0 failed, **0 needing a human.** Re-run with `--clean-install` after the merge. Check 3 reproduced the README's worked example verbatim from a pip-installed build — `AAAAUGUGG (9S) -> 302 aa, 906 nt, 4 fragments. Fidelity 0.828. QC PASS.` — which is what re-verified the costing removal on both sides at once, since it compares strings |
| falsification | 10 / 10 reverted defects caught. **Not re-run this pass** — no boundary test was added or changed, and nothing in `tests/test_policy_boundaries.py` was touched |
| notebooks | **Two.** Designer: 22 cells, 12 code cells run, 0 failed. Inventory: 20 cells, 10 code cells run, 0 failed. Both against a *pip-installed* build from `origin/main`, not the working tree |
| default behaviour | **unchanged** — every published number reproduces |

---

## 1. Your four corrections — all accepted, two independently verified

| your finding | status |
|---|---|
| "tighter than any published vendor guideline" is false | **confirmed by me**. I fetched Twist's page: *"global GC% of less than 25% or more than 65% and local GC windows (50 bp) of less than 35% or more than 65%."* Our band and window exactly. |
| one constant is insufficient | **confirmed by reading our own code**. `recoding.GC_BAND` was the validator's; `optimize_cds` kept its own `gc_bounds` default that `recode_inventory` never passed. |
| the experiment changed two parameters | accepted; your 2×2 reproduced and adopted |
| the reference has reusable exports | **confirmed by reading theirs**. `workflows.export_optimized_library` writes CSV/FASTA/Excel/GenBank; `compile_and_assemble_target` takes an already-optimised library. |

On the first: I had compared against the reference's *transcription* of Twist's rules — in a
document where I had myself written that the transcription should not be trusted for that
purpose. I wrote the caveat and did not follow it.

What survives: the band has **no recorded source in this repository**, Twist must not be
back-dated as its origin, it is Twist's *codon-optimisation* guidance while we order from IDT,
and we enforce the local half without the global 25–65% it is paired with.

---

## 2. The constraint profile — implemented, then twice found not to be

`src/clippr/synthesis_profile.py`. One versioned object carrying thresholds **and** their
provenance, consumed by every stage that judges a sequence.

**The three-way enforcement split you asked for:**

    hard              a breach is infeasible; the candidate is refused
    target            a breach is reported and delivered; an optimiser steers away
    vendor_judgement  not machine-checkable here; recorded so it is not mistaken for checked

An unlisted rule is **hard** — silence must not relax a constraint.

No profile can soften a forbidden site, and the way that is enforced changed under review. The
solver had always declined to move an enzyme site to its objective list, so a lenient profile
was simply ignored on that rule — but `solver_bounds()` still *advertised* `forbidden_sites`
among the soft rules, and the recorded binding carried it. An artefact could therefore state
that sites were advisory while every consumer had treated them as absolute. Declaring one soft
now raises at construction: refusing to build the policy is a stronger guarantee than every
consumer independently declining to honour it.

`substrate_problems` returns only hard breaches; `substrate_warnings` returns targets. They
never share a return value, so a soft warning cannot silently become approval and a target
breach cannot silently block a strict export.

**Three shipped profiles.** `clippr-strict-legacy` is the default with `source` recorded as
*unrecorded*; `twist-codon-optimisation-guidance` adds the global band and makes both advisory
as that page frames them; `broad-experimental` is measurement-only.

**The proof it reaches the solver, not just the validator:**

| profile | mean CAI | improved / unchanged |
|---|---:|---|
| `clippr-strict-legacy` | **0.618830** | 39 / 3 |
| `broad-experimental` | **0.819921** | 42 / 0 |

0.819921 is unreachable by patching the validator alone. That is why the regression asserts it
rather than asserting the profile object exists.

**Threading — and what "completed" turned out to mean.** My first pass wired the engine but
not the workflow layer, and I reported it as done. I then wired the workflow signatures and
reported *that* as done too. A second review found the second claim was the same error as the
first: `order_items_for` took a `profile=` and exported under the default, `library_search`
resolved a profile and passed none of its bounds, `joint_search` recorded no policy on the
front, and `select_interface` rechecked under whatever the default happened to be.

Four of that review's six findings were a consumer ignoring an argument it accepted. So the
useful statement is not that the entry points take a profile — they took one before — but that
each consumer is now pinned by a regression that fails when it stops reading it. Those live in
`tests/test_policy_boundaries.py`, and `validation/falsify_boundaries.py` reverts each repair
in a scratch worktree to prove the regressions are not vacuous.

The lesson is worth carrying: **a parameter in a signature is not evidence that the consumer
uses it**, and a grep for `profile=` cannot tell the two apart.

---

## 3. The four-arm policy benchmark — item 5 of your work order

`validation/experiments/policy_benchmark.py`, all 42 modules, matched seeds and budgets, each
arm judged on its own policy. Artefact `work/policy/policy_benchmark.json`.

| arm | mean CAI | improved | worse | hard breaches | advisory |
|---|---:|---:|---:|---:|---:|
| `strict_hard` (default) | 0.618830 | 39 | 0 | 0 | 0 |
| `broad_hard` | **0.819921** | 42 | 0 | 0 | 0 |
| `strict_soft` | **0.819921** | 42 | 0 | 0 | **40** |
| `broad_soft` | 0.819921 | 42 | 0 | 0 | 0 |

Rerun 2026-09-17. The first run of the two soft arms did not implement soft enforcement and
reported `strict_soft` at 0.657195 with 3 advisories; see section B.

**A — hard vs hard.** Broad ahead on 40 of 42, strict on 0, two tied. Mean +0.201091, largest
single gain +0.670973, **largest loss +0.000000**.

**B — hard vs soft at the same thresholds. Retracted once, now measured, and it reverses.**

The first run of this arm did not implement soft enforcement: `solver_bounds()` returned
identical arguments for `hard` and `target`, so both "soft" arms were a hard solver read by a
permissive validator. I reported its conclusion — *"soft enforcement has no independent value"* —
and that was wrong. Rerun 2026-09-17 against the repaired solver:

| arm | mean CAI | improved | advisory warnings | seconds |
|---|---:|---:|---:|---:|
| `strict_hard` (today's default) | 0.618830 | 39 / 42 | 0 | 2.4 |
| `broad_hard` | 0.819921 | 42 / 42 | 0 | 2.7 |
| **`strict_soft`** | **0.819921** | **42 / 42** | **40** | 6.2 |
| `broad_soft` | 0.819921 | 42 / 42 | 0 | 2.2 |

**Soft enforcement at the strict thresholds recovers the whole gain of widening the band** —
`strict_soft` and `broad_hard` agree on all 42 modules to six decimal places. The strict band
costs 0.201091 mean CAI as a wall and nothing as a gradient. The retracted version of this arm
scored 0.657195 and recovered three modules; that figure was the permissive validator, which is
why it coincided exactly with the one-line `GC_BAND` patch.

`broad_soft` ties `broad_hard` with zero warnings — the control, confirming the mechanism costs
nothing at thresholds that do not bind.

**This is a decision you should make, not one I have made.** `strict_soft` dominates on CAI and
on information retained, but it *delivers* sequences outside Twist's published local-GC guidance
where the hard policy refuses to. That is a guarantee traded for a documented warning, and it is
worth having only if someone reads the warnings. The default is unchanged pending that call.

Arm A (hard vs hard) is unaffected and stands: broad ahead on 40 of 42, mean +0.201091, largest
single gain +0.670973, largest loss +0.000000.
---

## 4. The degenerate third axis — closed

`repeat_burden` was identically zero, so a search declaring three objectives ranked by two.

**The obvious fix does not work.** The reference uses `repeat_k = 16` against our 20, which
invites the conclusion that the parameter was wrong. Measured:

| k | 20 | 16 | 12 | 10 | 8 |
|---|---:|---:|---:|---:|---:|
| modules with any duplicate | 0/42 | 0/42 | 0/42 | 0/42 | 5/42 |

Zero to k = 10. A repeated 10-mer inside a ~100 nt module is genuinely rare, so the **measure**
was wrong. A test pins this so nobody ships the k change as a fix.

`synthesis_fitness` replaces it — GC centrality, homopolymer headroom, internal repetition and
collection sharing, on the ordered substrate, each a *distance* so ranking moves when
everything comfortably passes.

| | old axis | new |
|---|---:|---:|
| deposited | 0 (degenerate) | mean 0.330694, spread 0.082140 |
| recoded | 0 (degenerate) | mean 0.262544, spread 0.085402 |

**Gate D's front grew 5 → 9 with 6 distinct fitness values. The recommendation is unchanged.**

**The weights were tested for load-bearing, not assumed harmless.**
`validation/experiments/fitness_weight_sensitivity.py`, six weightings including each term
taken to dominance:

| weighting | spread | front | recommendation |
|---|---:|---:|---|
| shipped | 0.085402 | 9 | unchanged |
| equal | 0.073520 | 9 | unchanged |
| gc only | 0.133522 | 5 | unchanged |
| homopolymer only | 0.154099 | 8 | unchanged |
| sharing only | 0.166598 | 6 | unchanged |
| **repetition only** | **0.000000** | 5 | unchanged |

The recommendation is identical under all six, so the weights are not load-bearing for the
choice. Front *membership* moves (5 to 9), so they do change which trade-offs are visible. The
last row is the old axis in isolation and is the only one that degenerates — the sweep
reproduces the original defect independently.

Two costs, reported rather than absorbed: the recoded inventory scores **lower** than the
deposited one, because recoding buys CAI by pushing GC to the band edge; and the axis initially
made the interface search 2× slower. That slowdown is fixed — `inventory_fitness` now computes
sharing from k-mer document frequency in one pass instead of rebuilding every neighbour's set
per module. Identical results (verified per-module to 1e-12), 133 ms → 82 ms, and the search is
back to 0.917 s.

---

## 5. Runtime — §12's dimension, measured twice

Neither system had instrumentation. The reference's 21 modules contain no timing at all;
CLIPPR had one figure for the 200-design corpus.

**The first version of this section was wrong about its own method.** Its script claimed in its
docstring to time stages "in isolation and again as a chain", and no chain existed: it called
the collection optimiser and discarded the inventory it returned, ran a search whose result
never reached an export, and reported a **sum of medians** — 3.270 s, a number no run ever took.
A sum of isolated stages also understates the workflow systematically, because the later stages
are cheap on the deposited kit and expensive on the recoded, optimised inventory they receive.

Rewritten and rerun, 2 repeats, `validation/experiments/runtime_profile.py`:

**Isolated** — every stage from the same fixed input, for comparing stages with each other:

| stage | median |
|---|---:|
| load deposited kit | 0.276 s |
| recode for host (42 modules, 4 seeds) | 2.173 s |
| compile 3 targets | 0.000 s |
| optimise collection (greedy, 42) | 0.788 s |
| build junction classes | 0.001 s |
| interface search (100 evaluations) | 1.149 s |
| build 42 order substrates | 0.005 s |
| check order eligibility | 0.005 s |
| *sum of the above* | *4.397 s — not a workflow time* |

**Chained** — each stage consuming the previous stage's actual output, timed end to end:

| route | end to end (median) | front | export |
|---|---:|---|---|
| default (`clippr-strict-legacy`) | **4.696 s** | 48 / 100 feasible, 9 on the front | 42 items, 0 failures |
| explicit `broad-experimental` | **6.260 s** | 85 / 100 feasible, 15 on the front | 42 items, 0 failures |

Two things worth noting. The broad route costs **33% more wall clock** and returns nearly twice
the feasible candidates and a front two-thirds larger — the extra time is spent on work the
strict route never had to do because the candidates were vetoed. And under the broad policy the
collection optimiser reports **no improvement**: the recoding already reached what it would have
found, so a stage that looks productive under one policy is redundant under another.

Every repeat is recorded, not a best-of. It deliberately does **not** time the reference:
running their optimiser under our budgets would compare different amounts of work and call it
speed.

---

## 6. Level-1 — narrowed from blocked to fillable

**The reference does not solve this.** Its `grasp_level1_reaction_overhangs` composes
`final_cassette (GGAG/CGCT) + ppr_outer (AGGT/TTCG) + block_joins`. Searching its whole package
shows the P2L2S2 linker and DYW domain appear **only in bundled GenBank data, never in reaction
composition** — and it emits `level1_fidelity` from that incomplete set regardless. On this
point CLIPPR is now stricter than the reference.

What their code *did* give us is the composition rule: four flanking overhangs, not two. Our
`DESTINATION_OVERHANGS["level1"] = ("GGAG","CGCT")` was never a rival to the derived
`AATG…TTCG` — it is the final cassette, and the real reaction contains both pairs.

`assembly_spec.level1_participants()` now makes the gap fillable: supply co-assembled ends
**with evidence naming their source** and the reaction scores, carrying its basis. Without
evidence it stays unscored and says why.

Demonstrated as a hypothesis, not a result: block subset + cassette ends + the three known
co-assembled ends gives a ten-overhang 19S set scoring **0.623247**, *below* level 0's
0.751683. If that participant list is right, level 1 is the bottleneck.

---

## 7. What is genuinely still open

**Decisions, not work:**

1. ~~**Merge `synthesis-policy-profile`.**~~ **Closed 2026-10-01.** Merged as a server-side
   fast-forward, `3100ff5..b6d0f68`, 9 commits. Verified in this order: full suite green
   *before* main was touched; fast-forward confirmed with
   `git merge-base --is-ancestor origin/main HEAD`; then the throwaway venv reinstalled from
   the **new** `main` and both notebooks re-run against it — designer 12 code cells, inventory
   10, 0 failed. The reviewer now opens a Colab notebook containing the six fixes from
   `3b737a9` and everything from this pass.
2. **Whether the default profile changes.** The benchmark is favourable and is **not** evidence
   of vendor acceptance, which is the only thing that matters for an order. Someone must own
   that risk.

**Blocked on the outside world:**

3. **Level-1 participants.** The mechanism exists and is now reachable — `level1_readiness`
   takes the list and scores it. What does not exist, on either side, is the linker/DYW
   context itself. This is blocked on a document, not on code.
4. **Vendor acceptance.** Needs a real submission to IDT's checker.
5. **External users.** Being arranged.

**Work, in the order I would take it:**

6. ~~`level1_participants()` is called by nothing.~~ **Closed.** `workflow.level1_readiness`
   exposes it: supply the linker, editing-domain, bridging and final-cassette ends plus a
   citation, and it either refuses with the roles it lacks or scores the reaction. Scoring it
   found a second defect on the way in — the first version passed *every contributed end*, so
   each internal junction competed with a duplicate of itself and returned 0.0039 where the
   distinct junctions give 0.9940. A quiet failure, because 0.0039 reads as a badly chosen
   overhang set rather than a counting error. Pinned by a regression.
7. **§12's matched comparisons — repeat/sharing now closed, and we lose it.** Measured with
   one implementation over both systems' coding spans (k = 12, reverse-complement folded):
   internal repetition is **0.000000 on both sides**, but cross-module sharing is **0.964147
   ours against 0.911273 theirs**, and their 42 designs draw on **1005 distinct 12-mers to our
   552**. Per module we share more on 33 of 42.

   **I got the mechanism wrong on the first write-up and the correction sharpens it.** I said
   nothing on our side sees collection sharing. In fact `library_search.optimise_library`
   weights a `sharing` term at **1.0, equal to adaptation**. What has none is
   `recoding.recode_inventory` — which produced the very inventory the benchmark compares. So
   the table above pitted our no-sharing artefact against their with-similarity one.

   The fair version: our collection optimiser's output (`inventory_greedy.json`) scores
   **0.953602 sharing over 621 distinct 12-mers** — real movement, and about a fifth of the
   distance to their 0.911273 / 1005. The gap survives the correction.

   Two plausible reasons, neither measured: their similarity penalty sits **inside each
   module's solve**, so every codon choice is scored against the collection as it is made,
   while ours ranks finished proposals over a bounded budget; and we recode every module toward
   the same optimum, so identical stretches get identical codons unless something pushes them
   apart. Moving the term into the solve, and giving `recode_inventory` one at all, are the two
   changes this points at — each trading adaptation for diversity at an unmeasured rate, so
   each needs its own experiment rather than being assumed an improvement.

   **Stage-specific fidelity remains unmatched**, and runtime is measured for us only — timing
   their optimiser under our budgets would compare different amounts of work.
8. **Four reference-derived ideas unadopted**, and one of them now has a number behind it.

   **Library similarity inside the per-module solve** moves to the top: the matched comparison
   in item 7 measures what it is worth on this kit. Their similarity penalty is applied as each
   codon is chosen; ours ranks finished proposals. That is the difference between 1005 and 621
   distinct 12-mers over the same 42 proteins.

   The other three are unchanged and unmeasured: a minimum codon-adaptiveness floor; scoped
   rescue codons; protocol provenance on the ligation matrices.

   None of the four should be adopted on the strength of the reference doing it. The similarity
   change in particular trades codon adaptation for diversity at a rate this project has not
   measured, and the four-arm benchmark is the template for how to settle that before shipping
   it — matched seeds, matched budgets, each arm judged on its own terms.
9. **The notebook cannot select a profile. Still open — I closed it, broke Colab, and
   reverted it.** This is the most useful thing in this document, so it is written out in full.

   I added an *Inventory 1b* form to `tools/notebook_inventory_cells.py` offering the three
   policies, threaded into the recode, the collection optimiser, the interface search and the
   export. It worked. The full suite passed, and `validation/run_notebook.py` reported every
   cell executing.

   It broke every Colab session. The cell imports `clippr.synthesis_profile`, which exists on
   this branch and **not on `main`** — and the notebook's own Setup cell installs with
   `pip install git+https://github.com/SyedZainAliShah/clippr.git`, which resolves to the
   default branch. So cell 1b raised `ModuleNotFoundError`, `policy` was never bound, and every
   downstream cell failed on `NameError`. `order_items_for` on `main` takes no `**kwargs`, so
   `profile=` was an immediate `TypeError` rather than a harmless extra.

   **Why the green suite did not catch it.** `validation/run_notebook.py` inserts the working
   tree's `src/` at the front of `sys.path`. It therefore tests whether the notebook agrees
   with *the code you are editing* — which is useful, and is structurally incapable of catching
   a notebook that depends on something the published package lacks. The test and the defect
   were in disjoint spaces.

   **The repair, and what is now in the repository.**
   `validation/notebook_on_published_build.py` runs the notebook's code cells against a
   virtualenv with CLIPPR pip-installed from GitHub, with the working tree kept off the path.
   Current result, measured rather than recalled: **21 cells ran, 0 failed** against the `main`
   build, and the same against this branch's build — so merging cannot break the notebook.

   The selector itself is reverted and unshipped. The shape it needs is a cell that detects
   whether the installed build has `synthesis_profile` and degrades to the defaults when it
   does not, rather than assuming it. A draft exists and is not in the repository, because it
   still has an escaping bug and has never been run against a `main` install.

   **The general rule this earns:** anything that ships to users who did not choose their
   package version has to degrade, not break — and it must be tested against the version they
   will actually get.
10. ~~Naming collision.~~ **Closed.** `plan_order`'s parameter is now `vendor_profile`. Every
    call site passed it positionally, so nothing broke. Note the notebook has a third meaning
    again in `enzyme_profile` — which restriction sites to avoid — so whoever finishes the
    selector in item 9 should name its control `synthesis_policy`, not a fourth `profile`.
11. ~~`apply_assignment` raises a bare `KeyError`.~~ **Closed.** It now refuses with a message
    naming the unknown class and listing the available ones, and separately refuses an overhang
    absent from a class's options — applying one would publish a design nobody scored.
12. ~~`runtime_profile.py` measures no chain.~~ **Closed and rerun.** Its docstring claimed
    stages were timed "again as a chain" and none was: it discarded the collection optimiser's
    inventory and fed no search result to an export, then reported a **sum of medians**. It now
    runs a real route where each stage consumes the previous one's output, under both a default
    and an explicit non-default policy. Measured end to end: **4.70 s** default, **6.26 s**
    broad-experimental, against a sum-of-isolated-medians of 4.40 s that no run ever took.

**Still open:**

13. **A default decision.** `strict-as-target` now matches the broad band's CAI exactly while
    keeping all 40 advisories, so the choice is no longer "conservative or optimised" but
    "refuse or record". It ships out-of-band sequences the strict default refuses to produce.
    That is a judgement about risk, and it is not mine to make.
14. **Region-based targeting. Scoped and reviewed, not built — `PROPOSAL_TO_CODEX.md`.** The
    reviewer recommended specifying a region of the target gene rather than a fixed
    9/14/19-mer. Still true: given a bare k-mer the tool cannot separate the intended binding
    site from a genuine off-target, which is what the reviewer asked when a hit landed in
    ORF1995.

    **The first version of that proposal was wrong in six ways and all six were mine.** Codex
    reviewed it and every finding reproduced. Two are worth carrying here because they are
    facts about this repository, not about the proposal:

      - **`load_transcripts()` is not a gene dictionary.** 109 annotated features, **100
        distinct names**; nine names occur twice. `psbA` sits at 48774 (minus strand) and
        138789 (plus strand) with **byte-identical** extracted sequence — the chloroplast
        inverted repeat. `Hit` carries neither feature identity nor genomic coordinate, so name
        plus local offset cannot tell two copies apart. Any occurrence accounting must key on
        feature identity. My proposal keyed on name, which silently forgave both copies and
        turned a 51% result into a 67% one.
      - **The CDS features are coding spans, not UTR-inclusive transcripts**, and `psbA`'s
        location is five joined segments. So `start + offset` is not a genomic coordinate, and
        a native-UTR workflow cannot be served from this cache at all — it needs a data source
        that does not yet exist here.

    What the review settled, against my proposal: the feature does **not** change the input
    model. `design_oneshot`, `design_library` and `parts.select` already accept a chosen
    sequence, so a selector sits upstream and hands one over — purely additive. And window
    selection applies to the **inventory** route too, which the proposal wrongly excluded:
    every window in a 200 nt sample was kit-available and every window produced a **different
    pick list**.
15. **Notebook UX beyond the cut.** The two routes are separated and the prose is down to 1,573
    words, but the designer notebook is still a linear scroll of Colab forms. Grouping inputs
    and outputs, and making detail collapsible, is unstarted — and depends on whether Step 1
    stays a sequence box, so it is held behind item 14.

    ~~**Two corrections are not held, because they are wrong in the repository today.**~~
    **Both shipped in `e1a09a7`.** The off-target cell's markdown had said an occurrence means
    "the protein binds there too", contradicting `offtarget.py`'s own docstring three paragraphs
    below it; and it printed `architecture_advice` beside a **transcript** verdict to argue the
    hit was "arithmetic", when that figure is a **genomic, both-strand** composition null —
    the pairing `report()` already warns about in prose, and which an earlier reviewer did read
    as a contradiction.

    `report()` changed with the markdown, which is why that commit touches `src/`: it printed
    "These are off-targets", so correcting only the cell would have left it contradicting itself
    one screen lower, in text the notebook prints rather than renders. Each occurrence is now
    "a candidate to review, NOT a demonstrated off-target". The cell also now states what was
    searched — 109 features, nine duplicated names, coding spans rather than UTR-inclusive
    transcripts — and what was not: the user's own construct and the rest of their library.

    Verified the way this repository requires, because `offtarget.py` is in the package Colab
    installs: 927 of 927 pass, and the notebook was run against a pip-installed build **three**
    times — before the edit, after the edit against the old build, and again after the merge
    against a build reinstalled from the new `main`, which is the only run that exercises the
    new `report()`. Designer 12 code cells, inventory 10, 0 failed each time.

    On the web app: the ordering argument in the proposal cited **4.70 s / 6.26 s** as the
    designer route. Those are `chained_routes` in `work/runtime/runtime_profile.json` — the
    **inventory** chain. The one-shot designer figures are a separate section,
    `oneshot_synthesis_route`: **0.5221 s** 9S, **14.1180 s** 14S, **204.9769 s** 19S, each
    n = 1 with screening disabled. The conclusion survives and hardens — 205 s cannot be a
    synchronous request — but anyone repeating the argument should quote the right table.

---

## 8. What none of this establishes

- **Vendor acceptance of anything.** No sequence from either system has been synthesised.
  Twist and IDT both score submissions with proprietary models.
- **Optimiser superiority in either direction.** Both share a mean-log-adaptiveness component;
  theirs also carries weighted GC, homopolymer, repeat and similarity penalties. Not the same
  objective, so no overall ranking follows. Two axes now *have* matched numbers and both are
  reported without spin: codon adaptation over the aligned span is 0.621004 ours against
  0.658672 theirs, and cross-module sharing is 0.964147 ours against 0.911273 theirs. Two axes
  are not a verdict on an optimiser, and on neither has any synthesis outcome been observed.
- **Whole-workflow parity.** A 10/10 release check is a statement about that checklist's scope.
- **That a higher CAI is a better protein.** `docs/objectives.md` §8 holds.
- **That any weighting of `synthesis_fitness` is correct.** The sweep measures sensitivity, not
  truth.

---

## 9. Where to look

| | |
|---|---|
| `PROPOSAL_TO_CODEX.md` | open issues 14 and 15 scoped, reviewed and corrected — read its revision note before quoting any figure from it |
| `docs/synthesis_policy.md` | the profile mechanism and the four-arm benchmark |
| `docs/gc_band_and_constraint_model.md` | the GC measurement, vendor sources and their limits |
| `docs/reference_implementation_findings.md` | everything established about GRASP Designer |
| `docs/objectives.md` §0, §4a | the scope rule; the composite axis and its sensitivity sweep |
| `FULL_WORKFLOW_HANDOFF.md` §14b | the three gaps closed overnight |
| `RESPONSE_TO_AUDIT.md` §12 | point-by-point reply to your review |
| `docs/verification_commands.md` | one invocation per operation |
| `tests/test_policy_boundaries.py` | the regressions pinning each of your six findings |
| `validation/falsify_boundaries.py` | proof those regressions are not vacuous — run it |
| `validation/notebook_on_published_build.py` | runs the notebook against a pip-installed build, which is the only test that can catch a notebook/package version split |
| `CLIPPR_FOR_THE_WET_LAB.md` | the briefing written for biologists — the metric glossary in Part 3 is the plainest statement of what every number means |

**Reproduction:** every experiment above writes a JSON artefact under `work/`, and every
command is in `docs/verification_commands.md`. The primary GenBank records resolve relative to
the repository or via `CLIPPR_GRASP_PRIMARY`; there are no absolute developer paths left.
