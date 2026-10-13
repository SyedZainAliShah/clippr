# Test sweep — both notebooks, all six hosts, 13 October 2026

Everything below was run on this checkout. Raw results are in the session scratchpad; the
numbers here are copied from them, not estimated.

## 1. Automated suite

**1229 tests · 1195 passed · 34 skipped · 0 failures · 0 errors · 499 s.**

33 of the 34 skips are one fixture: `tests/test_release_evidence_binding.py` refuses to
validate against `work/repair_release_20261012-225743`, because that bundle was built from
source `45303b4f…` and this checkout is `d358c4b3…`. That is the control working — it will
not check today's code against yesterday's record. **Those 33 controls are not running until
the release bundle is regenerated**, which is a hand-driven multi-phase process not scripted
in this repo. The 34th skip is one target that passes the final-artifact contract here.

## 2. Designer notebook — 6 hosts × 5 targets = 30 designs

**30 of 30 orderable. 0 errors. The ordered-fragment constraint applied in all 30.**

Targets: `AAAAUGUGG` (9S), `UUACACGUGCGUAC` (14S), `CUAUCACAUCACAUAAGCG` (19S),
`GGGCCGGCG` (9S, GC-rich), `AAUAAUAAU` (9S, AU-rich and repetitive).

| architecture | CDS | protein | fragments | predicted fidelity |
|---|---:|---:|---:|---:|
| 9S | 906 nt | 302 aa | 4 | 0.8282 |
| 14S | 1371 nt | 457 aa | 6 | 0.8254 |
| 19S | 1836 nt | 612 aa | 7 | 0.8282 |

Those figures are identical across all six hosts, which is the expected result and a useful
check: the host changes the codons, not the protein, the cut plan or the junctions. Fidelity
depends on the overhangs, so it does not move with the host either.

**Confirmed directly:** for the same target across three hosts the protein hash is identical
and the CDS hash is different in all three — `c_reinhardtii_nuclear` `8b195b94…`, `e_coli`
`eac8976a…`, `s_cerevisiae` `06f2a388…`. This is what makes `organism` the lever that produces
a genuinely different design.

**Runtime varies by host, a lot**, because the codon table decides how hard the GC band is to
satisfy: `e_coli` 0.8–3.6 s for every target; `c_reinhardtii_nuclear` 0.9 s for 9S but 286 s
for 19S; `c_reinhardtii_chloroplast` 68–165 s throughout. Budget for minutes, not seconds, on
19S targets in a live demo, and prefer `e_coli` if you need it to be instant.

## 3. Inventory notebook — all six hosts, run as a user would

Driven through the real notebook, not the library: **12 cells, 0 failures, for each of the six
hosts.** Each produced 42 assembly-ready substrates, 1 pool, every substrate inside the
synthesis contract, and a distinct inventory version.

| host | modules improved | mean CAI before → after |
|---|---:|---:|
| `c_reinhardtii_nuclear` | 39 of 42 | 0.2079 → 0.6183 |
| `c_reinhardtii_chloroplast` | 13 of 42 | 0.3478 → 0.4790 |
| `e_coli` | 42 of 42 | 0.5415 → 0.9416 |
| `s_cerevisiae` | 10 of 42 | 0.7007 → 0.7436 |
| `a_thaliana_nuclear` | 24 of 42 | 0.8073 → 0.8790 |
| `n_tabacum_chloroplast` | 10 of 42 | 0.6232 → 0.6870 |

CAI predicts how well codons match the host; it is not a measurement of expression, and a
large gain is not a proportional expression gain. A module reported *unchanged* is a delivered
module, not a gap.

**`candidate_seeds` matters.** At `seeds=2` three hosts reported one substrate outside the
synthesis contract; at the notebook default of `4`, all six reported none. Leave it at 4 or
higher.

## 4. The bug this sweep found, and fixed

Step 6 called `select_interface` without the host codon table. The front's binding guard
correctly refused — *the objective values were measured under one codon table and would have
been published against another* — and the cell then raised `KeyError: 'inventory'` reaching
into a result it never checked.

It reproduced on **five of six hosts**. `c_reinhardtii_nuclear` was the exception only because
it is the cached default, so the hashes matched by accident. Verified before and after on
every host: `select_interface` without the table refuses on 5 of 6; with the table it succeeds
on 6 of 6.

The library tests pass `codon_table=` explicitly in every call, so they could never have
caught the notebook omitting it. Two static gates now cover it, plus
`validation/run_notebook.py --notebook … --set organism=…`, which runs either notebook with
form values, and release gate 7, which now executes **both** notebooks rather than only the
designer.

A second defect surfaced while fixing the first: nulling any refused result would have hidden
a usable 42-module order behind one flagged substrate on 3 of 6 hosts. Tasks that refuse
*without writing* (`explore_interfaces`, `select_interface`) now stop the chain; tasks that
write a *partial* result (`recode_for_host`, `load_and_compile`, `order_items_for`) hand you
the rest and name what fell short.

## 5. Release gate

Gate 7 passes: `PPROPT_designer.ipynb` 15 cells 0 failed, `PPROPT_inventory.ipynb` 12 cells
0 failed. The full `release_check.py` was **not** re-run end to end; gate 7 was run directly.
