# What intron placement costs in assembly fidelity

*Measured 2026-10-16 on branch `feature/intron-advisory`. Reproduce with
`validation/experiments/intron_fidelity_tradeoff.py`; results land in `work/introns/`.*

---

## The question

Putting introns into a Golden Gate design adds junctions, and every junction competes in the
same reaction. The overhang at an intron junction is not free: it has to leave the intron's
own ordered part inside the GC band, and it has to sit between two guanines, which is Baier's
splicing rule. So the question is how much predicted assembly fidelity that costs.

PPROPT's baseline is **0.8282** for 9S and 19S and **0.8254** for 14S, and that number is
itself capped by the level 0 destination pair — `CTCA`/`CTCG` score 0.8282 on their own
(`docs/intron_feasibility.md` §5 has the mechanism).

## The first answer was wrong, and wrong in a dangerous direction

Taking the first usable site past each spacing target keeps selecting the **same four-base
overhang**, because every candidate sits between two Gs. Identical overhangs let fragments
swap places. Measured:

| | naive walk | baseline |
|---|---:|---:|
| 9S | 0.1795 | 0.8282 |
| 14S | 0.0307 | 0.8254 |
| 19S | 0.0026 | 0.8282 |

That placement was briefly what the advisory recommended. It would not have assembled. It
only surfaced because 14S and 19S were run — on 9S every ordered piece passed the synthesis
rules and the result looked fine.

## Making the overhang choosable

At a codon-aligned position the junction overhang can be chosen by synonymous substitution,
which is how PPROPT already picks its own cuts (`arelf.safe_overhangs`). Requiring the
position to be codon-aligned costs candidate sites but buys options at each:

| | sites with a fixed overhang | with a choosable one | total options |
|---|---:|---:|---:|
| 9S | 16 | 34 | 49 |
| 14S | 25 | 57 | 82 |
| 19S | 32 | 75 | — |

## The trade-off, and the lever that actually matters

Positions evenly spaced, overhangs solved with `overhangs.best_set`. Best found:

| architecture | introns | **level 0** | **level 1** |
|---|---:|---:|---:|
| 9S | 1 | 0.758 | **0.998** |
| 9S | 2 | — | 0.869 |
| 9S | 3 | 0.624 | 0.822 |
| 14S | 2 | 0.562 | 0.762 |
| 14S | 3 | 0.687 | **0.998** |
| 14S | 5 | — | 0.546 |
| 19S | 3 | 0.509 | 0.762 |

**The destination level dominates everything else.** At level 0 nothing exceeds 0.76, because
the backbone pair caps the whole reaction at 0.828 before any intron is added. At level 1 the
pair scores 1.000 alone, and 9S and 14S both reach **0.998** — better than the intron-free
design at level 0.

Three introns is Baier's architecture, and on 14S it is also the best-scoring count.

## What a dash means

"No valid overhang combination; all candidates conflict" — the solver's words, not a timeout.
On 9S at four introns the space is **eight** combinations and every one conflicts. So those
are real negatives.

**But they are negatives at those positions.** Positions were fixed by even spacing and only
the overhangs were optimised, because searching positions and overhangs together is
intractable: 19S at seven introns is choosing from ~75 sites, billions of combinations. A
smarter position search might beat these numbers. **Read the table as a floor, not a ceiling.**

## Where this leaves each architecture

- **9S — works.** One intron at level 1 reaches 0.998. Three reach 0.822.
- **14S — works.** Three introns at level 1 reach 0.998.
- **19S — not demonstrated.** Best found is 0.762, with three introns at level 1. No
  placement above that was found, and nothing at four or more introns produced a valid set at
  the positions tried.

## What this does not establish

- **Nothing here is measured.** Fidelity is a prediction from Pryor's mis-ligation data, not
  an assembly anyone ran. It says nothing about whether a spliced intron raises expression in
  a designed PPR, which no one has tested.
- **Level 1 is not free.** Moving the destination changes what the construct is and what it
  assembles into. This measures the overhang arithmetic, not whether level 1 suits the
  experiment.
- **Positions were not exhaustively searched**, as above.
- **One seed, one host, one intron.** Seed 42, *C. reinhardtii* nuclear, rbcS2i1 throughout.
  LHCBM1i2 is the stronger intron and remains untested against any of this.
