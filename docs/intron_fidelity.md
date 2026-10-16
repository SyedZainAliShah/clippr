# What intron placement costs in assembly fidelity

*Measured 2026-10-16 on branch `feature/intron-advisory`. Reproduce with
`validation/experiments/intron_fidelity_tradeoff.py`; results land in `work/introns/`.*

---

## Where it ends up

**Level 0, which is the case that matters**, because a PPR coding sequence is a reusable part
and gets assembled into a level 0 acceptor. MoInClo does the same: exons and level −1 intron
parts go into a level 0 backbone, which then feeds level 1 and 2.

| architecture | introns | fidelity | baseline | cost | |
|---|---:|---:|---:|---:|---|
| 9S | 3 | **0.803** | 0.828 | 0.026 | usable |
| 14S | 3 | **0.751** | 0.825 | 0.075 | usable |
| 19S | 4 | **0.686** | 0.828 | 0.143 | marginal, warned |

Three introns is Baier's own architecture. Introns always cost something here, and the cost
grows with the protein, but none of the three is closed.

**That 9S row is seed 42, and seed 42 is an outlier.** Across seeds 1, 7, 42, 99 and 2024 on
the same target, only seed 42 reaches 0.803; the other four land at exactly **0.7507**. Two of
the five take **two** introns rather than three, because a third would fall to 0.572 or 0.680.
Seed 42 also has the most candidate positions of the five, 65 against 54–57, and is the only
one with any usable route-A position at all.

So the measured claim is **two to three introns at 0.75–0.80**, and the 3-at-0.80 result holds
for one seed in five. The earlier framing overstated a lucky draw.

## The ceiling is the acceptor, not the introns

PPROPT's baseline is 0.8282, and that is already a cap. From `assembly_spec.py:15`:

> Level 0 is stored `("CTCA", "CGAG")`, but the strands that actually meet in the tube are
> **CTCA and CTCG** — which differ by a single nucleotide and score 0.794.

So the reaction is capped at about 0.83 before an intron exists. Introns then eat into what
is left.

Checked against a destination whose two ends are *not* near-identical — level 1's
`GGAG`/`AGCG` — the same placements reach **0.998** on 9S and 14S. That number is not a
recommendation to clone into level 1, which would make a transcription unit rather than a
reusable part. It is a measurement of where the ceiling actually comes from.

**The useful consequence: a level 0 acceptor whose fusion sites are better separated would
remove most of this cost while staying at level 0.** The notebook already accepts
`custom_destination`. Whether such a backbone is available is a wet-lab question.

## Decomposed: the introns are nearly free, the acceptor is not

For the 9S placement above, scoring the acceptor pair alone and then the whole reaction:

| acceptor | ends that meet | pair alone | with 3 introns |
|---|---|---:|---:|
| level −1 | `ACAT`/`ACAA` | 0.9364 | 0.9075 |
| **level 0** | `CTCA`/`CTCG` | **0.8282** | **0.8025** |
| level 1 | `GGAG`/`AGCG` | 1.0000 | 0.9691 |
| *(no acceptor)* | — | — | *0.9691* |

**The three intron junctions on their own score 0.9691.** They cost about 0.031. Everything
else — 0.17 at level 0 — is the acceptor's two ends being one base apart.

That reframes the question for the lab. It is not "can we afford introns"; it is "why are we
paying 0.17 for our backbone". PPROPT's own level −1 pair scores 0.9364, so a better-separated
pair is not exotic.

## The floor is deciding more than it should

Four of those five seeds score exactly 0.7507 against `FIDELITY_FLOOR = 0.75` — clearing it by
**seven ten-thousandths**. A floor of 0.751 would have dropped an intron from all four and
tripped the warning for each.

The threshold is a judgement call wearing the costume of a measurement. The result now carries
`margin_to_floor`, and a placement within 0.01 of the line is flagged rather than reported as a
clean pass, so the closeness is visible instead of buried in a constant.

## It works better outside Chlamydomonas, where the evidence is weaker

Same target and seed, the other two hosts with a spliceosome:

| host | CDS GC | route A clean | route B clean | 3 introns | evidence here? |
|---|---:|---:|---:|---:|---|
| *C. reinhardtii* nuclear | 0.606 | **1** of 65 | 16 of 65 | 0.803 | **yes** |
| *S. cerevisiae* | 0.402 | 28 of 51 | 43 of 51 | **0.828** | no |
| *A. thaliana* nuclear | 0.410 | 30 of 53 | 41 of 53 | **0.828** | no |

In both AT-rich hosts three introns cost **nothing** — 0.828 at one, two and three, identical
to the intron-free baseline. The trade-off this whole document is about simply disappears.

Chlamydomonas is the hard case, not the easy one: its GC-rich codon usage leaves one usable
embedded position in sixty-five, against roughly one in two elsewhere.

**And the caveat runs opposite to the numbers.** The hosts where placement is free are the
hosts where the enhancement is undocumented — rbcS2i1's 5.5-fold was measured in
*Chlamydomonas*, and Baier 2020 found most foreign introns failed in that alga until their
splice sites were adapted. Cheap to place is not the same as likely to work.

Two smaller findings from the same sweep. The **lower** GC bound never fails on its own in any
host: rbcS2i1 is GC 0.579, so a seam window is pulled *up* toward it, and in an AT-rich CDS
that jump is larger, not smaller (0.440 → 0.660 in yeast, against 0.640 → 0.700 in
Chlamydomonas). And a failure mode absent from Chlamydomonas appears in the AT-rich hosts:
homopolymer runs of 5, at two positions in yeast and three in Arabidopsis. Lowering GC trades
one rejection mode for another rather than removing rejection.

Every one of these junction failures is a **single-base miss** — on a 50-nt window the GC
quantum is 0.02, so 0.660 is 33/50 against a 0.65 ceiling.

## The stronger intron is not available

Baier 2020 ranked LHCBM1i2 (253 bp) highest of 33 native introns, about 8-fold against
rbcS2i1's 5.5. It **breaches the band on its own**: a 50-nt window at offset 14 reaches 0.680.
The breach is interior rather than at a wrapper seam, so no enzyme or pad choice moves it, and
an intron cannot be recoded. It is kept in the module as `LHCBM1I2` with a test, so the
question has an answer in the code.

rbcS2i1 is therefore not a preference. It is the only one of the two that fits.

## Two mistakes this went through, both caught by measurement

**Repeated overhangs.** The first placement took the first usable site past each spacing
target. Every candidate sits between two guanines, so that kept selecting the *same* four
bases, and identical overhangs let fragments swap. It scored 0.18 on 9S, 0.03 on 14S and
0.003 on 19S. It was briefly what the advisory recommended, and it would not have assembled.
It surfaced only because 14S and 19S were run; on 9S every ordered piece passed the synthesis
rules and the result looked fine.

**A fixed intron count.** The next version derived the count from the spacing rule —
`round(len(cds) / 240) - 1` — which is five introns on 14S and seven on 19S. Each extra
intron adds a junction to the same reaction and fidelity falls steeply, so that was picking
the worst end of the curve:

| | spacing-derived count | searched count |
|---|---|---|
| 14S | 5 introns, **0.517** | 3 introns, **0.751** |
| 19S | 7 introns, **0.307** | 4 introns, **0.686** |

The count is now searched and the whole curve reported. On 9S it reads
`{1: 0.828, 2: 0.828, 3: 0.803}` — one or two introns cost nothing at all, and the third
costs 0.026. Since every published series found more copies gave more expression, the
advisory takes the most that still clears `FIDELITY_FLOOR` (0.75) rather than the single
highest-scoring count.

## Making the overhang choosable

At a codon-aligned position the junction overhang can be picked by synonymous substitution,
the way cuts already are (`arelf.safe_overhangs`). Requiring codon alignment costs candidate
sites but buys options at each:

| | sites, overhang fixed | sites, overhang choosable |
|---|---:|---:|
| 9S | 16 | 34 |
| 14S | 25 | 57 |
| 19S | 32 | 75 |

Note this did **not** beat the fixed-overhang search at level 0 on 9S — 0.687 against 0.803 —
because requiring codon alignment removes positions that happened to carry good overhangs.
The two searches cover different spaces and neither dominates. The advisory uses the
fixed-overhang search, which scored better on all three.

## What a dash means in the raw results

"No valid overhang combination; all candidates conflict" — the solver's words, not a timeout.
On 9S at four introns the space is **eight** combinations and every one conflicts. Real
negatives.

But positions were fixed by even spacing in that sweep and only the overhangs optimised,
because searching both together is intractable: 19S at seven introns is choosing from ~75
sites. **Read the sweep table as a floor, not a ceiling.**

## What this does not establish

- **Nothing here is measured.** Fidelity is a prediction from Pryor's mis-ligation data, not
  an assembly anyone ran, and it says nothing about whether a spliced intron raises
  expression in a designed PPR. Nobody has tested that.
- **One seed, one host, one intron.** Seed 42, *C. reinhardtii* nuclear, rbcS2i1 throughout.
  LHCBM1i2 is the stronger intron in the literature and is untested against any of this.
- **Whether a better level 0 acceptor exists** for this lab, which is the single change that
  would matter most.
