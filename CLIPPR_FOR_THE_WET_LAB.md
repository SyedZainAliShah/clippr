# CLIPPR — everything you need to present it

A briefing for the dry lab lead, written to be read once before a meeting with biologists
who have not seen the software. It covers what CLIPPR does, the biology behind it, every
number it reports and what that number means, a walkthrough of the notebook, what we can
honestly claim today, and a pitch structure with the questions you will be asked.

Nothing here assumes the reader has used the tool.

---

# Part 1 — The pitch

## The one-sentence version

**CLIPPR designs the DNA for a custom RNA-binding protein, given only the RNA sequence you
want it to grab, and checks the design well enough that you can order it.**

## The three-minute version, in the order to say it

**1. The problem.** PPR proteins are the closest thing biology has to a programmable
RNA-binding module. Each repeat in the protein grips exactly one RNA base, and which base it
grips is set by just two amino acids. So in principle you can write down any RNA target and
read off the protein that binds it.

In practice you cannot order a protein. You have to order DNA, and the DNA is where it gets
hard:

- a 19-repeat PPR is roughly 2,000 bases of coding sequence
- it is *tandem repeats of one template*, so the DNA is highly repetitive and synthesis
  companies reject or mis-synthesise repetitive sequence
- it is too long for one oligo, so it has to be split into fragments and re-joined
- the joins are Golden Gate junctions, and if two junctions have similar sticky ends the
  parts assemble in the wrong order
- and all of that has to be true *at once*, in the same sequence

Doing this by hand is a day of work per target and mistakes are invisible until the cloning
fails.

**2. What CLIPPR is.** A design tool that takes the RNA target and returns fragments ready to
order, having checked all of those constraints simultaneously. It runs in a Colab notebook —
no install, no command line.

**3. What makes it worth using rather than doing by hand.** Three things:

- **It never reports success on something it has not checked.** Every sequence it exports is
  digested *in silico* and re-assembled to confirm the fragments actually produce the intended
  coding sequence. A design that fails is reported as a failure with the reason.
- **It shows the trade-offs instead of hiding them.** There is no single best design. Better
  codon usage often means worse GC content; better assembly fidelity often means worse codon
  usage. CLIPPR gives you the set of designs where you cannot improve one thing without
  worsening another, and you pick.
- **Every number it prints is labelled as predicted or measured.** It will tell you the
  predicted ligation fidelity; it will never tell you the protein works.

**4. The ask.** We need wet lab people to (a) sanity-check that the outputs are the shape you
actually want to order, (b) tell us which constraints matter at your bench and which are us
being over-cautious, and (c) eventually build one and tell us what happened.

## What *not* to claim in the room

Say these out loud before anyone asks — it buys credibility, and all three are true:

- **No sequence designed by this tool has been synthesised or tested.** Everything is
  computational.
- **CAI does not predict expression.** It predicts how well the codons match the host's
  preferences. It is a proxy, and a rough one.
- **Ligation fidelity is a prediction from a published enzyme dataset**, not a measurement of
  your reaction, your buffer or your thermocycler.

A biologist who hears a dry lab tool claim it "guarantees expression" stops listening. One who
hears clearly-labelled predictions engages.

---

# Part 2 — The biology, briefly

## The PPR code

A PPR protein is a chain of ~35-amino-acid repeats. Each repeat contacts one RNA base. Two
positions within the repeat — commonly numbered 5 and 35 — determine which base:

| RNA base you want bound | amino acid pair |
|---|---|
| A | **TN** |
| C | **NN** |
| G | **TD** |
| U | **ND** |

That is the whole code. A 14-base target needs 14 repeats, each carrying the pair for its base.
Everything else in the repeat is the same scaffold.

**Consequence for the DNA:** the protein is 14 copies of one sequence differing at two
positions. Translated back to DNA naively, that is 14 near-identical stretches — exactly what
synthesis companies refuse. Most of what CLIPPR does is choosing *synonymous codons* to break
that repetitiveness without changing a single amino acid.

## Architectures

CLIPPR supports three, and the architecture follows from the target length — it is not a
choice you make:

| target length | architecture | approximate protein |
|---|---|---|
| 9 bases | **9S** | ~300 aa |
| 14 bases | **14S** | ~470 aa |
| 19 bases | **19S** | ~660 aa |

## Golden Gate assembly, in one paragraph

Type IIS restriction enzymes cut *outside* their recognition site, so you can choose the
four-base sticky end they leave. Give every fragment ends that match only its intended
neighbour and the whole construct assembles in one tube, in order, with the enzyme sites
destroyed in the product. The enzymes CLIPPR works with:

| enzyme | recognition site |
|---|---|
| BsaI | `GGTCTC` |
| BbsI (= BpiI) | `GAAGAC` |
| SapI | `GCTCTTC` |

**The critical constraint:** these sites must not appear *anywhere else* in your sequence, or
the enzyme cuts your construct in the middle. Every codon choice CLIPPR makes is checked
against this.

## Assembly levels

The parts system builds up in stages. You will hear these numbers:

| level | enzyme | what it makes |
|---|---|---|
| level −1 | BsaI | the entry clone — how a module plasmid is *made* |
| level 0 | BbsI | releases one module for assembly — how a module is *used* |
| level 1 | BsaI | joins modules into a finished cassette |

---

# Part 3 — Every number CLIPPR reports, and what it means

This is the section to have open during the meeting.

## The two headline metrics

### CAI — Codon Adaptation Index

**Range 0 to 1. Higher is better. What it measures:** how closely the codons used match the
codons the host organism actually prefers, for the same protein.

**Why it matters:** a protein whose codons are rare in the host can translate slowly or stall.

**What it is not:** a prediction of expression level. Two sequences with the same CAI can
express very differently. Treat it as "nothing obviously wrong with the codon usage" rather
than "this will express well."

**Numbers you can quote:** recoding the 42-module kit for *Chlamydomonas* moves mean CAI from
**0.208 to 0.621** — the deposited sequences were not written for this host.

### Ligation fidelity

**Range 0 to 1. Higher is better. What it measures:** the predicted fraction of assemblies that
join in the correct order, given the set of four-base sticky ends in the reaction.

**How it is computed:** from published measurements of how often each overhang pair ligates,
including the mismatched pairs. If two junctions in your reaction have similar ends, they
compete, and fidelity falls.

**Why it matters:** this is the number that predicts whether your cloning works. A reaction at
0.99 is essentially clean; one at 0.75 means a meaningful fraction of colonies are wrong.

**It is a property of the whole reaction, not of one junction.** Adding a part changes the
fidelity of every other junction.

**Numbers you can quote:**

| reaction | predicted fidelity |
|---|---|
| level-0, six parts (BbsI) | **0.752** |
| 9S block joins (BsaI) | **0.998** |
| 14S block joins (BsaI) | **0.998** |
| 19S block joins (BsaI) | **0.996** |

A worked single design: a 9S target gives **302 aa, 906 nt, 4 fragments, fidelity 0.828**.

## Sequence quality metrics

### GC content

Measured two ways, and both matter:

- **Global GC** — across the whole sequence
- **Windowed GC** — in every sliding **50-base** window

A sequence can have a perfect global GC and still contain a 50-base stretch at 80% GC that
fails synthesis. The window is the one that actually bites.

**Thresholds CLIPPR uses for its QC verdict:**

| | warn below | warn above | fail below | fail above |
|---|---|---|---|---|
| global GC | 35% | 65% | 25% | 75% |
| any 50 nt window | 25% | 75% | 15% | 85% |

### Homopolymer run

The longest run of a single base (`AAAAAA`). Long runs cause polymerase slippage and
synthesis errors. **Warn at 6, fail at 9.**

### Repeat fraction

How much of the sequence is duplicated elsewhere in the same sequence. **Warn above 10%, fail
above 50%**, and the longest single repeat warns at 25 bases and fails at 50.

This is the metric the PPR architecture fights hardest — the protein *is* a repeat.

### QC verdict

Every design gets **PASS / WARNING / FAIL** from the thresholds above, plus every underlying
measurement, so you can rank designs on the numbers rather than on the label.

## Search and trade-off metrics

When CLIPPR explores alternative designs it scores each on three axes at once:

| axis | direction | meaning |
|---|---|---|
| **fidelity** | maximise | predicted assembly correctness, as above |
| **adaptation** | maximise | mean CAI across the whole inventory |
| **synthesis_fitness** | maximise | composite of GC centrality, homopolymer headroom, internal repetition and cross-module sharing |

It returns the **observed Pareto front**: every design where you cannot improve one axis
without losing another. Typically 9 to 15 designs. You pick one, and the tool records which
and why.

**Say "observed", not "optimal".** It is the best set among the designs actually evaluated,
not a proof that nothing better exists.

## The synthesis policy — the one real decision

Every stage that judges a sequence reads one policy. It sets the GC band, the homopolymer
limit, and — importantly — whether breaking a rule **refuses** a design or is **reported and
delivered**.

| policy | 50 nt GC window | a breach is |
|---|---|---|
| `clippr-strict-legacy` *(default)* | 0.35–0.65 | refused |
| `strict-as-target` | 0.35–0.65 | reported, sequence delivered |
| `broad-experimental` | 0.15–0.85 | refused |

**This matters and it is worth raising in the meeting.** Measured on the 42-module kit:

- strict band as a hard limit → mean CAI **0.619**
- strict band as a *target* → mean CAI **0.819**, with **40 advisory warnings**

Same thresholds. The difference is whether the tool refuses the sequence or hands it over with
a flag. Enforcing the band as guidance rather than a wall is worth **+0.20 mean CAI**, and you
still get told which 40 modules sit outside the conservative band.

**The question for the biologists in the room:** would you rather the tool refuse to give you a
sequence outside the conservative GC band, or give it to you with a warning attached? That is a
risk call about your synthesis vendor, and it is theirs to make, not ours.

## Ordering constraints

The export targets IDT oPools:

| | |
|---|---|
| oligo length | 40–350 nt |
| oligos per pool | 2–384 |

A design is split into fragments that fit, then packed into pools. Cost is estimated only when
a dated price list is supplied — otherwise it reports "unavailable" rather than guessing.

---

# Part 4 — The notebook, cell by cell

Open it in Colab. Nothing is installed locally. There are **two routes** and they answer
different questions.

## Route A — design a new protein from scratch

Use when you want a binder for a target and you do not have parts.

| cell | what it does | what to look at |
|---|---|---|
| **Setup** | installs CLIPPR | takes ~1 min, once |
| **Design parameters** | your RNA target, host organism, enzyme set | the only cell you normally edit |
| **Upload a codon table** *(optional)* | supply your own host codon usage | skip it and a *Chlamydomonas* table is used |
| **Design** | runs the whole design | the summary line: length, fragments, fidelity, QC |
| **Fragment table** | the actual oligos | this is what gets ordered |
| **Design audit** | why *this* design and not another | the constraint that was binding |
| **Download the design files** | FASTA / CSV / GenBank | take these to the bench |
| **Off-target check** | does this RNA sequence already occur in the host? | a hit here matters biologically |
| **Kit route** | can you build this from parts you already own? | cheaper if yes |
| **Explore the alternatives** | the trade-off front for this one target | pick and justify |
| **Design the whole library** | many targets at once | |
| **Library QC table** | QC verdict per member | scan for FAIL |
| **Homology** | DNA shared between library members | high sharing risks recombination |
| **Two-tier cross-talk** | will your binders bind each other's targets? | separation gates, affinity annotates |

## Route B — reuse a module kit

Use when you have the 42-module parts kit and want to *compile* targets from it rather than
synthesise new DNA. Far cheaper: the DNA is ordered once and reused.

| cell | what it does | what to look at |
|---|---|---|
| **Inventory 1 — load the deposited kit** | reads the module table | 42 modules; any structural problems |
| **Inventory 1b — choose the synthesis policy** | the decision from Part 3 | prints the policy and its version |
| **Inventory 2 — recode for your host** | rewrites every module's codons, **protein and both junctions frozen** | mean CAI before → after |
| **Inventory 3 — compile targets** | which targets these modules can spell | a target may simply not be buildable |
| **Inventory 4 — optimise the collection** | improves the kit *as a set*, not module by module | |
| **Inventory 5 — explore junction trade-offs** | the Pareto front over junction choices | front size, what each point trades |
| **Inventory 6 — commit to a choice** | rebuilds the inventory around your pick | **this changes the interface version** |
| **Inventory 7 — eligibility and pool plan** | can this be ordered, and in how many pools | |
| **Inventory 8 — download the package** | everything, with provenance | |

### One thing to flag at Inventory 6

Choosing a different junction assignment produces a **different interface version**. Modules
from the old version and the new version **cannot be mixed** in an assembly. The notebook says
so at the point of choosing rather than letting it be discovered at the bench.

### Runtime, so nobody thinks it has crashed

| | |
|---|---|
| whole inventory route (Route B), start to finish | **4.7 seconds** |
| same route under the permissive policy | **6.3 seconds** |
| a single 9S design | **0.5 seconds** |
| a single 14S design | **14 seconds** |
| a single 19S design | **3.4 minutes** |

The 19S design is slow because it is solving ~2,000 bases against GC, homopolymer and three
forbidden-site constraints simultaneously, and the cost grows steeply with target length —
0.5 s, 14 s, 205 s for 9, 14 and 19 bases. **It is working, not hung.** Say this before the
demo if you are going to run a 19S live; three minutes of a blank cell in front of an audience
is otherwise uncomfortable. Better: demo a 9S.

*Scope of those design figures:* codon design and fragmenting only, with the off-target genome
scan switched off. The notebook has off-target screening **on** by default, so a real run is
somewhat slower, and the very first run on a fresh machine also downloads the host genome
once and caches it.

---

# Part 5 — What we are achieving today

Everything here is measured and reproducible, not estimated.

## The pipeline is complete end to end

Target → PPR design → codon optimisation → fragmenting → junction selection → trade-off front →
selection → assembly-ready substrates → order package. No stubs in the main route.

## It checks itself

- **949 automated tests, all passing.**
- **A 10-point release check passes 10/10**, including building the package in a *clean, empty
  environment* and reproducing the documented example output character for character.
- Exported sequences are **digested and re-assembled in silico** and compared with the intended
  coding sequence — and deliberately corrupted controls are caught, so we know the check works
  rather than merely passing.
- 200 stored designs across all three architectures (100× 9S, 50× 14S, 50× 19S) reconstruct
  correctly, independently verified by code that shares nothing with the code it checks.

## The specific goals it currently meets

| goal | status |
|---|---|
| Design a binder for any 9/14/19-base RNA target | **working** |
| Never export a sequence that fails its own contract | **working** — refusals are reported with reasons |
| Recode a module kit for a new host without changing protein or junctions | **working** — CAI 0.208 → 0.621 |
| Predict assembly fidelity per reaction, not per junction | **working** — level-0 0.752, block joins 0.996–0.998 |
| Show trade-offs rather than a single answer | **working** — observed Pareto front, 9–15 designs |
| Make the synthesis rules explicit and versioned | **working** — three policies, recorded on every artefact |
| Check targets against the host's own genome | **working** |
| Check library members against each other | **working** — homology and two-tier cross-talk |
| Produce an orderable package with provenance | **working** — IDT oPools, pooled, hash-verified |
| Run without installing anything | **working** — Colab, 22/22 cells execute |

## Honest limits, stated plainly

| | |
|---|---|
| Nothing has been synthesised or tested at a bench | the entire project is computational |
| CAI does not predict expression | it predicts codon match, nothing more |
| Fidelity is predicted from published enzyme data | not a measurement of your reaction |
| The Pareto front is *observed*, not proven optimal | it is the best among designs evaluated |
| Level-1 reactions cannot be fully scored yet | the linker and editing-domain parts are not in any document we have — the tool refuses to score rather than guessing |
| Vendor acceptance is not established | passing our checks is not the same as a company agreeing to synthesise it |

That last table is not a weakness to hide. It is the reason the tool's positive claims are
worth believing.

---

# Part 6 — Questions you will be asked

**"How do I know the sequence it gives me is right?"**
Every exported sequence is cut and re-joined in software and compared with the intended coding
sequence. Deliberately broken sequences are fed in as controls and are caught. That is a check
on the checker, not just on the design.

**"What if I don't like the codon table it used?"**
Upload your own. One cell, one file.

**"Can it handle my target?"**
If it is 9, 14 or 19 bases, yes. Other lengths are refused with a message rather than silently
rounded — one PPR repeat binds one base, so a 12-base target is not a 12-repeat protein you can
half-build.

**"Why did it refuse my design?"**
It will tell you which constraint failed and by how much. The commonest is a 50-base GC window
outside the band. That is the policy decision in Part 3 — you may want the permissive setting.

**"Will this protein actually bind my RNA?"**
We do not know, and the tool does not claim to. The PPR code is well established, but predicted
recognition is not measured affinity. This is exactly what we want the wet lab to find out.

**"How much does it cost to order?"**
It estimates only when given a dated price list. Otherwise it reports the pool plan and says
cost is unavailable, rather than quoting a number that might be a year stale.

**"What do you need from us?"**
Three things: check the output is the shape you'd actually order; tell us which constraints are
real at your bench and which are us being over-cautious; and if you build one, tell us what
happened — good or bad.

---

## The closing line

*"It is a design tool that refuses to lie to you. Every number is labelled predicted or
measured, every refusal comes with a reason, and everything it claims to have checked, it has
actually run. What it cannot tell you is whether the protein works — and that is where we need
you."*
