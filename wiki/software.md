# Software

> **Draft — figures not filled.** Every number in this page is left as a placeholder of the
> form `[[FIG: what it is | source]]`. Grep for `[[FIG:` to find them all. They were left out
> deliberately so each one can be checked against its source before it ships.

## What PPR-opt is

PPR-opt designs synthetic PPR proteins that bind a chosen RNA sequence, and turns that design
into DNA somebody can actually order. You give it an RNA target; it gives you back a protein,
the DNA that encodes it, the fragments to order, and a verdict on whether those fragments are
fit to be ordered at all.

The last part is the part we care most about, and it is the reason this page is longer than
"we wrote a sequence generator". The software's job is not to produce a sequence. It is to
produce a sequence **together with an honest statement of whether it should be used.**

## Why this needs software

A PPR protein is built from repeats stacked in a row, each repeat gripping exactly one RNA
base. Two amino acids inside each repeat decide which base it grips, and that mapping is
published:

| RNA base | repeat grips it with |
|---|---|
| A | TN |
| C | NN |
| G | TD |
| U | ND |

So the protein design itself is mechanical — a lookup, not a search. If that were the whole
job, a spreadsheet would do.

The work is everything downstream. The same protein can be encoded by an enormous number of
DNA sequences, and they are not interchangeable: some are hard to synthesise, some fold back
on themselves, some repeat internally in ways that make assembly unreliable, and some break
the rules of the company that would manufacture them. Choosing among those spellings, and
then proving the choice is sound, is where the software earns its place.

## How the design works

**1 · Read the target.** You enter an RNA sequence. Its length picks the architecture —
9, 14 or 19 bases give a 9S, 14S or 19S protein of
[[FIG: residue counts per architecture | src/ppropt/constants.py]] residues respectively.

**2 · Pick the repeats.** Each base is looked up in the table above.

**3 · Write the protein.** The repeats are stacked into a full protein sequence with the
fixed scaffold ends the published design requires.

**4 · Encode it as DNA.** Each amino acid has several possible codons. We choose spellings
the host reads efficiently while avoiding patterns that cause trouble later — forbidden
restriction sites, long internal repeats, extreme local GC.

**5 · Cut it into orderable fragments.** A single long synthetic fragment is unreliable, so
the DNA is split into roughly [[FIG: typical fragment count for 9S | README.md]] pieces, each
carrying four-base sticky ends that match only their intended neighbour. The pieces therefore
assemble in exactly one order.

**Cut positions are chosen before the codons, not after.** This ordering matters: optimising
first would let the optimiser rewrite the very bases the junctions depend on.

**6 · Check what would actually be ordered.** Not the gene — the exact string that goes on the
order form, including the enzyme wrapper and entry overhangs bolted onto each end. This
distinction is the single most consequential design decision in the project, and section
*What it refuses* below explains why.

**7 · Decide.** Either an order sheet, or a refusal naming the rule that broke.

## Two routes

| Route | What it does | Needs synthesis? |
|---|---|---|
| **Inventory** | Builds from the deposited GRASP module kit you already hold. Recodes modules for your host, or compiles a target from modules that exist. | Only for recoded modules |
| **De novo** | Designs a brand-new gene for any target you choose. | Yes, entirely |

The inventory route is constrained by what the kit can spell; the de novo route is constrained
by what a vendor will make. They are different workflows with different answers, so they are
different notebooks.

## What it refuses, and why that is the point

Both routes produce an order sheet today. The de novo route refused every design until
2026-10-12, and the reason is the clearest illustration of what this software is for.

Our synthesis policy requires local GC content to stay inside a fixed band. The coding
sequence obeyed it. What you *order* is that coding sequence plus the GC-rich joining bases
added at each end — and on the panel the ordered fragment crossed the limit while the gene
inside it did not, because nothing was showing the optimiser the combined string.
[[FIG: panel designs refused before the ordered-fragment constraint, and the margin by which
the worst one failed | PANEL_WRAPPER_AWARE_2026-10-11.md]]

The software caught its own rule breaking and refused. We published the number rather than
widening the limit until it passed, then fixed the optimiser instead — see
*The wrapper-aware constraint* in [[model]]. The band is unchanged; what changed is the
sequence the constraint is evaluated on. [[FIG: panel designs orderable after the change, and
the codon-quality cost | PANEL_WRAPPER_AWARE_2026-10-11.md §1 and §3]]

A refusal is still not silent and not destructive. The files are still written, every filename and
key marked `REJECTED`, so a design can be inspected — but no downstream step can mistake it
for an order.

## Running it

Both notebooks run in Colab with no local install:

- **Designer** — design a PPR for a target of your choosing
- **Inventory** — build from the deposited kit

Neither notebook ships stored outputs, so the figures appear only once you run it. The setup
cell compares the installed commit against `main` and updates itself when the notebook is
newer than the copy you have, because a stale install fails deep in a later cell in a way that
reads like a broken notebook.

## What the software does not do

These are limits, not caveats to be read past:

- **Nothing here is measured in a cell.** Every number the software reports is a prediction
  or a software check. We have not tested whether any designed protein binds its RNA, or
  whether a host expresses it.
- **No vendor has accepted or refused anything.** We never submitted a sequence. Passing our
  policy is not vendor acceptance, and the refusals the panel used to record were our own
  policy's, not a vendor's.
- **The GC band is our own engineering choice**, settled as such on 2026-10-13. It equals
  Twist's published codon-optimisation guidance numerically, we enforce it harder than that
  document frames it, and no IDT document contains it — and IDT is who we order from.
- **Predicted fidelity is not assembly efficiency.** It comes from published ligation-count
  matrices, not from a reaction in anyone's hands.
- **The QC check is sequence complexity**, not a calibration against vendor outcomes.
- **The panel is synthetic.** The designs it covers were generated by this software from
  random sequences and frozen so that every experiment measures the same set. They are not
  biological targets.
- **Codon optimisation is the only expression lever here, and it is documented as
  insufficient for the default host.** A codon-optimised, intron-free gene expressed
  minimally in *Chlamydomonas*; three copies of *rbcS2* intron 1 in the same gene gave over
  100-fold the mRNA (Baier 2018). PPR-opt designs no introns. It reports candidate positions
  and says what it has not accounted for.

## Reproducibility

Every figure the project quotes is reproducible from a pinned environment, and the release
process re-derives rather than re-reads: each piece of recorded evidence states the artefacts
it was produced from, and the checker recomputes those identities rather than trusting the
record. A record that names no source, or whose stated hash no longer matches, fails.

[[FIG: test count, release-check count, and evidence bundle identity | RELEASE_STATUS, README]]
