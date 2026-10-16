# What PPROPT does

## The problem

Cells are full of RNA messages. We want to switch **one specific message** on or off without
touching the others. To do that you need a protein that grips only that message and ignores
everything else.

## The gripper

A **PPR protein** is built from repeats stacked in a row, each gripping one RNA letter — like
a zip tie cut to fit one cable and no other. Two amino acids inside each repeat decide which
letter it grips, and that mapping is known:

| RNA letter | repeat grips it with |
|---|---|
| A | TN |
| C | NN |
| G | TD |
| U | ND |

So the design is mechanical. Give us the RNA and we can say which repeats to stack.

## What the software does, step by step

**1. Read the target.** You type an RNA sequence, e.g. `AAAAUGUGG`. Its length picks the
size: 9, 14 or 19 letters → a 9S, 14S or 19S protein.

**2. Pick the repeats.** Look up each letter in the table above. Nine letters → nine repeats.

**3. Write the protein.** Stack the repeats into a full protein sequence, with the fixed
scaffold ends the published design needs.

**4. Turn the protein into DNA.** Every amino acid can be spelled by several three-letter DNA
*codons* — like synonyms. We pick spellings the host cell reads easily, while avoiding
patterns that cause trouble later.

**5. Cut it into fragments you can order.** Synthesis companies won't make one long piece reliably,
so we split the DNA into about four fragments, each with **sticky ends** — four unpaired
letters that only match their partner, like Lego studs, so the pieces snap together in one
order and no other.

**6. Check the pieces.** We check the exact string that would go on the order form, including
the fixed letters added at each end for joining. One check is **GC content** — how much of the
DNA is G or C. Too high and the synthesis machine struggles.

**7. Decide.** Either an order sheet, or a refusal naming the rule that broke. A refusal still
writes the files, labelled `REJECTED`, so you can inspect them.

## Two ways in

| Route | What it is | Status |
|---|---|---|
| **Kit route** | You already own the 42-part GRASP kit. We pick the parts and plan the pooling. Nothing new is made. | **Works today** — 42 items, one pool |
| **From-scratch route** | We design a brand-new gene for any target. | **Works today** — 200 of 200 saved designs orderable |

## Where we are now

Both routes produce an order sheet. The from-scratch route only started doing so on
2026-10-12, and the reason it used to refuse is worth keeping, because it is the most
instructive thing in the project.

Our rule says GC must stay between 0.35 and 0.65. The gene obeyed it. But what you order is
the gene fragment *plus* the GC-rich joining letters bolted on each end, and nothing was
checking that combined string while the codons were being chosen. On our example the gene sat
at 0.640 and the ordered piece reached **0.660** — over by 0.010. All 200 saved designs failed
that way, so the software refused all 200. We published the number rather than raising the
limit to make it pass.

## The fix, and what it cost

We changed one thing: show the optimiser the **whole ordered piece** instead of just the gene,
then let it pick different codon spellings for the same protein. The GC limit itself did not
move.

| | with the fix | without it |
|---|---:|---:|
| orderable | **200 of 200** | **0 of 200** |
| failing fragments | **0** | 277 |

Cost in codon quality: typically 0.003 on a 0–1 scale, worst 0.032, and thirteen designs got
slightly *better*. A second, independently written checker agreed on all 200.

The control matters as much as the result: with the one change removed, the software returned
**exactly** the old sequences — all 200, letter for letter. So the fix is what fixed it.

An earlier version of this run reported 8 failures and a reason for them. That was a bug in
our own test code, not a property of those targets. Once fixed, all 200 passed. The write-up
keeps the wrong claim on the record.

**It is switched on by default as of 2026-10-12.** Pass `wrapper_aware=False` to get the old
behaviour back, which is how the control above is run. Every number this document reports
comes from the version with it on; the earlier corpus is kept as
`work/independent/corpus200_prewrapper.jsonl` so the two can be compared directly.

## What we do not claim

- **No vendor said no.** We never submitted anything. Our own rule refused these, and the
  sequences may well be makeable.
- **Nothing biological is measured.** We have not tested whether any protein grips its RNA or
  whether a cell makes it. Every number is a prediction or a software check.
- **The GC limit is our own choice, not a vendor rule.** Nobody recorded where the number
  first came from; on 2026-10-13 we settled it as this project's own engineering decision. It
  happens to equal Twist's published codon-optimisation guidance, we enforce it harder than
  that document frames it, and no IDT document contains it at all — and IDT is who we order
  from. Defend it as conservatism, never as vendor acceptance.
- **The 200 designs are our own**, generated from random sequences by this software, frozen so
  every check measures the same set. They are not real biological targets.
- **Good codons are not enough for this host, and we only do codons.** In *Chlamydomonas*,
  a codon-optimised gene with no introns expressed barely at all; adding three copies of a
  short intron to the same gene raised mRNA more than 100-fold and protein about 8-fold
  (Baier 2018). Codon choice is the whole of what this software does for expression, so a
  design it calls good may still express poorly for a reason it does not model. The
  designer notebook now reports where introns could go; it does not put them there. See
  `docs/intron_feasibility.md`.
- **The kit route's pooling plan** is checked against vendor rules we copied by hand. That is
  not vendor approval.
