# Should PPR-opt insert introns?

*Written 2026-10-16, against commit `30cf908`. A feasibility assessment, not a decision and
not an implementation. Every measurement in §4 was taken from this checkout; every number in
§2 and §3 is from the cited literature and was measured by someone else, in someone else's
construct, in an organism that is sometimes not ours.*

> **Partly superseded, 2026-10-16.** This document established *whether* an intron can be
> ordered. It did not measure what placing one costs in assembly fidelity, and that turned out
> to dominate: the naive placement it implies scores 0.18 on 9S and 0.003 on 19S, because
> every candidate position sits between two guanines and the junction overhangs repeat.
> `docs/intron_fidelity.md` has the measurements, the corrected placement, and the finding
> that the acceptor vector rather than the introns is what caps the result. Read that one for
> the current numbers; this one for how the question was framed.

> **Independently verified 2026-10-16.** 37 claims were re-checked by a separate reader that
> did not write this document: every `file:line`, every measurement re-derived from scratch,
> and every published figure against the primary paper rather than a search result. 34 came
> back correct, two were wrong and one was over-reached. All three are corrected here: a
> citation off by one line, a table cell carrying the 19S value in the 14S column, and a claim
> that GRASP uses BpiI at level −1, which its paper does not state. Two framings that were
> technically true but would have misled a reader were also tightened. The corrections are in
> place; this note records that they were needed.

---

## 0. The short version

PPR-opt's default host is *Chlamydomonas reinhardtii*, and the published work on that organism
says plainly that **codon optimisation alone does not produce expression**. Codon
optimisation is the whole of what PPR-opt does for expression. The intervention that closes
the gap in the literature is intron insertion, and it is worth between 4- and 100-fold
depending on what is measured.

Three findings decide the feasibility:

1. **The canonical intron passes every PPR-opt synthesis rule unchanged.** I expected this to
   be the blocker. It is not. *rbcS2* intron 1 sits at GC 0.579 with no 50-nt window outside
   0.35–0.65, a longest homopolymer of exactly 4, and no BsaI, BbsI or SapI site.
2. **Insertion sites are abundant in PPR-opt's own output** — 65 of them in a 906 nt design,
   one every 13.9 nt — because the published placement rule rewards the GC-rich sequence that
   Chlamydomonas codon optimisation already produces.
3. **The binding constraint is oligo length, not biology.** IDT caps oPools at 350 nt and
   PPR-opt's fragments leave 53–86 nt of headroom. The full 145 nt intron does not fit inside
   a fragment on any architecture.

The honest verdict is **staged**: one thing worth recording today at zero risk, one advisory
feature worth building at low risk, and a real implementation that is expensive and should
not be attempted before the current release is behind us. §8 gives the detail.

---

## 1. Why the question exists at all

PPR-opt designs a PPR coding sequence, picks codons the host reads easily, splits the DNA into
orderable fragments and checks what would go on the order form. Nothing in it addresses
whether the resulting gene is *expressed*. That was a defensible scope while the project
claimed only orderable DNA — `src/ppropt/parts.py:42` says so explicitly: the route "does not
model the acceptor plasmids, the non-PPR elements of the transcriptional unit".

Two things make it less defensible now.

**The host is the hard case.** Nuclear transgene expression in *C. reinhardtii* is poor
enough to have its own literature. Transgenic DNA is selectively silenced by an epigenetic
pathway ([Nature Communications 2020](https://www.nature.com/articles/s41467-020-19983-4)),
which is why the field standardised on the UV-mutagenesis strains UVM4 and UVM11 — both
carrying mutations in the Sir2-type histone deacetylase *SRTA* that leave transgene chromatin
more open ([Neupert et al. 2009](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-313X.2008.03746.x)).
Even in those strains, nuclear expression remains below chloroplast expression and below
other eukaryotic platforms.

**The thing we are designing is big.** PPR-opt emits 302, 457 or 612 aa proteins — 906 to
1836 nt of CDS. That is squarely the size class the intron literature is about.

And the kit we build from has never been tested in this organism. The GRASP toolkit was
validated in *E. coli* cell-free lysate (NEBexpress), and its authors are careful to say the
assay "may not fully reflect protein activity in vivo"
([Nucleic Acids Research 2025, gkaf1169](https://academic.oup.com/nar/article/53/20/gkaf1169/8321212)).
No published study expresses a synthetic PPR in *Chlamydomonas* and tests RNA binding there.
The closest is nuclear expression of *natural* PPRs — maize PPR10 and Arabidopsis HCF107,
each fused to a transit peptide, activating chloroplast reporters carrying the matching 5′UTR
([Algal Research 2021](https://www.sciencedirect.com/science/article/abs/pii/S2211926421003544)).

So the expression question is not a detail someone else has already solved for us.

---

## 2. The evidence that introns matter

### 2.1 The result that should worry us

[Baier et al. 2018, *Nucleic Acids Research* 46(13):6909](https://academic.oup.com/nar/article/46/13/6909/5042046)
built a codon-optimised, intron-free 1,662 bp patchoulol synthase gene for *C. reinhardtii*.
It expressed minimally. The authors then placed three copies of *rbcS2* intron 1 into the same
CDS, roughly every 500 bp, adding about 435 bp (25% more DNA). Against the intron-free
control:

| measure | improvement |
|---|---|
| mRNA abundance | **more than 100-fold** |
| protein, mean fluorescence per cell | **about 8-fold** |
| protein titre achieved | 0.03–0.05% of total soluble protein |

mRNA stability measured with actinomycin D did not differ between constructs. Each stepwise
intron addition raised expression, and a fourth copy improved on the third.

**The control is the part that matters to us.** Their intron-free arm *was* codon-optimised.
Codon optimisation on its own was not sufficient. That is a direct statement about the
limits of what PPR-opt currently does.

### 2.2 How large the effect is, and which intron

[Baier et al. 2020, *PLoS Genetics* 16(8):e1008944](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008944)
screened 33 native and 13 non-native introns:

| intron | length | enhancement |
|---|---:|---:|
| LHCBM1i2 | 253 bp | **~8-fold** (best native tested) |
| RBCS2i1 | 145 bp | **~5.5-fold** |
| RBCS2i1, truncated ΔE3+4 | 115 bp | no loss |
| RBCS2i1, truncated ΔE1+3+4+5+6 | 70 bp | **~3.8-fold** |

There was no correlation between intron size and enhancement. Deletion analysis found no
required element: the authors describe the active sequence as "dispersed and partially
redundant", and neither multiple sequence alignment nor MEME motif analysis found a conserved
motif.

[The MoInClo system (PMC11925875)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11925875/)
combined two introns in one gene and reported **more than 16-fold** NanoLuc activity over the
intron-free control, with about 60-fold mRNA in the top pools. Fewer than 5% of intron-free
transformants gave any signal at all, against roughly one third for its two best
single-intron constructs — a third single-intron design managed only about 10% — so introns
change the *hit rate*, not only the ceiling.

The effect was first described by
[Lumbreras et al. 1998](https://onlinelibrary.wiley.com/doi/10.1046/j.1365-313X.1998.00145.x),
who reported that *rbcS2* intron 1 raised both transformation frequency and *ble* expression.
An enhancer element was proposed at the time; Baier 2020 showed there is none.

### 2.3 The placement rule, which is quantitative

Baier 2018 built a splicing reporter testing all 16 combinations of the nucleotides flanking
the insertion point. The result is the design rule, and it is sharp:

| flanking context | colony yield vs control |
|---|---:|
| **G \| G** | **132%** |
| purine-flanked (A/G) | 86–106% |
| pyrimidine-flanked | 15–16% |
| C \| C | **3%** |

So the rule is **NG/GN — put the intron between two guanines**. Nucleotides beyond the
immediate pair showed no correlation. The authors note, usefully for us, that high-GC
codon-optimised sequence contains many such sites.

Spacing follows the native architecture: exons of **≤500 bp**, against a *C. reinhardtii*
average exon length of **240 nt**. Worth noting that the paper's own best construct has a
583 bp exon, so its guideline is a guideline and not a threshold even to its authors. A construct with a 1,057 bp exon expressed worse than the
evenly-spaced design. MoInClo's own constructs used exons of 234–432 bp.

Position matters as well as spacing: enhancement is stronger closer to the transcription
start site, and removing the N-terminal introns lowered expression more than removing others.

### 2.4 The result that reframes the problem

[Nucleic Acids Research 2025, 53(4):gkaf097](https://academic.oup.com/nar/article/53/4/gkaf097/8039171)
made over 300,000 introns carrying **160 nt of fully random internal sequence** between
constant splice signals, integrated them singly into the genome of human HEK293T cells, and
analysed the **19,346** that met its minimum read-count threshold.

Nearly every one enhanced expression. The library averaged **8.5-fold mRNA and 8.9-fold
protein**, with about 95% of introns falling between roughly 6.5- and 11-fold. More than 98%
of GFP reads were spliced as expected despite the random interiors. The authors' conclusion
is that enhancement is "a general property of splicing, largely independent of intron
sequence", correlating only weakly with splicing efficiency (Spearman ρ = 0.072) and with
intronic poly-uridine content.

If that transferred to *Chlamydomonas*, it would change the engineering problem completely:
we would not need a particular natural intron, we could **design** one to satisfy PPR-opt's own
GC and homopolymer rules.

**It does not obviously transfer, and this is the most important caveat in this document.**
That study is human cells. Baier 2020 tested 13 non-native introns in *Chlamydomonas* and
found that most failed — the authors attribute this to foreign introns coordinating poorly
with an algal spliceosome evolved on native boundaries. Adapting the internal splice sites to
the algal consensus rescued 7 of the 13: an SV40 intron went from not splicing at all to
about 7-fold, and a synthetic murine intron exceeded RBCS2i1. Four gained nothing and two
were made worse. The Arabidopsis UBQ10 intron, a workhorse elsewhere, did not enhance at all.

So the algal spliceosome is pickier than the human one, and "any well-spliced intron works"
is a hypothesis for this host, not an established fact.

### 2.5 Mechanism, and the honest state of it

IME is not an enhancer effect. Baier 2020 showed that placing the intron outside the
transcribed region abolished the effect, so it must be transcribed and spliced, and concluded
the mechanism is post-transcriptional. Beyond that the mechanism is open. Candidate
explanations across the wider literature include exon junction complex deposition aiding
export, effects on transcript stability, and effects on translation efficiency
([Frontiers in Plant Science 2011](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2011.00098/full);
[Shaul 2017](https://www.sciencedirect.com/science/article/abs/pii/S1357272517301541)).
Plant studies point to U-rich motifs, which the 2025 random-intron work independently
recovers as poly-U content. A TREX-complex route is documented in human and viral systems but
has not been shown in plants.

For our purposes the mechanism being unknown is not disqualifying — the effect is
reproducible and large — but it does mean there is no model to design against beyond "splice
well, and put it in the right place".

### 2.6 There is precedent for a tool doing this

PPR-opt would not be inventing the category. The
[Wild Worm Codon Adapter](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8496300/)
is a published web tool that codon-optimises a transgene for non-*Caenorhabditis* nematodes
and then inserts introns automatically, dividing the optimised cDNA at three predicted
insertion sites spaced approximately equidistantly, with built-in GT…AG introns and an option
for user-supplied ones.

That is almost exactly the feature under discussion, for a different clade. I found no
general-purpose gene design tool that does it beyond that one; the mainstream codon
optimisers (IDT, GenScript, Eurofins GENEius) do not insert introns, though GENEius does
screen for *artificial* splice sites — the opposite concern.

---

## 3. Why it might not belong in PPR-opt

### 3.1 It is outside the declared scope

PPR-opt claims orderable DNA and refuses to claim expression. Adding introns is an expression
intervention, and `parts.py:42` already puts the non-PPR elements of the transcriptional unit
out of scope. The counter-argument is that an intron sits *inside* the PPR coding sequence,
not beside it — it is not a promoter or a terminator — but the boundary is genuinely arguable.

### 3.2 It would be a second unvalidated layer

Nothing in PPR-opt is validated against a bench result. Every figure it reports is a
prediction or a software check. Intron placement would add a second predictive layer on top
of the first, and the evidence for it is from other people's constructs.

### 3.3 It applies to at most three of six hosts, and is catastrophic in the other three

`e_coli`, `c_reinhardtii_chloroplast` and `n_tabacum_chloroplast` have no nuclear
spliceosome. A spliceosomal intron in those contexts is not merely inert — it is read as
coding sequence. Measured on the actual sequence:

- rbcS2i1 is **145 nt, and 145 mod 3 = 1**, so an unspliced copy **frameshifts** everything
  downstream;
- it carries a **stop codon in all three reading frames**, so the protein is truncated as
  well as frameshifted.

Any intron feature must refuse those three hosts outright rather than warn. This is the one
place where a quiet default would destroy a construct.

Of the remaining three, the evidence is *Chlamydomonas*-specific. Nothing supports carrying
rbcS2i1 into *S. cerevisiae* or *A. thaliana*, and the UBQ10 result above is a caution
against assuming introns travel.

### 3.4 The kit route cannot use it at all

`recoding.py:203` rebuilds a module as `dna[:start] + result["cds"] + dna[end:]`, where
`result["cds"]` comes from `optimize_cds` and is therefore exactly `3 × len(protein)`. An
intron inside a recoded span would be **silently deleted**, and the length guard at
`recoding.py:206` would then reject the candidate without explaining why. The deposited GRASP
modules are themselves intron-free fragments of one continuous ORF, spliced at four-base
overlaps.

So this is a de novo-route feature only.

---

## 4. What I measured against this checkout

### 4.1 The intron passes every synthesis rule, unchanged

*rbcS2* intron 1, NCBI X04472.1 bases 1307–1451, 145 bp, retrieved and scored against
`STRICT_LEGACY`:

| check | PPR-opt rule | rbcS2i1 | verdict |
|---|---|---|---|
| overall GC | — | 0.5793 | — |
| 50-nt window GC | 0.35–0.65, hard | min 0.50, max **0.64** | **0 of 96 windows breach** |
| longest homopolymer | ≤ 4 | **4** | passes, no margin |
| BsaI `GGTCTC` | none, either strand | 0 | passes |
| BbsI `GAAGAC` | none, either strand | 0 | passes |
| SapI `GCTCTTC` | none, either strand | 0 | passes |
| splice boundaries | — | starts `GT`, ends `AG` | canonical |

This is the headline, and it is worth stating plainly: **the canonical intron is drop-in
compatible with the synthesis policy we already enforce.** Note the two zero-margin results.
The worst window is 0.64 against a 0.65 ceiling, and the homopolymer is exactly at the limit.
Both would need re-checking in context, because a window straddling an exon/intron boundary
sees GC-rich coding sequence on one side — structurally the same problem the wrapper-aware
constraint already solves for the ordered fragment, and the same 0.01-margin failure that
used to refuse every de novo design.

### 4.2 Insertion sites are abundant

Applying Baier's G|G rule to a real 9S *C. reinhardtii* design produced by this checkout
(target `AAAAUGUGG`, seed 42):

- **65 G|G sites in 906 nt — one every 13.9 nt.**
- Under the 240 nt native exon spacing, three introns would fall at nt 249, 489 and 739.
- Under the 500 nt rule, one intron at nt 516.

Site scarcity is not a constraint. This is the mechanism Baier predicted: high-GC codon
optimisation generates the required context.

### 4.3 Phase is free

rbcS2i1 is natively a **phase-2** intron — 164 coding nt precede it (1143→1306), and
164 mod 3 = 2. Baier inserted it at NG/GN sites without regard to codon phase. Spliceosomal
introns occur in all three phases.

This matters because PPR-opt's cuts are codon-aligned by construction
(`assembly_spec.py:260`, established by measurement over 750 cut contexts). Phase-0 insertion
at a cut is therefore legitimate, and we are not fighting the existing geometry.

The existing cuts are not *already* at G|G — in the design above, cuts at nt 225, 453 and 678
have contexts `C|G`, `G|A` and `G|C`. None qualifies. But the codon optimiser chooses those
bases, so "make the cut context G|G" is a constraint it could satisfy, of exactly the same
class as the four-base overhang locks it already honours.

### 4.4 Oligo length is the binding constraint

IDT oPools cap at **350 nt** (`ordering.py:89`). Current oligos:

| arch | CDS | oligos | oligo length | headroom |
|---|---:|---:|---|---:|
| 9S | 906 nt | 4 | 257–264 nt | 86 nt |
| 14S | 1371 nt | 6 | 257–264 nt | 86 nt |
| 19S | 1836 nt | 7 | 293–297 nt | **53 nt** |

Adding an intron to a fragment's payload:

| intron | size | 9S | 14S | 19S |
|---|---:|---|---|---|
| rbcS2i1, full | 145 nt | 409 ✗ | 409 ✗ | 442 ✗ |
| rbcS2i1, truncated | 70 nt | 334 ✓ | 334 ✓ | 367 ✗ |
| LHCBM1i2, strongest | 253 nt | 517 ✗ | 517 ✗ | 550 ✗ |

**Only the truncated intron fits, and only on 9S and 14S.** This single table is what forces
the choice in §5.

A related finding worth recording on its own: **`design_oneshot` never checks oligo length at
all.** The 40–350 nt rule lives only in `ordering.check_eligibility`, reached through
`workflow.plan_order`. That is harmless today because `MAX_FRAGMENT_AA = 90` caps oligos at
about 306 nt, but it is an unguarded path, and anything that lengthens the DNA without
lengthening the protein would walk straight through it.

---

## 5. The two assembly routes

### Route A — intron inside the exon oligo

The intron becomes part of a fragment's payload. The wet lab orders the same number of
pieces and the assembly protocol is untouched.

**What it costs.** Only the 70 bp truncated intron fits, only on 9S and 14S, and it buys
3.8-fold rather than 5.5-fold. **19S cannot take an intron at any length.** A fixed ~70 nt
block would also appear in every library member, and `homology.diversify_library` cannot
break it, because `avoid_sequences` only constrains the mutable coding span.

There is a further unknown: the published deletion series does not give the retained sequence
of the 70 bp variant, so **its GC and homopolymer profile are unknown** and it may not pass
§4.1 at all. That would have to be obtained before Route A could be costed properly.

### Route B — intron as its own part

The intron is a separate module with four-base overhangs matching its insertion site, and the
exon fragments are cut where it goes. This is what MoInClo does: a level −1 intron module,
BbsI-flanked for the exon junctions, with outer BsaI sites for cloning.

**This sidesteps the 350 nt cap entirely.** Full 145 bp rbcS2i1 becomes available on every
architecture, and so does the stronger 253 bp LHCBM1i2.

**It fits the kit's existing architecture better than I expected.** GRASP assembles level 0
with **BpiI** and level 1 with **BsaI** — and PPR-opt already records
that BpiI is an isoschizomer of BbsI, same `GAAGAC` site (`policy.py:59`). So the two-enzyme
scheme MoInClo needs is the scheme GRASP already uses, and PPR-opt already models BbsI as a
wrapping enzyme (`substrates.py:134`). GRASP also places its own junctions in the ARELF
region at residues 20–24, chosen for codon diversity — the same reasoning PPR-opt applies.

**What it costs.** Three things, all real:

1. **Fidelity.** PPR-opt currently chooses every junction overhang freely to minimise
   mis-ligation. At an intron junction the overhang is fixed by the intron's ends. Pinning it
   will lower predicted fidelity, and by how much is precisely the kind of thing PPR-opt is
   built to measure — it should be quantified, not assumed small.
2. **Reachability.** The cut must land on a G|G site. `balanced_cuts` allows only ±4 residues
   (±12 nt) of movement (`arelf.py:274`), and G|G sites occur every ~13.9 nt, so a site is
   usually but not always within reach. The better fix is to constrain the optimiser to
   *create* G|G at the chosen cut rather than hunt for one.
3. **Enzyme policy.** `igem_rfc1000` forbids BbsI sites while this route uses BbsI for
   assembly. PPR-opt already supports BbsI as an assembly enzyme, so this is a policy question
   — which enzyme is "ours" at which level — not a blocker.

### Recommendation

**Route B.** Route A cannot serve 19S at all and depends on a 70 bp sequence whose compliance
is unverified. Route B matches a published working system and reuses an enzyme scheme GRASP
already uses.

*Corrected 2026-10-16:* an earlier version of this paragraph also claimed Route B was "the
only route that can carry the stronger intron". There is no route that carries it. LHCBM1i2
breaches the GC band 0.680 at offset 14, inside the intron rather than at a seam, so no
wrapping enzyme or pad moves it. rbcS2i1 is not a preference; it is the only one that fits.
See `docs/intron_fidelity.md`.

---

## 6. What it would cost to build

From a full survey of the codebase at `30cf908`. **Nothing in this repository models
non-coding DNA.** The CDS is built as `reverse_translate(protein)` (`codons.py:326`) and is
therefore exactly `3 × len(protein)` by construction.

**Four hard refusals.** Three are codon-multiple gates: `assembly.py:106`, which rejects any
CDS that is not a whole number of codons — and rbcS2i1 at 145 nt is not — plus
`products.py:334` and `products.py:379`. The fourth, `codons.py:369`, is a different thing and
worth naming separately: a locked-site bounds check that raises when a lock offset computed
against a longer sequence falls outside the optimiser's own. There is no whole-number-of-codons
refusal anywhere in `codons.py`.

**About fifteen silent-corruption sites** converting residue index to nucleotide offset by
multiplying or dividing by 3. These do not raise; they address the wrong bases. The ones that
matter most: the overhang locks at `design.py:422`, `design.py:476` and `search.py:172`; the
fragment boundary arithmetic at `assembly.py:109–137`; `arelf.py:46`; and **every PPR repeat
annotation in the GenBank export** at `export.py:79–95`, which would place every feature
after the intron at the wrong coordinates.

**Nine tests and six validation scripts** pin the invariant deliberately, including
`validation/integration.py:46` (`len(cds) != 3 * len(protein)`) and
`tests/test_assembly.py:131` (`test_non_codon_sequence_rejected`). These are guards, not
incidental assertions, and would need rewriting rather than deleting.

**One genuine seam exists.** `extra_constraints` (`codons.py:377`) is the caller-supplied
specification hook already used by `wrapper_gc.py`, and dnachisel 3.2.16 supports
`EnforceTranslation(location=)` and `AvoidChanges(location=)`. Per-exon translation with a
frozen intron is achievable within the library we already depend on.

**Two things nothing would catch.** `n_fragments` derives from protein length
(`design.py:398`), so an intron does not increase the fragment count — its whole length lands
in one oligo and nothing notices. And an intron is invisible to off-target screening:
`design.py:561` passes only the 9/14/19-nt target RNA to `offtarget.scan`, while a
host-derived intron is by definition host-identical sequence.

**It would invalidate the current evidence.** The pinned 200-design corpus and the release
bundle are both current. Changing how the CDS is constructed moves every number in them.

---

## 7. Risks that would need their own answers

**Repeated copies.** The spacing rule implies three or four identical intron copies in one
gene. Plant transformation guidance recommends avoiding identical stretches of roughly 30–50
bp between related sequences precisely because repeats drive recombination, deletion and
rearrangement. Against that, Baier 2018 reports using up to 15 copies in large fusions with
no observed alternative splicing or recombination. The tension is unresolved in the
literature; it is not ours to resolve, but it is ours to declare.

**Cryptic splicing.** Introducing strong splice signals into a synthetic CDS risks activating
cryptic sites elsewhere in it, and the consequence is context-dependent. The 2025 random-intron
work found cryptic sites absent beyond about 150 nt, implying a practical minimum intron
length near 70 nt — which happens to be the truncated variant. A real implementation should
screen the assembled transcript for cryptic donors and acceptors, which is a capability
PPR-opt does not currently have.

**Silencing is upstream of all of this.** Introns do not address the epigenetic pathway that
selectively silences transgenic DNA. Strain choice (UVM4/UVM11) and promoter choice sit
outside PPR-opt's scope and may matter more than anything PPR-opt can do to the sequence.

---

## 8. Verdict

Staged, because the three stages have very different risk.

**1 — Now, zero risk, no code.** Record two things in `WHAT_PPROPT_DOES.md`: that codon
optimisation alone is documented as insufficient for this host, and that rbcS2i1 passes our
synthesis policy unchanged. The first is a limitation of PPR-opt that is true whether or not
introns are ever built, and it belongs in the honest-limitations list next to "nothing
biological is measured". We should not wait for an implementation to admit it.

**2 — Next, low risk.** An advisory output: report the G|G sites and the suggested intron
positions alongside each design, with the exon-length rule applied, and change nothing about
the CDS. No published number moves, no invariant breaks, and the wet lab gets something
actionable. This is roughly the Wild Worm Codon Adapter's behaviour minus the insertion.

**3 — Later, expensive.** Route B proper, with a host gate that refuses `e_coli` and both
chloroplast hosts, oligo-length checking promoted into `design_oneshot`, and the
residue-to-nucleotide arithmetic routed through a single mapping function. This should not
start until the current release and corpus are no longer load-bearing.

---

## 9. What this does not decide

- **No expression was measured here, and adding this would measure none.** Every fold-change
  quoted is from another group's constructs. PPR-opt would be predicting that a published
  effect transfers to our gene, which is a weaker claim than the numbers look.
- **Whether LHCBM1i2 passes the synthesis policy.** It was not tested — only rbcS2i1 was. It
  is the stronger intron and the one Route B would most want.
- **Whether the 70 bp truncated rbcS2i1 passes the policy.** The published deletion series
  does not give the retained sequence.
- **Whether a designed-from-scratch intron would splice in this alga.** The 2025 random-intron
  result says sequence barely matters in human cells; Baier 2020 says foreign introns mostly
  fail in *Chlamydomonas* until their splice sites are adapted. Those are not reconcilable
  from the literature alone, and the question is worth an experiment.
- **Whether the wet lab can run a two-enzyme MoClo assembly**, and which enzyme is ours at
  which level.
- **Whether any of this matters more than strain and promoter choice**, which are outside
  PPR-opt entirely.

---

## 10. Sources

**Primary, *Chlamydomonas***

- Baier T, Wichmann J, Kruse O, Lauersen KJ (2018). Intron-containing algal transgenes mediate
  efficient recombinant gene expression in *Chlamydomonas reinhardtii*. *Nucleic Acids
  Research* 46(13):6909–6919.
  <https://academic.oup.com/nar/article/46/13/6909/5042046>
- Baier T, Jacobebbinghaus N, Einhaus A, Lauersen KJ, Kruse O (2020). Introns mediate
  post-transcriptional enhancement of nuclear gene expression in *Chlamydomonas reinhardtii*.
  *PLoS Genetics* 16(8):e1008944.
  <https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008944>
- Lumbreras V, Stevens DR, Purton S (1998). Efficient foreign gene expression in
  *Chlamydomonas reinhardtii* mediated by an endogenous intron. *The Plant Journal*.
  <https://onlinelibrary.wiley.com/doi/10.1046/j.1365-313X.1998.00145.x>
- MoInClo: a MoClo-mediated intron insertion system for enhanced transgene expression in
  *Chlamydomonas reinhardtii*.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11925875/>
- Neupert J, Karcher D, Bock R (2009). Generation of *Chlamydomonas* strains that efficiently
  express nuclear transgenes. *The Plant Journal*.
  <https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-313X.2008.03746.x>
- An epigenetic gene silencing pathway selectively acting on transgenic DNA in
  *Chlamydomonas*. *Nature Communications* 2020.
  <https://www.nature.com/articles/s41467-020-19983-4>

**Mechanism and general IME**

- Sequence-dependent and -independent effects of intron-mediated enhancement learned from
  thousands of random introns. *Nucleic Acids Research* 2025, 53(4):gkaf097.
  <https://academic.oup.com/nar/article/53/4/gkaf097/8039171>
- Shaul O (2017). How introns enhance gene expression.
  <https://www.sciencedirect.com/science/article/abs/pii/S1357272517301541>
- Intron-mediated enhancement review, *Frontiers in Plant Science* 2011.
  <https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2011.00098/full>

**PPR and the kit**

- GRASP: a modular toolkit for building synthetic pentatricopeptide repeat RNA-binding
  proteins. *Nucleic Acids Research* 2025, 53(20):gkaf1169.
  <https://academic.oup.com/nar/article/53/20/gkaf1169/8321212>
- Use of plant chloroplast RNA-binding proteins as orthogonal activators of chloroplast
  transgenes in *Chlamydomonas reinhardtii*. *Algal Research* 2021.
  <https://www.sciencedirect.com/science/article/abs/pii/S2211926421003544>
- Crozet P et al. (2018). Birth of a photosynthetic chassis: a MoClo toolkit for the green
  microalga *Chlamydomonas reinhardtii*. *ACS Synthetic Biology* 7:2074–2086.

**Tooling precedent**

- The Wild Worm Codon Adapter: automated codon adaptation and intron insertion for transgenes
  in non-*Caenorhabditis* nematodes.
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8496300/>

**Sequence**

- NCBI X04472.1 — *C. reinhardtii rbcS2* gene. Intron 1 at 1307..1451 (145 bp); intron 2 at
  1532..1860 (329 bp).
