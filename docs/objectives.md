# Objectives and constraints

Every quantity an optimiser may move, defined precisely enough that an independent
implementation computes the same number. Where a definition has a free parameter, the
parameter is named and enters the cache key.

**Hard constraints are not objectives.** A candidate violating one is infeasible and is never
rescued by a good score. This separation is stated here once and holds everywhere.

---

## 0. Scope: evaluate the contract on what leaves the building

Every constraint below is evaluated on the **artefact that is actually delivered** — the
sequence a vendor synthesises, the fragment a reaction consumes — never on a sub-span that is
convenient to check.

This rule is stated first because violating it produced the same defect four times in this
project, each time at a different layer and each time with the constraint technically
satisfied on the thing that was measured:

| measured | delivered | result |
|---|---|---|
| the optimiser's coding span | the whole module | GC 0.640 → **0.660** |
| the BsaI matrix | a BbsI level-0 reaction | a fictitious fidelity bottleneck |
| the bare insert | the released fragment | GC 0.640 → **0.660** |
| the bare insert, in the collection optimiser | the ordered substrate | 16 of 18 runs returned an inventory the export step refused |

A constraint checked one layer inside the deliverable is not a weaker check. It is a check of
something else.

**Fixing it four times in four places was the mistake.** The rule now has one implementation:
`substrates.substrate_problems(insert, block)` builds the sequence that would be ordered and
judges *that*. Recoding, the collection optimiser and the joint search all call it, and a
fifth acceptance path that does not is a defect regardless of what it checks instead.

## 1. Hard constraints

A candidate is feasible only if **all** of these hold. Each is checked on the final sequence
that would actually be ordered, not on an intermediate.

| constraint | statement |
|---|---|
| protein identity | the translation over the declared coding interval equals the intended protein, exactly |
| reading context | the coding interval is a whole number of codons and contains no internal stop |
| frozen interfaces | both four-base boundaries equal the version they must interoperate with |
| forbidden sites | no active recognition site from the enzyme profile, **on either strand**, outside its intended role |
| sequence limits | GC within band in every window; homopolymer run within limit |
| compatibility | junction roles permit every adjacency; no duplicated or self-complementary overhang where the assembly model forbids it |

An *intended* recognition site — the one wrapping an oligo so the enzyme can release it — is
handled by its role. An accidental internal occurrence is a violation. The two are
distinguished by position and role, never by hoping there is only one.

---

## 2. Codon adaptation (CAI)

The geometric mean of relative synonymous weights over the declared coding interval.

For codon *c* encoding amino acid *a*, with frequency *f(c)* in the supplied table:

    w(c) = f(c) / max{ f(c') : c' encodes a }

and over the codons actually used:

    CAI = exp( (1/n) * sum log w(c_i) )

**Declared treatment of the awkward cases**, because these are where independent
implementations diverge:

- **Methionine and tryptophan are excluded.** They have one codon each, so `w = 1` by
  construction; including them dilutes the score with values no optimiser can change.
- **Stop codons are excluded** from the coding-interval score.
- **A codon absent from the supplied table** is excluded from the product and counted in
  `missing_weights`, which is reported. It is not silently treated as weight 1.
- **If any included weight is zero**, CAI is 0 rather than `-inf`; the zero-weight codons are
  reported. A single zero must not make two otherwise different sequences incomparable.
- **The contributing interval and codon count are emitted with the score.** Two CAI values over
  different intervals are not comparable, and reporting the interval is what makes that
  visible rather than assumed.

Never average CAI across incomparable regions. A per-module mean and a whole-product value are
different statistics.

---

## 3. Within-sequence repetition

Parameterised by *k*, which is reported with every value and enters the cache key.

- `duplicated_kmers(s, k)` — the number of *k*-mer occurrences beyond the first for each
  distinct word: `len(windows) - len(set(windows))`.
- `distinct_repeated_words(s, k)` — how many distinct words occur more than once.

These are different numbers and are never used interchangeably. Orientation: forward strand
only unless `include_reverse_complement` is set, which is reported when it is.

Default `k = 20`, for continuity with the existing QC. Configurable.

**These do not make a usable search objective for this kit, and that was measured rather than
assumed.** Duplicated k-mers are identically zero over all 42 deposited and recoded modules at
k = 20, 16, 12 **and** 10, registering only at k = 8. A search declaring three objectives was
ranking by two. The third axis is now `synthesis_fitness` (§4a); the repetition figures are
still computed and reported as diagnostics, because a repetitive candidate is a real problem
even when this kit never produces one.

---

## 4. Between-module sharing

Measured over **distinct physical module records only**. Reusing one selected record across
several target products is the point of a reusable inventory; counting that reuse as sharing
would penalise the feature for working.

Two separate quantities, reported separately:

**Collection penalty.** With `c_i(w)` the number of occurrences of *k*-mer `w` in record *i*:

    P = sum over pairs i < j  of  sum over w  of  min( c_i(w), c_j(w) )

This is an engineering objective for steering a search. **It is not a recombination
probability** and carries no biological interpretation.

**Worst shared tract.** The length of the longest exact substring common to any two distinct
records. Reported independently.

**Lowering `P` does not necessarily lower the worst tract**, and no result may claim it does.
`P` aggregates many short coincidences; the worst tract is a single extreme. A search that
improves one can leave the other unchanged, and that has to be measured rather than inferred.

---

## 4a. Composite synthesis fitness

The third search objective, maximised, in [0, 1], computed on the **ordered substrate**.

Four terms, each a *distance* rather than a pass/fail so that ranking still moves when every
candidate comfortably satisfies the constraints:

| term | 1.0 means | why it is here |
|---|---|---|
| GC centrality | every window at the declared band's midpoint | spans 0.417-0.574 across our substrates |
| homopolymer headroom | no run longer than 1 | spans 3-4 today |
| internal repetition | every k-mer unique | degenerate now; fires on a genuinely bad candidate |
| collection sharing | no k-mer shared with the rest of the inventory | varies with the collection |

**GC centrality is measured against the declared band's centre, not a fixed 0.5.** A policy
permitting 0.15-0.85 is not asking for the same sequence as one permitting 0.35-0.65, and
scoring both toward 0.5 would impose a preference no profile declared. The **worst** window
sets the term, not the mean: one bad window is what a vendor's model reacts to.

**The weights are a declared engineering choice, not a measurement**, and are reported with
every score together with the components, because an aggregate whose parts cannot be inspected
is a number nobody can check.

**Degeneracy is reported.** `inventory_fitness` returns the spread and a `degenerate` flag. An
axis with no spread is not an axis, and the failure that produced this section must announce
itself rather than be discovered a second time.

---

## 5. Reaction fidelity

Predicted ligation fidelity of an overhang set, from the Pryor et al. 2020 mismatch data.

The pool for a directional ratio is the set **and its Watson–Crick pairs** — the denominator
sums over both strands. This is the definition the primary source uses; an earlier version
here summed only the reverse complements and produced a set/rc(set) asymmetry of 0.43, where
the corrected form gives 2.2×10⁻¹⁶.

Aggregation convention is versioned and recorded. The shipped convention is the geometric
mean of directional products.

**A reaction is scored only when its participant list is complete**, and every scored reaction
carries the overhang set it was scored on and where those ends came from. Fidelity is a
property of the whole competing set in one tube, so a score over a *subset* of a reaction is
not that reaction's fidelity and is never reported in the same field. Where a subset is worth
reporting it is carried separately and labelled a diagnostic.

This is not a formality. The level-1 block subset of the 19S product scores 0.996046; adding
only the co-assembled ends already known to be in that reaction gives 0.742067. A subset score
can sit on either side of the figure it would be mistaken for.

**A multi-reaction compilation reports every reaction's fidelity.** The default search summary
is the **minimum** across reactions, labelled as such. It is never an unlabelled product of
the reactions, because a product reads as a predicted overall yield and is not one.

Fidelity is a **prediction from published mismatch data**, not a measured yield.

---

## 6. Dominance and the observed Pareto set

Objective directions are declared before the search runs.

**A dominates B** iff A is no worse than B in every objective and strictly better in at least
one. Numerical comparison uses an absolute tolerance `1e-9` per objective; differences within
tolerance count as equal, so a floating-point wobble cannot manufacture dominance.

The result is the **observed** nondominated set — nondominated among the candidates actually
evaluated. It is called the globally optimal front only if the feasible space was exhaustively
enumerated, which for the real inventory it is not.

**Degeneracy is reported, not disguised.** If every observed candidate has the same fidelity,
the report says so rather than presenting a trade-off that does not exist.

**Shortlists are rescored under matched conditions.** A refined selected candidate must never
be compared against stale, cheaply evaluated alternatives without that difference being
marked.

---

## 7. Default recommendation policy

Stated in advance so the recommendation is a rule, not a preference applied after seeing
results:

1. Feasible candidates only.
2. Retain those within the documented fidelity tolerance of the best observed feasible
   fidelity.
3. Among those, select by the declared sequence-quality policy.
4. Break remaining ties deterministically.

The policy is displayed with the recommendation, and any other valid observed Pareto point can
be selected instead.

---

## 8. What none of these measure

No objective here predicts expression, binding affinity, specificity, recombination frequency,
assembly yield or manufacturability. They are sequence measurements and predictions from
published data. An improvement in any of them is an improvement in that measurement and
nothing more.
