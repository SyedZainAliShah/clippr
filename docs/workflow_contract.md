# Workflow contract

What each capability takes, returns, and refuses. Written before the code it governs, so the
implementation answers to this document rather than the document describing whatever got
built. Behaviour stated here is CLIPPR's own; where it follows a published fact rather than an
engineering choice, `methodology_sources.md` names the source.

Semantic contracts, not mandated class names. An existing CLIPPR type that already satisfies
one of these is the right implementation of it.

---

## 0. Coordinate conventions

Stated once, because most of the silent defects in this area are off-by-one.

- **Internal intervals are zero-based and half-open.** `(start, end)` covers `end - start`
  bases beginning at `start`. Every internal coordinate — cut positions, module spans,
  feature spans, coding intervals — uses this form.
- **GenBank export converts on the way out.** A written GenBank record is one-based inclusive;
  the conversion happens in the exporter and nowhere else.
- **An overhang is named by its coding-strand sequence, 5'→3'.** The four bases a Type IIS
  enzyme leaves on the top strand. The strand that physically anneals to it is the reverse
  complement, and code that needs the annealing partner computes it explicitly rather than
  assuming the stored value already is one. This distinction is the reason a destination pair
  stored as `("CTCA", "CGAG")` meets as `CTCA`/`CTCG` in the tube.
- **Shared junction bases belong to both neighbours and are counted once in the product.**
  Module spans in an assembled product therefore overlap by four bases. That overlap is the
  assembly, not an error. A product's length is `sum(len(module)) - 4 * (n_modules - 1)`.
- **A module's reading frame is the offset of its first whole codon**, in `{0, 1, 2}`,
  measured from the module's own first base. It is derived from the module's position in the
  product, never searched for by trying frames and picking the one without stop codons.

---

## 1. Design context

Everything that changes a result. Two runs with the same context and seed must produce the
same sequences; two runs with different contexts must never be compared without saying so.

| field | meaning |
|---|---|
| `host_label` | free-text organism name, for reporting |
| `expression_compartment` | where the **protein** is expressed (e.g. nuclear) |
| `target_compartment` | where the **RNA target** resides (e.g. chloroplast) |
| `genetic_code` | NCBI translation table id |
| `codon_table` | normalised frequencies, plus `effective_sha256` of the table as the optimiser received it |
| `enzyme_profile` | active forbidden recognition sites |
| `destination` | physical external ends, as coding sites, 5' then 3' |
| `ligation_matrix` | matrix name and fingerprint |
| `aggregation_version` | which fidelity aggregation convention applies |
| `synthesis_profile` | GC band, window, homopolymer limit, repeat *k* |
| `product_profile` | vendor product rules, when an order is in scope |
| `seed`, `budget` | search determinism and limits |

**`expression_compartment` and `target_compartment` are separate fields.** A PPR expressed in
the nucleus to bind a chloroplast transcript has different values for each, and collapsing
them would silently pick the wrong codon table or the wrong off-target genome.

**The current *Chlamydomonas* nuclear context is the benchmark, not a default for everything.**
Its GC band (0.35–0.65 over 50 nt) was chosen for that context. Applying it to another host is
a decision that must be made explicitly, not inherited.

## 2. Inventory

A versioned collection of modules that can be compiled into targets.

| field | meaning |
|---|---|
| `module_id` | stable identity, independent of sequence |
| `version` | which realisation of that identity this is |
| `dna` | the module's full insert sequence |
| `source_record` | supplying file identity and hash |
| `coding_interval` | half-open span of whole codons this module carries, in its own frame |
| `frame` | 0, 1 or 2 |
| `five_interface`, `three_interface` | the four boundary bases |
| `junction_roles` | which neighbours this module may sit beside |
| `context` | the design context it was built under |
| `constraints`, `objectives` | evaluated values, with the definitions that produced them |

**Deposited and redesigned records stay distinct.** A recoded module is a new version of the
same identity, never an overwrite of the deposited one.

**An interface assignment belongs to an inventory version.** Modules from inventories with
different interface assignments must not be mixed; a compile that is handed such a mixture
fails with both version identities named, rather than producing a product that cannot ligate.

**An unchanged module is a valid module.** If recoding finds no improvement, the incumbent is
the deliverable, recorded as unchanged — not a gap in coverage.

## 3. Compilation

Turning a target into an ordered set of modules and the product they yield.

| field | meaning |
|---|---|
| `target` | the original RNA target and a stable id |
| `inventory_version` | which inventory this was compiled from |
| `modules` | ordered `(module_id, version)` pairs |
| `stage_products` | the sub-assembly products, in reaction order |
| `junctions` | every join: left, right, overhang, coordinate |
| `protein` | the intended protein |
| `boundary_scope` | what the product does **not** include |
| `order_sequences` | wrapped sequences, when a synthesis route is in scope |

**`boundary_scope` is mandatory and is stated on every export.** A joined set of inserts is
not an expression construct: it carries no acceptor backbone and no vector context. A FASTA
travels without its JSON, so the scope is repeated in every record header.

**Failures name identities.** A missing module or an incompatible version reports exactly
which module and which version, never a bare "compile failed".

## 4. Search result

What any optimiser returns, whether it improved anything or not.

| field | meaning |
|---|---|
| `incumbent` | what was already there |
| `candidates` | every candidate evaluated, feasible or not |
| `feasible` | per candidate, with a **reason** when false |
| `objectives` | the objective vector per candidate |
| `nondominated` | the observed Pareto set |
| `recommended` | the single default, and the rule that selected it |
| `evaluations` | candidates fully evaluated — **not** cache hits |
| `timings`, `completion` | elapsed, and `complete` / `budget_exhausted` |

**A cache hit is not an evaluation.** They are counted separately, because a search reporting
"200 evaluations" that actually computed four is not reporting search effort.

**`budget_exhausted` never means "no improvement exists".** It means the search stopped. The
distinction is the difference between a measured negative result and an unexamined space.

**Returning the incumbent is a valid outcome.** An optimiser that finds no useful trade-off
and says so has worked correctly.

## 5. Order plan

| field | meaning |
|---|---|
| `mapping` | source design and module for every ordered sequence |
| `sequences` | the actual wrapped sequences, not the bare CDS |
| `lengths`, `quantities`, `scales` | as ordered |
| `pools` | pool membership |
| `eligibility` | per product profile: violations, warnings, unresolved rules |
| `price_version` | profile version, source and date |
| `estimate` | cost, or an explicit unavailable |

**Identical sequences from different orders are not silently merged.** Any deduplication
preserves every source mapping and the intended quantities; two designs that happen to share
a fragment still need both.

**No filler.** If a request cannot meet a profile's minimum count, the plan is "no eligible
plan" naming the constraint — never invented oligos to reach a threshold.

**A local eligibility pass is not vendor approval**, and a historical constant is never
presented as a current quotation.

**One feasibility test, shared by every step that accepts a candidate.**
`substrates.substrate_problems(insert, block)` builds the substrate that would be ordered and
judges *that*; recoding, collection acceptance, junction feasibility and export all call it.
A step that accepts a module the order plan would refuse is a contract violation, not a
difference of opinion between two checkers — the collection optimiser judging bare inserts
while export judged substrates is what made 16 of 18 search runs unshippable.

---

## 6. The six user-visible tasks

1. **Design CDSs for synthesis** from explicit RNA targets. *(exists)*
2. **Load an inventory and compile targets** — deposited or saved redesigned.
3. **Recode an inventory** for a supplied host context, interfaces fixed.
4. **Optimise an inventory as a collection**, optionally producing alternative versions.
5. **Explore junction/codon trade-offs** and select a compatible inventory.
6. **Check order eligibility, plan pools, export the package.**

Each is callable on its own. Selecting an alternative in task 5 changes the inventory and
interface version, and therefore requires recompilation — the interface says so at the point
of selection rather than leaving the user to discover it.

## 7. Evidence contract

Every validation artefact records: the command, input hashes, source and worker hashes,
dependencies, settings, requested identities, completed identities, failures, elapsed time and
output hashes.

Every finding states its claim, its falsifier and the command that settles it, and is labelled:

- **RAN** — executed here, result read from that run
- **READ** — inspected an artefact a previous run produced
- **ESTIMATE** — not measured
- **PROPOSED** — specified but not yet run

Expected values in tests come from declared rules and independent calculation, never from the
production serialiser or reconstructor. A checker that reads the value under test and agrees
with it has measured nothing.
