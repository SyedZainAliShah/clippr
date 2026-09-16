# Methodology sources

Where each fact and each method comes from, and — just as important — which parts are
CLIPPR's engineering choices rather than anything a source dictates.

Labels: **primary fact** (published or deposited, we did not choose it) · **documented
method** (standard mathematics or algorithm from the literature) · **engineering choice**
(ours; a different project could justifiably choose otherwise) · **vendor rule** (external,
versioned, expires).

---

## 1. Biological architecture — primary fact

**The PPR code.** Each ~31-residue repeat reads one RNA base through two specificity residues
at positions 5 and 35 (the 5th and last of the repeat as indexed here):

    TN → A      NN → C      TD → G      ND → U

Transcribed from the published code, not lifted from any implementation. The protein is a
deterministic function of the target, which is why CLIPPR generates designs rather than
searching a catalogue.

**Scaffold constants.** `N_TERMINAL` and `REPEAT_TEMPLATE` in `clippr.biology` derive from the
paper and were checked against the authors' deposited construct: for target `UUACACGUG`,
CLIPPR's 302-aa protein from index 2 equals the deposited sequence's first 300 residues.

**Deposited records.** Three GenBank files from the authors' repository at commit
`b569a8c499c87f0c825e4e5db237601dde99a1d5`, preserved under
`work/provenance/grasp_primary`. Independently re-fetched from the pinned URLs and confirmed
byte-identical. The authors' repository declares **no licence**: cite and read, never vendor.

**Module inserts.** Supplementary Table S1, user-supplied, sha256
`6324ede3e3863ced60f12cc822ebe68f530e9ebc6a5ebc45c490a8816930ef35`. Deliberately not shipped
— see `NOTICE.md`. Every code path that needs it degrades to an explicit "unavailable" rather
than a crash or a silent skip.

## 2. Ligation fidelity — primary fact

**Pryor et al. 2020**, ligase fidelity matrices, CC BY 4.0, shipped under
`src/clippr/data/matrices/`. `BsaI-HFv2` fingerprint `f5e4da3cacff1577`.

The pooling rule — a directional ratio's denominator covers the set **and its Watson–Crick
pairs** — is the source's definition. Getting this wrong is not a rounding difference: it
produced a 0.43 asymmetry between a set and its reverse complement where the correct form
gives 2.2×10⁻¹⁶.

The stage-matched conditions (BsaI-HFv2 at levels −1 and 1, BbsI-HF at level 0, 37↔16 °C) come
from the same source.

**Aggregation convention is an engineering choice.** Two defensible conventions exist —
geometric mean of directional products (shipped) versus per-junction pooled ratio. They
produce genuine rank reversals on constructed sets. On the real candidate pools the choice
changed selection in 0 of 200 cases, which is a measurement, not a proof that it never
matters. The convention is versioned and recorded with every result.

## 3. Assembly mechanics — primary fact

Type IIS enzymes, recognition and cut geometry:

| enzyme | site | cut | overhang |
|---|---|---|---|
| BsaI | `GGTCTC` | (1/5) | 4 nt, 5' |
| BbsI (= BpiI) | `GAAGAC` | (2/6) | 4 nt, 5' |
| SapI | `GCTCTTC` | (1/4) | 3 nt, 5' |

MoClo fusion sites `AATG` (carries its own initiator) and `AGGT` (N-terminal fusion) are
standard-defined.

Destination overhangs as coding sites: level −1 `("ACAT","TTGT")`, level 0 `("CTCA","CGAG")`,
level 1 `("GGAG","CGCT")` — read from the deposited vectors.

**Two different reactions are called "level 1".** The pair above is the acceptor vector's own,
used when the one-shot path assembles a fresh construct into it. A GRASP *block* assembly's
level-1 reaction joins assembled blocks and takes its ends from the compiled product
(`AATG`/`AGGT` … `TTCG`), derived in `validation/experiments/level1_geometry.py` from the 28
deposited block plasmids and confirmed against Table S1. The two describe different tubes and
neither figure contradicts the other.

**Independent validators restate these rather than importing them.** V5 writes out the enzyme
geometry and the destination pairs explicitly: taking them from `Bio.Restriction`, which the
production code also uses, would make a shared dependency look like independent corroboration.

## 4. Sequence design — documented method

**DNA Chisel** (Zulkower & Rosser 2020) provides the constraint-satisfaction engine —
`EnforceTranslation`, `EnforceSequence`, `UniquifyAllKmers`. Retained deliberately: it is a
published, maintained tool for exactly this problem, and reimplementing it would add risk
without adding capability.

**CAI** is the standard Sharp & Li formulation. The treatment of M/W, stops, missing weights
and zero frequencies is specified in `objectives.md` — the underlying statistic is standard;
those decisions are ours and are stated because independent implementations differ on them.

## 5. Multi-objective search — documented method

**Pareto dominance** is the textbook definition, restated in `objectives.md` with an explicit
numerical tolerance. Nothing here is novel; it is written down so an independent
implementation agrees.

**Simulated annealing** is Kirkpatrick et al. 1983: accept improving moves; accept a worsening
feasible move with probability `exp(-Δ/T)` under a declared cooling schedule.

**Every annealing parameter is an engineering choice** — schedule, temperatures, move mix,
iteration counts, scalarisation weights. They are chosen from CLIPPR's own measurements and
are not copied from any reference implementation's settings. Reported with every result so a
reader can see what was actually run.

The move generator and the acceptance test draw from a run-local RNG seeded per run. No
process-global random state, because global state makes a "reproducible" result depend on what
else ran first.

## 6. Vendor rules — vendor rule, versioned and expiring

Product profiles carry name, region, currency, **source URL and retrieval date**, length
bounds, count and scale limits, price tiers, modification rules, quantity assumptions, and the
rules that cannot be machine-checked.

Constraints:

- Rules and prices come from **current official vendor documentation**. Hardcoded numbers
  found in any other implementation are not authoritative and are not imported as fact.
- Where prices are unavailable, an explicitly supplied **dated** price table is accepted;
  otherwise cost is reported **unavailable** while eligibility and export still function.
- Monetary arithmetic is exact decimal. Tax, shipping and modification treatment is explicit.
- The existing `109.00 EUR` pool and `1.63 EUR/oligo` phosphorylation figures are **historical
  constants measured across the 200-design corpus** — a list price, not a quote, not fetched
  from a vendor. They are labelled historical and are never presented as current.
- A local eligibility pass is an implementation check. **It is not vendor approval** and no
  output may imply that it is.

## 7. Relationship to the reference implementation

GRASP Designer (`grasp-library-designer`, v0.1.15, commit `8882759ac267e79f`) is used as a
**behavioural reference and benchmark**, not as implementation material.

- Its reference imports and execution live in an isolated benchmark harness. Production and
  its normal tests run without the reference installed or its repository present.
- Benchmarking goes through public entry points and exported artefacts.
- Oracle output can be comparison evidence. It is never the sole expected answer in a
  correctness test — an expected value has to come from a declared rule.
- Reference functions, control flow, classes, notebooks, comments, internal structures,
  parameter bundles and test fixtures are not copied or translated. Renaming copied code
  would not be independence.
- **Reference implementation files may be read to understand the method** (authorised
  2026-09-16). What is taken is understanding, never code: anything adopted is reimplemented
  from the documented method and primary data, and any vendor rule learned this way is
  re-sourced from the vendor's own published page before it is relied on. Findings from
  reading are recorded with their source in `docs/reference_implementation_findings.md`.
  Until that date this line read "reference implementation files are not read to guide
  coding"; the change of posture is recorded rather than quietly applied.

**Stated plainly:** this project has had prior exposure to that repository. The process
described here is independent implementation using documented methods and primary data, with
the reference as a behavioural benchmark. It is **not** strict clean-room development, this
document does not describe it as such, and the engineering process alone does not resolve
every licensing question.

The earlier v0.1.5 wheel used in this project's oracle comparisons is a **different
reference** and does not contain the reusable-library workflow. Comparisons against it
establish nothing about reusable-library parity.

**The matched benchmark reads artefacts, not code.**
`validation/experiments/matched_oracle_benchmark.py` compares the ordered sequence each system
emits, reading the reference's preserved `optimized_library.csv` from disk. It neither imports
nor executes the reference, and its constraint checker is written inside the benchmark rather
than imported from either side — a checker borrowed from one of the two systems would grade
its owner on home rules. The one place the two rule sets are reconciled is the intended-site
exemption, which keys on the recognition sequence precisely because BbsI and BpiI are
isoschizomers of `GAAGAC` (§3).

## 8. Provenance of data versus provenance of implementation

Tracked separately, because they fail differently.

Publicly *accessible* sequences, figures and tables are not automatically redistributable.
Supplied-input routes stay in place where the licence is unclear. Any newly distributed input
carries a documented basis for distributing it.

## 9. What no source here supports

Nothing in this project has been validated at the bench. No binding was measured, no assembly
performed, no expression observed. Every fidelity figure is a prediction from the Pryor 2020
data; every sequence figure is a sequence measurement; every cost figure is an estimate from a
dated profile.
