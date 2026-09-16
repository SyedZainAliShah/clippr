# Handoff for independent verification

Everything below is reproducible from this checkout. Each claim names the command that
produces it and the artefact it writes, so nothing here needs to be taken on the word of the
implementation's own report generator.

Labels: **RAN** — executed here, result read from the run. **READ** — inspected an artefact a
previous run wrote. **ESTIMATE** — not measured.

---

## 1. Environment identity

| | |
|---|---|
| package | `clippr` 0.1.0, MIT |
| source fingerprint | covers 27 files — 25 `.py` plus `parts.json` and `scaffold.json` |
| ligation matrices | Pryor et al. 2020, `BsaI-HFv2` fingerprint `f5e4da3cacff1577` |
| dependencies | biopython 1.88, dnachisel 3.2.16, numpy 2.5.3, pandas 3.0.5, openpyxl 3.1.5, python_codon_tables 0.1.18 |
| reference | `grasp-library-designer` **v0.1.15**, commit `8882759ac267e79f`, its own venv |

The v0.1.5 wheel used by every earlier oracle comparison in this project is a **different
reference** and does not contain the reusable-library workflow. It cannot substitute.

## 2. Test suite — RAN

    python -m pytest -q

**627 collected, 627 passed, 0 failed.** The count is the collected total from
`--collect-only`, not a count of progress dots; an earlier report of "617" came from counting
decimal points in comparison output.

13 tests were added in this pass: 5 for the new FASTA/GenBank product exports, 4 for manifest
batch homogeneity, 2 for the manifest inside the library package, 2 for the widened source
fingerprint.

## 3. Pinned 200-design corpus — RAN

    python validation/independent/regenerate_corpus.py

| | |
|---|---|
| designs | 200 — 100 × 9S, 50 × 14S, 50 × 19S |
| `constraints_ok` | 200/200 |
| QC PASS | 200/200 |
| failures | 0 |
| wall | 19.6 min across all runs |
| planning cache | 5,214 hits / 232 misses |

Artefacts: `work/independent/corpus200_pinned.jsonl`, `corpus200_pinned_summary.json`.

**Its recorded identity is the one it actually ran under.** The journal header carries source
fingerprint `1e0a3c74…`, which is the 25-Python-file definition in force at the time of the
run. The fingerprint has since been widened to 27 files (adding `parts.json` and
`scaffold.json`), and this corpus is **not** relabelled with the new value: the artefact
states the identity it was produced under, which is the whole point of recording one. The
widened definition applies to runs made after the change.

The fingerprint still excludes the experiment worker and the target-generating validation
code, which live outside the package. Those are recorded separately in each validator's
`provenance` block rather than folded into the package digest.

The journal interleaves bookkeeping records with design rows, so the 200-design corpus is a
403-line file. `validation/independent/journal.py` filters the reserved keys; a checker that
reads every line blindly sees the run-identity and checkpoint records as malformed designs,
which is how V5 and V6 first failed against it.

## 4. Independent checks against that corpus — RAN

    python validation/independent/v5_digest_ligate.py
    python validation/independent/v6_constraints.py
    python validation/independent/v7_export_semantics.py

| check | result | independence |
|---|---|---|
| **V5** digest/ligate | 200/200 CDSs reconstructed; **5/5 corrupt controls detected** | reads only `oligo_sequence_5to3`; imports no CLIPPR assembly code |
| **V6** constraints | 200/200 satisfy every hard constraint; 0 disagreements with the designs' own verdicts; 200/200 have zero duplicated 20-mers | expectations written out from the contract and the published PPR code; imports no `clippr.qc` or `clippr.codons` |
| **V7** export semantics | 3 designs, 0 semantic problems, 3/3 malformed controls caught | reparses written files with stock parsers |

V6 rebuilds the **whole expected protein** from the scaffold template rather than decoding
only the specificity residues. The decode-only version accepted a Q→F change at residue 1.

## 5. M2 — joint full-CDS library selection — READ

**Verdict: DO NOT KEEP.** The preregistered gate required improvement on two of three panels;
it improved one. Preregistration frozen before any result existed, sha256
`ca5582455298044a008ce0d0ce939c94b4d795059a4e944a1572c788cd8939d9`.

The mechanism was measured, not assumed: 105 of 105 distinct-prefix pairs share ≥50 nt
(min 59, median 59, max 68). No selector restricted to these banks can beat the 48 nt and
56 nt diversified baselines. **The bank is the binding constraint, not the selector.**

Full detail and the per-panel table: `validation/experiments/m2_result.md`.

## 6. M3 — fixed-interface reusable inventory — RAN

    python validation/experiments/m3_inventory_pilot.py
    python validation/experiments/m3_validate_inventory.py
    python validation/experiments/m3_reference_comparison.py

**Pilot.** 21 modules, one declared host context, **1.7 s of a 600 s cap**.

| | min | median | max |
|---|---:|---:|---:|
| CAI before | 0.1784 | 0.2206 | 0.3033 |
| CAI after | 0.5076 | 0.6450 | 0.7199 |
| gain | +0.3009 | +0.4203 | +0.4956 |

21/21 improved; 17/21 are reused by more than one target and all 17 improved. Module QC
21/21 PASS, unchanged before and after.

**Assembly validation (E4).** 6/6 products assemble, every overlap verified, 0 failures, no
junction-spanning forbidden site, proteins identical to the originals'. The join is written
out in the validation script rather than taken from `products.reconstruct`, so a fault in the
production reconstructor cannot hide a fault in the inventory.

**The negative part, stated plainly.** Longest shared tract between two kit products:
**326 → 317 nt**. Recoding shared modules cannot reduce sharing between products that reuse
the same records. M3 is a reusable-inventory capability, not library diversification.

## 7. E6 — reference comparison — RAN

Reference v0.1.15 run in its own environment under a 3600 s cap; finished in **22.0 s** for
42 modules. It was given our *Chlamydomonas* table in its own `codon,aa,frequency` format.

**Span alignment.** The reference's CDS runs 0–2 codons longer because it includes the codon
spanning a junction, completed by the neighbouring module. Our protein is a contiguous
substring of theirs on **21/21** shared modules, so scores are taken over exactly that common
region, where both encode the same protein on 21/21.

| mean CAI on the common region | |
|---|---:|
| original deposited inventory | 0.2246 |
| ours | **0.6432** |
| reference | **0.6333** |

Per module: reference higher on 9, ours higher on 12, none tied.

**Reading: equivalent.** A 0.0099 mean difference with a 12–9 split is not a result in either
direction. One scorer was applied to both; neither system's self-reported metric was used,
because their `codon_score` and our `cai_after` are not the same statistic.

Declared mismatches, stated rather than aligned: GC band (0.25–0.65 global / 0.15–0.85
windowed vs our 0.35–0.65 over 50 nt), max homopolymer (3 vs 4), repeat *k* (16 vs 20), and
scope (42 modules from amino acids plus a coding mask vs the 21 needed, from DNA, interfaces
frozen). A looser GC band and a stricter homopolymer limit push in opposite directions.

What matched without adjustment: the forbidden-site set (SapI, BsaI, BpiI — BpiI is BbsI's
isoschizomer), genetic code 1, and the Pryor 2020 ligation data. Both inventories hold `AATG`
and `AGGT` variants, so `ppr_5prime_fusion_site` selects an assembly variant and does not
restrict the module-level comparison.

Its own log records four modules where it took a codon below its configured minimum relative
adaptiveness to deplete a SapI site. **Observed results under these inputs — not evidence
about the reference generally, and not evidence for us.**

**One failure here was ours.** The first attempt died on a `UnicodeEncodeError` because the
harness let a Windows subprocess inherit cp1252 while the reference printed an arrow. Recorded
as a harness fault, not a reference failure.

## 8. M3 gate decision (E5) — READ

**All five criteria met; the limited user-facing route is authorised.** Criterion 5 —
reproducibility, M1 and release preparation complete enough to protect the internal freeze —
was outstanding when the first four were settled and was closed later in the same pass by the
M1 exports, the manifest integration and the §10 release validation below. Full table:
`validation/experiments/m3_result.md`.

## 9. Integration work completed in this pass

- **M1 product exports.** `products.write_fasta` and `products.write_genbank`. Every FASTA
  header repeats the insert-only scope, because a FASTA travels without its JSON. The GenBank
  carries one feature per module (spans overlapping by the four junction bases, which is the
  assembly, not an off-by-one), one per junction, and a CDS **only** when a translation is
  supplied — no frame is invented. Verified: the CDS `/translation` is derivable from the
  bases under its own feature, and the final FASTA record equals the GenBank sequence.
- **Manifest in the result package.** `write_library` now writes `manifest.json`. A run where
  every target failed still gets a provenance record, because a directory of failures with no
  statement of what produced them is the case the manifest exists to prevent.
- **Batch homogeneity.** `manifest.build` reported the codon table, genetic code and enzyme
  profile once, read off the *first* design. It now refuses a heterogeneous batch and names
  the disagreement. Key order alone is not a disagreement.
- **Source fingerprint widened** from 25 `.py` files to 27, adding `parts.json` and
  `scaffold.json` — they determine a design as much as the code does. The falsifier is
  tested: adding a `.json` to the package changes the digest, and removing it restores it.

## 10. Release validation (§10)

    python validation/release_check.py --clean-install

Ten checks, each labelled RAN / READ / **NEEDS A HUMAN**. A check the script cannot evaluate
is never reported as passed. Result: `work/release/release_check.json`.

**10 of 10 passed.** Two of the ten failed on first run and both failures were in the checker,
not the project: `addopts = "-q"` in `pyproject.toml` meant passing `-q` again made it `-qq`,
which suppresses pytest's summary line entirely, and the `-k` filter for check 6 matched no
real test names. Both are recorded here because a checklist that goes green because the
checker was weaker than the claim is the exact failure this script exists to prevent.

Check 2 reports zero skipped tests, which is true **on this machine** because Supplementary
Table S1 is present. On a checkout without it the kit reconstruction tests skip with a stated
reason, and the check reports that count instead.

The clean-environment gate builds a wheel, installs it into a throwaway venv with a fresh
`CLIPPR_CACHE_DIR`, and runs the README's worked example **read out of the README itself**,
requiring its documented output verbatim. Confirmed: it reproduces

    AAAAUGUGG (9S) -> 302 aa, 906 nt, 4 fragments. Fidelity 0.828. QC PASS. 109.00 EUR (list price).

## 11. Findings from the orchestrator's verification, and what changed

The orchestrator reviewed the previous version of this handoff and found four validator
defects. All four were real; two were validators passing for the wrong reason, which is worse
than failing. Each is recorded with the probe that now holds it.

**V5 accepted incorrect products.** `cds in product` is containment, not correctness. Two
demonstrably wrong assemblies passed: one with its entry overhang changed from `CTCA` to
`AAAA`, and one with a `GGTCTC` inserted just inside the entry overhang. `digest()` also used
`find`/`rfind`, so extra recognition sites were invisible — an enzyme cuts at every site it
finds. Now: the product must **equal** `CTCA + CDS + CGAG` against the independently declared
destination; every site occurrence is counted and more than one pair is refused; the released
core is rechecked for further sites; and corrupt controls are judged by the same predicate as
real designs. **7/7 controls detected, corpus still 200/200 under the stricter gate.**

**V7 accepted missing and shifted annotations.** Stripping every GenBank feature gave zero
problems, because the loop iterated whatever existed and zero features means zero iterations.
Shifting the CDS three bases also gave zero problems, because no span was ever compared to an
expected coordinate. Now the feature layout is **derived** from the scaffold constants and
target length — 1 CDS + 2 + 3n + fragments, which reproduces the observed 34/51/67 exactly —
every promised feature must be present at its derived coordinates, the `/translation` must
match the bases under its own feature, and the two damaged files run through `inspect_exports`,
the same function real exports use. **5/5 controls detected.**

**Release checks were weaker than their labels.** `check_architectures()` would have accepted
a three-design smoke run. It now derives required coverage from the registered target list and
requires exact target identities, per-architecture counts, no duplicates and no failures. The
independent-agreement check read bare integers with nothing tying them to what produced them;
each validator now stamps a `provenance` block (corpus hash, its own file hash, package source
fingerprint) and a result that does not bind to the current tree is reported **STALE**, not
passed. V7's helper accepted an empty control dictionary; controls must now exist as well as
pass. Clean-install evidence records the wheel hash, package source, README hash and codon
table hash, drops inherited `PYTHONPATH`/`PYTHONHOME`, and verifies the import resolves
**inside** the throwaway environment.

**Budget accounting lost the interrupted item.** Checkpointing only between items meant a run
killed 30 s into its second item recorded 4 s, and the resume began with a falsely fresh
budget. Elapsed time is now persisted in a `finally`, covering every catchable exit, and each
checkpoint is labelled `exact` or `lower bound` — an abruptly killed process cannot write, so
its last value under-reports by up to one item. Fixing this exposed a second leak: the
returned total was read *after* the file closed and so always exceeded the persisted one, and
every resume inherited the difference. Both are now the same number. Three tests hold this.

**Bank identity ignored its source field.** A bank declaring `source_sha256='wrong-source'`
was compared as a matched run because only the scorer, matrix and codon table were checked.
`source_sha256` is now in that set, and the orchestrator's probe is refused.

Two corrections to this document's own claims: Table S1's hash was given as `1e0a3c74…`,
which is the package source fingerprint, not the file (`6324ede3…`); and "different sequence
every time" overstated a harness that counts modules receiving more than one distinct
sequence across five runs.

## 12. Known gaps

- **V1's matched-condition NEB Ligase Fidelity Viewer panel** has not been run. It needs a
  person at the NEB web tool; it is not automatable here.
- **V2 vendor screening** is untouched, and was declared optional and subordinate to
  correctness.
- Expansion of M3 beyond the 21 pilot modules to all 42 is a separate bounded extension and
  is not authorised by the M3 result.
- Nothing here is biologically validated. No binding was measured. Every fidelity number is a
  prediction from the Pryor 2020 data, and every sequence claim is a sequence measurement.
