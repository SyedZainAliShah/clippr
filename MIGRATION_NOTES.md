# Migration notes — the repair release

Behaviour that changed for callers. The earlier handoff claimed one behavioural break; that was
wrong, and the full list is below. Each entry says what used to happen, what happens now, and
what a caller should do.

## 1. A rejected design no longer produces normal order files

**Was:** `design_oneshot` wrote `oligo_csv`, `oligo_fasta`, `gene_fasta` and `genbank` for any
design that did not raise. `design_library` kept every such design, `LibraryResult.oligos()`
concatenated them, `write_library` exported them as `library_order_sheet.csv` and
`library_opool.csv`, and `design_for_synthesis` reported `ok=True`.

**Now:** every design carries a `final_artifacts` verdict judged on the sequence that would
actually be ordered. On a hard breach:

- `design_oneshot` writes the same four files under `diagnostic_oligo_csv`,
  `diagnostic_oligo_fasta`, `diagnostic_gene_fasta` and `diagnostic_genbank`, with `REJECTED`
  in every filename. No normal order-file key exists.
- `LibraryResult.oligos()` returns accepted designs only. `oligos(include_rejected=True)`
  returns everything with a `synthesis_ready` column.
- `write_library` writes `library_REJECTED_oligos.csv` and `library_REJECTED_reasons.json`
  instead of the order sheet.
- `design_for_synthesis` reports `ok=False` and lists each rejection in `failures`.

**If you consume these:** key off `result["synthesis_ready"]` or
`result["final_artifacts"]["status"]`, not the presence of files. A design with no recorded
assessment is treated as **not** ready; absent information is not readiness.

**Note:** all 200 saved corpus designs are rejected under the strict policy. This is a
correctly detected rejection against the project's own declared GC band, not a crash and not a
vendor refusal.

## 2. The joint inventory route refuses a non-standard genetic code

**Was:** `explore_interfaces` forwarded `genetic_code` to the search while completing the codon
table under a hardcoded `1`, so a non-standard code produced a front whose binding and whose
table disagreed. `select_interface` completed the table before reading the front's binding.

**Now:** both refuse any genetic code other than 1 with an explicit `WorkflowResult(ok=False)`
and write nothing. The code is not propagated through recoding, evaluation and export, so
honouring it partially would report objectives computed under a different code than declared.

**If you need another code:** use the one-shot route, whose own genetic-code support is
unchanged.

## 3. The weighted target-distance metric refuses targets longer than 10

**Was:** `metric="weighted"` accepted any length. Positions past 10 had weight zero, so two
14-mers differing only after position 10 scored weighted distance **0** against uniform **4**;
the 19-mer pair scored **0** against **9**.

**Now:** `_weights(length, "weighted")` raises `ValueError` above length 10 and states why. The
weights come from a ten-position editing assay; none exist for later positions and none are
invented.

**If you called it:** use `metric="uniform"`, which is the default and the supported route.
`distance()` also now validates its weight vector — a short vector used to compare only a
prefix and report it as the whole distance.

## 4. `level1_participants` never establishes participants from a flat role list

**Was:** `participants_established` could be `True` when end multiplicities looked like a
chain. That is a degree test, and an open path plus a disjoint cycle passes it.

**Now:** `topology` is always `"unverified"`, `participants_established` is always `False`, and
`scored` is `None` from this entry point. New keys: `topology`, `scored`, `malformed_ends`,
`end_multiplicities`, `degree_consistent_with_a_chain`. `workflow.level1_readiness` additionally
exposes `topology`, `why_incomplete` and `degree_consistent_with_a_chain` in `data`.

Level 0 reactions from `compile_target` are unaffected: their participants genuinely are known.

**If you relied on the score:** there is no replacement yet. A typed physical-fragment model is
what would restore it, and it does not exist.

## 5. `substrates.synthesis_findings` takes intended-site records

**Was:** `synthesis_findings(sequence, profile)` excused the wrapper's own sites by offset
alone and only for BbsI, so the one-shot BsaI route could not use it, and a forward site at the
reverse arm's offset was wrongly excused.

**Now:** `synthesis_findings(sequence, profile, *, intended_sites=None, enzymes=None)`.
Exemptions are `(enzyme, strand, offset)` records and all three must match. `intended_sites=None`
keeps the inventory's BbsI geometry, so existing callers are unchanged. `enzymes` supplies the
effective active set rather than always using the default profile.

`assembly.wrapper_intended_sites(payload, enzyme)` derives the records from the arms
`wrap_fragment` builds.

## 6. Target selections carry a reference content binding

**Was:** `Selection.reference` was a label. Two annotations can share an accession and differ.

**Now:** `Selection.reference_binding` holds a hash over the ordered feature payload screened —
index, name, kind, strand, start, normalised sequence, duplicates preserved — plus the
coordinate convention and a schema version. `targeting.verify_reference_binding(saved,
transcripts)` returns `match`, `mismatch`, `unverified` or `unreadable`.

Records written before this stay readable and are reported `unverified`. They cannot claim to
describe a current reference.

## 7. Screening and QC wording

No API change, but the text changed and tests asserting the old strings will fail:

- the universal "14-nt and 19-nt targets do not occur" is gone, with a counterexample in the
  regression;
- "there is no RNA for a PPR to bind at those loci" is gone;
- "not found in the host" is now "no exact match in the supplied reference";
- "occurrence is necessary for an off-target interaction" is gone — it implies absence excludes
  interaction, and a PPR tolerates mismatches;
- QC no longer claims to catch "sequences no vendor will make".

## 8. Validation harnesses

- `validation/notebook_on_published_build.py` now takes `--notebook`, `--workdir`, `--out` and
  repeatable `--input DEST=SRC`. It refuses a working directory it does not own, refuses staged
  inputs that escape it, refuses a wheel whose fingerprint differs from the source, fails any
  cell with a non-install magic, and fails a run that executed nothing.
- `validation/experiments/phasef_gate.py` gains `--keep-environment` and refuses to overwrite a
  previous run's retained environment or relocated package.
- `validation/release_check.py` gains `--evidence-dir` and `--refresh-state`, and checks 12 and
  13.
- `validation/independent/final_artifact_contract.py` is new.

---

# Second repair pass — 10 October 2026

The entries above cover the first pass. The independent review
`repair_review_20261010_revision2` found that the readiness boundary, the notebook gate and
the evidence binding were still incomplete. These are the further breaks.

## 9. One readiness predicate, and it reads the evidence

**Was:** two predicates disagreed. `library._ready` returned `synthesis_ready` when present
and otherwise fell back to `status != "fail"`, which admitted `not_assessed` and any
unrecognised status. `design_oneshot` chose the diagnostic export route only on
`status == "fail"` and then computed readiness from a stricter conjunction. A design whose
status was `not_assessed` was therefore reported unready **and** given a full set of normal
order files.

**Now:** `design.synthesis_readiness(result) -> (bool, reasons)` is the only predicate, and
both paths call it. It accepts `design.ORDERABLE_STATUSES` — `("pass", "warning")` — with
satisfied optimiser constraints and no hard findings. Everything else refuses, with a
reason: an absent assessment, an absent constraint result, an unknown status string, a hard
breach recorded beside a `pass`, and a `synthesis_ready` flag that disagrees with the record.
`design_oneshot` now also returns `not_ready_because`.

**If you consume these:** read `synthesis_ready` or call `synthesis_readiness`. Do not test
`final_artifacts["status"] != "fail"`; that is the defect. Note that a stale or hand-edited
`synthesis_ready=True` no longer wins over the record.

## 10. Saved summaries, QC tables and manifests state the synthesis verdict

**Was:** `qc_table()` had a column `qc` holding the coding span's verdict, so a rejected
design showed `PASS` with nothing beside it. `summary()` led with `QC 1 PASS`. `write_library`
derived `completion` from crashed designs alone, so a library in which every design was
correctly refused shipped a manifest saying `"complete"` with `failures: []`.

**Now:**

- the QC column is `cds_qc`, and `final_artifacts` and `synthesis_ready` sit beside it;
- `summary()` leads with `synthesis  N of M designs are orderable` and lists each
  `NOT ORDERABLE <target>: <reason>`;
- the manifest carries `synthesis_coverage` — `requested`, `accepted`, `rejected`,
  `not_assessed`, `failed` — and `rejections`, and `completion` is `"partial"` when anything
  was rejected.

**If you read the QC CSV:** the column name changed. The design dict's own `result["qc"]`
key is unchanged.

## 11. The notebook harness requires declared inputs and results

**Was:** `INSTALL_MAGICS` contained `%%capture` and `%load_ext`, and a cell containing any
permitted magic was skipped **whole**. A required cell beginning `%%capture` and then raising
reported "1 skipped" and PASS. Separately, the harness counted executed cells, so the
inventory notebook with Table S1 absent guarded every step on `deposited is None`, executed
ten cells, did no workflow, produced no package and reported `ok=True`.

**Now:**

- skipping is a property of a **line**. `SKIPPABLE_LINES` is `%pip`, `!pip`, `!git`,
  `%load_ext`; `SUPPORTED_CELL_MAGICS` is `%%capture`, honoured by executing the body
  underneath it. A cell is skipped only when nothing but permitted magic lines remains.
  Anything else — `%time`, `%%time`, `%matplotlib`, `%cd`, `!rm` — fails the cell.
  `classify_cell(src)` exposes this and is unit-tested;
- `--require-input DEST` refuses to start unless that input was staged, exit 2;
- `--expect EXPR` declares a postcondition evaluated in the notebook's namespace after the
  run. All must hold or the run fails. `ok` is now
  `failed == 0 and ran > 0 and not unmet`, and the record carries `expectations` and
  `postconditions_declared`.

**If you invoke it:** a run with no `--expect` can only report that cells executed, and
release check 13 now refuses such a record.

## 12. The independent corpus checker fails on missing data and short panels

**Was:** `return 1 if disagreements else 0`. A corpus of unreadable rows produced
`not_assessed` counts and exit 0 — agreement that nothing could be assessed read as success.

**Now:** any `not_assessed` design exits 1. `--expect-designs N` declares the panel this run
must cover and a shorter run exits 1. A **rejection** still exits 0: that is a verdict,
correctly reached, and must not be confused with a verifier failure.

## 13. Release check 13 is bound to the candidate release

**Was:** the two notebook fingerprints were compared only `if both exist`, so a record with
neither passed and a record reduced to `{"ok": true}` passed with no identity at all.
`comparisons_performed > 0` accepted one comparison for a two-hundred-design panel. Check 13
was also appended **after** the record was written, so `--refresh-state` rendered the State
table from an interim eleven-check file and the final write immediately invalidated it: the
table said "11 of 11" while a separate read-only check measured 12 of 12.

**Now:**

- a notebook record must carry `notebook_sha256`, both fingerprints, `interpreter`,
  `environment_prefix` and at least one declared postcondition, with none unmet;
- the installed fingerprint is compared against **this checkout**, not only against its
  neighbour in the same file;
- the corpus report must cover the declared 200 designs, compare every row it lists, and
  report no `not_assessed`;
- check 13 is computed **before** the record is written. Only check 11, which reads the
  record, runs after it, and the renderer excludes check 11, so the rendered tally stays
  true;
- `state_table.py --write`'s exit code is honoured. A renderer failure fails the release
  instead of being followed by a verification of whatever the file already said.

---

# Third repair pass — 10 October 2026

The review `repair_review_20261010_revision3` accepted the export predicate, the real
notebook execution and the state-table repair, and found three narrower issues: the
user-facing display still led with a green coding-span `PASS`, release check 13 checked
field presence rather than artifact identity, and the successor document overstated four
claims. These are the further breaks.

## 14. Every user-facing surface states the synthesis verdict

**Was:** the API refused to export a rejected design while every surface a person reads led
with the coding-span QC.

- the notebook result card took its badge text and its colour from `result["qc"]["status"]`,
  so a design the final-artifact check rejects showed a green **PASS** and the rejection
  appeared nowhere on the card;
- `DesignAudit.report()` printed `constraints satisfied` and `QC PASS` with no final verdict;
- `result["summary"]` — the README quickstart's whole printed output — ended `QC PASS.`;
- the fragment section was headed "The fragments to order" and called each row an
  "orderable piece";
- the download table labelled `diagnostic_*` keys with their raw key names.

**Now:**

- the badge is the readiness verdict, `READY TO ORDER` or `NOT ORDERABLE`, and a refused
  design gets a panel naming each reason and saying the files are not an order. The card
  carries `coding-span QC` and `final artifacts` as separate, scoped statistics, and the GC
  and repeat figures say `coding span`;
- `DesignAudit` gains `final_status`, `synthesis_ready` and `not_ready_because`, and
  `report()` leads with `SYNTHESIS  NOT ORDERABLE` and its reasons. The retained CDS result
  is kept and labelled `coding-span QC … (the CDS only, not the ordered sequences)`, and the
  optimiser line says `(the coding span)`;
- `summary` ends `Coding-span QC PASS. NOT ORDERABLE (final artifacts: fail).` for a refused
  design and `… Orderable under the declared software checks.` for an accepted one;
- the fragment section says a row is not an order; the download table labels diagnostic keys
  `DIAGNOSTIC ONLY - rejected …` and prefixes a refusal banner.

**If you depend on these strings:** `result["summary"]` changed, and so did the README
quickstart comment that release check 3 parses as its expected clean-install output — the two
are kept in step by `tests/test_display_verdict.py`. The audit's `QC` line is now
`coding-span QC`.

## 15. Release check 13 recomputes identities instead of reading fields

**Was:** the check required certain fields to be non-empty and never recomputed them. Each of
these passed: a notebook record whose `notebook_sha256` was sixty-four zeros; a staged Table
S1 whose recorded hash was zeros; a Gate F record whose installed-source fingerprint and
wheel hash were zeros; a corpus report whose `corpus_sha256` was zeros; a notebook record
with `expectations: []` and `postconditions_declared: 8`; and a corpus report whose two
hundred rows were the first row copied two hundred times. Only four record types were
required, so a bundle holding Gate F, two notebook reports and the corpus report passed with
no B–E or V5–V7 evidence at all.

**Now:**

- `REQUIRED_RECORDS` names eleven records — gates B–F, both notebooks, the corpus contract
  and V5–V7 — and a bundle missing any one of them fails;
- every stated identity is recomputed from the file it names: the notebook bytes, each staged
  input's source bytes, Gate F's retained wheel, the corpus file, and each validator script.
  A record that names no readable file to verify against also fails;
- `postconditions_declared` must equal the length of `expectations`, and `NOTEBOOK_CONTRACT`
  names the result each notebook must actually have asserted, so `--expect "1"` no longer
  satisfies the gate;
- the corpus report's row targets must match the declared panel as a multiset, repeats
  included, so a duplicated panel fails while a legitimately repeated target does not;
- the corpus report and V5–V7 records carry `provenance`, including the validator's own
  SHA-256, because the package fingerprint does not cover `validation/*.py`. Editing a
  checker invalidates results the previous version produced.

`validation/independent/v5_digest_ligate.py`, `v6_constraints.py` and
`v7_export_semantics.py` gain `--out`, so a bundle can hold its own copy; the shared
`work/independent/` path stays as the default. `check_independent_agreement` reads the bundle
when one is given.

## 16. The release record is frozen into the bundle

**Was:** every figure in a release document resolved to `work/release/release_check.json`,
which later runs rewrite. A renderer-failure control overwrote it with its own failing
record once.

**Now:** with `--evidence-dir`, the aggregate writes a copy to `<bundle>/release_check.json`
and its SHA-256 to `<bundle>/release_check.sha256`. The shared path remains as a convenience
and is explicitly mutable.

## 17. The notebook harness takes `--env`

**Was:** the harness set no environment, so the installed wheel resolved `cache_dir()` to a
per-user directory. The designer notebook consumed a codon table and a chloroplast genome
from `%LOCALAPPDATA%\clippr` while its record declared no inputs at all.

**Now:** `--env NAME=VALUE` is applied before the build is imported, and a relative value
resolves inside the working directory and may not escape it. Pointing `CLIPPR_CACHE_DIR` at a
staged cache makes every byte the run consumes a declared, hashed input. The record carries
`environment_overrides`.

**If you invoke the harness:** both documented notebook commands now stage their codon table,
and the designer also stages the genome. See section 8 of `RELEASE_STATUS_2026-10-10.md`.

## 18. V6's reported scope is named

**Was:** release check 5 printed "200/200 satisfy every hard constraint". V6 checks the
coding sequence; the final wrapped-oligo check rejects the same corpus. The unscoped sentence
read as both.

**Now:** "200/200 satisfy the V6 coding-sequence constraints". No behavioural change in V6.

---

# Fourth pass — 10 October 2026

## 19. `optimize_cds` accepts caller-supplied constraints

**Was:** the solver's constraint list was fixed. An experiment that needed the solver to see
something the shipped rules do not — the sequence that would actually be ordered, rather than
the coding span — had no way in without editing a default.

**Now:** `optimize_cds(..., extra_constraints=())` takes further DNA Chisel specifications
from the caller. The default is empty and **no shipped route passes it**;
`tests/test_extra_constraints.py` asserts that `extra_constraints` appears in no `src/` module
other than its own definition, so a shipped route starting to use it becomes a visible test
failure rather than a silent default change.

**If you call it:** an unsatisfiable extra constraint is reported as `constraints_ok=False`,
not dropped, and translation is still enforced — the hook cannot buy a constraint win with the
wrong protein.

The wrapper-aware pilot in `validation/experiments/wrapper_aware_pilot.py` is the only user.
See `PILOT_WRAPPER_AWARE_2026-10-10.md` for its result and its negative control.

## 20. The one-shot summary and the audit lead with the synthesis verdict

Covered in section 14. Repeated here because it changes two strings a caller may have pinned:

- `result["summary"]` now ends `Coding-span QC PASS. NOT ORDERABLE (final artifacts: fail).`
  for a refused design, or `… Orderable under the declared software checks.` for an accepted
  one. The README quickstart comment that release check 3 parses as its expected
  clean-install output changed with it, and `tests/test_display_verdict.py` keeps the two in
  step.
- `DesignAudit.report()` leads with `SYNTHESIS  NOT ORDERABLE` and its reasons, and its `QC`
  line is now `coding-span QC`.
