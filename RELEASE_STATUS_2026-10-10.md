# PPROPT release status — 10 October 2026, revision 2

**This document supersedes the release-status paragraphs in the frozen audit pack**
(`docs/*.md`, `RESPONSE_TO_AUDIT.md`, `PROPOSAL_TO_CODEX.md`) and revision 1 of this file.
Those files are left unchanged as the historical record. Where they disagree with this
document, this document is current.

Revision 1 overstated four things. The independent review was right about all four and they
are corrected in place, with the correction named: see §10.

Nothing here has been committed, pushed, published, ordered or sent to a third party. All
changes are local and uncommitted on `479690b53663c397926aa431e72f4965dfbe4998`.

**RAN** = executed for this document. **READ** = inspected a retained record.

---

## 1. What the software does today, in one paragraph

PPROPT turns a short target RNA into a synthetic PPR protein's coding sequence, splits it
into Golden Gate fragments, and checks the sequences that would actually be ordered. It has
two routes. The **inventory route** takes the deposited 42-module GRASP kit, recodes it for a
host, picks a junction interface and produces an order plan — this route works end to end and
is demonstrated in §4. The **one-shot de novo route** designs a fresh coding sequence; all 200
saved corpus designs were rejected on revalidation under this policy, and the demonstrated
default example is also rejected, because the optimiser constrains the coding span while the
ordered fragment also carries an entry overhang and an enzyme wrapper. That rejection is the
honest output, not a crash.

## 2. Identity of this release candidate

| | |
|---|---|
| checkout HEAD | `479690b53663c397926aa431e72f4965dfbe4998`, with uncommitted changes |
| package source fingerprint | `0f4d63540cd2c291e6c85185b8d2c58232b0d300f2df37060e9f465135bb14ed` |
| wheel SHA-256 | `db2a77070c628ff9e82ff95c7c63867604dff7a95da5789ad0b0cdb63a646e90` |
| evidence bundle | `work/repair_release_20261010-232032/` |
| designer notebook SHA-256 | `b92981317002794d…` |
| inventory notebook SHA-256 | `840105d98522a0fb…` |
| corpus checker SHA-256 | `c84873bef0494008…` |
| pinned corpus SHA-256 | `8e93d401be3b22cf…` |

Earlier bundles (`repair_release_20261009-221835`, `20261010-044742`, `20261010-074942`) are
preserved unchanged, including the `before/` snapshot of the pre-repair records.

**The bundle is self-contained.** It holds all eleven required records — gates B–F, both
notebook reports, the corpus contract, and V5–V7 — plus a frozen copy of the aggregate
record as `release_check.json`, whose hash
`d6e40803d2ce45b3ea2dba821060060da86f550f3fee8402d6c7b21f887af6f2` is beside it in
`release_check.sha256`. The shared paths
`work/release/` and `work/independent/` remain as conveniences and are **mutable**: a control
run overwrote `work/release/release_check.json` once, which is why the bundle now carries its
own copy.

## 3. Measured results

Every figure below resolves to a named record **inside the bundle**.

| Result | Value | Record in the bundle |
|---|---|---|
| full test suite, RAN | **1108 passed of 1108 collected, 0 failed** | `release_check.json`, check 1 |
| aggregate release check, RAN | **13 of 13 passed, 0 failed, 0 needing a human** | same file |
| clean-environment install, RAN | example reproduced its documented output verbatim | check 3 |
| Gate F, RAN | **PASS**, 12 artefacts re-verified after relocation | `phasef/gate_f.json` |
| Gates B, C, D, E, RAN | **PASS** | `gate_b.json`, `phase{c,d,e}/gate_*.json` |
| V5 digest/ligate, RAN | 200/200 reconstructed, 7/7 corrupt controls caught | `v5_digest_ligate.json` |
| V6 constraints, RAN | **200/200 satisfy the V6 coding-sequence constraints** | `v6_constraints.json` |
| V7 export semantics, RAN | 3 designs, 0 semantic problems, 5/5 controls caught | `v7_export_semantics.json` |
| designer notebook, RAN | 14 cells ran, 0 failed, **8 postconditions met**, 2 inputs staged | `notebook_designer.json` |
| inventory notebook, RAN | 11 cells ran, 0 failed, **8 postconditions met**, 2 inputs staged | `notebook_inventory.json` |
| corpus contract, RAN | 200 designs: 0 accepted, **200 rejected**, 0 not assessed; 200 comparisons, **0 disagreements** | `final_artifact_contract.json` |
| fragment breaches, RAN | 1050 fragments checked, **277 with a hard finding**, 0 unusable | same file |

**V6 and the corpus contract have different scopes and both results are true.** V6 checks the
**coding sequence** — its GC windows, homopolymers, blacklisted sites, translation and cut
geometry. The corpus contract checks the **final wrapped oligo**, which is the entry overhang
plus the coding slice plus the enzyme wrapper. So: 200/200 satisfy the V6 coding-sequence
constraints, and 200/200 fail the separate final-oligo policy check. Reading either number
without its scope is how a rejected design came to look acceptable.

## 4. The inventory route — what can be shown

RAN, in the inventory notebook on the wheel in §2, with Table S1 and the codon table staged:

```
42 assembly-ready substrates from selected-recommended@a6c3074572de546c
order plan : 42 oligos, eligible, 1 pool(s)
package    : 8 tasks, ok=True, 0 failure(s)
```

Honest wording for a talk or a wiki page:

> The software can take the deposited module kit, recode it for a host, choose a junction
> interface and produce a local order plan — 42 assembly-ready substrates in one pool. The
> eligibility result is a check against a vendor profile we transcribed ourselves. No order
> was placed, no vendor approved these sequences, and nothing here measures binding,
> expression or assembly yield.

Two limits to state alongside it:

- the vendor profile is **transcribed from published specifications on 2026-09-15** and only
  the 50 pmol scale is modelled. Confirm against the current page before ordering;
- the pool plan is **a feasible plan, not a proven cost optimum**, and the interface search
  reports an **observed front** within a recorded budget, not an exhaustive optimum.

## 5. The one-shot route — the rejection, and why

RAN, in the designer notebook on the installed wheel. The notebook's own postconditions
assert that the rejection is visible, not merely recorded:

```
AAAAUGUGG final-artifact status : fail
synthesis_ready                : False
hard  : AAAAUGUGG_F1: GC 0.660 outside (0.35, 0.65) in the window at 1
result card badge              : NOT ORDERABLE
exported as                    : diagnostic_oligo_csv  (no oligo_csv key)
```

Cause, traced: the optimiser enforces local GC on the **coding span**. The sequence that
would be ordered is the destination's entry overhang, plus the CDS slice, plus the enzyme
wrapper. On this design the CDS 50-nt windows run 0.540–0.640, inside the band, while the
ordered oligo reaches 0.660 at offset 1. The payload begins `CTCA`, which no CDS window
covered.

Scale: **277 of 1050 fragments** carry a breach; 268 worst windows at 0.66 and 9 at 0.68
against a 0.65 ceiling. **77 of them are F2–F7** (F2=5, F3=22, F4=24, F5=2, F6=12, F7=12), so
fixing only the outer gene boundaries cannot settle this.

**The threshold was not weakened and no optimiser default was changed.** The decision taken
is to keep the strict band and publish the counts. Whether a wrapper-aware optimiser can
recover a usable one-shot route is a separate, predeclared experiment that has **not** been
attempted.

Honest wording:

> Our current de novo designs fail our own strict final-sequence rule, because the optimiser
> does not yet account for every fragment's construction context. We publish the rejection
> rather than widening the limit.

Do not present the failed one-shot sequence as an accepted synthesis design.

## 6. What changed for anyone reading an older document

`MIGRATION_NOTES.md` is the full list. The five that change how a result should be read:

1. **A rejected or unassessed design no longer produces normal order files.** One predicate,
   `design.synthesis_readiness`, decides both the export route and the reported verdict. Only
   `pass` and `warning` with satisfied optimiser constraints are orderable; a missing
   assessment, an unknown status, and a `synthesis_ready` flag that disagrees with the record
   are all refusals. Rejected artefacts are still written, under `diagnostic_` keys with
   `REJECTED` filenames.
2. **Every user-facing surface states the synthesis verdict.** The notebook result card's
   badge is the readiness verdict, not the coding-span QC; it names the breach and says the
   files are not an order. `audit.report()` leads with `SYNTHESIS NOT ORDERABLE` and its
   reasons, and labels the CDS result `coding-span QC`. `result["summary"]` — the README
   quickstart's whole printed output — carries the verdict too.
3. **Saved summaries, QC tables and manifests state it as well.** The QC column is `cds_qc`,
   with `final_artifacts` and `synthesis_ready` beside it. A manifest carries
   `synthesis_coverage` and `rejections`, and `completion` is `partial` when anything was
   rejected.
4. **A notebook validation run must declare its inputs and its results.** The harness takes
   `--require-input`, `--expect` and `--env`; a run that executes cleanly but builds no result
   fails, and a run whose declared inputs were not staged refuses to start. Permitted magic
   lines are stripped from a cell rather than skipping the cell.
5. **`level1_participants` never establishes participants from a flat role list.** `topology`
   is always `"unverified"` and there is no score from that entry point. Level 0 reactions
   from `compile_target` are unaffected.

## 7. What is not assessed

- **Whether any vendor would refuse these sequences.** A failed engineering-policy check is
  not a vendor refusal. No order was submitted and no quote obtained.
- **Whether a wrapper-aware optimiser can recover a usable one-shot route.** Not attempted.
  There is **no positive final-export control under the strict policy**, so no successful de
  novo synthesis route should be claimed. Its runtime is unmeasured; nothing in this run
  predicts it.
- **Browser or Colab execution.** Only installed local notebook execution was verified.
- **Level-1 assembly context.** Participants are not established from a flat role list, so no
  level-1 score is reported.
- **The exact historical bytes behind the `c76dc3ca…` fingerprint.** The line-ending mechanism
  is demonstrated; the historical mix is not recovered. The defensible claim is the weaker
  one: a source edit is not needed to explain it.
- **Any biological property.** Nothing here measures binding, expression, regulation,
  specificity or assembly yield.

## 8. Reproduction

From the repository root. `$Run` is a **new** directory: Gate F refuses to overwrite a
previous run's retained environment, which is deliberate.

```powershell
Set-Location -LiteralPath 'C:\CodePractice\Drylab\ppropt'
$Py  = '.\.venv\Scripts\python.exe'
$Run = "work\release_$(Get-Date -Format yyyyMMdd-HHmmss)"
New-Item -ItemType Directory -Force $Run | Out-Null

& $Py validation\experiments\phaseb_gate.py --out "$Run\gate_b.json"
& $Py validation\experiments\phasec_gate.py --out "$Run\phasec"
& $Py validation\experiments\phased_gate.py --out "$Run\phased"
& $Py validation\experiments\phasee_gate.py --out "$Run\phasee"
& $Py validation\experiments\phasef_gate.py --out "$Run\phasef" --keep-environment
& $Py validation\independent\v5_digest_ligate.py     --out "$Run\v5_digest_ligate.json"
& $Py validation\independent\v6_constraints.py       --out "$Run\v6_constraints.json"
& $Py validation\independent\v7_export_semantics.py  --out "$Run\v7_export_semantics.json"
& $Py validation\independent\final_artifact_contract.py `
    --jsonl work\independent\corpus200_pinned.jsonl `
    --out "$Run\final_artifact_contract.json" --expect-designs 200
```

### The notebooks, and the data they consume

Both notebooks run on the wheel Gate F just built. `ipython` is installed afterwards, so it
is recorded as `notebook_runtime_dependency` and sits deliberately outside the clean-install
claim.

**Every substantive input is staged from a repository-relative path and hashed into the
record.** Without this the installed wheel resolves its cache to a per-user directory, and
the bytes a run consumed appear in no record — the designer notebook previously declared no
inputs at all while reading a codon table and a chloroplast genome from that cache. `--env`
points `PPROPT_CACHE_DIR` at a staged cache inside the working directory.

```powershell
$Env = "$Run\phasef\environment\env\Scripts\python.exe"
& $Env -m pip install -q ipython

& $Env validation\notebook_on_published_build.py `
    --notebook notebooks\PPROPT_designer.ipynb `
    --workdir "$Run\nb_designer" --out "$Run\notebook_designer.json" `
    --env "PPROPT_CACHE_DIR=cache" `
    --input "cache/codon_tables/kazusa_3055.json=data\codon_tables\kazusa_3055.json" `
    --input "cache/genomes/NC_005353.1.gb=data\genomes\NC_005353.1.gb" `
    --require-input "cache/codon_tables/kazusa_3055.json" `
    --require-input "cache/genomes/NC_005353.1.gb" `
    --expect "result['target_rna'] == 'AAAAUGUGG'" `
    --expect "len(result['oligos']) == 4" `
    --expect "result['final_artifacts']['status'] == 'fail'" `
    --expect "result['synthesis_ready'] is False" `
    --expect "'diagnostic_oligo_csv' in result['paths'] and 'oligo_csv' not in result['paths']" `
    --expect "bool(result['not_ready_because'])" `
    --expect "'NOT ORDERABLE' in result['audit'].report()" `
    --expect "'NOT ORDERABLE' in card"

& $Env validation\notebook_on_published_build.py `
    --notebook notebooks\PPROPT_inventory.ipynb `
    --workdir "$Run\nb_inventory" --out "$Run\notebook_inventory.json" `
    --env "PPROPT_CACHE_DIR=cache" `
    --input "Table S1.xlsx=data\grasp_supp\Table S1.xlsx" `
    --input "cache/codon_tables/kazusa_3055.json=data\codon_tables\kazusa_3055.json" `
    --require-input "Table S1.xlsx" `
    --require-input "cache/codon_tables/kazusa_3055.json" `
    --expect "deposited is not None and len(deposited) == 42" `
    --expect "order_items.ok is True" `
    --expect "'42 assembly-ready substrates' in order_items.summary" `
    --expect "order.ok is True" `
    --expect "'42 oligos' in order.summary and 'eligible' in order.summary" `
    --expect "summary['ok'] is True" `
    --expect "len(summary['failures']) == 0" `
    --expect "len(summary['tasks']) == 8"
```

**None of these three files is in the repository.** `data/` is gitignored, and the wheel
carries only the two Pryor fidelity matrices. The codon table and the genome are fetched once
from Kazusa and NCBI and cached; Table S1 is downloaded from the GRASP paper's supplementary
material. A fresh clone must obtain all three before the commands above will run.

The bytes this release consumed:

| Staged as | SHA-256 | Where it comes from |
|---|---|---|
| `Table S1.xlsx` | `6324ede3e3863ced60f12cc822ebe68f530e9ebc6a5ebc45c490a8816930ef35` | GRASP supplementary material, `10.1093/nar/gkaf1169` |
| `cache/codon_tables/kazusa_3055.json` | `a9d8f6fc88e15694be1473123ff3dc9cb7e501777c462f3d244ec7368fb78ab3` | Kazusa, via `table_from_kazusa(3055)` |
| `cache/genomes/NC_005353.1.gb` | `d39b1657055046b98e05d931f27aebe0a87b606a4b601950620d30af36ebe5ec` | NCBI, via `offtarget.scan()` |

Check your copies against these before concluding a difference is a code change. If a future
run needs an input not listed here, fetch it once, record the hash the run consumed, and stage
that file; do not let the run resolve it silently from a user cache.

Then the aggregate, which runs the suite itself, and the separate read-only state check that
is the acceptance gate:

```powershell
& $Py validation\release_check.py --evidence-dir $Run --refresh-state --clean-install
& $Py validation\state_table.py
```

The last command must exit 0 on its own. It reads the record the aggregate wrote; running it
separately is what catches a tally the aggregate invalidated after rendering it. The aggregate
also freezes its record into `$Run\release_check.json` with a `.sha256` beside it.

Do not run the full suite manually and then again inside the aggregate to obtain the same
number. Check 1 is the measurement.

## 9. R10 deliverables — where each one is, and what is still open

| R10 requirement | Where it is now | State |
|---|---|---|
| Technical: status semantics | §6 here; `MIGRATION_NOTES.md` 1, 9, 10; `design.ORDERABLE_STATUSES` and `design.synthesis_readiness` | **done** |
| Technical: schemas | `FINAL_ARTIFACT_SCHEMA`, `REFERENCE_BINDING_SCHEMA`, `manifest.SCHEMA_VERSION`; record schemas in `release_check.REQUIRED_RECORDS` and `NOTEBOOK_CONTRACT` | **done** |
| Technical: coordinates | `docs/workflow_contract.md` (unchanged and still valid); `targeting.verify_reference_binding` coordinate convention | **existing, valid** |
| Technical: cache | §8 here, "the notebooks, and the data they consume"; `paths.cache_dir` docstring | **done** |
| Wiki: problem and method | `docs/methodology_sources.md`, `docs/gc_band_and_constraint_model.md` | **existing, valid** |
| Wiki: results | §3, §4, §5 here | **done** |
| Wiki: limitations | §7 here | **done** |
| Wiki: next work | §5 here and the predeclared pilot in `REPAIR_RELEASE_HANDOFF.md` §3 | **done** |
| Saved demonstration | the two notebook records in the bundle, with their postconditions and staged inputs | **done** |
| Conversation script | §4 and §5 here, the two "honest wording" blocks | **done** |
| Publication checklist | §11 here | **partial** — the checklist exists; the official requirements are **not** verified, see below |
| Attribution and input rights | `NOTICE.md` | **done for the software** — see below |

**Correction.** Revision 1 of this file, and handoff revision 4, said the attribution and
input-rights statement was unwritten. That was wrong. `NOTICE.md` covers it: the two Pryor
matrices are redistributed under CC BY 4.0 with attribution; the GRASP module index is a
derived index that deliberately omits insert sequences; the unlicensed PPRmatcher score table
is **not** redistributed and must be supplied by the user. The supplementary tables, codon
tables and genomes in `data/` are gitignored and absent from the wheel, so nothing from the
paper is republished.

**Open, and not closed by this document:** the wiki-facing attribution text, which is a
separate deliverable from `NOTICE.md`, and verification of the current official publication
requirements. R10 is therefore **not** wholly closed. The frozen audit pack is referenced
rather than duplicated or edited.

## 10. What revision 1 of this file claimed, and the correction

| Revision 1 said | Correct statement |
|---|---|
| "every current design is rejected" | All 200 saved corpus designs were rejected on revalidation under this policy, and the demonstrated default one-shot example is also rejected. The corpus is the preserved historical panel, **not** newly generated by this release |
| "V6: 200/200 satisfy every hard constraint" | 200/200 satisfy the V6 **coding-sequence** constraints; 200/200 fail the separate **final-oligo** policy check. Both are true; see §3 |
| "Every figure is from the evidence directory" | Revision 1 listed shared `work/release` and `work/independent` records. The bundle now holds its own copy of every record, including a frozen `release_check.json` with its hash |
| reproduction commands, inputs | Revision 1's notebook commands declared no codon table and no genome, so the designer run consumed a per-user cache no record named. §8 now stages and hashes every substantive input |

## 11. Before any of this is published

- Revise the release-status paragraphs in `docs/*.md` from this document. They predate the
  repairs and describe a state that no longer holds.
- Write the wiki-facing attribution text. `NOTICE.md` already carries the software's
  third-party material and input rights; a wiki page is a separate deliverable.
- State the scope on every figure: a local check against a transcribed vendor profile is not
  vendor acceptance, a coding-sequence check is not a final-oligo check, and a passing
  software gate is not a biological result.
- Check the current official publication requirements at the time of publication. **This
  document does not certify them.**
