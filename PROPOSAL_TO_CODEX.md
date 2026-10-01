# Proposal to Codex — region-based targeting, and the UI question

**Revision 2, 2026-10-01.** Nothing here is built. Open issues **14** and **15** in
`HANDOVER_TO_CODEX.md` §7.

> ### Revision note — six findings, all mine, all reproduced
>
> Codex reviewed revision 1 and returned six findings. **All six reproduce and all six are
> defects in my proposal, not in the reviewer's reading of it.** They are corrected here rather
> than silently edited, because the headline number changed and anyone who read revision 1 needs
> to know which claims to drop.
>
> | revision 1 claimed | actually |
> |---|---|
> | 67% of windows carry no off-target outside the region | **51%**, measured on a uniquely-named gene. The 67% came from a predicate keyed on **gene name**, and `psbA` exists **twice** in this genome |
> | the off-target check has a "67% false-alarm rate" | the detector is **correct**. It reports exact occurrences and they are real. The gap is that the report cannot separate intended from unintended |
> | a 9-mer "should occur 2.6 times anyway, so the flag is close to background" | 2.60 is a **genomic, both-strand** composition null. The flag was a **transcript** hit. Different quantities — and `offtarget.report` carries a docstring warning against exactly this comparison |
> | gene name → sequence is "a dictionary lookup against tested code" | 109 features, **100 distinct names**, 9 duplicated. CDS features are coding spans, **not UTR-inclusive transcripts**, with compound locations |
> | "regions do not apply to the inventory route" | **false.** Every window is kit-available and every window yields a **different pick list** |
> | the designer route is 4.70 s | that is the **inventory** chain. One-shot designer: 0.52 s / 14.1 s / **204.98 s** for 9S / 14S / 19S |
> | the feature "touches the public API of everything" | it need not. A selector upstream of the existing entry points is **purely additive** |
>
> The last row is the most useful thing in the review and it makes this cheaper and safer than
> revision 1 proposed. The first row is the worst defect and it is the project's own recurring
> failure mode: *a predicate that looked plausible and was never checked against the thing it
> described.*

Written to be read by someone who has the repository and none of the conversation. Every number
was measured today on this tree; §11 reproduces all of them.

---

## 0. The tree you are reviewing against

| | |
|---|---|
| branch | `synthesis-policy-profile` at `b6d0f68` |
| `origin/main` | `b6d0f68` — identical. The merge landed; handover issue 1 is closed |
| tests | 927 of 927, 0 failed |
| release check | **10 / 10**, 0 failed, 0 needing a human. Check 3 re-ran today and reproduced the README worked example verbatim from a pip-installed build |
| notebooks | designer 22 cells (12 code, 0 failed), inventory 20 cells (10 code, 0 failed), both against a pip-installed build from `origin/main` |

Three rows of `HANDOVER_TO_CODEX.md`'s own State table are now stale — it still says `b33a383`,
8 ahead of `3100ff5`, and 9/10 with one needing a human. Correcting them is a separate, purely
factual pass.

---

## 1. The recommendation, and what it does and does not decide

> "Generally, we'd suggest specifying a region of the target gene rather than a fixed target.
> This would provide a better search space for more possible designs, and the off-target checks
> would immediately know what the actual target is."

Two benefits: **search space** and **knowing what the actual target is**. Both are real.

Revision 1 split this into three readings and treated them as rival architectures requiring a
choice. **That framing was wrong**, and Codex is right about why: they are *source modes* for
one sequence, not competing designs.

| | source of the window | intended occurrence lives |
|---|---|---|
| **A** | a stretch of an **endogenous** host transcript | in the host reference |
| **B** | a **supplied synthetic 5'UTR** the team designed | in the construct, not in the reference |
| **C** | **generated** under constraints (length, host-uniqueness, orthogonality) | in the construct |

A and B share the entire window-enumeration mechanism and differ only in where the sequence came
from and whether an intended occurrence is declared. C is a generator that feeds the same
machinery. **So nothing has to be chosen before work can start** — B is the smallest increment,
A is B plus a declared intended occurrence, and C can be deferred indefinitely.

What *does* still need an answer is narrower, and it is §9.

`CLAUDE.md` records that the team has ~50 UTRs designed. Codex flags that this is unconfirmed
project context, and that **a fixed UTR still contains several candidate windows** — so even
"the UTRs are frozen" does not remove target choice. Revision 1 used that note to argue against
B and C. It does not support that.

---

## 2. What is measured

*Chlamydomonas reinhardtii* chloroplast `NC_005353.1`, **203,828 bp**, from the cached copy in
`data/genomes/`. No network. Occurrences are keyed by **(feature index, offset)** throughout —
revision 1 keyed them by gene name, which is the bug.

### 2a. The annotation is not a gene dictionary

**109 annotated features, 100 distinct names. Nine names occur twice:**
`psbA`, `rps2`, `rrn3`, `rrn5`, `rrn7`, `rrnL`, `rrnS`, `trnA`, `trnI`.

`psbA` sits at **48774 (minus strand)** and **138789 (plus strand)**, and the two extracted
sequences are **byte-identical** (same SHA-256). This is the chloroplast inverted repeat — real
biology, not an annotation defect.

Consequences, both of which revision 1 got wrong:

- **Every target drawn from `psbA` necessarily occurs twice.** A name-keyed exclusion forgives
  both copies silently.
- A gene label does **not** identify a locus, and `Hit` carries neither feature identity nor
  genomic coordinate, so name + offset cannot tell the two apart.

### 2b. How much of a region is usable — corrected

The honest single-locus measurement uses a **uniquely-named** feature. `rbcL`, first 200 nt of
1,428:

| window | windows | clean for one intended occurrence |
|---:|---:|---:|
| 9 nt | 192 | **98 (51%)** |
| 14 nt | 187 | **187 (100%)** |
| 19 nt | 182 | **182 (100%)** |

And `psbA`, to show what the duplication does:

| window | clean for **one** intended occurrence | clean when **both copies** are declared intended |
|---:|---:|---:|
| 9 nt | **0 (0%)** | 128 (67%) |
| 14 nt | **0 (0%)** | 187 (100%) |
| 19 nt | **0 (0%)** | 182 (100%) |

Revision 1's 67% is exactly the psbA *declared-both* figure. It was not arbitrary; it was
computing a different quantity from the one it claimed.

**Declaring both psbA copies intended is biologically reasonable** — same gene, and repressing
it would want both. That is precisely the point: the tool must accept an **explicit set** of
intended occurrences, not one, and must record which set was declared. A scope that is not
declared cannot be audited, which is how revision 1 published 67% as a single-site number.

**Half of a region is still real freedom.** 98 usable 9-mers where today the user picks one by
hand, and complete freedom at 14 and 19.

### 2c. What the scanner does and does not establish

The demo default `AAAAUGUGG`: verdict **OCCURS IN A TRANSCRIPT**, 1 transcript occurrence in
**ORF1995**, 5 genomic. Two further numbers, named correctly this time:

| quantity | value |
|---|---:|
| genomic, both strands, **this** target's composition (`expected_by_chance`) | 2.60 |
| genomic, both strands, an **average** 9-mer (`architecture_advice`) | 3.57 |

**Neither is an expected count of transcript occurrences**, and the observation was a transcript
occurrence. Revision 1 set 2.60 against it and concluded "close to background". That comparison
is unsupported, and `offtarget.report` already warns against it in prose: the two were once
printed a few lines apart and a reviewer read the pair as a contradiction.

**And the gap is a classification gap, not detector error.** `scan` does exact substring
matching and is right every time it fires. Under reading A the intended site is in the reference,
so it fires on the intended site too — correctly — and the report cannot say which hit that was.
Calling that a "false-alarm rate", as revision 1 did, miscategorises a missing field as a wrong
answer. The user-visible cost is real and unchanged: the flag cannot be cleared without a manual
lookup, which is the work the tool exists to remove.

---

## 3. What I would build — additive, upstream, and smaller than revision 1 said

**The feature is a selector that sits in front of the existing API and hands it a sequence.**
`design_oneshot`, `design_library` and `parts.select` already accept a chosen sequence. A window
record carries provenance; its `.sequence` goes into today's entry points unchanged. No public
signature is replaced, so nothing that works today can regress.

Revision 1 claimed this "touches the public API of everything". It does not, and that claim was
the main reason it read as expensive.

### 3a. Correct the record first — no new capability

The three claims in §2 that revision 1 got wrong also exist in shipped prose, and fixing them is
independent of everything below. See §7.

### 3b. The occurrence model, before any ranking

This is the load-bearing piece and it must come first, because every later decision reads it.

1. **Key occurrences by feature identity**, never by name. Nine names collide today.
2. **Accept an explicit set of intended occurrences**, possibly several (both psbA copies), and
   record the set with the result.
3. **Classify, never discard.** Retain raw occurrences; add `intended` / `other`. A second
   occurrence *inside* the search interval must stay visible — excluding a whole interval hides
   it. Codex's synthetic control makes this concrete: the same 9-mer at offsets 0 and 9 leaves
   zero occurrences under interval exclusion and one under exact-site exclusion.
4. **`scan()` gains an optional intended-occurrence argument**, defaulting to `None` so today's
   behaviour is byte-identical.

### 3c. Provenance the record must carry

Source ID, source content hash, reference accession where applicable, coordinate convention,
the selected source interval, and the selected window. **A source-relative offset is not a
genomic coordinate** — `psbA` has five joined segments, so `start + offset` is wrong across joins
and wrong again on the minus strand. Do not compute one.

### 3d. Input: a supplied sequence first

Accept an RNA sequence in 5'→3' orientation, plus a length. That is reading B, and it needs no
annotation at all.

An annotation adapter can later offer extracted features **labelled as what they are**: the CDS
entries are coding spans, *not* UTR-inclusive transcripts, so for a native-UTR workflow the
cached annotation **does not contain the sequence the user wants**. Revision 1 missed this
entirely and it is the single biggest correction to the cost of reading A. Automatic native-UTR
retrieval is a separate capability blocked on choosing a data source.

### 3e. Selection: simple and inspectable

Rank on non-intended exact occurrences under a declared reference, plus kit availability. Permit
**"no acceptable window"** as an answer. Show ties; break them deterministically and do not
present the tie-break as biological optimisation.

Two terms revision 1 floated and I now propose to **exclude**: RNA-window GC is not the
coding-DNA synthesis GC that `synthesis_profile` governs, and overhang feasibility belongs to the
downstream design, not to a proxy. `HANDOVER_TO_CODEX.md` §4 is the precedent — a plausible
objective term was identically zero across 42 modules and ranked nothing.

### 3f. Architecture as an output — opt-in only

`ppr.architecture_of` refuses any length outside `{9S: 9, 14S: 14, 19S: 19}`. With a window the
tool can report *"no 9-mer here is unique, so 14 was used"*, which performs the reviewer's
"lengthening is the fix" advice instead of printing it. **Opt-in**: it changes how many modules
a user orders. Note from §2b that at 14 nt every window was clean in both genes measured.

---

## 4. The inventory route — revision 1 was wrong to exclude it

Revision 1 said *"regions do not apply, and pretending they do would put a control on a form that
cannot change the answer."* Measured on `rbcL`'s first 200 nt:

| window | windows | kit-available | **distinct pick lists** |
|---:|---:|---:|---:|
| 9 nt | 192 | 192 | **192** |
| 14 nt | 187 | 187 | **187** |
| 19 nt | 182 | 182 | **182** |

Every window changes the answer:

```
ATGGTTCCACAAACAGAAA -> ('C1', 'A5', 'D5', 'H5', 'E1')
TGGTTCCACAAACAGAAAC -> ('A1', 'H4', 'D5', 'G5', 'E1')
```

A fixed 42-module kit constrains *which modules exist*, not *which target you pick*. Window
selection should serve both routes through the same upstream operation.

One caution Codex is right to attach: availability did **not** discriminate here — 192 of 192
available. It must not be sold as a useful ranking term on this evidence. It is a filter that
happened not to bind.

---

## 5. Naming

`HANDOVER_TO_CODEX.md` item 10 records a collision where `profile` meant three things.
**`region` is already taken**: `ordering.ProductProfile.region` is a vendor shipping region
(`"EU"`). Proposed: **`target_window`** for the selected span, **`search_sequence`** for the
input, **never `region`**. Source origin must be explicit and must **not** default silently to
synthetic.

---

## 6. The Colab constraint — and the part of it I overstated

Both notebooks' Setup cell runs `pip install git+https://.../clippr.git`, resolving to the
**default branch**. Anything a cell imports must exist on `main` when a user opens Colab.
`validation/run_notebook.py` cannot catch a violation: it puts the working tree's `src/` first on
`sys.path`. That combination broke every Colab session once — `HANDOVER_TO_CODEX.md` §7 item 9.

Revision 1 called `CLIPPR_inventory.ipynb` cell 2 "the working pattern" for graceful
degradation. **It is a detection pattern, not a degradation pattern**: it prints a banner and
does not gate anything downstream. And `validation/notebook_on_published_build.py` counts cells
that complete without raising — it does not assert the requested output was produced. So "every
cell executed, 0 failed" is weaker evidence than revision 1 leaned on it for.

What is actually required:

- Keep the published-build check, and **additionally assert the selected operation's output
  exists**, recording the package revision.
- **Gate** dependent cells on the capability, not just banner it.
- If a build lacks window selection, **disable that mode with an explanation** — never silently
  reinterpret a supplied search sequence as an exact target, and never report it as processed.
- The package change reaches `main` before the notebook change, or in the same push.

---

## 7. Notebook prose that is wrong today, and needs no decision

These are corrections, not features. Three of them repair claims that revision 1 either made or
proposed to make.

1. **Remove the overstatements.** The off-target cell's markdown says an occurrence means "the
   protein binds there too" — the module's own docstring says occurrence is necessary and never
   sufficient. Revision 1 proposed replacing it with a *second* overstatement ("a hit here is a
   real off-target"). Both exceed exact substring matching. Codex's wording is close to right:

   > Exact matches in the scanned reference are shown below. If you supplied an intended
   > occurrence, it is listed separately. Other matches are candidates for review; this result
   > does not establish binding. The reference and annotation coverage are recorded with the
   > result.

2. **State the reference's scope.** 109 features, 9 duplicated names, CDS spans rather than
   UTR-inclusive transcripts. A user reading "not found in the host" should know what was
   searched. Also: for a supplied construct, scan the supplied library sequences too, or say
   plainly that library crosstalk was not screened — "no host match" is not a statement about the
   experiment's own sequences.
3. **Keep the flagged demo, labelled.** `AAAAUGUGG` flags on first run. Keep it as a deliberate
   occurrence example rather than swapping in a clean one, and add a clean example separately if
   wanted. Do **not** reinstate revision 1's "close to background" gloss (§2c).
4. **Move the library cell's three hardcoded targets.** Its markdown has to open by explaining
   they are "not derived from your target above", which says the cell is misplaced rather than
   badly worded.

Still held, because it depends on what Step 1 becomes: grouping and collapsing the 22-cell
scroll, 8 of which are `Optional`.

---

## 8. Notebook before web app — conclusion unchanged, evidence replaced

Revision 1 cited **4.70 s / 6.26 s** as the designer route. Those are `chained_routes` in
`work/runtime/runtime_profile.json` — the **inventory** recode → optimise → search → export
chain. The one-shot designer figures are a different section, `oneshot_synthesis_route`:

| architecture | saved one-shot duration |
|---|---:|
| 9S | 0.5221 s |
| 14S | 14.1180 s |
| 19S | **204.9769 s** |

Each **n = 1**, with off-target screening disabled because it fetches a genome. Read today, not
re-run. Neither bounds nor predicts arbitrary candidates.

**205 seconds for a 19S design makes the infrastructure point harder, not softer** — that cannot
be a synchronous web request under any design. But revision 1's conclusion rested on the wrong
numbers, and the honest version is: a web interface needs an async job model, and that is work
that buys nothing until the selection layer exists.

Separately, **screening is cheap and design is not.** Codex measured enumeration + occurrence
counting + kit selection across one length of a 200 nt input at a median **0.065–0.071 s**
(3 repeats per length, excluding imports, cache load, genomic scan, compilation, fidelity, de
novo design and export). So the runtime policy is forced: **screen every window, design one
selected candidate** — or an explicitly bounded shortlist. Never run full optimisation per
window. Record candidate count, shortlist budget, cache state and realised elapsed time.

Notebook-first stays the recommendation, but as a product choice about reaching the reviewer
quickly — not as something a 4-second runtime forced. A web interface can reuse the same
selection layer later.

---

## 9. What I need decided

Revision 1 asked for an A/B/C choice. That was a false choice (§1); these are what remain.

1. **Is the first increment reading B — "choose a window inside a sequence I supply"?** It needs
   no annotation, no intended-occurrence model, and no new data source. I recommend yes.
2. **Should §7's prose corrections ship now, separately?** They fix claims that are wrong in the
   repository today. I recommend yes, and would send them first.
3. **Is a native-gene UTR workflow wanted at all** (reading A), given that the cached annotation
   holds **coding spans, not UTRs**, so it would need a new data source? This is the question
   revision 1 should have asked instead of the A/B/C one.
4. **May architecture (9S/14S/19S) become an opt-in output?** It changes how many modules get
   ordered, so it is a wet-lab call.
5. **Naming** (§5): `target_window` / `search_sequence`, never `region`.

Defer: reading C until someone wants generated targets; native-UTR retrieval until its source is
named. Neither blocks increment B.

---

## 10. What this establishes, and what it does not

- **Not better binding, expression, or experimental success.** Nothing here predicts affinity. No
  sequence from this project has been synthesised.
- **Not a false-positive rate for the scanner.** It reports exact occurrences correctly. 51% and
  67% describe how many candidate windows carry additional occurrences under a *declared* scope,
  on one gene, one 200 nt interval, one genome.
- **Not single-locus specificity from length.** 14 and 19 nt were clean in both genes measured.
  That is two genes, and the inverted repeat shows how fast a plausible predicate goes wrong.
- **Not that kit availability is a useful ranking term.** 192 of 192 available — it did not
  discriminate.
- **Not that the 0.065 s screening figure predicts the feature's runtime.** It excludes most of
  the work.
- **Not that any region contains a usable window.** 0% of psbA's 9-mers qualify for a single
  declared site. "No acceptable window" must be an answer the tool can give.

---

## 11. Reproducing §2

Save the block below as `proposal_evidence.py` in the repository root. PowerShell 5.1:

```powershell
$env:PYTHONPATH = "src"; .venv\Scripts\python.exe .\proposal_evidence.py
```

Cached inputs; no network. Verified output:

```
203828 bp; 109 annotated features; 100 distinct names
duplicated names: ['psbA', 'rps2', 'rrn3', 'rrn5', 'rrn7', 'rrnL', 'rrnS', 'trnA', 'trnI']

AAAAUGUGG: OCCURS IN A TRANSCRIPT; transcript 1 ['ORF1995']; genomic 5
  genomic, both strands, this target's composition : 2.60
  genomic, both strands, an AVERAGE 9-mer          : 3.57
  neither is an expected number of TRANSCRIPT hits.

rbcL (unique name), first 200 nt of a 1428 nt feature
  9-mer: 192 windows | 98 clean for ONE intended site (51%) | 98 clean when every copy of rbcL in the window is declared intended (51%)
  14-mer: 187 windows | 187 clean for ONE intended site (100%) | 187 clean when every copy of rbcL in the window is declared intended (100%)
  19-mer: 182 windows | 182 clean for ONE intended site (100%) | 182 clean when every copy of rbcL in the window is declared intended (100%)

psbA (duplicated), first 200 nt of a 1059 nt feature
  9-mer: 192 windows | 0 clean for ONE intended site (0%) | 128 clean when every copy of psbA in the window is declared intended (67%)
  14-mer: 187 windows | 0 clean for ONE intended site (0%) | 187 clean when every copy of psbA in the window is declared intended (100%)
  19-mer: 182 windows | 0 clean for ONE intended site (0%) | 182 clean when every copy of psbA in the window is declared intended (100%)
```

```python
"""Corrected evidence. Occurrences are keyed by FEATURE INDEX, never by gene name."""
from collections import Counter
from clippr.offtarget import (load_genome, load_transcripts, scan,
                              find_in_transcripts, architecture_advice)

genome, trs = load_genome(), load_transcripts()
dup = {n: c for n, c in Counter(t.name for t in trs).items() if c > 1}
print(f"{len(genome)} bp; {len(trs)} annotated features; "
      f"{len(set(t.name for t in trs))} distinct names")
print(f"duplicated names: {sorted(dup)}")

def occurrences(target):
    """(feature index, offset) -- the identity a name-keyed hit throws away."""
    return [(j, h.position) for j, tr in enumerate(trs)
            for h in find_in_transcripts(target, [tr])]

# The demo default, with both quantities named correctly.
r = scan("AAAAUGUGG", genome, trs)
print()
print(f"AAAAUGUGG: {r['verdict']}; transcript {r['n_transcript']} {r['genes']}; "
      f"genomic {r['n_genomic']}")
print(f"  genomic, both strands, this target's composition : "
      f"{r['expected_by_chance']:.2f}")
print(f"  genomic, both strands, an AVERAGE 9-mer          : "
      f"{architecture_advice(genome)[9]:.2f}")
print("  neither is an expected number of TRANSCRIPT hits.")

for name in ("rbcL", "psbA"):
    pick = next(t for t in trs if t.name == name)
    pi = trs.index(pick)
    W = 200
    region = pick.sequence[:W]
    tag = "duplicated" if name in dup else "unique name"
    print()
    print(f"{name} ({tag}), first {W} nt of a {len(pick.sequence)} nt feature")
    for k in (9, 14, 19):
        n = one = declared = 0
        for i in range(len(region) - k + 1):
            occ = occurrences(region[i:i + k])
            n += 1
            one += not [o for o in occ if o != (pi, i)]
            declared += not [o for o in occ
                             if trs[o[0]].name != name or o[1] >= W - k + 1]
        print(f"  {k}-mer: {n} windows | {one} clean for ONE intended site "
              f"({100.0*one/n:.0f}%) | {declared} clean when every copy of "
              f"{name} in the window is declared intended ({100.0*declared/n:.0f}%)")
```

The kit figures in §4 come from `clippr.parts.select` over the same windows, counting
`PartsPlan.available` and distinct `PartsPlan.plates`. The runtime figures in §8 are read from
`work/runtime/runtime_profile.json` — `oneshot_synthesis_route`, not `chained_routes`.

Not committed, deliberately: nothing should land in `validation/` before the direction is
approved. If this proceeds it becomes `validation/experiments/region_targeting.py`, writes a JSON
artefact under `work/`, and gets a line in `docs/verification_commands.md`, which is this
repository's convention for every other measurement.

---

## 12. Acceptance tests this feature must pass

From Codex's review, and I would hold the work to them:

- duplicate names at distinct loci; repeated sites inside one interval
- zero, one and many acceptable windows
- coordinate convention and normalisation, including compound locations and minus strand
- explicit multi-length selection
- selected-window round trip into **both** the designer and inventory routes
- an exact-target run matching today's result under identical settings
- installed-build notebook checks **plus output assertions**
- the complete feature timed on a short and a long case, cold and warm cache
