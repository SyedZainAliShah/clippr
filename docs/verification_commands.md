# Verification commands

One copyable invocation per operation, with its required inputs and where its output lands.
Every validation program takes explicit `--` paths so an independent reviewer can write
somewhere else and never overwrite the executor's evidence.

Run from the repository root with `.venv\Scripts\python.exe`. Set `PYTHONUTF8=1` on Windows.

**Settings are part of the result.** Several figures reproduce only at a specific budget — the
defaults are smaller. Where that matters it is stated with the command, not left to be
discovered.

---

## Inventory

**Import and validate an inventory**

```bash
python -c "import sys;sys.path.insert(0,'src');from clippr import inventories as inv;lib=inv.load_deposited('data/grasp_supp/Table S1.xlsx');print(len(lib),lib.version);print(inv.validate(lib) or 'valid');print(inv.synthesis_profile_notes(lib))"
```

Needs Supplementary Table S1 (not shipped — see `NOTICE.md`). Prints module count, content
version, structural problems, and any module outside the synthesis profile.

**Recode all 42 modules with fixed interfaces**

```bash
python -c "import sys,json;sys.path.insert(0,'src');from clippr import workflow as w;t=json.load(open('data/codon_tables/kazusa_3055.json'));r=w.recode_for_host('data/grasp_supp/Table S1.xlsx',t,'out/recode',label='chlamy',seeds=8,wall_seconds=600);print(r.summary);print(r.artefacts)"
```

**`seeds=8` is required to reproduce the saved inventory.** Four attempts give a different
inventory (mean CAI 0.618830 rather than 0.621004) and differ from the saved one at ten
modules. Repeated four-attempt runs agree with each other; they simply are not the eight-
attempt result.

## Optimisation

**Optimise a library as a collection**

```bash
python -c "import sys,json;sys.path.insert(0,'src');from clippr import workflow as w;t=json.load(open('data/codon_tables/kazusa_3055.json'));r=w.optimise_collection('work/phaseb/inventory_recoded.json',t,'out/collection',mode='greedy',max_proposals=168,wall_seconds=600);print(r.summary)"
```

**`max_proposals=168` is required to reproduce the reported table.** Two different defaults
exist and must not be conflated: `optimise_library` defaults to **400**, while
`phasec_gate.py` defaults to **84**. Neither reproduces the table; pass 168 explicitly.

**The matched three-mode comparison, 18 runs**

```bash
python validation/experiments/phasec_gate.py --max-proposals 168 --out work/phasec
```

Record the seed, processing order and budget with any result taken from this. All three modes
draw proposals from a seeded generator; none is seed-free deterministic.

## Junction search

**Search junction/codon alternatives and export the observed front**

```bash
python -c "import sys,json;sys.path.insert(0,'src');from clippr import workflow as w;t=json.load(open('data/codon_tables/kazusa_3055.json'));r=w.explore_interfaces('work/phaseb/inventory_recoded.json',['AAAAUGUGG','UUACACGUGCGUAC','CUAUCACAUCACAUAAGCG'],t,'out/front',max_evaluations=100,wall_seconds=600);print(r.summary);print(r.data['recommended'])"
```

Searches **level-0 junction classes only**. Block joins belong to a level-1 reaction whose
block-plasmid context this package does not hold; changing them cannot affect any reaction
that can be scored.

**Commit to a candidate and rebuild from it**

```bash
python -c "import sys;sys.path.insert(0,'src');from clippr import workflow as w;s=w.select_interface('work/phaseb/inventory_recoded.json','out/front/interface_front.json','out/selected',choice='recommended');print(s.summary);c=w.load_and_compile(s.artefacts['inventory'],['AAAAUGUGG'],'out/recompiled');print(c.summary);i=w.order_items_for(s.artefacts['inventory'],'out/items');print(i.summary)"
```

`choice='incumbent'` is a valid path. Selecting anything else changes the interface version
and requires recompilation, which the result states.

## Compilation and products

**Compile targets from an inventory version**

```bash
python -c "import sys;sys.path.insert(0,'src');from clippr import workflow as w;r=w.load_and_compile('work/phaseb/inventory_recoded.json',['AAAAUGUGG','UUACACGUGCGUAC','CUAUCACAUCACAUAAGCG'],'out/compiled');print(r.summary);print(r.data['inventory_version'])"
```

**Validate complete exported products and annotations**

```bash
python validation/experiments/phaseb_gate.py --recoded work/phaseb/inventory_recoded.json --out work/phaseb
python validation/independent/v7_export_semantics.py
```

## Ordering

**Evaluate oPools eligibility, price and pool plans**

```bash
python -c "import sys,json;sys.path.insert(0,'src');from pathlib import Path;from clippr import workflow as w;from clippr.ordering import CURRENT_OPOOL_50PMOL;i=w.order_items_for('work/phaseb/inventory_recoded.json','out/items');items=json.loads(Path(i.artefacts['order_items']).read_text());r=w.plan_order(items,CURRENT_OPOOL_50PMOL,'out/order');print(r.summary);print(r.data['cost'])"
```

`CURRENT_OPOOL_50PMOL` carries the published rules (40–350 nt, 2–384 oligos per pool at
50 pmol) with **no prices**, so cost reports unavailable. `HISTORICAL_OPOOL` is kept for the
older figures and is labelled historical. Supply a dated profile through `load_profile` for an
estimate. No order is placed by any command here.

```bash
python validation/experiments/phasee_gate.py --out work/phasee
```

## Benchmark against the reference

**Run the pinned reference in its own environment**

```bash
python validation/experiments/m3_reference_comparison.py --repeats 5
```

Archives the previous comparison rather than deleting it, and writes each run separately under
`work/m3/reference/runs/`.

## Release

**Run release verification with exact artefact identities**

```bash
python validation/release_check.py --clean-install
```

Ten checks. A check the script cannot evaluate is reported as needing a human, never as passed.
`--clean-install` builds a wheel, installs it into a throwaway environment with `PYTHONPATH`
dropped, and requires the README example to reproduce its documented output verbatim.

**The full test suite, counted from collection**

```bash
python -m pytest --collect-only -q
python -m pytest
```

Count from the collection total and the summary line, never from progress dots.

---

## Gates

| Gate | Command | Output |
|---|---|---|
| B | `python validation/experiments/phaseb_gate.py` | `work/phaseb/gate_b.json` |
| C | `python validation/experiments/phasec_gate.py --max-proposals 168` | `work/phasec/gate_c.json` |
| D | `python validation/experiments/phased_gate.py` | `work/phased/gate_d.json` |
| E | `python validation/experiments/phasee_gate.py` | `work/phasee/gate_e.json` |
| F | `python validation/experiments/phasef_gate.py` | `work/phasef/gate_f.json` |
| oracle | `python validation/experiments/matched_oracle_benchmark.py` | `work/oracle/matched_benchmark.json` |
| level 1 | `python validation/experiments/level1_geometry.py` | `work/level1/geometry.json` |

Each accepts `--out` so a reviewer can write to their own directory.

The oracle benchmark reads the preserved reference output at
`work/m3/reference/output/optimized_library.csv` and the shipped inventory at
`work/phaseb/inventory_recoded.json`; it runs nothing in the reference's environment, so it
reproduces without that checkout being installed. Re-running the reference itself is
`validation/experiments/m3_reference_comparison.py`, which does need it.

## Supplied inputs

| Input | Needed for | If absent |
|---|---|---|
| Supplementary Table S1 | the deposited inventory | not shipped; `load_deposited` raises with a named reason |
| a saved inventory | every step after recoding | none — it carries its own sequences |
| codon table | recoding, optimisation | fetched once from Kazusa and cached, or supplied |
| chloroplast genome | host screening | fetched from NCBI on first use, or screening disabled |
| block-plasmid context | the level-1 block **geometry** | **not needed** — the released block fragment is the joined module inserts, so the ends come from the compiled product (`validation/experiments/level1_geometry.py`) |
| the level-1 reaction's **participant list** | scoring that reaction | **not available**; the deposited BsaI reaction co-assembles parts this package does not compile, so the reaction is reported unscorable and only a labelled PPR-only diagnostic is given |
