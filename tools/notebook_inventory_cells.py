"""The reusable-inventory section of the Colab notebook.

Kept apart from `build_notebook.py` because it is a self-contained route with its own inputs
and its own caveats, and because that file was already long enough that appending a second
workflow to it would bury both.

`cells(md, code)` takes the builder's own cell constructors, so formatting and form metadata
stay consistent with the rest of the notebook rather than being reimplemented here.
"""


def cells(md, code):
    """Every cell of the inventory route, in order."""
    out = []

    out.append(md("""
---

## Step 1 — load the deposited kit
"""))

    # Step 1 used to say "or use the upload cell near the top of this notebook" when the
    # file was missing. There was no such cell, so the error named a remedy that did not
    # exist -- worse than naming none. This is that cell.
    #
    # It prompts only when the file is absent. `files.upload()` blocks on a browser file
    # picker, which would hang the notebook gate in `notebook_on_published_build.py`; that
    # harness stages Table S1 itself, so the file is already there and this cell says so and
    # moves on.
    out.append(code('''
#@markdown Run this and choose **Supplementary Table S1** (`.xlsx`) from your computer. It
#@markdown holds the module insert sequences and is not distributed with this package --
#@markdown see `NOTICE.md`.
#@markdown
#@markdown The folder icon in the sidebar does the same job. Either way the file lives on
#@markdown this runtime only: Colab deletes it when the runtime recycles, so it is needed
#@markdown once per session, not once ever.
#@markdown
#@markdown Already have it here? This cell says so and does nothing.
from pathlib import Path as _Path

_WANTED = "Table S1.xlsx"

if _Path(_WANTED).is_file():
    print(f"{_WANTED} is already here -- nothing to upload")
else:
    try:
        from google.colab import files as _files
    except ImportError:
        print(f"not running in Colab. Put {_WANTED} next to this notebook, or set the")
        print("table_s1 field in Step 1 to wherever it already is.")
    else:
        for _name, _data in _files.upload().items():
            print(f"uploaded {_name}  ({len(_data):,} bytes)")
            if _name != _WANTED:
                print(f"   note: Step 1 looks for {_WANTED!r}. Either rename this file or")
                print(f"   set Step 1's table_s1 field to {_name!r}.")
''', title="Upload Supplementary Table S1"))

    out.append(code('''
#@markdown Upload **Supplementary Table S1** (`.xlsx`), then put its filename here.
table_s1 = "Table S1.xlsx"  #@param {type:"string"}

import json as _json
from pathlib import Path

from ppropt import inventories as inv
from ppropt import workflow as w

# Resolves the codon table from the host form above. It cannot borrow the name `table`:
# that name is used elsewhere for a dataframe
# to a display DataFrame, so every call here was handed a DataFrame and raised. The name is
# distinct so a later display cell cannot shadow it again.
from ppropt import constants as _C
from ppropt.design import _resolve_codon_table as _resolve


def done(result):
    """The result if the task succeeded, else None, with the reasons printed.

    For the two tasks that refuse *without writing an artefact*: `explore_interfaces` (an
    empty front has nothing to select) and `select_interface` (it returns before publishing
    anything the verifier rejected). Reaching into `.artefacts[...]` after one of those
    turned a designed refusal into `KeyError: 'inventory'`, which reads like a broken
    notebook rather than the explanation it actually is.

    Not for `recode_for_host`, `load_and_compile` or `order_items_for`. Those write their
    artefact either way and report `ok=False` for a *partial* result -- one infeasible
    module, one unspellable target, one substrate outside the synthesis contract. Nulling
    those would throw away 41 good modules because of the 42nd. They guard on the artefact
    being present instead, and print what fell short.
    """
    if result is None or result.ok:
        return result
    for failure in result.failures:
        print(f"   {failure.get('stage', '?')}: {failure.get('reason', failure)}")
    print()
    print("Nothing was written. Fix the above and re-run this cell; later steps will stay")
    print("skipped until it succeeds.")
    return None

_host_code = genetic_code_override or _C.ORGANISMS[organism][1]
host_codon_table = _resolve(codon_table_file or None, organism, _host_code)
print(f"codon table for {organism}, genetic code {_host_code}: "
      f"{len(host_codon_table)} amino acids")

deposited = None
if not Path(table_s1).is_file():
    print(f"{table_s1} not found. Run the upload cell just above this one, or use the")
    print("folder icon in the sidebar. Nothing below this step will run without it.")
else:
    deposited = inv.load_deposited(table_s1)
    print(f"loaded {len(deposited)} modules  |  version {deposited.version}")
    print(f"structural problems: {inv.validate(deposited) or 'none'}")

    notes = inv.synthesis_profile_notes(deposited)
    if notes:
        print()
        print("Modules outside PPROPT's synthesis profile, before you order anything:")
        for module_id, why in notes.items():
            print(f"   {module_id}: {'; '.join(why)}")
        print()
        print("That limit is PPROPT's engineering choice, not a published rule, so this")
        print("is a statement about fit -- not a defect in the deposited kit.")
''', title="Step 1 — load the deposited kit"))

    out.append(md("""
## Recode the kit for your host

Every module is rewritten to use your host's preferred codons while **the protein, the length
and both four-base interfaces stay exactly as they were**. A recoded module therefore drops
straight into an assembly built from unrecoded neighbours.

A module that cannot be improved keeps its original sequence and is reported as *unchanged*.
That is a delivered module, not a gap.
"""))

    out.append(code('''
label = "my-host"  #@param {type:"string"}
candidate_seeds = 4  #@param {type:"slider", min:1, max:8, step:1}
recode_seconds = 600  #@param {type:"integer"}

recoded = None
if deposited is not None:
    recoded = w.recode_for_host(table_s1, host_codon_table, "out/inventory", label=label,
                                seeds=candidate_seeds, wall_seconds=recode_seconds)
    print(recoded.summary)
    for failure in recoded.failures:
        print(f"   {failure['module']}: {failure['reason']}")

if recoded is not None and "inventory" in recoded.artefacts:
    print(f"saved to {recoded.artefacts['inventory']}")

    report = _json.loads(Path(recoded.artefacts["report"]).read_text(encoding="utf-8"))
    scored = [o for o in report["objectives"].values() if "cai_after" in o]
    if scored:
        before = sum(o["cai_before"] for o in scored) / len(scored)
        after = sum(o["cai_after"] for o in scored) / len(scored)
        print()
        print(f"mean CAI {before:.4f} -> {after:.4f} over {len(scored)} modules")
        print()
        print("CAI predicts how well codons match the host. It is not a measurement of")
        print("expression, and a large CAI gain is not a proportional expression gain.")
''', title="Step 2 — recode for your host, interfaces frozen"))

    out.append(md("""
## Compile targets from the inventory

Give it targets; it returns the module list, the assembled product and every junction. A
target the kit cannot spell is listed with the reason why, and the others still run.

The product is **joined module inserts only** — no acceptor backbone, so it is not an
expression construct.
"""))

    out.append(code('''
inventory_targets = "AAAAUGUGG, UUACACGUGCGUAC"  #@param {type:"string"}
fusion_site = "AATG"  #@param ["AATG", "AGGT"]

wanted = [t.strip().upper().replace("T", "U")
          for t in inventory_targets.split(",") if t.strip()]
compiled = None
if recoded is not None:
    compiled = w.load_and_compile(recoded.artefacts["inventory"], wanted,
                                  "out/compiled", fusion_site=fusion_site)
    print(compiled.summary)
    for failure in compiled.failures:
        print(f"   could not build {failure['target']}: {failure['reason']}")

# Not `compiled.ok`: one unspellable target makes this False while the others still
# compiled and the artefact still exists. The artefact's presence is the real condition.
if compiled is not None and "compiled_targets" in compiled.artefacts:
    products = _json.loads(
        Path(compiled.artefacts["compiled_targets"]).read_text(encoding="utf-8"))
    for name, product in products.items():
        print(f"   {name}: {product['product_nt']} nt, "
              f"{len(product['modules'])} modules, {product['reactions']} reaction(s)")
''', title="Step 3 — compile targets"))

    out.append(md("""
## Optimise the inventory as a collection

Recoding treats each module alone, so it cannot see what the modules share *with each other*.
This step can trade a little codon adaptation in one module for less shared sequence across
the set.

`greedy` is the default: it accepts only improvements, and over 18 matched runs it reaches
within 2% of what `anneal` achieves on shared sequence while making a third as many changes.
`anneal` explores further and gives up somewhat more codon adaptation for a small further
gain. `random` is a control — a method that cannot beat it has not been shown to work.

All three are **reproducible for a fixed seed and settings**, not seed-free deterministic:
the candidate replacements come from a seeded generator in every mode.
"""))

    out.append(code('''
mode = "greedy"  #@param ["greedy", "anneal", "random"]
max_proposals = 168  #@param {type:"integer"}

optimised = None
if recoded is not None:
    optimised = w.optimise_collection(recoded.artefacts["inventory"], host_codon_table,
                                      "out/collection", mode=mode,
                                      max_proposals=max_proposals, wall_seconds=600)
    print(optimised.summary)

if optimised is not None and "report" in optimised.artefacts:
    search = _json.loads(Path(optimised.artefacts["report"]).read_text(encoding="utf-8"))
    print(f"   adaptation {search['incumbent']['adaptation']} -> "
          f"{search['final']['adaptation']}")
    print(f"   shared k-mers per pair {search['incumbent']['sharing_per_pair']} -> "
          f"{search['final']['sharing_per_pair']}")
    print()
    print("The sharing term steers the search. It is not a recombination probability,")
    print("and lowering it does not necessarily shorten the worst shared tract.")
''', title="Step 4 — optimise the collection"))

    out.append(md("""
## Explore junction trade-offs

The kit's internal junction overhangs can be swapped for other four-base sequences encoding
the same protein. Different choices give different predicted assembly fidelity, so there is a
real trade-off to look at rather than a single answer.

⚠️ **Choosing an alternative changes the interface version.** Every target must be recompiled
from the selected inventory, and modules from the previous version cannot be mixed with it.

The result is the **observed** front — nondominated among the candidates actually evaluated.
A bounded search does not enumerate the space, so this is not the globally optimal front.
"""))

    out.append(code('''
max_evaluations = 100  #@param {type:"integer"}

front = None
if recoded is not None:
    front = w.explore_interfaces(recoded.artefacts["inventory"], wanted, host_codon_table,
                                 "out/interfaces", max_evaluations=max_evaluations,
                                 wall_seconds=600)
    print(front.summary)
    front = done(front)

if front is not None:
    print()
    print(f"incumbent  : {front.data['incumbent']}")
    print(f"recommended: {front.data['recommended']}")
    print()
    print(f"policy: {front.data['recommendation_policy']}")
    if front.data["fidelity_degenerate"]:
        print()
        print("Every candidate scored the same fidelity -- no trade-off to make here.")
    print()
    print(front.data["selecting_an_alternative"])
''', title="Step 5 — explore junction trade-offs"))

    out.append(md("""
## Commit to a choice

Exploring a front does nothing until you pick one. This saves the chosen candidate as a real
inventory and builds the order items from *it*, so what you order is what you chose.

The front is **bound** to the inventory it was measured against. Selecting it against a
different one is refused rather than silently reporting the wrong inventory's numbers.

`incumbent` is a perfectly good choice when the front offers nothing worth an interface change.
"""))

    out.append(code('''
choice = "recommended"  #@param ["recommended", "incumbent"]

selected = recompiled = order_items = None
if front is not None:
    selected = w.select_interface(recoded.artefacts["inventory"],
                                  front.artefacts["front"], "out/selected",
                                  choice=choice, codon_table=host_codon_table)
    print(selected.summary)
    selected = done(selected)

if selected is not None:
    recompiled = w.load_and_compile(selected.artefacts["inventory"], wanted,
                                    "out/recompiled")
    print(recompiled.summary)
    for failure in recompiled.failures:
        print(f"   could not build {failure['target']}: {failure['reason']}")

    # ok=False here means some substrate is outside the synthesis contract, not that the
    # order failed. All 42 are still written; the flagged ones are named so you can decide.
    order_items = w.order_items_for(selected.artefacts["inventory"], "out/order_items")
    print(order_items.summary)
    for failure in order_items.failures:
        print(f"   {failure['module']}: {failure['reason']}")
''', title="Step 6 — commit to a choice and rebuild"))

    out.append(md("""
## Plan the order

Eligibility, pooling and export against a written product profile.

The order items are **assembly-ready substrates**, not bare inserts: each is wrapped so BbsI
releases the level-0 fragment with the correct exposed ends, and each is digest-verified
before export. An A-module insert starts with its own `AATG` fusion site, while the substrate
must expose `CTCA` — ordering the insert would order something that cannot assemble.

**A local check is not vendor approval.** Eligibility and pooling are checked against a
written-down reading of the vendor's published rules; the vendor's current documentation is
always the authority, and no order is placed here.
"""))

    out.append(code('''
#@markdown Leave blank for the current published oPools rules. Supply your own profile
#@markdown JSON to check against a different product.
profile_file = ""  #@param {type:"string"}

from ppropt.ordering import CURRENT_OPOOL_50PMOL, load_profile

profile = load_profile(profile_file) if profile_file.strip() else CURRENT_OPOOL_50PMOL

# The order comes from the inventory you selected, never from whatever happens to be left in
# notebook memory. Taking sequences from the synthesis route's `lib` here meant a user could
# explore one inventory and order something unrelated to it.
order = None
if order_items is None or "order_items" not in order_items.artefacts:
    print("Run the cells above first -- the order is built from the selected inventory.")
else:
    items = _json.loads(Path(order_items.artefacts["order_items"]).read_text(
        encoding="utf-8"))
    print(f"ordering {len(items)} assembly-ready substrates from "
          f"{selected.data['to_version']}")
    print(f"exposed ends: {order_items.data['exposed_ends']}")
    order = w.plan_order(items, profile, "out/order")
    print(order.summary)
    for failure in order.failures:
        print(f"   {failure['stage']}: {failure['reason']}")

    print()
    for rule in order.data["eligibility"]["unresolved"]:
        print(f"   unresolved: {rule}")
''', title="Step 7 — eligibility and pool plan"))

    # Where to put the pipette. The order names 42 modules, and a list of plasmid
    # identifiers gets transcribed by hand into a picking order -- which is the step that
    # goes wrong.
    out.append(code('''
from ppropt.viz import kit_plate_svg

if order_items is None:
    print("Run the cells above first.")
else:
    import importlib.resources as _res
    _modules = _json.loads(_res.files("ppropt").joinpath("parts.json")
                           .read_text(encoding="utf-8"))["modules"]
    _wanted = [i["module"] for i in items]
    display(HTML(kit_plate_svg(
        _wanted, _modules,
        title=f"GRASP kit — {len(_wanted)} wells to pick",
        note=f"{len(_wanted)} modules for this build, from inventory "
             f"{selected.data['to_version']}")))
''', title="Step 7b — where to put the pipette"))

    out.append(md("""
## Take the package

One directory holding the manifest, every task's result, the artefacts and — importantly —
every failure. Artefact paths are recorded relative to the package, so it can be moved or
shared without carrying a machine-specific path along with it.
"""))

    out.append(code('''
done = [r for r in (recoded, compiled, optimised, front, selected, recompiled,
                    order_items, order) if r is not None]
if not done:
    print("Nothing to package yet -- run the cells above.")
else:
    package = w.write_package(done, "out/package",
                              inputs={"table_s1": table_s1, "targets": wanted})
    print(f"wrote {package}")
    summary = _json.loads(Path(package).read_text(encoding="utf-8"))
    print(f"   {len(summary['tasks'])} tasks, ok={summary['ok']}, "
          f"{len(summary['failures'])} failure(s)")
    print()
    print(summary["scope"])

    try:
        import shutil

        from google.colab import files
        shutil.make_archive("ppropt_package", "zip", "out")
        files.download("ppropt_package.zip")
    except ImportError:
        print()
        print("Not running in Colab -- the package is in out/.")
''', title="Step 8 — download the result package"))

    return out
