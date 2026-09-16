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

# The reusable-inventory route

Everything above designs a **new coding sequence per target**. This section is the other
route: take the deposited GRASP module kit, recode it once for your host, then *compile*
targets from it — the DNA is ordered once and reused.

The two routes answer different questions. Synthesis gives you any target the PPR code can
express. The kit gives you the targets its modules can spell, far more cheaply, because you
are assembling parts you already have.

**This route needs Supplementary Table S1.** It holds the module insert sequences and is not
distributed with this package — see `NOTICE.md`. Upload it below. Every step after the recode
works from a *saved inventory*, which carries its own sequences, so you need Table S1 once.
"""))

    out.append(code('''
#@markdown Upload **Supplementary Table S1** (`.xlsx`), then put its filename here.
table_s1 = "Table S1.xlsx"  #@param {type:"string"}

import json as _json
from pathlib import Path

from clippr import inventories as inv
from clippr import workflow as w

# This route resolves its own codon table, from the same organism and file chosen in the
# design form above. It cannot borrow the name `table`: the fragment-table cell rebinds that
# to a display DataFrame, so every call here was handed a DataFrame and raised. The name is
# distinct so a later display cell cannot shadow it again.
from clippr import constants as _C
from clippr.design import _resolve_codon_table as _resolve

_host_code = genetic_code_override or _C.ORGANISMS[organism][1]
host_codon_table = _resolve(codon_table_file or None, organism, _host_code)
print(f"codon table for {organism}, genetic code {_host_code}: "
      f"{len(host_codon_table)} amino acids")

deposited = None
if not Path(table_s1).is_file():
    print(f"{table_s1} not found. Upload it with the folder icon in the sidebar,")
    print("or use the upload cell near the top of this notebook.")
else:
    deposited = inv.load_deposited(table_s1)
    print(f"loaded {len(deposited)} modules  |  version {deposited.version}")
    print(f"structural problems: {inv.validate(deposited) or 'none'}")

    notes = inv.synthesis_profile_notes(deposited)
    if notes:
        print()
        print("Modules outside CLIPPR's synthesis profile, before you order anything:")
        for module_id, why in notes.items():
            print(f"   {module_id}: {'; '.join(why)}")
        print()
        print("That limit is CLIPPR's engineering choice, not a published rule, so this")
        print("is a statement about fit -- not a defect in the deposited kit.")
''', title="Inventory 1 - load the deposited kit"))

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
''', title="Inventory 2 - recode for your host, interfaces frozen"))

    out.append(md("""
## Compile targets from the inventory

Give it targets; it returns the module list, the assembled product and every junction. A
target the kit cannot spell is reported with its reason rather than raising.

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

    products = _json.loads(
        Path(compiled.artefacts["compiled_targets"]).read_text(encoding="utf-8"))
    for name, product in products.items():
        print(f"   {name}: {product['product_nt']} nt, "
              f"{len(product['modules'])} modules, {product['reactions']} reaction(s)")
''', title="Inventory 3 - compile targets"))

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
    search = _json.loads(Path(optimised.artefacts["report"]).read_text(encoding="utf-8"))
    print(f"   adaptation {search['incumbent']['adaptation']} -> "
          f"{search['final']['adaptation']}")
    print(f"   shared k-mers per pair {search['incumbent']['sharing_per_pair']} -> "
          f"{search['final']['sharing_per_pair']}")
    print()
    print("The sharing term steers the search. It is not a recombination probability,")
    print("and lowering it does not necessarily shorten the worst shared tract.")
''', title="Inventory 4 - optimise the collection"))

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
''', title="Inventory 5 - explore junction trade-offs"))

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
                                  choice=choice)
    print(selected.summary)

    recompiled = w.load_and_compile(selected.artefacts["inventory"], wanted,
                                    "out/recompiled")
    print(recompiled.summary)

    order_items = w.order_items_for(selected.artefacts["inventory"], "out/order_items")
    print(order_items.summary)
''', title="Inventory 6 - commit to a choice and rebuild"))

    out.append(md("""
## Plan the order

Eligibility, pooling and export against a written product profile.

The order items are **assembly-ready substrates**, not bare inserts: each is wrapped so BbsI
releases the level-0 fragment with the correct exposed ends, and each is digest-verified
before export. An A-module insert starts with its own `AATG` fusion site, while the substrate
must expose `CTCA` — ordering the insert would order something that cannot assemble.

**A local check is not vendor approval.** The shipped profile carries this project's
historical price constants and deliberately **refuses to price**, because an undated number
presented as current is worse than no number at all. Eligibility, pooling and export work
without one; supply a dated profile to get an estimate.
"""))

    out.append(code('''
#@markdown Leave blank for the current published oPools rules (unpriced).
#@markdown Supply a dated profile JSON to get a cost estimate.
profile_file = ""  #@param {type:"string"}

from clippr.ordering import CURRENT_OPOOL_50PMOL, load_profile

profile = load_profile(profile_file) if profile_file.strip() else CURRENT_OPOOL_50PMOL

# The order comes from the inventory you selected, never from whatever happens to be left in
# notebook memory. Taking sequences from the synthesis route's `lib` here meant a user could
# explore one inventory and order something unrelated to it.
order = None
if order_items is None:
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

    cost = order.data["cost"]
    if cost["available"]:
        print(f"   estimate {cost['total']} {cost['currency']} "
              f"(priced {cost['priced_on']}) -- an estimate, not a quote")
    else:
        print(f"   cost unavailable: {cost['reason']}")
    print()
    for rule in order.data["eligibility"]["unresolved"]:
        print(f"   unresolved: {rule}")
''', title="Inventory 7 - eligibility and pool plan"))

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
        shutil.make_archive("clippr_package", "zip", "out")
        files.download("clippr_package.zip")
    except ImportError:
        print()
        print("Not running in Colab -- the package is in out/.")
''', title="Inventory 8 - download the result package"))

    return out
