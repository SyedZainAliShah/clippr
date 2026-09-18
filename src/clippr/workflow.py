"""The six user-visible tasks, and one result package that ties them together.

Each task is callable on its own and writes persistent artefacts, so a long workflow can be
stopped, inspected and resumed rather than re-run. `docs/workflow_contract.md` §6 is the
contract this implements.

    1. design_for_synthesis   design complete CDSs from explicit RNA targets
    2. load_and_compile       load an inventory and compile targets from it
    3. recode_for_host        recode an inventory for a host context, interfaces fixed
    4. optimise_collection    optimise an inventory as a collection
    5. explore_interfaces     explore junction/codon trade-offs, return the observed front
    6. plan_order             check eligibility, plan pools, export the package
    7. level1_readiness       is a level-1 reaction fully specified, and how well does it ligate

**Selecting an alternative changes the inventory.** Task 5 returns candidates whose junction
assignments differ, and choosing one produces a *different interface version*. Products
compiled from the old version cannot be mixed with modules from the new one, so every result
that offers alternatives says so at the point of choosing rather than leaving it to be
discovered at the bench.

**A failed target stays in the report.** A package that lists only what worked has quietly
changed its own denominator, so failures are carried with their reasons — and a run where
every target failed still produces provenance, because "what produced this directory of
nothing" is exactly the question that then needs answering.
"""
from __future__ import annotations

import json
import shutil
from dataclasses import dataclass, field
from pathlib import Path

PACKAGE_SCHEMA_VERSION = 1


@dataclass
class WorkflowResult:
    """One task's output, its provenance, and whatever it could not do."""

    task: str
    ok: bool
    summary: str
    artefacts: dict[str, str] = field(default_factory=dict)
    data: dict = field(default_factory=dict)
    failures: list[dict] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {"task": self.task, "ok": self.ok, "summary": self.summary,
                "artefacts": self.artefacts, "data": self.data,
                "failures": self.failures}


def design_for_synthesis(targets, codon_table, outdir, **kwargs) -> WorkflowResult:
    """Task 1 — design complete CDSs for synthesis, one per target."""
    from .library import design_library, write_library

    out = Path(outdir)
    library = design_library(targets, outdir=None, codon_table=codon_table, **kwargs)
    artefacts = write_library(library, out) if library.designs else {}
    failures = [{"target": t, "reason": why} for t, why in sorted(library.failed.items())]
    return WorkflowResult(
        task="design_for_synthesis",
        ok=bool(library.designs) and not failures,
        summary=(f"{len(library.designs)} of {len(list(targets))} targets designed; "
                 f"{len(failures)} failed"),
        artefacts={k: str(v) for k, v in artefacts.items()},
        data={"designed": sorted(library.designs)}, failures=failures)


def load_and_compile(inventory_path, targets, outdir, fusion_site: str = "AATG"
                     ) -> WorkflowResult:
    """Task 2 — load a saved or deposited inventory and compile targets from it."""
    from . import inventories as inv

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = Path(inventory_path)
    library = (inv.load_deposited(path) if path.suffix.lower() in (".xlsx", ".csv")
               else inv.load(path))

    compiled, failures = {}, []
    for target in targets:
        try:
            got = inv.compile_target(library, target, fusion_site=fusion_site)
        except inv.IncompatibleInventory as exc:
            failures.append({"target": target, "reason": str(exc)})
            continue
        if not got.get("available"):
            failures.append({"target": target, "reason": got.get("reason", "unavailable")})
            continue
        compiled[target] = got

    written = out / "compiled_targets.json"
    written.write_text(json.dumps(compiled, indent=2, sort_keys=True) + "\n",
                       encoding="utf-8")
    notes = inv.synthesis_profile_notes(library)
    return WorkflowResult(
        task="load_and_compile", ok=bool(compiled) and not failures,
        summary=(f"{len(compiled)} of {len(list(targets))} targets compiled from "
                 f"{library.label}@{library.version}"),
        artefacts={"compiled_targets": str(written)},
        data={"inventory_version": library.version, "inventory_label": library.label,
              "interface_assignment": library.interface_assignment,
              "modules": len(library),
              "synthesis_profile_notes": notes},
        failures=failures)


def recode_for_host(inventory_path, codon_table, outdir, *, label: str,
                    profile=None, **kwargs) -> WorkflowResult:
    """Task 3 — recode an inventory for a host context, both interfaces frozen.

    `profile` is the synthesis policy, reaching both the solver and the validator.
    Default reproduces the shipped result exactly. It was previously reachable only
    through `**kwargs`, which meant a caller had to know it existed; it is explicit
    now, and its version is recorded on the result.
    """
    from . import inventories as inv
    from .recoding import recode_inventory

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = Path(inventory_path)
    library = (inv.load_deposited(path) if path.suffix.lower() in (".xlsx", ".csv")
               else inv.load(path))

    report = recode_inventory(library, codon_table, label=label, profile=profile,
                              **kwargs)
    saved = inv.save(report.inventory, out / f"inventory_{label}.json")
    written = out / "recoding_report.json"
    written.write_text(json.dumps(report.as_dict(), indent=2, sort_keys=True) + "\n",
                       encoding="utf-8")
    return WorkflowResult(
        task="recode_for_host", ok=not report.infeasible,
        summary=(f"{len(report.improved)} improved, {len(report.unchanged)} unchanged, "
                 f"{len(report.infeasible)} infeasible, {len(report.not_attempted)} not "
                 f"attempted"),
        artefacts={"inventory": str(saved), "report": str(written)},
        data={"from_version": library.version, "to_version": report.inventory.version,
              "synthesis_profile": report.profile,
              "budget_exhausted": report.budget_exhausted},
        failures=[{"module": m, "reason": "infeasible"} for m in report.infeasible])


def optimise_collection(inventory_path, codon_table, outdir, *, mode: str = "greedy",
                        profile=None, **kwargs) -> WorkflowResult:
    """Task 4 — optimise an inventory as a collection, not module by module."""
    from . import inventories as inv
    from .library_search import optimise_library

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    library = inv.load(Path(inventory_path))
    result = optimise_library(library, codon_table, mode=mode, profile=profile,
                              **kwargs)
    saved = inv.save(result.inventory, out / f"inventory_{mode}.json")
    written = out / f"search_{mode}.json"
    written.write_text(json.dumps(result.as_dict(), indent=2, sort_keys=True) + "\n",
                       encoding="utf-8")
    return WorkflowResult(
        task="optimise_collection", ok=True,
        summary=(f"{mode}: {'improved' if result.improved else 'returned the incumbent'}; "
                 f"{len(result.accepted)} module(s) accepted, "
                 f"{result.as_dict()['completion']}"),
        artefacts={"inventory": str(saved), "report": str(written)},
        data={"mode": mode, "improved": result.improved,
              "from_version": library.version, "to_version": result.inventory.version,
              "budget_exhausted": result.budget_exhausted})


def explore_interfaces(inventory_path, targets, codon_table, outdir,
                       profile=None, **kwargs) -> WorkflowResult:
    """Task 5 — explore junction/codon trade-offs and return the observed front.

    The result carries an explicit warning: taking any alternative changes the interface
    version, and every product must then be recompiled from the new inventory.
    """
    from . import inventories as inv, joint_search as js
    from .codons import complete_table

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    library = inv.load(Path(inventory_path))
    table = complete_table(codon_table, 1)
    result = js.search(library, list(targets), table, profile=profile, **kwargs)

    written = out / "interface_front.json"
    written.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                       encoding="utf-8")
    return WorkflowResult(
        task="explore_interfaces", ok=bool(result["observed_front"]),
        summary=(f"{result['evaluated']} evaluated, {result['feasible']} feasible, "
                 f"{len(result['observed_front'])} on the observed front; "
                 f"{result['completion']}"),
        artefacts={"front": str(written)},
        data={"incumbent": result["incumbent"]["objectives"],
              "recommended": (result["recommended"]["objectives"]
                              if result["recommended"] else None),
              "recommendation_policy": result["recommendation_policy"],
              "fidelity_degenerate": result["fidelity_degenerate"],
              "scope": result["scope"],
              "selecting_an_alternative": (
                  "changes the interface version; every target must be recompiled from the "
                  "selected inventory, and modules from the previous version cannot be "
                  "mixed with it")})


def _profile_name(profile) -> str:
    from .synthesis_profile import resolve

    return resolve(profile).name


def select_interface(inventory_path, front_path, outdir, *, choice: str = "recommended",
                     index: int | None = None, codon_table=None,
                     profile=None, allow_mismatch: bool = False) -> WorkflowResult:
    """Task 5b -- commit to one candidate from the front and save its inventory.

    Exploring a front and then ordering is only connected if the chosen candidate becomes a
    real, saved inventory. Without this step a user could read a recommendation and then reach
    an order built from something else entirely, which is exactly what happened: the notebook's
    order cell took its sequences from the synthesis route's `lib.oligos()`.

    `choice="incumbent"` is a valid path and is the honest default when the front offers
    nothing worth the interface change.
    """
    import json as _json

    from . import inventories as inv, joint_search as js

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    library = inv.load(Path(inventory_path))
    front = _json.loads(Path(front_path).read_text(encoding="utf-8"))

    # A front's objective values were measured against one inventory under one context. The
    # same assignment will happily rebuild against a *different* inventory -- the overhang
    # names exist in both -- and every reported number then describes the wrong DNA. This
    # once returned the recoded run's adaptation of 0.672749 for deposited sequences that
    # actually score 0.219935.
    table = _resolve_table(codon_table)
    # Scoring knobs come from the front that recorded them. Validating a k=10 front with the
    # default k=20 refused a perfectly valid search, and `allow_mismatch=True` is not the
    # answer -- that suppresses the check rather than honouring the context.
    binding = front.get("binding") or {}
    # The synthesis policy travels with the front for the same reason k does: a
    # candidate rechecked under a different policy is a candidate nobody measured.
    recorded_profile = binding.get("synthesis_profile")
    if recorded_profile:
        from .synthesis_profile import from_dict as _profile_from_dict

        profile = _profile_from_dict(recorded_profile)
    k = binding.get("k", 20)
    genetic_code = binding.get("genetic_code", 1)
    matrix = binding.get("matrix", "BsaI-HFv2")
    destination = binding.get("destination", "level0")
    mismatches = js.binding_mismatches(binding, library, table,
                                       targets=front.get("targets"),
                                       genetic_code=genetic_code, matrix=matrix,
                                       destination=destination, k=k)
    if mismatches and not allow_mismatch:
        return WorkflowResult(
            task="select_interface", ok=False,
            summary=("front does not belong to this inventory or context; refusing to "
                     "publish its objectives against these sequences"),
            failures=[{"stage": "binding", "reason": m} for m in mismatches])

    if choice == "incumbent":
        chosen = front["incumbent"]
    elif index is not None:
        chosen = front["observed_front"][index]
    else:
        chosen = front["recommended"]
    if chosen is None:
        return WorkflowResult(task="select_interface", ok=False,
                              summary="no candidate to select", failures=[
                                  {"stage": "selection", "reason": "front is empty"}])

    # Rebuild against the targets the front was computed over, so the classes match the
    # assignment being applied. A front that did not record them falls back to the three
    # architectures, which between them cover every junction class the kit has.
    classes = js.junction_classes(
        library, front.get("targets")
        or ["AAAAUGUGG", "UUACACGUGCGUAC", "CUAUCACAUCACAUAAGCG"])

    selected = js.apply_assignment(library, classes, chosen["assignment"],
                                   label=f"selected-{choice}")

    # Recomputed from the sequences actually delivered, never copied from the front. A value
    # carried across is a value nobody checked against what shipped. This runs *before* the
    # inventory is written: an earlier version saved first and returned no artefact on
    # failure, which left a rejected inventory on disk for the next step to pick up.
    verified = js.evaluate(library, classes, chosen["assignment"],
                           front.get("targets") or [], table, k=k, matrix=matrix,
                           destination=destination, profile=profile)
    if not verified.feasible:
        # Never publish an artefact the verifier just rejected. This returned ok=True
        # with objectives=None, which reads as a successful selection of a candidate
        # that cannot be built.
        return WorkflowResult(
            task="select_interface", ok=False,
            summary=(f"the selected candidate does not verify under the front's own "
                     f"policy ({_profile_name(profile)}); nothing was published"),
            data={"choice": choice, "index": index,
                  "synthesis_profile": _profile_name(profile)},
            failures=[{"stage": "verification", "reason": verified.reason}])

    saved = inv.save(selected, out / "inventory_selected.json")
    drift = None
    if verified.objectives and chosen.get("objectives"):
        drift = {name: [chosen["objectives"].get(name), verified.objectives.get(name)]
                 for name in verified.objectives
                 if chosen["objectives"].get(name) != verified.objectives.get(name)}
    record = out / "selection.json"
    record.write_text(_json.dumps({
        "choice": choice, "index": index,
        "assignment": chosen["assignment"],
        "objectives": verified.objectives,
        "objectives_as_recorded_in_the_front": chosen.get("objectives"),
        "objective_drift": drift,
        "binding_mismatches": mismatches,
        "synthesis_profile": _profile_name(profile),
        "from_inventory": library.version, "selected_inventory": selected.version,
        "interface_assignment": selected.interface_assignment,
        "requires_recompilation": True,
        "note": ("this is a different interface version; products compiled from the previous "
                 "inventory cannot be mixed with modules from this one"),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return WorkflowResult(
        task="select_interface", ok=True,
        summary=(f"selected the {choice}; inventory {library.version} -> "
                 f"{selected.version}, recompilation required"),
        artefacts={"inventory": str(saved), "selection": str(record)},
        data={"choice": choice, "from_version": library.version,
              "synthesis_profile": _profile_name(profile),
              "to_version": selected.version,
              "objectives": verified.objectives,
              "objectives_recomputed": True,
              "objective_drift": drift,
              "binding_mismatches": mismatches})


def level1_readiness(block_subset, outdir, *, roles=None, evidence: str = "",
                     matrix: str = "BsaI-HFv2") -> WorkflowResult:
    """Is a level-1 reaction fully specified, and if so, how well does it ligate?

    CLIPPR derives the **block join geometry** from the deposited kit, but a level-1 reaction
    also needs the linker, editing domain, bridging parts and final cassette, and the deposited
    material does not supply their ends. That gap is real and documented; what was missing was
    a way for a user who *has* those parts to say so and get an answer.

    Supply `roles` mapping each of `assembly_spec.LEVEL1_REQUIRED_ROLES` to the overhangs it
    contributes, plus `evidence` naming where they came from. An incomplete list is refused
    with the roles it lacks — never scored, because a fidelity number computed over a partial
    reaction describes a reaction nobody is running.

    The refusal is the common case today and is a successful outcome of this task in the sense
    that matters: it tells a user exactly what to go and find.
    """
    import json as _json

    from . import assembly_spec as spec
    from .overhangs import fidelity_components, set_fidelity

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    participants = spec.level1_participants(block_subset, roles=roles, evidence=evidence)

    scored = None
    if participants["participants_established"]:
        # The **distinct** junctions, not every contributed end. Each internal overhang is
        # supplied twice -- once by the part on each side -- and scoring the list as given
        # makes every junction compete against a perfect copy of itself. On the chain used in
        # the tests that is the difference between 0.0039 and 0.9940: not a smaller number,
        # a meaningless one.
        ends = participants["overhangs"]
        junctions = sorted(set(ends))
        forward, reverse = fidelity_components(junctions, matrix)
        scored = {"fidelity": set_fidelity(junctions, matrix),
                  "forward": forward, "reverse": reverse,
                  "matrix": matrix, "junctions": junctions,
                  "contributed_ends": len(ends),
                  "note": ("predicted from a published ligation table over the distinct "
                           "junctions; not a measurement of this reaction")}

    written = out / "level1_readiness.json"
    written.write_text(_json.dumps({**participants, "scored": scored}, indent=2,
                                   sort_keys=True) + "\n", encoding="utf-8")

    if scored:
        summary = (f"level-1 reaction fully specified: {len(scored['junctions'])} distinct "
                   f"junctions from {scored['contributed_ends']} contributed ends; "
                   f"predicted fidelity {scored['fidelity']:.4f} on {matrix}")
    else:
        summary = (f"level-1 reaction not fully specified: "
                   f"{participants['why_incomplete']}")

    return WorkflowResult(
        task="level1_readiness", ok=bool(scored), summary=summary,
        artefacts={"readiness": str(written)},
        data={"participants_established": participants["participants_established"],
              "roles_required": participants["roles_required"],
              "roles_missing": participants["roles_missing"],
              "evidence": participants["evidence"], "scored": scored},
        failures=([] if scored else
                  [{"stage": "participants", "reason": participants["why_incomplete"]}]))


def _resolve_table(codon_table):
    """The completed codon table, defaulting to the shipped Chlamydomonas one."""
    from .codons import complete_table

    if codon_table is None:
        path = (Path(__file__).resolve().parents[2] / "data" / "codon_tables"
                / "kazusa_3055.json")
        codon_table = json.loads(path.read_text(encoding="utf-8"))
    if codon_table and isinstance(next(iter(codon_table.values())), dict):
        return complete_table(codon_table, 1)
    return codon_table


def order_items_for(inventory_path, outdir, *, form: str = "assembly_ready",
                    profile=None) -> WorkflowResult:
    """Task 5c -- the order items for an inventory, each traceable to its module.

    `form="assembly_ready"` (the default) wraps each module so BbsI releases the level-0
    fragment with the correct exposed ends. **A bare insert is not an assembly substrate**: an
    A-module insert begins with its own `AATG`/`AGGT` fusion site, while the level-0 substrate
    must expose `CTCA`, which comes from the vector. Ordering the insert would order something
    that cannot assemble.

    `form="raw_inserts"` exports the bare inserts. That is a legitimate diagnostic and a vendor
    profile may well accept them; it is simply not the assembly substrate, and it is labelled.

    Every assembly-ready sequence is digested again before export, so a wrong wrapper fails
    here rather than reaching an order.
    """
    import json as _json

    from . import inventories as inv, substrates as sub

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    library = inv.load(Path(inventory_path))

    if form == "raw_inserts":
        items = [{"sequence": record.dna, "module": module_id,
                  "version": record.version, "design": library.label,
                  "inventory_version": library.version, "quantity": 1,
                  "form": "raw_insert",
                  "note": "bare insert; not an assembly substrate"}
                 for module_id, record in sorted(library.modules.items())]
        summary = f"{len(items)} raw inserts (diagnostic, not assembly-ready)"
        ends = None
    elif form == "assembly_ready":
        items = [{**s.as_dict(), "design": library.label,
                  "inventory_version": library.version, "form": "assembly_ready"}
                 for s in sub.substrates_for(library, profile)]
        failing = [i for i in items if i["synthesis_problems"]]
        summary = (f"{len(items)} assembly-ready substrates from "
                   f"{library.label}@{library.version}, each digest-verified; "
                   f"{len(failing)} outside the synthesis contract")
        ends = sorted({(i["five_overhang"], i["three_overhang"]) for i in items})
    else:
        raise ValueError(f"form {form!r} is not 'assembly_ready' or 'raw_inserts'")

    written = out / "order_items.json"
    written.write_text(_json.dumps(items, indent=2, sort_keys=True) + chr(10),
                       encoding="utf-8")
    # A substrate outside the declared band is not a successful order item. Reporting ok
    # while exporting one is how a 0.660 GC window reached a package marked ok=True.
    failures = [{"stage": "synthesis", "module": i["module"], "reason": problem}
                for i in items for problem in i.get("synthesis_problems", ())]
    return WorkflowResult(
        task="order_items_for", ok=bool(items) and not failures, summary=summary,
        artefacts={"order_items": str(written)},
        data={"inventory_version": library.version, "items": len(items), "form": form,
              "enzyme": "BbsI" if form == "assembly_ready" else None,
              "exposed_ends": ends,
              "synthesis_contract": ("evaluated on the wrapped sequence that would be "
                                     "ordered; the two intended BbsI sites are excluded "
                                     "by role")},
        failures=failures)


def plan_order(sequences_or_items, vendor_profile, outdir) -> WorkflowResult:
    """Task 6 — eligibility, pooling and an export package.

    `vendor_profile` is the **ordering** profile — pool sizes, oligo length limits, price — and
    is a different thing from the `profile=` every other task takes, which is the synthesis
    policy. The parameter was called `profile` here too, which put two unrelated meanings of
    one word on functions a user calls one after the other. Nothing misbehaved; it was a trap
    waiting for whoever passed the wrong one first, and there is no diagnostic that could
    distinguish them.
    """
    from .ordering import check_eligibility, estimate_cost, plan_pools, write_order_files

    items = [{"sequence": s} if isinstance(s, str) else dict(s)
             for s in sequences_or_items]
    sequences = [i["sequence"] for i in items]

    # Order-wide checks first (lengths, bases, whether there are enough oligos for one pool
    # at all), then pool, then check each resulting pool against the per-pool count limits.
    # Applying the per-pool maximum to the whole order refused every multi-pool request.
    eligibility = check_eligibility(sequences, vendor_profile)
    plan = plan_pools(items, vendor_profile) if eligibility["eligible"] else {
        "feasible": False, "reason": "order is not eligible under this vendor profile",
        "pools": []}

    pool_reports, pool_costs = [], []
    if plan["feasible"]:
        for pool in plan["pools"]:
            members = [m["sequence"] for m in pool["members"]]
            report = check_eligibility(members, vendor_profile, per_pool=True)
            pool_reports.append({"pool": pool["pool"], **report})
            pool_costs.append(estimate_cost(members, vendor_profile))
        if not all(r["eligible"] for r in pool_reports):
            plan = {"feasible": False, "pools": [],
                    "reason": "; ".join(v for r in pool_reports
                                        for v in r["violations"])}

    # Cost is the sum over the pools actually planned, not one call on the whole order.
    cost = _total_cost(pool_costs, vendor_profile) if plan["feasible"] else estimate_cost(
        sequences, vendor_profile)
    artefacts = write_order_files(plan, vendor_profile, outdir) if plan["feasible"] else {}

    failures = []
    if not eligibility["eligible"]:
        failures += [{"stage": "eligibility", "reason": v}
                     for v in eligibility["violations"]]
    if not plan["feasible"]:
        failures.append({"stage": "pooling", "reason": plan.get("reason", "infeasible")})

    return WorkflowResult(
        task="plan_order", ok=eligibility["eligible"] and plan["feasible"],
        summary=(f"{eligibility['oligos']} oligos, "
                 f"{'eligible' if eligibility['eligible'] else 'not eligible'}, "
                 f"{plan.get('pool_count', 0)} pool(s), cost "
                 f"{'estimated' if cost['available'] else 'unavailable'}"),
        artefacts={k: str(v) for k, v in artefacts.items()},
        data={"eligibility": eligibility, "cost": cost,
              "per_pool_eligibility": pool_reports,
              "pool_count": plan.get("pool_count", 0),
              "optimality": plan.get("optimality"),
              "vendor_approval": ("not obtained; these are local checks against a written "
                                  "vendor profile and no order was placed")},
        failures=failures)


def _digest(path: Path) -> str:
    """Content hash, so a relocated package can prove each artefact is the one recorded."""
    import hashlib

    if path.is_dir():
        digest = hashlib.sha256()
        for child in sorted(path.rglob("*")):
            if child.is_file():
                digest.update(child.relative_to(path).as_posix().encode())
                digest.update(child.read_bytes())
        return digest.hexdigest()
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _total_cost(pool_costs: list[dict], profile) -> dict:
    """Sum the per-pool estimates. A multi-pool order costs the sum of its pools."""
    from decimal import Decimal

    if not pool_costs:
        return {"available": False, "reason": "no pools planned",
                "currency": profile.currency, "price_status": profile.price_status}
    if not all(c["available"] for c in pool_costs):
        return next(c for c in pool_costs if not c["available"])
    total = sum((Decimal(c["total"]) for c in pool_costs), Decimal("0"))
    return {"available": True, "currency": profile.currency,
            "price_status": profile.price_status,
            "priced_on": profile.priced_on, "source_url": profile.source_url,
            "pools": len(pool_costs),
            "per_pool": [c["total"] for c in pool_costs],
            "total": str(total), "is_a_quote": False,
            "note": "sum over the planned pools; an estimate, not a quote"}


def write_package(results, outdir, *, inputs: dict | None = None) -> Path:
    """One coherent result package, written even when every task failed.

    A directory of failures with no statement of what produced them is the case this exists
    to prevent, so provenance is written first and unconditionally.
    """
    from . import manifest as _manifest

    out = Path(outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    missing: list[str] = []
    collisions: list[str] = []
    claimed: dict[str, str] = {}

    def portable(record: dict, position: int) -> dict:
        """Copy every artefact into the package under a per-instance path, and hash it.

        Two runs of the same task write the same default filename. Keying the destination on
        the task name alone meant the second copy silently overwrote the first, leaving a
        manifest whose first task pointed at the second task's result -- while still
        reporting ok=True. The position prefix gives each task instance its own directory,
        and any remaining clash is recorded rather than resolved by overwriting.
        """
        moved, digests = {}, {}
        for name, value in record.get("artefacts", {}).items():
            source = Path(value)
            if not source.exists():
                missing.append(f"{record['task']}.{name}: {value}")
                moved[name] = None
                continue
            try:
                moved[name] = source.resolve().relative_to(out).as_posix()
                digests[name] = _digest(source)
                continue
            except ValueError:
                pass
            target = out / "artefacts" / f"{position:02d}_{record['task']}" / source.name
            key = target.resolve().as_posix()
            if key in claimed and claimed[key] != str(source.resolve()):
                collisions.append(f"{record['task']}.{name} would overwrite "
                                  f"{claimed[key]}")
            claimed[key] = str(source.resolve())
            target.parent.mkdir(parents=True, exist_ok=True)
            if source.is_dir():
                shutil.copytree(source, target, dirs_exist_ok=True)
            else:
                shutil.copy2(source, target)
            moved[name] = target.resolve().relative_to(out).as_posix()
            digests[name] = _digest(target)
        return {**record, "artefacts": moved, "artefact_sha256": digests,
                "instance": position}

    payload = {
        "schema_version": PACKAGE_SCHEMA_VERSION,
        "software": _manifest._software(),
        "inputs": inputs or {},
        "artefact_paths": ("relative to this package directory; every named "
                           "artefact is copied in, so the package resolves "
                           "after being moved"),
        "tasks": [portable(r.as_dict(), i) for i, r in enumerate(results)],
        # A package that names an artefact it does not contain is not a complete success,
        # however well the tasks themselves went.
        "ok": all(r.ok for r in results) and not missing and not collisions,
        "missing_artefacts": missing,
        "artefact_collisions": collisions,
        "failures": ([{"task": r.task, **f} for r in results for f in r.failures]
                     + [{"task": "write_package", "stage": "artefact",
                         "reason": f"referenced but absent: {m}"} for m in missing]
                     + [{"task": "write_package", "stage": "artefact",
                         "reason": c} for c in collisions]),
        "scope": ("sequence design and predictions only; nothing here is biologically "
                  "validated and no order was placed"),
    }
    path = out / "result_package.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path
