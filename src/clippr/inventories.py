"""Versioned module inventories: load, validate, compile, save, reload.

Named `inventories` rather than `inventory` because `clippr.inventory()` is already
the kit *index* function re-exported from `parts`; two things called `inventory` in
one namespace is a collision a reader should not have to resolve.

An inventory is a collection of modules that can be compiled into targets, plus the identity
of the realisation those modules belong to. The deposited kit is one inventory; a recoded kit
is another; a kit with redesigned interfaces is a third. All three load, validate and compile
through this one contract, which is what makes them interchangeable at the call site.

**Versions, not overwrites.** A recoded module is a new *version* of the same module identity.
The deposited record is never replaced, so a result can always be traced back to what it came
from.

**An interface assignment belongs to an inventory.** Modules from inventories with different
interface assignments cannot be mixed: their four-base boundaries were chosen to interoperate
within one assignment, and a compile handed a mixture would emit a product that cannot ligate.
`compile_target` refuses such a mixture and names both versions rather than producing it.

**Frames are derived, never searched.** A module's reading frame follows from its position in
the assembled product and the fusion site that sets the product's frame. Trying all three
frames and keeping whichever has no stop codon is a heuristic dressed as a measurement, and it
succeeds often enough to be dangerous.

See `docs/workflow_contract.md` for the field-by-field contract this implements.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

from .parts import DEFAULT_FUSION_SITE, FUSION_SITES, inventory as deposited_index
from .parts import select

INTERFACE = 4
SCHEMA_VERSION = 1

#: Coding frame the product runs in, set by the 5' fusion site of the first module. Both
#: deposited sites are four-base MoClo CDS sites whose first base completes the codon upstream
#: of the insert, so both put the first whole codon at product position 1.
_FUSION_FRAME = {"AATG": 1, "AGGT": 1}


class InventoryInvalid(Exception):
    """The records do not form a usable inventory, and the reason names the records."""


class IncompatibleInventory(Exception):
    """Modules from different interface assignments cannot be compiled together."""


@dataclass(frozen=True)
class ModuleRecord:
    """One module, in one realisation."""

    module_id: str
    version: str
    dna: str
    frame: int
    block: str
    fifth: str | None
    last: str | None
    five_interface: str
    three_interface: str
    source_sha256: str
    derived_from: str | None = None
    unchanged: bool = True

    @property
    def coding_interval(self) -> tuple[int, int]:
        """Half-open span of the whole codons this module carries, in its own frame."""
        n = (len(self.dna) - self.frame) // 3
        return self.frame, self.frame + 3 * n

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.dna.encode()).hexdigest()

    def as_dict(self) -> dict:
        start, end = self.coding_interval
        return {"module_id": self.module_id, "version": self.version, "dna": self.dna,
                "frame": self.frame, "coding_interval": [start, end],
                "block": self.block, "fifth": self.fifth, "last": self.last,
                "five_interface": self.five_interface,
                "three_interface": self.three_interface,
                "source_sha256": self.source_sha256, "derived_from": self.derived_from,
                "unchanged": self.unchanged, "sha256": self.sha256}

    @classmethod
    def from_dict(cls, d: dict) -> "ModuleRecord":
        return cls(module_id=d["module_id"], version=d["version"], dna=d["dna"],
                   frame=d["frame"], block=d["block"], fifth=d["fifth"], last=d["last"],
                   five_interface=d["five_interface"],
                   three_interface=d["three_interface"],
                   source_sha256=d["source_sha256"], derived_from=d.get("derived_from"),
                   unchanged=d.get("unchanged", True))


@dataclass
class Inventory:
    """A complete, versioned set of modules and the assignment they interoperate under."""

    label: str
    interface_assignment: str
    context: dict
    modules: dict[str, ModuleRecord] = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)

    @property
    def version(self) -> str:
        """Content identity: the same modules in the same versions give the same value."""
        digest = hashlib.sha256()
        for module_id in sorted(self.modules):
            record = self.modules[module_id]
            digest.update(f"{module_id}\x1f{record.version}\x1f{record.sha256}\x1e".encode())
        digest.update(f"interfaces={self.interface_assignment}".encode())
        return digest.hexdigest()[:16]

    def __len__(self) -> int:
        return len(self.modules)

    def as_dict(self) -> dict:
        return {"schema_version": SCHEMA_VERSION, "label": self.label,
                "version": self.version,
                "interface_assignment": self.interface_assignment,
                "context": self.context, "provenance": self.provenance,
                "modules": {k: v.as_dict() for k, v in sorted(self.modules.items())}}


def _frames_from_architecture() -> dict[str, int]:
    """Every module's frame, derived from where the kit's own recipes place it.

    A module's frame is a property of its position class, not of its sequence. The witness
    recipes below visit every position every module can occupy -- one target shape is not
    enough, because a target of a single repeated base never selects the specificity variants
    that read the other three.

    A module found at two different frames is a contradiction, not something to resolve by
    picking one: raise it.
    """
    frames: dict[str, int] = {}
    for target, fusion_site in sorted(set(witness_targets().values())):
        plan = select(target, fusion_site=fusion_site)
        if not plan.available:
            continue
        product_frame = _FUSION_FRAME[fusion_site]
        offset = 0
        for module in plan.modules:
            frame = (product_frame - offset) % 3
            prior = frames.get(module.plasmid_id)
            if prior is not None and prior != frame:
                raise InventoryInvalid(
                    f"{module.plasmid_id} occupies frame {prior} in one recipe and {frame} "
                    f"in another; a module with two frames cannot be recoded in one of them")
            frames[module.plasmid_id] = frame
            offset += module.insert_length - INTERFACE
    return frames


def witness_targets() -> dict[str, tuple[str, str]]:
    """For every deposited module, a target and fusion site that actually uses it.

    The 200-target corpus exercises 40 of the 42 modules: the two `AGGT` start variants are
    unreachable through the default fusion site, because they are identical to their `AATG`
    twins in block, specificity residues and length. Coverage of 40 is not coverage of 42, and
    a witness for every module is what makes the difference checkable rather than assumed.
    """
    from .parts import BUILDABLE_LENGTHS, code_to_base

    base_for = code_to_base()
    witnesses: dict[str, tuple[str, str]] = {}
    for fusion_site in FUSION_SITES:
        for n in sorted(BUILDABLE_LENGTHS):
            # Vary the target so different specificity variants get selected.
            for base in sorted(set(base_for.values())):
                target = base * n
                plan = select(target, fusion_site=fusion_site)
                if not plan.available:
                    continue
                for module in plan.modules:
                    witnesses.setdefault(module.plasmid_id, (target, fusion_site))
    return witnesses


def load_deposited(table_path: str | Path, context: dict | None = None) -> Inventory:
    """The 42 deposited modules, cross-checked against the kit index.

    Both directions are checked. A record present in the table but absent from the index is
    as much a defect as the reverse: it means the two disagree about what the kit contains,
    and compiling from either alone would silently use one side's view.
    """
    from .products import load_inserts

    sequences, source = load_inserts(table_path)
    index = {m.plasmid_id: m for m in deposited_index()}

    missing = sorted(set(index) - set(sequences))
    extra = sorted(set(sequences) - set(index))
    if missing or extra:
        raise InventoryInvalid(
            f"the kit index and the supplied records disagree: {len(missing)} indexed "
            f"modules have no sequence {missing[:3]}, {len(extra)} supplied sequences are "
            f"not in the index {extra[:3]}")

    wrong_length = [pid for pid, module in index.items()
                    if len(sequences[pid]) != module.insert_length]
    if wrong_length:
        raise InventoryInvalid(
            f"{len(wrong_length)} supplied sequences do not match the indexed insert "
            f"length: {wrong_length[:3]}")

    frames = _frames_from_architecture()

    modules = {}
    for pid, module in sorted(index.items()):
        dna = sequences[pid]
        if pid not in frames:
            raise InventoryInvalid(
                f"{pid} has no witness context in any kit recipe, so its reading frame "
                f"cannot be derived; it cannot be recoded safely")
        modules[pid] = ModuleRecord(
            module_id=pid, version="deposited", dna=dna, frame=frames[pid],
            block=module.block, fifth=module.fifth, last=module.last,
            five_interface=dna[:INTERFACE], three_interface=dna[-INTERFACE:],
            source_sha256=source.sha256, derived_from=None, unchanged=True)

    return Inventory(
        label="deposited", interface_assignment="deposited",
        context=dict(context or {}), modules=modules,
        provenance={"source_path": str(source.path), "source_sha256": source.sha256,
                    "index_modules": len(index)})


def synthesis_profile_notes(inv: Inventory) -> dict[str, list[str]]:
    """Modules that do not fit CLIPPR's synthesis profile, whatever their provenance.

    Separate from `validate`, which is about structural integrity. This is about fit to a
    profile, and a deposited module can be structurally perfect and still fall outside it:
    `pPR-1_D_LD5T` carries `TTTTT` and fails our max-homopolymer-4 rule.

    That rule is an engineering choice of ours, not a published requirement, so this is not a
    defect in the deposited kit -- but someone about to order that module should be told
    before the order, not after.

    Judged on the substrate that would actually be ordered, through the one shared validator.
    Warning a user about the bare insert while the order step refuses the substrate would put
    this function one layer inside the thing it exists to warn about.
    """
    from .substrates import substrate_problems

    return {module_id: substrate_problems(record.dna, record.block)
            for module_id, record in sorted(inv.modules.items())
            if substrate_problems(record.dna, record.block)}


def validate(inv: Inventory) -> list[str]:
    """Every structural problem with an inventory, named. Empty means usable.

    Structural only. For fit to the synthesis profile -- which a valid inventory can fail
    without being malformed -- see `synthesis_profile_notes`.
    """
    problems = []
    index = {m.plasmid_id: m for m in deposited_index()}

    for pid, record in sorted(inv.modules.items()):
        if record.module_id != pid:
            problems.append(f"{pid} is filed under an identity it does not carry "
                            f"({record.module_id})")
        if set(record.dna) - set("ACGT"):
            problems.append(f"{pid} contains non-ACGT bases")
        if record.dna[:INTERFACE] != record.five_interface:
            problems.append(f"{pid} 5' interface field disagrees with its own sequence")
        if record.dna[-INTERFACE:] != record.three_interface:
            problems.append(f"{pid} 3' interface field disagrees with its own sequence")
        if record.frame not in (0, 1, 2):
            problems.append(f"{pid} has frame {record.frame}")
        known = index.get(pid)
        if known is not None and len(record.dna) != known.insert_length:
            problems.append(f"{pid} is {len(record.dna)} nt, the kit index says "
                            f"{known.insert_length}")

    absent = sorted(set(index) - set(inv.modules))
    if absent:
        problems.append(f"{len(absent)} indexed modules are absent from this inventory: "
                        f"{absent[:3]}")
    return problems


def _stage_summary(junctions: list[dict], level0_reactions: int,
                   flanks: tuple[str, str]) -> dict:
    """Each reaction with its stage, ends, enzyme and whether it can be scored at all.

    `flanks` are the compiled product's own first and last interfaces. A level-0 reaction
    takes its ends from the level-0 destination, which is a fixed vector pair; a **level-1
    reaction takes them from the product**, because the released block fragment is the joined
    module inserts and the block vector contributes no bases (`assembly_spec.LEVEL1_BLOCK_ENDS`).

    **A level-1 reaction is still not scored**, for a reason narrower than the one this code
    used to give. The block geometry is established; the reaction's *participant list* is not,
    because the deposited BsaI reaction co-assembles parts this package does not compile.
    Fidelity depends on the whole competing overhang set, so the PPR-only figure is carried as
    a named diagnostic in its own field and never as this reaction's fidelity.
    """
    from .assembly_spec import ASSEMBLY_STAGES

    grouped: dict[int, dict] = {}
    for junction in junctions:
        entry = grouped.setdefault(junction["reaction"], {
            "reaction": junction["reaction"], "stage": junction["stage"],
            "internal_overhangs": []})
        entry["internal_overhangs"].append(junction["overhang"])

    out = []
    for index in sorted(grouped):
        entry = grouped[index]
        spec = ASSEMBLY_STAGES[entry["stage"]]
        ends = tuple(flanks) if spec.get("flanks_from_product") else spec["destination"]
        known_ends = ends is not None
        complete = spec["context_available"] and known_ends
        row = {
            **entry,
            "enzyme": spec["enzyme"], "matrix": spec["matrix"],
            "destination": list(ends) if ends else None,
            "ends_from": "product" if spec.get("flanks_from_product") else "destination vector",
            "geometry_established": bool(spec.get("geometry_established", complete)),
            "participants_established": complete,
            "reaction_overhangs": ([ends[0]] + entry["internal_overhangs"] + [ends[1]])
                                  if complete else None,
            "scorable": complete,
            "note": None,
        }
        if not complete and known_ends and spec.get("geometry_established"):
            # The geometry is known and the subset is worth reporting -- in its own field,
            # named for what it is. Putting it in `reaction_overhangs` would make an
            # incomplete set look like the reaction.
            row["block_subset_overhangs"] = ([ends[0]] + entry["internal_overhangs"]
                                             + [ends[1]])
            row["note"] = (
                "the block release geometry is established, but the reaction's participants "
                "are not: the deposited reaction co-assembles parts this package does not "
                "compile. `block_subset_overhangs` is a PPR-only diagnostic, not this "
                "reaction's fidelity")
        elif not complete:
            row["note"] = "this reaction's ends are not determined by the supplied records"
        out.append(row)
    return {"reactions": out, "level0_reactions": level0_reactions,
            "level1_reactions": sum(1 for r in out if r["stage"] == "level1")}


def compile_target(inv: Inventory, target: str,
                   fusion_site: str = DEFAULT_FUSION_SITE) -> dict:
    """Compile one target from this inventory, or say exactly what is missing.

    Refuses rather than guesses: a module absent from the inventory, or one belonging to a
    different interface assignment, names itself in the error.
    """
    plan = select(target, fusion_site=fusion_site)
    if not plan.available:
        return {"target": target, "available": False, "reason": plan.reason}

    ids = [m.plasmid_id for m in plan.modules]
    absent = [pid for pid in ids if pid not in inv.modules]
    if absent:
        raise IncompatibleInventory(
            f"inventory {inv.label}@{inv.version} lacks {len(absent)} module(s) this "
            f"target needs: {absent}")

    records = [inv.modules[pid] for pid in ids]
    # Which reaction each junction belongs to. Golden Gate needs unique overhangs *within one
    # reaction*, not across the whole product: the kit reuses ACTC at both 1A->B and 2A->B,
    # and that is legal precisely because the plan splits them into different reactions.
    # Checking uniqueness product-wide would reject the deposited kit itself.
    # A junction belongs to the reaction that actually creates it. Joins *inside* a
    # sub-assembly are made at level 0 by BbsI; the join *between* two sub-assemblies is a
    # block join, made at level 1 by BsaI, and is not present in any level-0 tube.
    #
    # The previous version assigned the block join to the preceding level-0 reaction. That
    # put CTTC, GTGA and CACG into BbsI reactions they never take part in, and produced a
    # spurious "third reaction bottleneck" in the 19S architecture.
    reaction_of: dict[int, dict] = {}
    position = 0
    for index, sub in enumerate(plan.sub_assemblies):
        for _ in range(len(sub) - 1):
            reaction_of[position] = {"reaction": index, "stage": "level0"}
            position += 1
        if position < len(plan.modules) - 1:
            reaction_of[position] = {"reaction": len(plan.sub_assemblies),
                                     "stage": "level1"}
            position += 1

    records = [inv.modules[pid] for pid in ids]
    sequence = records[0].dna
    junctions = []
    for left, right in zip(records, records[1:]):
        if left.three_interface != right.five_interface:
            raise IncompatibleInventory(
                f"{left.module_id} ends {left.three_interface} but {right.module_id} starts "
                f"{right.five_interface}; they belong to different interface assignments")
        junctions.append({"left": left.module_id, "right": right.module_id,
                          "overhang": left.three_interface,
                          "start": len(sequence) - INTERFACE,
                          **reaction_of.get(len(junctions),
                                            {"reaction": 0, "stage": "level0"})})
        sequence += right.dna[INTERFACE:]

    frame = _FUSION_FRAME[fusion_site]
    coding = sequence[frame:]
    coding = coding[:len(coding) // 3 * 3]
    return {
        "target": target, "available": True,
        "inventory_version": inv.version, "inventory_label": inv.label,
        "interface_assignment": inv.interface_assignment,
        "fusion_site": fusion_site,
        "modules": [(r.module_id, r.version) for r in records],
        "product": sequence, "product_nt": len(sequence),
        "coding_interval": [frame, frame + len(coding)],
        "junctions": junctions,
        "reactions": len(plan.sub_assemblies),
        "stages": _stage_summary(junctions, len(plan.sub_assemblies),
                                 (records[0].five_interface, records[-1].three_interface)),
        "boundary_scope": ("joined module inserts only; no acceptor backbone supplied, so "
                           "this is not an expression construct"),
    }


def save(inv: Inventory, path: str | Path) -> Path:
    """Write an inventory so it can be reopened and compiled without regeneration."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(inv.as_dict(), indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    return out


def load(path: str | Path) -> Inventory:
    """Reopen a saved inventory. The recorded version must match what the records imply."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("schema_version") != SCHEMA_VERSION:
        raise InventoryInvalid(
            f"inventory schema {data.get('schema_version')}, this build reads "
            f"{SCHEMA_VERSION}")
    inv = Inventory(
        label=data["label"], interface_assignment=data["interface_assignment"],
        context=data.get("context", {}), provenance=data.get("provenance", {}),
        modules={k: ModuleRecord.from_dict(v) for k, v in data["modules"].items()})
    if inv.version != data["version"]:
        raise InventoryInvalid(
            f"saved version {data['version']} but the records hash to {inv.version}; the "
            f"file has been edited since it was written")
    return inv


def derive(inv: Inventory, label: str, replacements: dict[str, str],
           interface_assignment: str | None = None) -> Inventory:
    """A new inventory version with some modules replaced, the rest carried over unchanged.

    An unchanged module is a valid module. Carrying it forward marked `unchanged` keeps the
    denominator honest: coverage is over every module the inventory holds, not only the ones
    an optimiser happened to move.
    """
    modules = {}
    for pid, record in inv.modules.items():
        new_dna = replacements.get(pid)
        if new_dna is None or new_dna == record.dna:
            modules[pid] = ModuleRecord(**{**record.__dict__, "unchanged": True})
            continue
        modules[pid] = ModuleRecord(
            module_id=pid, version=label, dna=new_dna, frame=record.frame,
            block=record.block, fifth=record.fifth, last=record.last,
            five_interface=new_dna[:INTERFACE], three_interface=new_dna[-INTERFACE:],
            source_sha256=record.source_sha256,
            derived_from=f"{record.version}:{record.sha256[:16]}", unchanged=False)
    return Inventory(
        label=label,
        interface_assignment=interface_assignment or inv.interface_assignment,
        context=dict(inv.context), modules=modules,
        provenance={**inv.provenance, "derived_from_version": inv.version,
                    "derived_from_label": inv.label})
