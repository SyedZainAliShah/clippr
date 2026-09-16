"""Order-product profiles: eligibility, pooling and cost estimates, with their limits stated.

Local product support. **None of this is vendor approval**, and no function here places or
submits an order. A profile is a written-down reading of a vendor's published rules at a given
date; the vendor's current documentation is always the authority.

**Prices expire, so they carry a date or they are not offered.** A profile whose
`priced_on` is absent reports cost as *unavailable* while eligibility, pooling and export
continue to work — the useful half of this module does not depend on a number that goes stale.
A price found in some other tool's source is not a source; it is that tool's copy of one.

**Money is `Decimal`.** Binary floats cannot represent 0.01, and a tiered calculation that
accumulates float error produces an estimate whose last digits are noise.

**Eligibility is judged on the wrapped order sequence**, the thing that would actually be
synthesised — not on the bare coding sequence. They differ by the enzyme wrapper, which is
where several length limits bite.

**No filler, ever.** If a request cannot meet a profile's minimum oligo count, the answer is
"no eligible plan" naming the constraint. Inventing padding oligos to reach a threshold would
turn an honest refusal into an order the user did not ask for.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

CENT = Decimal("0.01")


@dataclass(frozen=True)
class PriceTier:
    """One tier: applies while the oligo count is within [min_oligos, max_oligos]."""

    min_oligos: int
    max_oligos: int | None
    price: Decimal
    per: str = "pool"                    # "pool" | "oligo" | "base"

    def cost(self, oligos: int, bases: int) -> Decimal:
        if self.per == "pool":
            return self.price
        if self.per == "oligo":
            return self.price * oligos
        if self.per == "base":
            return self.price * bases
        raise ValueError(f"unknown price basis {self.per!r}")


@dataclass(frozen=True)
class ProductProfile:
    """A vendor product's machine-checkable rules, plus the ones that are not.

    `unresolved_rules` is not decoration. Every vendor has requirements no local check can
    settle — secondary structure judgements, acceptable-use terms, region-specific
    restrictions. Listing them keeps "we checked what we could" distinct from "this will be
    accepted", which is the distinction a passing eligibility report would otherwise blur.
    """

    name: str
    vendor: str
    region: str
    currency: str
    min_oligo_nt: int
    max_oligo_nt: int
    min_oligos: int
    max_oligos: int
    scales: tuple[str, ...] = ()
    price_tiers: tuple[PriceTier, ...] = ()
    priced_on: str | None = None                 # ISO date the prices were read
    source_url: str | None = None
    price_status: str = "unavailable"            # "current" | "historical" | "unavailable"
    modifications: dict = field(default_factory=dict)
    quantity_assumptions: str = ""
    unresolved_rules: tuple[str, ...] = ()
    allowed_bases: str = "ACGT"

    @property
    def prices_usable(self) -> bool:
        return bool(self.price_tiers) and self.priced_on is not None

    def as_dict(self) -> dict:
        return {
            "name": self.name, "vendor": self.vendor, "region": self.region,
            "currency": self.currency, "price_status": self.price_status,
            "priced_on": self.priced_on, "source_url": self.source_url,
            "limits": {"oligo_nt": [self.min_oligo_nt, self.max_oligo_nt],
                       "oligos": [self.min_oligos, self.max_oligos],
                       "scales": list(self.scales)},
            "price_tiers": [{"min_oligos": t.min_oligos, "max_oligos": t.max_oligos,
                             "price": str(t.price), "per": t.per}
                            for t in self.price_tiers],
            "modifications": self.modifications,
            "quantity_assumptions": self.quantity_assumptions,
            "unresolved_rules": list(self.unresolved_rules),
        }


#: The figures this package has historically used, kept **as history**.
#:
#: 109.00 EUR per pool and 1.63 EUR per oligo for 5' phosphorylation were measured across the
#: 200-design corpus and documented in `export.py`. They are a list price from an
#: unrecorded date, not a quote and not fetched from a vendor. `price_status` says
#: `historical` so nothing downstream can present them as current, and `priced_on` is absent
#: so `prices_usable` is False: a caller wanting a number must supply a dated table.
HISTORICAL_OPOOL = ProductProfile(
    name="oPools DNA (historical constants)",
    vendor="IDT", region="EU", currency="EUR",
    min_oligo_nt=20, max_oligo_nt=350,
    min_oligos=2, max_oligos=384,
    scales=("50 pmol",),
    price_tiers=(PriceTier(min_oligos=1, max_oligos=None, price=Decimal("109.00")),),
    priced_on=None,
    source_url=None,
    price_status="historical",
    modifications={"5' phosphorylation": {"per_oligo": "1.63", "currency": "EUR"}},
    quantity_assumptions="one pool, 50 pmol per oligo, as measured across the corpus",
    unresolved_rules=(
        "acceptable-use and sequence-screening terms are not machine-checkable here",
        "secondary-structure and synthesis-difficulty judgements are the vendor's",
        "region, tax and shipping depend on the account placing the order",
        "the length and count limits here were not re-read from current vendor "
        "documentation and must be confirmed before ordering",
    ),
)


#: Current oPools rules, **separated from price provenance**.
#:
#: Rule provenance and price provenance are different things and expire differently. An old
#: length limit must not stand in for a current rule merely because pricing is switched off:
#: the historical profile's 20 nt minimum passes two 30-base oligos that the current published
#: minimum of 40 would reject, which is a false positive in the direction that matters.
#:
#: Values as read from IDT's published oPools specifications on 2026-09-15 (40-350 bases per
#: oligo; 2-384 oligos per pool at the 50 pmol scale). **Confirm against the current page
#: before ordering** -- this is a written-down reading, not a live query, and other scales
#: carry different count ranges. Prices are deliberately absent, so cost reports unavailable
#: while eligibility, pooling and export continue to work.
CURRENT_OPOOL_50PMOL = ProductProfile(
    name="oPools DNA, 50 pmol (current published rules, unpriced)",
    vendor="IDT", region="EU", currency="EUR",
    min_oligo_nt=40, max_oligo_nt=350,
    min_oligos=2, max_oligos=384,
    scales=("50 pmol",),
    price_tiers=(),
    priced_on=None,
    source_url=("https://www.idtdna.com/pages/products/custom-dna-rna/dna-oligos/"
                "custom-dna-oligos/opools-oligo-pools"),
    price_status="unavailable",
    modifications={},
    quantity_assumptions="50 pmol per oligo; other scales have different count ranges",
    unresolved_rules=(
        "rules transcribed from the published specifications on 2026-09-15; confirm against "
        "the current page before ordering",
        "only the 50 pmol scale is modelled; other scales carry different per-pool counts",
        "acceptable-use and sequence-screening terms are not machine-checkable here",
        "secondary-structure and synthesis-difficulty judgements are the vendor's",
        "region, tax and shipping depend on the account placing the order",
        "no price is supplied, so cost is unavailable by design",
    ),
)


def load_profile(path: str | Path) -> ProductProfile:
    """Read a profile from JSON, including a dated price table a caller supplies."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    tiers = tuple(PriceTier(min_oligos=t["min_oligos"], max_oligos=t.get("max_oligos"),
                            price=Decimal(str(t["price"])), per=t.get("per", "pool"))
                  for t in data.get("price_tiers", ()))
    return ProductProfile(
        name=data["name"], vendor=data["vendor"], region=data.get("region", ""),
        currency=data.get("currency", ""),
        min_oligo_nt=data["min_oligo_nt"], max_oligo_nt=data["max_oligo_nt"],
        min_oligos=data["min_oligos"], max_oligos=data["max_oligos"],
        scales=tuple(data.get("scales", ())), price_tiers=tiers,
        priced_on=data.get("priced_on"), source_url=data.get("source_url"),
        price_status=data.get("price_status", "unavailable"),
        modifications=data.get("modifications", {}),
        quantity_assumptions=data.get("quantity_assumptions", ""),
        unresolved_rules=tuple(data.get("unresolved_rules", ())),
        allowed_bases=data.get("allowed_bases", "ACGT"))


def check_eligibility(sequences: list[str], profile: ProductProfile, *,
                      per_pool: bool = False) -> dict:
    """Machine-checkable violations, warnings and unresolved rules, kept apart.

    They are different kinds of statement and merging them is how "eligible" comes to mean
    less than a reader assumes. A violation means this order cannot be placed as described; a
    warning means it can but something is worth looking at; an unresolved rule means no local
    check settles it either way.
    """
    violations, warnings = [], []

    if not sequences:
        violations.append("no sequences supplied")
        return {"eligible": False, "violations": violations, "warnings": warnings,
                "unresolved": list(profile.unresolved_rules),
                "eligible_scales": [], "oligos": 0, "total_bases": 0}

    count = len(sequences)
    # The count limits are **per pool**, not per order. Applying the per-pool maximum to the
    # whole order refused any request larger than one pool -- nine oligos against a four-per-
    # pool product failed outright, even though three pools of three is a valid plan. Pass
    # `per_pool=True` when the sequences are one pool; the default judges an order that will
    # be split.
    if per_pool:
        if count < profile.min_oligos:
            violations.append(
                f"{count} oligos is below the product minimum of {profile.min_oligos} "
                f"per pool; no filler is added to reach it")
        if count > profile.max_oligos:
            violations.append(f"{count} oligos exceeds the product maximum of "
                              f"{profile.max_oligos} per pool")
    elif count < profile.min_oligos:
        violations.append(
            f"{count} oligos cannot fill even one pool, which needs at least "
            f"{profile.min_oligos}; no filler is added to reach it")

    allowed = set(profile.allowed_bases)
    for index, sequence in enumerate(sequences):
        bad = sorted(set(sequence.upper()) - allowed)
        if bad:
            violations.append(f"oligo {index} contains {bad} outside "
                              f"{profile.allowed_bases}")
        length = len(sequence)
        if length < profile.min_oligo_nt:
            violations.append(f"oligo {index} is {length} nt, below the minimum "
                              f"{profile.min_oligo_nt}")
        if length > profile.max_oligo_nt:
            violations.append(f"oligo {index} is {length} nt, above the maximum "
                              f"{profile.max_oligo_nt}")

    duplicates = len(sequences) - len(set(sequences))
    if duplicates:
        warnings.append(f"{duplicates} duplicate sequence(s); quantities and source "
                        f"mappings are preserved rather than merged")

    return {"eligible": not violations, "violations": violations, "warnings": warnings,
            "unresolved": list(profile.unresolved_rules),
            "eligible_scales": list(profile.scales) if not violations else [],
            "oligos": count, "total_bases": sum(len(s) for s in sequences)}


def estimate_cost(sequences: list[str], profile: ProductProfile, *,
                  modifications: tuple[str, ...] = (),
                  tax_rate: Decimal | None = None,
                  shipping: Decimal | None = None) -> dict:
    """A cost estimate, or an explicit statement that no usable price exists.

    Tax and shipping are **never assumed**. They depend on the account and destination, so
    they are included only when a caller supplies them, and their absence is reported rather
    than silently treated as zero.
    """
    if not profile.prices_usable:
        return {"available": False,
                "reason": (f"profile '{profile.name}' has price_status "
                           f"'{profile.price_status}' and no priced_on date; supply a dated "
                           f"price table to obtain an estimate"),
                "currency": profile.currency, "price_status": profile.price_status}

    count = len(sequences)
    bases = sum(len(s) for s in sequences)
    tier = next((t for t in profile.price_tiers
                 if t.min_oligos <= count and (t.max_oligos is None
                                               or count <= t.max_oligos)), None)
    if tier is None:
        return {"available": False,
                "reason": f"no price tier covers {count} oligos",
                "currency": profile.currency, "price_status": profile.price_status}

    subtotal = tier.cost(count, bases)
    lines = [{"item": "synthesis", "amount": str(subtotal.quantize(CENT, ROUND_HALF_UP))}]

    for name in modifications:
        rule = profile.modifications.get(name)
        if rule is None:
            return {"available": False,
                    "reason": f"profile has no rule for modification {name!r}",
                    "currency": profile.currency, "price_status": profile.price_status}
        amount = Decimal(str(rule["per_oligo"])) * count
        subtotal += amount
        lines.append({"item": name, "amount": str(amount.quantize(CENT, ROUND_HALF_UP))})

    total = subtotal
    if shipping is not None:
        total += shipping
        lines.append({"item": "shipping", "amount": str(shipping.quantize(CENT,
                                                                          ROUND_HALF_UP))})
    if tax_rate is not None:
        tax = (total * tax_rate).quantize(CENT, ROUND_HALF_UP)
        total += tax
        lines.append({"item": f"tax at {tax_rate}", "amount": str(tax)})

    return {"available": True, "currency": profile.currency,
            "price_status": profile.price_status, "priced_on": profile.priced_on,
            "source_url": profile.source_url,
            "oligos": count, "total_bases": bases,
            "tier": {"min_oligos": tier.min_oligos, "max_oligos": tier.max_oligos,
                     "per": tier.per, "price": str(tier.price)},
            "lines": lines,
            "total": str(total.quantize(CENT, ROUND_HALF_UP)),
            "tax_included": tax_rate is not None,
            "shipping_included": shipping is not None,
            "is_a_quote": False,
            "note": ("estimate from a dated profile; not a quote, and not vendor "
                     "confirmation that the order will be accepted")}


def plan_pools(items: list[dict], profile: ProductProfile) -> dict:
    """Group order items into pools that satisfy the profile, preserving every mapping.

    `items` are dicts with at least `sequence`; anything else (`design`, `module`, `quantity`)
    is carried through untouched, because a pool plan that loses which design an oligo came
    from cannot be ordered against.

    A deterministic first-fit by descending length. **Not claimed to be cost-optimal**: it is
    a feasible plan, and `validation/experiments/phasee_gate.py` compares it against exhaustive
    solutions on small cases to show how far off it can be.
    """
    if not items:
        return {"feasible": False, "reason": "no items to pool", "pools": []}

    too_long = [i for i, item in enumerate(items)
                if len(item["sequence"]) > profile.max_oligo_nt]
    too_short = [i for i, item in enumerate(items)
                 if len(item["sequence"]) < profile.min_oligo_nt]
    if too_long or too_short:
        return {"feasible": False,
                "reason": (f"{len(too_long)} oligo(s) exceed {profile.max_oligo_nt} nt and "
                           f"{len(too_short)} fall below {profile.min_oligo_nt} nt; no pool "
                           f"can contain them"),
                "pools": []}

    # Use the fewest pools the maximum allows, then spread items as evenly as possible.
    #
    # Plain first-fit packs each pool to the maximum and leaves the remainder in the last
    # one: nine items at four per pool became 4+4+1, and that pool of one is below the
    # minimum, so the whole order was refused -- even though 3+3+3 is perfectly feasible.
    # The heuristic was manufacturing the infeasibility it then reported.
    pool_count = -(-len(items) // profile.max_oligos)
    order = sorted(range(len(items)), key=lambda i: (-len(items[i]["sequence"]), i))
    pools: list[list[int]] = [[] for _ in range(pool_count)]
    for position, index in enumerate(order):
        pools[position % pool_count].append(index)

    undersized = [n for n, pool in enumerate(pools) if len(pool) < profile.min_oligos]
    if undersized:
        return {"feasible": False,
                "reason": (f"{len(items)} oligo(s) cannot fill {pool_count} pool(s) to the "
                           f"minimum of {profile.min_oligos} each, and no filler is invented "
                           f"to reach it"),
                "pools": []}

    built = []
    for number, pool in enumerate(pools):
        members = [{**items[i], "pool": number} for i in sorted(pool)]
        built.append({
            "pool": number, "oligos": len(members),
            "total_bases": sum(len(m["sequence"]) for m in members),
            "unused_capacity": profile.max_oligos - len(members),
            "members": members,
        })
    return {"feasible": True, "pools": built, "pool_count": len(built),
            "method": ("fewest pools the maximum allows, filled evenly, ordered by "
                       "descending length"),
            "optimality": ("a feasible plan, not a proven cost optimum; compare against "
                           "exhaustive solutions on small cases before claiming otherwise")}


def write_order_files(plan: dict, profile: ProductProfile, outdir: str | Path) -> dict:
    """A vendor-shaped sequence file, a pool membership table and a local estimate report."""
    import csv

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {}

    sequences = out / "order_sequences.csv"
    with sequences.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Pool name", "Sequence"])
        for pool in plan.get("pools", []):
            for member in pool["members"]:
                writer.writerow([f"pool_{pool['pool']}", member["sequence"]])
    paths["order_sequences"] = str(sequences)

    membership = out / "pool_membership.csv"
    fields = ["pool", "design", "module", "quantity", "length", "sequence"]
    with membership.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for pool in plan.get("pools", []):
            for member in pool["members"]:
                writer.writerow({**member, "length": len(member["sequence"])})
    paths["pool_membership"] = str(membership)

    report = out / "order_estimate.json"
    report.write_text(json.dumps({"profile": profile.as_dict(), "plan": {
        "feasible": plan.get("feasible"), "pool_count": plan.get("pool_count"),
        "method": plan.get("method"), "optimality": plan.get("optimality"),
        "reason": plan.get("reason")}}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    paths["estimate_report"] = str(report)
    return paths
