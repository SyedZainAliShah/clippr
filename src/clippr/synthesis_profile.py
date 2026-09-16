"""One versioned synthesis policy, consumed by every stage that judges a sequence.

**Why this exists.** The GC band lived in two places. `recoding.GC_BAND` was the *validator's*
threshold; `codons.optimize_cds` carried its own `gc_bounds` default and `recode_inventory`
never passed one. Changing either alone changed which candidates were accepted without changing
which were searched for, and a review found the obvious repair -- "widening the band is one
line" -- was false for exactly that reason. A threshold that appears twice is not a constant; it
is a policy that has not been written down.

**What a profile carries, and why each field.**

A number alone cannot be audited. `0.65` does not say whether it came from a vendor page, which
product it applies to, when it was read, or whether breaching it stops an order or merely costs
a point. Every profile therefore records its `source`, its `product`, the `read_on` date, and
an `enforcement` map saying, per rule, which of three things it is:

    hard              a breach is infeasible; the candidate is refused
    target            a breach is reported and delivered; an optimiser steers away from it
    vendor_judgement  not machine-checkable here; recorded so it is not mistaken for checked

That three-way split is the point. Collapsing it is what made an advisory guideline into a
rejection criterion worth about 0.16 CAI.

**Provenance, stated exactly.** `STRICT_LEGACY` reproduces today's shipped behaviour. Its band
matches Twist's published codon-optimisation guidance exactly -- but this repository has never
recorded where its constant came from, so `source` says *unrecorded*, not *Twist*. Matching a
guideline after the fact is not being derived from it, and back-dating a citation would be the
same class of error as every other in this project's history.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, replace

#: How a rule is enforced. The middle value is the one this module exists to make expressible.
HARD = "hard"
TARGET = "target"
VENDOR_JUDGEMENT = "vendor_judgement"
ENFORCEMENT = (HARD, TARGET, VENDOR_JUDGEMENT)


@dataclass(frozen=True)
class SynthesisProfile:
    """A named, versioned, sourced synthesis policy.

    Frozen because a profile that can be mutated after a result is bound to it is not a
    version. Use `narrowed()` or `dataclasses.replace` to derive a new one, which gets a new
    `version`.
    """

    name: str
    source: str
    product: str
    read_on: str | None
    local_gc: tuple[float, float]
    window: int
    max_homopolymer: int
    global_gc: tuple[float, float] | None = None
    enforcement: dict[str, str] = field(default_factory=dict)
    unresolved: tuple[str, ...] = ()
    note: str = ""

    def __post_init__(self) -> None:
        for rule, how in self.enforcement.items():
            if how not in ENFORCEMENT:
                raise ValueError(f"{self.name}: rule {rule!r} has enforcement {how!r}; "
                                 f"expected one of {ENFORCEMENT}")
        lo, hi = self.local_gc
        if not 0.0 <= lo < hi <= 1.0:
            raise ValueError(f"{self.name}: local_gc {self.local_gc} is not a band in [0, 1]")
        if self.global_gc is not None:
            glo, ghi = self.global_gc
            if not 0.0 <= glo < ghi <= 1.0:
                raise ValueError(f"{self.name}: global_gc {self.global_gc} is not a band")
        if self.window < 1:
            raise ValueError(f"{self.name}: window must be positive")

    @property
    def version(self) -> str:
        """Content hash over every field that can change a verdict.

        `name` and `note` are included: two policies that score identically today but are
        documented differently are still different policies to cite.
        """
        payload = json.dumps(asdict(self), sort_keys=True, default=list)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]

    def how(self, rule: str) -> str:
        """Enforcement for one rule. Unlisted rules are hard, because silence must not relax."""
        return self.enforcement.get(rule, HARD)

    def is_hard(self, rule: str) -> bool:
        return self.how(rule) == HARD

    def solver_bounds(self) -> dict:
        """The keyword arguments `codons.optimize_cds` needs to search this policy's space.

        This is the half that was missing: the solver and the validator now take their limits
        from the same object, so they cannot drift apart again.
        """
        return {"gc_bounds": tuple(self.local_gc), "gc_window": self.window,
                "max_homopolymer": self.max_homopolymer}

    def as_dict(self) -> dict:
        return {**asdict(self), "version": self.version}

    def narrowed(self, **changes) -> "SynthesisProfile":
        """A derived profile. Renamed on purpose: a changed policy is not the same policy."""
        changes.setdefault("name", f"{self.name}+derived")
        return replace(self, **changes)


#: **The default, and today's shipped behaviour exactly.** Local 0.35-0.65 over 50 nt,
#: homopolymer <= 4, no global band, every rule hard. Reproduces every number in the current
#: handoff; changing it changes published results.
#:
#: `source` is *unrecorded* deliberately. The band happens to match Twist's published
#: codon-optimisation guidance (global 25-65%, local 35-65% over 50 bp, read 2026-09-16), but
#: this repository has never recorded that as its origin and it would be dishonest to claim it
#: retrospectively. See `docs/gc_band_and_constraint_model.md`.
STRICT_LEGACY = SynthesisProfile(
    name="clippr-strict-legacy",
    source="unrecorded in this repository; see docs/gc_band_and_constraint_model.md",
    product="not product-specific",
    read_on=None,
    local_gc=(0.35, 0.65),
    window=50,
    max_homopolymer=4,
    global_gc=None,
    enforcement={"local_gc": HARD, "homopolymer": HARD, "forbidden_sites": HARD},
    unresolved=(
        "no global GC band is enforced, though every vendor guideline read pairs one with a "
        "local window",
        "the constant's provenance is unknown; it matches Twist's published codon-optimisation "
        "guidance but was not recorded as derived from it",
        "vendor acceptance is multifactorial and scored by proprietary models; no profile here "
        "predicts it",
    ),
    note="today's shipped policy, preserved as the baseline so results stay reproducible",
)

#: Twist's published **codon-optimisation** guidance, read from their page on 2026-09-16:
#: "global GC% of less than 25% or more than 65% and local GC windows (50 bp) of less than 35%
#: or more than 65%".
#:
#: Two differences from `STRICT_LEGACY` that matter: it adds the global band that guidance
#: pairs with the window, and it treats both as **targets** rather than hard limits, because
#: that page frames them as steps taken during optimisation and separates them from its
#: sequence-acceptance rules. It is **not** an IDT oPools rule, and CLIPPR orders from IDT.
TWIST_CODON_GUIDANCE = SynthesisProfile(
    name="twist-codon-optimisation-guidance",
    source="https://www.twistbioscience.com/products/genes/resources",
    product="Twist gene synthesis, codon-optimisation guidance",
    read_on="2026-09-16",
    local_gc=(0.35, 0.65),
    window=50,
    max_homopolymer=4,
    global_gc=(0.25, 0.65),
    # Both bands are targets, because the page frames them together as steps taken during
    # optimisation. Declaring only the local one would leave the global one hard by default —
    # the safe fallback, but not what this source says.
    enforcement={"local_gc": TARGET, "global_gc": TARGET, "homopolymer": TARGET,
                 "forbidden_sites": HARD},
    unresolved=(
        "this is Twist's codon-optimisation guidance, not its sequence-acceptance rules, and "
        "not IDT's",
        "Twist's own homopolymer limit is a separate hard rule (< 14 nt); the value here is a "
        "CLIPPR target carried over from the strict profile",
        "final acceptance is ML-scored by the vendor and is not predicted here",
    ),
    note="opt-in; adds the global band and makes the local band advisory, as the page frames it",
)

#: A deliberately broad **experimental** policy, for measuring what the strict band costs.
#: Not a vendor rule and not for ordering: its bounds come from the 2x2 experiment in
#: `docs/gc_band_and_constraint_model.md`, where it produced a mean CAI of 0.819921 across all
#: 42 modules with every selected substrate passing its own checks.
BROAD_EXPERIMENTAL = SynthesisProfile(
    name="broad-experimental",
    source="CLIPPR experiment, not a vendor rule",
    product="none; measurement only",
    read_on="2026-09-16",
    local_gc=(0.15, 0.85),
    window=50,
    max_homopolymer=4,
    global_gc=None,
    enforcement={"local_gc": HARD, "homopolymer": HARD, "forbidden_sites": HARD},
    unresolved=("no vendor has been shown to accept sequences designed under this policy",),
    note="for measuring the cost of the strict band; never a default, never an order",
)

PROFILES: dict[str, SynthesisProfile] = {
    p.name: p for p in (STRICT_LEGACY, TWIST_CODON_GUIDANCE, BROAD_EXPERIMENTAL)
}

#: What every stage uses unless a caller passes something else. Changing this line changes
#: published results, which is why it is a named constant and not a literal.
DEFAULT_PROFILE = STRICT_LEGACY


def resolve(profile: "SynthesisProfile | str | None") -> SynthesisProfile:
    """Accept a profile, a registered name, or None for the default."""
    if profile is None:
        return DEFAULT_PROFILE
    if isinstance(profile, SynthesisProfile):
        return profile
    try:
        return PROFILES[profile]
    except KeyError:
        raise ValueError(
            f"unknown synthesis profile {profile!r}; known: {sorted(PROFILES)}"
        ) from None
