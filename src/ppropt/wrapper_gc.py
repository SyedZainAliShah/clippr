"""The GC constraint evaluated on the ordered fragment, not on the coding span.

A synthesis order is not the coding sequence. It is the coding sequence with an entry
overhang and an enzyme wrapper fixed to each end, and a GC window that straddles that
boundary exists only once the wrapper is attached. Constraining the coding span alone left
those windows unmeasured: on the pinned 200-design corpus every design satisfied the coding
constraint and every design broke the same rule on the sequence that would actually be sent.

This module is the constraint that closes that gap. It lived in
`validation/experiments/wrapper_aware_pilot.py` while it was an experiment; the pilot now
imports it from here, so the shipped route and the experiment cannot disagree about what is
being measured.

Scope, stated because it is easy to over-read: this constrains GC on the ordered string. It
says nothing about secondary structure, and nothing about whether a vendor will accept the
fragment -- the band it enforces is this project's own, and is stricter than any vendor rule
we have been able to verify first-hand.
"""
from __future__ import annotations


def _gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / max(1, len(seq))


def frozen_geometry(baseline: dict) -> list[dict]:
    """Each fragment's coding slice and the fixed DNA that flanks it when ordered.

    Derived from the baseline's own oligo table rather than rebuilt from the wrapper rules,
    so the pilot cannot disagree with production about where the payload sits. The slice
    must occur exactly once in the ordered string or the mapping is ambiguous and we refuse.
    """
    cds = baseline["cds"]
    out = []
    for _, row in baseline["oligos"].iterrows():
        full = str(row["oligo_sequence_5to3"])
        a, b = int(row["cds_start"]), int(row["cds_end"])
        chunk = cds[a:b]
        if full.count(chunk) != 1:
            raise ValueError(
                f"{row['fragment_id']}: its coding slice occurs {full.count(chunk)} times "
                f"in the ordered sequence, so the flanks cannot be identified")
        i = full.index(chunk)
        out.append({"fragment_id": str(row["fragment_id"]), "cds_start": a, "cds_end": b,
                    "left": full[:i], "right": full[i + len(chunk):]})
    return out


def make_constraint(fragments, policy):
    """A DNA Chisel constraint that scores the ordered fragment, not the coding span.

    Built here rather than in `src/` because this is an experiment. Promoting it would be a
    separate decision with its own evidence, and one passing example would not justify it.
    """
    import dnachisel as dc

    lo, hi = policy.local_gc
    window = policy.window

    class EnforceGCInOrderedFragments(dc.Specification):
        """Every GC window of every ordered fragment, flanks included."""

        def __init__(self, boost=1.0):
            self.boost = boost

        def evaluate(self, problem):
            seq = str(problem.sequence)
            breaches, worst = [], None
            for f in fragments:
                a, b = f["cds_start"], f["cds_end"]
                full = f["left"] + seq[a:b] + f["right"]
                offset = a - len(f["left"])        # full[i] is CDS base offset + i
                for i in range(0, max(1, len(full) - window + 1)):
                    gc = _gc(full[i:i + window])
                    if lo <= gc <= hi:
                        continue
                    if worst is None or abs(gc - (lo + hi) / 2) > abs(worst - (lo + hi) / 2):
                        worst = gc
                    # Clipped to the mutable coding span: the flanks cannot be changed, so
                    # pointing the solver at them would ask for the impossible.
                    start, end = max(a, offset + i), min(b, offset + i + window)
                    if end > start:
                        breaches.append((start, end))
            locations = [dc.Location(s, e) for s, e in sorted(set(breaches))]
            return dc.SpecEvaluation(
                self, problem, score=-len(locations), locations=locations,
                message=(f"{len(locations)} ordered-fragment window(s) outside "
                         # `if worst` dropped the figure at gc == 0.0, which is the most
                         # severe breach there is. The sibling fold constraint got this right.
                         f"({lo}, {hi})" + (f", worst {worst:.3f}"
                                            if worst is not None else "")))

        def localized(self, location, problem=None, with_righthand=True):
            # Evaluated globally: one codon change moves a window in a fragment whose
            # flanks sit outside this location, so a localised view would miss it.
            #
            # `with_righthand` is passed by some of DNA Chisel's resolution paths. Omitting
            # it from the signature raised TypeError inside the solver, which reported
            # constraints_ok=False and handed back the unoptimised sequence -- read as
            # eight targets where the band was unsatisfiable. It was this.
            return self

        def label_parameters(self):
            return [("window", str(window)), ("gc", f"{lo}-{hi}")]

    return EnforceGCInOrderedFragments()
