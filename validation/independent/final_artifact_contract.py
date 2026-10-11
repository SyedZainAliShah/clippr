"""Independently judge every saved design's ordered oligos, and compare with production.

The previous version of this file initialised `disagreements = []`, never appended to it, and
printed "disagreements with the production validator: 0". That number was the default value of
an untouched list, not the result of a comparison. It is computed here.

Three further defects it carried, each found by an adversarial control rather than by reading:

  * local GC and homopolymer limits were treated as hard regardless of how the effective
    profile enforces them, so a design that production correctly calls `warning` under
    `twist-codon-optimisation-guidance` was reported as rejected;
  * the global GC rule was not applied at all;
  * a missing oligo string was rebuilt with `wrap_fragment` and then described as validation
    of the stored ordered sequence, and a design with two missing fragments was still counted
    `accepted`.

Design status is exclusive: each design is accepted, rejected or not_assessed exactly once.
Fragment counts are reported separately and never mixed into the design denominator. A design
with any unusable fragment cannot be accepted.

Independence means the arithmetic here is its own. Declared constants are shared on purpose;
production *verdicts* are not consulted except to compare against.

    python validation/independent/final_artifact_contract.py \\
        --jsonl work/independent/corpus200_pinned.jsonl \\
        --out work/run/final_artifact_contract.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ppropt import constants as C                                     # noqa: E402
from ppropt.assembly import wrapper_intended_sites                    # noqa: E402
from ppropt.substrates import rc                                      # noqa: E402
from ppropt.synthesis_profile import resolve                          # noqa: E402


def _gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / max(1, len(seq))


def _worst_window(sequence: str, band, width: int):
    worst = None
    for start in range(0, max(1, len(sequence) - width + 1)):
        window = sequence[start:start + width]
        if len(window) < width:
            break
        gc = _gc(window)
        if band[0] <= gc <= band[1]:
            continue
        excess = max(band[0] - gc, gc - band[1])
        if worst is None or excess > worst[0]:
            worst = (excess, round(gc, 4), start)
    return worst


def _longest_run(sequence: str) -> int:
    run = longest = 1
    prev = ""
    for base in sequence:
        run = run + 1 if base == prev else 1
        prev = base
        longest = max(longest, run)
    return longest


def _independent_findings(sequence: str, intended, policy, enzyme_sites):
    """Hard and advisory findings, split the way this policy actually enforces them."""
    hard: list[str] = []
    advisory: list[str] = []

    def record(rule: str, message: str) -> None:
        (hard if policy.is_hard(rule) else advisory).append(message)

    allowed = set(intended)
    for name, site in enzyme_sites.items():
        for strand, pattern in (("forward", site), ("reverse", rc(site))):
            start = 0
            while (at := sequence.find(pattern, start)) >= 0:
                if (name, strand, at) not in allowed:
                    # Forbidden sites are hard in every profile.
                    hard.append(f"unintended {name} {strand} at {at}")
                start = at + 1

    if policy.global_gc is not None:
        overall = _gc(sequence)
        if not policy.global_gc[0] <= overall <= policy.global_gc[1]:
            record("global_gc",
                   f"global GC {overall:.3f} outside {tuple(policy.global_gc)}")

    worst = _worst_window(sequence, policy.local_gc, policy.window)
    if worst:
        record("local_gc",
               f"GC {worst[1]} outside {tuple(policy.local_gc)} at {worst[2]}")

    run = _longest_run(sequence)
    if run > policy.max_homopolymer:
        record("homopolymer", f"homopolymer run of {run} > {policy.max_homopolymer}")
    return hard, advisory


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jsonl", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--profile", default=None)
    ap.add_argument("--limit", type=int, default=None,
                    help="check only the first N designs (controls and smoke runs)")
    ap.add_argument("--expect-designs", type=int, default=None, metavar="N",
                    help=("the panel size this run is required to cover; a shorter run "
                          "fails. Without it, one comparison stood in for two hundred."))
    args = ap.parse_args()

    corpus = Path(args.jsonl)
    if not corpus.is_file():
        print(f"missing required corpus: {corpus}")
        return 2

    policy = resolve(args.profile)
    from Bio.Restriction import RestrictionBatch
    from ppropt.design import final_artifact_report          # production, for comparison
    import pandas as pd

    active = C.enzymes_for(C.DEFAULT_ENZYME_PROFILE)
    enzyme_sites = {n: str(RestrictionBatch([n]).get(n).site) for n in active}

    rows: list[dict] = []
    counts = {"accepted": 0, "rejected": 0, "not_assessed": 0}
    fragment_counts = {"checked": 0, "unusable": 0, "with_hard": 0}
    disagreements: list[dict] = []
    compared = 0

    for line in corpus.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        entry = json.loads(line)
        # The corpus interleaves journal records with design rows.
        if "__run__" in entry or "__run_identity__" in entry:
            continue
        if args.limit is not None and len(rows) >= args.limit:
            break
        target = entry.get("target") or entry.get("target_rna")
        oligos = entry.get("oligos") or []

        fragments, hard, advisory, unusable = [], [], [], []
        for o in oligos if isinstance(oligos, list) else []:
            if not isinstance(o, dict):
                unusable.append("a non-object entry in the oligo list")
                continue
            ident = o.get("order_fragment_id")
            sequence = (o.get("oligo_sequence_5to3") or "").strip()
            payload = (o.get("payload_5to3") or "").strip()
            enzyme = o.get("wrap_enzyme") or "BsaI"
            # Never reconstruct. A rebuilt sequence is not the exported bytes, and
            # checking it would not be a check of what was stored.
            if not sequence:
                unusable.append(f"{ident}: no stored oligo_sequence_5to3")
                continue
            if not payload or payload not in sequence:
                unusable.append(f"{ident}: payload absent or not inside the sequence")
                continue
            intended = wrapper_intended_sites(payload, enzyme)
            f_hard, f_adv = _independent_findings(sequence, intended, policy, enzyme_sites)
            fragment_counts["checked"] += 1
            if f_hard:
                fragment_counts["with_hard"] += 1
            fragments.append({"id": ident, "length_nt": len(sequence),
                              "sha256": hashlib.sha256(sequence.encode()).hexdigest()[:16],
                              "hard": f_hard, "advisory": f_adv})
            hard += [f"{ident}: {m}" for m in f_hard]
            advisory += [f"{ident}: {m}" for m in f_adv]

        fragment_counts["unusable"] += len(unusable)
        if not oligos or unusable:
            status = "not_assessed"
        elif hard:
            status = "rejected"
        else:
            status = "accepted"
        counts[status] += 1

        # The comparison this file previously claimed and never performed.
        production = None
        if fragments or unusable:
            frame = pd.DataFrame([o for o in oligos if isinstance(o, dict)])
            try:
                production = final_artifact_report(frame, profile=policy,
                                                   enzymes=active)["status"]
            except Exception as exc:                 # noqa: BLE001 - recorded, not hidden
                production = f"error: {type(exc).__name__}: {exc}"
            compared += 1
            expected = {"accepted": ("pass", "warning"), "rejected": ("fail",),
                        "not_assessed": ("not_assessed",)}[status]
            if production not in expected:
                disagreements.append({"target": target, "independent": status,
                                      "production": production})

        rows.append({"target": target, "status": status, "hard_count": len(hard),
                     "advisory_count": len(advisory), "unusable": unusable,
                     "production_status": production, "fragments": fragments})

    # The checker's own bytes and the package that produced the designs. The package
    # fingerprint does not cover validation/*.py, so editing this checker must invalidate
    # results the previous version produced rather than inheriting their verdict.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from journal import provenance

    report = {
        "provenance": provenance(__file__, args.jsonl),
        "scope": "wrapped_oligos, as stored; nothing is reconstructed",
        "policy": policy.name,
        "effective_rules": {
            "local_gc": list(policy.local_gc), "window": policy.window,
            "global_gc": list(policy.global_gc) if policy.global_gc else None,
            "max_homopolymer": policy.max_homopolymer,
            "forbidden_sites": sorted(active),
            "hard_rules": sorted(r for r in ("local_gc", "global_gc", "homopolymer")
                                 if policy.is_hard(r)),
        },
        "designs": len(rows),
        **counts,
        "fragments": fragment_counts,
        "comparisons_performed": compared,
        "disagreements_with_production": disagreements,
        "corpus_sha256": hashlib.sha256(corpus.read_bytes()).hexdigest(),
        "means": ("a rejected design failed the declared sequence rules on the sequence that "
                  "would be ordered; it is a correctly detected rejection, not a crash and "
                  "not an accepted design"),
        "rows": rows,
        "expected_designs": args.expect_designs,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"policy {policy.name}: local GC {tuple(policy.local_gc)} over {policy.window} nt, "
          f"homopolymer <= {policy.max_homopolymer}, "
          f"hard rules {report['effective_rules']['hard_rules']}")
    print(f"{len(rows)} designs: {counts['accepted']} accepted, {counts['rejected']} "
          f"rejected, {counts['not_assessed']} not assessed")
    print(f"fragments: {fragment_counts['checked']} checked, "
          f"{fragment_counts['with_hard']} with hard findings, "
          f"{fragment_counts['unusable']} unusable")
    print(f"comparisons performed: {compared}; "
          f"disagreements: {len(disagreements)}")
    for d in disagreements[:5]:
        print(f"  DISAGREEMENT {d['target']}: independent={d['independent']} "
              f"production={d['production']}")
    print(f"wrote {out}")
    if compared == 0:
        print("no comparison was performed; this is not an agreement")
        return 1
    # Agreeing that a design could not be assessed is not a completed validation. The
    # checker exited 0 on it, so a corpus of unreadable rows read as success. A rejection
    # is different: that is a verdict, correctly reached, and must stay 0.
    if counts["not_assessed"]:
        print(f"{counts['not_assessed']} design(s) carried no assessable final artifact; "
              f"required data is missing, so this run is incomplete")
        return 1
    if args.expect_designs is not None and len(rows) != args.expect_designs:
        print(f"expected {args.expect_designs} designs, checked {len(rows)}; a panel "
              f"smaller than the one declared is not the declared panel")
        return 1
    return 1 if disagreements else 0


if __name__ == "__main__":
    raise SystemExit(main())
