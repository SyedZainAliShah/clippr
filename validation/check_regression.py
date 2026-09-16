"""Drift gate: our own 200 designs against recorded hashes. Needs no external data.

**Why this exists.** The parity checks in this directory (`compare_assembly.py` and its
siblings) compare our output against a 200-design reference corpus. That corpus is 38 MB of
third-party program output and is not redistributed, so on a fresh checkout those checks
cannot run at all -- which left every headline number in the README unverifiable by anyone
who cloned the repository.

This check closes that gap for the part that can honestly be closed. It does **not** prove
parity with the reference; it proves that *this* package still produces exactly what it
produced when the fixture was recorded. That is a different and weaker claim, and it is
stated as such everywhere it appears.

What makes it work without shipping any sequence data is that the 200 target RNAs are
themselves generated from a fixed seed (`integration.corpus`), so the inputs are
reproducible from code. Only output hashes are stored -- a few tens of KB.

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\check_regression.py
    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\check_regression.py --write

`--write` re-records the fixture. Only do that when an output change is intended, and say
in the commit message which designs changed and why.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from clippr.design import design_oneshot          # noqa: E402
from integration import corpus                    # noqa: E402

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "design_hashes.json"

#: The seed every recorded design was produced under. Changing it invalidates the fixture.
SEED = 42


def digest(*parts: str) -> str:
    h = hashlib.sha256()
    for p in parts:
        h.update(p.encode())
        h.update(b"\x1f")                          # field separator, so fields cannot merge
    return h.hexdigest()[:16]


def fingerprint(result: dict) -> dict:
    """The outputs worth pinning, each hashed separately so a failure says which changed.

    Separate hashes rather than one, because "the oligos changed" and "the codon choice
    changed" are different bugs and a single combined digest would not distinguish them.
    """
    oligos = result["oligos"]
    return {
        "architecture": result["architecture"],
        "protein": digest(result["protein"]),
        "cds": digest(result["cds"]),
        "cuts": digest(",".join(str(c) for c in result["cuts"])),
        "overhangs": digest(",".join(result["junction_overhangs"])),
        "oligos": digest(*oligos["oligo_sequence_5to3"].tolist()),
        "n_fragments": len(oligos),
        "qc": result["qc"]["status"],
        "fidelity": f"{result['fidelity']:.6f}",
    }


def build(limit: int = 0) -> dict:
    targets = corpus()
    if limit:
        targets = targets[:limit]
    designs = {}
    for i, rna in enumerate(targets, 1):
        designs[rna] = fingerprint(design_oneshot(rna, seed=SEED))
        if i % 25 == 0:
            print(f"  ... {i}/{len(targets)}", flush=True)
    return {
        "seed": SEED,
        "n": len(designs),
        "note": "Self-regression only. Does NOT establish parity with any other "
                "implementation; see the module docstring.",
        "designs": designs,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true", help="re-record the fixture")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    if args.write:
        data = build(args.limit)
        FIXTURE.parent.mkdir(parents=True, exist_ok=True)
        FIXTURE.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
        print(f"\nwrote {FIXTURE.relative_to(ROOT)}  "
              f"({data['n']} designs, {FIXTURE.stat().st_size:,} bytes)")
        return

    if not FIXTURE.exists():
        sys.exit(f"no fixture at {FIXTURE}; run with --write to record one")

    want = json.loads(FIXTURE.read_text(encoding="utf-8"))
    if want["seed"] != SEED:
        sys.exit(f"fixture was recorded under seed {want['seed']}, this script uses {SEED}")

    targets = corpus()
    if args.limit:
        targets = targets[:args.limit]

    checked = drifted = 0
    missing: list[str] = []
    changes: dict[str, int] = {}
    detail: list[str] = []

    for i, rna in enumerate(targets, 1):
        expected = want["designs"].get(rna)
        if expected is None:
            missing.append(rna)
            continue
        got = fingerprint(design_oneshot(rna, seed=SEED))
        checked += 1
        differing = [k for k in expected if str(expected[k]) != str(got.get(k))]
        if differing:
            drifted += 1
            for k in differing:
                changes[k] = changes.get(k, 0) + 1
            if len(detail) < 5:
                k = differing[0]
                detail.append(f"{rna}: {k} {expected[k]} -> {got.get(k)}"
                              + (f" (+{len(differing) - 1} more fields)"
                                 if len(differing) > 1 else ""))
        if i % 25 == 0:
            print(f"  ... {i}/{len(targets)}", flush=True)

    print(f"\ndesigns checked          : {checked}")
    print(f"identical to fixture     : {checked - drifted}/{checked}")
    if missing:
        print(f"absent from fixture      : {len(missing)} (re-record with --write)")
    if changes:
        print(f"fields that drifted      : {changes}")
        print("\nfirst differences:")
        for d in detail:
            print("  " + d)

    ok = drifted == 0 and not missing
    print("\nProves reproducibility, not parity with any other implementation.")
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
