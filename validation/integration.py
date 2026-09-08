"""End-to-end gate: run the whole pipeline over the fixed corpus.

Passes only if every target designs without an unhandled exception, every result is
internally consistent, and a repeat run under the same seed is identical. This is the
check that the modules compose, which unit tests by construction cannot show.

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\integration.py [--limit N]
"""
from __future__ import annotations

import argparse
import random
import statistics
import sys
import time
import traceback
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from clippr.assembly import reassemble, split_cds
from clippr.design import design_oneshot

ARCH_LEN = {"9S": 9, "14S": 14, "19S": 19}


def corpus(n_per_arch=(100, 50, 50), seed=7) -> list[str]:
    """The fixed corpus, identical to the one compare_ppr.py uses. Never reseed."""
    rng = random.Random(seed)
    out = []
    for (_, length), n in zip(ARCH_LEN.items(), n_per_arch):
        out += ["".join(rng.choice("ACGU") for _ in range(length)) for _ in range(n)]
    return out


def check(result) -> list[str]:
    """Internal consistency of one design. Returns the problems found."""
    bad = []
    cds, protein = result["cds"], result["protein"]
    if str(Seq(cds).translate()) != protein:
        bad.append("translation")
    if len(cds) != 3 * len(protein):
        bad.append("length")
    for cut, oh in zip(result["cuts"], result["junction_overhangs"]):
        if cds[3 * cut - 4:3 * cut] != oh:
            bad.append("locked_overhang")
            break
    frags = result["fragments"]
    if len(frags) != len(result["oligos"]):
        bad.append("fragment_count")
    if reassemble(frags) != cds:
        bad.append("reassembly")
    for _, r in result["oligos"].iterrows():
        if r["payload_5to3"] not in r["oligo_sequence_5to3"]:
            bad.append("payload_not_in_oligo")
            break
    if not 0.0 <= result["fidelity"] <= 1.0:
        bad.append("fidelity_range")
    if result["qc"]["status"] not in ("PASS", "WARNING", "FAIL"):
        bad.append("qc_status")
    return bad


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    targets = corpus()
    if args.limit:
        targets = targets[:args.limit]

    ok = 0
    errors: list[str] = []
    inconsistent: Counter[str] = Counter()
    arch = Counter()
    qc = Counter()
    fid, secs = [], []
    first: dict[str, str] = {}

    for i, rna in enumerate(targets, 1):
        t0 = time.time()
        try:
            r = design_oneshot(rna, seed=42)
        except Exception as exc:
            errors.append(f"{rna}: {type(exc).__name__}: {exc}")
            if len(errors) == 1:
                traceback.print_exc()
            continue
        secs.append(time.time() - t0)
        problems = check(r)
        if problems:
            for p in problems:
                inconsistent[p] += 1
        else:
            ok += 1
        arch[r["architecture"]] += 1
        qc[r["qc"]["status"]] += 1
        fid.append(r["fidelity"])
        first[rna] = r["cds"]
        if i % 25 == 0:
            print(f"  ... {i}/{len(targets)}", flush=True)

    # determinism: re-run a sample under the same seed
    sample = targets[:5] + targets[100:103] + targets[150:153]
    same = sum(1 for t in sample if design_oneshot(t, seed=42)["cds"] == first.get(t))

    n = len(targets)
    print(f"\ntargets designed        : {n}")
    print(f"unhandled exceptions    : {len(errors)}")
    print(f"internally consistent   : {ok}/{n - len(errors)}")
    print(f"architectures           : {dict(arch)}")
    print(f"QC verdicts             : {dict(qc)}")
    if fid:
        print(f"fidelity                : {min(fid):.3f}-{max(fid):.3f} "
              f"(median {statistics.median(fid):.3f})")
    if secs:
        print(f"seconds per design      : {min(secs):.1f}-{max(secs):.1f} "
              f"(median {statistics.median(secs):.1f})")
    print(f"seed-reproducible       : {same}/{len(sample)} resampled targets identical")

    if inconsistent:
        print("\ninconsistencies:", dict(inconsistent))
    for e in errors[:5]:
        print("  " + e)

    passed = not errors and ok == n and same == len(sample)
    print("\n" + ("PASS" if passed else "FAIL"))
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
