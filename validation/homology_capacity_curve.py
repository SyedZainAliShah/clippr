"""How does achievable library homology grow with library size, and does a better solver help?

**Why this replaces the "84 encodings" headline.** That number is the saturation point of a
sampler under one arbitrary criterion (20-mer disjointness), and it invites two overclaims: that
84 is a proven maximum, and that a 50-member library therefore "cannot avoid homology". The
useful question is neither. It is:

    for N library members, what is the smallest achievable
    maximum exact shared tract between any two of them?

That is a curve, it needs no biological threshold to be meaningful, and it separates three things
a single number confuses: what the current encoder achieves, what a better assignment achieves,
and where the sequence space itself starts to bind.

**Scope, stated because it is easy to overread.** This models the *repeat body only* -- each
member is treated as its nine repeat copies, each drawn from a pool of valid synonymous encodings
of the 31-residue template. Real members also share a scaffold, assembly arms and regulatory
elements, so this is a lower bound on library homology, not the whole landscape. It is also a
combinatorial result about sequence, with no claim about recombination.

**Three strategies, per the "measure before you build" principle.** If the cheap one is already
close to the expensive one, the sophisticated solver is not worth writing.

    A  independent   each member picks encodings without looking at other members
                     -- the coordination the current per-member pass actually has
    B  global greedy each encoding is chosen to minimise conflict with everything placed so far
    C  minimax       B, then local repair swaps aimed at the worst pair specifically

    PYTHONUTF8=1 .venv\\Scripts\\python.exe validation\\homology_capacity_curve.py
"""
from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from clippr.biology import REPEAT_TEMPLATE          # noqa: E402
from clippr.homology import longest_shared          # noqa: E402

K = 20
REPEATS_PER_MEMBER = 9


def encoding_pool(size, seed=0):
    """Valid synonymous encodings of the repeat template, distinct but not disjoint."""
    from Bio import Restriction
    from Bio.Seq import Seq

    from clippr import constants as C
    from clippr.codons import complete_table, table_from_kazusa

    region = REPEAT_TEMPLATE.format(fifth="T", last="N")
    table = complete_table(table_from_kazusa(3055), 1)
    sites = [e for e in (getattr(Restriction, n, None)
                         for n in C.enzymes_for(C.DEFAULT_ENZYME_PROFILE)) if e is not None]
    opts = {aa: [c.replace("U", "T") for c in table[aa]] for aa in set(region)}

    rng = random.Random(seed)
    pool, seen = [], set()
    while len(pool) < size:
        dna = "".join(rng.choice(opts[aa]) for aa in region)
        if dna in seen:
            continue
        seen.add(dna)
        gc = (dna.count("G") + dna.count("C")) / len(dna)
        if not 0.35 <= gc <= 0.65:
            continue
        if any(b * 5 in dna for b in "ACGT"):
            continue
        if str(Seq(dna).translate()) != region:
            continue
        if any(e.search(Seq(dna)) for e in sites):
            continue
        pool.append(dna)
    return pool


def member_seq(encodings):
    return "".join(encodings)


def worst_pair(members):
    """Longest exact shared tract between any two members.

    Exact but quadratic in both member count and sequence length, so it is used for the final
    reported figure only. The search loops use `worst_pair_fast`, which cannot exceed it.
    """
    worst = 0
    for i in range(len(members)):
        for j in range(i + 1, len(members)):
            worst = max(worst, longest_shared(members[i], members[j])[0])
    return worst


def worst_pair_fast(members, lo=K, hi=200):
    """Longest shared tract, by binary search on shared substrings of a given length.

    Hashing every window of length L and asking whether two members share one answers "is the
    worst tract at least L" in linear time, and the answer is monotone in L, so a binary search
    finds the exact value in log steps. The dynamic-programming version is O(len^2) per pair and
    made a 20-member run intractable inside a repair loop.
    """
    def shared_at(length):
        seen: dict[str, int] = {}
        for idx, m in enumerate(members):
            for w in {m[i:i + length] for i in range(len(m) - length + 1)}:
                if seen.get(w, idx) != idx:
                    return True
                seen[w] = idx
        return False

    if not shared_at(lo):
        return lo - 1
    best = lo
    while lo <= hi:
        mid = (lo + hi) // 2
        if shared_at(mid):
            best, lo = mid, mid + 1
        else:
            hi = mid - 1
    return best


def kmers(s):
    return {s[i:i + K] for i in range(len(s) - K + 1)}


def build_independent(n, pool, rng):
    """A: each member samples on its own, blind to the others."""
    return [member_seq(rng.sample(pool, REPEATS_PER_MEMBER)) for _ in range(n)]


def build_global_greedy(n, pool, rng):
    """B: pick each encoding to add the fewest already-used windows."""
    used: set[str] = set()
    members = []
    for _ in range(n):
        chosen = []
        for _ in range(REPEATS_PER_MEMBER):
            sample = rng.sample(pool, min(40, len(pool)))
            best = min(sample, key=lambda e: len(kmers(e) & used))
            chosen.append(best)
            used |= kmers(best)
        members.append(member_seq(chosen))
    return members


def build_minimax(n, pool, rng, rounds=60):
    """C: B, then repair swaps aimed at whichever pair is currently worst."""
    members = build_global_greedy(n, pool, rng)
    parts = [[m[i:i + 93] for i in range(0, len(m), 93)] for m in members]
    for _ in range(rounds):
        cur = worst_pair_fast([member_seq(p) for p in parts])
        # locate the worst pair and try replacing one of its encodings
        target = None
        for i in range(len(parts)):
            for j in range(i + 1, len(parts)):
                if worst_pair_fast([member_seq(parts[i]), member_seq(parts[j])]) == cur:
                    target = (i, j)
                    break
            if target:
                break
        if not target:
            break
        i, slot = target[0], rng.randrange(REPEATS_PER_MEMBER)
        keep = parts[i][slot]
        parts[i][slot] = rng.choice(pool)
        if worst_pair_fast([member_seq(p) for p in parts]) > cur:
            parts[i][slot] = keep          # accept only improvements
    return [member_seq(p) for p in parts]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", type=int, default=400)
    ap.add_argument("--sizes", type=int, nargs="*", default=[6, 10, 20, 30, 50])
    args = ap.parse_args()

    # The pool must exceed the largest library's demand, or reuse is forced and the result
    # measures the pool rather than the sequence space. A 400-encoding pool against 50 members
    # (450 copies needed) produced exactly that artefact: the coordinated strategy scored
    # *worse* than the blind one, because reusing a whole 93-nt encoding creates a 93-nt tract.
    demand = max(args.sizes) * REPEATS_PER_MEMBER
    need = max(args.pool, int(demand * 1.5))
    if need > args.pool:
        print(f"raising pool {args.pool} -> {need}: {max(args.sizes)} members need {demand} "
              f"encodings, and a pool below that forces reuse")
    print(f"building a pool of {need} valid repeat-template encodings...", flush=True)
    pool = encoding_pool(need)
    print(f"pool ready: {len(pool)} encodings of {len(pool[0])} nt, "
          f"largest library needs {demand}\n")

    print("minimum achievable maximum exact shared tract (nt), by library size")
    print(f"{'members':>8}{'A independent':>16}{'B global greedy':>18}{'C minimax':>12}")
    print("-" * 56)
    rows = []
    for n in args.sizes:
        rng = random.Random(1)
        a = worst_pair_fast(build_independent(n, pool, rng))
        b = worst_pair_fast(build_global_greedy(n, pool, random.Random(1)))
        c = worst_pair_fast(build_minimax(n, pool, random.Random(1)))
        rows.append((n, a, b, c))
        print(f"{n:>8}{a:>16}{b:>18}{c:>12}", flush=True)

    print("\nEach member is modelled as its 9 repeat copies only. Real members also share a")
    print("scaffold, assembly arms and regulatory elements, so these are lower bounds on")
    print("library homology rather than the whole landscape. No recombination claim is made.")

    gains = [a - c for _, a, _, c in rows]
    if max(gains) <= 2:
        print("\nA, B and C agree within 2 nt at every size: coordination buys almost nothing")
        print("here, so a sophisticated joint solver is not worth building.")
    else:
        print(f"\nCoordination helps by up to {max(gains)} nt, so a joint solver is worth")
        print("building — the current independent approach leaves real capacity unused.")


if __name__ == "__main__":
    main()
