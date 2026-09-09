"""DNA shared between library members, and how to break it.

**A risk that only exists at library scale, and that we measured rather than assumed.** A PPR
gene is tandem repeats of one template, and every member of a library carries the same
scaffold. If two members share a long stretch of identical DNA, that stretch is a substrate for
homologous recombination between them once both are in the same cell, and it also makes the
pool harder to synthesise and to sequence unambiguously.

**What was measured** (six 9S designs, 906 nt each, default settings):

    every pair shared 166-269 of 887 distinct 20-mers
    every pair shared an exact stretch of at least 59 nt (median 59, max 84)
    15 of 15 pairs shared 50 nt or more
    varying the random seed per member barely helped: median 59 -> 62 nt

**The cause, located precisely.** Not poor codon diversification. Across 369 nt of *identical
protein* the longest shared DNA run was only 59 nt, so the optimiser is already diversifying
the repeat body well -- the `UniquifyAllKmers` objective forces internal uniqueness and that
incidentally separates members. Every shared stretch instead began at **position 0 in both
members** and decoded to `MQGGNSEEPRKSFDERPER...`: the fixed 23-residue N-terminal scaffold.

That region is protein-identical in every design by construction, sits at the very start where
no internal-repeat pressure applies, and receives the same codons every time because nothing
tells the optimiser that another member already used them.

**So it is fixable, and cheaply.** The scaffold's amino acids are fixed but its codons are
free, so each member can be given its own encoding of it. `diversify_library` does that and
reports the before/after, because a fix of this kind is only worth having if it is measured.
On the five designs above:

    longest shared stretch    84 nt -> 47 nt
    pairs over the 50 nt threshold    10 of 10 -> 0 of 10
    proteins still correct    5 of 5

**But the pairwise view alone would overstate the fix.** Fixing the worst *pair* says nothing
about blocks carried by most of the library, which are the risk that grows with library size, so
`block_profile` reports the library as a whole. On six designs:

    measure                          before   after
    k-mers shared by >1 member         1120     860
    ...present in EVERY member           56       0
    ...present in half or more           480     272
    widest block spans                6 of 6  5 of 6
    pairs over threshold                 15       1
    busiest member appears in             5       1
    members in no flagged pair            0       4

So the scaffold blocks are **eliminated** -- nothing is left in all six members -- but blocks in
five of six survive, and they decode to the *repeat template* (`GAGCTGTTCGACAAGATGCC` is
ELFDKMP..., `CAGAACGGCCGCATTGACGA` is QNGRID...).

**Whether that residue is avoidable was measured rather than assumed, and an earlier version of
this docstring got it wrong.** It claimed the reuse was "unavoidable" without establishing
anything, which is the kind of rationalisation this package is supposed to refuse.
`encoding_capacity` measures the supply instead: **84 mutually 20-mer-disjoint encodings of the
31-residue repeat template were obtained** under the codon table, GC band, enzyme set and
homopolymer limit. An unbiased random search and a greedy search that steers away from used
windows both reach 84, and the random search finds its last new encoding after 25,000 draws and
nothing in the following 775,000.

**That is a saturation point of these searches, not a proven maximum.** The exact maximum is a
set-packing problem and has not been solved. Write "84 encodings were obtained under the stated
constraints", never "the capacity is 84".

**A second measurement then overturned the conclusion drawn from the first, and this is the one
to quote.** Counting *disjoint* encodings asks the wrong question, because members do not need
disjoint encodings -- they need a short worst shared tract.
`validation/homology_capacity_curve.py` measures that directly: for N members, the smallest
achievable maximum exact shared tract, under three assignment strategies.

    members   A independent   B global greedy   C minimax
          6             101                26          26
         10             101                26          26
         20             102                29          26
         30             102                29          29
         50             113                36          35

**There is no wall at 50.** Coordinated assignment holds the worst shared tract to 26-36 nt at
every size tested, and the curve rises gently rather than breaking. The binding constraint was
never the sequence space: it is that assignment is made **per member, independently**, which
costs about 75 nt. So the earlier reading of the 84 figure -- "about 5-fold short, the criterion
cannot be satisfied at 50 members" -- was wrong, and it was wrong because 20-mer *disjointness* is
a far stricter requirement than a short worst tract.

Two consequences. `diversify_library` should coordinate assignment globally rather than
per-member; that is the single largest available improvement here. And the expensive solver is
**not** worth building: C beats B by at most 3 nt, so global greedy is enough.

One trap this exposed: the first run used a 400-encoding pool for a library needing 450, which
forced whole-encoding reuse and made the coordinated strategy score *worse* than the blind one
(194 nt against 186). The script now raises the pool above demand automatically. A capacity
experiment whose pool is smaller than its demand measures the pool.

**20-mer disjointness is our engineering criterion, not a biological threshold**, and the choice
of *k* dominates the answer: capacity measured 8 at k=12, 84 at k=20 and 2578 at k=40, because a
longer window is a *weaker* requirement -- sharing some 12-mer is near-inevitable, sharing a 40-mer
needs 40 consecutive identical bases. Since nothing calibrates any *k* to recombination
probability in this host, no single value should be treated as a safety threshold, and that
sensitivity is itself the argument against pretending one exists.

Finally, these six-member numbers move with library size -- the run leaves one pair at 60 nt where
the five-member run cleared all of them -- so measure the real library rather than quoting from
here.

**Two approaches that did not work, recorded so they are not retried.** Reacting to observed
sharing is unreliable. Forbidding a whole 59-nt stretch only obliges the optimiser to change
one base, leaving 58 nt shared. Forbidding every k-length window inside the stretch does better
but is **not monotone**: measured, the longest run went 84 -> 65 nt at a 20-nt ban window,
84 -> 65 at 30 nt, and 84 -> **86** at 25 nt. Iterating oscillated, 10 flagged pairs going to 5
on one pass and back to 7 on two. A fix that can make the metric worse is not a fix.

Two hypotheses were also tested and not supported. Varying the random seed alone barely helps
(median shared run 59 -> 62 nt). And **no positive relationship was found between target
similarity and DNA homology** -- the rank correlation between the longest run of shared target
positions and the longest shared DNA run came out at -0.38. That rests on 15 pairs from a
five-member set, so it is enough to say the expected relationship did not appear and not enough
to claim a negative one. The repeat body is genuinely well diversified either way: across 369 nt
of identical protein the longest shared DNA run was 59 nt.

**What this is not, stated plainly because the temptation runs the other way.** This module
measures *sequence identity between designs*. It does not establish that any of it is dangerous.

  * **There is no evidence here for a 50 nt danger threshold in *Chlamydomonas*.** The default is
    a rule of thumb from general practice, not a measured constant for this host, and nothing in
    this repository derives one. Never write "59 nt of homology is dangerous in Chlamydomonas".
  * A shared stretch is a **necessary substrate** for homologous recombination, never a
    prediction that recombination will occur. That depends on the host's repair pathways, the
    integration loci, copy number and expression context, none of which this module models.
  * **The risk depends on the physical library architecture, which is not ours to assume.** One
    construct per strain is a different situation from many constructs entering the same nuclear
    genome, which is different again from a pooled DNA mixture handled before transformation. The
    measurement is a design precaution across all three; it is a cellular claim in none of them.

The defensible sentence is: *independently designed library members acquired substantial
unintended DNA identity through a fixed shared scaffold, and synonymous redesign reduced the
longest shared tract from 84 to 47 nt.* Everything beyond that needs the bench.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

#: Shared-stretch length at which a pair is flagged. A rule of thumb, not a measured constant.
HR_THRESHOLD_NT = 50

#: k for the shared-k-mer count, matching the uniquifier's default so the two are comparable.
KMER = 20


def longest_shared(a: str, b: str) -> tuple[int, int, int]:
    """Longest exact common substring: (length, start in `a`, start in `b`).

    Dynamic programming over two rows. Coding sequences here are ~900 nt, so the O(len(a) *
    len(b)) cost is a few million operations per pair -- fine for a 50-member library, and
    worth the exactness over a k-mer approximation that could not report the true run length.
    """
    if not a or not b:
        return 0, 0, 0
    prev = [0] * (len(b) + 1)
    best = end_a = end_b = 0
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        ai = a[i - 1]
        for j in range(1, len(b) + 1):
            if ai == b[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best, end_a, end_b = cur[j], i, j
        prev = cur
    return best, end_a - best, end_b - best


def shared_kmers(a: str, b: str, k: int = KMER) -> int:
    """How many distinct k-mers the two sequences have in common."""
    if len(a) < k or len(b) < k:
        return 0
    ka = {a[i:i + k] for i in range(len(a) - k + 1)}
    return sum(1 for i in range(len(b) - k + 1) if b[i:i + k] in ka)


@dataclass(frozen=True)
class SharedStretch:
    """The longest identical DNA between two library members."""

    a: str
    b: str
    length: int
    start_a: int
    start_b: int
    sequence: str
    n_shared_kmers: int

    @property
    def flag(self) -> str:
        return "REVIEW" if self.length >= HR_THRESHOLD_NT else "ok"

    @property
    def at_sequence_start(self) -> bool:
        """Both copies at the very start means the shared scaffold, not a chance match."""
        return self.start_a == 0 and self.start_b == 0


def assess(cds_by_name: dict[str, str], *, k: int = KMER) -> list[SharedStretch]:
    """Every pair of members, longest shared stretch first."""
    names = list(cds_by_name)
    out = []
    for a, b in combinations(names, 2):
        sa, sb = cds_by_name[a], cds_by_name[b]
        n, ia, ib = longest_shared(sa, sb)
        out.append(SharedStretch(a=a, b=b, length=n, start_a=ia, start_b=ib,
                                 sequence=sa[ia:ia + n],
                                 n_shared_kmers=shared_kmers(sa, sb, k)))
    return sorted(out, key=lambda s: -s.length)


def repeat_blocks(cds_by_name: dict[str, str], *, k: int = KMER,
                  min_members: int = 2) -> list[tuple[str, int]]:
    """k-mers carried by several members at once, most widely shared first.

    The pairwise view can be misleading. Reducing the worst *pair* says nothing about whether a
    handful of blocks are still present in most of the library, and a k-mer in twenty members is
    a different problem from one in two: it is a repeat that scales with library size rather
    than a coincidence between two designs.
    """
    counts: dict[str, int] = {}
    for seq in cds_by_name.values():
        for km in {seq[i:i + k] for i in range(len(seq) - k + 1)}:
            counts[km] = counts.get(km, 0) + 1
    return sorted(((km, n) for km, n in counts.items() if n >= min_members),
                  key=lambda kv: -kv[1])


def block_profile(cds_by_name: dict[str, str], *, k: int = KMER) -> dict:
    """How much of the library is shared structure rather than pairwise coincidence."""
    n_members = len(cds_by_name)
    blocks = repeat_blocks(cds_by_name, k=k, min_members=2)
    pairs = assess(cds_by_name, k=k)
    over = [p for p in pairs if p.length >= HR_THRESHOLD_NT]
    degree: dict[str, int] = {n: 0 for n in cds_by_name}
    for p in over:
        degree[p.a] += 1
        degree[p.b] += 1
    return {
        "k": k,
        "n_members": n_members,
        "shared_kmers": len(blocks),
        "in_all_members": sum(1 for _, n in blocks if n == n_members),
        "in_half_or_more": sum(1 for _, n in blocks if n * 2 >= n_members),
        "max_multiplicity": blocks[0][1] if blocks else 0,
        "edges_over_threshold": len(over),
        "max_degree": max(degree.values(), default=0),
        "isolated_members": sum(1 for v in degree.values() if v == 0),
    }


def report(cds_by_name: dict[str, str], *, k: int = KMER, limit: int = 10) -> str:
    """A readable summary with the caveat attached to the numbers."""
    pairs = assess(cds_by_name, k=k)
    if not pairs:
        return "fewer than two members; nothing to compare"

    flagged = [p for p in pairs if p.flag == "REVIEW"]
    lines = [f"{len(cds_by_name)} members, {len(pairs)} pairs, "
             f"{len(flagged)} sharing {HR_THRESHOLD_NT} nt or more",
             "",
             f"{'member A':<16}{'member B':<16}{'shared nt':>10}{f'{k}-mers':>9}"
             f"{'at start':>10}  flag",
             "-" * 70]
    for p in pairs[:limit]:
        lines.append(f"{p.a:<16}{p.b:<16}{p.length:>10}{p.n_shared_kmers:>9}"
                     f"{str(p.at_sequence_start):>10}  {p.flag}")
    if len(pairs) > limit:
        lines.append(f"... {len(pairs) - limit} further pairs, all shorter")

    if flagged and all(p.at_sequence_start for p in flagged):
        lines += ["", "Every flagged stretch starts at position 0 in both members, which is "
                      "the fixed scaffold rather than a chance match. `diversify_library` "
                      "can break it: the scaffold's residues are fixed but its codons are "
                      "not."]

    # The pairwise view alone can mislead: fixing the worst pair says nothing about blocks
    # carried by most of the library, which are the risk that scales with library size.
    prof = block_profile(cds_by_name, k=k)
    lines += ["", f"library-wide structure at k={k}:",
              f"  {prof['shared_kmers']} k-mers appear in more than one member",
              f"  {prof['in_all_members']} appear in every member, "
              f"{prof['in_half_or_more']} in half or more",
              f"  widest block spans {prof['max_multiplicity']} of "
              f"{prof['n_members']} members",
              f"  {prof['edges_over_threshold']} pairs over threshold; busiest member is in "
              f"{prof['max_degree']} of them; {prof['isolated_members']} members in none"]
    lines += ["", f"A shared stretch is a necessary substrate for homologous recombination, "
                  f"not a prediction that it will happen. The {HR_THRESHOLD_NT} nt threshold "
                  f"is a rule of thumb, not a measured constant for this host."]
    return "\n".join(lines)


def scaffold_encoding(member: int, protein_prefix: str, codon_table, *,
                      genetic_code: int = 1, enzymes=None,
                      gc_bounds: tuple[float, float] = (0.35, 0.65),
                      tries: int = 200) -> str:
    """A distinct, valid synonymous encoding of a fixed protein prefix for one member.

    Deterministic in `member`, so a library is reproducible and a member keeps its encoding
    across reruns. Rejects any encoding that mistranslates, carries a blacklisted enzyme site
    or falls outside the GC band, and raises rather than returning something invalid.
    """
    import random

    from Bio import Restriction
    from Bio.Seq import Seq

    from . import constants as C
    from .codons import complete_table

    table = complete_table(codon_table, genetic_code)
    names = C.enzymes_for(C.DEFAULT_ENZYME_PROFILE if enzymes is None else enzymes)
    sites = [e for e in (getattr(Restriction, n, None) for n in names) if e is not None]
    options = {aa: [c.replace("U", "T") for c in table[aa]] for aa in set(protein_prefix)}

    rng = random.Random(member)
    for _ in range(tries):
        dna = "".join(rng.choice(options[aa]) for aa in protein_prefix)
        gc = (dna.count("G") + dna.count("C")) / len(dna)
        if not gc_bounds[0] <= gc <= gc_bounds[1]:
            continue
        if str(Seq(dna).translate()) != protein_prefix:
            continue
        if any(e.search(Seq(dna)) for e in sites):
            continue
        return dna
    raise ValueError(
        f"no valid synonymous encoding of the {len(protein_prefix)}-residue prefix found for "
        f"member {member} in {tries} tries under GC {gc_bounds} and sites {names}")


def encoding_capacity(protein_region: str, codon_table, *, k: int = KMER,
                      genetic_code: int = 1, enzymes=None,
                      gc_bounds: tuple[float, float] = (0.35, 0.65),
                      max_homopolymer: int = 4, draws: int = 200_000,
                      seed: int = 0) -> dict:
    """How many synonymous encodings of a protein region share no k-mer with each other?

    **This is the question that decides whether residual library homology is a bug or a bound.**
    A gene containing *r* copies of a repeat needs *r* mutually disjoint encodings for itself
    alone, so a library of *N* members needs roughly *r x N*. If supply exceeds demand, residual
    sharing means the encoder is failing to use available capacity; if supply falls short, the
    sharing is forced and no encoder can remove it.

    Searched two ways, because a plain random draw is a weak search whose acceptance curve
    flattens for efficiency reasons as much as for exhaustion: unbiased sampling, and a greedy
    pass that refuses codons which would recreate an already-used k-mer.

    **What the number is and is not.** Two searches agreeing is evidence that both have
    saturated, not a certificate that no larger set exists -- the exact maximum is a set-packing
    problem this does not solve. Report it as "N encodings were obtained under these
    constraints", never as "the capacity is N". It also depends heavily on `k`, and in the
    counter-intuitive direction: a longer window is a *weaker* disjointness requirement, so
    capacity rises with `k` (measured 8 at k=12, 84 at k=20, 2578 at k=40). Since no `k` here is
    calibrated to any biological process, none of these is a safety threshold.
    """
    import random

    from Bio import Restriction
    from Bio.Seq import Seq

    from . import constants as C
    from .codons import complete_table

    table = complete_table(codon_table, genetic_code)
    names = C.enzymes_for(C.DEFAULT_ENZYME_PROFILE if enzymes is None else enzymes)
    sites = [e for e in (getattr(Restriction, n, None) for n in names) if e is not None]
    options = {aa: [c.replace("U", "T") for c in table[aa]] for aa in set(protein_region)}

    def ok(dna: str) -> bool:
        gc = (dna.count("G") + dna.count("C")) / len(dna)
        if not gc_bounds[0] <= gc <= gc_bounds[1]:
            return False
        if any(b * (max_homopolymer + 1) in dna for b in "ACGT"):
            return False
        if str(Seq(dna).translate()) != protein_region:
            return False
        return not any(e.search(Seq(dna)) for e in sites)

    def windows(dna: str) -> set[str]:
        return {dna[i:i + k] for i in range(len(dna) - k + 1)}

    rng = random.Random(seed)
    found = last_hit = 0
    used: set[str] = set()
    for i in range(1, draws + 1):
        dna = "".join(rng.choice(options[aa]) for aa in protein_region)
        if ok(dna):
            w = windows(dna)
            if not (w & used):
                found += 1
                used |= w
                last_hit = i

    rng2 = random.Random(seed + 1)
    greedy, used2 = 0, set()
    for _ in range(max(1, draws // 4)):
        dna = ""
        for aa in protein_region:
            choices = options[aa][:]
            rng2.shuffle(choices)
            dna += next((c for c in choices
                         if len(dna + c) < k or (dna + c)[-k:] not in used2), choices[0])
        if ok(dna):
            w = windows(dna)
            if not (w & used2):
                greedy += 1
                used2 |= w

    return {
        "region_aa": len(protein_region),
        "k": k,
        "random_search": found,
        "greedy_search": greedy,
        "capacity": max(found, greedy),
        "searches_agree": found == greedy,
        "last_random_hit_at_draw": last_hit,
        "draws": draws,
    }


def _least_conflicting_scaffold(placed: list[str], protein_prefix: str, codon_table, *,
                                genetic_code: int, enzymes, candidates: int = 16,
                                k: int = KMER) -> str:
    """The candidate scaffold encoding sharing fewest k-mers with what is already placed.

    **Not used, and kept only to record a negative result.** The capacity curve
    (`validation/homology_capacity_curve.py`) showed coordinated assignment beating independent
    assignment by roughly 75 nt, so wiring this into `diversify_library` looked like the largest
    improvement available. Measured on the real pipeline over six targets, it made things
    **worse**: worst shared tract 60 -> 65 nt and flagged pairs 1 -> 6 of 15, against plain
    index assignment.

    The reason the prediction did not transfer is that the curve modelled a world where the
    assignment controls every repeat's encoding. Here it controls only the 69-nt scaffold of a
    906-nt sequence, and changing the locked prefix changes DNA Chisel's whole trajectory for
    the other 837 nt. Optimising the small part perturbs the large part more than it gains --
    the same effect that made reactive k-mer banning non-monotone.

    Left in place because "we tried coordinating and it was worse" is worth more to the next
    person than a silently absent idea. Any future attempt should control the body's encoding
    too, not just the prefix.
    """
    used: set[str] = set()
    for cds in placed:
        used |= {cds[i:i + k] for i in range(len(cds) - k + 1)}

    best, fewest = None, None
    for member in range(candidates):
        enc = scaffold_encoding(member, protein_prefix, codon_table,
                                genetic_code=genetic_code, enzymes=enzymes)
        clash = len({enc[i:i + k] for i in range(len(enc) - k + 1)} & used)
        if fewest is None or clash < fewest:
            best, fewest = enc, clash
        if clash == 0:
            break                      # nothing shared; no candidate can beat that
    return best


def diversify_library(targets, *, threshold: int = HR_THRESHOLD_NT, codon_table=None,
                      genetic_code: int = 1, retries: int = 3, on_progress=None,
                      **kwargs) -> dict:
    """Give every member its own encoding of the shared scaffold, then measure the effect.

    **Constructive, then verified.** Each member is assigned a distinct synonymous encoding of
    the fixed N-terminal scaffold and that encoding is locked. This is deterministic in the
    member's index, needs no search and cannot oscillate the way reactive banning did (see the
    module docstring for what failed). Measured over 50 encodings of the 23-residue scaffold:
    1225 pairs, longest shared run **median 9 nt, maximum 27 nt**, none above 30.

    Locking the scaffold changes the optimiser's whole trajectory, though, so the repeat body
    can coincidentally converge — measured once, nine pairs of ten improved while one went from
    84 nt to 181 nt. Each member is therefore re-rolled up to `retries` times while it still
    shares `threshold` or more with an already-finalised member, and a diversified member is
    **accepted only if it is no worse than the baseline it replaces**. That makes the whole
    procedure monotone by construction.

    Returns the designs plus a measured before/after, because an unmeasured fix is worth
    nothing. `kwargs` are passed to `design_oneshot`; cost is up to `2 + retries` designs per
    member.
    """
    from .biology import N_TERMINAL
    from .design import design_oneshot

    targets = [str(t).strip().upper() for t in targets if str(t).strip()]
    base_seed = kwargs.pop("seed", 42)
    baseline, current, encodings = {}, {}, {}
    for i, t in enumerate(targets, 1):
        base = design_oneshot(t, codon_table=codon_table, seed=base_seed, **kwargs)
        baseline[t] = base["cds"]
        # Assignment is by member index, deliberately. Choosing the scaffold that conflicts
        # least with the placed library sounds strictly better and was measured to be worse --
        # see `_least_conflicting_scaffold`, which is kept only to document that result.
        prefix = scaffold_encoding(
            i - 1, N_TERMINAL, base["codon_table"], genetic_code=base["genetic_code"],
            enzymes=base["enzyme_profile_effective"])
        encodings[t] = prefix
        # Locking the scaffold fixes the scaffold but leaves the repeat body to the optimiser,
        # and changing the locked prefix changes its whole search trajectory. Measured, that
        # occasionally makes one pair much worse: nine pairs of ten improved while one went
        # from 84 nt shared to 181 nt. So build the member, then keep re-rolling its seed
        # while it still shares too much with an already-finalised member, and accept only an
        # outcome no worse than the baseline. That makes the fix monotone by construction
        # rather than by luck, which is what the reactive k-mer ban failed to be.
        best_cds, best_worst = None, None
        for attempt in range(retries + 1):
            cand = design_oneshot(t, codon_table=codon_table, lock_prefix=prefix,
                                  seed=base_seed + 1000 * attempt + i, **kwargs)["cds"]
            worst = max((longest_shared(cand, done)[0] for done in current.values()),
                        default=0)
            if best_worst is None or worst < best_worst:
                best_cds, best_worst = cand, worst
            if best_worst < threshold:
                break
        baseline_worst = max((longest_shared(baseline[t], baseline[o])[0]
                              for o in current), default=0)
        # never accept a diversified member that is worse than what it replaced
        current[t] = best_cds if best_worst <= max(baseline_worst, threshold - 1) \
            else baseline[t]
        if on_progress:
            on_progress(i, len(targets), t)

    before, after = assess(baseline), assess(current)
    return {
        "cds": current,
        "baseline_cds": baseline,
        "scaffold_encodings": encodings,
        "before": before,
        "after": after,
        "max_before": max((p.length for p in before), default=0),
        "max_after": max((p.length for p in after), default=0),
        "flagged_before": sum(1 for p in before if p.length >= threshold),
        "flagged_after": sum(1 for p in after if p.length >= threshold),
        "n_pairs": len(before),
        "threshold": threshold,
    }
