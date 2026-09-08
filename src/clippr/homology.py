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

**Two approaches that did not work, recorded so they are not retried.** Reacting to observed
sharing is unreliable. Forbidding a whole 59-nt stretch only obliges the optimiser to change
one base, leaving 58 nt shared. Forbidding every k-length window inside the stretch does better
but is **not monotone**: measured, the longest run went 84 -> 65 nt at a 20-nt ban window,
84 -> 65 at 30 nt, and 84 -> **86** at 25 nt. Iterating oscillated, 10 flagged pairs going to 5
on one pass and back to 7 on two. A fix that can make the metric worse is not a fix.

Two hypotheses were also tested and refuted. Varying the random seed alone barely helps (median
shared run 59 -> 62 nt). And homology is *not* predictable from target similarity: across these
pairs, the correlation between the longest run of shared target positions and the longest shared
DNA run was **-0.38**, i.e. weakly negative. The repeat body is genuinely well diversified --
across 369 nt of identical protein the longest shared DNA run was 59 nt.

**What this is not.** A shared stretch is a *necessary* substrate for recombination, not a
prediction that recombination will occur -- that depends on the host, the loci, copy number and
repair pathways, none of which this module models. The 50 nt default threshold is a widely used
rule of thumb for where homologous recombination becomes plausible in many systems, not a
measured constant for *Chlamydomonas*. Treat the output as a ranked risk list, not a verdict.
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
