"""One call from a target RNA to an orderable, QC'd assembly plan.

The pipeline is `ppr` -> `arelf` -> `overhangs` -> `codons` -> `assembly` -> `qc` ->
`export`, and the ordering is forced by a dependency that is easy to get backwards: the
Golden Gate overhangs must be chosen *before* the coding sequence is optimised, because
they are then locked into it. Optimise first and the codon optimiser is free to change the
very bases the junctions depend on.

Cut selection is bounded rather than exhaustive. There are thousands of legal ways to
split a PPR into fragments, but predicted fidelity reaches its ceiling of 1.000 easily, so
the search tries the most balanced splits first and stops at the first proven optimum.
Balanced splits are tried first because fragment cost and synthesis risk both rise with
the longest fragment, not the average.
"""
from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from . import constants as C
from .arelf import (achievable_overhangs, balanced_cuts, realization_count,
                    safe_overhangs)
from .audit import DesignAudit, OverhangDecision
from .assembly import build_oligos, split_cds
from .codons import (optimize_cds, table_from_cds_fasta, table_from_csv,
                     table_from_genome, table_from_kazusa)
from .export import (opool_quote, write_fasta, write_gene_fasta, write_genbank,
                     write_oligo_csv)
from .overhangs import (best_set, enumerate_candidates, fidelity_components,
                        reaction_overhangs, set_fidelity)
from .ppr import describe
from .qc import synthesis_qc

#: Host and fragment-length choices live in `policy.py`, not here -- they are decisions
#: this project made, not facts about GRASP, and the split exists so that stays visible.
#:
#: `ORGANISMS` holds only the verified nuclear entry: taxid 3055 returns a
#: 21-amino-acid Chlamydomonas nuclear table with Leu CUG at 0.73. No Kazusa taxid for the
#: Chlamydomonas *chloroplast* table has been confirmed, and guessing one would silently
#: produce wrong DNA -- the whole failure this package is built to avoid. For chloroplast
#: work, pass `codon_table=` a CDS FASTA or CSV and `genetic_code=11`.
ORGANISMS = C.ORGANISMS
MAX_FRAGMENT_AA = C.MAX_FRAGMENT_AA
MIN_FRAGMENT_AA = C.MIN_FRAGMENT_AA

#: How many candidate splits to score for fidelity before settling.
SPLIT_BUDGET: int = 40

#: How many of those splits to run through codon optimisation before giving up. Each
#: attempt costs a full optimisation. The usual reason for a second attempt -- an overhang
#: that forces a blacklisted site -- is now prevented up front by `arelf.safe_overhangs`,
#: so this is a backstop for genuinely hard proteins, not the common path.
OPTIMISE_ATTEMPTS: int = 3


def _resolve_codon_table(codon_table, organism: str, genetic_code: int):
    """Accept a dict, a path to a CSV or CDS FASTA, or fall back to the named organism."""
    if isinstance(codon_table, Mapping):
        return dict(codon_table)
    if codon_table is not None:
        path = Path(codon_table)
        if not path.exists():
            raise FileNotFoundError(f"codon table not found: {path}")
        if path.suffix.lower() == ".csv":
            return table_from_csv(path)
        return table_from_cds_fasta(path, genetic_code=genetic_code)
    if organism not in ORGANISMS:
        raise ValueError(
            f"unknown organism {organism!r}; known: {sorted(ORGANISMS)}. For any other "
            f"host, pass codon_table= a CDS FASTA or CSV together with the right "
            f"genetic_code."
        )
    source, _ = ORGANISMS[organism]
    # An int is a Kazusa lookup; a str is a genome accession whose annotated genes are
    # counted. Organelles take the second path because Kazusa has no entry for them.
    return (table_from_genome(source, genetic_code) if isinstance(source, str)
            else table_from_kazusa(source))


def _plan_fragments(protein: str, n_fragments: int, destination, matrix: str,
                    enzyme_profile=C.DEFAULT_ENZYME_PROFILE):
    """Pick cut positions and their overhangs together, best balance first.

    Returns (cuts, overhangs, fidelity). Stops at the first split whose reaction reaches
    the metric's ceiling; otherwise returns the best seen within `SPLIT_BUDGET`.
    """
    splits = balanced_cuts(protein, n_fragments, MIN_FRAGMENT_AA, MAX_FRAGMENT_AA,
                           limit=SPLIT_BUDGET)
    if not splits:
        raise ValueError(
            f"no legal split of {len(protein)} residues into {n_fragments} fragments of "
            f"{MIN_FRAGMENT_AA}-{MAX_FRAGMENT_AA} aa")

    fixed = reaction_overhangs([], "level0") if destination is None else [
        destination[0], _rc(destination[1])]

    # Adding overhangs to a reaction can only take probability away, so the destination
    # pair scored alone is the ceiling any split can reach -- and it is reachable, when
    # the junctions cross-talk with nothing. Stopping there rather than at 1.000 matters:
    # the level 0 pair alone scores about 0.83, so a search waiting for 1.000 would never
    # exit early and would score every candidate split for nothing.
    ceiling = set_fidelity(fixed, matrix)

    out = []
    for cuts in splits:
        # safe_overhangs, not achievable_overhangs: an overhang that forces a blacklisted
        # site hands the codon optimiser an unsatisfiable problem, and no retry recovers.
        try:
            options = enumerate_candidates(
                [",".join(safe_overhangs(protein, c, enzymes=enzyme_profile))
                 for c in cuts])
        except ValueError:
            continue
        try:
            chosen, score = best_set(options, fixed=fixed, matrix=matrix)
        except ValueError:
            continue
        out.append((list(cuts), chosen, score))
        # Enough plans already at the ceiling: nothing later can score better, and each
        # extra split costs a full overhang enumeration. Scoring all of them regressed a
        # 19S design from 6 s to 190 s.
        if sum(1 for _, _, s in out if s >= ceiling - 1e-12) >= OPTIMISE_ATTEMPTS:
            break
    if not out:
        raise ValueError("no split produced a valid overhang set")
    # best fidelity first; ties keep the balanced order balanced_cuts produced
    out.sort(key=lambda x: -x[2])
    return out, ceiling


def _rc(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _overhang_decisions(protein, cuts, chosen, profile) -> list[OverhangDecision]:
    """Every candidate at every junction, and why it was or was not used.

    Three outcomes, and the distinction matters: *rejected* means no synonymous
    realization exists under the active profile -- the junction is impossible, not merely
    worse. *Considered* means it was viable and another was preferred.
    """
    out: list[OverhangDecision] = []
    for cut, pick in zip(cuts, chosen):
        feasible = set(safe_overhangs(protein, cut, enzymes=profile))
        for oh in achievable_overhangs(protein, cut):
            if oh not in feasible:
                out.append(OverhangDecision(
                    junction_cut=cut, sequence=oh, status="rejected",
                    reason=f"no locally feasible synonymous realization under {profile}",
                    local_realizations=0))
            elif oh == pick:
                out.append(OverhangDecision(
                    junction_cut=cut, sequence=oh, status="selected",
                    reason="highest predicted fidelity among feasible candidates",
                    local_realizations=realization_count(protein, cut, oh, profile)))
            else:
                out.append(OverhangDecision(
                    junction_cut=cut, sequence=oh, status="considered",
                    reason="feasible, but another candidate scored at least as well",
                    local_realizations=realization_count(protein, cut, oh, profile)))
    return out


def design_oneshot(
    target_rna: str,
    organism: str = "c_reinhardtii_nuclear",
    genetic_code: int | None = None,
    codon_table=None,
    n_fragments: int | None = None,
    destination: tuple[str, str] | None = None,
    enzyme: str = "BsaI",
    enzyme_profile: str | tuple[str, ...] = C.DEFAULT_ENZYME_PROFILE,
    extra_blacklist: tuple[str, ...] | str = (),
    matrix: str = "BsaI-HFv2",
    seed: int = 42,
    check_offtarget: bool = True,
    outdir: str | Path | None = None,
) -> dict:
    """Design a complete PPR binder for `target_rna`, ready to order.

    `architecture` is not a parameter: it follows from the target's length, since a PPR
    needs one repeat per base. Pass a 9, 14 or 19 base target to get 9S, 14S or 19S.

    Returns the protein, the coding sequence, the oligo table, the QC verdict, the
    predicted ligation fidelity, a cost estimate and -- when `outdir` is given -- the
    paths written.
    """
    d = describe(target_rna)
    protein, architecture = d["aa_sequence"], d["architecture"]

    if genetic_code is None:
        genetic_code = ORGANISMS.get(organism, (None, 1))[1]
    table = _resolve_codon_table(codon_table, organism, genetic_code)

    if n_fragments is None:
        n_fragments = max(1, -(-len(protein) // MAX_FRAGMENT_AA))

    # A profile says which sites this project always excludes. `extra_blacklist` is for
    # the site a particular experiment needs clear -- an EcoRI in a downstream vector, say
    # -- which is a property of that experiment, not of CLIPPR's policy, so it is added
    # here rather than becoming a new profile.
    if extra_blacklist:
        extra = (tuple(e.strip() for e in extra_blacklist.split(",") if e.strip())
                 if isinstance(extra_blacklist, str) else tuple(extra_blacklist))
        enzyme_profile = tuple(dict.fromkeys(C.enzymes_for(enzyme_profile) + extra))

    plans, _ceiling = _plan_fragments(protein, n_fragments, destination, matrix,
                                      enzyme_profile)

    # A locked overhang can create a blacklisted enzyme site that no synonymous change
    # can remove, because EnforceSequence has frozen those bases -- DNA Chisel then
    # reports NoSolutionError and returns a sequence that satisfies nothing. Measured on
    # the corpus, this hits 13 of 50 19S targets. The cure is a different split, so try
    # the next one rather than hand back a coding sequence that failed its constraints.
    attempt = fallback = None
    for cuts, overhangs, fidelity in plans[:OPTIMISE_ATTEMPTS]:
        locked = {3 * c - 4: o for c, o in zip(cuts, overhangs)}
        opt = optimize_cds(protein, locked_sites=locked, codon_table=table,
                           genetic_code=genetic_code, seed=seed,
                           enzymes=enzyme_profile)
        if opt["constraints_ok"]:
            attempt = (cuts, overhangs, fidelity, opt)
            break
        if fallback is None:
            fallback = (cuts, overhangs, fidelity, opt)
    if attempt is None:
        if fallback is None:
            raise ValueError(f"no viable design for {target_rna}")
        attempt = fallback
    cuts, overhangs, fidelity, opt = attempt
    cds = opt["cds"]

    decisions = _overhang_decisions(protein, cuts, overhangs, enzyme_profile)

    # split_cds re-derives the overhangs from the CDS and raises if the lock slipped
    fragments = split_cds(cds, cuts, overhangs, destination, enzyme)
    oligos = build_oligos(cds, cuts, overhangs, destination, enzyme, prefix=target_rna)
    qc = synthesis_qc(cds)
    cost = opool_quote(oligos)

    dest = destination or C.DESTINATION_OVERHANGS["level0"]
    reaction = reaction_overhangs(overhangs, "level0") if destination is None else [
        dest[0], *overhangs, _rc(dest[1])]

    paths: dict[str, str] = {}
    if outdir is not None:
        out = Path(outdir)
        stem = f"clippr_{target_rna}"
        paths = {
            "oligo_csv": str(write_oligo_csv(oligos, out / f"{stem}_oligos.csv")),
            "oligo_fasta": str(write_fasta(oligos, out / f"{stem}_oligos.fasta")),
            "gene_fasta": str(write_gene_fasta(cds, target_rna, out / f"{stem}_gene.fasta",
                                               architecture)),
            "genbank": str(write_genbank(cds, fragments, target_rna, out / f"{stem}.gb",
                                         genetic_code, enzyme, destination)),
        }

    # Does the target also occur where a PPR could bind it? A PPR binds RNA, so the tier
    # that matters is an annotated transcript in the sense orientation, not genomic DNA on
    # either strand. Measured over 200 random targets: 16% of 9-nt targets occur in a
    # transcript (48% occur somewhere in the DNA, which overstates it threefold); 14-nt
    # and 19-nt targets occur in neither.
    offtarget = None
    if check_offtarget:
        try:
            from .offtarget import scan as _scan
            offtarget = _scan(d["target_rna"])
        except (OSError, ValueError, ImportError) as exc:
            # Only genuine unavailability -- no network, no cached genome -- is tolerated.
            # A broad `except` here previously hid a signature mismatch behind a cheerful
            # "not checked", which is exactly how a silent no-op ships.
            offtarget = {"verdict": "not checked",
                         "error": f"{type(exc).__name__}: {exc}",
                         "n_transcript": None, "n_genomic": None, "genes": []}

    warnings = list(qc["warnings"])
    if not opt["constraints_ok"]:
        warnings.append("codon optimiser could not satisfy every constraint")
    if offtarget and offtarget.get("n_transcript"):
        genes = ", ".join(offtarget["genes"][:3])
        warnings.append(
            f"this target occurs inside {offtarget['n_transcript']} annotated host "
            f"transcript(s) ({genes}), in the sense orientation — RNA a PPR could bind "
            f"as well as your construct. A longer target is the reliable fix: 14-nt and "
            f"19-nt targets do not occur. Occurrence is necessary for an off-target "
            f"interaction, not sufficient; no binding affinity is predicted."
        )
    elif offtarget and offtarget.get("n_genomic"):
        warnings.append(
            f"this target occurs {offtarget['n_genomic']}x in host genomic DNA but in no "
            f"annotated transcript, so there is no RNA for a PPR to bind at those loci. "
            f"Worth noting, not acting on."
        )
    backbone = set_fidelity([dest[0], _rc(dest[1])], matrix)
    if backbone < 0.95:
        warnings.append(
            f"the destination pair {dest[0]}/{_rc(dest[1])} scores {backbone:.3f} on its "
            f"own, and the whole reaction scores {fidelity:.3f}. The backbone is the "
            f"limit here, not the junctions we chose."
        )

    summary = (
        f"{target_rna} ({architecture}) -> {len(protein)} aa, {len(cds)} nt, "
        f"{len(oligos)} fragments. Fidelity {fidelity:.3f}. QC {qc['status']}. "
        f"{cost['total_eur']:.2f} {cost['currency']} (list price)."
    )

    audit = DesignAudit(
        target_rna=d["target_rna"],
        architecture=architecture,
        enzyme_profile=enzyme_profile if isinstance(enzyme_profile, str) else "custom",
        enzymes_excluded=tuple(C.enzymes_for(enzyme_profile)),
        genetic_code=genetic_code,
        organism=organism,
        cuts=tuple(cuts),
        overhang_decisions=tuple(decisions),
        fidelity=fidelity,
        fidelity_components=fidelity_components(reaction, matrix),
        fidelity_ceiling=backbone,
        constraints_satisfied=bool(opt["constraints_ok"]),
        qc_status=qc["status"],
        findings=tuple(warnings),
    )

    return {
        "summary": summary,
        "audit": audit,
        "target_rna": d["target_rna"],
        "architecture": architecture,
        "ppr_code": d["ppr_code"],
        "protein": protein,
        "cds": cds,
        "cuts": cuts,
        "junction_overhangs": overhangs,
        "reaction_overhangs": reaction,
        "fragments": fragments,
        "oligos": oligos,
        "qc": qc,
        "offtarget": offtarget,
        "fidelity": fidelity,
        "cost": cost,
        "constraints_ok": opt["constraints_ok"],
        "warnings": warnings,
        "paths": paths,
    }
