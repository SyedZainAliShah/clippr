"""Target RNA -> PPR code -> binder protein sequence.

The PPR code is modular: each repeat contacts one RNA base, and the identity of that base
is set by two residues (the fifth and the last) of the repeat. Designing a binder for a
chosen sequence is therefore a lookup, not a search.
"""
from __future__ import annotations

from .constants import (
    ARCHITECTURES,
    ARELF_MOTIF,
    BASE_TO_CODE,
    CODE_TO_BASE,
    N_TERMINAL,
    REPEAT_TEMPLATE,
)

_VALID = set("ACGU")


def normalize_target(rna: str) -> str:
    """Uppercase, T->U, strip whitespace; validate alphabet and length.

    Raises ValueError for anything that is not a 9, 14 or 19 nt RNA over ACGU.
    """
    s = "".join(str(rna).split()).upper().replace("T", "U")
    if not s:
        raise ValueError("empty target sequence")
    bad = sorted(set(s) - _VALID)
    if bad:
        raise ValueError(f"target contains non-ACGU characters: {bad}")
    if len(s) not in set(ARCHITECTURES.values()):
        raise ValueError(
            f"target length must be one of {sorted(ARCHITECTURES.values())}, got {len(s)}"
        )
    return s


def architecture_of(rna: str) -> str:
    """'9S', '14S' or '19S' for a normalised target."""
    n = len(rna)
    for name, length in ARCHITECTURES.items():
        if length == n:
            return name
    raise ValueError(f"no architecture for length {n}")


def rna_to_code(rna: str) -> list[str]:
    """One code pair per base, 5'->3'. e.g. 'UUA' -> ['ND', 'ND', 'TN']."""
    seq = normalize_target(rna)
    return [BASE_TO_CODE[b.replace("U", "T")] for b in seq]


def code_to_rna(codes: list[str]) -> str:
    """Inverse of rna_to_code."""
    return "".join(CODE_TO_BASE[c] for c in codes).replace("T", "U")


def code_to_protein(codes: list[str]) -> str:
    """Assemble the binder: solvating helix followed by one repeat per code pair.

    A code pair is written (fifth, last), matching the deposited module naming
    convention L<last>5<fifth>.
    """
    parts = [N_TERMINAL]
    for c in codes:
        if c not in CODE_TO_BASE:
            raise ValueError(f"unknown PPR code pair {c!r}")
        parts.append(REPEAT_TEMPLATE.format(fifth=c[0], last=c[1]))
    return "".join(parts)


def describe(rna: str, include_start_codon: bool = True) -> dict:
    """Full description of the binder for a target RNA."""
    seq = normalize_target(rna)
    codes = rna_to_code(seq)
    protein = code_to_protein(codes)
    if not include_start_codon and protein.startswith("M"):
        protein = protein[1:]
    return {
        "target_rna": seq,
        "architecture": architecture_of(seq),
        "ppr_code": "".join(codes),
        "code_pairs": codes,
        "aa_sequence": protein,
        "aa_length": len(protein),
        "n_arelf": protein.count(ARELF_MOTIF),
    }
