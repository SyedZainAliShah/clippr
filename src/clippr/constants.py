"""Compatibility aggregator over the three constant layers.

**Prefer importing from the specific layer.** This module flattens all three into one
namespace, which is convenient but loses the distinction that matters most when reading
the package:

    from .biology import REPEAT_TEMPLATE        # published fact -- GRASP requires this
    from .assembly_spec import DESTINATION_OVERHANGS   # how the assembly works
    from .policy import ENZYME_PROFILES         # what CLIPPR chose

A contributor who sees `constants.ENZYME_PROFILES` beside `constants.REPEAT_TEMPLATE`
cannot tell that one is a published scaffold sequence and the other is this project's
compatibility preference. That confusion is how a project policy gets defended as though
it were biology -- and how excluding BsmBI, which no source requires, went unquestioned
long enough to make 13 designs infeasible.

Kept because existing modules and tests import `constants` as a single namespace.
"""
from __future__ import annotations

from .assembly_spec import (
    ARELF_OFFSETS,
    ARELF_OFFSETS_OBSERVED,
    CUT_IS_CODON_ALIGNED,
    DESTINATION_OVERHANGS,
)
from .biology import (
    ARCHITECTURES,
    ARELF_MOTIF,
    BASE_TO_CODE,
    CODE_TO_BASE,
    GENETIC_CODE_NUCLEAR,
    GENETIC_CODE_PLASTID,
    N_TERMINAL,
    PROVENANCE,
    REPEAT_TEMPLATE,
    SYNONYMOUS_CODONS,
)
from .policy import (
    CHLAMYDOMONAS_TAXID,
    DEFAULT_ENZYME_PROFILE,
    ENZYME_PROFILES,
    ENZYME_ROLES,
    MAX_FRAGMENT_AA,
    MIN_FRAGMENT_AA,
    ORGANISMS,
    REPEAT_K,
    enzymes_for,
)

#: Which layer each name comes from. Lets a reader -- or a test -- check that a value is
#: being treated as the kind of claim it actually is.
LAYERS: dict[str, tuple[str, ...]] = {
    "biology": ("CODE_TO_BASE", "BASE_TO_CODE", "N_TERMINAL", "REPEAT_TEMPLATE",
                "ARELF_MOTIF", "ARCHITECTURES", "SYNONYMOUS_CODONS", "PROVENANCE",
                "GENETIC_CODE_NUCLEAR", "GENETIC_CODE_PLASTID"),
    "assembly_spec": ("DESTINATION_OVERHANGS", "ARELF_OFFSETS", "ARELF_OFFSETS_OBSERVED",
                      "CUT_IS_CODON_ALIGNED"),
    "policy": ("ENZYME_ROLES", "ENZYME_PROFILES", "DEFAULT_ENZYME_PROFILE", "enzymes_for",
               "CHLAMYDOMONAS_TAXID", "ORGANISMS", "MAX_FRAGMENT_AA", "MIN_FRAGMENT_AA",
               "REPEAT_K"),
}

__all__ = [n for names in LAYERS.values() for n in names]
