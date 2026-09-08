"""CLIPPR — design tooling for synthetic PPR regulators.

iGEM Marburg 2026. Chloroplast Library for inducible PPR-mediated Post-transcriptional
Gene Expression Regulation.

Nothing in this package imports any third-party design pipeline. The corpus
implementation is used only under `validation/`, as a test oracle.
"""
from .arelf import (
    achievable_overhangs,
    balanced_cuts,
    candidate_cuts,
    realization_count,
    safe_overhangs,
    cut_context,
    find_motifs,
    offsets_from_motif,
    tunable_cuts,
)
from .assembly import (
    Fragment,
    build_oligos,
    enzyme_geometry,
    reassemble,
    split_cds,
    wrap_fragment,
)
from .codons import (
    BLACKLIST_ENZYMES,
    CHLAMYDOMONAS_TAXID,
    complete_table,
    optimize_cds,
    table_from_cds_fasta,
    table_from_csv,
    table_from_kazusa,
)
from .audit import DesignAudit, OverhangDecision
from .design import ORGANISMS, design_oneshot
from .policy import DEFAULT_ENZYME_PROFILE, ENZYME_PROFILES, ENZYME_ROLES, enzymes_for
from .qc import DEFAULT_THRESHOLDS, synthesis_qc
from .export import (
    opool_quote,
    write_fasta,
    write_gene_fasta,
    write_genbank,
    write_oligo_csv,
)
from .constants import (
    ARCHITECTURES,
    ARELF_MOTIF,
    BASE_TO_CODE,
    CODE_TO_BASE,
    DESTINATION_OVERHANGS,
    GENETIC_CODE_NUCLEAR,
    GENETIC_CODE_PLASTID,
    N_TERMINAL,
    PROVENANCE,
    REPEAT_TEMPLATE,
    SYNONYMOUS_CODONS,
)
from .overhangs import (
    best_set,
    enumerate_candidates,
    fidelity_components,
    load_matrix,
    reaction_overhangs,
    reverse_complement,
    set_fidelity,
    valid_set,
)
from .orthogonal import (
    POSITION_WEIGHTS,
    TargetSet,
    crosstalk_report,
    design_orthogonal_set,
    distance,
)
from .ppr import (
    architecture_of,
    code_to_protein,
    code_to_rna,
    describe,
    normalize_target,
    rna_to_code,
)

__version__ = "0.1.0"

__all__ = [
    "ARCHITECTURES", "ARELF_MOTIF", "BASE_TO_CODE", "CODE_TO_BASE",
    "DESTINATION_OVERHANGS", "GENETIC_CODE_NUCLEAR", "GENETIC_CODE_PLASTID",
    "N_TERMINAL", "PROVENANCE", "REPEAT_TEMPLATE", "SYNONYMOUS_CODONS",
    "achievable_overhangs", "candidate_cuts", "cut_context", "find_motifs",
    "offsets_from_motif", "tunable_cuts",
    "BLACKLIST_ENZYMES", "CHLAMYDOMONAS_TAXID", "complete_table", "optimize_cds",
    "table_from_cds_fasta", "table_from_csv", "table_from_kazusa",
    "Fragment", "build_oligos", "enzyme_geometry", "reassemble", "split_cds",
    "wrap_fragment",
    "opool_quote", "write_fasta", "write_gene_fasta", "write_genbank", "write_oligo_csv",
    "DEFAULT_THRESHOLDS", "synthesis_qc",
    "ORGANISMS", "design_oneshot", "balanced_cuts",
    "DesignAudit", "OverhangDecision", "realization_count", "safe_overhangs",
    "DEFAULT_ENZYME_PROFILE", "ENZYME_PROFILES", "ENZYME_ROLES", "enzymes_for",
    "best_set", "enumerate_candidates", "fidelity_components", "load_matrix",
    "reaction_overhangs", "reverse_complement", "set_fidelity", "valid_set",
    "POSITION_WEIGHTS", "TargetSet", "crosstalk_report", "design_orthogonal_set",
    "distance",
    "architecture_of", "code_to_protein", "code_to_rna", "describe",
    "normalize_target", "rna_to_code",
    "__version__",
]
