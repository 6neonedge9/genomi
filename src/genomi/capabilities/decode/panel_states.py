from __future__ import annotations

EMPTY_NATIVE_STATUSES = frozenset(
    {
        "requires_library_install",
        "source_unavailable",
        "out_of_scope_for_input",
        "skipped_missing_library",
        "skipped_tool_unavailable",
        "insufficient_overlap",
        "domain_id_required",
        "unknown_domain",
        "domain_out_of_scope_by_construction",
        "invalid_evidence_tier",
    }
)
EMPTY_PGX_STATUSES = EMPTY_NATIVE_STATUSES | frozenset(
    {
        "position_aware_pharmcat_export_required",
        "no_pharmcat_vcf_records",
        "active_genome_index_input_unavailable",
        "explicit_pharmcat_executable_unavailable",
        "no_pharmcat_artifacts",
    }
)
EMPTY_PRS_STATUSES = frozenset(
    {
        "requires_score_import",
        "requires_library_install",
        "out_of_scope_for_input",
        "source_unavailable",
    }
)
EMPTY_COVERAGE_STATES = frozenset({"in_scope_empty", "out_of_scope_for_input"})
