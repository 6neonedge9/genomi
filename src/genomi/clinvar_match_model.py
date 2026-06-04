from __future__ import annotations


MATCH_BASIS_EXACT_ALLELE = "exact_allele"
MATCH_BASIS_MULTIALLELIC_ALT = "multiallelic_alt"
MATCH_BASIS_CONSUMER_ARRAY_ALLELE_INFERENCE = "consumer_array_allele_inference"
MATCH_BASIS_LIFTOVER_EXACT_ALLELE = "liftover_exact_allele"
MATCH_BASIS_LIFTOVER_MULTIALLELIC_ALT = "liftover_multiallelic_alt"
MATCH_BASIS_VALUES = {
    MATCH_BASIS_EXACT_ALLELE,
    MATCH_BASIS_MULTIALLELIC_ALT,
    MATCH_BASIS_CONSUMER_ARRAY_ALLELE_INFERENCE,
    MATCH_BASIS_LIFTOVER_EXACT_ALLELE,
    MATCH_BASIS_LIFTOVER_MULTIALLELIC_ALT,
}


def match_basis_for_sample_mode(sample_mode: str, *, cross_build: bool = False) -> str:
    """Return the public ClinVar provenance basis for an AGI-staged sample mode."""

    if sample_mode not in {
        MATCH_BASIS_EXACT_ALLELE,
        MATCH_BASIS_MULTIALLELIC_ALT,
        MATCH_BASIS_CONSUMER_ARRAY_ALLELE_INFERENCE,
    }:
        raise ValueError(f"unknown ClinVar sample match mode: {sample_mode}")
    if not cross_build:
        return sample_mode
    if sample_mode == MATCH_BASIS_EXACT_ALLELE:
        return MATCH_BASIS_LIFTOVER_EXACT_ALLELE
    if sample_mode == MATCH_BASIS_MULTIALLELIC_ALT:
        return MATCH_BASIS_LIFTOVER_MULTIALLELIC_ALT
    return sample_mode


def match_basis_sql_expression(sample_mode_expression: str, *, cross_build: bool = False) -> str:
    """SQL expression equivalent of :func:`match_basis_for_sample_mode`."""

    if not cross_build:
        return sample_mode_expression
    return f"""
        case {sample_mode_expression}
            when '{MATCH_BASIS_EXACT_ALLELE}' then '{MATCH_BASIS_LIFTOVER_EXACT_ALLELE}'
            when '{MATCH_BASIS_MULTIALLELIC_ALT}' then '{MATCH_BASIS_LIFTOVER_MULTIALLELIC_ALT}'
            else {sample_mode_expression}
        end
    """


def evidence_scope_for_match_basis(match_basis: str) -> str:
    if match_basis == MATCH_BASIS_CONSUMER_ARRAY_ALLELE_INFERENCE:
        return "consumer_array_inferred_allele"
    if match_basis.startswith("liftover_"):
        return "liftover_sample_allele"
    if match_basis in {MATCH_BASIS_MULTIALLELIC_ALT, MATCH_BASIS_LIFTOVER_MULTIALLELIC_ALT}:
        return "selected_alternate_allele"
    return "sample_allele"
