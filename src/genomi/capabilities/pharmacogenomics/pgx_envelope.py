from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any

from ...evidence import envelope as _env


JsonObject = dict[str, Any]
ExtraStatusHandler = Callable[[JsonObject, JsonObject, JsonObject, JsonObject], JsonObject | None]
NextActionFactory = Mapping[str, Any] | Callable[[JsonObject], Mapping[str, Any]]


@dataclass(frozen=True)
class PublicPGxSourceEnvelopeSpec:
    operation: str
    library_id: str
    source_id: str
    observation_keys: tuple[str, ...]
    invalid_target_reason: str
    invalid_target_action: str
    missing_inputs: tuple[str, ...]
    invalid_target_guidance: str
    source_unavailable_reason: str
    alternate_operations: tuple[str, ...]
    no_match_status: str
    no_match_action: str
    no_match_guidance: tuple[str, ...]
    evidence_guidance: tuple[str, ...]
    evidence_next_action: NextActionFactory
    no_match_target_fields: tuple[str, ...] = ()


def build_medication_review_envelope(
    *,
    query: JsonObject,
    evidence_state: JsonObject,
    evidence_matrix_traceability: JsonObject,
    sample_context_requested: bool,
    clinical_context_requested: bool,
    unanswered_answer_components: Any,
    source_availability: JsonObject,
) -> JsonObject:
    sources = source_availability.get("sources") or []
    consulted = [
        str(item.get("source_id"))
        for item in sources
        if item.get("source_id") and item.get("availability") not in {"unavailable", "source_unavailable"}
    ]
    unavailable = [
        str(item.get("source_id"))
        for item in sources
        if item.get("availability") in {"unavailable", "source_unavailable"}
    ]
    libraries = [_medication_review_library_use(item) for item in sources if item.get("source_id")]
    coverage = _env._coverage(libraries=libraries, consulted_sources=consulted, unavailable_sources=unavailable)
    observations = {
        "source_evidence_count": evidence_state.get("source_evidence_count"),
        "sample_evidence_count": evidence_state.get("sample_evidence_count"),
        "evidence_matrix_item_count": evidence_matrix_traceability.get("item_count"),
        "unresolved_components": unanswered_answer_components,
    }
    personal_context = _env._personal_context(uses_personal_dna=bool(sample_context_requested))
    scope_payload = {
        "drug": query.get("drug"),
        "gene": query.get("gene"),
        "rsid": query.get("rsid"),
        "atc_code": query.get("atc_code"),
        "drugbank_id": query.get("drugbank_id"),
        "genome_build": query.get("genome_build"),
        "sample_context_requested": sample_context_requested,
        "clinical_context_requested": clinical_context_requested,
    }
    has_public = bool(evidence_state.get("has_public_pgx_evidence"))
    has_sample = bool(evidence_state.get("has_sample_evidence"))
    if source_availability.get("status") == "source_unavailable_no_evidence":
        return _env.not_assessed(
            operation="pharmacogenomics.review_medication",
            reason="All consulted PGx sources were unavailable.",
            query_scope=scope_payload,
            personal_context=personal_context,
            coverage=coverage,
            observations=observations,
        )
    if not has_public and not has_sample:
        return _env.empty_consulted_scope(
            operation="pharmacogenomics.review_medication",
            query_scope=scope_payload,
            personal_context=personal_context,
            coverage=coverage,
            observations=observations,
        )
    answer_readiness = _env.NEEDS_CLINICAL_CONFIRMATION if (clinical_context_requested and has_public) else _env.SCOPED_ANSWER_ONLY
    return _env.evidence_present(
        operation="pharmacogenomics.review_medication",
        query_scope=scope_payload,
        personal_context=personal_context,
        coverage=coverage,
        observations=observations,
        answer_readiness=answer_readiness,
    )


def build_public_pgx_source_envelope(
    result: JsonObject,
    spec: PublicPGxSourceEnvelopeSpec,
    *,
    extra_status_handler: ExtraStatusHandler | None = None,
) -> JsonObject:
    target = dict(result.get("query") or {})
    observations = _observations(result, spec)
    coverage = _coverage(result, spec)
    status = observations["status"]

    if status == "invalid_target":
        return _env.not_assessed(
            operation=spec.operation,
            reason=spec.invalid_target_reason,
            query_scope=target,
            coverage=coverage,
            observations=observations,
            next_actions=[
                {
                    "action": spec.invalid_target_action,
                    "missing_inputs": list(spec.missing_inputs),
                }
            ],
            guidance=[spec.invalid_target_guidance],
        )
    if status == "source_unavailable":
        return _env.not_assessed(
            operation=spec.operation,
            reason=spec.source_unavailable_reason,
            query_scope=target,
            coverage=coverage,
            observations=observations,
            next_actions=[
                {
                    "action": "use_alternate_pgx_source_or_retry",
                    "operations": list(spec.alternate_operations),
                }
            ],
            guidance=["source_unavailable:retry_or_use_other_pgx_sources"],
        )

    if extra_status_handler is not None:
        handled = extra_status_handler(result, target, coverage, observations)
        if handled is not None:
            return handled

    if status == spec.no_match_status:
        next_action: JsonObject = {
            "action": spec.no_match_action,
            "operations": list(spec.alternate_operations),
        }
        if spec.no_match_target_fields:
            next_action["target_fields"] = list(spec.no_match_target_fields)
        return _env.empty_consulted_scope(
            operation=spec.operation,
            query_scope=target,
            coverage=coverage,
            observations=observations,
            next_actions=[next_action],
            guidance=list(spec.no_match_guidance),
        )

    return _env.evidence_present(
        operation=spec.operation,
        query_scope=target,
        coverage=coverage,
        observations=observations,
        answer_readiness=_env.SCOPED_ANSWER_ONLY,
        next_actions=[_evidence_next_action(spec.evidence_next_action, result)],
        guidance=list(spec.evidence_guidance),
    )


def _observations(result: JsonObject, spec: PublicPGxSourceEnvelopeSpec) -> JsonObject:
    summary = dict(result.get("summary") or {})
    observations: JsonObject = {"status": str(result.get("status") or "")}
    for key in spec.observation_keys:
        observations[key] = summary.get(key, 0)
    return observations


def _coverage(result: JsonObject, spec: PublicPGxSourceEnvelopeSpec) -> JsonObject:
    status = str(result.get("status") or "")
    raw_calls = result.get("raw_calls") or []
    return _env._coverage(
        libraries=[
            {
                "library": spec.library_id,
                "state": "failed" if status == "source_unavailable" else "installed",
            }
        ],
        consulted_sources=[spec.source_id] if raw_calls and status != "source_unavailable" else [],
        unavailable_sources=[spec.source_id] if status == "source_unavailable" else [],
    )


def _medication_review_library_use(source: JsonObject) -> JsonObject:
    source_id = str(source.get("source_id") or "")
    library_id = {"fda_pgx": "fda-pgx"}.get(source_id, source_id)
    failed_states = {"unavailable", "source_unavailable"}
    return {
        "library": library_id,
        "state": "failed" if source.get("availability") in failed_states else "installed",
    }


def _evidence_next_action(action: NextActionFactory, result: JsonObject) -> JsonObject:
    if callable(action):
        return dict(action(result))
    return dict(action)
