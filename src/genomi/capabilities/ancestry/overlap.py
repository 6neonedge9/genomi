from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from ...active_genome_index.active_genome_index import ActiveGenomeIndexReader
from ...evidence import envelope as evidence_envelope
from . import policy, reference_panels, source_context

JsonObject = dict[str, Any]
# Overlap is graded purely as a fraction of the loaded panel. The
# 1000G-30x-GRCh38 panel is deliberately compact (~20k LD-pruned,
# MAF-filtered, ancestry-informative SNPs per genomi-ancestry-panel/
# docs/filters.md — "ancestry PCA on 3,202 samples stabilizes well below
# 10,000 informative markers"). Small-by-design ≠ low quality, so there
# is no absolute marker-count floor; what matters is how much of the
# chosen panel the sample's AGI actually covers.
HIGH_OVERLAP_FRACTION = policy.HIGH_OVERLAP_FRACTION
MODERATE_OVERLAP_FRACTION = policy.MODERATE_OVERLAP_FRACTION
LOW_OVERLAP_FRACTION = policy.LOW_OVERLAP_FRACTION


def check_sample_overlap(
    agi_reader: ActiveGenomeIndexReader,
    *,
    genome_build: str = "GRCh38",
    panel_root: str | Path | None = None,
) -> JsonObject:
    panel_or_missing = _load_panel_or_missing(genome_build, panel_root)
    if isinstance(panel_or_missing, dict) and panel_or_missing.get("status") == "panel_not_installed":
        return panel_or_missing
    panel = panel_or_missing
    genotype_context = collect_sample_genotypes(
        agi_reader,
        genome_build=genome_build,
        panel=panel,
    )
    sample_qc = genotype_context["sample_qc"]
    result = {
        "status": sample_qc["overlap_status"],
        "personal_context": {"uses_personal_dna": True},
        "reference_panel": _reference_panel_summary(panel),
        "sample_qc": sample_qc,
        "limitations": source_context.limitations(),
        "next_actions": _overlap_next_actions(sample_qc),
    }
    result["evidence_envelope"] = _overlap_envelope("ancestry.check_sample_overlap", result)
    return result


def _load_panel_or_missing(
    genome_build: str,
    panel_root: str | Path | None,
) -> JsonObject:
    """Load the panel that matches the sample's genome build.

    Returns the panel payload on success, or a ``panel_not_installed``
    envelope payload that the caller should propagate directly when the
    matching panel library is missing on disk.
    """

    normalized_build = _normalize_build(genome_build)
    try:
        return reference_panels.load_panel(normalized_build, panel_root)
    except FileNotFoundError:
        return _panel_not_installed_payload(
            genome_build=normalized_build,
        )


def _panel_not_installed_payload(*, genome_build: str) -> JsonObject:
    from ...runtime.libraries import manager

    library = source_context.panel_library_for_build(genome_build)
    status = manager.status(library)
    panel_id = source_context.panel_id_for_build(genome_build)
    note = (
        f"No ancestry panel is installed for the sample's {genome_build} build. "
        f"Install {library} to enable reference-panel projection for this sample."
    )
    sample_qc = {
        "genome_build": genome_build,
        "panel_marker_count": 0,
        "usable_marker_count": 0,
        "missing_marker_count": 0,
        "overlap_fraction": 0.0,
        "overlap_status": "panel_not_installed",
        "projection_allowed": False,
        "marker_overlap_quality": "unavailable",
        "required_library": library,
        "required_panel_id": panel_id,
        "install_command": status["install_command"],
        "note": note,
    }
    result = {
        "status": "panel_not_installed",
        "personal_context": {"uses_personal_dna": True},
        "reference_panel": {
            "panel_id": panel_id,
            "title": (
                source_context.PANEL_TITLE_GRCH38
                if genome_build == "GRCh38"
                else source_context.PANEL_TITLE_GRCH37
            ),
            "library": library,
            "genome_build": genome_build,
            "installed": False,
            "source_urls": source_context.source_urls(),
        },
        "sample_qc": sample_qc,
        "limitations": source_context.limitations(),
        "next_actions": [
            {
                "action": "install_library",
                "library": library,
                "install_command": status["install_command"],
                "reason": note,
            }
        ],
    }
    result["evidence_envelope"] = _overlap_envelope("ancestry.check_sample_overlap", result)
    return result


def collect_sample_genotypes(
    agi_reader: ActiveGenomeIndexReader,
    *,
    genome_build: str = "GRCh38",
    panel: JsonObject | None = None,
) -> JsonObject:
    normalized_build = _normalize_build(genome_build)
    panel_payload = panel or reference_panels.load_panel(normalized_build)
    markers = list(panel_payload["markers"])
    panel_marker_count = len(markers)

    # No readiness / incompleteness handling here: open_agi gated access
    # upstream (missing / incomplete -> active_genome_index_incomplete). A
    # variants_ready index proceeds; the dispatch chokepoint stamps
    # reference_pending.
    marker_results = _marker_dosage_results(agi_reader, markers)
    dosages: dict[str, float] = {}
    missing_marker_ids: list[str] = []
    missing_marker_reasons: dict[str, int] = {}
    missing_marker_examples: list[JsonObject] = []
    for marker, marker_result in zip(markers, marker_results):
        marker_id = str(marker["marker_id"])
        if marker_result.get("status") != "matched":
            missing_marker_ids.append(marker_id)
            reason = _ancestry_missing_reason(str(marker_result.get("reason") or "unusable_marker"))
            missing_marker_reasons[reason] = missing_marker_reasons.get(reason, 0) + 1
            if len(missing_marker_examples) < 10:
                detail = _ancestry_missing_detail(marker_result)
                missing_marker_examples.append(
                    {
                        "marker_id": marker_id,
                        "reason": reason,
                        **({"detail": detail} if detail else {}),
                    }
                )
            continue
        dosage = marker_result.get("dosage")
        if dosage is None:
            dosage = marker_result.get("effect_allele_dosage")
        if dosage is None or not math.isfinite(float(dosage)):
            missing_marker_ids.append(marker_id)
            missing_marker_reasons["nonfinite_dosage"] = missing_marker_reasons.get("nonfinite_dosage", 0) + 1
            continue
        dosages[marker_id] = float(dosage)

    usable_marker_ids = [str(marker["marker_id"]) for marker in markers if str(marker["marker_id"]) in dosages]
    usable_marker_count = len(usable_marker_ids)
    fraction = _overlap_fraction(usable_marker_count, panel_marker_count)
    sample_qc = _sample_qc(
        marker_count=panel_marker_count,
        usable_marker_count=usable_marker_count,
        missing_marker_count=len(missing_marker_ids),
        genome_build=normalized_build,
        overlap_status=_overlap_status(fraction),
        projection_allowed=fraction >= LOW_OVERLAP_FRACTION,
        marker_overlap_quality=_marker_overlap_quality(fraction),
        note=_overlap_note(fraction),
        missing_marker_reasons=missing_marker_reasons,
        missing_marker_examples=missing_marker_examples,
    )
    return {
        "sample_qc": sample_qc,
        "dosages": dosages,
        "usable_marker_ids": usable_marker_ids,
        "missing_marker_ids": missing_marker_ids,
        "missing_marker_reasons": missing_marker_reasons,
        "missing_marker_examples": missing_marker_examples,
    }


def _marker_dosage_results(agi_reader: ActiveGenomeIndexReader, markers: list[JsonObject]) -> list[JsonObject]:
    dosage_variants = [
        {
            "variant_index": index,
            "variant_id": marker.get("marker_id"),
            "rsid": marker.get("marker_id"),
            "chrom": marker["chrom"],
            "pos": int(marker["pos"]),
            "reference_allele": str(marker["ref"]).upper(),
            "effect_allele": str(marker["alt"]).upper(),
            "other_allele": str(marker["ref"]).upper(),
            "effect_weight": 1.0,
        }
        for index, marker in enumerate(markers)
    ]
    raw_results = agi_reader.dosage_for_variants(
        dosage_variants,
        skip_ambiguous_palindromic=False,
    )
    results_by_index = {
        int(result["variant_index"]): result
        for result in raw_results
        if result.get("variant_index") is not None
    }
    return [
        results_by_index.get(
            index,
            {
                "status": "missing",
                "reason": "no_record_at_locus",
                "variant_index": index,
                "variant_id": marker.get("marker_id"),
            },
        )
        for index, marker in enumerate(markers)
    ]


def _ancestry_missing_reason(reason: str) -> str:
    return {
        "no_call": "missing_genotype",
        "filter_fail": "filtered_record",
        "genotype_allele_outside_score_alleles": "genotype_allele_outside_panel_alleles",
        "score_allele_model_not_supported_by_array_genotype": "panel_allele_model_not_supported_by_array_genotype",
        "effect_allele_not_supported_by_array_genotype": "marker_allele_not_supported_by_array_genotype",
        "effect_allele_not_in_record": "alternate_allele_mismatch",
        "other_allele_not_in_record": "reference_allele_mismatch",
    }.get(reason, reason)


def _ancestry_missing_detail(marker_result: JsonObject) -> list[JsonObject]:
    record = marker_result.get("record")
    if not isinstance(record, dict):
        return []
    detail: JsonObject = {
        "reason": _ancestry_missing_reason(str(marker_result.get("reason") or "unusable_marker")),
    }
    basis = _record_basis(record)
    if basis:
        detail["basis"] = basis
    if record.get("observed_alleles"):
        detail["allele_bases"] = record["observed_alleles"]
    return [detail]


def _record_basis(record: JsonObject) -> str | None:
    record_kind = str(record.get("record_kind") or "")
    if record_kind.startswith("array_"):
        return "consumer_array"
    if record_kind == "reference_block":
        return "reference_block"
    if record_kind == "variant_call":
        return "exact_genotype"
    return None


def _sample_qc(
    *,
    marker_count: int,
    usable_marker_count: int,
    missing_marker_count: int,
    genome_build: str,
    overlap_status: str,
    projection_allowed: bool,
    marker_overlap_quality: str,
    note: str,
    missing_marker_reasons: dict[str, int] | None = None,
    missing_marker_examples: list[JsonObject] | None = None,
) -> JsonObject:
    return {
        "genome_build": genome_build,
        "supported_genome_builds": list(reference_panels.SUPPORTED_BUILDS),
        "panel_marker_count": marker_count,
        "usable_marker_count": usable_marker_count,
        "missing_marker_count": missing_marker_count,
        "missing_marker_reasons": dict(missing_marker_reasons or {}),
        "missing_marker_examples": list(missing_marker_examples or []),
        "overlap_fraction": usable_marker_count / marker_count if marker_count else 0.0,
        "overlap_status": overlap_status,
        "projection_allowed": projection_allowed,
        "marker_overlap_quality": marker_overlap_quality,
        "thresholds": policy.overlap_thresholds(),
        "note": note,
    }


def _overlap_fraction(usable_marker_count: int, panel_marker_count: int) -> float:
    return usable_marker_count / panel_marker_count if panel_marker_count else 0.0


def _overlap_status(fraction: float) -> str:
    return policy.overlap_status(fraction)


def _marker_overlap_quality(fraction: float) -> str:
    return policy.marker_overlap_quality(fraction)


def _overlap_note(fraction: float) -> str:
    return policy.overlap_note(fraction)


def _overlap_next_actions(sample_qc: JsonObject) -> list[JsonObject]:
    status = str(sample_qc.get("overlap_status") or "")
    if status == "panel_not_installed":
        return [
            {
                "action": "install_library",
                "library": sample_qc.get("required_library"),
                "install_command": sample_qc.get("install_command"),
                "reason": sample_qc.get("note"),
            }
        ]
    if status == "active_genome_index_incomplete":
        return [{"action": "parse_source", "operation": "genomi.parse_source"}]
    if status == "insufficient_overlap":
        return [{"action": "use_higher_overlap_index", "reason": sample_qc.get("note")}]
    return [{"action": "project_pca", "operation": "ancestry.project_pca"}]


def _reference_panel_summary(panel: JsonObject) -> JsonObject:
    manifest = panel.get("manifest") or {}
    stats = panel.get("stats") or {}
    genome_build = str(manifest.get("genome_build") or "GRCh38")
    try:
        title = source_context.panel_title_for_build(genome_build)
    except ValueError:
        title = str(manifest.get("title") or reference_panels.PANEL_TITLE)
    return {
        "panel_id": str(manifest.get("panel_id") or source_context.panel_id_for_build(genome_build)),
        "title": title,
        "library": str(manifest.get("library") or source_context.panel_library_for_build(genome_build)),
        "genome_build": genome_build,
        "sample_count": int(manifest.get("sample_count") or stats.get("sample_count") or len(panel.get("samples") or [])),
        "marker_count": int(manifest.get("marker_count") or stats.get("marker_count") or len(panel.get("markers") or [])),
        "component_count": int(manifest.get("component_count") or stats.get("component_count") or len(panel.get("component_names") or [])),
        "label_scope": "1000 Genomes reference-panel population labels",
        "source_urls": source_context.source_urls(),
    }


def _overlap_envelope(operation: str, result: JsonObject) -> JsonObject:
    status = str(result.get("status") or "")
    sample_qc = result.get("sample_qc") if isinstance(result.get("sample_qc"), dict) else {}
    reference_panel = result.get("reference_panel") if isinstance(result.get("reference_panel"), dict) else {}
    panel_id = str(reference_panel.get("panel_id") or reference_panels.PANEL_ID)
    panel_library = str(reference_panel.get("library") or reference_panels.PANEL_LIBRARY)
    panel_title = str(reference_panel.get("title") or reference_panels.PANEL_TITLE)
    library_state = "installed" if reference_panel.get("installed", True) else "missing"
    coverage = evidence_envelope._coverage(
        libraries=[{"library": panel_library, "state": library_state, "title": panel_title}],
        consulted_sources=["active_genome_index", panel_library],
    )
    observations = {
        "status": status,
        "panel_marker_count": sample_qc.get("panel_marker_count"),
        "usable_marker_count": sample_qc.get("usable_marker_count"),
        "missing_marker_count": sample_qc.get("missing_marker_count"),
        "overlap_fraction": sample_qc.get("overlap_fraction"),
        "marker_overlap_quality": sample_qc.get("marker_overlap_quality"),
        "projection_allowed": sample_qc.get("projection_allowed"),
    }
    if status == "completed":
        return evidence_envelope.evidence_present(
            operation=operation,
            query_scope={"method": "ancestry_marker_overlap", "reference_panel": panel_id},
            personal_context={"uses_personal_dna": True},
            coverage=coverage,
            observations=observations,
            answer_readiness=evidence_envelope.SCOPED_ANSWER_ONLY,
            next_actions=result.get("next_actions") if isinstance(result.get("next_actions"), list) else [],
            notes=[source_context.BOUNDARY_NOTE],
            guidance=["evidence_present:answer_as_marker_overlap_only"],
        )
    return evidence_envelope.not_assessed(
        operation=operation,
        reason=sample_qc.get("note") or status,
        query_scope={"method": "ancestry_marker_overlap", "reference_panel": panel_id},
        personal_context={"uses_personal_dna": True},
        coverage=coverage,
        observations=observations,
        next_actions=result.get("next_actions") if isinstance(result.get("next_actions"), list) else [],
        notes=[source_context.BOUNDARY_NOTE],
        guidance=["not_assessed:do_not_interpret_reference_similarity"],
    )


def _normalize_build(value: str | None) -> str:
    return policy.normalize_build(value, default="unknown")
