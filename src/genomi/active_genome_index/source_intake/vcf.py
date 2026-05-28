from __future__ import annotations

import contextlib
from pathlib import Path
from typing import Any

from ...runtime.paths import (
    run_evidence_db_path_for_source,
    run_evidence_dir_for_source,
    run_project_dir_for_source,
    run_reference_dir_for_source,
    run_work_dir_for_source,
    sample_slug_from_source,
    shared_evidence_db_path,
)
from ...runtime.static_dependencies import resolve_genome_build
from ..active_genome_index import (
    active_genome_index_readiness,
    create_active_genome_index,
    default_active_genome_index_path,
)
from ..canonical import build_canonical_bgzip
from .agi_store import SOURCE_PARSE_SCHEMA, JsonObject, _init_source_evidence_db
from .detection import SourceDetection


def _parse_vcf_active_genome_index(
    source_path: Path,
    *,
    detection: SourceDetection,
    evidence_db: str | Path | None,
    source_evidence_db: str | Path | None,
    shared_evidence_db: str | Path | None,
    genome_build: str,
    force: bool,
    max_records: int | None,
    parallel_workers: int | None,
) -> JsonObject:
    effective_build = resolve_genome_build(source_path, genome_build)
    project_dir = run_project_dir_for_source(source_path, source_format=detection.source_format)
    work_dir = run_work_dir_for_source(source_path, source_format=detection.source_format)
    evidence_dir = run_evidence_dir_for_source(source_path, source_format=detection.source_format)
    reference_dir = run_reference_dir_for_source(source_path, source_format=detection.source_format)
    db_path = Path(evidence_db) if evidence_db is not None else run_evidence_db_path_for_source(source_path, source_format=detection.source_format)
    shared_db = Path(shared_evidence_db) if shared_evidence_db is not None else shared_evidence_db_path()
    active_genome_index_path = default_active_genome_index_path(source_path)

    project_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    evidence_dir.mkdir(parents=True, exist_ok=True)
    reference_dir.mkdir(parents=True, exist_ok=True)
    _init_source_evidence_db(
        db_path,
        source_path,
        source_format=detection.source_format,
        source_evidence_db=source_evidence_db,
        shared_evidence_db=shared_db,
    )
    # The structured Active Genome Index is self-sufficient: once built, no
    # capability reopens the intake or the canonical to understand a genome.
    # If a complete, current-schema index already exists we skip touching the
    # source entirely (it may even be gone). Otherwise we materialize a
    # canonical bgzip ONLY to parse it into the index, then delete it to
    # reclaim disk — the index is the sole source of truth afterward.
    #
    # gVCFs are ~96% reference blocks, so we parse in two phases: Phase A
    # stores every variant and returns variants_ready (the whole interpretation
    # surface is live in minutes); Phase B appends the reference-block tail as a
    # detached background job. Only gVCFs get phased — a plain VCF has no
    # reference tail, and max_records / forced builds stay single-phase.
    readiness = active_genome_index_readiness(active_genome_index_path)
    two_phase = detection.source_format == "gvcf" and max_records is None
    reference_job: JsonObject | None = None
    if readiness.get("complete") and not force:
        active_genome_index_result = {
            "status": "cached",
            "active_genome_index_complete": True,
            "active_genome_index_path": str(active_genome_index_path),
        }
    elif two_phase and readiness.get("variants_ready") and not force:
        # Phase A already done from a prior call; just make sure the reference
        # tail is (still) being built instead of rebuilding variants.
        active_genome_index_result = {
            "status": "variants_ready",
            "active_genome_index_complete": False,
            "reference_pending": True,
            "active_genome_index_path": str(active_genome_index_path),
        }
        reference_job = _enqueue_reference_pass(active_genome_index_path, parallel_workers)
    else:
        canonical_result = build_canonical_bgzip(source_path, work_dir, force=force)
        canonical_path = Path(canonical_result["canonical_path"])
        active_genome_index_result = create_active_genome_index(
            canonical_path,
            active_genome_index_path,
            include_reference=True,
            max_records=max_records,
            parallel_workers=parallel_workers,
            reuse_existing=not force,
            defer_reference=two_phase,
        )
        # Drop the work-dir canonical; the index owns a per-index canonical
        # (metadata.vcf_path) that Phase B and every capability tool read from.
        for stale in (canonical_path, Path(str(canonical_path) + ".gzi")):
            with contextlib.suppress(FileNotFoundError):
                stale.unlink()
        if active_genome_index_result.get("status") == "variants_ready":
            reference_job = _enqueue_reference_pass(active_genome_index_path, parallel_workers)
    outputs: dict[str, Any] = {"active_genome_index_path": str(active_genome_index_path)}
    if reference_job is not None:
        outputs["reference_pass_job_id"] = reference_job.get("job_id")
        outputs["reference_pass_job_path"] = reference_job.get("job_path")
    return {
        "schema": SOURCE_PARSE_SCHEMA,
        "workflow_area": "active-genome-index",
        "status": "completed",
        "source": str(source_path),
        "vcf": str(source_path),
        "source_format": detection.source_format,
        "source_kind": detection.source_kind,
        "provider": detection.provider,
        "annotation_scope": "active_genome_index",
        "sample_slug": sample_slug_from_source(source_path, source_format=detection.source_format),
        "genome_build": effective_build,
        "evidence_db": str(db_path),
        "shared_evidence_db": str(shared_db),
        "project_dir": str(project_dir),
        "work_dir": str(work_dir),
        "evidence_dir": str(evidence_dir),
        "reference_dir": str(reference_dir),
        "outputs": outputs,
        "steps": [
            {
                "name": "build-active-genome-index",
                "result": active_genome_index_result,
                "reason": "The VCF/gVCF is digitized into an Active Genome Index for targeted sample lookup.",
            }
        ],
        "warnings": [],
        "semantics": [
            "The VCF/gVCF source is digitized into a local Active Genome Index.",
            "Targeted rsID, locus, region, and exact allele lookup can use the Active Genome Index.",
            "Public evidence libraries are materialized lazily by focused tools such as ClinVar, HPO, GenCC, panel, and region annotation operations.",
        ],
    }


def _enqueue_reference_pass(active_genome_index_path: Path, parallel_workers: int | None) -> JsonObject | None:
    """Launch Phase B (the reference-block tail) as a detached background job.

    Reuses the standard job machinery (job_id, heartbeat, dead-worker
    detection, check_background_job polling). start_operation_job dedups on
    operation+params, so a duplicate parse_source call attaches to the running
    reference job instead of starting a second one. Best-effort: if the job
    can't be launched the variant-ready index is still fully usable, so we
    swallow the error rather than fail the parse.
    """
    from ...runtime import background_jobs

    params: JsonObject = {"active_genome_index_path": str(active_genome_index_path)}
    if parallel_workers is not None:
        params["parallel_workers"] = parallel_workers
    try:
        return background_jobs.start_operation_job("active_genome_index.build_reference_pass", params)
    except Exception:
        return None
