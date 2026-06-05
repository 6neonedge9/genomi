from __future__ import annotations

from pathlib import Path
from typing import Any

from ....active_genome_index.genotype_qc import (
    assess_genotype_support,
    assess_genotype_support_from_agi,
    assess_region_callability,
    assess_region_callability_from_agi,
    assess_sample_qc,
    assess_sample_qc_from_agi,
)
from ....active_genome_index.active_genome_index import (
    ActiveGenomeIndexNeed,
    default_agi_path,
    active_genome_index_summary,
    open_reader,
)
from ....active_genome_index.vcf import parse_region
from ....evidence import (
    default_evidence_path,
    evidence_summary,
    query_clinvar,
    query_population_frequency,
)
from ....runtime.handoff import attach_evidence_context, evidence_context
from ....runtime.paths import sample_slug_from_source

from ._helpers import (
    WORKFLOW_AREA_ID,
    default_static_outputs,
    workflow_contract,
)


def run_static_sample_qc(
    vcf: str | Path,
    *,
    evidence_db: str | Path | None = None,
    agi_path: str | Path | None = None,
    output: str | Path | None = None,
    genome_build: str = "auto",
    scan_records: int = 1000,
) -> dict[str, Any]:
    return assess_sample_qc(
        vcf,
        agi_path=agi_path,
        evidence_db=evidence_db,
        output=output,
        genome_build=genome_build,
        scan_records=scan_records,
    )


def run_static_sample_qc_from_agi(
    agi_path: str | Path,
    *,
    evidence_db: str | Path | None = None,
    output: str | Path | None = None,
    genome_build: str = "auto",
    scan_records: int = 1000,
) -> dict[str, Any]:
    return assess_sample_qc_from_agi(
        agi_path,
        evidence_db=evidence_db,
        output=output,
        genome_build=genome_build,
        scan_records=scan_records,
    )


def run_static_genotype_support(
    vcf: str | Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    evidence_db: str | Path | None = None,
    agi_path: str | Path | None = None,
    output: str | Path | None = None,
    genome_build: str = "auto",
    reference_fasta: str | Path | None = None,
    min_depth: int = 10,
    min_genotype_quality: int = 20,
) -> dict[str, Any]:
    return assess_genotype_support(
        vcf,
        chrom,
        pos,
        ref,
        alt,
        agi_path=agi_path,
        evidence_db=evidence_db,
        output=output,
        genome_build=genome_build,
        reference_fasta=reference_fasta,
        min_depth=min_depth,
        min_genotype_quality=min_genotype_quality,
    )


def run_static_genotype_support_from_agi(
    agi_path: str | Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    evidence_db: str | Path | None = None,
    output: str | Path | None = None,
    genome_build: str = "auto",
    reference_fasta: str | Path | None = None,
    min_depth: int = 10,
    min_genotype_quality: int = 20,
) -> dict[str, Any]:
    return assess_genotype_support_from_agi(
        agi_path,
        chrom,
        pos,
        ref,
        alt,
        evidence_db=evidence_db,
        output=output,
        genome_build=genome_build,
        reference_fasta=reference_fasta,
        min_depth=min_depth,
        min_genotype_quality=min_genotype_quality,
    )


def run_static_callability(
    vcf: str | Path,
    region: str,
    *,
    evidence_db: str | Path | None = None,
    agi_path: str | Path | None = None,
    output: str | Path | None = None,
    genome_build: str = "auto",
    min_depth: int = 10,
    min_covered_fraction: float = 0.95,
    limit: int = 5000,
) -> dict[str, Any]:
    return assess_region_callability(
        vcf,
        region,
        agi_path=agi_path,
        evidence_db=evidence_db,
        output=output,
        genome_build=genome_build,
        min_depth=min_depth,
        min_covered_fraction=min_covered_fraction,
        limit=limit,
    )


def run_static_callability_from_agi(
    agi_path: str | Path,
    region: str,
    *,
    evidence_db: str | Path | None = None,
    output: str | Path | None = None,
    genome_build: str = "auto",
    min_depth: int = 10,
    min_covered_fraction: float = 0.95,
    limit: int = 5000,
) -> dict[str, Any]:
    return assess_region_callability_from_agi(
        agi_path,
        region,
        evidence_db=evidence_db,
        output=output,
        genome_build=genome_build,
        min_depth=min_depth,
        min_covered_fraction=min_covered_fraction,
        limit=limit,
    )


def query_static_variant(
    vcf: str | Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    agi_path: str | Path | None = None,
    pass_only: bool = True,
    limit: int = 50,
) -> dict[str, Any]:
    reader = _static_reader(vcf, agi_path, need=ActiveGenomeIndexNeed.VARIANT)
    records = reader.query_variant(chrom, pos, ref, alt, pass_only=pass_only, limit=limit)
    return {
        "workflow_area": WORKFLOW_AREA_ID,
        "query": {"type": "variant", "chrom": chrom, "pos": pos, "ref": ref, "alt": alt},
        "count": len(records),
        "records": records,
        "evidence_context": evidence_context(
            "research",
            reason="Variant query output is local sample context for source-backed interpretation.",
            commands=[
                "genomi call variant.gather_allele_context --params '{\"db\":\"<evidence.sqlite>\",\"matches\":\"<clinvar.matches.jsonl>\",\"chrom\":\"<chrom>\",\"pos\":123,\"ref\":\"<ref>\",\"alt\":\"<alt>\"}'",
            ],
        ),
    }


def query_static_region(
    vcf: str | Path,
    region: str,
    *,
    agi_path: str | Path | None = None,
    variants_only: bool = False,
    pass_only: bool = True,
    limit: int = 200,
) -> dict[str, Any]:
    chrom, start, end = parse_region(region)
    reader = _static_reader(vcf, agi_path, need=ActiveGenomeIndexNeed.VARIANT)
    records = reader.query_region(
        chrom,
        start,
        end,
        variants_only=variants_only,
        pass_only=pass_only,
        limit=limit,
    )
    return {
        "workflow_area": WORKFLOW_AREA_ID,
        "query": {"type": "region", "region": region},
        "count": len(records),
        "records": records,
        "evidence_context": evidence_context(
            "research",
            reason="Region query output is local sample context for target selection and source-backed interpretation.",
            commands=[
                "genomi call variant.gather_allele_context --params '{\"db\":\"<evidence.sqlite>\",\"matches\":\"<clinvar.matches.jsonl>\",\"chrom\":\"<chrom>\",\"pos\":123,\"ref\":\"<ref>\",\"alt\":\"<alt>\"}'",
            ],
        ),
    }


def query_static_rsid(
    vcf: str | Path,
    rsid: str,
    *,
    agi_path: str | Path | None = None,
    pass_only: bool = True,
    limit: int = 50,
) -> dict[str, Any]:
    reader = _static_reader(vcf, agi_path, need=ActiveGenomeIndexNeed.VARIANT)
    records = reader.query_rsid(rsid, pass_only=pass_only, limit=limit)
    return {
        "workflow_area": WORKFLOW_AREA_ID,
        "query": {"type": "rsid", "rsid": rsid},
        "count": len(records),
        "records": records,
        "evidence_context": evidence_context(
            "research",
            reason="rsID query output is local sample context for interpretation or claim-status assessment.",
            commands=[
                "genomi call variant.gather_allele_context --params '{\"db\":\"<evidence.sqlite>\",\"matches\":\"<clinvar.matches.jsonl>\",\"chrom\":\"<chrom>\",\"pos\":123,\"ref\":\"<ref>\",\"alt\":\"<alt>\"}'",
                "genomi call research.build_target_packet --params '{\"db\":\"<evidence.sqlite>\",\"target_type\":\"variant\",\"chrom\":\"<chrom>\",\"pos\":123,\"ref\":\"<ref>\",\"alt\":\"<alt>\"}'",
            ],
        ),
    }


def query_static_coverage(
    vcf: str | Path,
    region: str,
    *,
    agi_path: str | Path | None = None,
    limit: int = 200,
) -> dict[str, Any]:
    chrom, start, end = parse_region(region)
    reader = _static_reader(vcf, agi_path, need=ActiveGenomeIndexNeed.REFERENCE)
    payload = reader.coverage(chrom, start, end, limit=limit)
    payload["workflow_area"] = WORKFLOW_AREA_ID
    return attach_evidence_context(
        payload,
        "research",
        reason="Coverage/callability context is a static input for target-specific claim assessment.",
        commands=[
            "genomi call variant.gather_gene_context --params '{\"db\":\"<evidence.sqlite>\",\"matches\":\"<clinvar.matches.jsonl>\",\"gene\":\"<gene>\"}'",
        ],
    )


def _static_reader(
    vcf: str | Path,
    agi_path: str | Path | None,
    *,
    need: ActiveGenomeIndexNeed,
):
    path = Path(agi_path) if agi_path else default_agi_path(vcf)
    return open_reader(path, need=need)


def static_db_lookup(
    evidence_db: str | Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    genome_build: str = "GRCh38",
) -> dict[str, Any]:
    return {
        "workflow_area": WORKFLOW_AREA_ID,
        "clinvar": query_clinvar(evidence_db, chrom, pos, ref, alt, genome_build=genome_build),
        "population": query_population_frequency(evidence_db, chrom, pos, ref, alt, genome_build=genome_build),
        "evidence_context": evidence_context(
            "research",
            reason="Static DB lookup output is structured evidence for source-backed interpretation.",
            commands=[
                "genomi call variant.gather_allele_context --params '{\"db\":\"<evidence.sqlite>\",\"matches\":\"<clinvar.matches.jsonl>\",\"chrom\":\"<chrom>\",\"pos\":123,\"ref\":\"<ref>\",\"alt\":\"<alt>\"}'",
            ],
        ),
    }


def summarize_static_state(
    vcf: str | Path,
    *,
    evidence_db: str | Path | None = None,
    agi_path: str | Path | None = None,
) -> dict[str, Any]:
    db_path = Path(evidence_db) if evidence_db is not None else default_evidence_path(vcf)
    agi_path = Path(agi_path) if agi_path is not None else default_agi_path(vcf)
    return {
        "workflow_area": WORKFLOW_AREA_ID,
        "contract": workflow_contract(),
        "active_genome_index": active_genome_index_summary(agi_path) if agi_path.exists() else None,
        "evidence": evidence_summary(db_path) if db_path.exists() else None,
        "outputs": default_static_outputs(vcf),
        "evidence_context": evidence_context(
            "research",
            reason="Static state is summarized for user-intent target research.",
            commands=["genomi call research.build_target_packet --params '{\"db\":\"<evidence.sqlite>\",\"target_type\":\"gene\",\"gene\":\"<gene>\"}'"],
        ),
    }


def summarize_static_state_from_agi(
    agi_path: str | Path,
    *,
    evidence_db: str | Path | None = None,
) -> dict[str, Any]:
    agi_path = Path(agi_path)
    db_path = Path(evidence_db) if evidence_db is not None else agi_path.parent / "evidence.sqlite"
    active_summary = active_genome_index_summary(agi_path) if agi_path.exists() else None
    evidence = evidence_summary(db_path) if db_path.exists() else None
    _align_evidence_metadata_to_active_index(evidence, active_summary)
    return {
        "workflow_area": WORKFLOW_AREA_ID,
        "contract": workflow_contract(),
        "active_genome_index": active_summary,
        "evidence": evidence,
        "outputs": {
            "sample_qc": str(agi_path.parent / "sample-qc.json"),
            "clinvar_matches": str(agi_path.parent / "clinvar.matches.jsonl"),
            "clinvar_scan": str(agi_path.parent / "clinvar.candidates.json"),
        },
        "evidence_context": evidence_context(
            "research",
            reason="Static state is summarized for user-intent target research.",
            commands=["genomi call research.build_target_packet --params '{\"db\":\"<evidence.sqlite>\",\"target_type\":\"gene\",\"gene\":\"<gene>\"}'"],
        ),
    }


def _align_evidence_metadata_to_active_index(evidence: dict[str, Any] | None, active_summary: dict[str, Any] | None) -> None:
    if not isinstance(evidence, dict) or not isinstance(active_summary, dict):
        return
    metadata = evidence.get("metadata")
    active_metadata = active_summary.get("metadata")
    if not isinstance(metadata, dict) or not isinstance(active_metadata, dict):
        return
    source = active_metadata.get("source")
    source_format = active_metadata.get("source_format")
    if source and source_format:
        metadata["run_sample_slug"] = sample_slug_from_source(source, source_format=str(source_format))
