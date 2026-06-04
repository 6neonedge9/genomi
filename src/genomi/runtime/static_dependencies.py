"""Genome-build resolution utilities.

ClinVar/reference-FASTA materialization, freshness tracking, and the FASTA
gunzip+faidx transform now live in the central library manager
(``genomi.runtime.libraries``); only build inference/normalization remains here.
"""

from __future__ import annotations

from pathlib import Path

from ..active_genome_index.vcf import VcfHeader, read_header


def resolve_genome_build(vcf: str | Path, requested: str | None) -> str:
    requested_normalized = (requested or "auto").strip()
    if requested_normalized.lower() not in {"", "auto"}:
        return _normalize_genome_build(requested_normalized)
    return infer_genome_build_from_vcf(vcf) or "GRCh38"


def infer_genome_build_from_vcf(vcf: str | Path) -> str | None:
    header = _read_vcf_header_for_build_inference(vcf)
    if header is None:
        return None
    return _infer_genome_build_from_header(header)


def _read_vcf_header_for_build_inference(vcf: str | Path) -> VcfHeader | None:
    try:
        return read_header(vcf)
    except Exception:
        pass
    try:
        from ..active_genome_index.source_intake.text_io import open_genomic_binary

        meta: list[str] = []
        with open_genomic_binary(Path(vcf)) as handle:
            while True:
                line = handle.readline()
                if not line:
                    return None
                text = line.decode("utf-8", errors="replace").rstrip("\r\n")
                if text.startswith("##"):
                    meta.append(text)
                    continue
                if text.startswith("#CHROM"):
                    return VcfHeader(meta=meta, columns=text.split("\t"))
                return None
    except Exception:
        return None


def _infer_genome_build_from_header(header: VcfHeader) -> str | None:
    text = " ".join(
        value or ""
        for value in [
            header.first_meta_value("reference"),
            header.first_meta_value("referenceInfo"),
            header.first_meta_value("assembly"),
        ]
    ).lower()
    if any(token in text for token in ["grch37", "hg19", "g1k.37", "b37"]):
        return "GRCh37"
    if any(token in text for token in ["grch38", "hg38", "grch38.p", "b38"]):
        return "GRCh38"
    contigs = header.contigs()
    if contigs and any(contig.startswith("chr") for contig in contigs[:24]):
        return "GRCh38"
    return None


def _normalize_genome_build(value: str) -> str:
    normalized = value.strip().lower()
    if normalized in {"grch37", "hg19", "37"}:
        return "GRCh37"
    if normalized in {"grch38", "hg38", "38"}:
        return "GRCh38"
    raise ValueError(f"unsupported genome build for static dependencies: {value}")
