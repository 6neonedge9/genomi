from __future__ import annotations

from typing import Any

JsonObject = dict[str, Any]

ARRAY_FORMAT = "GT_ARRAY"
ARRAY_NO_CALLS = {"", ".", "--", "00", "NN"}
SUPPORTED_ARRAY_BASES = {"A", "C", "G", "T", "I", "D"}


def is_array_genotype_record(record: JsonObject) -> bool:
    value = record.get("format")
    if isinstance(value, list):
        return ARRAY_FORMAT in {str(item).upper() for item in value}
    return str(value or "").upper() == ARRAY_FORMAT


def array_genotype_bases(record: JsonObject) -> list[str] | None:
    if not is_array_genotype_record(record):
        return None
    if str(record.get("filter") or "") not in {"PASS", "."}:
        return None
    genotype = str(record.get("genotype") or record.get("sample") or "").strip().upper()
    if genotype in ARRAY_NO_CALLS:
        return None
    if any(separator in genotype for separator in ("/", "|", ",", ":", ";")):
        return None
    bases = list(genotype)
    if not bases or any(base not in SUPPORTED_ARRAY_BASES for base in bases):
        return None
    return bases


def called_genotype_tokens(genotype: Any) -> list[str]:
    text = str(genotype or "").strip().upper()
    if not text:
        return []
    if "/" in text or "|" in text:
        return [allele for allele in text.replace("|", "/").split("/") if allele not in {"", "."}]
    if len(text) in {1, 2} and all(base in SUPPORTED_ARRAY_BASES for base in text):
        return list(text)
    return [text]


def count_array_allele(
    record: JsonObject,
    *,
    target_allele: str,
    allowed_alleles: list[str] | tuple[str, ...] | set[str] | None = None,
) -> JsonObject:
    target = _single_array_allele(target_allele)
    if target is None:
        return {"status": "missing", "reason": "array_target_allele_not_single_base"}
    allowed: set[str] | None = None
    if allowed_alleles is not None:
        allowed = set()
        for allele in allowed_alleles:
            normalized = _single_array_allele(allele)
            if normalized is None:
                return {"status": "missing", "reason": "array_allele_model_not_single_base"}
            allowed.add(normalized)
        allowed.add(target)

    bases = array_genotype_bases(record)
    if bases is None:
        return {"status": "missing", "reason": "missing_genotype"}
    if allowed is not None and any(base not in allowed for base in bases):
        return {
            "status": "missing",
            "reason": "array_genotype_allele_outside_allowed_alleles",
            "allele_bases": bases,
        }
    return {
        "status": "matched",
        "allele_bases": bases,
        "dosage": float(sum(1 for base in bases if base == target)),
        "ploidy": len(bases),
    }


def _single_array_allele(value: str | None) -> str | None:
    allele = str(value or "").strip().upper()
    if len(allele) != 1 or allele not in SUPPORTED_ARRAY_BASES:
        return None
    return allele
