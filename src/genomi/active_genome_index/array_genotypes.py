from __future__ import annotations

import json
from typing import Any

from .record_kinds import RECORD_KIND_ARRAY_CALL, RECORD_KIND_ARRAY_NO_CALL

JsonObject = dict[str, Any]

ARRAY_NO_CALLS = {"", ".", "--", "00", "NN"}
SUPPORTED_ARRAY_BASES = {"A", "C", "G", "T", "I", "D"}


def is_array_genotype_record(record: JsonObject) -> bool:
    return record.get("record_kind") in {RECORD_KIND_ARRAY_CALL, RECORD_KIND_ARRAY_NO_CALL}


def array_genotype_bases(record: JsonObject) -> list[str] | None:
    if record.get("record_kind") != RECORD_KIND_ARRAY_CALL:
        return None
    if str(record.get("filter") or "") not in {"PASS", "."}:
        return None
    bases = _observed_alleles(record.get("observed_alleles"))
    if not bases:
        return None
    if any(base not in SUPPORTED_ARRAY_BASES for base in bases):
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


def _observed_alleles(value: Any) -> list[str]:
    raw = value
    if isinstance(raw, str):
        try:
            raw = json.loads(raw)
        except json.JSONDecodeError:
            return []
    if not isinstance(raw, (list, tuple)):
        return []
    return [str(base).strip().upper() for base in raw if str(base).strip()]
