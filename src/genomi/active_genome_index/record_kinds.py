from __future__ import annotations

import json
from typing import Any

JsonObject = dict[str, Any]

ARRAY_FORMAT = "GT_ARRAY"
ARRAY_NO_CALL_FILTER = "NO_CALL"
ARRAY_RECORD_KIND_VERSION = 1

RECORD_KIND_VARIANT_CALL = "variant_call"
RECORD_KIND_REFERENCE_BLOCK = "reference_block"
RECORD_KIND_NO_CALL = "no_call"
RECORD_KIND_ARRAY_CALL = "array_call"
RECORD_KIND_ARRAY_NO_CALL = "array_no_call"


def array_record_kind(*, is_called: bool) -> str:
    return RECORD_KIND_ARRAY_CALL if is_called else RECORD_KIND_ARRAY_NO_CALL


def array_record_info(*, source_format: str, is_called: bool) -> JsonObject:
    return {
        "source_format": source_format,
        "coordinate_semantics": "plus_strand_grch37",
        "record_kind": array_record_kind(is_called=is_called),
        "call_status": "called" if is_called else "no_call",
        "record_kind_version": ARRAY_RECORD_KIND_VERSION,
    }


def is_array_record_format(value: Any) -> bool:
    if isinstance(value, list):
        return ARRAY_FORMAT in {str(item).upper() for item in value}
    return str(value or "").upper() == ARRAY_FORMAT


def reference_block_sql(alias: str | None = None) -> str:
    prefix = f"{alias}." if alias else ""
    genotype = f"coalesce({prefix}genotype, '')"
    return (
        f"{prefix}is_variant = 0 "
        f"and coalesce({prefix}format, '') != '{ARRAY_FORMAT}' "
        f"and {genotype} != '' "
        f"and {genotype} not like '%.%'"
    )


def array_no_call_sql(alias: str | None = None) -> str:
    prefix = f"{alias}." if alias else ""
    return f"coalesce({prefix}format, '') = '{ARRAY_FORMAT}' and {prefix}filter = '{ARRAY_NO_CALL_FILTER}'"


def infer_record_kind(
    *,
    format_value: Any,
    is_variant: bool,
    filter_value: Any,
    info_raw: Any = None,
    genotype_value: Any = None,
) -> str:
    if not is_array_record_format(format_value):
        if _is_no_call_genotype(genotype_value):
            return RECORD_KIND_NO_CALL
        return RECORD_KIND_VARIANT_CALL if is_variant else RECORD_KIND_REFERENCE_BLOCK
    info_kind = _array_record_kind_from_info(info_raw)
    if info_kind in {RECORD_KIND_ARRAY_CALL, RECORD_KIND_ARRAY_NO_CALL}:
        return info_kind
    return RECORD_KIND_ARRAY_NO_CALL if str(filter_value or "") == ARRAY_NO_CALL_FILTER else RECORD_KIND_ARRAY_CALL


def is_reference_block_record(record: JsonObject) -> bool:
    kind = record.get("record_kind")
    if kind:
        return kind == RECORD_KIND_REFERENCE_BLOCK
    return (
        not bool(record.get("is_variant"))
        and not is_array_record_format(record.get("format"))
        and not _is_no_call_genotype(record.get("genotype"))
    )


def _is_no_call_genotype(value: Any) -> bool:
    genotype = str(value or "").strip()
    if not genotype:
        return True
    return "." in genotype


def _array_record_kind_from_info(info_raw: Any) -> str | None:
    if isinstance(info_raw, dict):
        kind = info_raw.get("record_kind")
        return str(kind) if kind is not None else None
    if not info_raw:
        return None
    try:
        parsed = json.loads(str(info_raw))
    except (TypeError, json.JSONDecodeError):
        return None
    if not isinstance(parsed, dict):
        return None
    kind = parsed.get("record_kind")
    return str(kind) if kind is not None else None
