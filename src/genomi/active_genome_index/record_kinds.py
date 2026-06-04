from __future__ import annotations

from typing import Any

JsonObject = dict[str, Any]

ARRAY_FORMAT = "GT_ARRAY"
ARRAY_NO_CALL_FILTER = "NO_CALL"

RECORD_KIND_VARIANT_CALL = "variant_call"
RECORD_KIND_REFERENCE_BLOCK = "reference_block"
RECORD_KIND_NO_CALL = "no_call"
RECORD_KIND_ARRAY_CALL = "array_call"
RECORD_KIND_ARRAY_NO_CALL = "array_no_call"


def array_record_kind(*, is_called: bool) -> str:
    return RECORD_KIND_ARRAY_CALL if is_called else RECORD_KIND_ARRAY_NO_CALL


def reference_block_sql(alias: str | None = None) -> str:
    prefix = f"{alias}." if alias else ""
    return f"{prefix}record_kind = '{RECORD_KIND_REFERENCE_BLOCK}'"


def array_no_call_sql(alias: str | None = None) -> str:
    prefix = f"{alias}." if alias else ""
    return f"{prefix}record_kind = '{RECORD_KIND_ARRAY_NO_CALL}'"


def is_reference_block_record(record: JsonObject) -> bool:
    return record.get("record_kind") == RECORD_KIND_REFERENCE_BLOCK


def _is_no_call_genotype(value: Any) -> bool:
    genotype = str(value or "").strip()
    if not genotype:
        return True
    return "." in genotype
