from __future__ import annotations

from typing import Any

PASSING_FILTER_VALUES = frozenset({"PASS", "."})


def passing_filter_sql(column: str = "filter") -> str:
    return f"{column} in ('PASS', '.')"


def is_passing_filter(value: Any) -> bool:
    return str(value or "") in PASSING_FILTER_VALUES
