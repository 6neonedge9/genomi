from __future__ import annotations

import json
import re
from typing import Any

JsonObject = dict[str, Any]


def observed_alleles_from_record(record: JsonObject) -> list[str]:
    """Return called allele bases from AGI observation fields or numeric GT."""
    observed = _observed_alleles(record.get("observed_alleles"))
    if observed:
        return observed
    genotype = str(record.get("genotype") or "").strip().upper()
    if not re.fullmatch(r"[0-9.]+([/|][0-9.]+)*", genotype):
        return []
    ref = str(record.get("ref") or "").strip().upper()
    alts = [item.strip().upper() for item in str(record.get("alt") or "").split(",") if item.strip()]
    alleles: list[str] = []
    for token in re.split(r"[/|]", genotype):
        if token in {"", "."}:
            continue
        try:
            index = int(token)
        except ValueError:
            continue
        if index == 0 and ref:
            alleles.append(ref)
        elif 0 < index <= len(alts):
            alleles.append(alts[index - 1])
    return alleles


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
