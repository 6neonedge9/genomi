from __future__ import annotations

import json
import sqlite3
from unittest import mock

from genomi.capabilities.prs import scorer as prs_scorer


def memory_prs_index() -> sqlite3.Connection:
    connection = sqlite3.connect(":memory:")
    connection.row_factory = sqlite3.Row
    connection.executescript(
        """
        create table records (
            chrom text not null,
            chrom_sort integer not null,
            pos integer not null,
            end integer not null,
            ref text not null,
            alt text not null,
            filter text not null,
            format text,
            genotype text not null,
            record_kind text not null,
            observed_alleles text,
            offset integer not null,
            sample_index integer not null
        );
        create table spans (
            chrom text not null,
            chrom_sort integer not null,
            pos integer not null,
            end integer not null,
            offset integer not null,
            sample_index integer not null
        );
        """
    )
    return connection


def insert_prs_record(
    connection: sqlite3.Connection,
    *,
    pos: int,
    ref: str,
    alt: str,
    genotype: str,
) -> None:
    record_kind, observed_alleles = vcf_record_observation(ref=ref, alt=alt, genotype=genotype)
    _insert_record(
        connection,
        pos=pos,
        ref=ref,
        alt=alt,
        genotype=genotype,
        record_kind=record_kind,
        observed_alleles=observed_alleles,
        format_value="GT",
        filter_value="PASS",
    )


def insert_array_prs_record(
    connection: sqlite3.Connection,
    *,
    pos: int,
    genotype: str,
    filter_value: str = "PASS",
) -> None:
    is_called = filter_value == "PASS" and genotype not in {"", ".", "--", "00", "NN"}
    _insert_record(
        connection,
        pos=pos,
        ref=".",
        alt=".",
        genotype=genotype,
        record_kind="array_call" if is_called else "array_no_call",
        observed_alleles=list(genotype) if is_called else None,
        format_value="GT_ARRAY",
        filter_value=filter_value,
    )


def vcf_record_observation(*, ref: str, alt: str, genotype: str) -> tuple[str, list[str] | None]:
    tokens = [token for token in genotype.replace("|", "/").split("/") if token]
    if not tokens or any(token == "." for token in tokens):
        return "no_call", None
    alts = [] if alt == "." else alt.split(",")
    observed: list[str] = []
    for token in tokens:
        if token == "0":
            observed.append(ref)
            continue
        try:
            observed.append(alts[int(token) - 1])
        except (IndexError, ValueError):
            return "no_call", None
    record_kind = "variant_call" if any(token != "0" for token in tokens) and alts else "reference_block"
    return record_kind, observed


def score_variant(*, pos: int, effect_allele: str, other_allele: str) -> dict[str, object]:
    return {
        "variant_index": 0,
        "variant_id": f"1:{pos}:{other_allele}:{effect_allele}",
        "rsid": "rs-test",
        "chrom": "1",
        "pos": pos,
        "effect_allele": effect_allele,
        "other_allele": other_allele,
        "effect_weight": 1.0,
        "harmonized": True,
        "palindromic": False,
    }


def tiny_thresholds(*, min_variants: int = 1, min_fraction: float = 0.10):
    return mock.patch.multiple(
        prs_scorer,
        MIN_SCORE_VARIANTS=min_variants,
        MIN_OVERLAP_FRACTION=min_fraction,
        MODERATE_OVERLAP_FRACTION=0.50,
        HIGH_OVERLAP_FRACTION=0.90,
    )


def _insert_record(
    connection: sqlite3.Connection,
    *,
    pos: int,
    ref: str,
    alt: str,
    genotype: str,
    record_kind: str,
    observed_alleles: list[str] | None,
    format_value: str,
    filter_value: str,
) -> None:
    connection.execute(
        """
        insert into records(
            chrom, chrom_sort, pos, end, ref, alt, filter, format, genotype,
            record_kind, observed_alleles, offset, sample_index
        )
        values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            "1",
            1,
            pos,
            pos,
            ref,
            alt,
            filter_value,
            format_value,
            genotype,
            record_kind,
            json.dumps(observed_alleles, sort_keys=True) if observed_alleles is not None else None,
            pos,
            0,
        ),
    )
