from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any

from ..clinvar_match_model import (
    MATCH_BASIS_CONSUMER_ARRAY_ALLELE_INFERENCE,
    MATCH_BASIS_EXACT_ALLELE,
    MATCH_BASIS_MULTIALLELIC_ALT,
)
from .record_kinds import RECORD_KIND_ARRAY_CALL, RECORD_KIND_VARIANT_CALL
from .reader import ActiveGenomeIndexReader


JsonObject = dict[str, Any]


def stage_clinvar_match_records(
    reader: ActiveGenomeIndexReader,
    connection: sqlite3.Connection,
) -> JsonObject:
    """Materialize the AGI rows ClinVar matching may read into temp tables."""

    reader.ensure_ready()
    alias = "_agi_clinvar_records"
    reader.attach_to(connection, alias)
    try:
        _ensure_ready_for_clinvar_match(connection, alias, reader.agi_path)
        connection.executescript(
            """
            drop table if exists temp.selected_active_genome_index_agi_records;
            drop table if exists temp.selected_active_genome_index_records;
            create temp table selected_active_genome_index_agi_records (
                record_rowid integer not null,
                chrom text not null,
                chrom_sort integer,
                pos integer not null,
                rsid text,
                ref text,
                alt text,
                qual text,
                filter text,
                info text,
                sample_index integer,
                sample_name text,
                format text,
                genotype text,
                depth integer,
                genotype_quality integer,
                record_kind text,
                observed_alleles text
            );
            create temp table selected_active_genome_index_records (
                record_rowid integer not null,
                chrom text not null,
                chrom_sort integer,
                pos integer not null,
                rsid text,
                ref text,
                alt text,
                qual text,
                filter text,
                info text,
                sample_index integer,
                sample_name text,
                format text,
                genotype text,
                depth integer,
                genotype_quality integer,
                record_kind text,
                observed_alleles text,
                clinvar_batch_id text not null,
                clinvar_match_mode text not null,
                clinvar_match_alt text not null
            );
            create index selected_active_genome_index_agi_records_order_idx
                on selected_active_genome_index_agi_records(chrom_sort, pos, record_rowid, sample_index);
            create index selected_active_genome_index_records_locus_idx
                on selected_active_genome_index_records(chrom, pos);
            """
        )
        sql = f"""
            insert into temp.selected_active_genome_index_agi_records (
                record_rowid, chrom, chrom_sort, pos, rsid, ref, alt, qual, filter,
                info, sample_index, sample_name, format, genotype, depth,
                genotype_quality, record_kind, observed_alleles
            )
            select rowid, chrom, chrom_sort, pos, rsid, ref, alt, qual, filter,
                   info, sample_index, sample_name, format, genotype, depth,
                   genotype_quality, record_kind, observed_alleles
            from {alias}.records
            where record_kind in ('{RECORD_KIND_VARIANT_CALL}', '{RECORD_KIND_ARRAY_CALL}')
        """
        connection.execute(sql)
        connection.execute(
            """
            insert into temp.selected_active_genome_index_records (
                record_rowid, chrom, chrom_sort, pos, rsid, ref, alt, qual, filter,
                info, sample_index, sample_name, format, genotype, depth,
                genotype_quality, record_kind, observed_alleles,
                clinvar_batch_id, clinvar_match_mode, clinvar_match_alt
            )
            select distinct
                   r.record_rowid, r.chrom, r.chrom_sort, r.pos, r.rsid, r.ref, r.alt,
                   r.qual, r.filter, r.info, r.sample_index, r.sample_name, r.format,
                   r.genotype, r.depth, r.genotype_quality, r.record_kind, r.observed_alleles,
                   case
                       when instr(r.alt, ',') > 0
                           then cast(r.record_rowid as text) || ':' || upper(observed.value)
                       else cast(r.record_rowid as text)
                   end as clinvar_batch_id,
                   case
                       when instr(r.alt, ',') > 0 then ?
                       else ?
                   end as clinvar_match_mode,
                   upper(observed.value) as clinvar_match_alt
            from temp.selected_active_genome_index_agi_records as r
            join json_each(r.observed_alleles) as observed
            where r.record_kind = ?
              and r.alt not in ('', '.')
              and upper(observed.value) <> upper(r.ref)
              and instr(',' || upper(r.alt) || ',', ',' || upper(observed.value) || ',') > 0
            """,
            (
                MATCH_BASIS_MULTIALLELIC_ALT,
                MATCH_BASIS_EXACT_ALLELE,
                RECORD_KIND_VARIANT_CALL,
            ),
        )
        connection.execute(
            """
            insert into temp.selected_active_genome_index_records (
                record_rowid, chrom, chrom_sort, pos, rsid, ref, alt, qual, filter,
                info, sample_index, sample_name, format, genotype, depth,
                genotype_quality, record_kind, observed_alleles,
                clinvar_batch_id, clinvar_match_mode, clinvar_match_alt
            )
            select distinct
                   r.record_rowid, r.chrom, r.chrom_sort, r.pos, r.rsid, r.ref, r.alt,
                   r.qual, r.filter, r.info, r.sample_index, r.sample_name, r.format,
                   r.genotype, r.depth, r.genotype_quality, r.record_kind, r.observed_alleles,
                   cast(r.record_rowid as text) || ':' || upper(observed.value) || ':array'
                       as clinvar_batch_id,
                   ? as clinvar_match_mode,
                   upper(observed.value) as clinvar_match_alt
            from temp.selected_active_genome_index_agi_records as r
            join json_each(r.observed_alleles) as observed
            where r.record_kind = ?
              and upper(observed.value) in ('A', 'C', 'G', 'T')
            """,
            (
                MATCH_BASIS_CONSUMER_ARRAY_ALLELE_INFERENCE,
                RECORD_KIND_ARRAY_CALL,
            ),
        )
        connection.commit()
        return {"source_format": _source_format(connection, alias)}
    finally:
        if connection.in_transaction:
            connection.rollback()
        connection.execute(f"detach database {alias}")


def _ensure_ready_for_clinvar_match(connection: sqlite3.Connection, alias: str, agi_path: Path) -> None:
    stats_count = connection.execute(f"select count(*) from {alias}.stats").fetchone()[0]
    index_names = {
        str(row["name"])
        for row in connection.execute(
            f"""
            select name
            from {alias}.sqlite_master
            where type = 'index' and tbl_name = 'records'
            """
        )
    }
    required_indexes = {"records_export_idx", "records_variant_idx"}
    missing_indexes = sorted(required_indexes - index_names)
    if stats_count == 0 or missing_indexes:
        details = []
        if stats_count == 0:
            details.append("missing stats rows")
        if missing_indexes:
            details.append(f"missing query indexes: {', '.join(missing_indexes)}")
        raise RuntimeError(
            f"Active Genome Index is incomplete for ClinVar refresh ({agi_path}): "
            f"{'; '.join(details)}. Rebuild the Active Genome Index from the source genome file once."
        )


def _source_format(connection: sqlite3.Connection, alias: str) -> str | None:
    try:
        row = connection.execute(f"select value from {alias}.metadata where key = 'source_format'").fetchone()
    except sqlite3.Error:
        return None
    if row is None:
        return None
    try:
        parsed = json.loads(str(row["value"]))
    except (TypeError, json.JSONDecodeError):
        parsed = row["value"]
    return str(parsed) if parsed else None
