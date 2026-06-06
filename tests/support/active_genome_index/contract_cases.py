from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from genomi.active_genome_index.source_intake.arrays import SUPPORTED_CONSUMER_ARRAY_FORMATS

from tests.support.active_genome_index.contract_fixtures import (
    EXPECTED_CLINVAR_MATCHED_ALLELES,
    LOCUS_MODEL,
)


@dataclass(frozen=True)
class SourceContractCase:
    case_id: str
    expected_format: str
    writer: Callable[[Path], Path]
    parse_overrides: dict[str, object] | None = None

    @property
    def is_consumer_array(self) -> bool:
        return self.expected_format in SUPPORTED_CONSUMER_ARRAY_FORMATS

    @property
    def expected_record_stats(self) -> dict[str, int]:
        if self.is_consumer_array:
            return {
                "total_records": len(LOCUS_MODEL),
                "variant_records": 0,
                "reference_records": 0,
                "pass_records": len(LOCUS_MODEL),
                "fail_records": 0,
            }
        if self.expected_format == "gvcf":
            return {
                "total_records": len(LOCUS_MODEL) + 1,
                "variant_records": 2,
                "reference_records": 3,
                "pass_records": len(LOCUS_MODEL) + 1,
                "fail_records": 0,
            }
        return {
            "total_records": len(LOCUS_MODEL),
            "variant_records": 2,
            "reference_records": 2,
            "pass_records": len(LOCUS_MODEL),
            "fail_records": 0,
        }

    @property
    def expected_callability_for_called_site(self) -> str:
        if self.is_consumer_array:
            return "unknown_no_reference_blocks"
        return "unknown_missing_depth"

    @property
    def expected_callability_for_unrepresented_site(self) -> str:
        if self.is_consumer_array:
            return "unknown_no_reference_blocks"
        return "not_callable"

    @property
    def expected_clinvar_scanned_records(self) -> int:
        if self.is_consumer_array:
            return len(LOCUS_MODEL)
        return EXPECTED_CLINVAR_MATCHED_ALLELES

    @property
    def expected_clinvar_queried_alleles(self) -> int:
        if self.is_consumer_array:
            return sum(len(set(str(locus["bases"]))) for locus in LOCUS_MODEL)
        return EXPECTED_CLINVAR_MATCHED_ALLELES

    @property
    def expected_source_kind(self) -> str:
        return {
            "bam": "alignment_reads",
            "fastq": "paired_reads_input",
        }.get(
            self.expected_format,
            "consumer_genotype_array" if self.is_consumer_array else "variant_callset",
        )


@dataclass(frozen=True)
class LocusContract:
    rsid: str
    chrom: str
    pos: int
    ref: str
    alt: str
    expected_alt_observed: bool
    expected_alt_count: int | None
    expected_zygosity: str


LOCUS_CONTRACTS = (
    LocusContract("rs900000001", "1", 100, "A", "C", True, 2, "homozygous_alternate"),
    LocusContract("rs900000002", "1", 200, "T", "G", True, 1, "heterozygous"),
    LocusContract("rs900000003", "1", 300, "A", "G", False, 0, "reference_or_other_alternate"),
    LocusContract("rs900000004", "1", 400, "C", "T", False, 0, "reference_or_other_alternate"),
)

UNREPRESENTED_LOCUS = LocusContract(
    rsid="rs900000099",
    chrom="1",
    pos=500,
    ref="A",
    alt="G",
    expected_alt_observed=False,
    expected_alt_count=None,
    expected_zygosity="unknown",
)
