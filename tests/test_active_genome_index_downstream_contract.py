from __future__ import annotations

import os
import tempfile
from pathlib import Path
from unittest import mock

from genomi.active_genome_index.active_genome_index import active_genome_index_readiness
from genomi.evidence import import_clinvar_vcf
from genomi.operations import call_operation

from _active_genome_index_contract_fixtures import (
    EXPECTED_CLINVAR_MATCHED_ALLELES,
    EXPECTED_RAW_SCORE,
    LOCUS_MODEL,
    ActiveGenomeIndexContractFixtureMixin,
)
from _genomi_runtime_helpers import GenomiRuntimeTestCase


class ActiveGenomeIndexDownstreamContractTests(
    ActiveGenomeIndexContractFixtureMixin,
    GenomiRuntimeTestCase,
):
    """PGP-HMS-shaped fake sources must feed every coordinate consumer.

    Public PGP-HMS downloads are used as the source of truth for wrappers,
    member names, comments, and columns. The genotype rows here are synthetic
    and deliberately tiny.
    """

    def test_pgp_hms_shaped_supported_sources_feed_coordinate_consumers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            previous = os.getcwd()
            os.chdir(tmp)
            try:
                scoring_file = self._write_scoring_file(Path("PGSAGI001_hmPOS_GRCh37.txt"))
                imported_score = call_operation(
                    "prs.import_scoring_file",
                    {
                        "pgs_id": "PGSAGI001",
                        "scoring_file": str(scoring_file),
                        "genome_build": "GRCh37",
                        "force": True,
                    },
                )
                self.assertEqual(imported_score["status"], "completed")
                self._install_contract_ancestry_panel()

                clinvar_db = Path("contract-clinvar.sqlite")
                clinvar_vcf = self._write_clinvar_fixture(Path("contract.clinvar.vcf"))
                import_clinvar_vcf(clinvar_vcf, clinvar_db, source_version="contract-fixture", genome_build="GRCh37")

                for case_id, expected_format, writer in self._source_cases():
                    with self.subTest(source=case_id):
                        source = writer(Path(case_id))
                        self._assert_source_feeds_coordinate_consumers(
                            source,
                            expected_format=expected_format,
                            case_id=case_id,
                            imported_score=imported_score,
                            clinvar_db=clinvar_db,
                        )

                reference = self._write_reference_fasta(Path("contract-reference.fa"))
                with self.subTest(source="bam"):
                    bam = self._write_bam_source(Path("Nebula_Genomics_BAM_format.bam"))
                    with self._mock_derived_vcf_materialization():
                        self._assert_source_feeds_coordinate_consumers(
                            bam,
                            expected_format="bam",
                            case_id="bam",
                            imported_score=imported_score,
                            clinvar_db=clinvar_db,
                            parse_overrides={"reference_fasta": str(reference)},
                        )

                with self.subTest(source="fastq"):
                    fastq = self._write_fastq_sources(Path("60820188475559_SA_L001_R1_001.fastq.gz"))
                    with self._mock_derived_vcf_materialization(), mock.patch(
                        "genomi.active_genome_index.source_intake.sequencing.align_fastq_to_bam",
                        side_effect=self._fake_align_fastq_to_bam,
                    ):
                        self._assert_source_feeds_coordinate_consumers(
                            fastq,
                            expected_format="fastq",
                            case_id="fastq",
                            imported_score=imported_score,
                            clinvar_db=clinvar_db,
                            parse_overrides={"reference_fasta": str(reference)},
                        )
            finally:
                os.chdir(previous)

    def _assert_source_feeds_coordinate_consumers(
        self,
        source: Path,
        *,
        expected_format: str,
        case_id: str,
        imported_score: dict[str, object],
        clinvar_db: Path,
        parse_overrides: dict[str, object] | None = None,
    ) -> None:
        parse_params = {"source": str(source), "genome_build": "GRCh37", "force": True}
        parse_params.update(parse_overrides or {})
        parsed = call_operation("genomi.parse_source", parse_params)

        self.assertEqual(parsed["status"], "completed")
        self.assertEqual(parsed["source_format"], expected_format)
        readiness = active_genome_index_readiness(parsed["outputs"]["active_genome_index_path"])
        self.assertTrue(readiness["complete"], readiness)
        self.assertEqual(readiness["missing_objects"], [])

        variant = call_operation("variant.resolve", {"rsid": "rsagi2", "genome_build": "GRCh37"})
        self.assertEqual(variant["sample_context"]["count"], 1, variant)
        self.assertEqual(variant["sample_context"]["matches"][0]["genotype"], self._expected_genotype_for_source(expected_format, 1))

        summary = call_operation("active_genome_index.summarize")
        self.assertTrue(summary["active_genome_index"]["active_genome_index_readiness"]["complete"], summary)
        self.assertGreaterEqual(summary["active_genome_index"]["stats"]["total_records"], len(LOCUS_MODEL))

        callset_qc = call_operation(
            "active_genome_index.classify_callset_qc",
            {"genome_build": "GRCh37", "scan_records": 100},
        )
        self.assertEqual(callset_qc["status"], "completed", callset_qc)
        self.assertGreaterEqual(callset_qc["summary"]["total_records"], len(LOCUS_MODEL))

        callability = call_operation(
            "active_genome_index.classify_region_callability",
            {
                "region": "1:100-100",
                "genome_build": "GRCh37",
                "min_covered_fraction": 0.1,
            },
        )
        self.assertEqual(callability["status"], "completed", callability)
        self.assertIn(
            callability["callability_status"],
            {
                "callable",
                "not_callable",
                "unknown",
                "unknown_missing_depth",
                "unknown_no_coverage",
                "unknown_no_reference_blocks",
            },
        )

        with self._tiny_prs_thresholds():
            overlap_result = call_operation(
                "prs.check_score_overlap",
                {
                    "score_dir": imported_score["score_cache"]["score_dir"],
                    "genome_build": "GRCh37",
                },
            )
            prs_result = call_operation(
                "prs.calculate_score",
                {
                    "score_dir": imported_score["score_cache"]["score_dir"],
                    "genome_build": "GRCh37",
                },
            )
        self.assertEqual(overlap_result["status"], "score_ready", overlap_result)
        self.assertEqual(overlap_result["sample_qc"]["matched_variant_count"], len(LOCUS_MODEL))
        self.assertEqual(prs_result["status"], "completed", prs_result)
        self.assertEqual(prs_result["sample_qc"]["matched_variant_count"], len(LOCUS_MODEL))
        self.assertEqual(prs_result["sample_qc"]["missing_variant_count"], 0)
        self.assertAlmostEqual(prs_result["score_result"]["raw_weighted_score"], EXPECTED_RAW_SCORE)

        ancestry_result = call_operation("ancestry.check_sample_overlap", {"genome_build": "GRCh37"})
        self.assertEqual(ancestry_result["status"], "completed", ancestry_result)
        self.assertEqual(ancestry_result["sample_qc"]["usable_marker_count"], len(LOCUS_MODEL))
        self.assertEqual(ancestry_result["sample_qc"]["missing_marker_count"], 0)

        matches_path = Path(f"{case_id}.clinvar.matches.jsonl")
        clinvar_result = call_operation(
            "clinvar.match_variants",
            {
                "db": str(clinvar_db),
                "output": str(matches_path),
                "genome_build": "GRCh37",
                "force": True,
            },
        )
        self.assertEqual(clinvar_result["status"], "completed", clinvar_result)
        self.assertEqual(clinvar_result["stats"]["matched_alleles"], EXPECTED_CLINVAR_MATCHED_ALLELES)
        self.assertEqual(clinvar_result["stats"]["written_records"], EXPECTED_CLINVAR_MATCHED_ALLELES)
        self._assert_clinvar_payloads_are_real_alleles(matches_path, expected_format=expected_format)

        scanned = call_operation(
            "clinvar.scan_candidates",
            {
                "matches": str(matches_path),
                "db": str(clinvar_db),
                "output": str(Path(f"{case_id}.clinvar.candidates.json")),
                "genome_build": "GRCh37",
                "force": True,
            },
        )
        self.assertEqual(scanned["status"], "completed", scanned)
        candidates_by_pos = {int(candidate["variant"]["pos"]): candidate for candidate in scanned["candidate_inventory"]}
        self.assertIn(200, candidates_by_pos)
        self.assertIn("heterozygous_p_lp_context_needed", candidates_by_pos[200]["buckets"])
        if expected_format in {"vcf", "gvcf", "bam", "fastq"}:
            self.assertEqual(candidates_by_pos[200]["match_provenance"]["primary_match_basis"], "exact_allele")
        else:
            self.assertEqual(
                candidates_by_pos[200]["match_provenance"]["primary_match_basis"],
                "consumer_array_allele_inference",
            )
            self.assertEqual(candidates_by_pos[200]["variant"]["source_record_ref"], "N")
            self.assertEqual(candidates_by_pos[200]["variant"]["source_record_format"], "GT_ARRAY")

        observed_support = self._genotype_support(chrom="1", pos=100, ref="A", alt="C")
        self.assertEqual(observed_support["sample_observation"]["target_alt_observed"], True)
        self.assertEqual(observed_support["sample_observation"]["alt_allele_count"], 2)

        heterozygous_support = self._genotype_support(chrom="1", pos=200, ref="T", alt="G")
        self.assertEqual(heterozygous_support["sample_observation"]["target_alt_observed"], True)
        self.assertEqual(heterozygous_support["sample_observation"]["alt_allele_count"], 1)
        self.assertEqual(heterozygous_support["sample_observation"]["zygosity"], "heterozygous")

        absent_support = self._genotype_support(chrom="1", pos=300, ref="A", alt="G")
        self.assertEqual(absent_support["support_status"], "not_observed")
        self.assertEqual(absent_support["sample_observation"]["target_alt_observed"], False)
        self.assertEqual(absent_support["sample_observation"]["alt_allele_count"], 0)
