from __future__ import annotations

import unittest
from pathlib import Path
from unittest import mock

from genomi.active_genome_index.active_genome_index import create_active_genome_index, default_agi_path
from genomi.capabilities.prs import harmonize as prs_harmonize
from genomi.capabilities.prs import scorer as prs_scorer
from genomi.operations import call_operation
from genomi.runtime import context as runtime_context

from tests.support.capabilities.prs_contract import PolygenicScoreTestBase


class PolygenicScoreCalculationTests(PolygenicScoreTestBase):
    def test_calibration_uses_only_supplied_parameters(self) -> None:
        scoring_file = self._write_scoring_file()
        call_operation("prs.import_scoring_file", {"pgs_id": "PGS900001", "scoring_file": str(scoring_file)})
        vcf = self._write_indexed_vcf("sample_calibrated.vcf")
        self._select_approved_agi(vcf)

        with self._tiny_thresholds():
            result = call_operation(
                "prs.calculate_score",
                {"pgs_id": "PGS900001", "score_mean": 1.0, "score_sd": 0.5},
            )

        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["score_result"]["calibration"]["status"], "standardized_from_supplied_parameters")
        self.assertAlmostEqual(result["score_result"]["calibration"]["z_score"], 2.0)
        self.assertIn("user-supplied", result["score_result"]["calibration"]["meaning"])

    def test_cross_build_score_without_liftover_chains_prompts_install(self) -> None:
        # GRCh38 score against a GRCh37 sample, but liftover-chains library is
        # not installed in the tmp GENOMI_HOME — the runtime must surface the
        # liftover-chains install prompt rather than silently producing the
        # wrong result.
        scoring_file = self._write_scoring_file()
        imported = call_operation(
            "prs.import_scoring_file",
            {"pgs_id": "PGS900001", "scoring_file": str(scoring_file), "genome_build": "GRCh38"},
        )
        vcf = self._write_indexed_vcf("sample_grch37.vcf")
        runtime_context.set_active_agi_from_source(
            vcf,
            status="parsed",
            agi_path=default_agi_path(vcf),
            genome_build="GRCh37",
        )
        runtime_context.approve_agi_access(reason="test approved Active Genome Index access")

        result = call_operation("prs.check_score_overlap", {"score_dir": imported["score_cache"]["score_dir"]})

        self.assertEqual(result["status"], "requires_library_install")
        self.assertEqual(result["missing_library"]["library"], "liftover-chains")
        self.assertEqual(result["score_genome_build"], "GRCh38")
        self.assertEqual(result["sample_genome_build"], "GRCh37")
        self.assertEqual(result["polygenic_score"]["pgs_id"], "PGS900001")

    def test_cross_build_score_with_chains_but_missing_pyliftover_prompts_install(self) -> None:
        scoring_file = self._write_scoring_file()
        imported = call_operation(
            "prs.import_scoring_file",
            {"pgs_id": "PGS900001", "scoring_file": str(scoring_file), "genome_build": "GRCh38"},
        )
        vcf = self._write_indexed_vcf("sample_grch37_missing_pyliftover.vcf")
        runtime_context.set_active_agi_from_source(
            vcf,
            status="parsed",
            agi_path=default_agi_path(vcf),
            genome_build="GRCh37",
        )
        runtime_context.approve_agi_access(reason="test approved Active Genome Index access")
        self._write_fake_liftover_chains()

        with mock.patch("genomi.runtime.liftover.importlib.import_module", side_effect=ImportError("missing pyliftover")):
            result = call_operation("prs.check_score_overlap", {"score_dir": imported["score_cache"]["score_dir"]})

        self.assertEqual(result["status"], "requires_library_install")
        self.assertEqual(result["reason"], "missing_python_dependency")
        self.assertEqual(result["missing_library"]["library"], "pyliftover")
        self.assertEqual(result["score_genome_build"], "GRCh38")
        self.assertEqual(result["sample_genome_build"], "GRCh37")
        self.assertEqual(result["polygenic_score"]["pgs_id"], "PGS900001")
        self.assertEqual(result["evidence_envelope"]["finding_state"], "blocked_missing_library")
        self.assertEqual(result["evidence_envelope"]["coverage"]["libraries"][0]["library"], "pyliftover")

    def test_cross_build_score_with_liftover_chains_lifts_variants(self) -> None:
        # Same scenario but with the real UCSC chain files linked into the
        # test GENOMI_HOME. The scoring file declares GRCh38 coordinates for
        # APOE rs429358 and rs7412; the AGI on GRCh37 carries those SNPs at
        # their GRCh37 coordinates. The runtime must lift the score variants
        # onto GRCh37 and match them in the Active Genome Index, completing the calculation.
        if not self._link_real_liftover_chains():
            self.skipTest("liftover setup not available on this host")
        from genomi.capabilities.prs import harmonize as prs_harmonize

        prs_harmonize.get_liftover.cache_clear()

        scoring_file = Path(self._home_tmp.name) / "PGS900099_hmPOS_GRCh38.txt"
        # rs429358: GRCh38 chr19:44908684 -> GRCh37 chr19:45411941
        # rs7412:   GRCh38 chr19:44908822 -> GRCh37 chr19:45412079
        scoring_file.write_text(
            "\n".join(
                [
                    "#pgs_id=PGS900099",
                    "hm_chr\thm_pos\trsID\teffect_allele\tother_allele\teffect_weight",
                    "19\t44908684\trs429358\tC\tT\t0.5",
                    "19\t44908822\trs7412\tT\tC\t1.0",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        imported = call_operation(
            "prs.import_scoring_file",
            {"pgs_id": "PGS900099", "scoring_file": str(scoring_file), "genome_build": "GRCh38"},
        )

        # Build a tiny GRCh37 AGI that carries both APOE SNPs at their
        # GRCh37 coordinates so the lifted score variants find matches.
        vcf = Path(self._home_tmp.name) / "sample_apoe_grch37.vcf"
        vcf.write_text(
            "\n".join(
                [
                    "##fileformat=VCFv4.2",
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                    "19\t45411941\trs429358\tT\tC\t.\tPASS\t.\tGT\t0/1",
                    "19\t45412079\trs7412\tC\tT\t.\tPASS\t.\tGT\t0/0",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        create_active_genome_index(vcf, parallel_workers=1, reuse_existing=False)
        runtime_context.set_active_agi_from_source(
            vcf,
            status="parsed",
            agi_path=default_agi_path(vcf),
            genome_build="GRCh37",
        )
        runtime_context.approve_agi_access(reason="test approved Active Genome Index access")

        with self._tiny_thresholds(min_variants=1, min_fraction=0.10):
            result = call_operation(
                "prs.calculate_score",
                {"score_dir": imported["score_cache"]["score_dir"]},
            )

        self.assertEqual(result["status"], "completed", result)
        self.assertEqual(result["sample_qc"]["genome_build"], "GRCh37")
        self.assertEqual(result["sample_qc"]["score_genome_build"], "GRCh38")
        liftover = result["sample_qc"]["liftover"]
        self.assertEqual(liftover["source_build"], "GRCh38")
        self.assertEqual(liftover["target_build"], "GRCh37")
        self.assertEqual(liftover["lifted_variant_count"], 2)
        self.assertEqual(liftover["dropped_variant_count"], 0)
        self.assertEqual(result["sample_qc"]["matched_variant_count"], 2)

    def test_cross_build_liftover_drops_are_excluded_in_variant_accounting(self) -> None:
        scoring_file = self._write_scoring_file()
        imported = call_operation(
            "prs.import_scoring_file",
            {"pgs_id": "PGS900001", "scoring_file": str(scoring_file), "genome_build": "GRCh38"},
        )
        vcf = self._write_indexed_vcf("sample_grch37_liftover_drops.vcf")
        runtime_context.set_active_agi_from_source(
            vcf,
            status="parsed",
            agi_path=default_agi_path(vcf),
            genome_build="GRCh37",
        )
        runtime_context.approve_agi_access(reason="test approved Active Genome Index access")

        class FakeLifter:
            def lift_position_full(self, chrom: str, pos: int) -> tuple[str, int, str] | None:
                if pos == 200:
                    return None
                if pos == 300:
                    return str(chrom), pos, "-"
                return str(chrom), pos, "+"

        with (
            mock.patch.object(prs_scorer, "liftover_preflight", return_value={"status": "available"}),
            mock.patch.object(prs_harmonize, "get_liftover", return_value=FakeLifter()),
            self._tiny_thresholds(min_variants=1, min_fraction=0.10),
        ):
            result = call_operation(
                "prs.calculate_score",
                {"score_dir": imported["score_cache"]["score_dir"]},
            )

        self.assertEqual(result["status"], "completed", result)
        sample_qc = result["sample_qc"]
        self.assertEqual(sample_qc["score_variant_count"], 4)
        self.assertEqual(sample_qc["matched_variant_count"], 2)
        self.assertEqual(sample_qc["missing_variant_count"], 0)
        self.assertEqual(sample_qc["excluded_variant_count"], 2)
        self.assertEqual(sample_qc["accounted_variant_count"], 4)
        self.assertEqual(sample_qc["unaccounted_variant_count"], 0)
        self.assertEqual(sample_qc["overaccounted_variant_count"], 0)
        self.assertTrue(sample_qc["accounting_complete"])
        self.assertEqual(sample_qc["excluded_reasons"]["liftover_unmapped"], 1)
        self.assertEqual(sample_qc["excluded_reasons"]["liftover_strand_flipped"], 1)
        self.assertEqual(sample_qc["liftover"]["dropped_variant_count"], 2)
        self.assertEqual(sample_qc["liftover"]["dropped_reasons"], {"unmapped": 1, "strand_flipped": 1})
        accounting = result["variant_accounting"]
        self.assertEqual(accounting["accounted_variant_count"], 4)
        self.assertEqual(accounting["excluded_count"], 2)
        excluded_reasons = {item["reason"] for item in accounting["excluded_examples"]}
        self.assertEqual(excluded_reasons, {"liftover_unmapped", "liftover_strand_flipped"})

    def test_low_overlap_blocks_default_score_calculation(self) -> None:
        scoring_file = self._write_scoring_file()
        call_operation("prs.import_scoring_file", {"pgs_id": "PGS900001", "scoring_file": str(scoring_file)})
        vcf = self._write_indexed_vcf("sample_partial.vcf", include_positions={100})
        self._select_approved_agi(vcf)

        with self._tiny_thresholds(min_variants=2, min_fraction=0.75):
            result = call_operation("prs.calculate_score", {"pgs_id": "PGS900001"})

        self.assertEqual(result["status"], "insufficient_overlap")
        self.assertIsNone(result["score_result"])
        self.assertFalse(result["sample_qc"]["calculation_allowed"])
        self.assertLessEqual(
            set(result["sample_qc"]),
            {
                "genome_build",
                "score_genome_build",
                "score_variant_count",
                "matched_variant_count",
                "missing_variant_count",
                "excluded_variant_count",
                "accounted_variant_count",
                "unaccounted_variant_count",
                "overaccounted_variant_count",
                "accounting_complete",
                "overlap_fraction",
                "overlap_status",
                "calculation_allowed",
                "overlap_quality",
                "missing_reasons",
                "excluded_reasons",
                "note",
                "liftover",
            },
        )
        self.assertEqual(result["evidence_envelope"]["coverage"]["consulted_sources"], ["local_active_genome_index", "local_prs_score_cache"])
        self.assertEqual(result["evidence_envelope"]["observations"]["matched_variant_count"], 1)


if __name__ == "__main__":
    unittest.main()
