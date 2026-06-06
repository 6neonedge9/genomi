from __future__ import annotations

import tempfile
from pathlib import Path

from genomi.capabilities.clinvar.static_annotation import (
    run_static_callability,
    run_static_genotype_support,
    run_static_sample_qc,
)
from genomi.evidence import (
    evidence_summary,
    gather_variant_evidence,
    import_clinvar_vcf,
    match_clinvar_variants,
)
from tests.support.capabilities.external_layers import (
    TINY_CLINVAR,
    TINY_VCF,
    EvidenceImportTestBase,
)


class ExternalGenotypeContextTests(EvidenceImportTestBase):
    def test_sample_qc_records_evidence_boundaries_in_private_db(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"
            index = Path(tmp) / "active-genome-index.sqlite"
            output = Path(tmp) / "sample-qc.json"

            result = run_static_sample_qc(TINY_VCF, evidence_db=db, agi_path=index, output=output)
            summary = evidence_summary(db)

            self.assertEqual(result["status"], "completed")
            self.assertEqual(result["agi_intake_source_path"], str(TINY_VCF))
            self.assertEqual(result["input_type"], "callset_with_reference_blocks")
            self.assertTrue(result["has_reference_blocks"])
            self.assertTrue(result["absence_claims_allowed_by_default"])
            self.assertIn("evidence_boundaries", result["evidence_boundaries"])
            self.assertEqual(summary["tables"]["sample_qc"], 1)
            self.assertTrue(output.exists())

    def test_genotype_support_classifies_supported_and_weak_calls(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"
            index = Path(tmp) / "active-genome-index.sqlite"

            supported = run_static_genotype_support(
                TINY_VCF,
                "1",
                10250,
                "A",
                "C",
                evidence_db=db,
                agi_path=index,
            )
            self.assertEqual(supported["support_status"], "supported")
            self.assertEqual(supported["agi_intake_source_path"], str(TINY_VCF))
            self.assertIn("genotype_support_supported", supported["accepted_report_evidence_classes"])

            weak_vcf = Path(tmp) / "weak.vcf"
            weak_vcf.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "##reference=GRCh38",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                        "1\t200\t.\tA\tG\t.\tLowQual\t.\tGT:DP:GQ\t0/1:4:9",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            weak = run_static_genotype_support(
                weak_vcf,
                "1",
                200,
                "A",
                "G",
                evidence_db=db,
                agi_path=Path(tmp) / "weak-active-genome-index.sqlite",
            )

            self.assertEqual(weak["support_status"], "weak")
            self.assertEqual(weak["evidence_class"], "genotype_support_weak")
            self.assertIn("limitation context", weak["evidence_boundaries"]["evidence_boundaries"][0])
            self.assertEqual(evidence_summary(db)["tables"]["genotype_support"], 2)

    def test_genotype_support_resolves_interior_gvcf_reference_block_with_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"
            index = Path(tmp) / "active-genome-index.sqlite"
            fasta = Path(tmp) / "ref.fa"
            gvcf = Path(tmp) / "sample.g.vcf"
            fasta.write_text(">1\nAAAAACCCCCGGGGGTTTTT\n", encoding="utf-8")
            Path(f"{fasta}.fai").write_text("1\t20\t3\t20\t21\n", encoding="utf-8")
            gvcf.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "##reference=GRCh38",
                        "##contig=<ID=1,length=20>",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                        "1\t5\t.\tA\t.\t.\tPASS\tEND=15\tGT:DP:GQ\t0/0:35:0",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            result = run_static_genotype_support(
                gvcf,
                "1",
                10,
                "C",
                "T",
                evidence_db=db,
                agi_path=index,
                reference_fasta=fasta,
            )

            self.assertEqual(result["support_status"], "not_observed")
            self.assertEqual(result["sample_observation"]["observed_genotype"], "C/C")
            self.assertTrue(result["sample_observation"]["reference_call_supported"])
            self.assertEqual(result["sample_observation"]["matched_by"], "reference_block")
            self.assertIn("reference_inference_or_assay_completeness", result["accepted_report_evidence_classes"])

    def test_genotype_support_preserves_sample_specific_records(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"
            index = Path(tmp) / "active-genome-index.sqlite"
            vcf = Path(tmp) / "multi.vcf"
            vcf.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "##reference=GRCh38",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\tSample1",
                        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t0/0:20:60\t0/1:50:80",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            result = run_static_genotype_support(
                vcf,
                "1",
                100,
                "A",
                "G",
                evidence_db=db,
                agi_path=index,
            )

            self.assertEqual(result["support_status"], "supported")
            self.assertEqual(result["sample_observation"]["observed_genotype"], "A/G")
            self.assertEqual(result["sample_observation"]["source_record"]["sample_name"], "Sample1")

    def test_genotype_support_projects_simple_complex_record_to_target_base(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"
            index = Path(tmp) / "active-genome-index.sqlite"
            vcf = Path(tmp) / "sample.vcf"
            vcf.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "##reference=GRCh38",
                        "##contig=<ID=1,length=200>",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                        "1\t100\t.\tAC\tAT\t.\tPASS\t.\tGT:DP:GQ\t0/1:35:80",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            result = run_static_genotype_support(
                vcf,
                "1",
                101,
                "C",
                "T",
                evidence_db=db,
                agi_path=index,
            )

            self.assertEqual(result["support_status"], "supported")
            self.assertEqual(result["sample_observation"]["observed_genotype"], "C/T")
            self.assertEqual(result["sample_observation"]["matched_by"], "overlapping_variant_projection")
            self.assertEqual(result["sample_observation"]["record_type"], "complex_projection")

    def test_genotype_support_projects_deletion_allele_to_dash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"
            index = Path(tmp) / "active-genome-index.sqlite"
            vcf = Path(tmp) / "sample.vcf"
            vcf.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "##reference=GRCh38",
                        "##contig=<ID=1,length=200>",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                        "1\t100\t.\tAC\tA\t.\tPASS\t.\tGT:DP:GQ\t0/1:35:80",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            result = run_static_genotype_support(
                vcf,
                "1",
                101,
                "C",
                "-",
                evidence_db=db,
                agi_path=index,
            )

            self.assertEqual(result["support_status"], "supported")
            self.assertEqual(result["sample_observation"]["observed_genotype"], "C/-")
            self.assertEqual(result["sample_observation"]["matched_by"], "overlapping_variant_projection")

    def test_gather_variant_consumes_private_genotype_support_context(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"
            matches = Path(tmp) / "matches.jsonl"
            index = Path(tmp) / "active-genome-index.sqlite"
            import_clinvar_vcf(TINY_CLINVAR, db, source_version="fixture")
            match_clinvar_variants(TINY_VCF, db, matches)
            run_static_sample_qc(TINY_VCF, evidence_db=db, agi_path=index, output=Path(tmp) / "sample-qc.json")
            run_static_genotype_support(
                TINY_VCF,
                "1",
                10257,
                "A",
                "C",
                evidence_db=db,
                agi_path=index,
                min_depth=100,
            )

            result = gather_variant_evidence(db, "1", 10257, "A", "C", matches_path=matches)

            private_context = result["private_sample_context"]
            self.assertEqual(private_context["sample_qc"]["count"], 1)
            self.assertEqual(private_context["genotype_support"]["latest"]["support_status"], "weak")
            self.assertEqual(result["evidence_options"][0]["available_operation"], "active_genome_index.classify_genotype_support")

    def test_callability_requires_reference_blocks_for_negative_claims(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"
            callable_result = run_static_callability(
                TINY_VCF,
                "1:10001-10249",
                evidence_db=db,
                agi_path=Path(tmp) / "active-genome-index.sqlite",
            )

            self.assertEqual(callable_result["callability_status"], "callable")
            self.assertEqual(callable_result["agi_intake_source_path"], str(TINY_VCF))
            self.assertTrue(callable_result["can_support_negative_or_reference_claim"])
            self.assertIn(
                "reference_inference_or_assay_completeness",
                callable_result["accepted_report_evidence_classes"],
            )

            variant_only_vcf = Path(tmp) / "variant-only.vcf"
            variant_only_vcf.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "##reference=GRCh38",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                        "1\t10250\t.\tA\tC\t.\tPASS\t.\tGT:DP:GQ\t0/1:50:99",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            unknown = run_static_callability(
                variant_only_vcf,
                "1:10251-10251",
                evidence_db=db,
                agi_path=Path(tmp) / "variant-only-active-genome-index.sqlite",
            )

            self.assertEqual(unknown["callability_status"], "unknown_no_reference_blocks")
            self.assertFalse(unknown["can_support_negative_or_reference_claim"])
            self.assertEqual(evidence_summary(db)["tables"]["region_callability"], 2)
