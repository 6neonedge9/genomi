from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from genomi.capabilities.clinvar.static_annotation import (
    default_static_outputs,
)
from genomi.capabilities.clinvar.static_annotation import workflow_contract as static_contract
from genomi.capabilities.research.intent_research import query_reviewed_research, record_reviewed_research
from genomi.capabilities.research.intent_research import workflow_contract as research_contract
from genomi.evidence import (
    fetch_gene_evidence,
    gather_variant_evidence,
    import_clinvar_vcf,
    match_clinvar_variants,
    query_research_findings,
    record_research_findings,
    search_research_findings,
)
from genomi.evidence.sources import evidence_source_catalog
from tests.support.capabilities.external_layers import (
    TINY_CLINVAR,
    TINY_VCF,
    EvidenceImportTestBase,
)


class ExternalResearchContextTests(EvidenceImportTestBase):
    def test_research_findings_are_recorded_and_returned_with_gather_actions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"
            matches = Path(tmp) / "matches.jsonl"
            import_clinvar_vcf(TINY_CLINVAR, db, source_version="fixture")
            match_clinvar_variants(TINY_VCF, db, matches)

            recorded = record_research_findings(
                db,
                {
                    "findings": [
                        {
                            "target": {"type": "variant", "chrom": "1", "pos": 10250, "ref": "A", "alt": "C"},
                            "source": {
                                "title": "Example Variant Source",
                                "url": "https://example.test/variant",
                                "type": "source",
                                "accessed_at": "2026-05-05T00:00:00+00:00",
                            },
                            "finding": {
                                "text": "Original short finding text from the source.",
                                "summary": "Variant source says the fact.",
                                "type": "variant_assertion",
                            },
                            "searched_query": "1 10250 A C",
                        },
                        {
                            "target": {"type": "gene", "gene": "GENE1"},
                            "source": {
                                "title": "Example Gene Source",
                                "url": "https://example.test/gene",
                                "accessed_at": "2026-05-05T00:00:00+00:00",
                            },
                            "finding": {
                                "text": "Original gene finding text from the source.",
                                "summary": "Gene source says the fact.",
                            },
                        },
                    ]
                },
            )

            self.assertEqual(recorded["inserted_findings"], 2)
            variant_evidence = gather_variant_evidence(db, "1", 10250, "A", "C", matches_path=matches)
            self.assertEqual(variant_evidence["research_evidence"]["exact_variant"]["count"], 1)
            self.assertEqual(
                variant_evidence["research_evidence"]["exact_variant"]["records"][0]["source"]["url"],
                "https://example.test/variant",
            )
            self.assertEqual(variant_evidence["research_evidence"]["genes"]["GENE1"]["count"], 1)

            gene_evidence = fetch_gene_evidence("GENE1", db, matches_path=matches)
            self.assertEqual(gene_evidence["research_evidence"]["count"], 1)
            self.assertEqual(
                gene_evidence["research_evidence"]["records"][0]["finding"]["text"],
                "Original gene finding text from the source.",
            )

    def test_research_findings_support_drug_condition_and_topic_targets(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"

            recorded = record_research_findings(
                db,
                {
                    "findings": [
                        {
                            "target": {"type": "drug", "drug": "Warfarin"},
                            "source": {
                                "title": "Example CPIC Warfarin",
                                "url": "https://example.test/cpic-warfarin",
                                "type": "pgx_guideline",
                                "accessed_at": "2026-05-07T00:00:00+00:00",
                            },
                            "finding": {
                                "text": "Short warfarin pharmacogenomic finding.",
                                "summary": "Warfarin source context.",
                                "type": "pharmacogenomic_guideline",
                            },
                        },
                        {
                            "target": {"type": "condition", "condition": "Alpha-1 antitrypsin deficiency"},
                            "source": {
                                "title": "Example Condition Source",
                                "url": "https://example.test/aatd",
                                "accessed_at": "2026-05-07T00:00:00+00:00",
                            },
                            "finding": {"text": "Short condition finding."},
                        },
                        {
                            "target": {"type": "topic", "topic": "smoking related genetic risk"},
                            "source": {
                                "title": "Example Topic Source",
                                "url": "https://example.test/smoking-topic",
                                "accessed_at": "2026-05-07T00:00:00+00:00",
                            },
                            "finding": {"text": "Short topic finding."},
                        },
                    ]
                },
            )

            self.assertEqual(recorded["inserted_findings"], 3)
            self.assertEqual(
                {step["target_type"] for step in recorded["evidence_options"]},
                {"drug", "condition", "topic"},
            )

            drug = query_research_findings(db, "drug", drug="warfarin")
            self.assertEqual(drug["count"], 1)
            self.assertEqual(drug["query"]["drug"], "warfarin")
            self.assertEqual(drug["records"][0]["target"]["drug"], "Warfarin")
            self.assertEqual(drug["records"][0]["finding"]["type"], "pharmacogenomic_guideline")

            condition = query_research_findings(db, "condition", condition="Alpha-1 antitrypsin deficiency")
            self.assertEqual(condition["count"], 1)
            self.assertEqual(condition["records"][0]["target"]["condition"], "Alpha-1 antitrypsin deficiency")

            topic = query_research_findings(db, "topic", topic="smoking related genetic risk")
            self.assertEqual(topic["count"], 1)
            self.assertEqual(topic["records"][0]["target"]["topic"], "smoking related genetic risk")

            search = search_research_findings(db, "smoking risk")
            self.assertEqual(search["count"], 1)
            self.assertEqual(search["records"][0]["target"]["type"], "topic")

            semantic_search = search_research_findings(
                db,
                "blood thinner after stent",
                semantic_context={
                    "raw_query": "blood thinner after stent",
                    "host_expansions": ["warfarin"],
                    "host_entities": [{"text": "warfarin", "type": "drug"}],
                },
            )
            self.assertEqual(semantic_search["count"], 1)
            self.assertEqual(semantic_search["records"][0]["target"]["drug"], "Warfarin")
            self.assertIn(
                "warfarin",
                {item["text"] for item in semantic_search["semantic_context"]["term_matches"]},
            )

    def test_source_catalog_exposes_more_databases_and_storage_contract(self) -> None:
        catalog = evidence_source_catalog(target_type="drug")
        source_ids = {source["source_id"] for source in catalog["sources"]}

        self.assertIn("cpic", source_ids)
        self.assertIn("pharmgkb", source_ids)
        self.assertIn("fda_pharmacogenomics", source_ids)
        self.assertIn("fda_pharmacogenetic_associations", source_ids)
        self.assertIn("drug", catalog["storage_contract"]["record_research_target_types"])

    def test_intent_research_scope_separates_shared_and_private_findings(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "evidence.sqlite"
            shared_db = Path(tmp) / "shared.sqlite"
            shared = {
                "target": {"type": "gene", "gene": "GENE1"},
                "source": {
                    "title": "Shared Source",
                    "url": "https://example.test/shared-gene",
                    "accessed_at": "2026-05-07T00:00:00+00:00",
                },
                "finding": {"text": "Shared source finding.", "type": "clinical_review"},
            }
            private = {
                "target": {"type": "gene", "gene": "GENE1"},
                "source": {
                    "title": "Private Source",
                    "url": "https://example.test/private-gene",
                    "accessed_at": "2026-05-07T00:00:00+00:00",
                },
                "finding": {"text": "Private user-specific combination finding."},
            }

            shared_record = record_reviewed_research(db, shared, scope="shared", shared_evidence_db=shared_db)
            private_record = record_reviewed_research(db, private, scope="private", shared_evidence_db=shared_db)

            all_rows = query_reviewed_research(db, "gene", gene="GENE1")
            shared_rows = query_reviewed_research(db, "gene", gene="GENE1", scope="shared")
            private_rows = query_reviewed_research(db, "gene", gene="GENE1", scope="private")

            self.assertEqual(all_rows["count"], 2)
            self.assertEqual(shared_rows["count"], 1)
            self.assertEqual(private_rows["count"], 1)
            self.assertEqual(shared_rows["records"][0]["scope"], "shared")
            self.assertEqual(private_rows["records"][0]["scope"], "private")
            self.assertEqual(shared_record["shared_sync"]["status"], "completed")
            self.assertEqual(private_record["shared_sync"]["status"], "private_not_synced")
            shared_db_rows = query_reviewed_research(shared_db, "gene", gene="GENE1")
            self.assertEqual(shared_db_rows["count"], 1)
            self.assertEqual(shared_db_rows["records"][0]["scope"], "shared")

    def test_workflow_contracts_are_explicit(self) -> None:
        self.assertEqual(static_contract()["id"], "static")
        self.assertEqual(research_contract()["id"], "research")
        self.assertIn("local parsing, database import, and deterministic evidence checks", static_contract()["purpose"])
        self.assertEqual(
            static_contract()["primary_outputs"],
            [
                "run project layout",
                "Active Genome Index",
                "source-format metadata",
                "sequencing-derived sample QC and genotype/callability support rows",
                "consumer-array rsID/locus observations when supplied",
                "library-scoped ClinVar exact-match JSONL when ClinVar matching is requested",
                "target-scoped candidate inventory when requested",
                "canonical shared evidence DB for reusable static rows",
                "per-run user evidence DB for sample-specific context",
                "SQLite evidence rows for lazily materialized public sources",
            ],
        )
        self.assertFalse({"panel_json", "panel_markdown"} & set(default_static_outputs(TINY_VCF)))
        self.assertIn("shared/private", research_contract()["purpose"])


if __name__ == "__main__":
    unittest.main()
