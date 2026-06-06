from __future__ import annotations

from pathlib import Path
from unittest import mock

from genomi.capabilities.prs import pgs_catalog as prs_pgs_catalog
from genomi.operations import OperationError, call_operation, list_operations
from genomi.runtime import context as runtime_context

from tests.support.capabilities.prs_contract import PolygenicScoreTestBase


class PolygenicScoreDiscoveryTests(PolygenicScoreTestBase):
    def test_public_tools_do_not_require_personal_approval(self) -> None:
        source_context = call_operation("prs.build_source_context")
        imported = call_operation("prs.list_imported_scores")

        self.assertEqual(source_context["status"], "completed")
        self.assertIn("raw weighted score", " ".join(source_context["method_boundaries"]["does"]).lower())
        self.assertEqual(imported["status"], "completed")
        self.assertEqual(imported["score_count"], 0)

    def test_search_scores_uses_host_semantic_terms_without_hardcoded_synonyms(self) -> None:
        rows = [
            self._pgs_metadata_row(
                pgs_id="PGS001987",
                name="portability-PLR_M_less_hair",
                reported_trait="Hair/balding pattern",
                mapped_trait_labels="balding measurement",
                mapped_trait_ids="EFO_0007825",
                variant_count="23692",
            ),
            self._pgs_metadata_row(
                pgs_id="PGS900010",
                name="lipids-ldl",
                reported_trait="LDL cholesterol",
                mapped_trait_labels="low density lipoprotein cholesterol measurement",
                mapped_trait_ids="EFO_0004611",
            ),
            self._pgs_metadata_row(
                pgs_id="PGS900011",
                name="hair-color",
                reported_trait="Hair color",
                mapped_trait_labels="hair color measurement",
                mapped_trait_ids="EFO_0007824",
            ),
        ]

        with mock.patch.object(prs_pgs_catalog, "_fetch_score_metadata_rows", return_value=rows):
            result = call_operation(
                "prs.search_scores",
                {
                    "query": "will I go bald",
                    "limit": 3,
                    "semantic_context": {
                        "raw_query": "will I go bald",
                        "host_expansions": ["male pattern baldness", "androgenetic alopecia", "hair loss"],
                        "host_entities": [
                            {"text": "androgenetic alopecia", "type": "trait_or_condition"}
                        ],
                    },
                },
            )

        self.assertEqual(result["status"], "completed")
        self.assertIn(result["retrieval"]["model"], {"hybrid_bm25_rrf_v1", "persistent_sqlite_fts5_bm25_rrf_v1"})
        self.assertFalse(result["retrieval"]["semantic_query_model"]["hardcoded_synonyms"])
        self.assertEqual(result["results"][0]["pgs_id"], "PGS001987")
        self.assertEqual(result["results"][0]["mapped_trait_ids"], "EFO_0007825")
        semantic_context = result["semantic_context"]
        self.assertEqual(
            {
                "raw_query",
                "host_expansions",
                "host_entities",
                "term_matches",
                "term_misses",
                "ignored_hints",
                "retrieval_streams",
                "retrieval_boundary",
            },
            set(semantic_context),
        )
        matches = {item["text"] for item in semantic_context["term_matches"]}
        misses = {item["text"] for item in semantic_context["term_misses"]}
        self.assertIn("male pattern baldness", matches)
        self.assertIn("androgenetic alopecia", misses)
        self.assertTrue(all(item["status"] == "hit" for item in semantic_context["term_matches"]))
        self.assertTrue(all(item["status"] in {"miss", "ignored_for_exact_identifier"} for item in semantic_context["term_misses"]))

    def test_search_scores_filters_by_efo_trait_id(self) -> None:
        rows = [
            self._pgs_metadata_row(
                pgs_id="PGS001987",
                reported_trait="Hair/balding pattern",
                mapped_trait_labels="balding measurement",
                mapped_trait_ids="EFO_0007825",
            ),
            self._pgs_metadata_row(
                pgs_id="PGS900010",
                reported_trait="LDL cholesterol",
                mapped_trait_labels="low density lipoprotein cholesterol measurement",
                mapped_trait_ids="EFO_0004611",
            ),
        ]

        with mock.patch.object(prs_pgs_catalog, "_fetch_score_metadata_rows", return_value=rows):
            result = call_operation("prs.search_scores", {"efo_id": "EFO:0007825", "limit": 5})

        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["summary"]["matched_count"], 1)
        self.assertEqual([item["pgs_id"] for item in result["results"]], ["PGS001987"])

    def test_search_scores_refreshes_public_retrieval_index(self) -> None:
        rows = [
            self._pgs_metadata_row(
                pgs_id="PGS001987",
                reported_trait="Hair/balding pattern",
                mapped_trait_labels="balding measurement",
                mapped_trait_ids="EFO_0007825",
            )
        ]

        with mock.patch.object(prs_pgs_catalog, "_fetch_score_metadata_rows", return_value=rows):
            call_operation("prs.search_scores", {"query": "balding", "limit": 1})

        listed = call_operation("genomi.search_indexes", {"source": "pgs_scores", "query": "balding"})
        self.assertEqual(listed["search_results"][0]["source"], "pgs_scores")
        self.assertEqual(listed["search_results"][0]["hits"][0]["doc_id"], "PGS001987")

    def test_search_scores_requires_installed_metadata_library(self) -> None:
        result = call_operation("prs.search_scores", {"query": "balding", "limit": 1})

        self.assertEqual(result["status"], "requires_library_install")
        self.assertFalse(result["tool_will_work"])
        self.assertEqual(result["missing_library"]["library"], "pgs-catalog-score-metadata")
        self.assertEqual(result["operation"], "prs.search_scores")
        self.assertIn("install_command", result["ask_user"])

    def test_private_tools_require_approval_for_existing_active_context(self) -> None:
        vcf = Path(self._home_tmp.name) / "sample.vcf"
        runtime_context.set_active_agi_from_source(
            vcf,
            status="parsed",
            agi_path=vcf.with_suffix(".sqlite"),
            genome_build="GRCh38",
        )

        with self.assertRaises(OperationError) as raised:
            call_operation("prs.calculate_score", {"pgs_id": "PGS900001"})
        self.assertEqual(raised.exception.code, "active_genome_index_approval_required")

    def test_active_agi_returns_requires_score_import_with_defaults(self) -> None:
        vcf = self._write_indexed_vcf("sample_requires_import.vcf")
        self._select_approved_agi(vcf)

        result = call_operation("prs.calculate_score", {"pgs_id": "PGS900001"})

        self.assertEqual(result["status"], "requires_score_import")
        self.assertTrue(result["personal_context"]["uses_personal_dna"])
        self.assertEqual(result["missing_library"]["library"], "PGS900001")
        self.assertEqual(result["missing_library"]["status"], "not_installed")
        self.assertIn("genomi call prs.import_scoring_file", result["ask_user"]["install_command"])
        self.assertIn("PGS900001", result["ask_user"]["question"])
        envelope = result["evidence_envelope"]
        self.assertEqual(envelope["finding_state"], "blocked_missing_library")
        self.assertEqual(envelope["answer_readiness"], "needs_user_install")
        self.assertEqual(envelope["coverage"]["libraries"][0]["library"], "PGS900001")
        self.assertIn("genomi call prs.import_scoring_file", envelope["coverage"]["libraries"][0]["install_command"])
        defaults = {item["parameter"]: item for item in result["defaults_applied"]}
        self.assertEqual(defaults["genome_build"]["value"], "GRCh38")
        self.assertTrue(defaults["skip_ambiguous_palindromic"]["value"])

    def test_discovery_registers_all_prs_handlers(self) -> None:
        tools = {tool["name"]: tool for tool in list_operations(capability="polygenic-score")}

        self.assertEqual(
            set(tools),
            {
                "prs.search_scores",
                "prs.fetch_score_metadata",
                "prs.import_scoring_file",
                "prs.list_imported_scores",
                "prs.check_score_overlap",
                "prs.calculate_score",
                "prs.build_source_context",
            },
        )
        self.assertEqual(tools["prs.calculate_score"]["annotations"]["discoveryRole"], "entry_tool")
        self.assertEqual(tools["prs.calculate_score"]["annotations"]["privacyScope"], "local_private_prs_score")
        self.assertEqual(tools["prs.calculate_score"]["annotations"]["agiNeed"], "reference")
        self.assertIn("pgs_catalog_ftp", tools["prs.import_scoring_file"]["annotations"]["externalIO"])
