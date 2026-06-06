from __future__ import annotations

from unittest import mock

from genomi.evidence import init_evidence_db
from genomi.operations import call_operation

from tests.support.active_genome_index.contract_cases import LocusContract, SourceContractCase
from tests.support.matrix.capability_contract import MatrixCaseContext
from tests.support.matrix.source_support_operations import SupportsAssertions


def assert_external_source_operations(
    testcase: SupportsAssertions,
    contract: SourceContractCase,
    ctx: MatrixCaseContext,
    *,
    allele_locus: LocusContract,
) -> set[str]:
    seen_operations: set[str] = set()
    _assert_gnomad_fetch(testcase, ctx, allele_locus, seen_operations)
    _assert_pgx_fetchers(testcase, seen_operations)
    _assert_gwas_fetcher(testcase, seen_operations)
    _assert_functional_genomics_fetchers(testcase, ctx, seen_operations)
    _assert_phenotype_fetchers(testcase, contract, seen_operations)
    return seen_operations


def _assert_gnomad_fetch(
    testcase: SupportsAssertions,
    ctx: MatrixCaseContext,
    allele_locus: LocusContract,
    seen_operations: set[str],
) -> None:
    run_db = ctx.tmp_path / "gnomad.run.sqlite"
    shared_db = ctx.tmp_path / "gnomad.shared.sqlite"
    init_evidence_db(run_db)
    variant = {
        "variant_id": f"{allele_locus.chrom}-{allele_locus.pos}-{allele_locus.ref}-{allele_locus.alt}",
        "rsids": [allele_locus.rsid],
        "chrom": allele_locus.chrom,
        "pos": allele_locus.pos,
        "ref": allele_locus.ref,
        "alt": allele_locus.alt,
        "exome": {"ac": 2, "an": 20, "af": 0.1, "homozygote_count": 0, "populations": []},
        "genome": None,
    }
    with mock.patch("genomi.evidence._post_graphql", return_value={"data": {"variant": variant}}):
        result = call_operation(
            "gnomad.fetch_population_frequency",
            {
                "db": str(run_db),
                "shared_db": str(shared_db),
                "chrom": allele_locus.chrom,
                "pos": allele_locus.pos,
                "ref": allele_locus.ref,
                "alt": allele_locus.alt,
            },
        )
    seen_operations.add("gnomad.fetch_population_frequency")
    testcase.assertEqual(result["status"], "completed", result)
    testcase.assertEqual(result["population_frequency"]["count"], 1)


def _assert_pgx_fetchers(testcase: SupportsAssertions, seen_operations: set[str]) -> None:
    with mock.patch(
        "genomi.operations.registry.handlers_pgx.clinpgx.lookup_clinpgx",
        return_value={"status": "completed", "summary": {"guideline_annotation_count": 1}, "guideline_annotations": [{}]},
    ):
        clinpgx = call_operation("pharmacogenomics.fetch_clinpgx", {"drug": "clopidogrel", "gene": "CYP2C19"})
    seen_operations.add("pharmacogenomics.fetch_clinpgx")
    testcase.assertEqual(clinpgx["status"], "completed")

    with mock.patch(
        "genomi.operations.registry.handlers_pgx.fda_pgx.lookup_fda_pgx",
        return_value={"status": "completed", "summary": {"association_count": 1}, "associations": [{}]},
    ):
        fda = call_operation("pharmacogenomics.fetch_fda_labels", {"drug": "clopidogrel", "gene": "CYP2C19"})
    seen_operations.add("pharmacogenomics.fetch_fda_labels")
    testcase.assertEqual(fda["status"], "completed")

    with mock.patch(
        "genomi.operations.registry.handlers_pgx.pgxdb.lookup_pgxdb",
        return_value={"status": "completed", "summary": {"pgx_record_count": 1}, "pgx_records": [{}]},
    ):
        pgxdb = call_operation("pharmacogenomics.fetch_pgxdb", {"drug": "infliximab", "rsid": "rs1061622"})
    seen_operations.add("pharmacogenomics.fetch_pgxdb")
    testcase.assertEqual(pgxdb["status"], "completed")


def _assert_gwas_fetcher(testcase: SupportsAssertions, seen_operations: set[str]) -> None:
    with mock.patch(
        "genomi.capabilities.research.intent_research.compare_gwas_variant_evidence",
        return_value={
            "status": "completed",
            "top_observed_candidate": "rs900000002",
            "association_records": [{"variant": "rs900000002"}],
        },
    ):
        result = call_operation("gwas.compare_variant_associations", {"phenotype": "LDL cholesterol", "variants": ["rs900000002"]})
    seen_operations.add("gwas.compare_variant_associations")
    testcase.assertEqual(result["status"], "completed")
    testcase.assertEqual(result["top_observed_candidate"], "rs900000002")


def _assert_functional_genomics_fetchers(
    testcase: SupportsAssertions,
    ctx: MatrixCaseContext,
    seen_operations: set[str],
) -> None:
    effect_path = ctx.write_text(
        "depmap-effect.csv",
        "ModelID,EGFR (1956),MYC (4609)\nACH-000001,-1.25,-0.18\n",
    )
    model_path = ctx.write_text("depmap-model.csv", "ModelID,CellLineName,CCLEName\nACH-000001,A549,A549_LUNG\n")
    records = call_operation(
        "functional_genomics.retrieve_perturbation_records",
        {
            "context": "A549 CRISPR dependency",
            "genes": ["EGFR", "MYC"],
            "cell_line": "A549",
            "perturbation": "CRISPR knockout",
            "phenotype": "dependency",
            "perturbation_sources": ["depmap"],
            "depmap_gene_effect_url": effect_path,
            "depmap_model_url": model_path,
        },
    )
    seen_operations.add("functional_genomics.retrieve_perturbation_records")
    testcase.assertEqual(records["coverage_state"], "data_returned", records)

    with mock.patch(
        "genomi.operations.registry.geo.query_geo_datasets",
        return_value={
            "status": "geo_metadata_found",
            "coverage_state": "metadata_only",
            "summary": {"record_count": 0, "geo_hit_count": 1},
            "geo_hits": [{"accession": "GSE12345"}],
            "download_candidates": [],
            "source_records": [],
            "direct_perturbation_source_records": [],
            "records_by_gene": {},
        },
    ):
        geo = call_operation("functional_genomics.query_geo", {"context": "GSE12345 A549 CRISPR dependency", "accession": "GSE12345"})
    seen_operations.add("functional_genomics.query_geo")
    testcase.assertEqual(geo["status"], "geo_metadata_found")


def _assert_phenotype_fetchers(
    testcase: SupportsAssertions,
    contract: SourceContractCase,
    seen_operations: set[str],
) -> None:
    with mock.patch(
        "genomi.operations.registry.handlers_evidence_phenotype.gene_identification.retrieve_trait_gene_records",
        return_value={
            "status": "trait_gene_records_found",
            "coverage_state": "data_returned",
            "gene_records": [{"gene": "PCSK9", "direct_record_count": 1}],
            "source_records": [{"source_id": "opentargets"}],
        },
    ):
        trait = call_operation("phenotype.retrieve_trait_gene_records", {"trait": "LDL cholesterol", "genes": ["PCSK9"]})
    seen_operations.add("phenotype.retrieve_trait_gene_records")
    testcase.assertEqual(trait["status"], "trait_gene_records_found")

    with mock.patch(
        "genomi.operations.registry.handlers_evidence_phenotype.targets.retrieve_disease_clinical_drug_targets",
        return_value={
            "status": "completed",
            "coverage_state": "data_returned",
            "targets": [{"gene": "ADRB2", "disease": "asthma"}],
            "source_records": [{"source_id": "opentargets"}],
        },
    ):
        disease = call_operation("phenotype.retrieve_disease_drug_targets", {"disease": "asthma", "genes": ["ADRB2"]})
    seen_operations.add("phenotype.retrieve_disease_drug_targets")
    testcase.assertEqual(disease["status"], "completed")
    testcase.assertTrue(contract.expected_format)
