from __future__ import annotations

from genomi.operations import call_operation
from genomi.runtime import background_jobs

from tests.support.active_genome_index.contract_cases import LocusContract, SourceContractCase
from tests.support.matrix.source_support_operations import SupportsAssertions


def assert_stateful_runtime_operations(
    testcase: SupportsAssertions,
    contract: SourceContractCase,
    _source: object,
    _parsed: dict[str, object],
    *,
    allele_locus: LocusContract,
) -> set[str]:
    seen_operations: set[str] = {"genomi.parse_source"}
    _assert_runtime_metadata_operations(testcase, seen_operations)
    _assert_dispatcher_and_background_job(testcase, allele_locus, seen_operations)
    _assert_journal_operations(testcase, contract, seen_operations)
    _assert_research_operations(testcase, contract, seen_operations)
    return seen_operations


def _assert_runtime_metadata_operations(testcase: SupportsAssertions, seen_operations: set[str]) -> None:
    context = call_operation("genomi.describe_context")
    seen_operations.add("genomi.describe_context")
    testcase.assertIn("has_active_genome_index", context)

    libraries = call_operation("genomi.check_libraries", {"libraries": ["clinvar-grch38"]})
    seen_operations.add("genomi.check_libraries")
    testcase.assertEqual(libraries["libraries"][0]["library"], "clinvar-grch38")
    testcase.assertIn("install_command", libraries["libraries"][0])

    resources = call_operation("genomi.list_resources")
    seen_operations.add("genomi.list_resources")
    testcase.assertTrue(resources["resource_groups"])

    indexed = call_operation("genomi.search_indexes", {"query": "CYP2C19", "limit": 1})
    seen_operations.add("genomi.search_indexes")
    testcase.assertEqual(indexed["status"], "completed")

    profile = call_operation("genomi.set_response_profile", {"profile": "literate"})
    seen_operations.add("genomi.set_response_profile")
    testcase.assertEqual(profile["status"], "completed")


def _assert_dispatcher_and_background_job(
    testcase: SupportsAssertions,
    allele_locus: LocusContract,
    seen_operations: set[str],
) -> None:
    invoked = call_operation(
        "genomi.invoke",
        {
            "tool": "variant.resolve",
            "params": {
                "query": f"chr{allele_locus.chrom}:{allele_locus.pos}:{allele_locus.ref}:{allele_locus.alt}",
                "genome_build": "GRCh37",
            },
        },
    )
    seen_operations.add("genomi.invoke")
    testcase.assertEqual(invoked["dispatched_tool"], "variant.resolve")

    job_id = f"source-matrix-{allele_locus.rsid}"
    raw_result = call_operation("genomi.list_resources")
    job_path = background_jobs.jobs_dir() / f"{job_id}.json"
    background_jobs.write_job(
        job_path,
        {
            "job_id": job_id,
            "operation": "genomi.list_resources",
            "params": {},
            "params_digest": background_jobs.operation_params_digest("genomi.list_resources", {}),
            "status": "completed",
            "created_at": "2026-05-20T00:00:00+00:00",
            "started_at": "2026-05-20T00:00:00+00:00",
            "finished_at": "2026-05-20T00:00:01+00:00",
            "pid": 123,
            "result": raw_result,
        },
    )
    checked = call_operation("genomi.check_background_job", {"job_id": job_id})
    seen_operations.add("genomi.check_background_job")
    testcase.assertEqual(checked["status"], "completed")


def _assert_journal_operations(
    testcase: SupportsAssertions,
    contract: SourceContractCase,
    seen_operations: set[str],
) -> None:
    appended = call_operation(
        "journal.append_entry",
        {
            "entry_type": "observation",
            "content": f"{contract.expected_format} source matrix runtime observation",
            "tags": ["source-matrix", contract.expected_format],
        },
    )
    seen_operations.add("journal.append_entry")
    testcase.assertEqual(appended["status"], "completed")

    searched = call_operation("journal.search_entries", {"query": "source matrix", "limit": 5})
    seen_operations.add("journal.search_entries")
    testcase.assertTrue(searched["entries"])

    summarized = call_operation("journal.summarize", {"limit": 5})
    seen_operations.add("journal.summarize")
    testcase.assertEqual(summarized["status"], "completed")

    exported = call_operation("journal.export_memory", {})
    seen_operations.add("journal.export_memory")
    testcase.assertEqual(exported["status"], "completed")


def _assert_research_operations(
    testcase: SupportsAssertions,
    contract: SourceContractCase,
    seen_operations: set[str],
) -> None:
    listed = call_operation("research.list_sources", {"target_type": "drug"})
    seen_operations.add("research.list_sources")
    testcase.assertEqual(listed["filters"]["target_type"], "drug")
    testcase.assertGreater(listed["summary"]["source_count"], 0)

    topic = f"source-matrix-{contract.expected_format}"
    payload = {
        "target": {"type": "topic", "topic": topic},
        "source": {
            "source_id": "source-matrix",
            "title": "Source matrix fixture",
            "type": "test_fixture",
            "url": "https://example.test/source-matrix",
            "accessed_at": "2026-05-20",
        },
        "finding": {
            "type": "source_matrix_runtime_fixture",
            "text": f"{contract.expected_format} matrix reviewed research fixture.",
            "summary": f"{contract.expected_format} matrix fixture.",
        },
        "captured_by": "source-matrix-runtime-test",
    }
    recorded = call_operation("research.record", {"payload": payload, "scope": "shared"})
    seen_operations.add("research.record")
    testcase.assertEqual(recorded["status"], "completed")

    queried = call_operation("research.query", {"target_type": "topic", "topic": topic})
    seen_operations.add("research.query")
    testcase.assertEqual(queried["count"], 1)

    searched = call_operation("research.search", {"query": topic, "limit": 5})
    seen_operations.add("research.search")
    testcase.assertGreater(searched["count"], 0)

    packet = call_operation("research.build_target_packet", {"target_type": "topic", "topic": topic})
    seen_operations.add("research.build_target_packet")
    testcase.assertEqual(packet["target"]["target_type"], "topic")
    testcase.assertEqual(packet["target"]["topic"], topic)
