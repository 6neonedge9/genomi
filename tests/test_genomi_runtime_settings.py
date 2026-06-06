from __future__ import annotations

import json
import tempfile
from pathlib import Path
from unittest import mock

from genomi.interfaces.mcp import handle_request
from genomi.operations import OperationError, call_operation

from tests.support.runtime.genomi import (
    GenomiRuntimeTestCase,
)


class GenomiRuntimeSettingsTests(GenomiRuntimeTestCase):
    def test_mcp_parse_default_disclosure_hides_artifact_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "sample.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                "10\t94761900\trs4244285\tG\tA\t50\tPASS\t.\tGT:DP:GQ\t0/1:31:99\n",
                encoding="utf-8",
            )
            response = handle_request(
                {
                    "jsonrpc": "2.0",
                    "id": 4,
                    "method": "tools/call",
                    "params": {
                        "name": "genomi.parse_source",
                        "arguments": {"source": str(vcf), "genome_build": "GRCh38"},
                    },
                }
            )

        self.assertIsNotNone(response)
        assert response is not None
        text = response["result"]["content"][0]["text"]
        payload = json.loads(text)
        self.assertNotIn("disclosure", payload)
        self.assertEqual(payload["steps"][0]["result"]["stats"]["total_records"], 1)
        self.assertNotIn("outputs", payload)
        self.assertNotIn("project_dir", payload)
        self.assertNotIn(str(vcf.resolve(strict=False)), text)

    def test_parse_blocker_preserves_missing_libraries_without_marking_digitized(self) -> None:
        source = self.genomi_home / "sample_R1_001.fastq.gz"
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_bytes(b"stub")
        agi_path = self.genomi_home / "blocked" / "active-genome-index.sqlite"

        with mock.patch(
            "genomi.operations.registry.handlers_admin.source_intake.parse_source",
            return_value={
                "status": "requires_library_install",
                "source": str(source),
                "source_format": "fastq",
                "source_kind": "sequencing_reads",
                "sample_slug": "fastq-sha256-blocked",
                "genome_build": "GRCh38",
                "missing_libraries": [{"binary": "samtools", "install_library": "samtools"}],
                "message": "samtools is required",
                "outputs": {"agi_path": str(agi_path)},
                "steps": [
                    {
                        "name": "align-fastq-to-bam",
                        "status": "requires_library_install",
                        "result": {
                            "status": "requires_library_install",
                            "missing_libraries": [{"binary": "samtools", "install_library": "samtools"}],
                            "message": "samtools is required",
                        },
                    }
                ],
            },
        ):
            parsed = call_operation("genomi.parse_source", {"source": str(source)})

        self.assertEqual(parsed["status"], "requires_library_install")
        self.assertEqual(parsed["missing_libraries"][0]["binary"], "samtools")
        self.assertEqual(parsed["steps"][0]["result"]["missing_libraries"][0]["install_library"], "samtools")
        active = parsed["active_genome_index"]
        self.assertEqual(active["status"], "requires_library_install")
        self.assertFalse(active["digitized"])
        self.assertEqual(active["active_genome_index_readiness"]["status"], "missing")
        self.assertFalse(active["active_genome_index_readiness"]["complete"])

    def test_describe_context_surfaces_active_response_profile_default(self) -> None:
        context = call_operation("genomi.describe_context")
        profile = context.get("active_response_profile")
        self.assertIsInstance(profile, dict)
        self.assertEqual(profile["id"], "eli5")
        self.assertEqual(profile["source"], "default")
        self.assertTrue(profile["guidance"].strip())
        self.assertTrue(profile["label"].strip())

    def test_set_response_profile_persists_and_surfaces(self) -> None:
        from genomi.runtime.host_response import host_response_profiles

        result = call_operation("genomi.set_response_profile", {"profile": "literate"})
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["active_response_profile"]["id"], "literate")
        self.assertEqual(result["active_response_profile"]["source"], "explicit")

        catalog = host_response_profiles()
        literate_entry = next(
            profile
            for profile in catalog["profiles"]
            if isinstance(profile, dict) and profile.get("id") == "literate"
        )

        context = call_operation("genomi.describe_context")
        profile = context["active_response_profile"]
        self.assertEqual(profile["id"], "literate")
        self.assertEqual(profile["source"], "explicit")
        self.assertEqual(profile["guidance"], literate_entry["guidance"])
        self.assertEqual(profile["label"], literate_entry["label"])

    def test_set_response_profile_rejects_invalid_id(self) -> None:
        with self.assertRaises(OperationError) as raised:
            call_operation("genomi.set_response_profile", {"profile": "nope"})
        self.assertEqual(raised.exception.code, "invalid_response_profile")
if __name__ == "__main__":
    import unittest

    unittest.main()
