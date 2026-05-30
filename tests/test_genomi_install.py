from __future__ import annotations

from unittest import mock

from genomi.active_genome_index.active_genome_index import (
    connect_existing,
    create_active_genome_index,
)
from genomi.active_genome_index._agi_schema import _upsert_metadata
from genomi.interfaces.cli import build_parser
from genomi.operations import call_operation
from genomi.runtime import context as runtime_context

from _genomi_runtime_helpers import GenomiRuntimeTestCase


def _register_stale_genome(home, *, stored_schema: int = 1):
    """Build a real index, mark its stored schema stale, and register it."""
    vcf = home / "stale.vcf"
    vcf.parent.mkdir(parents=True, exist_ok=True)
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    index = vcf.with_suffix(".sqlite")
    create_active_genome_index(vcf, index)
    with connect_existing(index) as connection:
        _upsert_metadata(connection, "schema_version", stored_schema)
        connection.commit()
    runtime_context.set_active_genome_index(
        vcf, status="parsed", active_genome_index_path=index, genome_build="GRCh38"
    )
    return vcf, index


class GenomiInstallTests(GenomiRuntimeTestCase):
    def test_bare_install_defaults_to_everything(self) -> None:
        # `genomi install` / `genomi update` with no flags updates everything:
        # libraries default to 'everything' (runtime + reparse-stale are set by
        # the command handler, not via flags). setup-only stays available opt-in.
        args = build_parser().parse_args(["install"])
        self.assertEqual(args.libraries, "everything")
        self.assertFalse(args.force)

    def test_update_is_a_cli_alias_of_install_not_a_separate_command(self) -> None:
        install = build_parser().parse_args(["install", "--libraries", "everything"])
        update = build_parser().parse_args(["update", "--libraries", "everything"])
        # Same handler, same defaults — a true alias, not a duplicate command.
        self.assertIs(update.func, install.func)
        self.assertEqual(build_parser().parse_args(["update"]).libraries, "everything")

    def test_cli_install_requests_everything(self) -> None:
        # The command front door asks the operation to update everything; the
        # operation's own defaults stay conservative (tested separately).
        captured: dict[str, object] = {}

        def _capture(operation: str, params: dict[str, object]) -> dict[str, object]:
            captured["operation"] = operation
            captured["params"] = params
            return {"status": "completed"}

        args = build_parser().parse_args(["update"])
        with mock.patch("genomi.interfaces.cli.call_operation", _capture):
            args.func(args)

        self.assertEqual(captured["operation"], "genomi.install")
        params = captured["params"]
        self.assertEqual(params["libraries"], "everything")
        self.assertTrue(params["update_runtime"])
        self.assertTrue(params["reparse_stale"])

    def test_no_separate_update_tool_and_install_description_covers_update(self) -> None:
        from genomi.operations.registry.table import OPERATIONS

        by_name = {op.name: op for op in OPERATIONS}
        # The alias lives in wording, not a duplicate MCP tool.
        self.assertNotIn("genomi.update", by_name)
        description = by_name["genomi.install"].tool_definition()["description"]
        self.assertIn("update", description.lower())
        # A bare update call must be valid: libraries is not a required field.
        self.assertNotIn(
            "libraries",
            by_name["genomi.install"].tool_definition()["inputSchema"].get("required", []),
        )

    def test_setup_only_install_persists_response_profile_without_context_disclosure(self) -> None:
        result = call_operation("genomi.install", {"libraries": "setup-only", "response_profile": "expert"})

        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["install"]["status"], "skipped")
        self.assertEqual(result["active_response_profile"]["id"], "expert")
        self.assertEqual(result["install_scope"]["updates"][0], "genomi_home_setup")
        # A bare operation call does not request a runtime git pull.
        self.assertEqual(result["runtime_update"]["status"], "not_requested")

    def test_library_inventory_points_to_genomi_install_command(self) -> None:
        result = call_operation("genomi.install", {"libraries": "setup-only"})
        commands = [item["install_command"] for item in result["library_inventory"]["libraries"]]

        self.assertTrue(commands)
        self.assertTrue(all(command.startswith("genomi install --libraries ") for command in commands))

    def test_bare_operation_call_does_not_reparse_or_pull(self) -> None:
        # The operation stays conservative by default (only the CLI front door
        # turns everything on): no git pull, no reparse scan.
        _register_stale_genome(self.genomi_home)
        result = call_operation("genomi.install", {"libraries": "setup-only"})
        self.assertIsNone(result["reparse"])
        self.assertEqual(result["runtime_update"]["status"], "not_requested")

    def test_reparse_stale_launches_background_job_per_genome(self) -> None:
        vcf, _index = _register_stale_genome(self.genomi_home)

        launched: list[tuple[str, dict]] = []

        def _fake_start(operation: str, params: dict) -> dict:
            launched.append((operation, params))
            return {"job_id": "job-x", "job_path": "/tmp/job-x.json"}

        with mock.patch(
            "genomi.runtime.background_jobs.start_operation_job", side_effect=_fake_start
        ):
            result = call_operation(
                "genomi.install", {"libraries": "setup-only", "reparse_stale": True}
            )

        reparse = result["reparse"]
        self.assertEqual(reparse["stale"], 1)
        self.assertEqual(len(reparse["launched"]), 1)
        self.assertEqual(reparse["launched"][0]["source"], str(vcf))
        self.assertEqual(launched, [("genomi.parse_source", {"source": str(vcf), "force": True})])

    def test_reparse_skips_genome_whose_source_is_gone(self) -> None:
        vcf, _index = _register_stale_genome(self.genomi_home)
        vcf.unlink()  # source no longer available — cannot rebuild

        with mock.patch("genomi.runtime.background_jobs.start_operation_job") as start:
            result = call_operation(
                "genomi.install", {"libraries": "setup-only", "reparse_stale": True}
            )

        start.assert_not_called()
        reparse = result["reparse"]
        self.assertEqual(reparse["stale"], 1)
        self.assertEqual(reparse["launched"], [])
        self.assertEqual(reparse["skipped"][0]["reason"], "source_unavailable")

    def test_skip_env_suppresses_runtime_git_pull(self) -> None:
        # The gate's whole reason for existing: a non-git distribution sets it so
        # `genomi update` never tries to git pull.
        with mock.patch.dict("os.environ", {"GENOMI_SKIP_RUNTIME_GIT_PULL": "1"}):
            result = call_operation(
                "genomi.install", {"libraries": "setup-only", "update_runtime": True}
            )
        runtime_update = result["runtime_update"]
        self.assertEqual(runtime_update["status"], "skipped")
        self.assertFalse(runtime_update["restart_required"])
        self.assertIn("GENOMI_SKIP_RUNTIME_GIT_PULL", runtime_update["message"])

    def test_legacy_runtime_update_env_still_suppresses_git_pull(self) -> None:
        # Backward-compat: an existing install that set the retired command env
        # must not suddenly start pulling.
        with mock.patch.dict("os.environ", {"GENOMI_RUNTIME_UPDATE": "anything"}):
            result = call_operation(
                "genomi.install", {"libraries": "setup-only", "update_runtime": True}
            )
        self.assertEqual(result["runtime_update"]["status"], "skipped")
        self.assertIn("GENOMI_RUNTIME_UPDATE", result["runtime_update"]["message"])
