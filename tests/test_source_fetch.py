from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from genomi.runtime import source_fetch
from genomi.runtime.external import write_manifest


class SourceFetchTests(unittest.TestCase):
    def test_cached_refresh_uses_freshness_timeout(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "library.dat"
            output.write_text("cached", encoding="utf-8")
            manifest = output.with_name(output.name + ".genomi-manifest.json")
            expected = {
                "library": "example",
                "source_url": "https://example.test/library.dat",
                "output": str(output),
            }
            write_manifest(manifest, {**expected, "etag": "previous"})

            with mock.patch.object(
                source_fetch,
                "conditional_fetch",
                return_value={"status": "up_to_date", "validators": {"etag": "previous"}, "reason": "not_modified"},
            ) as conditional:
                result = source_fetch.refresh_or_download(
                    expected["source_url"],
                    output,
                    manifest,
                    expected=expected,
                    force=False,
                    refresh=True,
                    timeout=120,
                )

        self.assertEqual(result["status"], "up_to_date")
        self.assertEqual(result["freshness"], "not_modified")
        self.assertEqual(conditional.call_args.kwargs["timeout"], source_fetch.DEFAULT_FRESHNESS_TIMEOUT_SECONDS)

    def test_missing_file_uses_download_timeout(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "library.dat"
            manifest = output.with_name(output.name + ".genomi-manifest.json")
            expected = {
                "library": "example",
                "source_url": "https://example.test/library.dat",
                "output": str(output),
            }

            def _download(_url, target, *, timeout, user_agent=None):
                Path(target).write_text("downloaded", encoding="utf-8")
                return {"etag": "new"}

            with mock.patch.object(source_fetch, "download", side_effect=_download) as download:
                result = source_fetch.refresh_or_download(
                    expected["source_url"],
                    output,
                    manifest,
                    expected=expected,
                    force=False,
                    refresh=True,
                    timeout=120,
                )

        self.assertEqual(result["status"], "downloaded")
        self.assertEqual(download.call_args.kwargs["timeout"], 120)

    def test_existing_file_without_manifest_is_downloaded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "library.dat"
            output.write_text("untracked cached bytes", encoding="utf-8")
            manifest = output.with_name(output.name + ".genomi-manifest.json")
            expected = {
                "library": "example",
                "source_url": "https://example.test/library.dat",
                "output": str(output),
            }

            def _download(_url, target, *, timeout, user_agent=None):
                Path(target).write_text("downloaded", encoding="utf-8")
                return {"etag": "new"}

            with mock.patch.object(source_fetch, "download", side_effect=_download) as download:
                result = source_fetch.refresh_or_download(
                    expected["source_url"],
                    output,
                    manifest,
                    expected=expected,
                    force=False,
                    refresh=False,
                )

            self.assertEqual(result["status"], "downloaded")
            self.assertEqual(output.read_text(encoding="utf-8"), "downloaded")
            self.assertEqual(download.call_count, 1)


if __name__ == "__main__":
    unittest.main()
