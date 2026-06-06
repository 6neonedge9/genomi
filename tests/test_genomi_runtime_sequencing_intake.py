from __future__ import annotations

import gzip
import os
import tempfile
import zipfile
from pathlib import Path
from unittest import mock

from genomi.operations import call_operation

from tests.support.runtime.genomi import GenomiRuntimeTestCase


def _hidden_path_leaks(payload: object, *paths: Path) -> list[str]:
    hidden = {
        value
        for path in paths
        for value in (str(path), str(path.expanduser().resolve(strict=False)))
    }
    leaks: list[str] = []

    def visit(value: object, location: str) -> None:
        if isinstance(value, dict):
            for key, item in value.items():
                visit(item, f"{location}.{key}")
            return
        if isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, f"{location}[{index}]")
            return
        if isinstance(value, str) and any(path in value for path in hidden):
            leaks.append(location)

    visit(payload, "$")
    return leaks


class GenomiRuntimeSequencingIntakeTests(GenomiRuntimeTestCase):
    def test_fastq_detect_source_recognizes_paired_inputs(self) -> None:
        from genomi.active_genome_index.alignment import detect_paired_fastq
        from genomi.active_genome_index.source_intake import detect_source

        with tempfile.TemporaryDirectory() as tmp:
            previous = os.getcwd()
            os.chdir(tmp)
            try:
                for r1_name, r2_name in [
                    ("sample_R1_001.fastq", "sample_R2_001.fastq"),
                    ("foo_1.fastq.gz", "foo_2.fastq.gz"),
                ]:
                    Path(r1_name).write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
                    Path(r2_name).write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
                    detection = detect_source(r1_name)
                    self.assertEqual(detection.source_format, "fastq")
                    self.assertEqual(detection.source_kind, "paired_reads_input")
                    pair = detect_paired_fastq(Path(r1_name))
                    self.assertIsNotNone(pair)
                    assert pair is not None
                    self.assertEqual(pair[1].name, r2_name)
                    Path(r1_name).unlink()
                    Path(r2_name).unlink()
            finally:
                os.chdir(previous)

    def test_fastq_aligner_pick_uses_median_read_length(self) -> None:
        from genomi.active_genome_index.alignment import pick_aligner_for_reads, sniff_fastq_read_length

        with tempfile.TemporaryDirectory() as tmp:
            previous = os.getcwd()
            os.chdir(tmp)
            try:
                short = Path("short_R1.fastq")
                short_seq = "ACGT" * 37 + "AC"  # 150 bp
                short.write_text(f"@s\n{short_seq}\n+\n{'I' * len(short_seq)}\n", encoding="utf-8")
                self.assertEqual(sniff_fastq_read_length(short), 150)
                self.assertEqual(pick_aligner_for_reads(sniff_fastq_read_length(short)), "bwa-mem2")

                long_path = Path("long_R1.fastq")
                long_seq = "ACGT" * 200  # 800 bp
                long_path.write_text(f"@l\n{long_seq}\n+\n{'I' * len(long_seq)}\n", encoding="utf-8")
                self.assertEqual(sniff_fastq_read_length(long_path), 800)
                self.assertEqual(pick_aligner_for_reads(sniff_fastq_read_length(long_path)), "minimap2")
            finally:
                os.chdir(previous)

    def test_fastq_parse_returns_requires_library_install_when_aligner_missing(self) -> None:
        from genomi.active_genome_index import alignment, source_intake

        with tempfile.TemporaryDirectory() as tmp:
            previous = os.getcwd()
            os.chdir(tmp)
            try:
                r1 = Path("sample_R1_001.fastq.gz")
                r2 = Path("sample_R2_001.fastq.gz")
                short_seq = "ACGT" * 37 + "AC"  # 150 bp
                record = f"@r\n{short_seq}\n+\n{'I' * len(short_seq)}\n".encode("utf-8")
                with gzip.open(r1, "wb") as handle:
                    handle.write(record)
                with gzip.open(r2, "wb") as handle:
                    handle.write(record)
                reference = Path("reference.fa")
                reference.write_text(">chr1\n" + "A" * 200 + "\n", encoding="utf-8")

                with (
                    mock.patch.object(alignment, "resolve_aligner_binary", return_value=None),
                    mock.patch.object(alignment.shutil, "which", return_value=None),
                ):
                    result = source_intake.parse_source(
                        r1,
                        reference_fasta=reference,
                        auto_reference_fasta=False,
                    )

                self.assertEqual(result["status"], "requires_library_install")
                self.assertEqual(result["source_format"], "fastq")
                binaries = {entry["binary"] for entry in result["missing_libraries"]}
                self.assertIn("bwa-mem2", binaries)
                self.assertIn("samtools", binaries)
                libs = {entry["install_library"] for entry in result["missing_libraries"]}
                self.assertIn("bwa-mem2-binary", libs)
            finally:
                os.chdir(previous)

    def test_fastq_parse_materializes_paired_reads_from_zip_archive_pair(self) -> None:
        from genomi.active_genome_index import source_intake

        with tempfile.TemporaryDirectory() as tmp:
            previous = os.getcwd()
            os.chdir(tmp)
            try:
                archive = Path("fastq-pair.zip")
                short_seq = "ACGT" * 37 + "AC"
                r1_record = f"@r1\n{short_seq}\n+\n{'I' * len(short_seq)}\n".encode("utf-8")
                r2_record = (
                    f"@r2\n{short_seq}\n+\n{'I' * len(short_seq)}\n"
                    f"@r2b\n{short_seq}\n+\n{'I' * len(short_seq)}\n"
                ).encode("utf-8")
                with zipfile.ZipFile(archive, "w") as bundle:
                    bundle.writestr("reads/sample_R1_001.fastq.gz", gzip.compress(r1_record))
                    bundle.writestr("reads/sample_R2_001.fastq.gz", gzip.compress(r2_record))
                reference = Path("reference.fa")
                reference.write_text(">chr1\n" + "A" * 200 + "\n", encoding="utf-8")
                detection = source_intake.detect_source(archive)
                self.assertEqual(detection.source_format, "fastq")
                self.assertEqual(detection.member_name, "reads/sample_R1_001.fastq.gz")

                with mock.patch(
                    "genomi.active_genome_index.source_intake.sequencing.align_fastq_to_bam",
                    return_value={
                        "status": "requires_library_install",
                        "missing_libraries": [{"binary": "bwa-mem2", "install_library": "bwa-mem2-binary"}],
                        "message": "stubbed missing aligner",
                    },
                ) as align:
                    result = source_intake.parse_source(
                        archive,
                        reference_fasta=reference,
                        auto_reference_fasta=False,
                    )

                self.assertEqual(result["status"], "requires_library_install")
                self.assertEqual(result["source_format"], "fastq")
                r1_arg, r2_arg = (Path(value) for value in align.call_args.args[:2])
                self.assertNotEqual(r1_arg.parent, archive.parent)
                self.assertEqual(r1_arg.read_text(encoding="utf-8"), r1_record.decode())
                self.assertEqual(r2_arg.read_text(encoding="utf-8"), r2_record.decode())
                self.assertEqual(Path(result["fastq"]["r1"]), r1_arg)
                self.assertEqual(Path(result["fastq"]["r2"]), r2_arg)
            finally:
                os.chdir(previous)

    def test_fastq_parse_raises_when_r2_sibling_missing(self) -> None:
        from genomi.active_genome_index import source_intake

        with tempfile.TemporaryDirectory() as tmp:
            previous = os.getcwd()
            os.chdir(tmp)
            try:
                r1 = Path("orphan_R1_001.fastq")
                r1.write_text("@r\nACGT\n+\nIIII\n", encoding="utf-8")
                reference = Path("reference.fa")
                reference.write_text(">chr1\nA\n", encoding="utf-8")

                with self.assertRaises(ValueError) as ctx:
                    source_intake.parse_source(
                        r1,
                        reference_fasta=reference,
                        auto_reference_fasta=False,
                    )
                self.assertIn("paired-end R1", str(ctx.exception))
            finally:
                os.chdir(previous)

    def test_active_genome_index_parse_accepts_bam_by_materializing_derived_vcf(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            previous = os.getcwd()
            os.chdir(tmp)
            try:
                bam = Path("sample.bam")
                bam.write_bytes(b"BAM\x01")
                reference = Path("reference.fa")
                reference.write_text(">chr1\n" + "A" * 200 + "\n", encoding="utf-8")

                def fake_materialize_bam_variant_vcf(
                    bam_path: Path,
                    reference_fasta: Path,
                    output_vcf: Path,
                    *,
                    force: bool = False,
                ) -> dict[str, object]:
                    del bam_path, reference_fasta, force
                    Path(output_vcf).write_text(
                        "##fileformat=VCFv4.2\n"
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
                        "1\t100\trs555\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t0/1:31:99\n",
                        encoding="utf-8",
                    )
                    return {
                        "status": "completed",
                        "output": str(output_vcf),
                        "manifest_path": str(Path(f"{output_vcf}.genomi-manifest.json")),
                    }

                with (
                    mock.patch("genomi.active_genome_index.source_intake.sequencing.infer_genome_build_from_bam", return_value="GRCh38"),
                    mock.patch(
                        "genomi.active_genome_index.source_intake.sequencing.materialize_bam_variant_vcf",
                        side_effect=fake_materialize_bam_variant_vcf,
                    ),
                ):
                    self.approve_access()
                    parsed = call_operation(
                        "genomi.parse_source",
                        {"source": str(bam), "reference_fasta": str(reference)},
                    )

                self.assertEqual(parsed["status"], "completed")
                self.assertEqual(parsed["source_format"], "bam")
                self.assertEqual(parsed["source_kind"], "alignment_reads")
                self.assertIn("active_genome_index", parsed)
                self.assertEqual([step["name"] for step in parsed["steps"]], ["init-source", "materialize-variants-from-bam", "build-active-genome-index-from-derived-vcf"])
                self.assertEqual(set(parsed["outputs"]), {"agi_path", "bam_variant_call_manifest", "derived_vcf"})
                self.assertTrue(Path(parsed["outputs"]["derived_vcf"]).exists())
                self.assertIsNone(parsed.get("genotype_reference_fasta"))
                self.assertIsNone(parsed.get("evidence_summary"))
                self.assertEqual(_hidden_path_leaks(parsed, bam), [])

                current = call_operation("genomi.describe_context")
                self.assertTrue(current["has_active_genome_index"])
                self.assertEqual(current["active_genome_index"]["agi_source_format"], "bam")
                self.assertTrue(current["active_genome_index"]["digitized"])
                self.assertEqual(_hidden_path_leaks(current, bam), [])
            finally:
                os.chdir(previous)

    def test_active_genome_index_parse_accepts_zipped_bam_member(self) -> None:
        import pysam

        with tempfile.TemporaryDirectory() as tmp:
            previous = os.getcwd()
            os.chdir(tmp)
            try:
                archive = Path("alignment.zip")
                bam = Path("sample.bam")
                header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "1", "LN": 1000}]}
                with pysam.AlignmentFile(str(bam), "wb", header=header):
                    pass
                original_bam_bytes = bam.read_bytes()
                with zipfile.ZipFile(archive, "w") as bundle:
                    bundle.write(bam, "nested/sample.bam")
                reference = Path("reference.fa")
                reference.write_text(">chr1\n" + "A" * 200 + "\n", encoding="utf-8")

                seen_bam_paths: list[Path] = []

                def fake_materialize_bam_variant_vcf(
                    bam_path: Path,
                    reference_fasta: Path,
                    output_vcf: Path,
                    *,
                    force: bool = False,
                ) -> dict[str, object]:
                    del reference_fasta, force
                    seen_bam_paths.append(Path(bam_path))
                    self.assertEqual(Path(bam_path).read_bytes(), original_bam_bytes)
                    Path(output_vcf).write_text(
                        "##fileformat=VCFv4.2\n"
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
                        "1\t100\trs555\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t0/1:31:99\n",
                        encoding="utf-8",
                    )
                    return {
                        "status": "completed",
                        "output": str(output_vcf),
                        "manifest_path": str(Path(f"{output_vcf}.genomi-manifest.json")),
                    }

                with (
                    mock.patch("genomi.active_genome_index.source_intake.sequencing.infer_genome_build_from_bam", return_value="GRCh38"),
                    mock.patch(
                        "genomi.active_genome_index.source_intake.sequencing.materialize_bam_variant_vcf",
                        side_effect=fake_materialize_bam_variant_vcf,
                    ),
                ):
                    self.approve_access()
                    parsed = call_operation(
                        "genomi.parse_source",
                        {"source": str(archive), "reference_fasta": str(reference)},
                    )

                self.assertEqual(parsed["status"], "completed")
                self.assertEqual(parsed["source_format"], "bam")
                self.assertIn("materialize-bam-archive-member", [step["name"] for step in parsed["steps"]])
                self.assertEqual(len(seen_bam_paths), 1)
                self.assertNotEqual(seen_bam_paths[0], archive)
                lookup = call_operation("variant.resolve", {"rsid": "rs555"})
                self.assertEqual(lookup["sample_context"]["count"], 1)
            finally:
                os.chdir(previous)

if __name__ == "__main__":
    import unittest

    unittest.main()
