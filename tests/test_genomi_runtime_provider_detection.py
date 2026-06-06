from __future__ import annotations

import tempfile
from pathlib import Path

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


class GenomiRuntimeProviderDetectionTests(GenomiRuntimeTestCase):
    def test_detect_source_tags_vcf_provider_for_sequencingdotcom(self) -> None:
        from genomi.active_genome_index.source_intake import detect_source

        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "sample.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "##source=Sequencing.com (30x WGS)\n"
                "##dataAnalysisProvider=Sequencing.com\n"
                "##reference=GRCh38.p13\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNG1ABCDEFG\n",
                encoding="utf-8",
            )
            detection = detect_source(vcf)
            self.assertEqual(detection.source_format, "vcf")
            self.assertEqual(detection.provider, "sequencingdotcom")
            self.assertEqual(detection.reference_build, "GRCh38")

    def test_parse_source_surfaces_vcf_provider(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "sample.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "##source=Sequencing.com (30x WGS)\n"
                "##dataAnalysisProvider=Sequencing.com\n"
                "##reference=GRCh38.p13\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNG1ABCDEFG\n"
                "1\t100\trs900000001\tA\tC\t50\tPASS\t.\tGT\t0/1\n",
                encoding="utf-8",
            )

            parsed = call_operation("genomi.parse_source", {"source": str(vcf)})

        self.assertEqual(parsed["status"], "completed")
        self.assertEqual(parsed["provider"], "sequencingdotcom")
        self.assertEqual(parsed["active_genome_index"]["agi_source_provider"], "sequencingdotcom")

    def test_detect_source_tags_vcf_provider_for_dantelabs(self) -> None:
        from genomi.active_genome_index.source_intake import detect_source

        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "DTC7U778.raw.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                '##DRAGENCommandLine=<ID=dragen,Version="SW: 05.121.645.4.0.3">\n'
                "##reference=file:///references/grch37/reference.bin\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDTC7U778\n",
                encoding="utf-8",
            )
            detection = detect_source(vcf)
            self.assertEqual(detection.source_format, "vcf")
            self.assertEqual(detection.provider, "dantelabs")
            self.assertEqual(detection.reference_build, "GRCh37")

    def test_detect_source_tags_vcf_provider_for_nebula_via_sample_id(self) -> None:
        from genomi.active_genome_index.source_intake import detect_source

        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "NG176JZTG8.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "##reference=file:///mnt/ssd/MegaBOLT_scheduler/reference/hg38.fa\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNG176JZTG8\n",
                encoding="utf-8",
            )
            detection = detect_source(vcf)
            self.assertEqual(detection.source_format, "vcf")
            self.assertEqual(detection.provider, "nebula")
            self.assertEqual(detection.reference_build, "GRCh38")
