from __future__ import annotations

import unittest
from importlib import resources as importlib_resources
from pathlib import Path

from genomi.capabilities.pharmacogenomics import pgx_requirements, pgx_star


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src" / "genomi"
MIGRATION_PROMISE_TERMS = (
    "former single-module",
    "former single-file",
    "former monolithic",
    "keeps working unchanged",
    "keep working unchanged",
    "continue to work unchanged",
    "existing import paths continue",
    "preserves the public surface",
    "preserves the full public surface",
)


class ImplementationCatalogTests(unittest.TestCase):
    def test_pgx_marker_definitions_and_requirements_are_packaged_data(self) -> None:
        marker_resource = importlib_resources.files("genomi.capabilities.pharmacogenomics").joinpath("data").joinpath("star_marker_definitions.json")
        requirement_resource = importlib_resources.files("genomi.capabilities.pharmacogenomics").joinpath("data").joinpath("gene_requirements.json")
        self.assertTrue(marker_resource.is_file())
        self.assertTrue(requirement_resource.is_file())

        marker_catalog = pgx_star.marker_definition_catalog()
        requirement_catalog = pgx_requirements.gene_requirements_catalog()
        self.assertTrue(marker_catalog["curation_scope"].strip())
        self.assertTrue(marker_catalog["definition_sets"])
        self.assertTrue(requirement_catalog["curation_scope"].strip())
        self.assertTrue(requirement_catalog["named_allele_matcher_genes"])
        self.assertEqual(pgx_star.implemented_marker_definition_genes(), ["CYP2C19"])
        self.assertIn("CYP2D6", requirement_catalog["outside_call_genes"])

    def test_source_docs_do_not_publish_migration_promises(self) -> None:
        violations: list[str] = []
        for path in sorted(SRC_ROOT.rglob("*.py")):
            text = path.read_text(encoding="utf-8").lower()
            for term in MIGRATION_PROMISE_TERMS:
                if term in text:
                    violations.append(f"{path.relative_to(REPO_ROOT)}: {term}")
        self.assertEqual(violations, [])


if __name__ == "__main__":
    unittest.main()
