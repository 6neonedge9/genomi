from __future__ import annotations

import unittest
from importlib import resources as importlib_resources
from pathlib import Path

from genomi.capabilities.pharmacogenomics import pgx_requirements, pgx_star

REPO_ROOT = Path(__file__).resolve().parents[1]


class ImplementationCatalogTests(unittest.TestCase):
    def test_marker_panel_infrastructure_is_removed(self) -> None:
        # The marker_panel.* operations + panel.py module + builtin_panels.json
        # were removed when the nutrigenomics capability landed. The curated
        # marker context that lived there now lives in
        # src/genomi/capabilities/nutrigenomics/catalog.py.
        self.assertFalse((REPO_ROOT / "src/genomi/active_genome_index/panel.py").exists())
        self.assertFalse((REPO_ROOT / "src/genomi/active_genome_index/panels").exists())
        self.assertFalse((REPO_ROOT / "skills/panels").exists())

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

if __name__ == "__main__":
    unittest.main()
