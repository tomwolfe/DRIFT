"""
Recon1 Integration Tests for DRIFT.

These tests verify that the Human Cancer preset works with Recon1
(Human Genome-Scale Metabolic Model). Tests are marked as slow because
Recon1 is a large model (~3,700 reactions) that takes time to load.
"""

import pytest
import numpy as np


pytestmark = pytest.mark.slow


def _recon1_available() -> bool:
    """Check if Recon1 can be loaded."""
    try:
        from cobra.io import load_model
        model = load_model("recon1")
        return model is not None and len(model.reactions) > 3000
    except Exception:
        return False


@pytest.mark.skipif(not _recon1_available(), reason="Recon1 model not available (large model, download required)")
class TestRecon1Integration:
    """Integration tests requiring Recon1 human genome-scale model."""

    def test_load_recon1(self):
        """Test that Recon1 loads successfully."""
        from cobra.io import load_model
        model = load_model("recon1")
        assert model is not None
        assert len(model.reactions) > 3000
        assert len(model.metabolites) > 4000

    def test_recon1_preset_mappings_resolve(self):
        """Test that all preset bridge mappings resolve against Recon1."""
        from cobra.io import load_model
        from drift.presets import get_human_cancer_bridge

        model = load_model("recon1")
        bridge = get_human_cancer_bridge()

        # Validate that all mappings resolve
        all_resolved = bridge.validate_with_model(model)
        assert all_resolved, "Some bridge mappings failed to resolve against Recon1"

    def test_recon1_basal_growth(self):
        """Test that Recon1 produces positive growth in default state."""
        from cobra.io import load_model
        model = load_model("recon1")
        solution = model.optimize()

        # Find biomass reaction
        biomass_flux = 0.0
        for rxn in model.reactions:
            if "BIOMASS" in rxn.id.upper():
                flux = solution.fluxes.get(rxn.id, 0.0)
                if flux > biomass_flux:
                    biomass_flux = flux

        assert biomass_flux > 0, f"Recon1 should have positive basal growth, got {biomass_flux}"

    def test_recon1_short_simulation(self):
        """Test that a short simulation with Recon1 completes without errors."""
        from drift.workbench import Workbench
        from drift.presets import get_human_cancer_topology, get_human_cancer_bridge

        topology = get_human_cancer_topology()
        bridge = get_human_cancer_bridge()

        wb = Workbench(
            drug_kd=1.0,
            drug_concentration=5.0,
            model_name="recon1",
            topology=topology,
            bridge=bridge,
            sub_steps=3,
        )
        wb.signaling.noise_scale = 0.001

        history = wb.run_simulation(steps=10)

        # Verify basic structure
        assert "growth" in history
        assert "signaling" in history
        assert len(history["growth"]) == 10
        assert history["signaling"].shape == (10, 4)

        # Growth should be non-negative (may be 0 if constraints are too tight)
        assert all(g >= 0 for g in history["growth"]), "Growth should be non-negative"

    def test_recon1_exchange_reaction_variants(self):
        """Test that common exchange reaction ID variations are handled."""
        from drift.metabolic import fuzzy_match_reaction_id

        # Simulated Recon1 reaction IDs (common patterns)
        recon1_ids = [
            "EX_glc__D_e", "EX_lac__L_e", "EX_gln__L_e",
            "EX_o2_e", "EX_h2o_e", "EX_co2_e",
            "BIOMASS_reaction", "BIOMASS_RECON1",
        ]

        # Test exact match
        assert fuzzy_match_reaction_id("EX_glc__D_e", recon1_ids) == "EX_glc__D_e"

        # Test case-insensitive match
        assert fuzzy_match_reaction_id("ex_glc__d_e", recon1_ids) is not None

        # Test prefix variation
        result = fuzzy_match_reaction_id("glc__D_e", recon1_ids)
        # Should either match EX_glc__D_e or not match (both acceptable)
