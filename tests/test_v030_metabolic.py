import unittest
import numpy as np
from drift.metabolic import MetabolicBridge, BridgeBuilder, DFBASolver, sigmoidal_mapping, michaelis_menten_mapping, _is_common_pattern_match, fuzzy_match_reaction_id, DesyncError
from unittest.mock import MagicMock, patch

class TestV030Metabolic(unittest.TestCase):
    def test_sigmoidal_mapping(self):
        self.assertAlmostEqual(sigmoidal_mapping(0.5), 0.5)

    def test_mm_mapping(self):
        self.assertEqual(michaelis_menten_mapping(1.0, km=0.5), 1.0)
        self.assertEqual(michaelis_menten_mapping(0.0), 0.0)

    def test_pattern_match(self):
        self.assertTrue(_is_common_pattern_match("abc", "EX_abc"))
        self.assertFalse(_is_common_pattern_match("abc", "xyz"))

    def test_bridge_builder_fluent(self):
        bridge = BridgeBuilder().set_flux_unit("custom").set_strict_mapping(True).build()
        self.assertEqual(bridge.flux_unit, "custom")
        self.assertTrue(bridge.strict_mapping)

    def test_bridge_get_feedback(self):
        bridge = MetabolicBridge(species_names=["S1"])
        bridge.reverse_mappings = [{"flux_id": "F1", "species_name": "S1", "influence": "negative", "weight": 1.0}]
        fb = bridge.get_feedback({"F1": 10.0, "growth": 0.2})
        self.assertAlmostEqual(fb[0], 0.0)

    def test_dfba_solver_proxy(self):
        solver = DFBASolver(model_name="textbook")
        if not hasattr(solver, "proxy_params"):
            solver.proxy_params = {"base_growth": 0.2, "yield_coefficient": 0.02}
        res = solver._solve_proxy({"R1": -10.0}, scalings={"R1": 0.5})
        self.assertEqual(res["status"], "proxied")

    def test_dfba_solver_validation(self):
        solver = DFBASolver(model_name="textbook")
        if not solver.headless:
            solver.validate_model()
            m = MagicMock(); m.objective = None
            orig = solver.model; solver.model = m
            with self.assertRaises(RuntimeError): solver.validate_model()
            solver.model = orig

    def test_bridge_desync(self):
        bridge = MetabolicBridge(species_names=["A"])
        with self.assertRaises(DesyncError): bridge.get_constraints(np.array([0.5, 0.5]))

if __name__ == "__main__":
    unittest.main()
