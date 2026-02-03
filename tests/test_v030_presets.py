import unittest
from drift.presets import get_diabetes_topology, get_diabetes_bridge, get_human_cancer_topology, get_human_cancer_bridge

class TestV030Presets(unittest.TestCase):
    def test_diabetes_preset(self):
        topo = get_diabetes_topology()
        bridge = get_diabetes_bridge()
        self.assertIn("AS160", topo.species)
        self.assertTrue(any(m["protein_name"] == "AS160" for m in bridge.mappings))

    def test_cancer_preset(self):
        topo = get_human_cancer_topology()
        bridge = get_human_cancer_bridge()
        self.assertIn("AMPK", topo.species)

if __name__ == "__main__":
    unittest.main()
