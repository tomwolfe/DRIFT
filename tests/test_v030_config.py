import unittest
from drift.config import SimulationConfig

class TestV030Config(unittest.TestCase):
    def test_web_schema(self):
        config = SimulationConfig(drug_kd=0.1)
        schema = config.to_web_schema()
        self.assertEqual(schema["schema_version"], "0.3.0")
        params = {p["id"]: p for p in schema["parameters"]}
        self.assertEqual(params["drug_kd"]["default"], 0.1)

if __name__ == "__main__":
    unittest.main()
