import unittest
import os
import json
import numpy as np
from drift.topology import Topology, get_default_topology, get_complex_topology, pi3k_akt_mtor_drift, complex_signaling_drift, drift_model
from unittest.mock import MagicMock, patch

class TestV030Topology(unittest.TestCase):
    def setUp(self):
        self.temp_sbml = "test_model_topo.xml"
        sbml_content = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">
  <model id="test_model" name="Test Model">
    <listOfSpecies>
      <species id="S1" compartment="c" initialConcentration="1.0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="S2" compartment="c" initialConcentration="0.0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
    <listOfParameters>
      <parameter id="k1" value="0.1" constant="true"/>
    </listOfParameters>
    <listOfReactions>
      <reaction id="r1" reversible="false" fast="false">
        <listOfReactants>
          <speciesReference species="S1" stoichiometry="1" constant="true"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="S2" stoichiometry="1" constant="true"/>
        </listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><ci> k1 </ci><ci> S1 </ci></apply>
          </math>
        </kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>"""
        with open(self.temp_sbml, "w") as f:
            f.write(sbml_content)

    def tearDown(self):
        if os.path.exists(self.temp_sbml):
            os.remove(self.temp_sbml)

    def test_sbml_loading(self):
        topology = Topology.from_sbml(self.temp_sbml)
        self.assertEqual(topology.name, "Test Model")
        self.assertIn("S1", topology.species)

    def test_default_topologies(self):
        topo = get_default_topology()
        self.assertEqual(len(topo.species), 3)
        complex_topo = get_complex_topology()
        self.assertEqual(len(complex_topo.species), 4)

    def test_pi3k_drift_variants(self):
        state = np.full(3, 0.5)
        params7 = np.array([0.1, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5])
        d7 = pi3k_akt_mtor_drift(state, params7)
        self.assertEqual(len(d7), 3)

    def test_topology_from_json(self):
        json_path = "test_topo.json"
        with open(json_path, "w") as f:
            json.dump({"species": ["A"], "parameters": {"k": 0.1}, "name": "J"}, f)
        try:
            topo = Topology.from_json(json_path)
            self.assertEqual(topo.name, "J")
        finally:
            if os.path.exists(json_path): os.remove(json_path)

    def test_sbml_qual_jit_fail_fallback(self):
        model = MagicMock(); qual_ext = MagicMock()
        s1 = MagicMock(); s1.getId.return_value = "S1"
        qual_ext.getListOfQualitativeSpecies.return_value = [s1]
        qual_ext.getListOfTransitions.return_value = []
        with patch("builtins.compile", side_effect=Exception("JIT FAIL")):
            topo = Topology._from_sbml_qual(model, qual_ext)
            self.assertIsNotNone(topo.drift_fn)

    def test_topology_rr_mock(self):
        doc = MagicMock(); doc.getNumErrors.return_value = 0
        model = MagicMock(); model.getPlugin.return_value = None
        model.getListOfSpecies.return_value = [MagicMock(getId=lambda: "S1")]
        model.getListOfParameters.return_value = []
        doc.getModel.return_value = model
        mock_rr_mod = MagicMock()
        with patch.dict("sys.modules", {"roadrunner": mock_rr_mod}):
            mock_rr_mod.RoadRunner = MagicMock()
            rr_inst = mock_rr_mod.RoadRunner.return_value
            rr_inst.model.getGlobalParameterIds.return_value = ["k1"]
            rr_inst.model.getGlobalParameterValues.return_value = [0.1]
            with patch("libsbml.SBMLReader.readSBML", return_value=doc):
                topo = Topology.from_sbml("fake.xml")
                if topo.drift_fn:
                    topo.drift_fn(np.array([0.5]), None)
                    self.assertTrue(rr_inst.reset.called)

if __name__ == "__main__":
    unittest.main()
