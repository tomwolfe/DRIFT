import unittest
import os
import numpy as np
from drift.topology import Topology

class TestSBMLIntegration(unittest.TestCase):
    def setUp(self):
        self.temp_sbml = "test_model_integration.xml"
        # Minimal SBML model for testing
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
            <apply>
              <times/>
              <ci> k1 </ci>
              <ci> S1 </ci>
            </apply>
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

    def test_sbml_loading_success(self):
        """Verify that SBML models can be loaded into a Topology object."""
        topology = Topology.from_sbml(self.temp_sbml)
        self.assertEqual(topology.name, "Test Model")
        self.assertIn("S1", topology.species)
        self.assertIn("S2", topology.species)
        self.assertEqual(topology.parameters["k1"], 0.1)

    def test_sbml_with_inhibited_species(self):
        """Verify that inhibited species can be specified during SBML loading."""
        topology = Topology.from_sbml(self.temp_sbml, inhibited_species="S1")
        self.assertEqual(topology.inhibited_species, "S1")

if __name__ == "__main__":
    unittest.main()
