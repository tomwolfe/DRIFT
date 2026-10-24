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

    def test_parameter_range_validation_flags_out_of_range(self):
        topo = Topology(
            species=["S1"],
            parameters={"ok": 1.0, "neg": -0.1, "big": 11.0},
        )
        # The > 10 branch and the < 0 branch are separate verdicts; a flipped
        # bool literal in either would have reported an implausible value as
        # valid.
        self.assertEqual(
            topo.validate_parameter_ranges(),
            {"ok": True, "neg": False, "big": False},
        )

    def test_pi3k_akt_mtor_drift_matches_equation(self):
        state = np.full(3, 0.5)
        params7 = np.array([0.1, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5])
        d = pi3k_akt_mtor_drift(state, params7)
        # effective_pi3k = 0.5 * (1 - 0.5); dpi3k = 0.1 - 0.1 * 0.5
        expected = np.array([0.05, 0.0125, 0.075])
        np.testing.assert_allclose(d, expected, rtol=0, atol=1e-12)

    def test_complex_signaling_drift_matches_equation(self):
        state = np.full(4, 0.5)
        params = np.array(
            [0.1, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1, 0.0, 0.0, 0.0, 0.0]
        )
        d = complex_signaling_drift(state, params)
        expected = np.array([
            0.1 / 3.5 - 0.05,   # mtor_feedback = 1 / (1 + 5 * 0.5)
            0.125 - 0.05,
            0.125 / 6.0 - 0.05,  # ampk_inhibition = 1 / (1 + 10 * 0.5)
            -0.05,
        ])
        np.testing.assert_allclose(d, expected, rtol=0, atol=1e-12)

    def test_sbml_qual_compiles_the_generated_drift(self):
        model = MagicMock()
        qual_ext = MagicMock()
        s1 = MagicMock()
        s1.getId.return_value = "S1"
        qual_ext.getListOfQualitativeSpecies.return_value = [s1]
        out = MagicMock()
        out.getQualitativeSpecies.return_value = "S1"
        trans = MagicMock()
        trans.getListOfOutputs.return_value = [out]
        default_term = MagicMock()
        default_term.getResultLevel.return_value = 0.0
        ft = MagicMock()
        ft.getMath.return_value = MagicMock()
        ft.getResultLevel.return_value = 1.0
        terms = MagicMock()
        terms.getDefaultTerm.return_value = default_term
        terms.__iter__ = lambda self: iter([ft])
        trans.getListOfFunctionTerms.return_value = terms
        qual_ext.getListOfTransitions.return_value = [trans]
        with patch("libsbml.formulaToL3String", return_value="S1"):
            topo = Topology._from_sbml_qual(model, qual_ext)
        # A non-default term is emitted with the `if` prefix.  If that prefix
        # is corrupted the generated source does not compile and the code
        # silently falls back to the interpreted drift, which works but is not
        # the compiled artefact this path promises.  The compiled function is
        # the one that carries the model-name stamp.
        self.assertTrue(hasattr(topo.drift_fn, "_drift_model_name"))

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
                    rates = topo.drift_fn(np.array([0.5]), None)
                    self.assertTrue(rr_inst.reset.called)
                    # The drift function's contract is to RETURN the rates, not
                    # merely to reset the RoadRunner instance on the way.  The
                    # reset-only assertion above let a `return None` mutant
                    # survive, so pin the value too.
                    self.assertIs(rates,
                                  rr_inst.getRatesOfChange.return_value)

if __name__ == "__main__":
    unittest.main()
