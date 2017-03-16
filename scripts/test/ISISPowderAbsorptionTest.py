from __future__ import (absolute_import, division, print_function)

import mantid.simpleapi as mantid
import unittest
from six import iterkeys

from isis_powder.routines import absorb_corrections


class ISISPowderAbsorptionTest(unittest.TestCase):

    def test_sample_is_set_correctly(self):
        sample_properties = {
            "cylinder_sample_height": 4.0,
            "cylinder_sample_radius": 0.25,
            "cylinder_position": [0., 0., 0.],
            "chemical_formula": "V"
        }

        ws = mantid.CreateSampleWorkspace(Function='Flat background', NumBanks=1, BankPixelWidth=1, XMax=10, BinWidth=1)
        ws = absorb_corrections.run_cylinder_absorb_corrections(ws_to_correct=ws, multiple_scattering=False,
                                                                config_dict=sample_properties)

        self.assertAlmostEqual(ws.dataY(0)[2], 1.16864808, delta=1e-8)
        self.assertAlmostEqual(ws.dataY(0)[5], 1.16872761, delta=1e-8)
        self.assertAlmostEqual(ws.dataY(0)[9], 1.16883365, delta=1e-8)

    def test_missing_property_is_detected(self):
        sample_properties = {
            "cylinder_sample_height": 4.0,
            "cylinder_sample_radius": 0.25,
            "cylinder_position": [0., 0., 0.],
            "chemical_formula": "V"
        }

        ws = mantid.CreateSampleWorkspace(Function='Flat background', NumBanks=1, BankPixelWidth=1, XMax=2, BinWidth=1)

        # Test each key one at a time
        for blacklisted_key in iterkeys(sample_properties):
            # Force python to make a shallow copy
            modified_dict = sample_properties.copy()
            modified_dict.pop(blacklisted_key)

            # Check that is raises an error
            with self.assertRaisesRegexp(KeyError, "The following key was not found in the advanced configuration"):
                ws = absorb_corrections.run_cylinder_absorb_corrections(ws_to_correct=ws, multiple_scattering=False,
                                                                        config_dict=modified_dict)

            # Then check the error actually has the key name in it
            with self.assertRaisesRegexp(KeyError, blacklisted_key):
                ws = absorb_corrections.run_cylinder_absorb_corrections(ws_to_correct=ws, multiple_scattering=False,
                                                                        config_dict=modified_dict)

if __name__ == "__main__":
    unittest.main()
