import unittest

from Muon.GUI.Common.muon_pair import MuonPair
from Muon.GUI.Common.muon_workspace import MuonWorkspace


class MuonPairTest(unittest.TestCase):
    """
    The MuonPair object encapsulates the information the describes a pair:

    - Its name
    - The name of the two groups that are combined to form the pair
    - The value of alpha
    - A workspace (Optional)

    It is intended to be lightweight, but with some type checking and basic logic, intended
    to prevent obvious misuse. The class is used by the MuonAnalysis interface and so these
    tests should be amended with care.
    """

    def test_that_cannot_initialize_MuonPair_without_name(self):
        with self.assertRaises(TypeError):
            MuonPair()

    def test_that_MuonPair_name_is_set_correctly(self):
        pair = MuonPair(pair_name="pair1")

        self.assertEqual(pair.name, "pair1")

    def test_that_cannot_set_new_name_on_pair(self):
        pair = MuonPair(pair_name="pair1")

        with self.assertRaises(AttributeError):
            pair.name = "new_name"
        self.assertEqual(pair.name, "pair1")

    def test_that_can_get_and_set_group1_name(self):
        pair = MuonPair(pair_name="pair1")

        pair.forward_group = "group1"

        self.assertEqual(pair.forward_group, "group1")

    def test_that_can_get_and_set_group2_name(self):
        pair = MuonPair(pair_name="pair1")

        pair.backward_group = "group2"

        self.assertEqual(pair.backward_group, "group2")

    def test_that_can_only_set_workspace_if_MuonWorkspace_object(self):
        pair = MuonPair(pair_name="pair1")
        self.assertIsNone(pair.workspace)

        pair.workspace = MuonWorkspace([1, 2, 3])
        self.assertIsNotNone(pair.workspace)

    def test_that_AttributeError_thrown_if_setting_workspace_to_non_MuonWorkspace_object(self):
        pair = MuonPair(pair_name="pair1")

        self.assertIsNone(pair.workspace)
        with self.assertRaises(AttributeError):
            pair.workspace = [1, 2, 3]
        self.assertIsNone(pair.workspace)

    def test_that_can_set_and_get_float_value_for_alpha(self):
        pair = MuonPair(pair_name="pair1")

        pair.alpha = 2.0

        self.assertEqual(pair.alpha, 2.0)

    def test_that_can_set_string_value_for_alpha(self):
        pair = MuonPair(pair_name="pair1")

        pair.alpha = "2.0"

        self.assertEqual(pair.alpha, 2.0)

    def test_that_cannot_add_negative_alpha(self):
        pair = MuonPair(pair_name="pair1")

        with self.assertRaises(AttributeError):
            pair.alpha = -1.0
            pair.alpha = "-1.0"


if __name__ == '__main__':
    unittest.main(buffer=False, verbosity=2)
