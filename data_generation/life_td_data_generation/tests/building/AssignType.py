import unittest
import numpy as np
from astropy.table import Table

# Assume assign_type is imported from the module
from building import assign_type


class TestAssignTypeFunctionAstropyTable(unittest.TestCase):
    def test_assign_type_masked_constant(self):
        """Test when the value in 'type_2' is a MaskedConstant."""
        cat = Table({
            'type_1': ['TypeA', 'TypeB'],
            'type_2': [np.ma.masked, 'TypeC'],
            'type': [None, None],
        })
        result = assign_type(cat, 0)
        self.assertEqual(result, 'TypeA')
        self.assertEqual(cat['type'][0], 'TypeA')

    def test_assign_type_none_string(self):
        """Test when the value in 'type_2' is the string 'None'."""
        cat = Table({
            'type_1': ['TypeA', 'TypeB'],
            'type_2': ['None', 'TypeC'],
            'type': [None, None],
        })
        result = assign_type(cat, 0)
        self.assertEqual(result, 'TypeA')
        self.assertEqual(cat['type'][0], 'TypeA')

    def test_assign_type_valid_type_2(self):
        """Test when the value in 'type_2' is valid and not 'None' or MaskedConstant."""
        cat = Table({
            'type_1': ['TypeA', 'TypeB'],
            'type_2': ['TypeD', 'TypeC'],
            'type': [None, None],
        })
        result = assign_type(cat, 0)
        self.assertEqual(result, 'TypeD')
        self.assertEqual(cat['type'][0], 'TypeD')

    def test_assign_type_edge_case_empty(self):
        """Test an edge case with an empty table."""
        cat = Table({
            'type_1': [],
            'type_2': [],
            'type': [],
        })
        with self.assertRaises(IndexError):
            assign_type(cat, 0)


if __name__ == '__main__':
    unittest.main()
