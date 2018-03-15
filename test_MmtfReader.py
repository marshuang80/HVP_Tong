#!/usr/bin/env python
'''
docstring
'''

import unittest
from anchors_generator import parse_j_genes


class Test(unittest.TestCase):

    def setUp(self):
        self.a = 1

    def test_a(self):
        self.assertTrue(self.a == 1)

    def tearDown(self):
        print("done")


if __name__ == '__main__':
    unittest.main()
