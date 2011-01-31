#!/usr/bin/env python

from taolib.CoreLib.Algorithm.BWT import *
import unittest

class TestBWT(unittest.TestCase):
    def setUp(self):
        self.s = "test is terrible!"
        self.b = "!tseilttr breies \x00"

    def test_bwt(self):
        b = bwt(self.s)
        #print b
        self.assertEqual(b,self.b)

    def test_ibwt(self):
        s = ibwt(self.b)
        #print s
        self.assertEqual(s,self.s)

if __name__ == '__main__':
    unittest.main()
