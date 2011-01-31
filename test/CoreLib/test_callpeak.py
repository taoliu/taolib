#!/usr/bin/env python

from taolib.CoreLib.Parser import *
import unittest

class TestCallPeak(unittest.TestCase):
    def setUp(self):
        self.fhd = open("chr22.score.wig")

    def test_callpeak (self):
        wio = WiggleIO.WiggleIO(self.fhd)
        wtrack = wio.build_wigtrack()
        wpeaks = wtrack.call_peaks(cutoff=10,min_window=300,max_gap=50)
        print wpeaks.tobed()


if __name__ == '__main__':
    unittest.main()
