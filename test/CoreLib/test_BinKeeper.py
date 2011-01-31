#!/usr/bin/env python

import os
from Cistrome.CoreLib.BinKeeper import *
from Cistrome.CoreLib.Parser import WiggleIO
from Cistrome.CoreLib.Parser import BedIO
from time import time
import unittest

class TestBinKeeper(unittest.TestCase):
    def setUp(self):
        pass
    
#     def test_load (self):
#         print "\ntest BinKeeper..."        
#         t0=time()        
#         bed = BedIO.parse_BED(file("sample.bed"))
#         print "load bed time:",time()-t0,"s for ",bed.total()," regions"
# #         t0=time()
#         w = WiggleIO.WiggleIO("sample.wig")
# #         wigtrack = w.build_wigtrack()
# #         print "load wig time:",time()-t0,"s for ",wigtrack.total," points"
#         t0=time()
#         binkeeper = w.build_binKeeper()
#         print "load binkeeper time:",time()-t0,"s for the same number of points"
# #         t0=time()
# #         a= bed.extract_wiggle_values(wigtrack)
# #         print "extract wig from bed time:",time()-t0,"s"
# #         print "total %d, first 3 result: %.3f %.3f %.3f" % (len(a),a[0],a[1],a[2]),
# #         print "and last 3 result: %.3f %.3f %.3f" % (a[-3],a[-2],a[-1])        
#         t0=time()
#         b= bed.extract_binkeepers(binkeeper)
#         print "extract binkeeper from bed time:",time()-t0,"s"        
#         print "total %d, first 3 result: %.3f %.3f %.3f" % (len(b),b[0],b[1],b[2]),
#         print "and last 3 result: %.3f %.3f %.3f" % (b[-3],b[-2],b[-1])        

    def test_dbbinkeeper (self):
        print "\ntest DBBinKeeper..."
        t0=time()        
        bed = BedIO.parse_BED(file("sample.bed"))
        print "load bed time:",time()-t0,"s for ",bed.total()," regions"
        t0=time()
#         print "build DB binkeeper..."
#         w = WiggleIO.WiggleIO("sample.wig")
#         if os.path.exists("test"):
#             os.rmdir("test")

#         binkeeper = w.build_DBBinKeeper("test",templatedb="template.db")
#         print "build DB binkeeper time:",time()-t0,"s for the same number of points"
        t0=time()
        b= bed.extract_DBBinKeeperI("test")
        print "extract DBBinKeeper from bed time:",time()-t0,"s"
        print "total %d, first 3 result: %.3f %.3f %.3f" % (len(b),b[0],b[1],b[2]),
        print "and last 3 result: %.3f %.3f %.3f" % (b[-3],b[-2],b[-1])               

if __name__ == '__main__':
    unittest.main()
