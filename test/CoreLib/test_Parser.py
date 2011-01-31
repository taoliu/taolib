#!/usr/bin/env python

from taolib.CoreLib.Parser import *
import unittest
import tempfile

test_wiggle_content_correct = """track type=wiggle_0 name='' description=''
browse unused options...
variableStep chrom=chrIII span=10
1   1.0
11  2.0
21  1.0
31  2.2
41  3.0
430984 43.0984
variableStep chrom=chrV span=10
9177243 1.1
1   1.0
11  2.0
21  1.0
31  2.2
41  3.0
variableStep chrom=chrV span=10
1   1.0
11  2.0
21  1.0
31  2.2
41  3.0
variableStep chrom=chrI span=10
1   1.0
11  2.0
21  1.0
31  2.2
41  3.0
variableStep chrom=chrII span=0
1   1.0
11  2.0
21  1.0
31  2.2
41  3.0
fixedStep chrom=chrX start=100 step=10 span=0
1.0
2.0
1.0
2.2
3.0
5.0
2.0
"""

class TestWiggleIO(unittest.TestCase):
   def setUp(self):
      self.fhd = open("tmp.wig","w")
      self.fhd.write(test_wiggle_content_correct)
      self.fhd.close()
      self.fhd = open("tmp.wig")

   def test_wiggleParser(self):
      wio = WiggleIO.WiggleIO(self.fhd)
      wtrack = wio.build_wigtrack()
      wtrack.sort()
      chroms = wtrack.get_chr_names()
      for chrom in chroms:
         d = wtrack.get_data_by_chr (chrom)
         print chrom,d

# class TestTagIO(unittest.TestCase):
#     def setUp(self):
#         self.lmultifhd = file('test_1_eland_multi.txt')
#         self.rmultifhd = file('test_2_eland_multi.txt')
#         self.lmultinum4EMP = 24
#         self.PEEMP_wiggle = """track type=wiggle_0 name="tag list" description=""
# variableStep chrom=chrIII.fa span=0 strand=0
# 430984
# variableStep chrom=chrV.fa span=0 strand=0
# 9177243
# variableStep chrom=chrV.fa span=0 strand=1
# 15381282
# variableStep chrom=chrI.fa span=0 strand=0
# 3865077
# 6406321
# 7695304
# variableStep chrom=chrII.fa span=0 strand=0
# 8664810
# 9680739
# 11486537
# 14009203
# 14037935
# variableStep chrom=chrX.fa span=0 strand=0
# 1598625
# 14024566
# variableStep chrom=chrX.fa span=0 strand=1
# 7135311
# """


#     def test_ELANDMultiParser(self):
#         EMP = TagIO.ELANDMultiParser()
#         fw = EMP.build_fwtrack(self.lmultifhd)
#         self.assertEqual(fw.total,self.lmultinum4EMP)

#     def test_PairEndELANDMultiParser(self):
#         PEEMP = TagIO.PairEndELANDMultiParser()
#         fw = PEEMP.build_fwtrack(self.lmultifhd,self.rmultifhd)
#         fw.sort()
#         self.assertEqual(str(fw),self.PEEMP_wiggle)

# class TestTagPeakOverlap(unittest.TestCase):
#     def setUp(self):
#         self.tagfhd = file('tags.bed')
#         self.peakfhd = file('peak.bed')
#         self.vresult = {'chr2:14...46': '\t0...25+\t10...35+\t40...65+\t20...45-\t30...55-', 'chr1:24...56': '\t0...25+\t10...35+\t40...65+\t20...45-\t30...55-\t50...75-', 'chr1:100...200': '\t80...105+\t100...125+\t90...115-'}
#         self.nresult = 3

#     def test_Overlap(self):
#         print "\nTest OverLap..."
#         tagtrack =  TagIO.BEDParser().build_fwtrack(self.tagfhd)
#         tagtrack.fw=25

#         peaktrack = BedIO.parse_BED(self.peakfhd)

#         (n,v) = peaktrack.overlap_with_FWTrackI(tagtrack)
#         self.assertEqual(n,self.nresult)
#         self.assertEqual(v,self.vresult)
#         print "Done"


if __name__ == '__main__':
    unittest.main()
