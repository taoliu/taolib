#!/usr/bin/env python

import os
from Cistrome.CoreLib.DB import *
import unittest
import tempfile

class TestDBRawChIP(unittest.TestCase):
    def setUp(self):
        self.result_id1 = (1, u'Oct4, Sox2, Klf4, c-Myc,Nanog, Dax1, Rex1, Zpf281, and Nac1', u'mm7', u'mouse ES cells', u'biotin mediated ChIP-chip', u'Affy Mouse Promoter 1.0R', u'Stuart H Orkin lab', 3, u'/cluster/homes/zcoban/1/CEL', u'/cluster/homes/zcoban/1/Factorname_TAG/Factorname_Cell.bed')
        self.select_result = (u'hg18', u'FoxA1')
        self.processedfileno15 = (u'/cluster/homes/zcoban/ChIP/Affy/FoxA1/TAG/FoxA1_MCF7.bed')
        self.samplendd = 'CSM11\nCTCF\nHomo sapiens\nhg18\nHeLa\nAffymetrix GeneChip Human Tiling 1.0R Array set\ncompleted\n/cluster/homes/zcoban/11/CEL/CTCF_raw.tar.gz\n/cluster/homes/zcoban/11/TAG/CTCF_peak.xls\n/cluster/homes/zcoban/11/TAG/CTCF_wiggle.tar.gz\nLaboratory of Chromosome Structure and Inheritance, Graduate School of Bioscience and Biotechnology, Tokyo Institute of Technology, Department of Biological Sciences\nCohesin mediates transcriptional insulation by CCCTC-binding factor\nG2 phase\nGPL6131:GSE9613:GSM243234,GSM243235,GSM243236,GSM243237,GSM243238,GSM243239,GSM243240,GSM243241,GSM243242,GSM243243,GSM243244,GSM243245,GSM243246,GSM243247\n18235444\nChIP-chip\nUsing MAT, default, p-value cutoff 1e-3\nZcoban\n'

    def test_create (self):
        db = dbRawChIP.DBRawChIP(file=':memory:')
        db.create_table()

    def test_add (self):
        db = dbRawChIP.DBRawChIP(file=':memory:')
        db.create_table()
        db.add_chip(1, u'Oct4, Sox2, Klf4, c-Myc,Nanog, Dax1, Rex1, Zpf281, and Nac1', u'mm7', u'mouse ES cells', u'biotin mediated ChIP-chip', u'Affy Mouse Promoter 1.0R', u'Stuart H Orkin lab', 3, u'/cluster/homes/zcoban/1/CEL', u'/cluster/homes/zcoban/1/Factorname_TAG/Factorname_Cell.bed')
        db.add_chip(15,'FoxA1','hg18','MCF7','','Affy WG','Brown lab','23456','/cluster/homes/zcoban/Affywg_p63_ActD_Molcel06_Struhl/CEL/GSE6132_RAW','/cluster/homes/zcoban/3/TAG/p63_Act1.bed')
        c = db.execute('select * from chips where id == 1')
        for row in c:
            self.assertEqual(row[:-1],self.result_id1)
        c = db.execute('select organism,factor from chips where id == 15')
        for row in c:
            self.assertEqual(row,self.select_result)

    def test_add_table (self):
        db = dbRawChIP.DBRawChIP(file=':memory:')
        db.create_table()
        db.add_from_table (open('datacollection.csv','r'))
        c = db.execute('select * from chips where id == 1')
        for row in c:
            self.assertEqual(row[:-1],self.result_id1)
        c = db.execute('select organism,factor from chips where id == 15')
        for row in c:
            self.assertEqual(row,self.select_result)
        
#     def test_id2processedpath (self):
#         db = dbRawChIP.DBRawChIP(file=':memory:')
#         db.create_table()
#         db.add_from_table (open('datacollection.csv','r'))
#         self.assertEqual(db.id2processedpath(tfile,15),self.processedfileno15)

    def test_NDDparse (self):
        d = NDD.parse(file('sample.ndd','r'))
        self.assertEqual(str(d),self.samplendd)

    def test_dbRawData (self):
        if os.path.exists('tmpdbrd'):
            os.remove('tmpdbrd')
        db = dbRawChIP.DBRawData(file='tmpdbrd')
        db.init_tables()
        db.add_NDDfile('sample.ndd')
        c = db.execute('select * from dataset where Factorid like "%CTCF%"')
        for i in c:
            print i
        
if __name__ == '__main__':
    unittest.main()
