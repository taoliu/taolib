# Time-stamp: <2008-10-01 01:16:31 Tao Liu>

"""Module Description: BinKeeper for Wiggle-like tracks.

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import sys
import re
from bisect import insort,bisect_left,bisect_right
from array import array
import sqlite3
# ------------------------------------
# constants
# ------------------------------------
# to determine the byte size
if array('H',[1]).itemsize == 2:
    BYTE2 = 'H'
else:
    raise Exception("BYTE2 type cannot be determined!")

if array('I',[1]).itemsize == 4:
    BYTE4 = 'I'
elif array('L',[1]).itemsize == 4:
    BYTE4 = 'L'
else:
    raise Exception("BYTE4 type cannot be determined!")

if array('f',[1]).itemsize == 4:
    FBYTE4 = 'f'
elif array('d',[1]).itemsize == 4:
    FBYTE4 = 'd'
else:
    raise Exception("BYTE4 type cannot be determined!")

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------
class DBBinKeeperI:
    """DBBinKeeper is a binKeeper way to store data in a SQL DB. A
    DBBinKeeper is only for a single chromosome.
    
    """
    def __init__ (self,file=None,chromosome="chrNA",chromosomesize=1e9,bin=8):
        """initialize the connection.
        
        Note, only 'file' parameter is available. Once it is set, all the
        data will be saved in that file.

        Parameters:
        bin : size of bin in Kilo Basepair
        chromosomesize : size of chromosome, default is 1G

        If file is ':memory:', all data will be saved in memory.
        """
        self.conn = sqlite3.connect(file)
        self.chromosome = chromosome
        self.chromosomesize = chromosomesize
        self.binsize = bin*1024
        self.binnumber = int(chromosomesize/self.binsize)+1
        self.c = self.conn.cursor()

    def init_tables (self):
        #c = self.conn.cursor()
        for i in xrange(1,self.binnumber+1):
            self.c.execute("create table B%d (pos long, score float)" % i)
        self.conn.commit()
        return True

    def load_binkeeper (self, bk):
        """Load a CoreLib/BinKeeper object.

        Run this function will delete all original data.
        """
        #c = self.conn.cursor()        
        for i in xrange(1,self.binnumber+1):
            self.c.execute("drop table B%d" % i)
            self.c.execute("create table B%d (pos long, score float)" % i)
            for (p,v) in bk.cage[i]:
                self.c.execute("insert into B%d values (%d,%.6f)" % (i,p,v))
        self.conn.commit()
        return True

    def p2table (self, p ):
        """Return the table number for a position.
        
        """
        return int(p/self.binsize)+1

    def add (self, p, value):
        """self.conn.commit() needed!
        """
        #c = self.conn.cursor()        
        table = int(p/self.binsize)+1
        self.c.execute("insert into B%d values (%d,%.6f)" % (table,p,value))

    def __pp2list ( self, p1, p2):
        assert p1<=p2
        #c = self.conn.cursor()
        table1 = self.p2table(p1)
        table2 = self.p2table(p2)        
        t = [array(BYTE4,[]),array(FBYTE4,[])]
        t0a = t[0].append
        t1a = t[1].append
        for i in xrange(table1,table2+1):
            for (p,v) in self.c.execute("select * from B%d" % (i)).fetchall():
                t0a(p)
                t1a(v)
        return t

    def pp2p (self, p1, p2):
        """Give the position list between two given positions.

        Parameters:
        p1 : start position
        p2 : end position
        Return Value:
        list of positions between p1 and p2.
        """
        (ps,vs) = self.__pp2list(p1,p2)
        p1_in_list = bisect_left(ps,p1)
        p2_in_list = bisect_right(ps,p2)
        return ps[p1_in_list:p2_in_list]

    def pp2v (self, p1, p2):
        """Give the value list between two given positions.

        Parameters:
        p1 : start position
        p2 : end position
        Return Value:
        list of values whose positions are between p1 and p2.
        """
        table1 = self.p2table(p1)
        table2 = self.p2table(p2)        
        d = []
        da=d.append
        #t = [array(BYTE4,[]),array(FBYTE4,[])]
        #t0a = t[0].append
        #t1a = t[1].append
        tabnames = map(lambda x:"B"+str(x),range(table1,table2+1))
        #for i in xrange(table1,table2+1):
        #print tabname
        #if tabname.count(",")>2:
        #    raise Exception("huh")
        for i in self.c.executemany("select * from %s",tabnames).fetchall():
            for n in xrange(0,len(i),2):
                p=i[n]
                v=i[n+1]
                if p1<=p and p<=p2:
                    da(v)
                

        #(ps,vs) = self.__pp2list(p1,p2)
        #p1_in_list = bisect_left(ps,p1)
        #p2_in_list = bisect_right(ps,p2)
        #return vs[p1_in_list:p2_in_list]
        return d

    def pp2pv (self, p1, p2):
        """Give the (position,value) list between two given positions.

        Parameters:
        p1 : start position
        p2 : end position
        Return Value:
        list of (position,value) between p1 and p2.
        """
        (ps,vs) = self.__pp2list(p1,p2)
        p1_in_list = bisect_left(ps,p1)
        p2_in_list = bisect_right(ps,p2)
        return zip(ps[p1_in_list:p2_in_list],vs[p1_in_list:p2_in_list])
        
