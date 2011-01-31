#!/usr/bin/env python
# Time-stamp: <2008-02-04 13:20:05 Tao Liu>

"""Module Description

Copyright (c) 2007 Tao Liu <taoliu@jimmy.harvard.edu>

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

import os
import sys
import re
from optparse import OptionParser

# ------------------------------------
# constants
# ------------------------------------
MIN_DIST = 50
MAX_DIST = 500
# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------
class Pos:
    def __init__ (self,chr,s):
        self.chr = chr
        self.start =s

class Pair:
    def __init__ (self,n):
        self.n = n
        self.left = []
        self.right = []

    def addleft (self,pos):
        self.left.append(pos)

    def addright (self,pos):
        self.right.append(pos)

    def pair(self):
        n = 0
        for rp in self.right:
            for lp in self.left:
                if lp.chr == rp.chr:
                    dist = rp.start - lp.start
                    if dist<MAX_DIST and dist>MIN_DIST:
                        n+=1
        return n

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog [options]"
    description = "Analyze the mapping result from xMAN, report the ratio for unique mapping\nIt's the step #3 after sample_seq.py (#1) and xMAN(#2) of the whole pipeline"
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-i","--ifile",dest="ifile",type="string",
                         help="input xMAN mapping file")
    optparser.add_option("-o","--ofile",dest="ofile",type="string",
                         help="output file")
    optparser.add_option("-p","--pair",dest="pair",action="store_true",
                         help="Whether or not to parse the pair-end mapping result") 
    (options,args) = optparser.parse_args()
    
    # ... you have alarge list of positions
    if not options.ifile or not optiions.ofile:
        optparser.print_help()
        sys.exit(1)
    ifhd = open(options.ifile,"r")
    ofhd = open(options.ofile,"w")

    col_tagno = 4
    col_chr = 2
    col_start = 3

    pairs = {}

    for l in ifhd.readlines():
        if l.startswith("#"): continue
        fields = l.split("\t")
        #print fields
        chr = fields[col_chr]
        start = int(fields[col_start])
        tagno = int(fields[col_tagno])
        right = False
        if tagno % 2 ==0:
            tagno-=1
            right = True
        if not pairs.has_key(tagno):
            pairs[tagno]=Pair(tagno)
        if chr == "Nomatch":
            continue
        if right:
            pairs[tagno].addright(Pos(chr,start))
        else:
            pairs[tagno].addleft(Pos(chr,start))

    ns = pairs.keys()
    ns.sort()
    total_unique_pairs = 0
    total_pairs = len(ns)
    for n in ns:
        p = pairs[n].pair()
        ofhd.write( "%d\t%d\n" % (n,p))
        if p == 1:
            total_unique_pairs += 1
    ofhd.write( "total: %d\nmapped: %d\nratio: %.2f%%\n" % (total_pairs,total_unique_pairs,float(total_unique_pairs)/total_pairs*100) )
    ofhd.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        msgl(_("\n;-) See you!"))
        sys.exit(0)
