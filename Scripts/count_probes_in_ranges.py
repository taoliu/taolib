#!/usr/bin/env python
# Time-stamp: <2008-11-04 15:37:29 Tao Liu>

"""Module Description

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

import os
import sys
import re
import logging
from random import sample
from optparse import OptionParser
from Cistrome.CoreLib.Parser import XLSIO,WiggleIO
from Cistrome.CoreLib.BasicStat.Func import * 
# ------------------------------------
# constants
# ------------------------------------

logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )
# ------------------------------------
# Misc functions
# ------------------------------------

error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info


# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog [options]"
    description = "Summarize ChIP-chip experiment. Calculate how many probes fall into certain score ranges '(-inf,-1],(-1,0],(0,+1],(+1,+inf)'."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-x","--wig1",dest="wig1",type="string",
                         help="wiggle file for #1 replicate")
    optparser.add_option("-o","--ofile",dest="ofile",
                         help="output file")
    (options,args) = optparser.parse_args()
    if not options.wig1 or not options.ofile:
        optparser.print_help()
        sys.exit(1)

    ofhd = open(options.ofile,"w")

    info("#1 read wiggle track from peak file")
    wigtrack1 = WiggleIO.WiggleIO(options.wig1).build_wigtrack()
    info("#1 finish reading wiggle files")

    info("#2 count probes in peaks")
    chroms = wigtrack1.get_chr_names()
    ofhd.write("chr\t(-inf,-1]\tpercentage\t(-1,0]\tpercentage\t(0,+1]\tpercentage\t(+1,+inf)\tpercentage\n")
    for chrom in chroms:
        wchr = wigtrack1.get_data_by_chr(chrom)[1]
        total = len(wchr)
        # -inf .. -1, -1 .. 0 0.. +1 +1 .. +inf
        a0=0
        a1=0
        a2=0
        a3=0
        for i in wchr:
            if i <= -1:
                a0+=1
            elif i<=0:
                a1+=1
            elif i<=1:
                a2+=1
            else:
                a3+=1
        
        p0=float(a0)/total*100
        p1=float(a1)/total*100
        p2=float(a2)/total*100
        p3=float(a3)/total*100    
        
        info(" chromosome %s: total: %d, (,-1]: %.2f, (-1,0]: %d, (0,1]: %.2f, (1,): %.2f" % (chrom,total,p0,p1,p2,p3))
        ofhd.write("%s\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n" % (chrom,a0,p0,a1,p1,a2,p2,a3,p3))

    ofhd.close()
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
