#!/usr/bin/env python
# Time-stamp: <2008-11-04 15:34:42 Tao Liu>

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
    description = "Summarize ChIP-chip experiment. Calculate how many probes are included in peak regions."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-p","--peak1",dest="peak1",type="string",
                         help="peak file in xls format for #1 replicate")
    optparser.add_option("-x","--wig1",dest="wig1",type="string",
                         help="wiggle file for #1 replicate")
    optparser.add_option("-o","--ofile",dest="ofile",
                         help="output file")
    optparser.add_option("-f","--format",dest="format",type="string",
                         help="ma2c, mat or macs, default: ma2c",default="ma2c")
    (options,args) = optparser.parse_args()
    if not options.peak1 or not options.wig1 or not options.ofile:
        optparser.print_help()
        sys.exit(1)
    format = options.format.lower()
    if format == 'ma2c':
        xlsparser = XLSIO.parse_pMA2C_xls
    elif format == 'mat':
        xlsparser = XLSIO.parse_MAT_xls
    elif format == 'macs':
        xlsparser = XLSIO.parse_MACS_xls
    else:
        print "unrecognized format: %s" % (format)
        sys.exit(1)

    ofhd = open(options.ofile,"w")

    info("#1 read peaks from first replicate")
    peaks1 = XLSIO.parse_pMA2C_xls(open(options.peak1,"r"))
    info("#1 read wiggle track from peak file")
    wigtrack1 = WiggleIO.WiggleIO(options.wig1).build_wigtrack()
    info("#1 finish reading wiggle files")

    info("#2 count probes in peaks")
    counts = peaks1.extract_wiggle_values_by_chrom(wigtrack1,func=len)

    info("#3 output")

    chroms = wigtrack1.get_chr_names()
    ofhd.write("chr\ttotal\tin_peak\tpercentage\n")
    for chrom in chroms:
        total_probe_chr = len(wigtrack1.get_data_by_chr(chrom)[0])
        if counts.has_key(chrom):
            peak_probe_chr = sum(counts[chrom])
        else:
            peak_probe_chr = 0
        info(" chromosome %s: total %d, in peak %d" % (chrom,total_probe_chr,peak_probe_chr))
        ofhd.write("%s\t%d\t%d\t%.2f\n" % (chrom,total_probe_chr,
                                           peak_probe_chr,
                                           100.0*peak_probe_chr/total_probe_chr))
    ofhd.close()
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
