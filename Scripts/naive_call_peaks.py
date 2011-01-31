#!/usr/bin/env python
# Time-stamp: <2010-04-28 22:59:47 Tao Liu>

"""Description: 

Copyright (c) 2009 Tao Liu <taoliu@jimmy.harvard.edu>

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
from optparse import OptionParser
from taolib.CoreLib.Parser import *
from taolib.CoreLib.BasicStat.Prob import normal_cdf_inv

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
    description = "Call peaks from Wiggle file."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-i","--wfile",dest="wfile",type="string",
                         help="input wiggle file. *REQUIRED")
    optparser.add_option("-o","--bfile",dest="bfile",type="string",
                         help="output bed file. *REQUIRED")
    optparser.add_option("-c","--cutoff",dest="cutoff",type="float",
                         help="pvalue cutoff. default: 1e-5",default=1e-5)
    optparser.add_option("-l","--min-length",dest="minlen",type="int",
                         help="minimum length of peak, default: 300",default=300)
    optparser.add_option("-g","--maxgap",dest="maxgap",type="int",
                         help="maximum gap between significant points in a peak, default: 50",default=50)
    (options,args) = optparser.parse_args()

    if not options.wfile or not options.bfile or not options.cutoff:
        optparser.print_help()
        sys.exit()
    
    f = options.wfile
    if not os.path.isfile(f):
        error("%s is not valid!" % f)
        sys.exit(1)

    try:
        fhd = open(f)
    except:
        error("Can't read %s" % f)
        sys.exit(1)

    try:
        bfhd = open(options.bfile,"w")
    except:
        error("Can't open %s to write" % options.bfile)
        sys.exit(1)

    info("open wiggle file...")
    wio = WiggleIO.WiggleIO(fhd)
    info("construct wiggle track object...")
    wtrack = wio.build_wigtrack()

    scorecutoff = float(options.cutoff)
    info("scorecutoff is %f" % scorecutoff)
    info("call peaks...")
    wpeaks = wtrack.call_peaks(cutoff=scorecutoff,min_length=options.minlen,max_gap=options.maxgap)
    info("write to bed file...")
    bfhd.write(wpeaks.tobed())
    bfhd.close()
    info("finished")

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
