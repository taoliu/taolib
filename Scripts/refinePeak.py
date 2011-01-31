#!/usr/bin/env python
# Time-stamp: <2009-01-08 13:21:16 Tao Liu>

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
from optparse import OptionParser

from Cistrome.CoreLib.Parser import WiggleIO
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
    description = "Refine ChIP-chip peak sub-summit, boundaries and find classes of binding events."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-i","--ifile",dest="ifile",type="string",
                         help="input wig file from ChIP-chip")
    optparser.add_option("-b","--bw",dest="bw",type="int",
                         help="band width, default: 100",default=100)
    optparser.add_option("--peak",dest="plow",type="int",
                         help="lowest score for peak, default: 1",default=1)
    optparser.add_option("--valley",dest="vup",type="int",
                         help="highest score for valley, default: 0.2",default=0.2)
    optparser.add_option("-o","--ofile",dest="ofile",
                         help="output wig file") 
    (options,args) = optparser.parse_args()
    if not options.ifile or not options.ofile:
        optparser.print_help()
        sys.exit(1)
    else:
        iwig = options.ifile
        owig = options.ofile

    info("read wiggle file")
    wigtrack = WiggleIO.WiggleIO(iwig).build_wigtrack()
    info("find all local summits")
    peaks = wigtrack.find_peaks(bw=options.bw)
    info("filter peaks using cutoff")
    peaks = peaks.filter_score(options.plow)
    info("find all local valleys")
    valleys = wigtrack.find_valleys(bw=options.bw)
    info("filter valleys using cutoff")
    valleys = valleys.filter_score_below(options.vup)

    info("write result")
    peaks.write_wig(file(owig,"w"),"refine")
    valleys.write_wig(file(owig+"valley","w"),"refine")    
    

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
