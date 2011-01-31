#!/usr/bin/env python
# Time-stamp: <2011-01-31 16:56:46 Tao Liu>

"""Script description: HMM conception modeling the Histone Marks K4me3 and K36me3 in worm.

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
from optparse import OptionParser

import logging

from taolib.CoreLib.Parser import *
from taolib.CoreLib.BasicStat.Prob import *
from taolib.CoreLib.FeatIO import *

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

states=('active','inactive')
observations=([],[])            # K4 and K36
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
    description = "blah blah blah"
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-i","--k4wig",dest="k4wig",type="string",
                         help="input wig for k4me3")
    optparser.add_option("-j","--k36wig",dest="k36wig",type="string",
                         help="input wig for k36me3")
    optparser.add_option("-o","--ofile",dest="ofile",
                         help="output file") 
    (options,args) = optparser.parse_args()

    
    f_k4wig = options.k4wig
    f_k36wig = options.k36wig
    if not os.path.isfile(f_k4wig) or not os.path.isfile(f_k36wig):
        error("wiggle files are not valid!")
        sys.exit(1)

    try:
        fhd_k4 = open(f_k4wig)
        fhd_k36 = open(f_k36wig)
    except:
        error("Can't read wiggle files")
        sys.exit(1)

    info("open treatment wiggle file...")
    k4io = WiggleIO.WiggleIO(fhd_k4)
    info("construct treatment wiggle track object...")
    ttrack = tio.build_wigtrack()
    tsum = ttrack.summary()[0]
    lambda_bg = tsum/options.gsize*50
    info("background average value: %.2f" % lambda_bg)

    


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
