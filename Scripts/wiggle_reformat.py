#!/usr/bin/env python
# Time-stamp: <2009-06-02 14:00:00 Tao Liu>

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
from taolib.CoreLib.Parser import WiggleIO

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
    description = "Sort and format wiggle file."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-i","--input",dest="input",type="string",
                         help="wiggle file for input")
    optparser.add_option("-o","--output",dest="output",
                         help="wiggle output file.")
    (options,args) = optparser.parse_args()
    if not options.input or not options.output:
        optparser.print_help()
        sys.exit(1)

    info("#1 read wiggle track from first wiggle file")
    wiginput = WiggleIO.WiggleIO(options.input).build_wigtrack()
    info("#1 finish reading wiggle file")
    info("#2 sorting")
    wiginput.sort()
    info("#3 write wiggle file")
    ofhd = open(options.output,"w")
    wiginput.write_wig(ofhd,name="reformatted wiggle file")
    info("#4 finished!")

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
