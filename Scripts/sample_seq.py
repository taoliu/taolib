#!/usr/bin/env python
# Time-stamp: <2008-02-03 00:39:36 Tao Liu>

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
from Bio import SeqIO
from Bio.Seq import reverse_complement as rc
from random import randint
from glob import glob
from optparse import OptionParser

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

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
    optparser.add_option("-i","--ifile",dest="ifile",type="string",
                         help="input files, if you give a pattern for files, please use \" to surround the pattern string")
    optparser.add_option("-w","--window",dest="window",type="int",
                         default=1000,help="window to extract a sample sequence, default:1000")
    optparser.add_option("-s","--sample",dest="sample",type="int",
                         default=36,help="size of a sample sequence,default=36")
    optparser.add_option("-f","--frag",dest="frag",type="int",
                         help="if pair-end calculation is needed, give the fragment size")
    optparser.add_option("-o","--ofile",dest="ofile",
                         help="output file") 
    (options,args) = optparser.parse_args()
    if not options.ifile or not options.window or not options.sample or not options.ofile:
        optparser.print_help()
        sys.exit(1)
    

    ohd = open(options.ofile,"w")
    
    files = glob(options.ifile)
    if not files:
        sys.stderr.write("no file found: %s\n" % (options.ifile))
        sys.exit(1)

    m = 0
    for f in files:
        sys.stdout.write("%s ... \n" % (f))
        sys.stdout.flush()
        fahd = open(f,"r")
        record = SeqIO.parse(fahd,"fasta")
        for i in record:
            n = len(i.seq)/1000
            for l in range(0,n-1):
                s = l*1000+randint(0,999)
                if options.frag:
                    slice = i.seq[s:s+options.frag].tostring().upper()
                    if slice.find("N") == -1:
                        m+=1
                        ohd.write("%s\t%d\n" % (slice[:options.sample],m))
                        m+=1
                        ohd.write("%s\t%d\n" % (rc(slice[-1*options.sample:]),m))
                        #ohd.write("> slice%d_left\n%s\n" % (m,slice[:options.sample]))
                        #ohd.write("> slice%d_right\n%s\n" % (m,rc(slice[-1*options.sample:])))
                    else:
                        slice = i.seq[s:s+options.sample].tostring().upper()
                        if slice.find("N") == -1:
                            m+=1
                            ohd.write("%s\t%d\n" % (slice,m))
                            #ohd.write("> slice%d\n%s\n" % (m,slice))
        fahd.close()
    ohd.close()
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        msgl(_("\n;-) See you!"))
        sys.exit(0)
