#!/usr/bin/env python
# Time-stamp: <2010-12-16 12:53:17 Tao Liu>

"""Description: Draw correlation plot for many wiggle files.

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
import urllib2
import tempfile
import gzip
import random
# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

UCSC_chrom_URL = r'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/chromInfo.txt.gz'

# ------------------------------------
# Misc functions
# ------------------------------------

error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info

def get_chrom_length ( dbname ):
    # first try to find a local file called dbname
    try:
        fhd = open(dbname,"r")
        chrom_len = {}
        for l in fhd:
            fs = l.split()
            try:
                chrom_len[fs[0]] = int(fs[1])
            except:
                pass
        fhd.close()
    except:
        # get chromosome length from UCSC db download page
        f = urllib2.urlopen(UCSC_chrom_URL % (dbname))
        tmpfname = tempfile.mkstemp(prefix="randomRegions")[1]
        tmpf = open(tmpfname,'w')
        tmpf.write(f.read())                # write file content to temp file
        tmpf.close()
        f.close
        # read it
        fhd = gzip.open(tmpfname,'r')
        chrom_len = {}
        for l in fhd:
            fs = l.split()
            try:
                chrom_len[fs[0]] = int(fs[1])
            except:
                pass
        fhd.close()
        os.unlink(tmpfname)
    return chrom_len

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog [-d chromInfo] [-w width] [-o out.bed] <number> ..."
    description = "Generate <number> random regions in the given genome."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-w",type="int",dest="wd",help="The width of each random region in bps, default: 1",default=1)    
    optparser.add_option("-d",type="str",dest="dbname",help="UCSC db name for the assembly. It can be an actual chromInfo.txt file downloaded from UCSC, or the UCSC database name like hg19 or mm9 (require internet connection). Default: hg19",default="hg19")
    optparser.add_option("-o",type="str",dest="ofile",help="Output BED file. If not set, use stdout.")
    
    (options,args) = optparser.parse_args()

    if not args:
        optparser.print_help()
        sys.exit(1)

    if not args[0].isdigit():
        error("<number> must be integer, but you provide '%s'" % args[0])
        sys.exit(1)
    number = int(args[0])

    # get chromosome length info from UCSC
    try:
        chrom_len = get_chrom_length(options.dbname)
        chrom_names = chrom_len.keys()
    except:
        error("Error to get the chromosome information!")
        sys.exit(1)
    wd = options.wd

    if options.ofile:
        ofhd = open(options.ofile,"w")
    else:
        ofhd = sys.stdout

    total_len = sum(chrom_len.values())

    # build bins for chromosome lengths
    chrom_bins = {}
    p=0
    for c in chrom_names:
        s=p
        p+=chrom_len[c]
        e=p
        chrom_bins[c]=(s,e)             # record the start and end of chrom bin
    print chrom_bins

    rr = random.randrange
    for i in range(number):
        # pick a random position in the whole genome
        p = rr(total_len)
        for c in chrom_bins:
            (s,e)=chrom_bins[c]
            # find the chromosome name
            if p>=s and p<e:
                ofhd.write("%s\t%d\t%d\n" % (c,p,p+wd))
                break
    ofhd.close()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
