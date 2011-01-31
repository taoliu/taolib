#!/usr/bin/env python
# Time-stamp: <2011-01-31 17:12:01 Tao Liu>

"""Description: Convert wiggle file to a bedGraph file with fixed bins.

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
import subprocess
from taolib.CoreLib.Parser import WiggleIO
from taolib.CoreLib.BasicStat.Func import * 

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
        tmpfname = tempfile.mkstemp(prefix="qcmanychip")[1]
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
    usage = "usage: %prog [options] <wiggle files> ..."
    description = "convert wiggle file to a bedGraph file with fixed bins"
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-d","--db",type="str",dest="dbname",help="UCSC db name for the assembly. Default: ce4",default="ce4")
    optparser.add_option("-s","--step",dest="step",type="int",
                         help="sampling step in kbps. default: 200, minimal: 100",default=200)
    optparser.add_option("-m","--method",dest="method",type="string",default="median",
                         help="method to process the paired two sets of data in the sampling step. Choices are 'median', 'mean', and 'sample' (just take one point out of a data set). Default: median")
    optparser.add_option("-l","--wig-label",dest="wiglabel",type="string",action="append",
                         help="the wiggle file labels in the figure. No space is allowed. This option should be used same times as wiggle files, and please input them in the same order as -w option. default: will use the wiggle file filename as labels.")
    
    (options,wigfiles) = optparser.parse_args()

    method = options.method.lower()
    if method == 'median':
        medfunc = median
    elif method == 'mean':
        medfunc = mean
    elif method  == 'sample':
        medfunc = lambda u: u[-1]
    else:
        print "unrecognized method: %s" % (method)
        sys.exit(1)

    wigfilenum = len(wigfiles)

    # must provide >=1 wiggle files
    if wigfilenum < 1:
        optparser.print_help()
        sys.exit(1)

    # wig labels
    if options.wiglabel and len(options.wiglabel) == wigfilenum:
        wiglabel = options.wiglabel
    else:                               # or use the filename
        wiglabel = map(lambda x:os.path.basename(x),wigfiles)
        

    if options.step < 100:
        error("Step can not be lower than 100bps!")
        sys.exit(1)

    # check the files
    for f in wigfiles:
        if not os.path.isfile(f):
            error("%s is not valid!" % f)
            sys.exit(1)
        
    wigfhds = map(open,wigfiles)        # file handlers for wiggle files

    info("number of wiggle files: %d" % wigfilenum)
    # get chromosome length info from UCSC
    info("connect to UCSC to get chromosome length information")
    try:
        chrom_len = get_chrom_length(options.dbname)
    except:
        error("Error!")
        sys.exit(1)

    # get the common chromosome list:
    chromsdict = {}
    for wigfhd in wigfhds:
        for l in wigfhd:
            if l.find("chrom=") != -1:
                c = re.search("chrom=(\w+)",l).group(1)
                chromsdict[c] = chromsdict.setdefault(c,0)+1
    chroms = []
    for c in chromsdict.keys():
        if chromsdict[c]==wigfilenum:
            chroms.append(c)

    info("common chromosomes are %s..." % ",".join(chroms))
    step = options.step
    # for each wig file, sample...
    for i in range(len(wigfhds)):
        datafilename = wiglabel[i]+".data"
        datafhd = open(datafilename,"w")
        wigfhd = wigfhds[i]
        wigfhd.seek(0)                  # reset
        info("read wiggle track from wiggle file #%d" % (i+1))
        wig = WiggleIO.WiggleIO(wigfhd).build_binKeeper(chromLenDict=chrom_len,binsize=step)
        p = []                          # size of p is number of cages in binkeeper
        ap = p.append
        for chrom in chroms:
            step_region = []
            step_region_a = step_region.append
            try:
                end_of_chrom = chrom_len[chrom]
            except:
                warn("chromosome %s can not be found in UCSC... skip..." % chrom)
                continue
            bkeeper = wig[chrom] # binkeeper object
            j = 0
            for (parray,varray) in bkeeper.cage: # for each cage in binkeeper, there are [position array] and [value array].
                j = j+1
                if varray:
                    datafhd.write("%s\t%d\t%d\t%f\n" % (chrom,(j-1)*step,j*step,medfunc(varray)))
                else:
                    datafhd.write("%s\t%d\t%d\tNA\n" % (chrom,(j-1)*step,j*step))
        wigfhd.close()
        datafhd.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
