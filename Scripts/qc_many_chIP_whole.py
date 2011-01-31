#!/usr/bin/env python
# Time-stamp: <2009-08-11 18:45:24 Tao Liu>

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
        chrom_len[fs[0]] = int(fs[1])
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
    usage = "usage: %prog <-r rfile> [options] <wiggle files> ..."
    description = "Draw correlation plot for many wiggle files. Based on qc_chIP_whole.py"
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-d","--db",type="str",dest="dbname",help="UCSC db name for the assembly. Default: ce4",default="ce4")
    optparser.add_option("-r","--rfile",dest="rfile",
                         help="R output file. If not set, do not save R file.")
    optparser.add_option("-s","--step",dest="step",type="int",
                         help="sampling step in kbps. default: 100, minimal: 1",default=100)
    optparser.add_option("-z","--imgsize",dest="imgsize",type="int",
                         help="image size. default: 10, minimal: 10",default=10)    
    optparser.add_option("-m","--method",dest="method",type="string",default="median",
                         help="method to process the paired two sets of data in the sampling step. Choices are 'median', 'mean', and 'sample' (just take one point out of a data set). Default: median")
    
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

    # must provide >=2 wiggle files
    if wigfilenum < 2 or not options.rfile:
        optparser.print_help()
        sys.exit(1)

    if options.step < 1:
        error("Step can not be lower than 1!")
        sys.exit(1)
    if options.imgsize < 10:
        error("Image size can not be lower than 10!")
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

    if options.rfile:
        rfhd = open(options.rfile,"w")
        rfhd.write('''
library("geneplotter")  ## from BioConductor
require("RColorBrewer") ## from CRAN
''')

    # for each wig file, sample...
    for i in range(len(wigfhds)):
        wigfhd = wigfhds[i]
        wigfhd.seek(0)                  # reset
        info("read wiggle track from wiggle file #%d" % i+1)
        wig = WiggleIO.WiggleIO(wigfhd).build_binKeeper(chromLenDict=chrom_len,bin=options.step)
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
            for (parray,varray) in bkeeper.cage: # for each cage in binkeeper, there are [position array] and [value array].
                if varray:
                    ap(medfunc(varray))
                else:
                    ap(None)
        info("write values to r file")
        rfhd.write("p%d <- c(" % i )
        if p[0]:
            rfhd.write("%f" % p[0])
        else:
            rfhd.write("NA")
        for v in p[1:]:
            if v:
                rfhd.write(",%f" % v)
            else:
                rfhd.write(",NA")
        rfhd.write(")\n")
        
    rfhd.write("m <- cbind(p0")
    for i in range(wigfilenum-1):
        rfhd.write(",p%d" % (i+1))
    rfhd.write(''')
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,use="complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}
''')

        
    rfhd.write("bitmap(\"%s.bmp\",width=%d,height=%d)\n" % (options.rfile,options.imgsize,options.imgsize))
    labels = ",".join(map(lambda x:"\""+os.path.basename(x)+"\"",wigfiles))
    rfhd.write('''
pairs(m, lower.panel=function(...) {par(new=TRUE);smoothScatter(...)}, upper.panel=panel.cor, labels=c(%s))
''' % (labels))
    rfhd.write("dev.off()\n")
    rfhd.close()
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
