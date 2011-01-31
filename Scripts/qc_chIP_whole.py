#!/usr/bin/env python
# Time-stamp: <2009-07-24 16:22:56 Tao Liu>

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
    description = "QC for two replicates of ChIP-chip experiment."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-x","--wig1",dest="wig1",type="string",
                         help="wiggle file for #1 replicate")
    optparser.add_option("-y","--wig2",dest="wig2",type="string",
                         help="wiggle file for #2 replicate")
    optparser.add_option("-r","--rfile",dest="rfile",
                         help="R output file. If not set, do not save R file.")
    optparser.add_option("-s","--step",dest="step",type="int",
                         help="sampling step, two sets of the data points in the region will be paried. default: 100",default=100)
    optparser.add_option("-f","--format",dest="format",type="string",
                         help="ma2c, mat or macs, default: ma2c",default="ma2c")
    optparser.add_option("-m","--method",dest="method",type="string",default="sample",
                         help="method to process the paired two sets of data in the sampling step. Choices are 'median', 'mean','sum', and 'sample' (just take one point out of a data set). Default: sample")
    (options,args) = optparser.parse_args()
    if not options.wig1 or not options.wig2:
        optparser.print_help()
        sys.exit(1)
    method = options.method.lower()
    if method == 'median':
        medfunc = median
    elif method == 'mean':
        medfunc = mean
    elif method == 'sum':
        medfunc = sum
    elif method  == 'sample':
        medfunc = lambda u: u[-1]
    else:
        print "unrecognized method: %s" % (method)
        sys.exit(1)

    info("#1 read wiggle track from first wiggle file")
    wig1 = WiggleIO.WiggleIO(options.wig1).build_wigtrack()
    info("#1 read wiggle track from second wiggle file")    
    wig2 = WiggleIO.WiggleIO(options.wig2).build_wigtrack()
    info("#1 finish reading wiggle files")

    info("#4 whole data set level")

    p1 = []
    p2 = []
    ap1 = p1.append
    ap2 = p2.append
    for chrom in wig1.get_chr_names():
        count = 0
        sum_region1 = []
        sum_region2 = []
        (p1chr,w1chr) = wig1.get_data_by_chr(chrom)
        try:
            (p2chr,w2chr) = wig2.get_data_by_chr(chrom)
        except:
            info("no chromosome %s in %s!" % (chrom, options.wig2))
            continue
        j = 0
        l1chr = len(w1chr)
        l2chr = len(w2chr)
        for i in range(l1chr):
            count += 1
            while j < l2chr:
                if p1chr[i] == p2chr[j]:
                    if count >= options.step:
                        count = 0
                        sum_region1.append(w1chr[i])
                        sum_region2.append(w2chr[j])
                        ap1(medfunc(sum_region1))
                        ap2(medfunc(sum_region2))
                        sum_region1 = []
                        sum_region2 = []
                    else:
                        sum_region1.append(w1chr[i])
                        sum_region2.append(w2chr[j])
                    j+=1
                    break
                elif p1chr[i] < p2chr[j]:
                    break
                else:
                    j+=1
            if j >= l2chr:
                count=0
                break

    print len(p1),len(p2)
    if options.rfile:
        rfhd = open(options.rfile,"w")
        rfhd.write('''
library("geneplotter")  ## from BioConductor
require("RColorBrewer") ## from CRAN
''')
        rfhd.write("p1 <- c(%s)\n" % (",".join(map(str,p1)))   )
        rfhd.write("p2 <- c(%s)\n" % (",".join(map(str,p2)))   )

    p1 = centering(p1)
    p2 = centering(p2)
    cor = sum(map(lambda x:x[0]*x[1],zip(p1,p2)))/(len(p1)-1)/std(p1)/std(p2)

    p1name = os.path.basename(options.wig1.rsplit(".wig",2)[0])
    p2name = os.path.basename(options.wig2.rsplit(".wig",2)[0])


    if options.rfile:
        rfhd.write("bitmap(\"%s.bmp\",width=10,height=10)\n" % options.rfile)
        rfhd.write("smoothScatter(p1,p2,main=\"cor=%.2f\",xlab=\"%s\",ylab=\"%s\")\n" % (cor,p1name,p2name))
        rfhd.write("dev.off()\n")
        rfhd.close()
    
    info("#4 overall whole data set level cor: %f" % cor)
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
