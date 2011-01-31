#!/usr/bin/env python
# Time-stamp: <2010-05-27 14:31:30 Tao Liu>

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
from taolib.CoreLib.Parser import XLSIO,WiggleIO,BedIO
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
    optparser.add_option("-p","--peak1",dest="peak1",type="string",
                         help="peak file in bed or xls format for #1 replicate")
    optparser.add_option("-q","--peak2",dest="peak2",type="string",
                         help="peak file in bed of xls format for #2 replicate")
    optparser.add_option("-x","--wig1",dest="wig1",type="string",
                         help="wiggle file for #1 replicate")
    optparser.add_option("-y","--wig2",dest="wig2",type="string",
                         help="wiggle file for #2 replicate")
    optparser.add_option("-r","--rfile",dest="rfile",
                         help="R output file")
    optparser.add_option("-f","--format",dest="format",type="string",
                         help="bed, ma2c, mat or macs, default: ma2c",default="bed")
    optparser.add_option("-s","--step",dest="step",type="int",
                         help="number of steps to calculate cor based on some score ranges, default: 5",default=5)
    optparser.add_option("-m","--method",dest="method",type="string",default="sample",
                         help="method to process the paired two sets of data in the sampling step. Choices are 'median', 'mean','sum', and 'sample' (just take one point out of a data set). Default: sample")
    (options,args) = optparser.parse_args()
    if not options.peak1 or not options.peak2 or not options.wig1 or not options.wig2 or not options.rfile:
        optparser.print_help()
        sys.exit(1)
    peakformat = options.format.lower()
    if peakformat == "bed":
        peakparser = BedIO.parse_BED
    elif peakformat == 'ma2c':
        peakparser = XLSIO.parse_pMA2C_xls
    elif peakformat == 'mat':
        peakparser = XLSIO.parse_MAT_xls
    elif peakformat == 'macs':
        peakparser = XLSIO.parse_MACS_xls
    else:
        print "unrecognized format: %s" % (peakformat)
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


#    ofhd = open(options.ofile,"w")
    rfhd = open(options.rfile,"w")

    info("#1 read peaks from first replicate")
    peaks1 = peakparser(open(options.peak1,"r"))
    info("#1 read peaks from second replicate")    
    peaks2 = peakparser(open(options.peak2,"r"))
    info("#1 read wiggle track from first wiggle file")
    bk1 = WiggleIO.WiggleIO(options.wig1).build_binKeeper()
    info("#1 read wiggle track from second wiggle file")    
    bk2 = WiggleIO.WiggleIO(options.wig2).build_binKeeper()
    info("#1 finish reading wiggle files")

    info("#2 overall correlation coefficient")
    info("#2 merge replicates")
    union =peaks1.dup()
    union.add_other_peakI(peaks2)
    union.sort()
    union.merge_overlap()
    info("#2 num of union regions: %d" % union.total())
    info("#2 extract scores in first replicate")
    pscore1 = union.extract_binkeepers_pv(bk1)
    info("#2 extract scores in second replicate")    
    pscore2 = union.extract_binkeepers_pv(bk2)
    #info("#2 num of points %d and %d" % (len(pscore1),len(pscore2)))

    union_chroms = set(pscore1.keys()+pscore2.keys())

    p1 = []
    p2 = []
    ap1 = p1.append
    ap2 = p2.append
    for chrom in union_chroms:
        count = 0
        try:
            pv1chr = reduce(lambda x,y:x+y,pscore1[chrom])
            pv2chr = reduce(lambda x,y:x+y,pscore2[chrom])
        except:
            next
        #(p1chr,w1chr) = wig1.get_data_by_chr(chrom)
        #(p2chr,w2chr) = wig2.get_data_by_chr(chrom)
        j = 0
        l1chr = len(pv1chr)
        l2chr = len(pv2chr)
        sum_region1 = []
        sum_region2 = []
        for i in range(l1chr):
            count += 1
            while j < l2chr:
                if pv1chr[i][0] == pv2chr[j][0]:
                    if count >= options.step:
                        count = 0
                        sum_region1.append(pv1chr[i][1])
                        sum_region2.append(pv2chr[j][1])
                        try:
                            ap1(medfunc(sum_region1))
                            ap2(medfunc(sum_region2))
                        except:
                            print sum_region1
                            print sum_region2
                        sum_region1 = []
                        sum_region2 = []
                    else:
                        sum_region1.append(pv1chr[i][1])
                        sum_region2.append(pv2chr[j][1])
                    j+=1
                    break
                elif pv1chr[i][0] < pv2chr[j][0]:
                    break
                else:
                    j+=1
            if j >= l2chr:
                count=0
                break

    p1name = os.path.basename(options.wig1.rsplit(".wig",2)[0])
    p2name = os.path.basename(options.wig2.rsplit(".wig",2)[0])
    rfhd.write('''library("geneplotter")  ## from BioConductor
    require("RColorBrewer") ## from CRAN
    ''')
    rfhd.write("p1 <- c(%s)\n" % (",".join(map(str,p1)))   )
    rfhd.write("p2 <- c(%s)\n" % (",".join(map(str,p2)))   )
    rfhd.write("pdf(\"%s.pdf\",width=10,height=10)\n" % options.rfile)

    info("#2 centering pscore1")
    pscore1 = centering(p1)
    info("#2 centering pscore2")    
    pscore2 = centering(p2)
    # overall correlation coefficient
    cor = sum(map(lambda x:x[0]*x[1],zip(pscore1,pscore2)))/(len(pscore1)-1)/std(pscore1)/std(pscore2)

    rfhd.write("smoothScatter(p1,p2,main=\"cor=%.2f\",xlab=\"%s\",ylab=\"%s\")\n" % (cor,p1name,p2name))
    #rfhd.write("plot(p1,p2,cex=0.5,main=\"cor=%.2f\",xlab=\"%s\",ylab=\"%s\")\n" % (cor,p1name,p2name))
    rfhd.write("dev.off()\n")
    rfhd.close()


        
def cor4region(wig1,wig2,lower,upper):
    p1 = []
    p2 = []
    ap1 = p1.append
    ap2 = p2.append
    for chrom in wig1.get_chr_names():
        (p1chr,w1chr) = wig1.get_data_by_chr(chrom)
        (p2chr,w2chr) = wig2.get_data_by_chr(chrom)
        j = 0
        for i in xrange(len(w1chr)):
            while 1:
                if p1chr[i] == p2chr[j]:
                    if ((w1chr[i] >=lower and w1chr[i] < upper) or
                        (w2chr[j] >=lower and w2chr[j] < upper)):
                        ap1(w1chr[i])
                        ap2(w2chr[j])
                        j+=1
                    break
                elif p1chr[i] < p2chr[j]:
                    break
                else:
                    j+=1
    zp12 = zip(p1,p2)
    if len(zp12)>1000:
       zp12=sample(zp12,1000)
       p1 = map(lambda x:x[0],zp12)
       p2 = map(lambda x:x[1],zp12)        
    p1 = centering(p1)
    p2 = centering(p2)
    cor = sum(map(lambda x:x[0]*x[1],zip(p1,p2)))/(len(p1)-1)/std(p1)/std(p2)
    return (cor,p1,p2)
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
