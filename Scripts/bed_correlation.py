#!/usr/bin/env python
# Time-stamp: <2015-01-14 13:48:42 Tao Liu>

"""Module Description

Copyright (c) 2008, 2015 Tao Liu <tliu4@buffalo.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: tliu4@buffalo.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import re
from subprocess import Popen, PIPE
#from math import sqrt
#from scipy import stats

from optparse import OptionParser

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

#hypercdf = stats.hypergeom.cdf
#std = stats.std
def check_Rscript_BEDTools ():
    try:
        p = Popen(["bedtools","coverage","-h"],stdout=PIPE,stderr=PIPE)
    except OSError as (errno,detail):
        sys.stderr.write("'bedtools coverage' can not be accessed through the command line!\n")
        sys.stderr.write("OS gave me this error message: [Errno %d]: %s\n" % (errno,detail))
        sys.stderr.write("Please download BEDtools from <http://code.google.com/p/bedtools/>.\n")
        sys.exit(1)
    try:
        p = Popen(["Rscript","--version"],stdout=PIPE,stderr=PIPE)
    except OSError as (errno,detail):
        sys.stderr.write("'Rscript' can not be accessed through the command line!\n")
        sys.stderr.write("OS gave me this error message: [Errno %d]: %s\n" % (errno,detail))
        sys.stderr.write("Please make sure you installed R(>2.9).\n")
        sys.exit(1)

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    check_Rscript_BEDTools()
    usage = "usage: %prog [options]"
    description = "Calculate bed correlation using BEDtools, based on basepair overlap comparing to random genome background. The one-side hypergeometric tests will be applied to decide whether to accept alternative hypothesis of enrichment or depletion. 'bedtools' and 'R' are both required. Note, please make sure the chromosome naming convientions of the twoe files are the same."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-i","--ifile",dest="ifile",type="string",
                         help="BED file 1")
    optparser.add_option("-j","--jfile",dest="jfile",type="string",
                         help="BED file 2")
    optparser.add_option("-g","--gsize",dest="gsize",type="float",
                         help="genome size",default=90000000)
    (options,args) = optparser.parse_args()

    if not options.ifile or not options.jfile:
        print "Need 2 BED files!"
        sys.exit(1)

    gsize = int(options.gsize)
    
    p = Popen(["bedtools coverage","-a",options.ifile,"-b",options.jfile,"-hist"],stdout=PIPE)
    q = Popen(["bedtools coverage","-b",options.ifile,"-a",options.jfile,"-hist"],stdout=PIPE)
    coverage2 = [0,0]
    for i in p.stdout:
        if i.startswith('all'):
            (tmp,depth,cover,len2,percentage)=i.strip().split()
            if int(depth)>0:
                coverage2[1] += int(cover)
                coverage2[0] =  int(len2)

    coverage1 = [0,0]
    for i in q.stdout:
        if i.startswith('all'):
            (tmp,depth,cover,len1,percentage)=i.strip().split()
            if int(depth)>0:
                coverage1[1] += int(cover)
                coverage1[0] =  int(len1)

    if coverage1[1] != coverage2[1]:
        print "Overlapped of 1 with 2 is different of 2 with 1:"
        print coverage1
        print coverage2
        print "Please run mergeBed on both of them first."
        sys.exit(1)

    #print coverage1
    #print coverage2
    common = float(coverage1[1])
    only1 = float(coverage1[0]-coverage1[1])
    only2 = float(coverage2[0]-coverage2[1])
    #print common,only1,only2

    #X = [1]*(common+only1)+[0]*only2
    #Y = [1]*(common+only2)+[0]*only1
    # r = sxy / (sx * sy)
    # sxy = ( sum(Xi*Yi) - ((sum(Xi)*sum(Yi))/n) ) / (n - 1)

    #sxy = ( common - (common+only1)*(common+only2)/(common+only1+only2) ) / (common+only1+only2-1)
    #sx = sqrt((common+only1)/(common+only1+only2)*(1- (common+only1)/(common+only1+only2)))
    #sy = sqrt((common+only2)/(common+only1+only2)*(1- (common+only2)/(common+only1+only2)))
    #print sx,sy,sxy
    #print "cor: %.4f" % (sx*sy/sxy)

    # hyper p-value for enrichment of overlap against random
    #p = hypercdf(common,gsize,gsize-common-only1,common+only2)
    r = "phyper(%d,%d,%d,%d,lower.tail=F,log.p=T)/log(10)" % (common,common+only1,gsize-common-only1,common+only2)
    p = Popen(["Rscript","-e",r],stdout=PIPE,stderr=PIPE)
    for i in p.stdout:
        if i.startswith("[1]"):
            i=i.strip()
            print "hypergeometric -log10(p-value) (upper tail CDF -- enriched against random) is: %.2f" % (abs(-1*float(i[i.find("[1]")+3:])))
            p.terminate()

    r = "phyper(%d,%d,%d,%d,lower.tail=T,log.p=T)/log(10)" % (common,common+only1,gsize-common-only1,common+only2)
    p = Popen(["Rscript","-e",r],stdout=PIPE,stderr=PIPE)
    for i in p.stdout:
        if i.startswith("[1]"):
            i=i.strip()
            print "hypergeometric -log10(p-value) (lower tail CDF -- depleted against random) is: %.2f" % (abs(-1*float(i[i.find("[1]")+3:])))
            p.terminate()


    
    #p2 = hypercdf(common,gsize,gsize-common-only2,common+only1)    
    #print "hyper pvalue of overlap comparing to genome random: %e" % p
    #print "hyper pvalue of overlap comparing to genome random: %e" % p2

    #p = hypercdf(common,common+only1+only2,only1,common+only2)
    #p2 = hypercdf(common,common+only1+only2,only2,common+only1) 
    #print "phyper(%d,%d,%d,%d,lower.tail=F,log.p=T)/log(10)" % (common,common+only1,only2,common+only2)   
    #print "hyper pvalue of overlap comparing to union regions random: %e" % p
    #print "hyper pvalue of overlap comparing to union regions random: %e" % p2        


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
