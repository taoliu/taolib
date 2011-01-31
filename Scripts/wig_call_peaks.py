#!/usr/bin/env python
# Time-stamp: <2009-10-21 17:36:38 Tao Liu>

"""Description: 

Copyright (c) 2009 Tao Liu <taoliu@jimmy.harvard.edu>

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
from taolib.CoreLib.Parser import *
from taolib.CoreLib.BasicStat.Prob import normal_cdf_inv

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
    description = "Call peaks from Wiggle file."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-i","--wfile",dest="wfile",type="string",
                         help="input wiggle file. *REQUIRED")
    optparser.add_option("-o","--bfile",dest="bfile",type="string",
                         help="output bed file. *REQUIRED")
    optparser.add_option("-c","--cutoff",dest="cutoff",type="float",
                         help="pvalue cutoff. default: 1e-5",default=1e-5)
    optparser.add_option("-l","--min-length",dest="minlen",type="int",
                         help="minimum length of peak, default: 300",default=300)
    optparser.add_option("-g","--maxgap",dest="maxgap",type="int",
                         help="maximum gap between significant points in a peak, default: 50",default=50)
    optparser.add_option("--skip-normal",dest="normalize",action="store_false",
                         help="Skip normalization the data before calling peaks",default=True)
    optparser.add_option("--skip-null",dest="nullmodel",action="store_false",
                         help="Skip sampling the data to use null (background) data to calculate mean and std",default=True)
    optparser.add_option("--sample-percent",dest="samplepercent",type="int",default=10,
                         help="Percentage to sample the data to build null model. Default: 10")
    
#     optparser.add_option("-m","--model",dest="model",type="string",
#                          help="statistic model for the score distribution. Choices are"+
#                          "'norm','binomal',and 'poisson', the default is 'norm'.", default="norm"
#                          )
    (options,args) = optparser.parse_args()

    if not options.wfile or not options.bfile or not options.cutoff:
        optparser.print_help()
        sys.exit()
    
#     model = options.model.lower()
#     if model == "norm":
#         pass
#     elif model == "binomal":
#         pass
#     elif model == "poisson":
#         pass
#     else:
#         error("Unrecognized model: %s" % (model))
#         sys.exit(1)

    f = options.wfile
    if not os.path.isfile(f):
        error("%s is not valid!" % f)
        sys.exit(1)

    try:
        fhd = open(f)
    except:
        error("Can't read %s" % f)
        sys.exit(1)

    try:
        bfhd = open(options.bfile,"w")
    except:
        error("Can't open %s to write" % options.bfile)
        sys.exit(1)

    info("open wiggle file...")
    wio = WiggleIO.WiggleIO(fhd)
    info("construct wiggle track object...")
    wtrack = wio.build_wigtrack()
    if options.normalize:
        info("Normalizing...")
        (sum_v,max_v,min_v,mean_v,std_v) = wtrack.normalize(null=options.nullmodel,sample_percent=options.samplepercent)
        if options.nullmodel:
            info("before normalization (NULL model):")
            info("mean of data: %.2f" % mean_v)
            info("std of data: %.2f" % std_v)
            info("after normalization (NULL model):")
        else:
            info("before normalization:")
            info("mean of data: %.2f" % mean_v)
            info("std of data: %.2f" % std_v)
            info("after normalization")
    if options.nullmodel:
        (sum_v,max_v,min_v,mean_v,std_v) = wtrack.null_model_summary(sample=options.samplepercent)
    else:
        (sum_v,max_v,min_v,mean_v,std_v) = wtrack.summary()
    info("mean of data: %.2f" % mean_v)
    info("std of data: %.2f" % std_v)
        
    info("Calculate score cutoff using pvalue cutoff...")
    scorecutoff = normal_cdf_inv(options.cutoff,mu=mean_v,sigma2=std_v,lower=False)
    info("scorecutoff is %f" % scorecutoff)
    info("call peaks...")
    wpeaks = wtrack.call_peaks(cutoff=scorecutoff,min_length=options.minlen,max_gap=options.maxgap)
    info("write to bed file...")
    bfhd.write(wpeaks.tobed())
    bfhd.close()
    info("finished")

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
