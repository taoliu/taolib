#!/usr/bin/env python
# Time-stamp: <2008-03-17 03:47:15 Tao Liu>

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
from optparse import OptionParser
from glob import glob
import Cistrome
from Cistrome.TabIO import TabFile

# ------------------------------------
# constants
# ------------------------------------
GENOME_SIZE = {"mm8":2644077689L,
               "hg18":3080419480L}

# ------------------------------------
# Misc functions
# ------------------------------------
def log( msg ):
    sys.stderr.write( msg )

def fact (n):
    return reduce(lambda a,b:a*(b+1),range(n),1)
# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog [options]"
    description = "Slide window based on BED file, then plot the distribution of hits in window along the whole genome."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-s","--species",dest="species",type="string",
                         help="species, must be \"mm8\" or \"hg18\"")    
    optparser.add_option("-i","--ifile",dest="ifile",type="string",
                         help="input BED file, e.g. the ChIP regions")
    optparser.add_option("-w","--window",dest="window",type="int",
                         default=1000,help="window length, default:1000")
    optparser.add_option("-g","--gap",dest="gap",type="int",
                         default=0,help="gap length between nearby windows, must be no less than 0, default:0")
    optparser.add_option("-o","--ofile",dest="ofile",
                         help="output file") 
    (options,args) = optparser.parse_args()

    if not options.ifile or not options.ofile:
        optparser.print_help()
        sys.exit(1)

    ofhd = open(options.ofile,"w")

    if options.species == "hg18":
        chromosomes_len = {                     # store start and number of file-pos for every chromosome in bin file
            "chr1":247249719,"chr2":242951149,"chr3":199501827,
            "chr4":191273063,"chr5":180857866,"chr6":170899992,
            "chr7":158821424,"chr8":146274826,"chr9":140273252,
            "chr10":135374737,"chr11":134452384,"chr12":132349534,
            "chr13":114142980,"chr14":106368585,"chr15":100338915,
            "chr16":88827254,"chr17":78774742,"chr18":76117153,
            "chr19":63811651,"chr20":62435964,"chr21":46944323,
            "chr22":49691432,"chrX":154913754,"chrY":57772954,
            "chrM" :16571
            }
        genome_len = 3080436051L
    elif options.species == "mm8":
        chromosomes_len = {                     # store start and number of file-pos for every chromosome in bin file
            "chr1" :197069962,"chr2" :181976762,"chr3" :159872112,
            "chr4" :155029701,"chr5" :152003063,"chr6" :149525685,
            "chr7" :145134094,"chr8" :132085098,"chr9" :124000669,
            "chr10":129959148,"chr11":121798632,"chr12":120463159,
            "chr13":120614378,"chr14":123978870,"chr15":103492577,
            "chr16":98252459, "chr17":95177420, "chr18":90736837,
            "chr19":61321190, "chrX" :165556469,"chrY" :16029404,
            "chrM" :16299
            }
        genome_len = 2644093988L
    else:
        log ("Species must be \"mm8\" for mouse or \"hg18\" for human!\n")
        sys.exit(1)


    log("#1. Read BED file: %s\n" % (options.ifile))
    f = TabFile(options.ifile)
    # add tags
    track = f.build_fwtrack()
    track_chrs = track.get_chr_names()
            

    log("#2. Sliding...\n")
    total_num = 0
    range_num = 0
    for k in chromosomes_len.keys():
        log(" on %s..." % (k))
        if not track_chrs.count(k):
            log(" skipped!\n")
            continue
        else:
            log(" computing...\n")
        chr_end = chromosomes_len[k]
        cs = 0                     # current position on chromosome
        ce = cs+options.window
        
        track_k = track.generate_rangeI_by_chr(k)
        track_k_range = track_k.next()
        while ce <= chr_end:
            # do something
            range_num += 1
            num = 0
            while True:
                try:
                    if ce < track_k_range.end:
                        ofhd.write( "%d\n" % (num) )
                        cs = options.gap+ce
                        ce = cs+options.window
                        break
                    else:
                        if cs < track_k_range.start:
                            num += 1
                            total_num+=1
                        track_k_range = track_k.next()
                        
                    #ofhd.write( "%d\n" % (num) )
                except StopIteration: # end of tag list
                    ofhd.write( "%d\n" % (num) )
                    cs = options.gap+ce
                    ce = cs+options.window
                    #ofhd.write( "0\n")
                    break
    ofhd.close()

    sys.stdout.write( "tags: %d\ntotal: %d\nlambda: %f\n" % (total_num,range_num,float(total_num)/range_num) )
    log("Over!\n")

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\n;-) See you!")
        sys.exit(0)
