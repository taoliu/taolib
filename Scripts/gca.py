#!/usr/bin/env python

"""Module Description

Copyright (c) 2009 H. Gene Shin <shin@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  H. Gene Shin
@contact: shin@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import os
import sys
import re
import logging
import operator
import itertools
from array import array
from optparse import OptionParser
from Cistrome.Assoc import *

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
error   = logging.critical        # function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info


def prepare_optparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    
    """
    
    usage = "usage: %prog [options]"
    description = "GCA -- Gene-centered Annotation"
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-b","--bed",dest="bed",type="string",
                         help="BED file of ChIP regions.")
    optparser.add_option("-g","--gdb",dest="gdb",type="string",
                         help="Gene annotation table. This can be a sqlite3 local db file, BED file or genome version of UCSC. The BED file must have an extension of '.bed'")
    optparser.add_option("--span", dest="span", type="int",\
                         help="Span in search of ChIP regions from TSS and TTS, DEFAULT=3000bp", default=3000)
    optparser.add_option("--name",dest="name",\
                         help="Experiment name. This will be used to name the output file. If an experiment name is not given, input BED file name will be used instead.")      
    
    return optparser


def opt_validate (optparser):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    
    # if gdb not given, print help, either BED or WIG must be given 
    if not options.gdb or not options.bed:
        optparser.print_help()
        sys.exit(0)
        
    if not os.path.isfile(options.gdb):
        if inout.MYSQL:
            options.gdbtype = 'ucsc'
            options.Host = "genome-mysql.cse.ucsc.edu"
            options.User = "genome"
        else:
            error('MySQLDb must be installed to use UCSC gene annotation table')
            sys.exit(0)
    else:
        if options.gdb.endswith('.bed'): # when bed input
            options.gdbtype = 'localbed'
        else:
            options.gdbtype = 'localdb'
            options.Host = None
            options.User = None
        
    if not os.path.isfile(options.bed):
        error("No such file: %" %options.bed)
        sys.exit(0)
        
    if not options.name:
        options.name=os.path.split(options.bed)[-1].rsplit('.bed',2)[0]

    return options


def convert_BED2GeneTable(Bed):
    """Convert a BED-formatted gene annotation table into a GeneTable object
    
    Note that the Bed must strictly follow the BED format not to generate error conversion error.
    """

    # make a GeneTable object
    GeneT = inout.GeneTable()
    GeneT.reset()
    # set the column names of the gene table
    GeneT.columns=('name','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts', 'exonEnds')
    Bed_columns = ('name', 'strand', 'start','end','thickStart','thickEnd', 'blockCount','blockSizes','blockStarts')
    for chrom in Bed.get_chroms():
        GeneT[chrom] = {}
        # copy the columns into the GeneTable
        for gc, bc in itertools.izip(GeneT.columns, Bed_columns):
            GeneT[chrom][gc] = Bed[chrom][bc]
                
    return GeneT    
    
# ------------------------------------
# Main function
# ------------------------------------

def main():
    
    # read the options and validate them
    options=opt_validate(prepare_optparser())
    
    jobcount=1
    info("#%d read the gene table..." %jobcount)
    if options.gdbtype == 'ucsc' or options.gdbtype == 'localdb':
        GeneT = inout.GeneTable()
        GeneT.read(Host = options.Host, User= options.User, Db=options.gdb, annotation='refGene')
    elif options.gdbtype == 'localbed':
        GBed = inout.Bed()
        GBed.read(options.gdb)
        GeneT = convert_BED2GeneTable(GBed)
    GeneT.sort()
    jobcount+=1
    
    info("#%d read the ChIP region BED file..." %jobcount)
    ChIP = inout.Bed()
    ChIP.read(options.bed)
    jobcount+=1
    
    info("#%d perform gene-centered annotation" %jobcount)
    GAnnotator=annotator.GeneAnnotator()
    GAnnotator.annotate(GeneT, ChIP, u=options.span, d=options.span)
    GAnnotator.map.set_name(options.name)
    GAnnotator.map.write()
    jobcount+=1
    
    info("# Gene-centered annotation is done! Check %s" %(options.name+'.xls'))
    

if __name__=="__main__":
    try:
        main()
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0)
