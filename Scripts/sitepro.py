#!/usr/bin/env python
# 

"""Module Description

Copyright (c) 2008 H. Gene Shin <shin@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

site.py gives an average enrichment profile of given regions of interest (eg, binding sites or motifs).

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
import string
import logging
import re
import itertools
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
error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info

# ------------------------------------
# classes
# ------------------------------------
class WigProfilerwBed:
    """WigProfiler for regions in Bed"""
    
    def __init__(self, span=1000, step=20, dir=False):
        """Constructor"""
        
        # parameters
        self.span = span
        self.step = step
        self.dir = dir
        
        # output values
        self.names = []
        self.profiles = []
        self.xaxis = []
        
        
    def set_params(self, span, step, dir):
        """Set parameters
        
        Parameters:
        1. wig: an Wig object in Cistrome.Assoc.inout
        2. bed: a Bed object in Cistrome.Assoc.inout
        3. span: span size from the center of each region
        4. dir: if set True, the start and end of '-' will be regarded as the end and start.
        """

        self.span=span
        self.step=step
        self.dir=dir
        
        
    def capture_regions(self, wig, bed):
        """Capture the regions that will be profiled"""

        wigint = self.estimate_wig_interval(wig)    # estimate the wig interval for where2 
        step = self.step        
        hf = step/2
        binned = []
        chroms = list(set(wig.get_chroms()).intersection(set(bed.get_chroms())))
        chroms = corelib.sort_chroms(chroms)
        
        for chrom in chroms:
            x=wig[chrom][0]
            y=wig[chrom][1]
            
            # if the direction is considered, include strand info from the BED
            if self.dir:
                try:
                    bediter=itertools.izip(bed[chrom]['start'],bed[chrom]['end'], bed[chrom]['strand'])
                except KeyError:
                    dummy=['' for i in xrange(len(bed[chrom]['start']))]
                    bediter=itertools.izip(bed[chrom]['start'],bed[chrom]['end'], dummy)
            else:
                dummy=['' for i in xrange(len(bed[chrom]['start']))]
                bediter=itertools.izip(bed[chrom]['start'],bed[chrom]['end'], dummy)
            
            init = 0  
            for begin, cease, strand in bediter:
                # get the center, right edge and left edge for search
                center = (begin+cease) /2
                left = center -self.span - hf
                right = center + self.span + hf
                n = (right-left)/step
                
                # get the region from the wig, binning
                start, end= corelib.where2(left, right, x[init:], wigint)
                regionx = x[init+start:init+end]
                regiony = y[init+start:init+end]
                
                if regionx and regiony:
                    bins = corelib.linspace(left, right, n+1)
                    this = corelib.binxy(bins, regionx, regiony, binfunc='middle', NaN=False)
                else:
                    this = [0] * n
                
                # if -, reverse
                if strand == '-':
                    this.reverse()
                    
                # save the binned signal
                binned.append(this)
                    
                # update initial search point   
                init+=start
                
        return binned
    
    
    def estimate_wig_interval(self, wig):
        """Estimate the interval between two consecutive points. 
        This method is exactly the same as estimate_wig_interval function in inout.py.
        
        This methods select randomly 10 regions in each chromosome and take the median of two consecutive intervals.
        """   
        
        chroms = wig.get_chroms()
        if not chroms: return None
    
        n_random_positions = 10
        intervals = []
        for chrom in chroms:
            len_this_chr = len(wig[chrom][0])
            a = corelib.randints(0, len_this_chr - 2, 2 * n_random_positions)   # we need at least two element array to get difference
            a.sort()
            starts = [a[i] for i in xrange(len(a)) if i%2 == 0]
            ends = [a[i] + 2 for i in xrange(len(a)) if i%2 == 1]# we need at least two element array to get difference

            for start, end in itertools.izip(starts, ends):
                intervals.append(corelib.median(corelib.diff(wig[chrom][0][start:end])))
             
        return corelib.median(intervals)
        
        
    def get_breaks(self, start, end):
        """Return breaks for bins
        
        Parameters:
        1. start: the start value
        2. end: the end value. This end value is included in the resulting breaks.
        
        """
        step = self.step
        n = (end-start) / step + 1
        breaks = map(lambda x: int(round(x)), corelib.linspace(start, end, n))
        
        return breaks
       
   
    def profile(self, wig, bed):
        """Wrapper function of WigProfilewBed. 
        Through this function, the user can set the parameters for profiling and get 
        a list of profiles of the regions.
        
        """
        
        start = -1 * self.span
        end = self.span
        return self.get_breaks(start, end), self.capture_regions(wig, bed)
   
   
# ------------------------------------
# functions
# ------------------------------------
def prepare_optparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    
    """
    usage = "usage: %prog <-w wig -b bed> [options]"
    description = "sitepro.py -- Average profile around given genomic sites"
    
    # option processor
    optparser = OptionParser(version="%prog 0.5",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    
    optparser.add_option("-w","--wig",dest="wig",type="string",
                         help="input WIG file. WARNING: only variableStep WIG is accepted.")
    optparser.add_option("-b","--bed",dest="bed",type="string",\
                         help="BED file of regions of interest. (eg, binding sites or motif locations)")
    optparser.add_option("--span",dest="span",type="int",\
                         help="Span from the center of each BED region in both directions(+/-) (eg, [c - span, c + span], where c is the center of a region), default:1000 bp", default=1000)   
    optparser.add_option("--pf-res", dest="pf_res", type="int",\
                          help="Profiling resolution, default: 50 bp", default=50) 
    optparser.add_option("--dir",action="store_true",dest="dir",\
                         help="If set, the direction (+/-) is considered in profiling. If no strand info given in the BED, this option is ignored.",default=False)
    optparser.add_option("--dump",action="store_true",dest="dump",\
                         help="If set, profiles are dumped as a TXT file",default=False)
    optparser.add_option("--name",dest="name",type="string",
                         help="Name of this run")
      
    return optparser


def opt_validate (optparser):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    
    # input BED file and GDB must be given 
    if not (options.wig and options.bed):
        optparser.print_help()
        sys.exit(1)

    if not os.path.isfile(options.wig):
        error('Check -w (--wig). No such file exists: %s' %options.wig)
        sys.exit(1)
        
    if not os.path.isfile(options.bed):
        error('Check -b (--bed). No such file exists:%s' %options.bed)
        sys.exit(1)
            
     # get namename
    if not options.name:
        options.name=os.path.split(options.wig)[-1].rsplit('.wig',2)[0]
        
    return options


def draw_siteprofiles(sitebreaks, avg_siteprof):
    """Return a R script that draws the average profile on the given sites"""
    
    #comment
    R.comment('')
    R.comment('Draw wig profiles around binding sites')
    R.comment('')
    
    #rscript=R.par(mar=[4,3.8,5,4],oma=[4,2,4,2],mfrow=[3,1])
    minlen = min([len(sitebreaks), len(avg_siteprof)])
    rscript=inout.draw_single_profile(sitebreaks[:minlen],avg_siteprof[:minlen],col=["red"],main='Average Profile around the Center of Sites',xlab='Relative Distance from the Center (bp)',ylab='Average Profile',ylim=[],v=0)
    
    return rscript

def dump(chrom, sites, siteprofs):
    """Dump the sites and their profiles in a long string
    
    """
    
    starts = sites['start']
    ends = sites['end']
    
    txt = ''
    for start, end, siteprof in itertools.izip(starts, ends, siteprofs):
        s = map(str, siteprof)
        txt += "%s\t%d\t%d\t%s\n" %(chrom, start, end, ','.join(s)+',')
    
    return txt
    
    
# ------------------------------------
# Main function
# ------------------------------------
def main():
    
    # read the options and validate them
    options=opt_validate(prepare_optparser())
    
    # reading a gene annotation table   
    jobcount=1
    
    # read regions of interest (bed file)
    info("#%d read the bed file of regions of interest..." %jobcount)
    Sites=inout.Bed()
    Sites.read(options.bed)
    chroms_bed=Sites.get_chroms()
    jobcount+=1
        
    # Doing profiling
    chrom=''
    chrcount=1
    swig=inout.Wig()
    FIRST=True

    for line in open(options.wig,'r').xreadlines():
        if not line: continue
        # read a chromosome
        if re.search(r'track',line): 
            try:
                description=re.search(r'description="(\w+)"\s',line).group(1)
            except AttributeError:
                pass
            continue
        if re.search(r'chrom=(\S+)\s',line):
            newchrom=re.search(r'chrom=(\S+)\s',line).group(1)
            try:
                newchrom=inout.standard_chroms[newchrom]
            except KeyError:
                pass
            continue
        l=line.strip().split()

        # the beginning
        if chrom=='' and chrom!=newchrom:
            # if the chromosome is not in gene table, continue
            chrom=newchrom
            
            # if this chromosome is not in chromosome list of the BED, just ignore it
            if chrom in chroms_bed:
            
                info("#%d-%d run wig profiling of %s..." %(jobcount,chrcount,chrom))
                input=inout.Wig()
                input.add_line(chrom,l)
                chrcount+=1

        elif chrom!='' and chrom!=newchrom:    # new chromosome
            
            if chrom in chroms_bed:                
            # wig profiling for given regions of interest
                
                profwbed=WigProfilerwBed(span=options.span, step=options.pf_res, dir=options.dir)
                sitebreaks, siteprofs=profwbed.profile(input, Sites)
                avg_sp, spcount=corelib.mean_col_by_col(siteprofs, counts=True)
                
                if not FIRST:    # if not first chromosome    
                    # average site profiles
                    avg_siteprof,avg_spcount=corelib.weight_mean_col_by_col([avg_siteprof,avg_sp],[avg_spcount,spcount],counts=True)
                    
                    # if --dump, dump the profiles along with their corresponding BED regions
                    if options.dump:
                        dfhd.write(dump(chrom, Sites[chrom], siteprofs))
                    
                    del avg_sp,spcount
                
                else:   # if first chromosome
                    avg_siteprof=avg_sp
                    avg_spcount=spcount
                    
                    # open a txt file to dump profiles if --dump is set
                    if options.dump:
                        dfhd = open(options.name+'.txt', 'w')
                        dfhd.write(dump(chrom, Sites[chrom], siteprofs))
                    
                    FIRST=False
                
                # delete unnucessary variables to maximize usable memory
                del siteprofs   
                
            # set chrom to the new chromosome
            chrom=newchrom
            if chrom in chroms_bed:
                info("#%d-%d run wig profiling of %s..." %(jobcount,chrcount,chrom))
                input=inout.Wig()
                input.add_line(chrom,l)
                chrcount+=1
        else:    # in the middle of chromosome
            if chrom in chroms_bed:
                input.add_line(chrom,l)

    if chrom in chroms_bed:
    # last chromosome
            
        # wig profiling for given regions of interest
        profwbed=WigProfilerwBed(span=options.span, step=options.pf_res, dir=options.dir)
        sitebreaks, siteprofs=profwbed.profile(input, Sites)
        avg_sp,spcount=corelib.mean_col_by_col(siteprofs,counts=True)
        del input             
        
        if not FIRST:    # the first chromosome profiling
            avg_siteprof,avg_spcount=corelib.weight_mean_col_by_col([avg_siteprof,avg_sp],[avg_spcount,spcount],counts=True)
            
            # if --dump, dump the profiles along with their corresponding BED regions
            if options.dump:
                dfhd.write(dump(chrom, Sites[chrom], siteprofs))
            
            del avg_sp,spcount
                    
        else:
            # average site profiles
            avg_siteprof=avg_sp
            avg_spcount=spcount
            
            # if --dump, dump the profiles along with their corresponding BED regions
            if options.dump:
                dfhd = open(options.name+'.txt', 'w')
                dfhd.write(dump(chrom, Sites[chrom], siteprofs))
            
        # delete unnecessary variables
        del siteprofs
     
    if options.dump:
        dfhd.close()
    jobcount+=1

    # write the R script
    info('#%d append R script of wig profiling...' %jobcount)
    ofhd=open(options.name+'.R','w')
    pdfname=options.name+'.pdf'
    rscript=R.pdf(pdfname,height=6, width=8.5)         
    rscript+=draw_siteprofiles(sitebreaks,avg_siteprof)
    ofhd.write(rscript)    # write wig profiling
    ofhd.write(R.devoff())
    ofhd.close()
   
    jobcount+=1
    
    if options.dump:
        info ('#... cong! Run R on %s and see %s!' %(options.name+'.R', options.name+'.txt'))
    else:    
        info ('#... cong! Run R on %s!' % (options.name+'.R'))


# program running
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0)
