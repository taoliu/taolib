#!/usr/bin/env python
# 

"""Module Description

Copyright (c) 2008 H. Gene Shin <shin@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

gpro.py was made to provide enrichment profiles given a gene annotation 
table from UCSC (or local sqlite3 copy), a wig file of enrichment 
(variableStep format), and a bed file of interest. Unless the bed file 
is given, profiling is performed on specific gennomic loci such as 
promoters, immediate downstreams, and gene body.

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
import operator
from optparse import OptionParser
from Cistrome.Assoc import *
import Cistrome.Assoc.inout as inout

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
    """WigProfiler for regions indicated by Bed"""
    
    def __init__(self,wig=None,bed=None,span=1000,step=20):
        """Constructor"""
        
        # parameters
        self.span=span
        self.wig=wig
        self.bed=bed
        self.step=step
        
        # output values
        self.names=[]
        self.profiles=[]
        self.xaxis=[]
        
    def set_params(self,wig,bed,span,step):
        """Set parameters
        
        Parameters:
        1. wig: an Wig object in Cistrome.Assoc.inout
        2. bed: a Bed object in Cistrome.Assoc.inout
        3. span: span size from the center of each region
        """
        self.wig=wig
        self.bed=bed
        self.span=span
        self.step=step
        
    def capture_regions(self):
        """"""

        hf=self.step/2
        chrom=''
        regions=[]
        chroms=list(set(self.wig.get_chroms()).intersection(set(self.bed.get_chroms())))
        for chrom in chroms:
            init=0
            x=self.wig[chrom][0]
            y=self.wig[chrom][1]
            for begin,cease in itertools.izip(self.bed[chrom]['start'],self.bed[chrom]['end']):
                # get the center, right edge and left edge for search
                center=(begin+cease)/2
                left=center-self.span-hf
                right=center+self.span+hf
                
                # get the region from the wig, binning
                start,end=corelib.where(left,right,x[init:])
                bins=range(left,right+hf,self.step)
                region=corelib.binxy(bins,x[init+start:init+end],y[init+start:init+end])
                init+=start
                
                # if region exists, store
                if region:
                    regions.append(region)
        
        return regions
    
    def get_breaks(self):
        """Return the breaks of the bins"""
        
        hf=self.step/2
        return range(-1*self.span,self.span+hf,self.step)
    
    
    def profile(self,wig,bed,span,step):
        """Wrapper function of WigProfilewBed. 
        Through this function, the user can set the parameters for profiling and get 
        a list of profiles of the regions.
        
        """
        
        self.set_params(wig,bed,span,step)
        
        return self.get_breaks(), self.capture_regions()
   
# ------------------------------------
# functions
# ------------------------------------
def prepare_optparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    
    """
    usage = "usage: %prog <-w wig -g gdb> [options]"
    description = "gpro.py --Enrichment Signal Profiling"
    
    # option processor
    optparser = OptionParser(version="%prog 0.5",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    
    optparser.add_option("-w","--wig",dest="wig",type="string",
                         help="input WIG file")
    optparser.add_option("-g","--gene-db",dest="gdb",type="string",
                         help="Gene annotation database file from Cistrome website")
    optparser.add_option("-b","--bed",dest="bed",type="string",\
                         help="Bed file of regions of interest (eg, binding sites or motif locations). If no bed file is given, profiling only goes on several predefined types of genomic loci such as promoters, immediate downstreams, and gene bodies")
    optparser.add_option("--rel-dist",dest="rel_dist",type="int",
                         help="Relative distance to TSS/TTS for gene profiling, default: 3000 bp", default=3000)
    optparser.add_option("--span",dest="span",type="int",\
                         help="Span from the center of a region of interest. This option only works with the bed file of regions of interest, default:1000 bp",default=1000)
    optparser.add_option("--metagene-size",dest="metagene_size",type="int",
                         help="Metagene size, If not given, the median length of genes will be used", default=None)                      
    optparser.add_option("--pf-res", dest="pf_res", type="int",\
                          help="Profiling resolution, default: 50 bp", default=50) 
    optparser.add_option("--exon-lowers",dest="elowers",type="str",\
                         help="Lower limits of the length (bp) of exons to profile in meta-exon plots. A series of numbers speparated by comma(s) can be given (eg, --exon-lowers=100,300,500). If not given, length values at 10%, 35%, and 65% of exon lengths will be used.", default=None)
    optparser.add_option("--exon-uppers",dest="eupper",type="str",\
                         help="upper limits of the length (bp) of introns to profile in meta-intron plots. A series of numbers speparated by comma(s) can be given (eg, --exon-lowers=300,500,1000). If not given, length values at 35%, 65%, and 90% of intron length will be used.",default=None)
    optparser.add_option("--intron-lowers",dest="ilowers",type="str",\
                         help="Lower limits of the length (bp) of introns to profile in meta-intron plots. A series of numbers speparated by comma(s) can be given (eg, --exon-lowers=100,300,500). If not given, length values at 10%, 35%, and 65% of intron length will be used.", default=None)
    optparser.add_option("--intron-uppers",dest="iupper",type="str",\
                         help="upper limits of the length (bp) of exons to profile in meta-exon plots. A series of numbers speparated by comma(s) can be given (eg, --exon-lowers=300,500,1000). If not given, length values at 35%, 65%, and 90% of exon lengths will be used.",default=None)
    optparser.add_option("--gn-groups",dest="gn_groups",type="string",\
                         help="Gene-group file names (eg, top10.txt,bottom10.txt), default: all genes", default=None) 
    optparser.add_option("--gn-group-names", dest="gn_names",type="string",\
                         help="Gene-group names for profiling (eg, top 10%,bottom 10%), default: None",default=None)
    optparser.add_option("--no-refseq", action="store_true", dest="name2",\
                         help="Whether gene names (eg, 'name2' in refGene of UCSC) or RefSeq accession IDs (eg, 'name' in refGene of UCSC) are used in --gn-groups or not. This flag is meaningful only if --gn-groups is set.",default=False)
    optparser.add_option("--name",dest="name",type="string",
                         help="Name of this run. ")
    optparser.add_option("--verbose",dest="verbose",type="string",
                         help="Name of a directory. if set, save verbose information in the direcotry.", default=None)
      
    return optparser


def opt_validate (optparser):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    
    # input BED file and GDB must be given 
    if not (options.wig and options.gdb):
        optparser.print_help()
        sys.exit(1)

    if not os.path.isfile(options.wig):
        error('Check -w (--wig). No such file exists' %options.wig)
        sys.exit(1)
        
    # get gdb lower case
    HAVELOCALGDB=os.path.isfile(options.gdb)
    if HAVELOCALGDB:
        options.Host=None
        options.User=None
        options.Db=options.gdb.lower()
    else:
        if inout.MYSQL:
            options.Host="genome-mysql.cse.ucsc.edu"
            options.User="genome"
            options.Db=os.path.split(options.gdb)[-1].lower()
        else:
            error('MySQLdb package needs to be installed to use UCSC or a local sqlite3 db file must exist.')
            error('Check -g (--gdb). No such file or species exists: %s' %options.gdb)
            sys.exit(1)
                
    # regions of interest
    if options.bed:
        if not os.path.isfile(options.bed):
            error('Check -b (--bed). No such file exists: %s' %options.ebed)
            sys.exit(1)
            
    # get the median length of genes in case metagene size is not given
    if options.metagene_size:
        options.catexon_size = int(round(0.5 * options.metagene_size))
        options.catintron_size = options.catexon_size
    
     # get namename
    if not options.name:
        options.name=os.path.split(options.wig)[-1].rsplit('.wig',2)[0]
    
    #check if name2 is going to be used instead of name
    if options.name2:
        options.which='name2'
    else:
        options.which='name'
                    
    # check the gene group files    
    if options.gn_groups:
        parsed=options.gn_groups.split(',')
        for p in parsed:
            if not os.path.isfile(p):
                error('Check --gn-groups. No such files exist: %s' %p)
                sys.exit(1)
        options.gn_groups=parsed
        
        if options.gn_names:
            parsed_names=options.gn_names.split(',')
            if len(parsed_names)!=len(options.gn_groups):
                error('There must be the same number of group names as gene groups')
                sys.exit(1)
            options.gn_names=parsed_names
        else:
            options.gn_names=[]
            for i in range(len(options.gn_groups)):
                options.gn_names.append('Group %d' %(i+1))
    
    # exon length and intron lengths
    # consider exons with lengths of [lower lim (%), upper lim (%)]
    # consider introns with lengths of [lower lim (%), upper lim (%)]
#    options.epercentlim=[50, 90]    # in percent
#    options.ipercentlim=[50, 90]    # in percent
#    options.minexonlen = options.pf_res*2   # minimum exon length to consider
#    options.minintronlen = options.pf_res*2   # minimum exon length to consider
#    options.metaexonsize = options.pf_res*10 # the meta exon and intron length for plotting
#    options.metaintronsize = options.pf_res*10 # the meta exon and intron length for plotting

    #fixed metagene size
    options.metagene_size = 3000
    options.catexon_size = options.metagene_size/2
    options.catintron_size = options.metagene_size/2
    options.epercentlowers = [10, 35, 65]
    options.epercentuppers = [35, 65, 90]
    options.ipercentlowers = [10, 35, 65]
    options.ipercentuppers = [35, 65, 90]

#    options.elowers = [30]
#    options.euppers = [99]
#    options.ilowers = [30]
#    options.iuppers = [99]
    
    # number of wig overview points
    options.nwigpts=2000
    
    
    return options

def draw_siteprofiles(sitebreaks,avg_siteprof):
    """Return a R script that draws the average profile on the given sites"""
    
    #comment
    R.comment('')
    R.comment('Draw wig profiles on binding sites')
    R.comment('')
    
    rscript=R.par(mar=[4,3.8,5,4],oma=[4,2,4,2],mfrow=[3,1])
    rscript+=inout.draw_single_profile(sitebreaks,avg_siteprof,col=["red"],main='Average Enrichment on the Sites',xlab='Relative Distance wrt the Centers',ylab='Average Enrichment',ylim=[],v=0)
    
    return rscript


def get_min_max_lengths(GeneT, epercentlimit, ipercentlimit, minexonlen, minintronlen):
    """Return the minimum and maximum lengths of exons or introns in the gene annotation table
    
    Parameters:
    1. GeneT: gene annotation table
    2. epercentlimit: [lower limit for exon length, upper limit for exon length] to consider. lower limit and upper limit must be percentage values (0-100)
    3. ipercentlimit: [lower limit for intron length, upper limit for exon lentht] to consider. 
    4. minexonlen: mininum exon length to consider
    5. minintronlen: minimum intron length to consider
    """
    
    exonLens, intronLens =GeneT.get_exon_intron_lengths()
    exon_enough_long, intron_enough_long = corelib.get_certain_part(exonLens, percentlimit=epercentlimit), corelib.get_certain_part(intronLens, percentlimit=ipercentlimit)
    exonLenLim, intronLenLim = (max(minexonlen, exon_enough_long[0]), exon_enough_long[-1]), (max(minintronlen, intron_enough_long[0]), intron_enough_long[-1])
    #exonLenLim = (500, 1500)
    #intronLenLim = (500, 1500)
    return exonLenLim, intronLenLim



# ------------------------------------
# Main function
# ------------------------------------
def main():
    
    # read the options and validate them
    options=opt_validate(prepare_optparser())
    # reading a gene annotation table   
    jobcount=1
    info("#%d read the gene table from db..." %jobcount)
    
    # read the gene table
    GeneT=inout.GeneTable()
    GeneT.read(Host=options.Host,User=options.User,Db=options.Db,annotation='refGene',which=options.which)    
    GeneT.sort()
    chroms_GeneT=GeneT.get_chroms()
    
    # determine the metagenesize, concated exon size, concatenated intron size
    gLens = GeneT.get_gene_lens()
    catexonLens, catintronLens =GeneT.get_cat_exon_intron_lens()
    medgLen = corelib.median(gLens)
#    medcatexonLen, medcatintronLen = corelib.median(catexonLens), corelib.median(catintronLens)
    if options.metagene_size:
#        catexonRat, catintronRat = 1.0 * medcatexonLen/medgLen, 1.0 * medcatintronLen/medgLen
#        options.catexon_size, options.catintron_size = int(round(catexonRat * options.metagene_size)), int(round(catintronRat * options.metagene_size))
        options.catexon_size, options.catintron_size = options.metagene_size/2, options.metagene_size/2
    else:    
        options.metagene_size = medgLen
#        options.catexon_size, options.catintron_size = medcatexonLen, medcatintronLen
        options.catexon_size, options.catintron_size = medgLen/2, medgLen/2
    
    # delete unnecessary variables
    del gLens, catexonLens, catintronLens
    
    # get the exon lengths and intron lengths
    if not (options.elowers and options.euppers and options.ilowers and options.iuppers):
        exonLens, intronLens = GeneT.get_exon_intron_lens()
        options.elowers, medexonsizes, options.euppers = corelib.get_boundaries_medians(exonLens, lowers = options.epercentlowers, uppers = options.epercentuppers)
        options.ilowers, medintronsizes, options.iuppers = corelib.get_boundaries_medians(intronLens, lowers = options.ipercentlowers, uppers = options.ipercentuppers)
    n_ranges = len(options.elowers)
    options.metaexonsizes = map(max, medexonsizes, [2*options.pf_res]*n_ranges)
    options.metaintronsizes = map(max, medintronsizes, [2*options.pf_res]*n_ranges)
    jobcount+=1
    
    # determine the metagenesize, concated exon size, concatenated intron size
#    gLens = GeneT.get_gene_lens()
#    catexonLens, catintronLens =GeneT.get_cat_exon_intron_lens()
#    medgLen = corelib.median(gLens)
#    medcatexonLen, medcatintronLen = corelib.median(catexonLens), corelib.median(catintronLens)
#    if options.metagene_size:
#        catexonRat, catintronRat = 1.0 * medcatexonLen/medgLen, 1.0 * medcatintronLen/medgLen
#        options.catexon_size, options.catintron_size = int(round(catexonRat * options.metagene_size)), int(round(catintronRat * options.metagene_size))
#    else:    
#        options.metagene_size = medgLen
#        options.catexon_size, options.catintron_size = medcatexonLen, medcatintronLen
#    
#    # delete unnecessary variables
#    del gLens, catexonLens, catintronLens
    
#    # get the exon lengths and intron lengths
#    if not (options.elowers and options.euppers and options.ilowers and options.iuppers):
#        exonLens, intronLens = GeneT.get_exon_intron_lens()
#        options.elowers, medexonsizes, options.euppers = corelib.get_boundaries_medians(exonLens, lowers = options.epercentlowers, uppers = options.epercentuppers)
#        options.ilowers, medintronsizes, options.iuppers = corelib.get_boundaries_medians(intronLens, lowers = options.ipercentlowers, uppers = options.ipercentuppers)
#    n_ranges = len(options.elowers)
##    options.metaexonsizes = map(max, options.metaexonsizes, [2*options.pf_res]*n_ranges)
##    options.metaintronsizes = map(max, options.metaintronsizes, [2*options.pf_res]*n_ranges)
#    options.metaexonsizes = [10*options.pf_res] * n_ranges
#    options.metaintronsizes = [10*options.pf_res] * n_ranges

    # the minimum value should be larger than the resolution
#    options.elowers = map(max, options.elowers, [2*options.pf_res]*n_ranges)
#    options.euppers = map(max, options.euppers, [2*options.pf_res]*n_ranges)
#    options.ilowers = map(max, options.ilowers, [2*options.pf_res]*n_ranges)
#    options.iuppers = map(max, options.iuppers, [2*options.pf_res]*n_ranges)
#    options.metaexonsizes = map(max, options.metaexonsizes, [2*options.pf_res]*n_ranges)
#    options.metaintronsizes = map(max, options.metaintronsizes, [2*options.pf_res]*n_ranges)
    jobcount+=1
    
    # read reginos of interest (bed file)
    if options.bed:
        info("#%d read the bed file of regions of interest..." %jobcount)
        Sites=inout.Bed()
        Sites.read(options.bed)
        jobcount+=1
        
    # if gene groups are give
    if options.gn_groups:
        subsets=inout.read_gene_subsets(options.gn_groups)
        
    chrom=''
    chrcount=1
    prof=profiler.WigProfiler2(rel_dist=options.rel_dist, step=options.pf_res, metagenesize=options.metagene_size, \
                               catexonsize=options.catexon_size, catintronsize=options.catintron_size, metaexonsizes=options.metaexonsizes, \
                               metaintronsizes=options.metaintronsizes, elowers=options.elowers, euppers=options.euppers, ilowers=options.ilowers, \
                               iuppers=options.iuppers)
    
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
            if chrom in chroms_GeneT: # only if the new chromosome is in the chroms of gene table, a wig object is initiated.
                
                info("#%d-%d run wig profiling of %s..." %(jobcount,chrcount,chrom))
                input=inout.Wig()
                input.add_line(chrom,l)
                chrcount+=1

        elif chrom!='' and chrom!=newchrom:    # new chromosome
            if chrom in chroms_GeneT:
                    
                # wig overview
                #
                ws=sampler.WigSampler()
                resolution=input.chrom_length(chrom)/options.nwigpts
                swig+=sampler.fillupwig(ws.sample(input,resolution=resolution),resolution=resolution)
                                                             
                # wig profiling
                profiles = prof.profile(input, GeneT)
#                names,breaks,upstreams,downstreams,metagene_breaks,metagenes,metacatexon_breaks,metacatexons,metacatintron_breaks,metacatintrons,metaexon_breaks,metaexons,metaintron_breaks,metaintrons=\
#                prof.profile(input,GeneT,rel_pos=options.rel_dist,metagenesize=options.metagene_size,step=options.pf_res,which=options.which,exonratio=0.5,emask=options.emask,imask=options.imask,\
#                             elenlim=exonLenLim,ilenlim=intronLenLim,metaexonsize=options.metaexonsize,metaintronsize=options.metaintronsize)
                
                n_metas = len(profiles['exons']) 
                # get average of this chromosome
                avg_up,upcount=corelib.mean_col_by_col(profiles['upstreams'], counts=True)
                avg_down,downcount=corelib.mean_col_by_col(profiles['downstreams'], counts=True)
                avg_mg,genecount=corelib.mean_col_by_col(profiles['genes'], counts=True)
                avg_mce,cexoncount=corelib.mean_col_by_col(profiles['catexons'], counts=True)
                avg_mci,cintroncount=corelib.mean_col_by_col(profiles['catintrons'],counts=True)
                
#                avg_me,exoncount=corelib.mean_col_by_col(corelib.extend_list_series(profiles['exons']), counts=True)
#                avg_mi,introncount=corelib.mean_col_by_col(corelib.extend_list_series(profiles['introns']), counts=True)
                outs = map(corelib.mean_col_by_col, map(corelib.extend_list_series, profiles['exons']), [True] * len(profiles['exons']))
                avg_me = map(operator.itemgetter(0), outs)
                exoncount = map(operator.itemgetter(1), outs)
                outs = map(corelib.mean_col_by_col, map(corelib.extend_list_series, profiles['introns']), [True] * len(profiles['introns']))
                avg_mi = map(operator.itemgetter(0), outs)
                introncount = map(operator.itemgetter(1), outs)
                del outs
                
                # wig profiling for given regions of interest
                if options.bed:
                    profwbed=WigProfilerwBed()
                    sitebreaks,siteprofs=profwbed.profile(input,Sites,options.span,options.pf_res)
                
                    # get the average of site profiles
                    avg_sp,spcount=corelib.mean_col_by_col(siteprofs,counts=True)
                
                if not FIRST:    # if not first chromosome
                    avg_upstream,avg_upcount=corelib.weight_mean_col_by_col([avg_upstream,avg_up],[avg_upcount,upcount],counts=True)
                    avg_downstream,avg_downcount=corelib.weight_mean_col_by_col([avg_downstream,avg_down],[avg_downcount,upcount],counts=True)
                    avg_metagene,avg_genecount=corelib.weight_mean_col_by_col([avg_metagene,avg_mg],[avg_genecount,genecount],counts=True)
                    avg_metacatexon,avg_cexoncount=corelib.weight_mean_col_by_col([avg_metacatexon,avg_mce],[avg_cexoncount,cexoncount],counts=True)
                    avg_metacatintron,avg_cintroncount=corelib.weight_mean_col_by_col([avg_metacatintron,avg_mci],[avg_cintroncount,cintroncount],counts=True)
#                    avg_metaexon,avg_exoncount=corelib.weight_mean_col_by_col([avg_metaexon,avg_me],[avg_exoncount,exoncount],counts=True)
#                    avg_metaintron,avg_introncount=corelib.weight_mean_col_by_col([avg_metaintron,avg_mi],[avg_introncount,introncount],counts=True)
                    outs = map(corelib.weight_mean_col_by_col, map(lambda x, y: [x, y], avg_metaexon, avg_me),  map(lambda x, y: [x, y], avg_exoncount, exoncount), [True] * n_metas)
                    avg_metaexon = map(operator.itemgetter(0), outs)
                    avg_exoncount = map(operator.itemgetter(1), outs)
                    outs = map(corelib.weight_mean_col_by_col, map(lambda x, y: [x, y], avg_metaintron, avg_mi), map(lambda x, y: [x, y], avg_introncount, introncount), [True] * n_metas)
                    avg_metaintron = map(operator.itemgetter(0), outs)
                    avg_introncount = map(operator.itemgetter(1), outs)
                    del outs
                    
                    if options.bed:
                        # average site profiles
                        avg_siteprof,avg_spcount=corelib.weight_mean_col_by_col([avg_siteprof,avg_sp],[avg_spcount,spcount],counts=True)
                        del avg_sp,spcount
                    
                    del avg_up,avg_down,avg_mg,avg_mce,avg_mci,avg_me,avg_mi,upcount,downcount,genecount,cexoncount,cintroncount,exoncount,introncount
                
                    if options.gn_groups:    # when gene sub-gropus are given
                        
                        if options.which=='name2':
                            ixs, subsets = profiler.get_gene_indicies(profiles['name2'], subsets)
                        else:
                            ixs, subsets = profiler.get_gene_indicies(profiles['name'], subsets)
            
#                        avg_ups, upcs, avg_downs, downcs, avg_mgs, gcs, avg_mces, cecs, avg_mcis, cics, avg_mes, ecs, avg_mis, ics=\
#                        profiler.select_profiles_chr_by_chr(ixs, profiles['upstreams'], profiles['downstreams'], profiles['genes'], profiles['catexons'], profiles['catintrons'], profiles['exons'], profiles['introns'])
                        avg_ups, upcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['upstreams'])
                        avg_downs, downcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['downstreams'])
                        avg_mgs, gcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['genes'])
                        avg_mces, cecs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catexons'])
                        avg_mcis, cics = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catintrons'])
#                        outs = map(profiler.select_take_average_profiles_chr_by_chr, [ixs] *n_metas, map(corelib.extend_list_series, profiles['exons']))
#                        avg_mes, ecs = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                        outs = map(profiler.select_take_average_profiles_chr_by_chr, [ixs] *n_metas, map(corelib.extend_list_series, profiles['introns']))
#                        avg_mis, ics = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                        del outs
                        avg_mes, ecs = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['exons'])
                        avg_mis, ics = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['introns'])
                        
                        avg_upstreams,avg_upcounts=profiler.weight_mean_profiles_chr_by_chr(avg_upstreams,avg_upcounts,avg_ups,upcs)
                        avg_downstreams,avg_downcounts=profiler.weight_mean_profiles_chr_by_chr(avg_downstreams,avg_downcounts,avg_downs,downcs)
                        avg_metagenes,avg_genecounts=profiler.weight_mean_profiles_chr_by_chr(avg_metagenes,avg_genecounts,avg_mgs,gcs)
                        avg_metacatexons,avg_cexoncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metacatexons,avg_cexoncounts,avg_mces,cecs)
                        avg_metacatintrons,avg_cintroncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metacatintrons,avg_cintroncounts,avg_mcis,cics)
#                        avg_metaexons,avg_exoncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metaexons,avg_exoncounts,avg_mes,ecs)
#                        avg_metaintrons,avg_introncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metaintrons,avg_introncounts,avg_mis,ics)
#                        outs = map(profiler.weight_mean_profiles_chr_by_chr, avg_metaexons, avg_exoncounts, avg_mes, ecs)
#                        avg_metaexons, avg_exoncounts = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                        outs = map(profiler.weight_mean_profiles_chr_by_chr, avg_metaintrons, avg_introncounts, avg_mis, ics)
#                        avg_metaintrons, avg_introncounts = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                        del outs
                        avg_metaexons, avg_exoncounts = profiler.weight_mean_profiles_chr_by_chr_meta(avg_metaexons, avg_exoncounts, avg_mes, ecs)
                        avg_metaintrons, avg_introncounts = profiler.weight_mean_profiles_chr_by_chr_meta(avg_metaintrons, avg_introncounts, avg_mis, ics)
                        del avg_ups,avg_downs,avg_mgs,avg_mces,avg_mcis,avg_mes,avg_mis,upcs,downcs,gcs,cecs,cics,ecs,ics
                
                else:   # if first chromosome
                    avg_upstream = avg_up
                    avg_downstream = avg_down
                    avg_metagene = avg_mg
                    avg_metacatexon = avg_mce
                    avg_metacatintron = avg_mci
                    avg_metaexon = avg_me
                    avg_metaintron = avg_mi
                    avg_upcount = upcount
                    avg_downcount = downcount
                    avg_genecount = genecount
                    avg_cexoncount = cexoncount
                    avg_cintroncount = cintroncount
                    avg_exoncount=exoncount
                    avg_introncount=introncount
                    
                    
                    if options.bed:
                        # average site profiles
                        avg_siteprof=avg_sp
                        avg_spcount=spcount
            
                    if options.gn_groups:
                        if options.which=='name2':
                            ixs, subsets = profiler.get_gene_indicies(profiles['name2'], subsets)
                        else:
                            ixs, subsets = profiler.get_gene_indicies(profiles['name'], subsets)
#                        avg_upstreams, avg_upcounts, avg_downstreams, avg_downcounts, avg_metagenes, avg_genecounts, avg_metacatexons, avg_cexoncounts, avg_metacatintrons, avg_cintroncounts, avg_metaexons, avg_exoncounts, avg_metaintrons, avg_introncounts =\
#                        profiler.select_profiles_chr_by_chr(ixs, profiles['upstreams'], profiles['downstreams'], profiles['genes'], profiles['catexons'], profiles['catintrons'], profiles['exons'], profiles['introns'])
                        avg_upstreams, avg_upcounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['upstreams'])
                        avg_downstreams, avg_downcounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['downstreams'])
                        avg_metagenes, avg_genecounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['genes'])
                        avg_metacatexons, avg_cexoncounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catexons'])
                        avg_metacatintrons, avg_cintroncounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catintrons'])
#                        outs = map(profiler.select_take_average_profiles_chr_by_chr, [ixs] *n_metas, map(corelib.extend_list_series, profiles['exons']))
#                        avg_metaexons, avg_exoncounts = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                        outs = map(profiler.select_take_average_profiles_chr_by_chr, [ixs] *n_metas,  map(corelib.extend_list_series, profiles['introns']))
#                        avg_metaintrons, avg_introncounts = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                        del outs
                        avg_metaexons, avg_exoncounts = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['exons'])
                        avg_metaintrons, avg_introncounts = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['introns'])
                        
                    FIRST=False
                
                del profiles
                
                # delete unnucessary variables to maximize usable memory

                if options.bed:
                    del siteprofs   
                
            # set chrom to the new chromosome
            chrom=newchrom
            if chrom in chroms_GeneT:    # only if the new chromosome is in the chroms of gene table, a wig object is initiated.
               
                info("#%d-%d run wig profiling of %s..." %(jobcount,chrcount,chrom))
                input=inout.Wig()
                input.add_line(chrom,l)
                chrcount+=1
        else:    # in the middle of chromosome
            if chrom in chroms_GeneT:   # only if the new chromosome is in the chroms of gene table, the wig object is updated.
                input.add_line(chrom,l)

    # last chromosome
    if chrom in chroms_GeneT:
            
        # wig overview
        ws=sampler.WigSampler()
        resolution=input.chrom_length(chrom)/options.nwigpts
        swig+=sampler.fillupwig(ws.sample(input,resolution=resolution),resolution=resolution)
            
        # profiling
        profiles = prof.profile(input, GeneT)
#        names,breaks,upstreams,downstreams,metagene_breaks,metagenes,metacatexon_breaks,metacatexons,metacatintron_breaks,metacatintrons,metaexon_breaks,metaexons,metaintron_breaks,metaintrons=\
#        prof.profile(input,GeneT,rel_pos=options.rel_dist,metagenesize=options.metagene_size,step=options.pf_res,which=options.which,exonratio=0.5,emask=options.emask,imask=options.imask,\
#                     elenlim=exonLenLim,ilenlim=intronLenLim,metaexonsize=options.metaexonsize,metaintronsize=options.metaintronsize)
        
        
        n_metas = len(profiles['exons'])  
        # get average of this chromosome
        avg_up,upcount = corelib.mean_col_by_col(profiles['upstreams'], counts=True)
        avg_down,downcount = corelib.mean_col_by_col(profiles['downstreams'], counts=True)
        avg_mg,genecount = corelib.mean_col_by_col(profiles['genes'], counts=True)
        avg_mce,cexoncount = corelib.mean_col_by_col(profiles['catexons'], counts=True)
        avg_mci,cintroncount = corelib.mean_col_by_col(profiles['catintrons'],counts=True)
#        avg_me,exoncount = corelib.mean_col_by_col(corelib.extend_list_series(profiles['exons']), counts=True)
#        avg_mi,introncount = corelib.mean_col_by_col(corelib.extend_list_series(profiles['introns']), counts=True)
        outs = map(corelib.mean_col_by_col, map(corelib.extend_list_series, profiles['exons']), [True] * n_metas)
        avg_me = map(operator.itemgetter(0), outs)
        exoncount = map(operator.itemgetter(1), outs)
        outs = map(corelib.mean_col_by_col, map(corelib.extend_list_series, profiles['introns']), [True] * n_metas)
        avg_mi = map(operator.itemgetter(0), outs)
        introncount = map(operator.itemgetter(1), outs)
        
        # wig profiling for given regions of interest
        if options.bed:
            profwbed=WigProfilerwBed()
            sitebreaks,siteprofs=profwbed.profile(input,Sites,options.span,options.pf_res)
    
            # get the average of site profiles
            avg_sp,spcount=corelib.mean_col_by_col(siteprofs,counts=True)
        
        del input             
        if not FIRST:    # the first chromosome profiling
            avg_upstream,avg_upcount=corelib.weight_mean_col_by_col([avg_upstream,avg_up],[avg_upcount,upcount],counts=True)
            avg_downstream,avg_downcount=corelib.weight_mean_col_by_col([avg_downstream,avg_down],[avg_downcount,upcount],counts=True)
            avg_metagene,avg_genecount=corelib.weight_mean_col_by_col([avg_metagene,avg_mg],[avg_genecount,genecount],counts=True)
            avg_metacatexon,avg_cexoncount=corelib.weight_mean_col_by_col([avg_metacatexon,avg_mce],[avg_cexoncount,cexoncount],counts=True)
            avg_metacatintron,avg_cintroncount=corelib.weight_mean_col_by_col([avg_metacatintron,avg_mci],[avg_cintroncount,cintroncount],counts=True)
#            avg_metaexon,avg_exoncount=corelib.weight_mean_col_by_col([avg_metaexon,avg_me],[avg_exoncount,exoncount],counts=True)
#            avg_metaintron,avg_introncount=corelib.weight_mean_col_by_col([avg_metaintron,avg_mi],[avg_introncount,introncount],counts=True)
            
            outs = map(corelib.weight_mean_col_by_col, map(lambda x, y: [x, y], avg_metaexon, avg_me),  map(lambda x, y: [x, y], avg_exoncount, exoncount), [True] * n_metas)
            avg_metaexon = map(operator.itemgetter(0), outs)
            avg_exoncount = map(operator.itemgetter(1), outs)
            outs = map(corelib.weight_mean_col_by_col, map(lambda x, y: [x, y], avg_metaintron, avg_mi), map(lambda x, y: [x, y], avg_introncount, introncount), [True] * n_metas)
            avg_metaintron = map(operator.itemgetter(0), outs)
            avg_introncount = map(operator.itemgetter(1), outs)
            del outs
                    
            if options.bed:
                # average site profiles
                avg_siteprof,avg_spcount=corelib.weight_mean_col_by_col([avg_siteprof,avg_sp],[avg_spcount,spcount],counts=True)
                del avg_sp,spcount
            
            del avg_up,avg_down,avg_mg,avg_mce,avg_mci,avg_me,avg_mi,upcount,downcount,genecount,cexoncount,cintroncount,exoncount,introncount
        
            if options.gn_groups:
                if options.which=='name2':
                    ixs, subsets = profiler.get_gene_indicies(profiles['name2'], subsets)
                else:
                    ixs, subsets = profiler.get_gene_indicies(profiles['name'], subsets)
                    
                # take an average of each profile (upstream, downstream, gene, exon, intron, cat-exon, cat-intron
#                avg_ups, upcs, avg_downs, downcs, avg_mgs, gcs, avg_mces, cecs, avg_mcis, cics, avg_mes, ecs, avg_mis, ics=\
#                profiler.select_profiles_chr_by_chr(ixs, profiles['upstreams'], profiles['downstreams'], profiles['genes'], profiles['catexons'], profiles['catintrons'], profiles['exons'], profiles['introns'])
                avg_ups, upcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['upstreams'])
                avg_downs, downcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['downstreams'])
                avg_mgs, gcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['genes'])
                avg_mces, cecs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catexons'])
                avg_mcis, cics = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catintrons'])
#                outs = map(profiler.select_take_average_profiles_chr_by_chr, [ixs] *n_metas, map(corelib.extend_list_series, profiles['exons']))
#                avg_mes, ecs = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                outs = map(profiler.select_take_average_profiles_chr_by_chr, [ixs] *n_metas, map(corelib.extend_list_series, profiles['introns']))
#                avg_mis, ics = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                del outs
                avg_mes, ecs = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['exons'])
                avg_mis, ics = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['introns'])
                    
                   
                avg_upstreams,avg_upcounts=profiler.weight_mean_profiles_chr_by_chr(avg_upstreams,avg_upcounts,avg_ups,upcs)
                avg_downstreams,avg_downcounts=profiler.weight_mean_profiles_chr_by_chr(avg_downstreams,avg_downcounts,avg_downs,downcs)
                avg_metagenes,avg_genecounts=profiler.weight_mean_profiles_chr_by_chr(avg_metagenes,avg_genecounts,avg_mgs,gcs)
                avg_metacatexons,avg_cexoncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metacatexons,avg_cexoncounts,avg_mces,cecs)
                avg_metacatintrons,avg_cintroncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metacatintrons,avg_cintroncounts,avg_mcis,cics)
#                avg_metaexons,avg_exoncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metaexons,avg_exoncounts,avg_mes,ecs)
#                avg_metaintrons,avg_introncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metaintrons,avg_introncounts,avg_mis,ics)
#                outs = map(profiler.weight_mean_profiles_chr_by_chr, avg_metaexons, avg_exoncounts, avg_mes, ecs)
#                avg_metaexons, avg_exoncounts = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                outs = map(profiler.weight_mean_profiles_chr_by_chr, avg_metaintrons, avg_introncounts, avg_mis, ics)
#                avg_metaintrons, avg_introncounts = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                del outs
                avg_metaexons, avg_exoncounts = profiler.weight_mean_profiles_chr_by_chr_meta(avg_metaexons, avg_exoncounts, avg_mes, ecs)
                avg_metaintrons, avg_introncounts = profiler.weight_mean_profiles_chr_by_chr_meta(avg_metaintrons, avg_introncounts, avg_mis, ics)
               
                del avg_ups,avg_downs,avg_mgs,avg_mces,avg_mcis,avg_mes,avg_mis,upcs,downcs,gcs,cecs,cics,ecs,ics
        else:
            avg_upstream=avg_up
            avg_downstream=avg_down
            avg_metagene=avg_mg
            avg_metacatexon=avg_mce
            avg_metacatintron=avg_mci
            avg_metaexon=avg_me
            avg_metaintron=avg_mi
            avg_upcount=upcount
            avg_downcount=downcount
            avg_genecount=genecount
            avg_cexoncount=cexoncount
            avg_cintroncount=cintroncount
            avg_exoncount=exoncount
            avg_introncount=introncount
            
            if options.bed:
                # average site profiles
                avg_siteprof=avg_sp
                avg_spcount=spcount
            
            if options.gn_groups:
                if options.which=='name2':
                    ixs, subsets = profiler.get_gene_indicies(profiles['name2'], subsets)
                else:
                    ixs, subsets = profiler.get_gene_indicies(profiles['name'], subsets)
#                avg_upstreams, avg_upcounts, avg_downstreams, avg_downcounts, avg_metagenes, avg_genecounts, avg_metacatexons, avg_cexoncounts, avg_metacatintrons, avg_cintroncounts, avg_metaexons, avg_exoncounts, avg_metaintrons, avg_introncounts =\
#                profiler.select_profiles_chr_by_chr(ixs, profiles['upstreams'], profiles['downstreams'], profiles['genes'], profiles['catexons'], profiles['catintrons'], profiles['exons'], profiles['introns'])
                avg_upstreams, avg_upcounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['upstreams'])
                avg_downstreams, avg_downcounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['downstreams'])
                avg_metagenes, avg_genecounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['genes'])
                avg_metacatexons, avg_cexoncounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catexons'])
                avg_metacatintrons, avg_cintroncounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catintrons'])
#                outs = map(profiler.select_take_average_profiles_chr_by_chr, [ixs] *n_metas, map(corelib.extend_list_series, profiles['exons']))
#                avg_metaexons, avg_exoncounts = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                outs = map(profiler.select_take_average_profiles_chr_by_chr, [ixs] *n_metas,  map(corelib.extend_list_series, profiles['introns']))
#                avg_metaintrons, avg_introncounts = map(operator.itemgetter(0), outs), map(operator.itemgetter(1), outs)
#                del outs
                avg_metaexons, avg_exoncounts = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['exons'])
                avg_metaintrons, avg_introncounts = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['introns'])
                
        # delete unnecessary variables
        if options.bed:
            del siteprofs
    jobcount+=1

    # writing R script
    info('#%d write an R script of wig profiling...' %jobcount)
    ofhd=open(options.name+'.R','w')
    pdfname=options.name+'.pdf'
    rscript=R.pdf(pdfname,height=11.5,width=8.5)
    
    #
    # write a R script of drawing wig overview
    #
    rscript+=inout.draw_wig_heatmap_overview(swig)
    
    #
    # write a R script of drawing profiles near genes
    #
    
    # get breaks for the plots
    breaks = profiles['breaks']
    metagene_breaks = profiles['genebreaks']
    metacatexon_breaks = profiles['catexonbreaks']
    metacatintron_breaks = profiles['catintronbreaks']
    metaexon_breaks = profiles['exonbreaks']
    metaintron_breaks = profiles['intronbreaks']
    
    #
    # when gene groups of interest are given and not given
    #
    if options.gn_groups:           
        # append the profiles of all genes
        avg_upstreams.append(avg_upstream)
        avg_downstreams.append(avg_downstream)
        avg_metagenes.append(avg_metagene)
        avg_metacatexons.append(avg_metacatexon)
        avg_metacatintrons.append(avg_metacatintron)
        map(lambda x, y: x.append(y), avg_metaexons, avg_metaexon)
        map(lambda x, y: x.append(y), avg_metaintrons, avg_metaintron)  
        rscript+=inout.draw_profile_plots(breaks,avg_upstreams,avg_downstreams,metagene_breaks,avg_metagenes,metacatexon_breaks,avg_metacatexons,metacatintron_breaks,avg_metacatintrons,metagene_breaks_lim=[-1000,1000],legends=options.gn_names)
        rscript+=inout.draw_exon_intron_profile_plots(metaexon_breaks, avg_metaexons,metaintron_breaks,avg_metaintrons, options.elowers, options.euppers, options.ilowers, options.iuppers, legends=options.gn_names)
    else:                 
        rscript+=inout.draw_profile_plot(breaks,avg_upstream,avg_downstream,metagene_breaks,avg_metagene,metacatexon_breaks,avg_metacatexon,metacatintron_breaks,avg_metacatintron,metagene_breaks_lim=[-1000,1000])
        if options.bed:
            rscript+=draw_siteprofiles(sitebreaks,avg_siteprof)
        rscript+=inout.draw_exon_intron_profile_plot(metaexon_breaks,avg_metaexon,metaintron_breaks,avg_metaintron, options.elowers, options.euppers, options.ilowers, options.iuppers)
        
        
    ofhd.write(rscript)    # write wig profiling
    ofhd.write(R.devoff())
    ofhd.close()
    jobcount+=1
    
    info ('#... cong! Run R on %s!' % (options.name+'.R'))


     
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0)
