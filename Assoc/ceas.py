#!/usr/bin/env python

"""Module Description

Copyright (c) 2008 H. Gene Shin <shin@jimmy.harvard.edu>

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
from optparse import OptionParser
from Cistrome.Assoc import *
from Cistrome.Assoc.inout import MYSQL

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
# Main function
# ------------------------------------
def main():
    
    # read the options and validate them
    options=opt_validate(prepare_optparser())

    # CEAS run
    # read the gene annotation table
    jobcount=1
    info("#%d read the gene table..." %jobcount)
        
    # read
    GeneT=inout.GeneTable()
    GeneT.read(Host=options.Host,User=options.User,Db=options.Db,annotation='refGene',which=options.which)    
    GeneT.sort()
    chroms_GeneT=GeneT.get_chroms()
    jobcount+=1

    # read ChIP regions
    info("#%d read the bed file of ChIP regions..." %jobcount)
    Cbed=inout.Bed()
    Cbed.read(options.bed)
    Csampler=sampler.ChIPSampler()
    ChIP=Csampler.sample(Cbed,resolution=options.chip_res)
    del Cbed
    jobcount+=1
    
    # read regions of interest if it is given
    if options.ebed:
        info("#%d read the bed file of regions of interest..." %jobcount)
        roi=inout.Bed()
        roi.read(options.ebed)
        jobcount+=1
    else: roi=None

    # if wig profiling is not being run.
    if not options.bg:
        
        # iterate through chromosomes of the gene table
        info("#%d read the pre-computed genome bg annotation..." %jobcount)
        GenomeBGS=tables.SummaryGBG(name='GenomeBGS')
        GenomeBGS.readdb(Db=options.gdb)
        GP=_interpoloate_gbg(gdb,options.promoter,options.bipromoter,options.downstream)
        chroms_bg=GP.get_chroms()
        
        # if any regions of interest are given
        if options.ebed:
            GP=_get_bgroi(GP,GenomeBGS,roi=roi,bg_res=options.bg_res)
        
        # annotate ChIP regions
        info('#%d annotate the ChIP regions...' %jobcount)
        Annot=annotator.Annotator()
        ChIPA=Annot.annotate(genome_coordinates=ChIP,gene_table=GeneT,roi=roi,prom=options.promoter,biprom=options.bipromoter,down=options.downstream,gene_div=(3,5))
        CS,CP=Annot.summarize(ChIPA)
        # make the table complete with missing chromsomes, if there are
        annotator.make_table_complete(CS,chroms_bg)
        annotator.make_table_complete(CP,chroms_bg)
        # get the pvalues
        CPval=annotator.estimate_pvals(GP,CS,CP)
        jobcount+=1

        # open outfile 
        info('#%d write a R script of CEAS...' %jobcount)
        ofhd=open(options.name+'.R','w')
        pdfname=options.name+'.pdf'
        # the first part of CEAS R script. Because wig profiling is not run, just terminate
        rscript=R.pdf(pdfname,height=11.5,width=8.5)   
        rscript+=inout.draw_CEAS(GP,CP,CPval,bg_res=options.bg_res,chip_res=options.chip_res,prom=options.promoter,biprom=options.bipromoter,down=options.downstream,gene_div=(3,5))      
        ofhd.write(rscript)    # write CEAS
    
    # when wig profiling is running
    if options.pf:
        
        if options.bg:
            GenomeBGS=tables.Summary()
        
        # if gene groups are give
        if options.gn_groups:
            subsets=inout.read_gene_subsets(options.gn_groups)
        
        chrom=''
        chrcount=1
        prof=profiler.WigProfiler()
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
            if re.search(r'chrom=(\w+)\s',line):
                newchrom=re.search(r'chrom=(\w+)\s',line).group(1)
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
                    info("#%d-%d work on %s..." %(jobcount,chrcount,chrom))
                    input=inout.Wig()
                    input.add_line(chrom,l)
                    chrcount+=1
            elif chrom!='' and chrom!=newchrom:    # new chromosome
                if chrom in chroms_GeneT:
                    # do genome BG annotation
                    if options.bg:
                        Sampler=sampler.GenomeSampler()
                        Annotator=annotator.Annotator()
                        GA=Annotator.annotate(Sampler.sample(input,resolution=options.bg_res),GeneT,roi=roi,prom=options.promoter,biprom=options.bipromoter,down=options.downstream,gene_div=(3,5))
                        tempS,tempP=Annotator.summarize(GA)
                        GenomeBGS.add_row(chrom,tempS.get_row(chrom))
                                                
                    # wig profiling
                    names,breaks,upstreams,downstreams,metagene_breaks,metagenes,metaexon_breaks,metaexons,metaintron_breaks,metaintrons=prof.profile(input,GeneT,rel_pos=options.rel_dist,metagenesize=options.metagene_size,step=options.pf_res,which=options.which,exonratio=0.5,emask=options.emask,imask=options.imask)
            
                    # get average of this chromosome
                    avg_up,upcount=corelib.mean_col_by_col(upstreams,counts=True)
                    avg_down,downcount=corelib.mean_col_by_col(downstreams,counts=True)
                    avg_mg,genecount=corelib.mean_col_by_col(metagenes,counts=True)
                    avg_me,exoncount=corelib.mean_col_by_col(metaexons,counts=True)
                    avg_mi,introncount=corelib.mean_col_by_col(metaintrons,counts=True)
            
                    if not FIRST:    # if not first chromosome
                        avg_upstream,avg_upcount=corelib.weight_mean_col_by_col([avg_upstream,avg_up],[avg_upcount,upcount],counts=True)
                        avg_downstream,avg_downcount=corelib.weight_mean_col_by_col([avg_downstream,avg_down],[avg_downcount,upcount],counts=True)
                        avg_metagene,avg_genecount=corelib.weight_mean_col_by_col([avg_metagene,avg_mg],[avg_genecount,genecount],counts=True)
                        avg_metaexon,avg_exoncount=corelib.weight_mean_col_by_col([avg_metaexon,avg_me],[avg_exoncount,exoncount],counts=True)
                        avg_metaintron,avg_introncount=corelib.weight_mean_col_by_col([avg_metaintron,avg_mi],[avg_introncount,introncount],counts=True)
                        del avg_up,avg_down,avg_mg,avg_me,avg_mi,upcount,downcount,genecount,exoncount,introncount
                
                        if options.gn_groups:    # when gene sub-gropus are given
                            ixs,subsets=profiler.get_gene_indicies(names,subsets)
                            avg_ups,upcs,avg_downs,downcs,avg_mgs,gcs,avg_mes,ecs,avg_mis,ics=profiler.select_profiles_chr_by_chr(ixs,upstreams,downstreams,metagenes,metaexons,metaintrons)
                            avg_upstreams,avg_upcounts=profiler.weight_mean_profiles_chr_by_chr(avg_upstreams,avg_upcounts,avg_ups,upcs)
                            avg_downstreams,avg_downcounts=profiler.weight_mean_profiles_chr_by_chr(avg_downstreams,avg_downcounts,avg_downs,downcs)
                            avg_metagenes,avg_genecounts=profiler.weight_mean_profiles_chr_by_chr(avg_metagenes,avg_genecounts,avg_mgs,gcs)
                            avg_metaexons,avg_exoncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metaexons,avg_exoncounts,avg_mes,ecs)
                            avg_metaintrons,avg_introncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metaintrons,avg_introncounts,avg_mis,ics)
                            del avg_ups,avg_downs,avg_mgs,avg_mes,avg_mis,upcs,downcs,gcs,ecs,ics
                
                    else:   # if first chromosome
                        avg_upstream=avg_up
                        avg_downstream=avg_down
                        avg_metagene=avg_mg
                        avg_metaexon=avg_me
                        avg_metaintron=avg_mi
                        avg_upcount=upcount
                        avg_downcount=downcount
                        avg_genecount=genecount
                        avg_exoncount=exoncount
                        avg_introncount=introncount
            
                        if options.gn_groups:
                            ixs,subsets=profiler.get_gene_indicies(names,subsets)
                            avg_upstreams,avg_upcounts,avg_downstreams,avg_downcounts,avg_metagenes,avg_genecounts,avg_metaexons,avg_exoncounts,avg_metaintrons,avg_introncounts=profiler.select_profiles_chr_by_chr(ixs,upstreams,downstreams,metagenes,metaexons,metaintrons)
                        FIRST=False
                
                    del upstreams,downstreams,metagenes,metaexons,metaintrons    
                
                # set chrom to the new chromosome
                chrom=newchrom
                if chrom in chroms_GeneT:    # only if the new chromosome is in the chroms of gene table, a wig object is initiated.
                    info("#%d-%d work on %s..." %(jobcount,chrcount,chrom))
                    input=inout.Wig()
                    input.add_line(chrom,l)
                    chrcount+=1
            else:    # in the middle of chromosome
                if chrom in chroms_GeneT:   # only if the new chromosome is in the chroms of gene table, the wig object is updated.
                    input.add_line(chrom,l)
                        
    # do profiling for the last chromosome 
        if chrom in chroms_GeneT:
            
            if options.bg:
                Sampler=sampler.GenomeSampler()
                Annotator=annotator.Annotator()
                GA=Annotator.annotate(Sampler.sample(input,resolution=options.bg_res),GeneT,roi=roi,prom=options.promoter,biprom=options.bipromoter,down=options.downstream,gene_div=(3,5))
                tempS,tempP=Annotator.summarize(GA)
                GenomeBGS.add_row(chrom,tempS.get_row(chrom))
                GenomeBGS.summarize()
                GP=GenomeBGS.get_p()
                
                if options.ebed:
                    GP=_get_bgroi(GP,GenomeBGS,roi=roi,bg_res=options.bg_res)
            
            # profiling
            names,breaks,upstreams,downstreams,metagene_breaks,metagenes,metaexon_breaks,metaexons,metaintron_breaks,metaintrons=prof.profile(input,GeneT,rel_pos=options.rel_dist,metagenesize=options.metagene_size,step=options.pf_res,which=options.which,exonratio=0.5,emask=options.emask,imask=options.imask)
            del input 
            # get average of this chromosome
            avg_up,upcount=corelib.mean_col_by_col(upstreams,counts=True)
            avg_down,downcount=corelib.mean_col_by_col(downstreams,counts=True)
            avg_mg,genecount=corelib.mean_col_by_col(metagenes,counts=True)
            avg_me,exoncount=corelib.mean_col_by_col(metaexons,counts=True)
            avg_mi,introncount=corelib.mean_col_by_col(metaintrons,counts=True)
            
            if not FIRST:    # the first chromosome profiling
                avg_upstream,avg_upcount=corelib.weight_mean_col_by_col([avg_upstream,avg_up],[avg_upcount,upcount],counts=True)
                avg_downstream,avg_downcount=corelib.weight_mean_col_by_col([avg_downstream,avg_down],[avg_downcount,upcount],counts=True)
                avg_metagene,avg_genecount=corelib.weight_mean_col_by_col([avg_metagene,avg_mg],[avg_genecount,genecount],counts=True)
                avg_metaexon,avg_exoncount=corelib.weight_mean_col_by_col([avg_metaexon,avg_me],[avg_exoncount,exoncount],counts=True)
                avg_metaintron,avg_introncount=corelib.weight_mean_col_by_col([avg_metaintron,avg_mi],[avg_introncount,introncount],counts=True)
                del avg_up,avg_down,avg_mg,avg_me,avg_mi,upcount,downcount,genecount,exoncount,introncount
        
                if options.gn_groups:
                    ixs,subsets=profiler.get_gene_indicies(names,subsets)
                    avg_ups,upcs,avg_downs,downcs,avg_mgs,gcs,avg_mes,ecs,avg_mis,ics=profiler.select_profiles_chr_by_chr(ixs,upstreams,downstreams,metagenes,metaexons,metaintrons)
                    avg_upstreams,avg_upcounts=profiler.weight_mean_profiles_chr_by_chr(avg_upstreams,avg_upcounts,avg_ups,upcs)
                    avg_downstreams,avg_downcounts=profiler.weight_mean_profiles_chr_by_chr(avg_downstreams,avg_downcounts,avg_downs,downcs)
                    avg_metagenes,avg_genecounts=profiler.weight_mean_profiles_chr_by_chr(avg_metagenes,avg_genecounts,avg_mgs,gcs)
                    avg_metaexons,avg_exoncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metaexons,avg_exoncounts,avg_mes,ecs)
                    avg_metaintrons,avg_introncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metaintrons,avg_introncounts,avg_mis,ics)
                    del avg_ups,avg_downs,avg_mgs,avg_mes,avg_mis,upcs,downcs,gcs,ecs,ics
            else:
                avg_upstream=avg_up
                avg_downstream=avg_down
                avg_metagene=avg_mg
                avg_metaexon=avg_me
                avg_metaintron=avg_mi
                avg_upcount=upcount
                avg_downcount=downcount
                avg_genecount=genecount
                avg_exoncount=exoncount
                avg_introncount=introncount
                if options.gn_groups:
                    ixs,subsets=profiler.get_gene_indicies(names,subsets)
                    avg_upstreams,avg_upcounts,avg_downstreams,avg_downcounts,avg_metagenes,avg_genecounts,avg_metaexons,avg_exoncounts,avg_metaintrons,avg_introncounts=profiler.select_profiles_chr_by_chr(ixs,upstreams,downstreams,metagenes,metaexons,metaintrons)
    
            del upstreams,downstreams,metagenes,metaexons,metaintrons
        jobcount+=1
        
        if options.bg:
            info('#%d annotate ChIP regions...' %jobcount)
            Annot=annotator.Annotator()
            ChIPA=Annot.annotate(genome_coordinates=ChIP,gene_table=GeneT,roi=roi,prom=options.promoter,biprom=options.bipromoter,down=options.downstream,gene_div=(3,5))
            CS,CP=Annot.summarize(ChIPA)
            CPval=annotator.estimate_pvals(GP,CS,CP)
            jobcount+=1    
            
            info('#%d write R script of CEAS and wig profiling...' %jobcount)
            ofhd=open(options.name+'.R','w')
            pdfname=options.name+'.pdf'
            # the first part of CEAS R script. Because wig profiling is not run, just terminate
            rscript=R.pdf(pdfname,height=11.5,width=8.5)
            rscript+=inout.draw_CEAS(GP,CP,CPval,bg_res=options.bg_res,chip_res=options.chip_res,prom=options.promoter,biprom=options.bipromoter,down=options.downstream,gene_div=(3,5)) 
            ofhd.write(rscript)    # writing CEAS
        else:
            info('#%d append R script of wig profiling...' %jobcount)
        
        # write R script
        if options.gn_groups:
            # append the profiles of all genes
            avg_upstreams.append(avg_upstream)
            avg_downstreams.append(avg_downstream)
            avg_metagenes.append(avg_metagene)
            avg_metaexons.append(avg_metaexon)
            avg_metaintrons.append(avg_metaintron)  
                
            rscript=inout.draw_profile_plots(breaks,avg_upstreams,avg_downstreams,metagene_breaks,avg_metagenes,metaexon_breaks,avg_metaexons,metaintron_breaks,avg_metaintrons,metagene_breaks_lim=[-1000,1000],legends=options.gn_names)
        else:              
            rscript=inout.draw_profile_plot(breaks,avg_upstream,avg_downstream,metagene_breaks,avg_metagene,metaexon_breaks,avg_metaexon,metaintron_breaks,avg_metaintron,metagene_breaks_lim=[-1000,1000])
        ofhd.write(rscript)    # write wig profiling
        
    ofhd.write(R.devoff())   
    ofhd.close()
    
    info ('#... cong! Run R on %s!' %(options.name+'.R'))


# ------------------------------------
# functions
# ------------------------------------
  
def prepare_optparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    
    """
    
    usage = "usage: %prog <-b bed -g gdb> [options]"
    description = "CEAS -- Cis-regulatory Element Annotation System"
    
    optparser = OptionParser(version="%prog 1.0",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-b","--bed",dest="bed",type="string",
                         help="BED file of ChIP regions.")
    optparser.add_option("-w","--wig",dest="wig",type="string",
                         help="WIGGLE file for either wig profiling or genome background re-annotation. WARNING: A WIGGLE file must be given for wig profiling.")
    optparser.add_option("-e","--ebed",dest="ebed",type="string",
                         help="BED file of extra regions of interest (eg, non-coding regions)")
    optparser.add_option("-g","--gene-db",dest="gdb",type="string",
                         help="Gene annotation table (a local sqlite3 db file provided by CEAS or species name in UCSC). CEAS searches the designated directory for the the db file. If not find, CEAS looks up UCSC for the table. WARNING: When using UCSC, MySQLdb package must be installed.")
    optparser.add_option("--bg",action="store_true",dest="bg",\
                         help="Run genome BG annotation. WARNING: This flag is effective only if a wig file is given through -w (--wig). Otherwise, ignored.",default=False)
    optparser.add_option("--name",dest="name",\
                         help="Experiment name. This will be used to name the output file. If an experiment name is not given, input BED file name will be used instead.")      
    optparser.add_option("--chip-res",dest="chip_res",type="int",
                         help="ChIP annotation resolution, DEFAULT: 600bp. WARNING: Value less than 600bp turns to be 600 p", default=600)    
    optparser.add_option("--promoter",dest="promoter",type="int",
                         help="Promoter size for annotation, DEFAULT: 3000bp", default=3000)    
    optparser.add_option("--bipromoter",dest="bipromoter",type="int",
                         help="Bidirectional-promoter size for annotation, DEFAULT: 5000bp", default=5000)  
    optparser.add_option("--downstream",dest="downstream",type="int",
                         help="Downstream size for annotation, DEFAULT: 3000bp", default=3000)     
    optparser.add_option("--pf-res", dest="pf_res", type="int",\
                          help="Wig profiling resolution, DEFAULT: 50bp. WARNING: Value smaller than the wig step (resolution) may cause aliasing error.", default=50) 
    optparser.add_option("--rel-dist",dest="rel_dist",type="int",
                         help="Relative distance to TSS/TTS in wig profiling, DEFAULT: 3000bp", default=3000)
    optparser.add_option("--metagene-size",dest="metagene_size",type="int",
                         help="Normalized gene length in wig profiling, DEFAULT: 3000bp. Every gene is normalized to have this length.", default=3000)   
    optparser.add_option("--gn-groups",dest="gn_groups",type="string",\
                         help="Gene-groups of particular interest in wig profiling. Each gene group file must have gene names in the 1st column. The file names are separated by commas w/ no space (eg, --gn-groups=top10.txt,bottom10.txt)") 
    optparser.add_option("--gn-group-names", dest="gn_names",type="string",\
                         help="The names of the gene groups in --gn-groups. The gene group names are separated by commas. (eg, --gn-group-names='top 10%,bottom 10%'). These group names appear in the legends of the wig profiling plots. If no group names given, the groups are represented as 'Group 1, Group2,...Group n'.")
    optparser.add_option("--alt-gn", action="store_true", dest="name2",\
                         help="Whether alternative gene names (eg, 'name2' in refGene of UCSC) are used in --gn-groups or not. This flag is meaningful only if --gn-groups is set.",default=False)
    optparser.add_option("--verbose",dest="verbose",type="string",
                         help="Name of a directory. if set, save verbose information in the direcotry.")
    
    return optparser


def opt_validate (optparser):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    
    # input BED file and GDB must be given 
    if not (options.bed and options.gdb):
        optparser.print_help()
        sys.exit(1)

    # get gdb lower case
    HAVELOCALGDB=os.path.isfile(options.gdb)
    if HAVELOCALGDB:
        options.Host=None
        options.User=None
        options.Db=options.gdb.lower()
    else:
        if MYSQL:
            options.Host="genome-mysql.cse.ucsc.edu"
            options.User="genome"
            options.Db=os.path.split(options.gdb)[-1].lower()
        else:
            error('MySQLdb package needs to be installed to use UCSC or a local sqlite3 db file must exist.')
            error('Check -g (--gdb). No such file or species: %s' %options.gdb)
            sys.exit(1)
        
    
    # bg background resolution is set to 100    
    options.bg_res=100
    
    # wig file
    if options.wig:
        HAVEWIG=os.path.isfile(options.wig)
    else: HAVEWIG=False
        
    REBG=False
    if HAVEWIG and options.bg:
        REBG=True
    elif not HAVELOCALGDB and HAVEWIG and not options.bg:
        error('Genome BG annotation must be run when no pre-computed BG annotation exists. Set --bg.')
        sys.exit(1)
    elif not HAVELOCALGDB and not HAVEWIG and options.bg:
        if not options.wig:
            error('Genome BG annotation must be run when no pre-computed BG annotation exists. Give a WIGGLE file through -w (--wig)')
        else:
            error('Genome BG annotation must be run when no pre-computed BG annotation exists. Check -w (--wig). No such file: %s' %options.wig)
        sys.exit(1)
    elif not HAVELOCALGDB and not HAVEWIG and not options.bg:
        if not options.wig:
            error('Genome BG annotation must be run when no pre-computed BG annotation exists. set --bg and give a WIGGLE file through -w (--wig)')
        else:
            error('Genome BG annotation must be run when no pre-computed BG annotation exists. Set --bg and check -w (--wig). No such file: %s' %options.wig)
        sys.exit(1)
    options.bg=REBG
    
    # if a WIGGLE file given, do profiling
    if HAVEWIG:
        options.pf=True
    else: options.pf=False
        
    # non-coding regions
    if options.ebed:
        if not os.path.isfile(options.ebed):
            error('Check -e (--ebed). No such file: %s' %options.ebed)
            sys.exit(1)
    
     # get namename
    if not options.name:
        options.name=os.path.split(options.bed)[-1].rsplit('.bed',2)[0]
    
    # the minimum ChIP annotation resolution is 600
    options.chip_res=max(600,options.chip_res)
    
    #check if name2 is going to be used instead of name
    if options.name2:
        options.which='name2'
    else:
        options.which='name'
                    
    # check the gene group files    
    if options.pf and options.gn_groups:
        parsed=options.gn_groups.split(',')
        for p in parsed:
            if not os.path.isfile(p):
                error('Check --gn-groups. No such file: %s' %p)
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
    
    # emask and imask -- potentially can be added later
    options.emask=0
    options.imask=0
    
    return options

def _get_bgroi(GenomeBGP,GenomeBGS,roi,bg_res=100):
    """Get the background annotation for regions of interest given through -e (or --ebed) option
    
    Parameters:
    1. GenomeGBP: a P object (see inout.py) of genome background annotation. This will be modified by this function and returned.
    2. GenomeGBS: a SummaryGBG object (see inout.py) of genome background annotation
    2. roi: a Bed object of regions of interest
    3. bg_res: genome background annotation resolution (default=100bp)
    
    """
    
    # sampler
    Sampler=sampler.ChIPSampler()
    roisamp=Sampler.sample(bed=roi,resolution=bg_res)
    
    chroms=set(GenomeBGP.get_chroms()).intersection(roi.keys())
    bgroi={}
    whole=0
    for chrom in chroms:
        num_this_chr=len(roi[chrom])
        bgroi[chrom]=num_this_chr
        whole+=num_this_chr
    bgroi['whole']=whole
    
    for chrom in bgroi.keys():
        try:
            GenomeBGP[chrom]['roi']=1.0*bgroi[chrom]/GenomeBGS[chrom]['Ns']
        except ZeroDivisionError:
            pass
        except KeyError:
            pass
    
    return GenomeBGP

def _interpoloate_gbg(gdb,promoter,bipromoter,downstream):
    """In using the pre-computed genome bg model, this function performs linear interpolation of 
    genome-wide enrichments of promoter, bidirectional promoter, and downstream.
    
    Parameters:
    1. gdb: sqlite3 db file. This file must have GenomeBGS and GenomeBGP tables
    2. promoter: promoter length given through options.promoter
    3. bipromoter: bidirectional promoter length given through options.bipromoter
    4. downstream: downstream length given through options.downstream
    
    Return
    GP: a P object (see tables.py). This object contains genome bg annotation
    
    """
    
    GenomeBGP=tables.PGBG(name='GenomeBGP',numprom=11,numbiprom=21,numdown=11)
    GenomeBGP.readdb(Db=gdb)
    
    GP=tables.P(name='GP')
    # the given promoter, bipromoter, and downstream lengths
    new['promoter'] = [promoter/3, 2*promoter/3, promoter]
    new['bipromoter'] = [bipromoeter/2, bipromoter]
    new['downstream'] = [downstream/3, 2*downstream/3, downstream]
    
    # the model promoter, bipromoter, and downstream lengths
    mod['promoter']=[0, 500]+seq(fr=1000,to=10000,by=1000)
    mod['bipromoter']=[0, 500]+seq(fr=1000,to=20000,by=1000)
    mod['downstream']=[0, 500]+seq(fr=1000,to=10000,by=1000)
    
    for chrom in GenomeBGP.get_chroms():
        GP.init_table(chrom)
        for column in GenomeBGP.columns[1:]:
            if column!='promoter' and column!='bipromoter' and column!='downstream':
                GP[chrom][column]=GenomeBGP[chrom][column]
            else:
                vals=[0.0]+GenomeBGP[chrom][column]
                interpol=[]
                for x in new[column]:
                    i=corelib.findbin(x, mod[column])
                    interpol.append(corelib.lininterpol([modproms[i],vals[i]], [modproms[i+1],vals[i+1]],p))
                
                GP[chrom][column]=interpol
                
    return GP
                    
                    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0)
