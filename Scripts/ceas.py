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
    if options.gdbtype == 'ucsc' or options.gdbtype == 'localdb':
        GeneT = inout.GeneTable()
        GeneT.read(Host = options.Host, User= options.User, Db=options.gdb, annotation='refGene', which=options.which)
    elif options.gdbtype == 'localbed':
        GBed = inout.Bed()
        GBed.read(options.gdb)
        GeneT = convert_BED2GeneTable(GBed)
        del GBed
    GeneT.sort()
    chroms_GeneT=GeneT.get_chroms()
    
    # determine the metagenesize, concated exon size, concatenated intron size
    if options.metagene_size:
        options.catexon_size, options.catintron_size = corelib.find_nearest_multiple(options.metagene_size/2, options.pf_res), corelib.find_nearest_multiple(options.metagene_size/2, options.pf_res)
    else:    
        options.metagene_size = return_med_gene(GeneT, options.pf_res)
        options.catexon_size, options.catintron_size = return_med_catexon_catintron(GeneT, options.pf_res)
    
    # get the exon and intron sizes to consider in the average profiling
    if options.exin_pf:     # only when exon and intron profiling is on
        options = determine_exon_intron_sizes(GeneT, options)    
    jobcount+=1

    if options.ceas:
        # read ChIP regions
        info("#%d read the bed file of ChIP regions..." %jobcount)
        ChIP=inout.Bed()
        ChIP.read(options.bed)
        ChIP.sort()
        Csampler=sampler.ChIPSampler()
        ChIPSamp=Csampler.sample(ChIP,resolution=options.chip_res)
        jobcount+=1
    
        # do gene-centered annotation
        info('#%d perform gene-centered annotation...' %jobcount)
        GAnnotator=annotator.GeneAnnotator()
        GAnnotator.annotate(GeneT, ChIP, u=options.span, d=options.span)
        GAnnotator.map.set_name(options.name)
        GAnnotator.map.write()
        jobcount+=1
    
    # read regions of interest if it is given
    if options.ceas and options.ebed:
        info("#%d read the bed file of regions of interest..." %jobcount)
        roi=inout.Bed()
        roi.read(options.ebed)
        jobcount+=1
    else: roi=None

    # if background annotation is not being run.
    if options.ceas and not options.bg:
        
        # iterate through chromosomes of the gene table
        info("#%d read the pre-computed genome bg annotation..." %jobcount)
        GenomeBGS=tables.SummaryGBG(name='GenomeBGS')
        GenomeBGS.readdb(Db=options.gdb)
        GP=_interpoloate_gbg(options.gdb,options.promoter,options.bipromoter,options.downstream)
        chroms_bg=GP.get_chroms()
        jobcount+=1
        
        # if any regions of interest are given
        if options.ebed:
            GP=_get_bgroi(GP,GenomeBGS,roi=roi,bg_res=options.bg_res)
        
        # annotate ChIP regions
        info('#%d perform ChIP region annotation...' %jobcount)
        Annot=annotator.Annotator()
        ChIPA=Annot.annotate(genome_coordinates=ChIPSamp,gene_table=GeneT,roi=roi,prom=options.promoter,biprom=options.bipromoter,down=options.downstream,gene_div=(3,5))
        CS,CP=Annot.summarize(ChIPA)
        CES, CEP = Annot.obtain_distribution_of_sites(ChIPA)
        # make the table complete with missing chromsomes, if there are
        annotator.make_table_complete(CS,chroms_bg)
        annotator.make_table_complete(CP,chroms_bg)
        # get the pvalues
        CPval=annotator.estimate_pvals(GP,CS,CP)
        jobcount+=1
        
        # open outfile 
        info('#%d write a R script of ChIP region annotation...' %jobcount)
        ofhd=open(options.name+'.R','w')
        pdfname=options.name+'.pdf'
        # the first part of CEAS R script. Because wig profiling is not run, just terminate
        rscript=R.pdf(pdfname,height=11.5,width=8.5)   
        rscript+=inout.draw_CEAS(GP,CP,CPval,bg_res=options.bg_res,chip_res=options.chip_res,prom=options.promoter,biprom=options.bipromoter,down=options.downstream,gene_div=(3,5))      
        ofhd.write(rscript)    # write CEAS
        rscript = inout.draw_pie_distribution_of_elements(CEP, prom=options.promoter, down=options.downstream)
        ofhd.write(rscript)
        jobcount+=1
        
    # when wig profiling is running
    if options.pf:
        
        if options.bg:
            GenomeBGS=tables.Summary()
        
        # if gene groups are give
        if options.gn_groups:
            subsets=inout.read_gene_subsets(options.gn_groups)
        
        chrom=''
        chrcount=1
        prof=profiler.WigProfiler2(rel_dist=options.rel_dist, step=options.pf_res, metagenesize=options.metagene_size, \
                                   catexonsize=options.catexon_size, catintronsize=options.catintron_size, metaexonsizes=options.metaexonsizes, \
                                   metaintronsizes=options.metaintronsizes, elowers=options.elowers, euppers=options.euppers, ilowers=options.ilowers, \
                                   iuppers=options.iuppers)
        wigsize = {}
        swig=inout.Wig()
        ws = sampler.WigSamplerFast()
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
                    if options.bg:
                        info("#%d-%d run wig profiling and genome bg annotation of %s..." %(jobcount,chrcount,chrom))
                    else:
                        info("#%d-%d run wig profiling of %s..." %(jobcount,chrcount,chrom))
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
                    
                    # update ChIP regions with WIG scores
                    profiler.scan_scores_in_wig(ChIP, input)                    
                    #swig += tempswig
                    #del tempswig
#                    ws=sampler.WigSampler()
#                    resolution=input.chrom_length(chrom)/options.nwigpts
#                    swig+=sampler.fillupwig(ws.sample(input,resolution=resolution),resolution=resolution)
                    
                    try:
                        wigsize[chrom] = (input[chrom][0][0], input[chrom][0][-1])
                    except IndexError:
                        wigsize[chrom] = (0, 0)
                                                                 
                    # wig profiling
                    profiles = prof.profile(input, GeneT)
        
                    # get average of this chromosome
                    avg_up,upcount=corelib.mean_col_by_col(profiles['upstreams'], counts=True)
                    avg_down,downcount=corelib.mean_col_by_col(profiles['downstreams'], counts=True)
                    avg_mg,genecount=corelib.mean_col_by_col(profiles['genes'], counts=True)
                    avg_mce,cexoncount=corelib.mean_col_by_col(profiles['catexons'], counts=True)
                    avg_mci,cintroncount=corelib.mean_col_by_col(profiles['catintrons'],counts=True)
                    if options.exin_pf:
                        n_metas = len(profiles['exons'])
                        outs = map(corelib.mean_col_by_col, map(corelib.extend_list_series, profiles['exons']), [True] * len(profiles['exons']))
                        avg_me = map(operator.itemgetter(0), outs)
                        exoncount = map(operator.itemgetter(1), outs)
                        outs = map(corelib.mean_col_by_col, map(corelib.extend_list_series, profiles['introns']), [True] * len(profiles['introns']))
                        avg_mi = map(operator.itemgetter(0), outs)
                        introncount = map(operator.itemgetter(1), outs)
                        del outs
            
                    if not FIRST:    # if not first chromosome
                        avg_upstream,avg_upcount=corelib.weight_mean_col_by_col([avg_upstream,avg_up],[avg_upcount,upcount],counts=True)
                        avg_downstream,avg_downcount=corelib.weight_mean_col_by_col([avg_downstream,avg_down],[avg_downcount,upcount],counts=True)
                        avg_metagene,avg_genecount=corelib.weight_mean_col_by_col([avg_metagene,avg_mg],[avg_genecount,genecount],counts=True)
                        avg_metacatexon,avg_cexoncount=corelib.weight_mean_col_by_col([avg_metacatexon,avg_mce],[avg_cexoncount,cexoncount],counts=True)
                        avg_metacatintron,avg_cintroncount=corelib.weight_mean_col_by_col([avg_metacatintron,avg_mci],[avg_cintroncount,cintroncount],counts=True)
                        if options.exin_pf:
                            outs = map(corelib.weight_mean_col_by_col, map(lambda x, y: [x, y], avg_metaexon, avg_me),  map(lambda x, y: [x, y], avg_exoncount, exoncount), [True] * n_metas)
                            avg_metaexon = map(operator.itemgetter(0), outs)
                            avg_exoncount = map(operator.itemgetter(1), outs)
                            outs = map(corelib.weight_mean_col_by_col, map(lambda x, y: [x, y], avg_metaintron, avg_mi), map(lambda x, y: [x, y], avg_introncount, introncount), [True] * n_metas)
                            avg_metaintron = map(operator.itemgetter(0), outs)
                            avg_introncount = map(operator.itemgetter(1), outs)
                            del outs, avg_me, avg_mi, exoncount, introncount
                        del avg_up,avg_down,avg_mg,avg_mce,avg_mci,upcount,downcount,genecount,cexoncount,cintroncount
                
                        if options.gn_groups:    # when gene sub-gropus are given
                            if options.which=='name2':
                                ixs, subsets = profiler.get_gene_indicies(profiles['name2'], subsets)
                            else:
                                ixs, subsets = profiler.get_gene_indicies(profiles['name'], subsets)
            
                            avg_ups, upcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['upstreams'])
                            avg_downs, downcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['downstreams'])
                            avg_mgs, gcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['genes'])
                            avg_mces, cecs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catexons'])
                            avg_mcis, cics = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catintrons'])
                            if options.exin_pf:
                                avg_mes, ecs = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['exons'])
                                avg_mis, ics = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['introns'])
                        
                            avg_upstreams,avg_upcounts=profiler.weight_mean_profiles_chr_by_chr(avg_upstreams,avg_upcounts,avg_ups,upcs)
                            avg_downstreams,avg_downcounts=profiler.weight_mean_profiles_chr_by_chr(avg_downstreams,avg_downcounts,avg_downs,downcs)
                            avg_metagenes,avg_genecounts=profiler.weight_mean_profiles_chr_by_chr(avg_metagenes,avg_genecounts,avg_mgs,gcs)
                            avg_metacatexons,avg_cexoncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metacatexons,avg_cexoncounts,avg_mces,cecs)
                            avg_metacatintrons,avg_cintroncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metacatintrons,avg_cintroncounts,avg_mcis,cics)
                            if options.exin_pf:
                                avg_metaexons, avg_exoncounts = profiler.weight_mean_profiles_chr_by_chr_meta(avg_metaexons, avg_exoncounts, avg_mes, ecs)
                                avg_metaintrons, avg_introncounts = profiler.weight_mean_profiles_chr_by_chr_meta(avg_metaintrons, avg_introncounts, avg_mis, ics)
                                del avg_mes,avg_mis, ecs,ics
                            del avg_ups,avg_downs,avg_mgs,avg_mces,avg_mcis,upcs,downcs,gcs,cecs,cics
                
                    else:   # if first chromosome
                        avg_upstream = avg_up
                        avg_downstream = avg_down
                        avg_metagene = avg_mg
                        avg_metacatexon = avg_mce
                        avg_metacatintron = avg_mci
                        avg_upcount = upcount
                        avg_downcount = downcount
                        avg_genecount = genecount
                        avg_cexoncount = cexoncount
                        avg_cintroncount = cintroncount
                        if options.exin_pf:
                            avg_metaexon = avg_me
                            avg_metaintron = avg_mi
                            avg_exoncount=exoncount
                            avg_introncount=introncount
                       
                        if options.gn_groups:
                            if options.which=='name2':
                                ixs, subsets = profiler.get_gene_indicies(profiles['name2'], subsets)
                            else:
                                ixs, subsets = profiler.get_gene_indicies(profiles['name'], subsets)

                            avg_upstreams, avg_upcounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['upstreams'])
                            avg_downstreams, avg_downcounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['downstreams'])
                            avg_metagenes, avg_genecounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['genes'])
                            avg_metacatexons, avg_cexoncounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catexons'])
                            avg_metacatintrons, avg_cintroncounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catintrons'])
                            if options.exin_pf:
                                avg_metaexons, avg_exoncounts = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['exons'])
                                avg_metaintrons, avg_introncounts = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['introns'])
                        
                        FIRST=False
                
                    del profiles  
                
                # set chrom to the new chromosome
                chrom=newchrom
                if chrom in chroms_GeneT:    # only if the new chromosome is in the chroms of gene table, a wig object is initiated.
                    if options.bg:
                        info("#%d-%d run wig profiling and genome bg annotation of %s..." %(jobcount,chrcount,chrom))
                    else:
                        info("#%d-%d run wig profiling of %s..." %(jobcount,chrcount,chrom))
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
                
                # if extra bed file exists, get the bg statistics of the regions
                if options.ceas and options.ebed:
                    GP=_get_bgroi(GP,GenomeBGS,roi=roi,bg_res=options.bg_res)
            
            # update ChIP regions with wig scores
            profiler.scan_scores_in_wig(ChIP, input)
            #swig += tempswig
            #del tempswig
#            ws=sampler.WigSampler()
#            resolution=input.chrom_length(chrom)/options.nwigpts
#            swig+=sampler.fillupwig(ws.sample(input,resolution=resolution),resolution=resolution)
            
            try:
                wigsize[chrom] = (input[chrom][0][0], input[chrom][0][-1])
            except IndexError:
                wigsize[chrom] = (0, 0)
                                
            # profiling
            profiles = prof.profile(input, GeneT)
            del input 
            
            # get average of this chromosome
            avg_up,upcount = corelib.mean_col_by_col(profiles['upstreams'], counts=True)
            avg_down,downcount = corelib.mean_col_by_col(profiles['downstreams'], counts=True)
            avg_mg,genecount = corelib.mean_col_by_col(profiles['genes'], counts=True)
            avg_mce,cexoncount = corelib.mean_col_by_col(profiles['catexons'], counts=True)
            avg_mci,cintroncount = corelib.mean_col_by_col(profiles['catintrons'],counts=True)
            if options.exin_pf:
                n_metas = len(profiles['exons'])  
                outs = map(corelib.mean_col_by_col, map(corelib.extend_list_series, profiles['exons']), [True] * len(profiles['exons']))
                avg_me = map(operator.itemgetter(0), outs)
                exoncount = map(operator.itemgetter(1), outs)
                outs = map(corelib.mean_col_by_col, map(corelib.extend_list_series, profiles['introns']), [True] * len(profiles['introns']))
                avg_mi = map(operator.itemgetter(0), outs)
                introncount = map(operator.itemgetter(1), outs)
            
            if not FIRST:    # the first chromosome profiling
                avg_upstream,avg_upcount=corelib.weight_mean_col_by_col([avg_upstream,avg_up],[avg_upcount,upcount],counts=True)
                avg_downstream,avg_downcount=corelib.weight_mean_col_by_col([avg_downstream,avg_down],[avg_downcount,upcount],counts=True)
                avg_metagene,avg_genecount=corelib.weight_mean_col_by_col([avg_metagene,avg_mg],[avg_genecount,genecount],counts=True)
                avg_metacatexon,avg_cexoncount=corelib.weight_mean_col_by_col([avg_metacatexon,avg_mce],[avg_cexoncount,cexoncount],counts=True)
                avg_metacatintron,avg_cintroncount=corelib.weight_mean_col_by_col([avg_metacatintron,avg_mci],[avg_cintroncount,cintroncount],counts=True)            
                if options.exin_pf:
                    outs = map(corelib.weight_mean_col_by_col, map(lambda x, y: [x, y], avg_metaexon, avg_me),  map(lambda x, y: [x, y], avg_exoncount, exoncount), [True] * n_metas)
                    avg_metaexon = map(operator.itemgetter(0), outs)
                    avg_exoncount = map(operator.itemgetter(1), outs)
                    outs = map(corelib.weight_mean_col_by_col, map(lambda x, y: [x, y], avg_metaintron, avg_mi), map(lambda x, y: [x, y], avg_introncount, introncount), [True] * n_metas)
                    avg_metaintron = map(operator.itemgetter(0), outs)
                    avg_introncount = map(operator.itemgetter(1), outs)
                    del outs,avg_me,avg_mi,exoncount,introncount
                del avg_up,avg_down,avg_mg,avg_mce,avg_mci,upcount,downcount,genecount,cexoncount,cintroncount
        
                if options.gn_groups:
                    if options.which=='name2':
                        ixs, subsets = profiler.get_gene_indicies(profiles['name2'], subsets)
                    else:
                        ixs, subsets = profiler.get_gene_indicies(profiles['name'], subsets)
                    
                # take an average of each profile (upstream, downstream, gene, exon, intron, cat-exon, cat-intron
                    avg_ups, upcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['upstreams'])
                    avg_downs, downcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['downstreams'])
                    avg_mgs, gcs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['genes'])
                    avg_mces, cecs = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catexons'])
                    avg_mcis, cics = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catintrons'])
                    if options.exin_pf:
                        avg_mes, ecs = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['exons'])
                        avg_mis, ics = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['introns'])
    
                    avg_upstreams,avg_upcounts=profiler.weight_mean_profiles_chr_by_chr(avg_upstreams,avg_upcounts,avg_ups,upcs)
                    avg_downstreams,avg_downcounts=profiler.weight_mean_profiles_chr_by_chr(avg_downstreams,avg_downcounts,avg_downs,downcs)
                    avg_metagenes,avg_genecounts=profiler.weight_mean_profiles_chr_by_chr(avg_metagenes,avg_genecounts,avg_mgs,gcs)
                    avg_metacatexons,avg_cexoncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metacatexons,avg_cexoncounts,avg_mces,cecs)
                    avg_metacatintrons,avg_cintroncounts=profiler.weight_mean_profiles_chr_by_chr(avg_metacatintrons,avg_cintroncounts,avg_mcis,cics)
                    if options.exin_pf:
                        avg_metaexons, avg_exoncounts = profiler.weight_mean_profiles_chr_by_chr_meta(avg_metaexons, avg_exoncounts, avg_mes, ecs)
                        avg_metaintrons, avg_introncounts = profiler.weight_mean_profiles_chr_by_chr_meta(avg_metaintrons, avg_introncounts, avg_mis, ics)
                        del avg_mes,avg_mis,ecs,ics
                        
                    del avg_ups,avg_downs,avg_mgs,avg_mces,avg_mcis,upcs,downcs,gcs,cecs,cics
            else:
                avg_upstream=avg_up
                avg_downstream=avg_down
                avg_metagene=avg_mg
                avg_metacatexon=avg_mce
                avg_metacatintron=avg_mci
                avg_upcount=upcount
                avg_downcount=downcount
                avg_genecount=genecount
                avg_cexoncount=cexoncount
                avg_cintroncount=cintroncount
                if options.exin_pf:
                    avg_metaexon=avg_me
                    avg_metaintron=avg_mi
                    avg_exoncount=exoncount
                    avg_introncount=introncount
                
                if options.gn_groups:
                    if options.which=='name2':
                        ixs, subsets = profiler.get_gene_indicies(profiles['name2'], subsets)
                    else:
                        ixs, subsets = profiler.get_gene_indicies(profiles['name'], subsets)
                    
                    avg_upstreams, avg_upcounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['upstreams'])
                    avg_downstreams, avg_downcounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['downstreams'])
                    avg_metagenes, avg_genecounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['genes'])
                    avg_metacatexons, avg_cexoncounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catexons'])
                    avg_metacatintrons, avg_cintroncounts = profiler.select_take_average_profiles_chr_by_chr(ixs, profiles['catintrons'])
                    if options.exin_pf:
                        avg_metaexons, avg_exoncounts = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['exons'])
                        avg_metaintrons, avg_introncounts = profiler.select_take_average_profiles_chr_by_chr_meta(ixs, profiles['introns'])
        jobcount+=1
        
        if options.bg:
            info('#%d perform ChIP region annotation...' %jobcount)
            Annot=annotator.Annotator()
            ChIPA=Annot.annotate(genome_coordinates=ChIPSamp,gene_table=GeneT,roi=roi,prom=options.promoter,biprom=options.bipromoter,down=options.downstream,gene_div=(3,5))
            CS,CP=Annot.summarize(ChIPA)
            CPval=annotator.estimate_pvals(GP,CS,CP)
            CES, CEP = Annot.obtain_distribution_of_sites(ChIPA)
            jobcount+=1
            
            info('#%d write R script of ChIP region annotation and wig profiling...' %jobcount)
            ofhd=open(options.name+'.R','w')
            pdfname=options.name+'.pdf'
            # the first part of CEAS R script. Because wig profiling is not run, just terminate
            rscript=R.pdf(pdfname,height=11.5,width=8.5)
            rscript+=inout.draw_CEAS(GP,CP,CPval,bg_res=options.bg_res,chip_res=options.chip_res,prom=options.promoter,biprom=options.bipromoter,down=options.downstream,gene_div=(3,5)) 
            ofhd.write(rscript)    # writing CEAS
            rscript = inout.draw_pie_distribution_of_elements(CEP, prom=options.promoter, down=options.downstream)
            ofhd.write(rscript)
            jobcount+=1
        
        # 
        # writing R script for wig profiling
        #
        if options.ceas:
            info('#%d append an R script of wig profiling...' %jobcount)
        else:
            info('#%d write an R script of wig profiling...' %jobcount)
            ofhd = open(options.name+'.R', 'w')
            pdfname = options.name + '.pdf'
            rscript = R.pdf(pdfname, height=11.5, width=8.5)
            ofhd.write(rscript)
        #
        # write a R script of drawing wig overview
        
        ### NOTE THAT this option should be run only if ChIP regions are given
        #
        #print wigsize
        #print ChIP
        rscript = inout.draw_ChIP_over_genome_mono_col(ChIP, wigsize, n_peaks=options.n_peaks)
        #rscript = inout.draw_wig_heatmap_overview(swig, ylim=None, chrominfo=None, colormap='GBR256')
        #print rscript
        ofhd.write(rscript)
        
        #
        # write a R script of drawing profiles near genes
        #
    
        # get breaks for the plots
        breaks = profiles['breaks']
        metagene_breaks = profiles['genebreaks']
        metacatexon_breaks = profiles['catexonbreaks']
        metacatintron_breaks = profiles['catintronbreaks']
        if options.exin_pf:
            metaexon_breaks = profiles['exonbreaks']
            metaintron_breaks = profiles['intronbreaks']
         
        # write R script
        if options.gn_groups:       # when multiple gene groups are given
            # append the profiles of all genes
            avg_upstreams.append(avg_upstream)
            avg_downstreams.append(avg_downstream)
            avg_metagenes.append(avg_metagene)
            avg_metacatexons.append(avg_metacatexon)
            avg_metacatintrons.append(avg_metacatintron)
            rscript=inout.draw_profile_plots(breaks,avg_upstreams,avg_downstreams,metagene_breaks,avg_metagenes,metacatexon_breaks,avg_metacatexons,metacatintron_breaks,avg_metacatintrons,metagene_breaks_lim=[-1000,1000],legends=options.gn_names)
            if options.exin_pf:
                map(lambda x, y: x.append(y), avg_metaexons, avg_metaexon)
                map(lambda x, y: x.append(y), avg_metaintrons, avg_metaintron)
                rscript+=inout.draw_exon_intron_profile_plots(metaexon_breaks, avg_metaexons,metaintron_breaks,avg_metaintrons, options.elowers, options.euppers, options.ilowers, options.iuppers, legends=options.gn_names)
        else:                       # only when a single master profiling is obtianed
            rscript=inout.draw_profile_plot(breaks,avg_upstream,avg_downstream,metagene_breaks,avg_metagene,metacatexon_breaks,avg_metacatexon,metacatintron_breaks,avg_metacatintron,metagene_breaks_lim=[-1000,1000])
            if options.exin_pf:
                rscript+=inout.draw_exon_intron_profile_plot(metaexon_breaks,avg_metaexon,metaintron_breaks,avg_metaintron, options.elowers, options.euppers, options.ilowers, options.iuppers)
    
    # write to the file and close it
    ofhd.write(rscript)    # write wig profiling
    ofhd.write(R.devoff())
    ofhd.close()
    jobcount+=1
    
    info ('#... cong! Run R on %s!' %(options.name+'.R'))


# ------------------------------------
# functions
# ------------------------------------
  
def prepare_optparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    
    """
    
    usage = "usage: %prog [options]"
    description = "CEAS -- Cis-regulatory Element Annotation System"
    
    optparser = OptionParser(version="%prog 1.0",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-b","--bed",dest="bed",type="string",
                         help="BED file of ChIP regions.")
    optparser.add_option("-w","--wig",dest="wig",type="string",
                         help="WIGGLE file for either wig profiling or genome background annotation. WARNING: --bg flag must be set for genome background annotation.")
    optparser.add_option("-e","--ebed",dest="ebed",type="string",
                         help="BED file of extra regions of interest (eg, non-coding regions)")
    optparser.add_option("-g","--gene-db",dest="gdb",type="string",
                         help="Gene annotation table (a local sqlite3 db file provided by CEAS or species name in UCSC). CEAS searches the designated directory for the the db file. If not find, it looks up UCSC for the table. WARNING: When using UCSC, MySQLdb package must be installed and a WIGGLE file must be given for genome bg annotation.")
    optparser.add_option("--span", dest="span", type="int",\
                         help="Span for gene-centered annotation. ChIP regions within this range from TSS and TTS are considered in the gene-centered annotation, DEFAULT=3000bp", default=3000)
    optparser.add_option("--bg",action="store_true",dest="bg",\
                         help="Run genome BG annotation. WARNING: This flag is effective only if a WIGGLE file is given through -w (--wig). Otherwise, ignored.",default=False)
    optparser.add_option("--name",dest="name",\
                         help="Experiment name. This will be used to name the output file. If an experiment name is not given, input BED file name will be used instead.")      
    optparser.add_option("--promoter",dest="promoter",type="int",
                         help="Promoter size for annotation, DEFAULT: 3000bp. WARNING: Value > 10000bp are automatically set to 10000bp.", default=3000)    
    optparser.add_option("--bipromoter",dest="bipromoter",type="int",
                         help="Bidirectional-promoter size for annotation, DEFAULT: 5000bp. WARNING: Value > 20000bp are automatically set to 20000bp.", default=5000)  
    optparser.add_option("--downstream",dest="downstream",type="int",
                         help="Downstream size for annotation, DEFAULT: 3000bp. WARNING: Value > 10000bp are automatically set to 10000bp.", default=3000)     
    optparser.add_option("--pf-res", dest="pf_res", type="int",\
                          help="Wig profiling resolution, DEFAULT: 50bp. WARNING: Value smaller than the wig step (resolution) may cause aliasing error.", default=50) 
    optparser.add_option("--rel-dist",dest="rel_dist",type="int",
                         help="Relative distance to TSS/TTS in wig profiling, DEFAULT: 3000bp", default=3000)   
#    optparser.add_option("--exon-lowers",dest="elowers",type="str",\
#                         help="Lower limits of the length (bp) of exons to profile in meta-exon plots. A series of numbers speparated by comma(s) can be given (eg, --exon-lowers=100,300,500). If not given, values at 10%, 35%, and 65% of exon lengths will be used.", default=None)
#    optparser.add_option("--exon-uppers",dest="euppers",type="str",\
#                         help="upper limits of the length (bp) of exons to profile in meta-exon plots. A series of numbers speparated by comma(s) can be given (eg, --exon-uppers=300,500,1000). If not given, values at 35%, 65%, and 90% of exon lengths will be used.",default=None)
#    optparser.add_option("--intron-lowers",dest="ilowers",type="str",\
#                         help="Lower limits of the length (bp) of introns to profile in meta-intron plots. A series of numbers speparated by comma(s) can be given (eg, --intron-lowers=100,300,500). If not given, values at 10%, 35%, and 65% of intron length will be used.", default=None)
#    optparser.add_option("--intron-uppers",dest="iuppers",type="str",\
#                         help="upper limits of the length (bp) of introns to profile in meta-intron plots. A series of numbers speparated by comma(s) can be given (eg, --intron-uppers=300,500,1000). If not given, values at 35%, 65%, and 90% of intron lengths will be used.",default=None)
    optparser.add_option("--exin-pf", action="store_true", dest="exin_pf",\
                         help="Switch to turn on average profiling for exons and introns", default=False)
    optparser.add_option("--gn-groups",dest="gn_groups",type="string",\
                         help="Gene-groups of particular interest in wig profiling. Each gene group file must have gene names in the 1st column. The file names are separated by commas w/ no space (eg, --gn-groups=top10.txt,bottom10.txt)") 
    optparser.add_option("--gn-group-names", dest="gn_names",type="string",\
                         help="The names of the gene groups in --gn-groups. The gene group names are separated by commas. (eg, --gn-group-names='top 10%,bottom 10%'). These group names appear in the legends of the wig profiling plots. If no group names given, the groups are represented as 'Group 1, Group2,...Group n'.")
    optparser.add_option("--no-refseq", action="store_true", dest="name2",\
                         help="Whether gene names (eg, 'name2' in refGene of UCSC) or RefSeq accession IDs (eg, 'name' in refGene of UCSC) are used in --gn-groups or not. This flag is meaningful only if --gn-groups is set.",default=False)
    
    return optparser


def opt_validate (optparser):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    
    # if gdb not given, print help, either BED or WIG must be given 
    if not options.gdb or not (options.bed or options.wig):
        optparser.print_help()
        sys.exit(1)

    HAVELOCALGDB=os.path.isfile(options.gdb)
    if not HAVELOCALGDB:
        if inout.MYSQL:
            options.gdbtype = 'ucsc'
            options.Host = "genome-mysql.cse.ucsc.edu"
            options.User = "genome"
        else:
            error('MySQLDb must be installed to use UCSC gene annotation table')
            sys.exit(0)
    else:
        if options.gdb.endswith('.bed'): # when bed input
            if options.name2:
                error('--no-refseq does not work with a BED formatted gene annotation table (-g or --gdb)')
                sys.exit(0)
            options.gdbtype = 'localbed'
        else:
            options.gdbtype = 'localdb'
            options.Host = None
            options.User = None
            
#    # get gdb lower case
#    HAVELOCALGDB=os.path.isfile(options.gdb)
#    if HAVELOCALGDB:
#        options.Host=None
#        options.User=None
#        options.Db=options.gdb
#    else:
#        if MYSQL:
#            options.Host="genome-mysql.cse.ucsc.edu"
#            options.User="genome"
#            options.Db=os.path.split(options.gdb)[-1].lower()
#        else:
#            error('MySQLdb package needs to be installed to use UCSC or a local sqlite3 db file must exist.')
#            error('Check -g (--gdb). No such file: %s' %options.gdb)
#            sys.exit(1)
    
    # bed file
    if options.bed:
        HAVEBED = os.path.isfile(options.bed)
        if not HAVEBED:
            error('Check -b (--bed). No such bed file: %s' %options.bed)
            sys.exit(1)
    else: HAVEBED = False
    
    # bg background resolution is set to 100    
    options.bg_res=100
    
    # wig file
    if options.wig:
        HAVEWIG=os.path.isfile(options.wig)
    else: HAVEWIG=False
    
    REBG=False
    if HAVEWIG and HAVEBED and options.bg:
        REBG=True
    elif not HAVELOCALGDB and HAVEWIG and HAVEBED and not options.bg:
        error('Genome BG annotation must be run when no pre-computed BG annotation exists. Set --bg.')
        sys.exit(1)
    elif not HAVELOCALGDB and not HAVEWIG and HAVEBED and options.bg:
        if not options.wig:
            error('Genome BG annotation must be run when no pre-computed BG annotation exists. Give a WIGGLE file through -w (--wig)')
        else:
            error('Genome BG annotation must be run when no pre-computed BG annotation exists. Check -w (--wig). No such file: %s' %options.wig)
        sys.exit(1)
    elif not HAVELOCALGDB and not HAVEWIG and HAVEBED and not options.bg:
        if not options.wig:
            error('Genome BG annotation must be run when no pre-computed BG annotation exists. set --bg and give a WIGGLE file through -w (--wig)')
        else:
            error('Genome BG annotation must be run when no pre-computed BG annotation exists. Set --bg and check -w (--wig). No such file: %s' %options.wig)
        sys.exit(1)
    options.bg=REBG
    
    # only when a bed file is given, annotation will run
    options.ceas=HAVEBED   
    
    # if a WIGGLE file given, do profiling
    options.pf = HAVEWIG
        
    # non-coding regions
    if options.ebed:
        if not os.path.isfile(options.ebed):
            error('Check -e (--ebed). No such file: %s' %options.ebed)
            sys.exit(1)
    
     # get namename
    if not options.name:
        if HAVEBED:
            options.name=os.path.split(options.bed)[-1].rsplit('.bed',2)[0]
        elif HAVEWIG:
            options.name=os.path.split(options.wig)[-1].rsplit('.wig',2)[0]
    
    # the minimum ChIP annotation resolution is 1000
    options.chip_res = 1000
    
    # the maximum promoter, bidirectional promoter and downstream
    options.promoter=min(10000, options.promoter)
    options.bipromoter=min(20000, options.bipromoter)
    options.downstream=min(10000, options.downstream)
    
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
    
    # profiling resolution
    options.pf_res = max(1, options.pf_res)
    # relative distance
    options.rel_dist = max(options.rel_dist, options.pf_res)
    # metagene_size
    options.metagene_size = 3000
    options.catexon_size = options.metagene_size/2
    options.catintron_size = options.metagene_size/2
    
    # exon length and intron lengths 
    # consider exons with lengths of [lower lim (%), upper lim (%)]
    # consider introns with lengths of [lower lim (%), upper lim (%)]
    options.epercentlowers=[10, 35, 65]    # in percent
    options.epercentuppers=[35, 65, 90]    # in percent
    options.ipercentlowers=[10, 35, 65]    # in percent
    options.ipercentuppers=[35, 65, 90]    # in percent
    
    options.elowers = None
    options.euppers = None
    options.ilowers = None
    options.iuppers = None
    
    options.metaexonsizes = None
    options.metaintronsizes = None
    #parse and convert to an integer list
#    if options.elowers:
#        options.elowers = map(int, options.elowers.split(','))
#    
#    if options.euppers:
#        options.euppers = map(int, options.euppers.split(','))
#        
#    if options.ilowers:
#        options.ilowers = map(int, options.ilowers.split(','))
#    
#    if options.iuppers:
#        options.iuppers = map(int, options.iuppers.split(','))
#        
#    # check if elowers and euppers have the same number of elements
#    if options.elowers and options.euppers:
#        if len(options.elowers)!=len(options.euppers):
#            error('ELOWERS and EUPPERS must have the same number of elements')
#            sys.exit(0)
#        else:
#            for elower, eupper in itertools.izip(options.elowers, options.euppers):
#                if eupper <= elower:
#                    error('EUPPER value must be larger than its corresponding ELOWER value')
#                    sys.exit(0)
#    elif (not options.elowers and options.euppers) or (options.elowers and not options.euppers): 
#        error('ELOWERS and EUPPERS must be given together')
#        sys.exit(0)
#    
#    # check if ilowers and iuppers have the same number of elements
#    if options.ilowers and options.iuppers:
#        if len(option.ilowers)!=len(options.iuppers):
#            error('ILOWERS and IUPPERS must have the same number of elements')
#            sys.exit(0)
#        else:
#            for ilower, iupper in itertools.izip(options.ilowers, options.iuppers):
#                if iupper <= ilower:
#                    error('IUPPER value must be larger than its corresponding ILOWER value')
#                    sys.exit(0)
#    elif (not options.ilowers and options.iuppers) or (options.ilowers and not options.iuppers):
#        error('ILOWERS and IUPPERS must be given together')
#        sys.exit(0)
    
    # number of wig overview resolution
    options.n_peaks = 10000
    
        
    return options


def get_min_max_lengths(GeneT, epercentlimit, ipercentlimit, minexonlen, minintronlen):
    """Return the minimum and maximum lengths of exons or introns in the gene annotation table
    
    Parameters:
    1. GeneT: gene annotation table
    2. epercentlimit: [lower limit for exon length, upper limit for exon length] to consider. lower limit and upper limit must be percentage values (0-100)
    3. ipercentlimit: [lower limit for intron length, upper limit for exon lentht] to consider. 
    4. minexonlen: mininum exon length to consider
    5. minintronlen: minimum intron length to consider
    """
    
    exonLens, intronLens=GeneT.get_exon_intron_lens()
    exon_enough_long, intron_enough_long = corelib.get_certain_part(exonLens, percentlimit=epercentlimit), corelib.get_certain_part(intronLens, percentlimit=ipercentlimit)
    exonLenLim, intronLenLim = [max(minexonlen, exon_enough_long[0]), exon_enough_long[-1]], [max(minintronlen, intron_enough_long[0]), intron_enough_long[-1]]

    return exonLenLim, intronLenLim


def return_med_gene(GeneT, pf_res):
    """Return the median gene length of the given gene annotation table.
    
    The nearest multiple of pf_res to the median gene length will be returned.
    
    Parameters:
    1. GeneT: the gene annotation table
    2. pf_res: the profiling resolution
    
    """
    
    gLens = GeneT.get_gene_lens()
    medgLen = corelib.median(gLens)
    
    return corelib.find_nearest_multiple(medgLen, pf_res)


def return_med_catexon_catintron(GeneT, pf_res):
    """Return the concatenated exons and introns sizes for scaling
    
    Parameters:
    1. GeneT: the gene annotation table to consider
    2. pf_res: profiling resolution
    
    """
        
    catexonLens, catintronLens =GeneT.get_cat_exon_intron_lens()
    medcatexonLen, medcatintronLen = corelib.median(catexonLens), corelib.median(catintronLens)
    medcatexonLen, medcatintronLen = corelib.find_nearest_multiple(medcatexonLen, pf_res), corelib.find_nearest_multiple(medcatintronLen, pf_res)

    return medcatexonLen, medcatintronLen

def return_med_exons_introns(GeneT, eplowers, epuppers, iplowers, ipuppers):
    """Return the median exon lengths and intron lengths
    
    Parameters:
    1. GeneT: genome annotation table
    2. eplowers: a list of numbers (1-100) that indicate lower percentage limits for exon length
    3. epuppers: a list of numbers (1-100) that indicate upper percentage limits for exon length
    4. iplowers: a list of numbers (1-100) that indicate lower percentage limits for intron length
    5. ipuppers: a list of numbers (1-100) that indicate upper percentage limits for intron length
    
    """
    
    exonLens, intronLens = GeneT.get_exon_intron_lens()
    elowers, medexonsizes, euppers = corelib.get_boundaries_medians(exonLens, lowers = eplowers, uppers = epuppers)
    ilowers, medintronsizes, iuppers = corelib.get_boundaries_medians(intronLens, lowers = iplowers, uppers = ipuppers)
    
    return elowers, medexonsizes, euppers, ilowers, medintronsizes, iuppers


def determine_exon_intron_sizes(GeneT, options):
    """Determine the exon and intron sizes to consider in exon-intron average profiling"""
    
    # get the exon lengths and intron lengths when either elowers or ilowers is not given.
    if not (options.elowers and options.ilowers):
        elowers, medexonsizes, euppers, ilowers, medintronsizes, iuppers = return_med_exons_introns(GeneT, options.epercentlowers, options.epercentuppers, options.ipercentlowers, options.ipercentuppers)
    
    # update the options with the exonsizes and intronsizes
    if options.elowers and not options.ilowers:
        options.ilowers = ilowers
        options.iuppers = iuppers
        medexonsizes = map(lambda x,y: (x+y)/2, options.elowers, options.euppers)
    elif not options.elowers and options.ilowers:
        options.elowers = elowers
        options.euppers = euppers
        medintronsizes  = map(lambda x,y: (x+y)/2, options.ilowers, options.iuppers)
    elif not options.elowers and not options.ilowers:
        options.elowers = elowers
        options.euppers = euppers
        options.ilowers = ilowers
        options.iuppers = iuppers
    else:
        medexonsizes = map(lambda x,y: (x+y)/2, options.elowers, options.euppers)
        medintronsizes = map(lambda x,y: (x+y)/2, options.ilowers, options.iuppers)
    
    # approximate the exonsizes to the nearest the multiple of the profiling resolution
    n_ranges_e = len(options.elowers)
    n_ranges_i = len(options.ilowers)
#    options.elowers = map(corelib.find_nearest_multiple, options.elowers, [options.pf_res]*n_ranges_e)
#    options.euppers = map(corelib.find_nearest_multiple, options.euppers, [options.pf_res]*n_ranges_e)
#    options.ilowers = map(corelib.find_nearest_multiple, options.ilowers, [options.pf_res]*n_ranges_i)
#    options.iuppers = map(corelib.find_nearest_multiple, options.iuppers, [options.pf_res]*n_ranges_i)
    options.metaexonsizes = map(corelib.find_nearest_multiple, map(max, medexonsizes, [2*options.pf_res]*n_ranges_e), [options.pf_res]*n_ranges_e)
    options.metaintronsizes = map(corelib.find_nearest_multiple, map(max, medintronsizes, [2*options.pf_res]*n_ranges_i), [options.pf_res]*n_ranges_i)

#    options.metaexonsizes = medexonsizes
#    options.metaintronsizes = medintronsizes
    
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
    new={}
    new['promoter'] = [promoter/3, 2*promoter/3, promoter]
    new['bipromoter'] = [bipromoter/2, bipromoter]
    new['downstream'] = [downstream/3, 2*downstream/3, downstream]
    
    # the model promoter, bipromoter, and downstream lengths
    mod={}
    mod['promoter']=[0, 500]+corelib.seq(fr=1000,to=10000,by=1000)
    mod['bipromoter']=[0, 500]+corelib.seq(fr=1000,to=20000,by=1000)
    mod['downstream']=[0, 500]+corelib.seq(fr=1000,to=10000,by=1000)
    
    for chrom in GenomeBGP.get_chroms():
        GP.init_table(chrom)
        for column in GenomeBGP.columns[1:]:
            if column not in ['promoter','bipromoter','downstream']:
                GP[chrom][column]=GenomeBGP[chrom][column]
            else:
                vals=[0.0]+GenomeBGP[chrom][column]
                interpol=[]
                for x in new[column]:
                    i=corelib.findbin(x, mod[column])
                    interpol.append(corelib.lininterpol([mod[column][i],vals[i]], [mod[column][i+1],vals[i+1]],x))
                
                GP[chrom][column]=interpol
                
    return GP
                    
                    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0)
