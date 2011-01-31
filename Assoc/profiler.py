

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
# Python modules
# ------------------------------------
import sys,operator,warnings
from math import *
from array import *
import copy
import logging

# ------------------------------------
# My own Python modules
# ------------------------------------
from Cistrome.Assoc.inout import *
from Cistrome.Assoc.corelib import *

#-------------------------------------
# classes
#-------------------------------------  
    
class WigProfiler2:
    """Class WigProfiler"""
    
    def __init__(self, rel_dist, step, metagenesize, catexonsize, catintronsize, metaexonsizes, metaintronsizes, elowers, euppers, ilowers, iuppers):
        """Null constructor"""
        
        self.wigint = 10 
        self.step = step
        self.rel_dist = rel_dist
        self.metagenesize = metagenesize
        self.catexonsize = catexonsize
        self.catintronsize = catintronsize
        self.metaexonsize = metaexonsizes
        self.metaintronsize = metaintronsizes
        self.elowers = elowers
        self.euppers = euppers
        self.ilowers = ilowers
        self.iuppers = iuppers
        self.wingsize = 1000
        self.binfunc = 'middle'
    
    
    def set_parameters(self, rel_dist, step, metagenesize, catexonsize, catintronsize, metaexonsizes, metaintronsizes, elowers, euppers, ilowers, iuppers):
        """Set parameters for profiling
        """
        
        self.rel_dist = rel_dist
        self.step = step
        self.metagenesize = metagenesize
        self.catexonsize = catexonsize
        self.catintronsize = catintronsize
        self.metaexonsizes = metaexonsizes
        self.metaintronsizes = metaintronsizes
        self.elowers = elowers
        self.euppers = euppers
        self.ilowers = ilowers
        self.iuppers = iuppers
        self.wingsize = 1000
        self.binfunc = 'middle'
        
        
    def profile(self, wig, gene_table):
        """Profile signal on and around genes"""
        
        # get the current profiling parameters
        step = self.step
        rel_dist = self.rel_dist
        metagenesize = self.metagenesize
        catexonsize = self.catexonsize
        catintronsize = self.catintronsize
        metaexonsizes = self.metaexonsize
        metaintronsizes = self.metaintronsize
        elowers = self.elowers
        euppers = self.euppers
        ilowers = self.ilowers
        iuppers = self.iuppers
        wingsize = min(self.wingsize, self.rel_dist)
        
        # get the chromosomes shared by wig and gene_table
        chroms = list(set(gene_table.get_chroms()).intersection(set(wig.get_chroms())))
        
        hf = 1.0*step/2 # half the step size
        
        # initialize profiles
        profiles = {}
        # iterate through chromosomes
        for chrom in chroms:
            curr_chrom = gene_table[chrom]
            curr_wig = wig[chrom]
            n_genes = len(gene_table[chrom]['txStart'])
            
            # copy the gene names
            try:
                profiles['name'].extend(curr_chrom['name'])
            except KeyError:
                profiles['name'] = curr_chrom['name']
            
            if curr_chrom.has_key('name2'):
                try: 
                    profiles['name2'].extend(curr_chrom['name2'])
                except KeyError:
                    profiles['name2'] = curr_chrom['name2']
            
            # obtain the start and end positions of genes and their surrounding areas
            rstarts = map(lambda txStart: txStart - rel_dist - hf, curr_chrom['txStart'])
            rends = map(lambda txEnd: txEnd + rel_dist + hf, curr_chrom['txEnd'])
            
            generegions = self.extract_regions(curr_wig, starts=rstarts, ends=rends)
            
            ustarts = map(self.return_upstart, curr_chrom['txStart'], curr_chrom['txEnd'], curr_chrom['strand'], [rel_dist] * n_genes)
            uends = map(self.return_upend, curr_chrom['txStart'], curr_chrom['txEnd'], curr_chrom['strand'], [rel_dist] * n_genes)
            
            gstarts = map(lambda txStart: txStart - step - hf, curr_chrom['txStart'])
            gends =  map(lambda txEnd: txEnd + step + hf, curr_chrom['txEnd'])
            
            dstarts = map(self.return_downstart, curr_chrom['txStart'], curr_chrom['txEnd'], curr_chrom['strand'], [rel_dist] * n_genes)
            dends = map(self.return_downend, curr_chrom['txStart'], curr_chrom['txEnd'], curr_chrom['strand'], [rel_dist] * n_genes)
            
            upstreams, genes, downstreams = self.extract_up_gene_down_from_regions(curr_wig, generegions, ustarts, gstarts, dstarts, uends, gends, dends)
            
            # get the exons 
            estarts = map(lambda exonStarts: map(lambda exonStart: exonStart - hf, exonStarts), curr_chrom['exonStarts'])
            eends = map(lambda exonEnds: map(lambda exonEnd: exonEnd + hf, exonEnds), curr_chrom['exonEnds'])
            exons = self.extract_regions_from_given_genes(curr_wig, genes, starts=estarts, ends=eends)
            
            # get the introns
            istarts = map(lambda exonEnds: map(lambda exonEnd: exonEnd - hf, exonEnds[:-1]), curr_chrom['exonEnds'])
            iends = map(lambda exonStarts: map(lambda exonStart: exonStart + hf, exonStarts[1:]), curr_chrom['exonStarts'])
            introns = self.extract_regions_from_given_genes(curr_wig, genes, starts=istarts, ends=iends)
    
        
            ### do binning for those regions
            ### From here on, downstream, upstreams, genes are binned signals
            # do binning for the upstreams: warning! here upstreams is a nested list of binned signals
            upstreams = self.do_binning(curr_wig, upstreams, starts=ustarts, ends=uends, strands=curr_chrom['strand'])
    
            # do binning for the downstreams
            downstreams = self.do_binning(curr_wig, downstreams, starts=dstarts, ends=dends, strands=curr_chrom['strand'])
    
            # do binning for the genes, including scaling too
            ### NOTE: genes, exons, introns are now nested lists of binned signals
            genes = self.do_binning_gene(curr_wig, genes, starts=gstarts, ends=gends, strands=curr_chrom['strand'], metasize=metagenesize)
            catexons = self.do_binning_cat(curr_wig, exons, starts=estarts, ends=eends, strands=curr_chrom['strand'], metasize=catexonsize)
            catintrons = self.do_binning_cat(curr_wig, introns, starts=istarts, ends=iends, strands=curr_chrom['strand'], metasize=catintronsize)
            
            # save
            try:
                profiles['upstreams'].extend(upstreams)
                profiles['downstreams'].extend(downstreams)
                profiles['genes'].extend(genes)
                profiles['catexons'].extend(catexons)
                profiles['catintrons'].extend(catintrons)
            except KeyError:        # when the first chromosome
                profiles['upstreams'] = upstreams
                profiles['downstreams'] = downstreams
                profiles['genes'] = genes
                profiles['catexons'] = catexons
                profiles['catintrons'] = catintrons
            
            # delete variables
            del upstreams, downstreams, genes, catexons, catintrons
            
            # exon and intron profiling when the switch is on
            if elowers and euppers and ilowers and iuppers:
                exons = self.do_binning_meta(curr_wig, exons, starts=estarts, ends=eends, strands=curr_chrom['strand'], metasizes=metaexonsizes, lowers=elowers, uppers=euppers)
                introns = self.do_binning_meta(curr_wig, introns, starts=istarts, ends=iends, strands=curr_chrom['strand'], metasizes=metaintronsizes, lowers=ilowers, uppers=iuppers)
                
                # exons = [A, B, C], where A, B, and C represent exon groups by length. Each of them has a structure of [[E1-1, E1-2, E1-3], [E2-1, E2-2], ...], 
                try:
                    map(lambda x, y: x.extend(y), profiles['exons'], exons)
                    map(lambda x, y: x.extend(y), profiles['introns'], introns)
                except KeyError:    # when the first chromosome
                    profiles['exons'] = exons
                    profiles['introns'] = introns
                
                # delete variables
                del exons, introns
            
            # delete variables that were used
            del gstarts, gends, estarts, eends, istarts, iends, ustarts, uends, dstarts, dends
            
        # get breaks
        profiles['breaks'] = self.get_breaks(-1*rel_dist, rel_dist)
        profiles['genebreaks'] = self.get_breaks(0, metagenesize)
        profiles['catexonbreaks'] = self.get_breaks(0, catexonsize, bp=False)
        profiles['catintronbreaks'] = self.get_breaks(0, catintronsize, bp=False)
        
        if elowers and euppers and ilowers and iuppers:    
            profiles['exonbreaks'] = map(lambda size:  self.get_breaks(0, size, bp=False), metaexonsizes)
            profiles['intronbreaks'] = map(lambda size: self.get_breaks(0, size, bp=False), metaintronsizes)
            
        return  profiles
        
            
    def return_upstart(self, start, end, strand, span):
        """Return the start location of the upstream
        """
        
        hf = 1.0*self.step/2
        if strand=='+':
            return start - span - hf
        else:
            return end - span - hf
        
        
    def return_upend(self, start, end, strand, span):
        """Return the end location of the upstream
        """
        
        hf = 1.0*self.step/2
        if strand=='+':
            return start + span + hf
        else:
            return end + span + hf
        
        
    def return_downstart(self, start, end, strand, span):
        """Return the start location of the downstream
        """
        
        hf = 1.0*self.step/2
        if strand=='+':
            return end - span - hf
        else:
            return start - span - hf
        
        
    def return_downend(self, start, end, strand, span):
        """Return the end location of the downstream
        """
        
        hf = 1.0*self.step/2
        if strand=='+':
            return end + span + hf
        else:
            return start + span + hf
        
                
    def estimate_wig_interval(self, wig):
        """Estimate median interval between points next to each other. This method is exactly the same as estimate_wig_interval function in inout.py.
        """   
        
        chroms = wig.get_chroms()
        if not chroms: return None
    
        n_random_positions = 10
        intervals = []
        for chrom in chroms:
            len_this_chr = len(wig[chrom][0])
            a = randints(0, len_this_chr - 2, 2 * n_random_positions)   # we need at least two element array to get difference
            a.sort()
            starts = [a[i] for i in xrange(len(a)) if i%2 == 0]
            ends = [a[i] + 2 for i in xrange(len(a)) if i%2 == 1]# we need at least two element array to get difference

            for start, end in itertools.izip(starts, ends):
                intervals.append(median(diff(wig[chrom][0][start:end])))
             
        self.wigint = median(intervals)


    def extract_regions(self, curr_wig, starts, ends):
        """Extract regions of interest from the Wig object
        """
        
        # regions = [[start_index0, end_index0], [start_index1, end_index1], ...]
        regions = []
        init = 0
        for start, end in itertools.izip(starts, ends):
            s, e = where(start, end, curr_wig[0][init:]) # use where2 in corelib.py
            regions.append([init+s, init+e])    
            init += s
        
        return regions
    
    
    def extract_up_gene_down_from_regions(self, curr_wig, regions, ustarts, gstarts, dstarts, uends, gends, dends):
        """Extract upstreams, genes, and downstreams from the regions"""
        
        upstreams = []
        genes = []
        downstreams = []
        
        for region, ustart, gstart, dstart, uend, gend, dend in itertools.izip(regions, ustarts, gstarts, dstarts, uends, gends, dends):
#            init = 0
#            s, e = where(ustart, uend, curr_wig[0][region[0]+init:region[1]])
#            upstreams.append([region[0]+init+s, region[0]+init+e])
#            init = s
#            
#            s, e = where(gstart, gend, curr_wig[0][region[0]+init:region[1]])
#            genes.append([region[0]+init+s, region[0]+init+e])
#            
#            s, e = where(dstart, dend, curr_wig[0][region[0]+init:region[1]])
#            downstreams.append([region[0]+init+s, region[0]+init+e])

            s, e = where(ustart, uend, curr_wig[0][region[0]:region[1]])
            upstreams.append([region[0]+s, region[0]+e])
            
            s, e = where(gstart, gend, curr_wig[0][region[0]:region[1]])
            genes.append([region[0]+s, region[0]+e])
            
            s, e = where(dstart, dend, curr_wig[0][region[0]:region[1]])
            downstreams.append([region[0]+s, region[0]+e])
    
        return upstreams, genes, downstreams
            
    
    def extract_regions_from_given_genes(self, curr_wig, genes, starts, ends):
        """extract regions such as exons or introns given Wig and genes"""
        
        regions = []
        for gene, start, end in itertools.izip(genes, starts, ends):
            this_regions = []            
            init = 0
            slice = curr_wig[0][gene[0]:gene[1]]
            for s, e in itertools.izip(start, end):
                s2, e2 = where(s, e, slice[init:])
                this_regions.append([gene[0]+init+s2, gene[0]+init+e2])
                init += s2
                
            regions.append(this_regions)
            
        return regions


    def get_breaks(self, start, end, bp=True):
        """Return meta gene or exon or intron breaks
        
        Parameters:
        1. start: the start value
        2. end: the end value. This end value is included in the resulting breaks.
        3. bp: True - breaks in bp, False - breaks in percentage
        
        """
        step = self.step
        n=(end-start)/step+1
        if bp:
            meta_breaks = map(lambda x: int(round(x)), linspace(start, end, n))
        else:
            meta_breaks = linspace(0, 100, n)
        
        return meta_breaks
    
    
    def do_binning(self, curr_wig, regions, starts, ends, strands):
        """Do binning for the regions"""
        
        step = self.step
        binned = []
        for region, start, end, strand in itertools.izip(regions, starts, ends, strands):
            
            # if region[0] == region[1], [] will be return
            x = curr_wig[0][region[0]:region[1]]
            y = curr_wig[1][region[0]:region[1]]
            
            if x and y: # if there are some signals in x and y, otherwise, just put null []
                n = int(end-start)/step
                bins = linspace(start, end, n + 1)
                this = binxy_equibin(bins, x, y, binfunc=self.binfunc, NaN=False)
                if strand != '+': this.reverse()
                binned.append(this)
            else:
                binned.append([])
                
        return binned


    def do_binning_gene(self, curr_wig, regions, starts, ends, strands, metasize):
        """Do binning for the regions"""
        
        step = self.step
        m = (metasize+step)/step
        binned = []
        for region, start, end, strand in itertools.izip(regions, starts, ends, strands):
            
            # if region[0] == region[1], [] will be return
            x = curr_wig[0][region[0]:region[1]]
            y = curr_wig[1][region[0]:region[1]]
            
            if x and y: # if there are some signals in x and y, otherwise, just put null []                                   
                n = int(end-start)/step
                bins = linspace(start, end, n + 1)
                this = binxy_equibin(bins, x, y, binfunc=self.binfunc, NaN=False)
                if strand != '+': this.reverse()
                if not this: binned.append([])
                if n-2 > m:
                    this = self.decimate(bins[1:-1], this[1:-1], m)
                elif n-2 < m:
                    this = self.interpol(bins[1:-1], this[1:-1], m, this[0])
                else:
                    this = this[1:-1]
                binned.append(this)
            else:
                binned.append([])
                
        return binned
    
    
    def do_binning_cat(self, curr_wig, cats, starts, ends, strands, metasize):
        """do binning for the concatenated exons or introns"""

        step = self.step
        wigint = self.wigint
        m = (metasize+step)/step
        binned = []
        for cat, start, end, strand in itertools.izip(cats, starts, ends, strands):
            if cat and start and end: 
                x = []
                y = []
                gaps = cumsum([0] + map(lambda a, b: a-b, start[1:], end[:-1]))    # the cumulative gaps between exons or between introns: using cumsum in corelib.py
                for c, gap in itertools.izip(cat, gaps):
                    slice = curr_wig[0][c[0]:c[1]]
                    if slice:
                        try:
                            d = max(wigint - (slice[0] - gap - x[-1]), 0)          # put at least wig interval between exons or (between introns)
                        except IndexError:  # if x is empty 
                            d = 0
                        x += map(lambda a: a - gap + d, slice)
                    else:
                        x += []
                    y += curr_wig[1][c[0]:c[1]]
            
                if x and y: # only if x and y exist, otherwise, put []
                    adj_end = map(lambda a, offset: a - offset, end, gaps)
                    n = int(adj_end[-1]-start[0])/step
                    bins = linspace(start[0], adj_end[-1], n+1)
                    this = binxy_equibin(bins, x, y, binfunc=self.binfunc, NaN=False)
                    if strand != '+': this.reverse()
                    if n > m:
                        this = self.decimate(bins, this, m)
                    elif n < m:
                        this = self.interpol(bins, this, m, this[0])
                    binned.append(this)
                else:
                    binned.append([])
            else:
                binned.append([])
                          
        return binned
            
            
    def do_binning_meta(self, curr_wig, metas, starts, ends, strands, metasizes, lowers, uppers):
        """Do binning for the exons and introns"""
        
        step = self.step
        ms = map(lambda m: (m+step)/step, metasizes)
        binned = map(lambda x: [], lowers)
        for meta, start, end, strand in itertools.izip(metas, starts, ends, strands):
            this_binned = map(lambda x: [], lowers)            
            for m, s, e in itertools.izip(meta, start, end):
                x = curr_wig[0][m[0]:m[1]]
                y = curr_wig[1][m[0]:m[1]]
                
                if x and y: # only when x and y exist, do binning. Otherwise, just put []
                    realsize = e - s - step
                    n = int(e - s)/step    
                    for m, lower, upper, this in itertools.izip(ms, lowers, uppers, this_binned):
                        if realsize >=lower and realsize < upper and n > 0:
                            bins = linspace(s, e, n+1)
                            th = binxy_equibin(bins, x, y, binfunc=self.binfunc, NaN=False)
                            if strand != '+': th.reverse()
                            if n > m:
                                th = self.decimate(bins, th, m)
                            elif n < m:
                                th = self.interpol(bins, th, m, th[0])
                            this.append(th)
                            break
                        
            # append binned    
            for b, this in itertools.izip(binned, this_binned):
                b.append(this)
                
        return binned

    
    def interpol(self, bins, binned, meta_n, beforestart):
        """Interpolate (in case the number of bins) to fit to meta"""
        
        hf = 1.0*self.step/2
        x = map(lambda b: b + hf, bins[:-1])
        metabins = linspace(bins[0], bins[-1], meta_n + 1)
        metabinned = binxy_equibin(metabins, x, binned, binfunc=self.binfunc, NaN=True) # binning first and fill out NaNs with their prior values
        prior = beforestart
        interpolated = []
        for mb in metabinned:
            if mb==mb: 
                prior = mb
                interpolated.append(mb)
            else: 
                interpolated.append(prior)
        
        return interpolated
    
        
    def decimate(self, bins, binned, meta_n):
        """Decimate (reduce the number of bins to fit to newn"""
        
        hf = 1.0*self.step/2
        x = map(lambda b: b + hf, bins[:-1])
        metabins = linspace(bins[0], bins[-1], meta_n + 1)
        
        return binxy_equibin(metabins, x, binned, binfunc=self.binfunc, NaN=False)
                      
                    
    def append_up_down_to_genes(self, upstreams, genes, downstreams):
        """Append upstream and downstream to the genes. Upstreams, genes, and downstreams are all binned data
        """
        
        wingsize = min(self.wingsize, self.rel_dist)
        rel_dist = self.rel_dist
        step = self.step
        # get the breaks
        breaks = self.get_breaks(-1*rel_dist, rel_dist, bp=True)
        su, eu = where(-1*wingsize, -1*step, breaks, step)
        sd, ed = where(step, wingsize, breaks, step)
        
        return map(lambda upstream, gene, downstream: upstream[su:eu] + gene + downstream[sd:ed], upstreams, genes, downstreams)


#-------------------------------------
# functions
#-------------------------------------  
def enumerate_genes(genes):
    """Enumerate genes"""
    
    schema={}
    for i,x in enumerate(genes): 
        if schema.has_key(x):
            schema[x].append(i)
        else:
            schema[x] = [i]
    
    return schema

def get_gene_indicies(genes,subsets):
    """Return the indicies of genes belonging to subsets"""
    
    # get the table
    schema=enumerate_genes(genes)
    
    ixs=[]
    missing_genes=[]
    for sub in subsets:
        ix=[]
        missing_gene=[]
        for s in sub:
            try:
                ix.extend(schema[s])
            except KeyError,e:
#                warnings.warn("%s does not exist in the gene annotation table" %e)
                missing_gene.append(s)
                pass
        ixs.append(ix)
        missing_genes.append(missing_gene)
    
    return ixs,missing_genes

def select_profiles(ix,profiles):
    """Select profiles based on the indicies returned by get_gene_indicies"""
    
    subprofiles=[]
    for ix in ixs:
        subprofiles.append(array_extractor(profiles,ix))
    
    return subprofiles


def select_take_average_profiles_chr_by_chr(ixs, profiles):
    """Select profiles using the sub indexes"""
    
    subavg = []
    subc = []
    for ix in ixs:
        try:
            a, c = mean_col_by_col(array_extractor(profiles, ix), counts = True)
            subavg.append(a)
            subc.append(c)
        except: # in case ix is empty
            subavg.append([])
            subc.append([])
        
    return subavg, subc


def select_take_average_profiles_chr_by_chr_meta(ixs, metaprofiles):
    """Do like select_take_average_profiles_chr_by_chr for metaexons and metaintrons"""
    
    subavgs = map(lambda x: [], metaprofiles)
    subcs =copy.deepcopy(subavgs)
    
    for subavg, subc, metaprofile in itertools.izip(subavgs, subcs, metaprofiles):
        for ix in ixs:
            try:
                a, c = mean_col_by_col(extend_list_series(array_extractor(metaprofile, ix)), counts = True)
                subavg.append(a)
                subc.append(c)
            except: # in case ix is empty
                subavg.append([])
                subc.append([])
            
    return subavgs, subcs

    
def select_profiles_chr_by_chr(ixs, upstreams, downstreams, metagenes, metacatexons, metacatintrons, metaexons, metaintrons):
    """Select profiles according to subset ixs"""
    
    # initialize the output values
    avg_upstreams=[]
    avg_downstreams=[]
    avg_metagenes=[]
    avg_metacatexons=[]
    avg_metacatintrons=[]
    avg_metaexons=[]
    avg_metaintrons=[]
    avg_upcounts=[]
    avg_downcounts=[]
    avg_genecounts=[]
    avg_cexoncounts=[]
    avg_cintroncounts=[]
    avg_exoncounts=[]
    avg_introncounts=[]
    
    for ix in ixs:
        # take the averages
        #try:
        # upstream
        a, c=mean_col_by_col(array_extractor(upstreams,ix),counts=True)
        avg_upstreams.append(a)
        avg_upcounts.append(c)
        
        #downstream
        a, c=mean_col_by_col(array_extractor(downstreams,ix),counts=True)
        avg_downstreams.append(a)
        avg_downcounts.append(c)
        
        #metagene
        a, c=mean_col_by_col(array_extractor(metagenes,ix),counts=True)
        avg_metagenes.append(a)
        avg_genecounts.append(c)
        
        #meta concatenated exons
        a, c=mean_col_by_col(array_extractor(metacatexons,ix),counts=True)
        avg_metacatexons.append(a)
        avg_cexoncounts.append(c)
        
        #meta concatenated introns
        a, c=mean_col_by_col(array_extractor(metacatintrons,ix),counts=True)
        avg_metacatintrons.append(a)
        avg_cintroncounts.append(c)
        
        #meta exons
        a, c=mean_col_by_col(extend_list_series(array_extractor(metaexons,ix)),counts=True)
        avg_metaexons.append(a)
        avg_exoncounts.append(c)
        
        #meta introns
        a, ct=mean_col_by_col(extend_list_series(array_extractor(metaintrons,ix)),counts=True)
        avg_metaintrons.append(a)
        avg_introncounts.append(c)
                
    return avg_upstreams, avg_upcounts, avg_downstreams, avg_downcounts, avg_metagenes, avg_genecounts, avg_metacatexons, avg_cexoncounts, avg_metacatintrons, avg_cintroncounts, avg_metaexons, avg_exoncounts, avg_metaintrons, avg_introncounts

def weight_mean_profiles_chr_by_chr(avgs, avgcounts, avg, avgcount):
    """Update sub set profiles in a chromosome by chromsome manner"""
    
    for i in xrange(0,len(avgs)):
        try:
            avgs[i],avgcounts[i]=weight_mean_col_by_col([avgs[i],avg[i]],[avgcounts[i],avgcount[i]],counts=True)
        except IndexError: # sometimes, avg is empty. 
            pass
        
    return avgs,avgcounts


def weight_mean_profiles_chr_by_chr_meta(avgs, avgcounts, avg, avgcount):
    """Do like wieght_mean_profiles_chr_by_chr_meta
    
    avgs = [A, B, C], where A, B, and C represent three groups of exons (or introns) according to their lengths.
    A = [G1, G2], where G1, and G2 present gene groups (for example, top 10% and bottom 10% of expressed genes)
    """
    
    for av, avcount, a, ac in itertools.izip(avgs, avgcounts, avg, avgcount):    
        for i in xrange(0, len(av)):
            try:
                av[i], avcount[i] = weight_mean_col_by_col([av[i], a[i]], [avcount[i], ac[i]], counts=True)
            except IndexError:
                pass
        
    return avgs, avgcounts


def return_genetables_from_gene_lists(gene_lists,Host='',User='',Db='',annotation='refGene',which='name',table_names=[]):
    """Return gene tables given gene lists"""
    
    gts=[]
    for i in range(0,len(gene_lists)):
        gt=GeneTable()
        gt.read(Host=Host,User=User,Db=Db,annotation=annotation,which=which,where=tuple(gene_lists[i]))
        try:
            gt.set_name(table_names[i])
        except IndexError:
            gt.set_name('Group %d' %i)
            
        gts.append(gt)
    
    return gts

def profile_tag_w_multiple_gene_sets(tag,gene_tables):
    """Return multiple profiles of give gene list in a single plot"""
        
    # tag profiler
    tp=TagProfiler()
    n_sets=len(gene_list)
    
    profiles={'avg_upstreams':[],'avg_downstreams':[],'avg_metagenes':[],'avg_metaexons':[],'avg_metaintrons':[],\
              'keys':[],'breaks':[],'metagene_breaks':[],'metaexon_breaks':[],'metaintron_breaks':[],'legends':[]}
    
    for gt in gene_tables:
        # get the profile
        keys,breaks,upstreams,downstreams,metagene_breaks,metagenes,metaexon_breaks,metaexons,metaintron_breaks,metaintrons=tp.do_individual_profiling(tag,gt,rel_pos=3000,metagenesize=3000,n=100,which='name')
        avg_upstream,avg_downstream,avg_metagene,avg_metaexon,avg_metaintron=tp.take_average_profiles(upstreams,downstreams,metagenes,metaexons,metaintrons)
        profiles['avg_upstreams'].append(avg_upstream)
        profiles['avg_downstreams'].append(avg_downstream)
        profiles['avg_metagenes'].append(avg_metagene)
        profiles['avg_metaexons'].append(avg_metaexon)
        profiles['avg_metaintrons'].append(avg_metaintron)
        profiles['keys'].append(keys)
        if not profiles['breaks']: profiles['breaks']+=breaks
        if not profiles['metagene_breaks']: profiles['metagene_breaks']+=[100.0*i/(len(avg_metagene)-1) for i in xrange(0,len(avg_metagene))]
        if not profiles['metaexon_breaks']: profiles['metaexon_breaks']+=[100.0*i/(len(avg_metaexon)-1) for i in xrange(0,len(avg_metaexon))]
        if not profiles['metaintron_breaks']: profiles['metaintron_breaks']+=[100.0*i/(len(avg_metaintron)-1) for i in xrange(0,len(avg_metaintron))]
    
    return profiles

def uppack_profiles(profiles):
    """Unpack multiple profiles"""
    
    breaks=profiles['breaks']
    avg_upstreams=profiles['avg_upstreams']
    avg_downstreams=profiles['avg_downstreams']
    metagene_breaks=profiles['metagene_breaks']
    avg_metagenes=profiles['avg_metagenes']
    metaexon_breaks=profiles['metaexon_breaks']
    avg_metaexons=profiles['avg_metaexons']
    metaintron_breaks=profiles['metaintron_breaks']
    avg_metaintrons=profiles['avg_metaintrons']
    keys=profiles['keys']
    legends=profiles['legends']
    
    
    return keys,breaks,upstreams,downstreams,metagene_breaks,metagenes,metaexon_breaks,metaexons,metaintron_breaks,metaintrons

def get_genelen_tagnum(gene_table,tag):
    """Return the length and number of tags within each gene"""
    
    names,chroms,tss,tts,exons,exone,strands=extract_gene_info(gene_table,tag,which='name')
    info=zip(names,chroms,tss,tts,exons,exone,strands)
    
    old_chrom=''
    for nm,chrom,ts,tt,exs,exe,strand in info:
        if chrom!=old_chrom:
            init=0
            old_chrom=chrom
        
        start,end=where(ts,tt,tag[chrom]['start'][init:])
        init=start
        print "%s\t%d\t%d\t%s\t%d" %(chrom,ts,tt,nm,end-start)

def extract_gene_info(gene_table,tag,which='name'):
    """Extract TSS,TTS and strand from the gene table"""
        
    # get the chromosomes
    chrs=set(gene_table.get_chroms()).intersection(tag.get_chroms())
    chrs=list(chrs)
        
    names=[]
    chroms=[]
    tss=[]
    tts=[]
    exons=[]
    exone=[]
    strand=[]
    for chr in chrs:
        names+=gene_table[chr][which]
        tss+=gene_table[chr]['txStart']
        tts+=gene_table[chr]['txEnd']
        strand+=gene_table[chr]['strand']
        chroms+=[chr]*len(gene_table[chr]['txStart'])
        exons+=gene_table[chr]['exonStarts']
        exone+=gene_table[chr]['exonEnds']

    return names,chroms,tss,tts,exons,exone,strand


def scan_scores_in_wig_w_options(bed, wig, score_type='max'):
    """Scan WIG to determine the score for each BED region
    
    Parameters:
    bed: a Bed object (in Cistrome.Assoc.inout.py)
    wig: a Wig object (in Cistrome.Assoc.inout.py)
    score_type: 'max' - set the score as the maximum value within a Bed region
                'min' - set the score as the minimum value within a Bed region
                'med' - set the score as the median value within a Bed region
                'mean' - set the score as the mean value within a Bed region
                
    WARNING! This function does not have a return value. the input bed will be modified with scores.
    """
    
    standard_chroms = {'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
    chroms = bed.get_chroms()
    chroms_wig = wig.get_chroms()
    if not chroms or not chroms_wig: return None
    
    for chrom in chroms:
        
        # get the standard chromosome name
        try:
            standardchrom=standard_chroms[chrom]
        except KeyError:
            standardchrom=chrom
            
        try:
            len_chrom = len(bed[chrom]['start'])
        except KeyError:    # might be an empty Bed
            continue
        
        if not bed[chrom].has_key('name'):  # if no names exist, put dummy names of 0s
            bed[chrom]['name'] = [0] *len_chrom
        
        if chrom in chroms_wig: # when chrom is in the wig
            init = 0
            scores = array('d', [])
            this_wig = wig[chrom]
            this_bed = bed[chrom]
            for start, end in itertools.izip(this_bed['start'], this_bed['end']):
                left, right = where(start, end, this_wig[0][init:])
                
                # if no region is found, just put 0 and go on
                if right == left: 
                    scores.append(0.0)
                    continue
                
                if score_type == 'max':
                    scores.append(max(this_wig[1][init+left:init+right]))
                elif score_type == 'min':
                    scores.append(min(this_wig[1][init+left:init+right]))
                elif score_type == 'med':
                    scores.append(median(this_wig[1][init+left:init+right]))
                elif score_type == 'mean':
                    scores.append(sum(this_wig[1][init+left:init+right])/(right-left))
                
                # update the init
                init += right
            
            # update scores
            bed[chrom]['score'] = scores[:]
    
    
def scan_scores_in_wig(bed, wig):
    """Scan WIG to determine the score for each BED region
    
    Parameters:
    bed: a Bed object (in Cistrome.Assoc.inout.py)
    wig: a Wig object (in Cistrome.Assoc.inout.py)
                
    Return:
    newbed: a new Bed object with scores
    """
    standard_chroms = {'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
    chroms = bed.get_chroms()
    chroms_wig = wig.get_chroms()
    if not chroms or not chroms_wig: return
    
    for chrom in chroms:
        
        # get the standard chromosome name
        try:
            standardchrom=standard_chroms[chrom]
        except KeyError:
            standardchrom=chrom
        
        try:
            len_chrom = len(bed[chrom]['start'])
        except KeyError:    # might be an empty Bed
            continue
        
        if not bed[chrom].has_key('name'):  # if no names exist, put dummy names of 0s
            bed[chrom]['name'] = [0] *len_chrom
        
        if not bed[chrom].has_key('score'): # if no scores exist, put dummy values of 0
              bed[chrom]['score'] = array('d', [0.0]*len_chrom) 
            
        if chrom in chroms_wig: # when chrom is in the wig
            init = 0
            bed[chrom]['score'] = array('d', [0.0]*len_chrom) 
            this_wig = wig[chrom]
            this_bed = bed[chrom]
            j = 0                                       # initialize jth ChIP region
            flag = False                                # flag, True: within a ChIP region, False: no
            for cor, val in itertools.izip(this_wig[0], this_wig[1]):
                
                # get the maximum value within a ChIP region as the score for the region
                # simple state machine: within a ChIP region or not
                if j < len_chrom:
                    if cor >= this_bed['start'][j] and cor < this_bed['end'][j]:
                        if flag:
                            maxval = max(maxval, val)
                        else:
                            flag = True
                            maxval = val    
                    elif cor >= this_bed['end'][j]:
                        if flag:
                            bed[chrom]['score'][j] = maxval
                            flag = False
                            j += 1
                        else:
                            j += 1
                        
                        # if the current coordinate falls in the next ChIP region
                        # if the ChIP regions have already been explored, just pass
                        try:
                            if cor >= this_bed['start'][j] and cor < this_bed['end'][j]:
                                flag = True
                                maxval = val
                        except IndexError:
                            pass     


def scan_scores_in_wig_and_sample_wig(bed, wig, wig_sample_resol):
    """Scan WIG to determine the score for each BED region
    
    Parameters:
    bed: a Bed object (in Cistrome.Assoc.inout.py)
    wig: a Wig object (in Cistrome.Assoc.inout.py)
    wig_sample_resol: resolution for wig sampling
                
    Return:
    newbed: a new Bed object with scores
    swig: sampled Wig object
    """
    standard_chroms = {'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
    chroms = bed.get_chroms()
    chroms_wig = wig.get_chroms()
    if not chroms or not chroms_wig: return
    
    swig = Wig()
    newbed = bed
    for chrom in chroms:
        
        # get the standard chromosome name
        try:
            standardchrom=standard_chroms[chrom]
        except KeyError:
            standardchrom=chrom
        
        try:
            len_chrom = len(newbed[chrom]['start'])
        except KeyError:    # might be an empty Bed
            continue
        
        if not newbed[chrom].has_key('name'):  # if no names exist, put dummy names of 0s
            bed[chrom]['name'] = [0] *len_chrom
        
        if not newbed[chrom].has_key('score'): # if no scores exist, put dummy values of 0
            newbed[chrom]['score'] = [0.0]*len_chrom   
            
        if chrom in chroms_wig: # when chrom is in the wig
            init = 0
            scores = []
            this_wig = wig[chrom]
            this_bed = newbed[chrom]
            j = 0                                       # initialize jth ChIP region
            flag = False                                # flag, True: within a ChIP region, False: no
            samp=[array('l',[]),array('d',[])]          # initialize the sampled wig
            for cor, val in itertools.izip(this_wig[0], this_wig[1]):
                
                # get the maximum value within a ChIP region as the score for the region
                # simple state machine: within a ChIP region or not
                if j < len_chrom:
                    if cor >= this_bed['start'][j] and cor < this_bed['end'][j]:
                        if flag:
                            maxval = max(maxval, val)
                        else:
                            flag = True
                            maxval = val    
                    elif cor >= this_bed['end'][j]:
                        if flag:
                            this_bed['score'][j] = maxval
                            flag = False
                            j += 1
                        else:
                            j += 1
                        
                        # if the current coordinate falls in the next ChIP region
                        # if the ChIP regions have already been explored, just pass
                        try:
                            if cor >= this_bed['start'][j] and cor < this_bed['end'][j]:
                                flag = True
                                maxval = val
                        except IndexError:
                            pass     
                
                # wig sampling
                sc = (cor/wig_sample_resol)*wig_sample_resol+1
                if len(samp[0])==0:
                    samp[0].append(sc)
                    samp[1].append(val)
                    continue
                if sc!=samp[0][-1]:
                    samp[0].append(sc)
                    samp[1].append(val)
            
            #added to sampWig only if there are some point(s)
            if samp[0]: swig.wig[standardchrom]=samp
        
    return newbed, swig
            