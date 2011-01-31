# Time-stamp: <2010-05-24 11:01:21 Tao Liu>

"""Module for PeakI, WigTrackI, FWTrackI and TrackI classes

Copyright (c) 2007,2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import os
import re
import sys
import time
import logging
from taolib.CoreLib.BasicStat.Func import mean,median,std
from array import array
from math import sqrt
from random import sample as random_sample
# ------------------------------------
# constants
# ------------------------------------
__version__ = "TabIO $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "PeakIO, WigTrackI, FWTrackI and TrackI classes"

# ------------------------------------
# Misc functions
# ------------------------------------

# to determine the byte size
if array('H',[1]).itemsize == 2:
    BYTE2 = 'H'
else:
    raise Exception("BYTE2 type cannot be determined!")

if array('I',[1]).itemsize == 4:
    BYTE4 = 'I'
elif array('L',[1]).itemsize == 4:
    BYTE4 = 'L'
else:
    raise Exception("BYTE4 type cannot be determined!")

if array('f',[1]).itemsize == 4:
    FBYTE4 = 'f'
elif array('d',[1]).itemsize == 4:
    FBYTE4 = 'd'
else:
    raise Exception("BYTE4 type cannot be determined!")

# ------------------------------------
# Classes
# ------------------------------------

class PeakIO:
    """IO for peak region information.

    This class can hold peak information from MAT/MA2C/MACS, and
    provide some extra functions:

    1. filtering functions, filter_pvalue/score/fdr/fold are functions
    to filter the peaks according to pvalue/score/fdr/fold range given
    by the user.

    2. overlapping functions

    """
    def __init__ (self, comment=""):
        """Initialization function.

        comment: you can add any comments to the peakIO object like
        whether or not it is from a ChIP-chip or a ChIP-seq
        experiments.
        """
        self.peaks = {}
        self.comment = comment

    def dup (self):
        """return a duplicate peakI.
        """
        r = PeakIO(comment=self.comment)
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[]
            new_peaks[chrom].extend(peaks[chrom])
        r.peaks = new_peaks
        return r
        
    
    def add (self, chromosome, start, end, summit=None, 
             score=None, total_p=None, 
             pvalue=None, fold_enrichment=None, fdr=None):
        """Use this function to add items to PeakIO object.

        items: (peak start,peak end, peak length, peak summit, peak
        score, number of tags/probes in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type

        Parameters:
        1. chromosome
        2. start
        3. end
        4. summit: the highest position for the peak region
        5. score:  the score for peak region
        6. total_p: total points in peak region. For ChIP-seq, it's
        how many tags in the region; for ChIP-chip, it's the number
        of probes.
        7. pvalue: -10*log(10,p-value) for peak region
        8. fold_enrichment: fold enrichment for the region
        9. fdr: False Discovery Rate for the region
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append((start,end,end-start,summit,
                                       score,total_p,
                                       pvalue,fold_enrichment,fdr))

    def filter_pvalue (self, pvalue_cut_low, pvalue_cut_up=None ):
        """Filter peaks in a given pvalue range.

        Note, pvalue is actually -10*log(10,pvalue)

        If pvalue_cut_low and pvalue_cut_up is assigned, the peaks with pvalue in [pvalue_cut_low,pvalue_cut_up).
        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if pvalue_cut_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[6] >= pvalue_cut_low and p[6]<pvalue_cut_up]
                if not new_peaks[chrom]: del new_peaks[chrom]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[6] >= pvalue_cut_low]
                if not new_peaks[chrom]: del new_peaks[chrom]
        self.peaks = new_peaks

    def filter_score (self, score_low, score_up=None ):
        """Filter peaks in a given score range.

        If score_low and score_up is assigned, the peaks with score in [score_low,score_up).
        
        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if score_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[4] >= score_low and p[4]<score_up]
                if not new_peaks[chrom]: del new_peaks[chrom]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[4] >= score_low]
                if not new_peaks[chrom]: del new_peaks[chrom]
        self.peaks = new_peaks

    def filter_fold (self, fold_low, fold_up=None ):
        """Filter peaks in a given fold enrichment range.

        If fold_low and fold_up is assigned, the peaks with fold in [fold_low,fold_up)
        
        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if fold_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[7] >= fold_low and p[7]<fold_up]
                if not new_peaks[chrom]: del new_peaks[chrom]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[7] >= fold_low]
                if not new_peaks[chrom]: del new_peaks[chrom]
        self.peaks = new_peaks

    def filter_fdr (self, fdr_up, fdr_low=None ):
        """Filter peaks in a given FDR range.

        If fdr_low and fdr_up is assigned, the peaks with fold in (fdr_low,fdr_up]. Otherwise, return the peaks with FDR lower or equal to fdr_up.
        
        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if fdr_low:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[8] > fdr_low and p[8]<=fdr_up]
                if not new_peaks[chrom]: del new_peaks[chrom]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[8] <= fdr_up]
                if not new_peaks[chrom]: del new_peaks[chrom]
        self.peaks = new_peaks

    def total (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x
    
    def ave_fold_enrichment (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        t = 0
        for chrom in chrs:
            x += len(peaks[chrom])
            for p in peaks[chrom]:
                t+=p[7]
        return t/x

    def max_fold_enrichment (self):
        """Return the maximum fold enrichment.
        
        """
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            if peaks[chrom]:
                m = max([i[7] for i in peaks[chrom]])
                if m>x:
                    x=m
        return x
        
    def max_min_score (self):
        """Return the maximum and minimum score.
        
        """
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        x = 0
        y = 100000
        for chrom in chrs:
            if peaks[chrom]:
                m = max([i[4] for i in peaks[chrom]])
                if m>x:
                    x=m
                m = min([i[4] for i in peaks[chrom]])
                if m<y:
                    y=m
        return (x,y)
    
    def tobed (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                text+= "%s\t%d\t%d\n" % (chrom,peak[0],peak[1])
        return text

    def towig (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                text+= "%s\t%d\t%d\n" % (peak[0],peak[1])
        return text

    def length (self):
        chrs = self.peaks.keys()
        chrs.sort()
        l = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                l+= peak[1]-peak[0]
        return l
        
        
    def init_from_dict (self, data):
        """Initialize the data from a dictionary. Improper assignment
        will damage the data structure.
        
        """
        self.peaks = {}
        chrs = data.keys()
        chrs.sort()
        for chrom in chrs:
            self.peaks[chrom]=[]
            a = self.peaks[chrom].append
            for i in data[chrom]:
                a(i)

    def sort (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            peaks[chrom].sort(lambda x,y:cmp(x[0],y[0]))


    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.peaks.keys())
        return l

    def add_other_peakI (self,peak2):
        """Add another peakI object to self. 
        
        """
        peaks1 = self.peaks
        peaks2 = peak2.peaks
        chrs1 = peaks1.keys()
        chrs1.extend(peaks2.keys())
        chrs = set(chrs1)
        for chrom in chrs:
            if not peaks1.has_key(chrom):
                peaks1[chrom]=[]
            if peaks2.has_key(chrom):
                peaks1[chrom].extend(peaks2[chrom])
        
    def merge_overlap ( self ):
        """peak_candidates[chrom] = [(peak_start,peak_end,peak_length,peak_summit,peak_height,number_cpr_tags)...]

        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[]
            n_append = new_peaks[chrom].append
            prev_peak = None
            peaks_chr = peaks[chrom]
            for i in xrange(len(peaks_chr)):
                if not prev_peak:
                    prev_peak = peaks_chr[i]
                    continue
                else:
                    if peaks_chr[i][0] <= prev_peak[1]:
                        s_new_peak = prev_peak[0]
                        e_new_peak = max(peaks_chr[i][1],prev_peak[1])
                        l_new_peak = e_new_peak-s_new_peak
                        if peaks_chr[i][4] > prev_peak[4]:
                            summit_new_peak = peaks_chr[i][3]
                            h_new_peak = peaks_chr[i][4]
                        else:
                            summit_new_peak = prev_peak[3]
                            h_new_peak = prev_peak[4]
                        prev_peak = (s_new_peak,e_new_peak,l_new_peak,summit_new_peak,h_new_peak,peaks_chr[i][5]+prev_peak[5])
                    else:
                        n_append(prev_peak)
                        prev_peak = peaks_chr[i]
            if prev_peak:
                n_append(prev_peak)
        del peaks
        self.peaks = new_peaks
        return True

    def overlap_with_other_peaks (self, peaks2, cover=0):
        """Peaks2 is a PeakIO object or dictionary with can be
        initialzed as a PeakIO. check __init__ for PeakIO for detail.

        return how many peaks are intersected by peaks2 by percentage
        coverage on peaks2(if 50%, cover = 0.5).
        """
        peaks1 = self.peaks
        if isinstance(peaks2,PeakIO):
            peaks2 = peaks2.peaks
        total_num = 0
        chrs1 = peaks1.keys()
        chrs2 = peaks2.keys()
        for k in chrs1:
            if not chrs2.count(k):
                continue
            rl1_k = iter(peaks1[k])
            rl2_k = iter(peaks2[k])
            tmp_n = False
            try:
                r1 = rl1_k.next()
                r2 = rl2_k.next()
                while (True):
                    if r2[0] < r1[1] and r1[0] < r2[1]:
                        a = sorted([r1[0],r1[1],r2[0],r2[1]])
                        if float(a[2]-a[1]+1)/r2[2] > cover:
                            if not tmp_n:
                                total_num+=1
                                tmp_n = True
                    if r1[1] < r2[1]:
                        r1 = rl1_k.next()
                        tmp_n = False
                    else:
                        r2 = rl2_k.next()
            except StopIteration:
                continue
        return total_num

    def overlap_with_FWTrackI (self, fwt, cover=0, twostrand=True):
        """
        
        """
        self_chrs = list(self.get_chr_names())
        if twostrand == False:
            fwt.self_merge_plus_minus_ranges_w_duplicates() # merge + and - strand
        fwt_chrs = list(fwt.get_chr_names())
        w = fwt.fw
        v = {}
        total_num = 0
        for k in self_chrs:
            if not fwt_chrs.count(k):
                continue
            # + strand fwt
            peak_k = iter(self.peaks[k])
            fwt_k = iter(fwt.get_ranges_by_chr(k)[0])
            tmp_n = False
            try:
                p1 = peak_k.next()
                f2b = fwt_k.next()
                f2e = f2b+w
                pt = k+":"+str(p1[0])+'...'+str(p1[1])
                tt = ""
                while (True):
                    if f2b < p1[1] and p1[0] < f2e:
                        a = sorted([p1[0],p1[1],f2b,f2e])
                        if float(a[2]-a[1]+1)/w > cover:
                            if not tmp_n:
                                total_num+=1
                                tmp_n = True
                        
                        tt += "\t"+str(f2b)+'...'+str(f2e)+"+"
                    if p1[1] < f2e:
                        if tt.find('\t') != -1:
                            v[pt]=tt
                        #print t
                        p1 = peak_k.next()
                        pt = k+":"+str(p1[0])+'...'+str(p1[1])
                        tt = ""
                        tmp_n = False
                    else:
                        f2b = fwt_k.next()
                        f2e = f2b + w
            except StopIteration:
                if tt and tt.find('\t') != -1:
                    v[pt]=tt
            # - strand fwt
            peak_k = iter(self.peaks[k])
            fwt_k = iter(fwt.get_ranges_by_chr(k)[1])
            tmp_n = False
            try:
                p1 = peak_k.next()
                f2e = fwt_k.next()
                f2b = f2e-w

                pt = k+":"+str(p1[0])+'...'+str(p1[1])
                tt = ""
                while (True):
                    if f2b < p1[1] and p1[0] < f2e:
                        a = sorted([p1[0],p1[1],f2b,f2e])
                        if float(a[2]-a[1]+1)/w > cover:
                            if not tmp_n:
                                total_num+=1
                                tmp_n = True
                
                        tt += "\t"+str(f2b)+'...'+str(f2e)+"-"
                    if p1[1] < f2e:
                        if tt.find('\t') != -1:
                            if v.has_key(pt):
                                v[pt]+=tt
                            else:
                                v[pt]=tt
                        #print t
                        p1 = peak_k.next()
                        pt = k+":"+str(p1[0])+'...'+str(p1[1])
                        tt = ""
                        tmp_n = False
                    else:
                        f2e = fwt_k.next()
                        f2b = f2e - w
            except StopIteration:
                if tt and tt.find('\t') != -1:
                    if v.has_key(pt):
                        v[pt]+=tt
                    else:
                        v[pt]=tt

                continue


        return (len(v.keys()),v)
        
    def extract_binkeepers (self,bks,func=max, NA=True):
        """Extract values from BinKeeper Object.

        Parameters:
        bks : BinKeeper dictionary
        func: a function which can operate list, default: max
        NA  : if True, missing value will be set to None, otherwise, skip it.
        Return:
        list of fun()ed values
        """
        d = []
        da = d.append
        chrs = set(self.peaks.keys())
        for chrom in chrs:
            if not bks.has_key(chrom):
                continue
            for cpeak in self.peaks[chrom]:
                try:
                    da (func(bks[chrom].pp2v(cpeak[0],cpeak[1])))
                except:
                    if NA:
                        da (None)
        return d

    def extract_binkeepers_pv (self,bks,NA=True):
        """Extract positions and values from BinKeeper Object.

        Parameters:
        bks : BinKeeper dictionary
        NA  : if True, missing value will be set to None, otherwise, skip it.
        Return:
        dictionary of (pos,value) for each chromosome
        """
        d = {}
        chrs = set(self.peaks.keys())
        for chrom in chrs:
            d[chrom]=[]
            da = d[chrom].append
            if not bks.has_key(chrom):
                continue
            for cpeak in self.peaks[chrom]:
                try:
                    da (bks[chrom].pp2pv(cpeak[0],cpeak[1]))
                except:
                    if NA:
                        da (None)
        return d



    def extract_wiggle_values (self,wigtrackI,func=max):
        """Extract values from wigtrackI.

        Parameters:
        wigtrackI : WigTrackI object
        func: a function which can operate list, default: max
        Return:
        list of fun()ed values
        """
        values = []
        vappend = values.append
        self.sort()
        chrs = set(self.peaks.keys())
        for chrom in chrs:
            cpeak = iter(self.peaks[chrom]).next
            cwig  = wigtrackI.get_data_by_chr(chrom)
            if not cwig:
                continue
            cpos = cwig[0]
            lpos = len(cpos)
            cpos = iter(cpos).next
            cvalue = iter(cwig[1]).next
            tmp=[]
            try:
                cur_peak = cpeak()
                for i in xrange(lpos):
                    p= cpos()
                    v= cvalue()
                    if p < cur_peak[0]:
                        continue
                    elif cur_peak[1] < p:
                        cur_peak = cpeak()
                        if tmp:
                            vappend(func(tmp))
                        tmp=[]
                    else:
                        tmp.append(v)
            except:
                if tmp:
                    vappend(func(tmp))
                    tmp = []
                continue
            if tmp:
                vappend(func(tmp))
        return values        

    def extract_wiggle_pv (self,wigtrackI):
        """Extract values from wigtrackI.

        Parameters:
        wigtrackI : WigTrackI object
        func: a function which can operate list, default: max
        Return:
        list of relative positions, and values
        """
        values = []
        vappend = values.append
        self.sort()
        chrs = set(self.peaks.keys())
        for chrom in chrs:
            cpeak = iter(self.peaks[chrom]).next
            cwig  = wigtrackI.get_data_by_chr(chrom)
            if not cwig:
                continue
            cpos = cwig[0]
            lpos = len(cpos)
            cpos = iter(cpos).next
            cvalue = iter(cwig[1]).next
            tmp=[]
            try:
                cur_peak = cpeak()
                for i in xrange(lpos):
                    p= cpos()
                    v= cvalue()
                    if p < cur_peak[0]:
                        continue
                    elif cur_peak[1] < p:
                        cur_peak = cpeak()
                        if tmp:
                            vappend(tmp)
                        tmp=[]
                    else:
                        tmp.append((p-cur_peak[0],v))
            except:
                if tmp:
                    vappend(tmp)
                    tmp = []
                continue
            if tmp:
                vappend(tmp)
        return values        

    def extract_wiggle_values_by_chrom (self,wigtrackI,func=max):
        """Extract values from wigtrackI.

        Parameters:
        wigtrackI : WigTrackI object
        func: a function which can operate list, default: max
        Return:
        Result will be stored in a dictionary whose keys are
        chromosome names.
        """
        values = {}
        self.sort()
        chrs = set(self.peaks.keys())
        for chrom in chrs:
            values[chrom]=[]
            vappend = values[chrom].append
            cpeak = iter(self.peaks[chrom]).next
            cwig  = wigtrackI.get_data_by_chr(chrom)
            if not cwig:
                continue
            cpos = cwig[0]
            lpos = len(cpos)
            cpos = iter(cpos).next
            cvalue = iter(cwig[1]).next
            tmp=[]
            try:
                cur_peak = cpeak()
                for i in xrange(lpos):
                    p= cpos()
                    v= cvalue()
                    if p < cur_peak[0]:
                        continue
                    elif cur_peak[1] < p:
                        cur_peak = cpeak()
                        #print tmp
                        #print max(tmp)
                        vappend(func(tmp))
                        tmp=[]
                        #vappend(v)
                    else:
                        tmp.append(v)
                        #vappend(v)
            except:
                if tmp:
                    #print max(tmp)
                    vappend(func(tmp))
                continue
        return values        

class WigTrackI:
    """Designed only for wig files generated by MACS/pMA2C/MAT(future
    version). The limitation is 'span' parameter for every track must
    be the same.
    
    """
    def __init__ (self):
        self.__data = {}
        self.span=0
        self.maxvalue =-10000
        self.minvalue = 10000

    def add_loc (self,chromosome,pos,value):
        if not self.__data.has_key(chromosome):
            self.__data[chromosome] = [array(BYTE4,[]),array(FBYTE4,[])] # for (pos,value)
        self.__data[chromosome][0].append(pos)
        self.__data[chromosome][1].append(value)
        if value > self.maxvalue:
            self.maxvalue = value
        if value < self.minvalue:
            self.minvalue = value
        
    def sort (self):
        """Naive sorting for tags. After sorting, counts are massed
        up.

        """
        for k in self.__data.keys():
            (p,v) = self.__data[k]
            pv = zip(p,v)
            pv = sorted(pv)
            self.__data[k] = [array(BYTE4,[]),array(FBYTE4,[])]
            pappend = self.__data[k][0].append
            vappend = self.__data[k][1].append
            for (tp,tv) in pv:
                pappend(tp)
                vappend(tv)

    def get_data_by_chr (self, chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([pos],[value])
        """
        if self.__data.has_key(chromosome):
            return self.__data[chromosome]
        else:
            return None
            #raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.__data.keys())
        return l

    def write_wig (self, fhd, name, shift=0):
        """Write all data to fhd in Wiggle Format.

        shift will be used to shift the coordinates. default: 0
        """
        fhd.write("track type=wiggle_0 name=\"%s\" description=\"%s\"\n" % (name,name))
        chrs = self.get_chr_names()
        for chrom in chrs:
            fhd.write("variableStep chrom=%s span=%d\n" % (chrom,self.span))
            (p,s) = self.__data[chrom]
            pnext = iter(p).next
            snext = iter(s).next            
            for i in xrange(len(p)):
                pos = pnext()
                score = snext()
                fhd.write("%d\t%.4f\n" % (pos+shift,score))                

    def filter_score (self, cutoff=0):
        """Filter using a score cutoff. Return a new WigTrackI.
        
        """
        ret = WigTrackI()
        ret.span = self.span
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]

            ret.__data[chrom] = [array(BYTE4,[]),array(FBYTE4,[])]
            (np,ns) = ret.__data[chrom]
            npa = np.append
            nsa = ns.append

            pnext = iter(p).next
            snext = iter(s).next            
            for i in xrange(len(p)):
                pos = pnext()
                score = snext()
                if score > cutoff:
                    npa(pos)
                    nsa(score)
        return ret

    def filter_score_below (self, cutoff=0):
        """Keep points below a score cutoff. Return a new WigTrackI.
        
        """
        ret = WigTrackI()
        ret.span = self.span
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]

            ret.__data[chrom] = [array(BYTE4,[]),array(FBYTE4,[])]
            (np,ns) = ret.__data[chrom]
            npa = np.append
            nsa = ns.append

            pnext = iter(p).next
            snext = iter(s).next            
            for i in xrange(len(p)):
                pos = pnext()
                score = snext()
                if score < cutoff:
                    npa(pos)
                    nsa(score)
        return ret

    def write_gff (self, fhd, shift=0, source=".", feature="."):
        """Write all data to fhd in GFF format.

        shift will be used to shift the coordinates. default: 0
        """
        assert isinstance(fhd,file)
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                fhd.write(
                    "\t".join( (chrom,source,feature,
                                str(pi-shift),str(pi-shift+self.span-1),
                                str(si),'+','.'
                                ) )+"\n"
                    )

    def write_bed (self, fhd):
        """Write all data to fhd in BED format.
        
        """
        pass

    def remove_redundant (self):
        """Remove redundant position, keep the highest score.
        
        """
        chrs = set(self.__data.keys())
        ndata = {}
        for chrom in chrs:
            ndata[chrom] = [array(BYTE4,[]),array(FBYTE4,[])] # for (pos,value)
            nd_p_append = ndata[chrom][0].append
            nd_s_append = ndata[chrom][1].append
            (p,s) = self.__data[chrom]
            prev_p = None
            prev_s = None
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                if not prev_p:
                    prev_p = pi
                    prev_s = si
                else:
                    if pi == prev_p:
                        if si>prev_s:
                            prev_s = si
                    else:
                       nd_p_append (prev_p)
                       nd_s_append (prev_s)
            nd_p_append (prev_p)
            nd_s_append (prev_s)
        del self.__data
        self.__data = ndata

    def find_peaks (self, bw=None):
        """A naive peak finding algorithm to find all local maximum
        points.

        bw : Bandwidth. Two local maximums will compete if their
        distance is less than bw. The higher one will be kept, or if
        they are equal, the last point will be kept. Default is 4*span.

        Return type:
        
        find_peak will return a new WigTrackI object with position
        and score for each local maximum(peaks).
        """
        if not bw:
            bw = 4*self.span
        bw = max(1,bw)
        ret = WigTrackI()
        ret.span = 1
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            prev_peak_p = -10000
            prev_peak_s = -10000
            prev_peak_summits = []
            prev_p = -10000
            prev_s = -10000
            increase = False
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                if si > prev_s:
                    # increase
                    increase = True
                elif si < prev_s:
                    # decrease
                    if increase:
                        # prev_p is a summit
                        #print prev_p,prev_s,pi,si
                        if prev_p-prev_peak_p < bw:
                            # two summits are near
                            if prev_s > prev_peak_s:
                                # new summit is high
                                prev_peak_s = prev_s
                                prev_peak_p = prev_p
                        else:
                            # it's a new summit
                            #print chrom,prev_peak_p,prev_peak_s
                            try:
                                ret.add_loc(chrom,prev_peak_p,prev_peak_s)
                            except OverflowError:
                                pass
                            prev_peak_s = prev_s
                            prev_peak_p = prev_p
                            
                    increase = False
                prev_p = pi
                prev_s = si
            try:
                ret.add_loc(chrom,prev_peak_p,prev_peak_s)            
            except OverflowError:
                pass
        return ret

    def find_valleys (self, bw=None):
        """A naive peak finding algorithm to find all local minimum
        points.

        bw : Bandwidth. Two local maximums will compete if their
        distance is less than bw. The higher one will be kept, or if
        they are equal, the last point will be kept. Default is 4*span.

        Return type:
        
        find_peak will return a new WigTrackI object with position
        and score for each local maximum(peaks).
        """
        if not bw:
            bw = 4*self.span
        bw = max(1,bw)
        ret = WigTrackI()
        ret.span = 1
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            prev_peak_p = -10000
            prev_peak_s = -10000
            prev_peak_summits = []
            prev_p = -10000
            prev_s = -10000
            decrease = False
            for i in xrange(len(p)):
                pi = p[i]
                si = s[i]
                if si < prev_s:
                    # decrease
                    decrease = True
                elif si > prev_s:
                    # increase
                    if decrease:
                        # prev_p is a valley
                        #print prev_p,prev_s,pi,si
                        if prev_p-prev_peak_p < bw:
                            # two summits are near
                            if prev_s < prev_peak_s:
                                # new summit is lower
                                prev_peak_s = prev_s
                                prev_peak_p = prev_p
                        else:
                            # it's a new valley
                            #print chrom,prev_peak_p,prev_peak_s
                            try:
                                ret.add_loc(chrom,prev_peak_p,prev_peak_s)
                            except OverflowError:
                                pass
                            prev_peak_s = prev_s
                            prev_peak_p = prev_p
                            
                    decrease = False
                prev_p = pi
                prev_s = si
            try:
                ret.add_loc(chrom,prev_peak_p,prev_peak_s)            
            except OverflowError:
                pass
        return ret

    def summary (self):
        """Calculate the sum, max, min, mean, and std. Return a tuple for (sum, max, min, mean, std).
        
        """
        n_v = 0
        sum_v = 0
        max_v = -100000
        min_v = 100000
        for (p,v) in self.__data.values():
            sum_v += sum(v)
            n_v += len(v)
            max_v = max(max(v),max_v)
            min_v = min(min(v),min_v)
        mean_v = float(sum_v)/n_v
        variance = 0.0
        for (p,v) in self.__data.values():
            for vv in v:
                tmp = vv-mean_v
                variance += tmp*tmp
        variance /= float(n_v-1)
        std_v = sqrt(variance)
        return (sum_v, max_v, min_v, mean_v, std_v)

    def null_model_summary (self, sample=10):
        """Calculate the sum, max, min, mean, and std. Return a tuple for (sum, max, min, mean, std).

        This is for NULL model which is a symetric normal distribution
        abased on sample of the whole data set.
        """
        # sample data step
        data_step = int(100/sample)

        null_list = array(FBYTE4,[])
        na = null_list.append
        for (p,v) in self.__data.values():
            i = 0
            for vv in v:
                i+=1
                if i==data_step:
                    na(vv)
                    i=0
                
        
        sum_v = sum(null_list)
        mean_v = sum_v/float(len(null_list))

        null_list_len = len(null_list)
        null_list = sorted(null_list)
        median_index1 = (null_list_len - 1) / 2
        median_index2 = null_list_len / 2
        median_v = (null_list[median_index1]+null_list[median_index2])/2.0

        # make right part of nullList

        for i in xrange(null_list_len/2):
            null_list[null_list_len-i-1] = 2* median_v - null_list[i]
        
        std_v = std(null_list)

        return (sum_v,max(null_list),min(null_list),median_v,std_v)

    def normalize (self,null=False,sample_percent=10):
        """Normalize values centered at 0 and variance as 1.

        If null is True, it will use the null list to calculate mean and std.
        When null is True, sample_percent will be passed to null_model to sample the data.
        """
        if null:
            (sum_v,max_v,min_v,mean_v,std_v) = self.null_model_summary(sample=sample_percent)
        else:
            (sum_v,max_v,min_v,mean_v,std_v) = self.summary()
        for (p,v) in self.__data.values():
            for i in range(len(v)):
                v[i] = float(v[i]-mean_v)/std_v
        return (sum_v, max_v, min_v, mean_v, std_v)
                

    def call_peaks (self, cutoff=1, up_limit=1e310, min_length=200, max_gap=50):
        """This function try to find some region within which, scores
        are continuously higher than a cutoff.

        cutoff:  cutoff of value, default 1
        min_length :  minimum peak length, default 200
        gap   :  maximum gap to merge nearby peaks
        """
        chrs = self.get_chr_names()
        peaks = PeakIO()                      # dictionary to save peaks
        for chrom in chrs:
            (ps,ss) = self.get_data_by_chr(chrom)
            psn = iter(ps).next         # assign the next function to a virable to speed up
            ssn = iter(ss).next
            x = 0
            while True:
                # find the first point above cutoff
                try:
                    p = psn()
                    s = ssn()
                except:
                    break
                x += 1                  # index for the next point
                if s >= cutoff and s<=up_limit:
                    peak_content = [(p,s),]
                    break               # found the first point above cutoff

            for i in range(x,len(ps)):
                p = psn()
                s = ssn()
                if s < cutoff or s > up_limit:
                    continue
                # for points above cutoff
                if p - peak_content[-1][0] <= max_gap:
                    peak_content.append((p,s))
                else:
                    # a new peak
                    peak_length = peak_content[-1][0]-peak_content[0][0]+self.span
                    if peak_length >= min_length:
                        summit = None
                        summit_score = None
                        for (m,n) in peak_content:
                            if not summit_score or summit_score < n:
                                summit = m
                                summit_score = n
                        peaks.add(chrom,peak_content[0][0],peak_content[-1][0]+self.span,
                                  summit=summit-peak_content[0][0],score=summit_score)
                    peak_content = [(p,s),]
        return peaks

    def total (self):
        t = 0
        chrs = set(self.__data.keys())
        for chrom in chrs:
            (p,s) = self.__data[chrom]
            t += len(p)
        return t
    
class FWTrackI:
    """Fixed Width Ranges along the whole genome (commonly with the
    same annotation type), which are stored in a dict.

    Ranges are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.

    """
    def __init__ (self,fw=0,anno=""):
        """fw is the fixed-width for all ranges.

        some member values you can access:
        1. self.fw : the fixed-width for all the positions
        2. self.total: total tags in this track
        3. self.total_unique: total unique tags in this track,
        available after self_merge_overlap()
        4. self.annotation: the annotation
        """
        self.fw = fw
        self.__ranges = {}
        self.__counts = {}    # store the counts if a position
        self.__well_merged = False
        self.total = 0                  # total tags
        self.total_unique = 0		# total unique tags
        self.annotation = anno          # need to be figured out

    def add_loc (self, chromosome, fiveendpos, strand):
        """Add a range to the list according to the sequence name.
        
        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, neg for neg strand
        strand     -- 0: plus, 1: minus
        """
        if not self.__ranges.has_key(chromosome):
            self.__ranges[chromosome] = [array(BYTE4,[]),array(BYTE4,[])] # for (+strand, -strand)
            self.__counts[chromosome] = [array(BYTE2,[]),array(BYTE2,[])] # for (+,-)
        self.__ranges[chromosome][strand].append(fiveendpos)
        self.__counts[chromosome][strand].append(1)
        self.total+=1

    def get_ranges_by_chr (self, chromosome):
        """Return array of locations by chromosome.

        The return value is a tuple:
        ([array for plus strand],[array for minus strand])
        """
        if self.__ranges.has_key(chromosome):
            return self.__ranges[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_counts_by_chr (self, chromosome):
        """Return array of counts by chromosome.

        The return value is a tuple:
        ([array for plus strand],[array for minus strand])
        """
        if self.__counts.has_key(chromosome):
            return self.__counts[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.__ranges.keys())
        return l

    def length (self):
        """Total covered length = total number of tags * width of tag		

        """
        return self.total*self.fw

    def sort (self):
        """Naive sorting for tags. After sorting, counts are massed
        up.

        Note: counts are massed up, so they will be set to 1 automatically.
        """
        for k in self.__ranges.keys():
            self.__ranges[k][0] = array(BYTE4,sorted(self.__ranges[k][0]))
            self.__ranges[k][1] = array(BYTE4,sorted(self.__ranges[k][1]))
            self.__counts[k][0] = array(BYTE2,[1])*len(self.__ranges[k][0])
            self.__counts[k][1] = array(BYTE2,[1])*len(self.__ranges[k][1])
            

    def merge_overlap (self, sec=None):
        """merge the SAME positions. Record the duplicate number in
        self.__counts{}.

        If sec is non-None, the second FWTrackI will be added before
        merging.
        
        *Note: different with the merge_overlap() in TrackI class,
        which merges the overlapped ranges.
        *Note: plus and minus strands are processed separately.
        """
        self.total = 0
        self.total_unique = 0

        if sec:
            if not isinstance(sec,FWTrackI):
                raise Exception("sec must be FWTrackI object!")
            else:
                secchrs= sec.get_chr_names()
                for k in self.__ranges.keys(): # for each chromosome
                    secr = sec.get_ranges_by_chr(k)
                    secc = sec.get_counts_by_chr(k)
                    if k in secchrs:
                        self.__ranges[k][0] += secr[0]
                        self.__counts[k][1] += secr[1]
                        self.__ranges[k][0] += secc[0]
                        self.__counts[k][1] += secc[1]                        

        for k in self.__ranges.keys(): # for each chromosome
            # + strand
            plus = sorted(self.__ranges[k][0])
            if len(plus) <1:
                logging.warning("NO records for chromosome %s, plus strand!" % (k))
                new_plus = []
                new_plus_c = []
            else:
                plus_c = self.__counts[k][0]
                (new_plus,new_plus_c) = (array(BYTE4,[plus[0]]),array(BYTE2,[1]))
            
                pappend = new_plus.append
                pcappend = new_plus_c.append
                n = 0                # the position in new list
                for p in plus[1:]:
                    if p == new_plus[n]:
                        try:
                            new_plus_c[n]+=1
                        except OverflowError:
                            logging.warning("> 32767 + strand tags mapped to position %d on chromosome %s!" % (p,k))
                            new_plus_c[n]=32767

                    else:
                        pappend(p)
                        pcappend(1)
                        n += 1
                self.total_unique +=  len(new_plus)
                self.total += sum(new_plus_c)
            # - strand
            minus = sorted(self.__ranges[k][1])
            if len(minus) <1:
                logging.warning("NO records for chromosome %s, minus strand!" % (k))
                new_minus = []
                new_minus_c = []
            else:
                minus_c = self.__counts[k][0]
                (new_minus,new_minus_c) = (array(BYTE4,[minus[0]]),array(BYTE2,[1]))
            
                mappend = new_minus.append
                mcappend = new_minus_c.append
                n = 0                # the position in new list
                for p in minus[1:]:
                    if p == new_minus[n]:
                        try:
                            new_minus_c[n]+=1
                        except OverflowError:
                            logging.warning("> 32767 - strand tags mapped to position %d on chromosome %s!" % (p,k))
                            new_minus_c[n]=32767
                    else:
                        mappend(p)
                        mcappend(1)
                        n += 1
                self.total_unique +=  len(new_minus)
                self.total += sum(new_minus_c)

            self.__ranges[k]=[new_plus,new_minus]
            self.__counts[k]=[new_plus_c,new_minus_c]
            self.__well_merged = True
		
    def self_merge_plus_minus_ranges_w_duplicates (self):
        """Merge minus strand ranges to plus strand. The duplications
        on a single strand is erased. But if the same position is on
        both pos and neg strand, keep them both.
        
        Side effect: Reset the counts. self.total_unique is set to
        None, use self.total instead.
        """
        self.total_unique = None
        self.total = 0
        for chrom in self.__ranges.keys():
            (plus_tags,minus_tags) = self.__ranges[chrom]
            new_plus_tags = array(BYTE4,[])
            #reset counts
            self.__counts[chrom][0] = array(BYTE2,[])
            self.__counts[chrom][1] = array(BYTE2,[])
            ip = 0
            im = 0
            lenp = len(plus_tags)
            lenm = len(minus_tags)
            while ip < lenp and im < lenm:
                if plus_tags[ip] < minus_tags[im]:
                    new_plus_tags.append(plus_tags[ip])
                    ip += 1
                else:
                    new_plus_tags.append(minus_tags[im])
                    im += 1
            if im < lenm:
                # add rest of minus tags
                new_plus_tags.extend(minus_tags[im:])
            if ip < lenp:
                # add rest of plus tags
                new_plus_tags.extend(plus_tags[ip:])

            self.__ranges[chrom] = [new_plus_tags,[]]
            self.total += len(new_plus_tags)

    def sample (self, percent):
        """Sample the tags for a given percentage.
        
        Side effect: self.total_unique is set to None, and counts are unset.
        """
        self.total = 0
        self.total_unique = None
        for key in self.__ranges.keys():
            num = int(len(self.__ranges[key][0])*percent)
            self.__ranges[key][0]=array(BYTE4,sorted(random_sample(self.__ranges[key][0],num)))
            num = int(len(self.__ranges[key][1])*percent)
            self.__ranges[key][1]=array(BYTE4,sorted(random_sample(self.__ranges[key][1],num)))
            self.total += len(self.__ranges[key][0]) + len(self.__ranges[key][1])
            self.__counts[key] = [[],[]]
            
    def __str__ (self):
        return self.__to_wiggle()
        
    def __to_wiggle (self):
        t = "track type=wiggle_0 name=\"tag list\" description=\"%s\"\n" % (self.annotation)
        for key in self.__ranges.keys():
            if self.__ranges[key][0]:
                t += "variableStep chrom=%s span=%d strand=0\n" % (key,self.fw)
                for i in self.__ranges[key][0]:
                    t += str(i)+"\n"
            if self.__ranges[key][1]:
                t += "variableStep chrom=%s span=%d strand=1\n" % (key,self.fw)
                for i in self.__ranges[key][1]:
                    t += str(i)+"\n"
        return t

class TrackI:
    """Track Class, container of regions with certain chromosome,
    start and end position. Region is strandless. For Track for
    regions with strand, use TrackII. For maximum performance of
    searching, feed TrackI to BinKeeper then use BinKeeper functions.

    Basic operations:
    1. add, add a region to TrackI
    2. init_from_dict, initialize trackI with a existing dictionary
    3. sort, sort regions in ever chromosome. It must be called before merge_overlap.
    4. clean_overlap, merge overlapped regions. It muse follow a self.sort() call.

    Operation provided by this class:
    1. dup,   duplicate self
    2. merge, merge another TrackI instance.

    """
    def __init__ (self, comment=""):
        """Initialization function.

        The instance is EMPTY after initialization. You need to add
        something.
        """
        self.regions = {}
        self.comment = comment

    def dup (self):
        """Return a duplicate TrackI.
        
        """
        r = TrackI(comment=self.comment)
        regions = self.regions
        new_regions = {}
        chrs = regions.keys()
        chrs.sort()
        for chrom in chrs:
            new_regions[chrom]=[]
            new_regions[chrom].extend(regions[chrom])
        r.regions = new_regions
        return r
    
    def add (self, chromosome, start, end):
        """Use this function to add region to TrackI instance.

        items: (peak start,peak end, peak length, peak summit, peak
        score, number of tags/probes in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type

        :parameter chromosome: string
        :parameter start: integer
        :parameter end: integer
        """
        if not self.regions.has_key(chromosome):
            self.regions[chromosome]=[]
        self.regions[chromosome].append((start,end,end-start))

    def total (self):
        """Return the total number of regions in TrackI instance.
        
        """
        regions = self.regions
        chrs = regions.keys()
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(regions[chrom])
        return x
    
    def tobed (self):
        """Return all region in BED format text string.

        You can save the return value to a BED file.        
        """
        text = ""
        chrs = self.regions.keys()
        chrs.sort()
        for chrom in chrs:
            for region in self.regions[chrom]:
                text+= "%s\t%d\t%d\n" % (chrom,region[0],region[1])
        return text

    def length (self):
        """Return the total lengt of TrackI object.

        Summarize all regions, regardless of overlap.
        """
        chrs = self.regions.keys()
        chrs.sort()
        l = 0
        for chrom in chrs:
            for region in self.regions[chrom]:
                l+= region[1]-region[0]
        return l
        
        
    def init_from_dict (self, data):
        """Initialize the data from a dictionary. Improper assignment
        will damage the data structure. Not recommended, use
        self.add() instead.

        :parameter data: a dictionary using chromosome names as keys,
        and list of tuples of (start, end position) as values.
        """
        self.regions = {}
        chrs = data.keys()
        chrs.sort()
        for chrom in chrs:
            self.regions[chrom]=[]
            a = self.regions[chrom].append
            for i in data[chrom]:
                a(i)

    def sort (self):
        """Sort regions by their start positions for each chromosome.

        """
        regions = self.regions
        chrs = regions.keys()
        chrs.sort()
        for chrom in chrs:
            regions[chrom].sort(lambda x,y:cmp(x[0],y[0]))


    def get_chr_names (self):
        """Return all the chromosome names stored.
        
        """
        l = set(self.regions.keys())
        return l

    def merge (self,track2):
        """Merge another TrackI object to self. 

        After merging, the instance need to be sorted!
        """
        regions1 = self.regions
        regions2 = track2.regions
        chrs1 = regions1.keys()
        chrs1.extend(regions2.keys())
        chrs = set(chrs1)
        for chrom in chrs:
            if not regions1.has_key(chrom):
                regions1[chrom]=[]
            if regions2.has_key(chrom):
                regions1[chrom].extend(regions2[chrom])
        
    def clean_overlap ( self ):
        """Merge the overlap regions for each chromosome.

        This function can only be called after self.sort()! Otherwise,
        the data will be messed up.
        """
        regions = self.regions
        new_regions = {}
        chrs = regions.keys()
        chrs.sort()
        for chrom in chrs:
            new_regions[chrom]=[]
            n_append = new_regions[chrom].append
            prev_region = None
            regions_chr = regions[chrom]
            for i in xrange(len(regions_chr)):
                if not prev_region:
                    prev_region = regions_chr[i]
                    continue
                else:
                    if regions_chr[i][0] <= prev_region[1]:
                        s_new_region = prev_region[0]
                        e_new_region = max(regions_chr[i][1],prev_region[1])
                        l_new_region = e_new_region-s_new_region
                        prev_region = (s_new_region,e_new_region)
                    else:
                        n_append(prev_region)
                        prev_region = regions_chr[i]
            if prev_region:
                n_append(prev_region)
        del regions
        self.regions = new_regions
        return True

    def intersect_with_other (self, track2, cover=0):
        """track2 is a TrackI instance or dictionary with can be
        initialzed as a TrackI.

        return how many regions are intersected by track2 by percentage
        coverage on track2(if cover = 0.5, percentage >50%).
        """
        regions1 = self.regions
        if isinstance(track2,TrackI):
            regions2 = track2.regions
        else:
            regions = track2
        total_num = 0
        chrs1 = regions1.keys()
        chrs2 = regions2.keys()
        for k in chrs1:
            if not chrs2.count(k):
                continue
            rl1_k = iter(regions1[k])
            rl2_k = iter(regions2[k])
            tmp_n = False
            try:
                r1 = rl1_k.next()
                r2 = rl2_k.next()
                while (True):
                    if r2[0] < r1[1] and r1[0] < r2[1]:
                        a = sorted([r1[0],r1[1],r2[0],r2[1]])
                        if float(a[2]-a[1]+1)/r2[2] > cover:
                            if not tmp_n:
                                total_num+=1
                                tmp_n = True
                    if r1[1] < r2[1]:
                        r1 = rl1_k.next()
                        tmp_n = False
                    else:
                        r2 = rl2_k.next()
            except StopIteration:
                continue
        return total_num

    def extract_binkeepers (self,bks,func=max, NA=True):
        """Extract values from BinKeeper Object within contained regions.

        Parameters:
        :parameter bks : BinKeeper dictionary
        :parameter func: a function which can operate list, default: max
        :parameter NA  : if True, missing value will be set to None, otherwise, skip it.
        :retval: list of fun()ed values
        """
        d = []
        da = d.append
        chrs = set(self.regions.keys())
        for chrom in chrs:
            if not bks.has_key(chrom):
                continue
            for cregion in self.regions[chrom]:
                try:
                    da (func(bks[chrom].pp2v(cregion[0],cregion[1])))
                except:
                    if NA:
                        da (None)
        return d

#     def extract_DBBinKeeperI (self,dirname,func=max, NA=True):
#         """Extract values from BinKeeper Object.

#         Parameters:
#         dirname : DBBinKeeper directory
#         func: a function which can operate list, default: max
#         NA  : if True, missing value will be set to None, otherwise, skip it.
#         Return:
#         list of fun()ed values
#         """
#         d = []
#         da = d.append
#         chrs = set(self.peaks.keys())
#         for chrom in chrs:
#             dbfile = os.path.join(dirname,chrom+".db")
#             if not os.path.exists(dbfile):
#                 continue
#             t0=time.time()
#             dbbk = DBBinKeeperI(dbfile,chromosome=chrom,bin=8,chromosomesize=250000000)
#             print time.time()-t0
#             t0=time.time()
#             for cpeak in self.peaks[chrom]:
#                 l = dbbk.pp2v(cpeak[0],cpeak[1])
#                 #print "l:",l
#                 try:
#                     da (func(l))
#                 except:
#                     if NA:
#                         da (None)
#             print time.time()-t0
#         #print d
#         return d

    def extract_binkeepers_pv (self,bks,NA=True):
        """Extract positions and values from BinKeeper Object.

        :parameter bks : BinKeeper dictionary
        :parameter NA  : if True, missing value will be set to None, otherwise, skip it.
        :retval: dictionary of (pos,value) for each chromosome
        """
        d = {}
        chrs = set(self.regions.keys())
        for chrom in chrs:
            d[chrom]=[]
            da = d[chrom].append
            if not bks.has_key(chrom):
                continue
            for cregion in self.regions[chrom]:
                try:
                    da (bks[chrom].pp2pv(cregion[0],cregion[1]))
                except:
                    if NA:
                        da (None)
        return d

    def extract_wiggle_values (self,wigtrackI,func=max):
        """Extract values from wigtrackI.
        
        :parameter wigtrackI : WigTrackI object
        :parameter func: a function which can operate list, default: max
        :retval: list of fun()ed values
        """
        values = []
        vappend = values.append
        self.sort()
        chrs = set(self.regions.keys())
        for chrom in chrs:
            cregion = iter(self.regions[chrom]).next
            cwig  = wigtrackI.get_data_by_chr(chrom)
            if not cwig:
                continue
            cpos = cwig[0]
            lpos = len(cpos)
            cpos = iter(cpos).next
            cvalue = iter(cwig[1]).next
            tmp=[]
            try:
                cur_region = cregion()
                for i in xrange(lpos):
                    p= cpos()
                    v= cvalue()
                    if p < cur_region[0]:
                        continue
                    elif cur_region[1] < p:
                        cur_region = cregion()
                        if tmp:
                            vappend(func(tmp))
                        tmp=[]
                    else:
                        tmp.append(v)
            except:
                if tmp:
                    vappend(func(tmp))
                    tmp = []
                continue
            if tmp:
                vappend(func(tmp))
        return values        

    def extract_wiggle_values_by_chrom (self,wigtrackI,func=max):
        """Extract values from wigtrackI.

        :parameter wigtrackI : WigTrackI object
        :parameter func: a function which can operate list, default: max
        :retval: Result will be stored in a dictionary whose keys are
        chromosome names.
        """
        values = {}
        self.sort()
        chrs = set(self.regions.keys())
        for chrom in chrs:
            values[chrom]=[]
            vappend = values[chrom].append
            cregion = iter(self.regions[chrom]).next
            cwig  = wigtrackI.get_data_by_chr(chrom)
            if not cwig:
                continue
            cpos = cwig[0]
            lpos = len(cpos)
            cpos = iter(cpos).next
            cvalue = iter(cwig[1]).next
            tmp=[]
            try:
                cur_region = cregion()
                for i in xrange(lpos):
                    p= cpos()
                    v= cvalue()
                    if p < cur_region[0]:
                        continue
                    elif cur_region[1] < p:
                        cur_region = cregion()
                        #print tmp
                        #print max(tmp)
                        vappend(func(tmp))
                        tmp=[]
                        #vappend(v)
                    else:
                        tmp.append(v)
                        #vappend(v)
            except:
                if tmp:
                    #print max(tmp)
                    vappend(func(tmp))
                continue
        return values        
