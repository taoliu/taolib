
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
import sys
from array import *
from bisect import *

# ------------------------------------
# My own Python modules
# ------------------------------------
from Cistrome.Assoc.inout import *

#-------------------------------------
# classes
#-------------------------------------  

class Sampler:
    """Abstract class Sampler"""
    
    def __init__(self,name=''):
        self.__name__=name
         # for C elegans
        self.standard_chroms={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}

    def sample(self):
        """Doing sampling"""
        
        pass
    
    def get_name(self):
        """Return the name of the sampler"""
        
        return self.__name__
    
    def set_name(self,name):
        """Set the name of the sampler"""
        
        self.__name__=name
        

class GenomeSampler(Sampler):
    
    def __init__(self,name=''):
        """Constructor"""
        
        Sampler.__init__(self,name)
        
    def sample(self,wig=None,resolution=100):
        """Doing sampling"""
        
        try:
            chroms=wig.get_chroms()
        except AttributeError:
            raise Exception("Argument 'wig' must be given")
        
        coordinates={}
        for chrom in chroms:
            try:
                standardchrom=self.standard_chroms[chrom]
            except KeyError:
                standardchrom=chrom
            wigcoord=wig[chrom][0]
            coordinates[standardchrom]=[]
            for wc in wigcoord:
                coordinate=(int(round(1.0*wc/resolution)))*resolution+1
                
                if not coordinates[standardchrom] or coordinate!=coordinates[standardchrom][-1]: 
                    coordinates[standardchrom].append(coordinate)
        
        return coordinates
  
        
class ChIPSampler(Sampler):
    """Class ChIPSampler"""
    
    def __init__(self,name=''):
        """Constructor"""
        
        Sampler.__init__(self,name='')
        # for C elegans
        self.standard_chroms={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
   
    def sample(self,bed=None,resolution=600):
        """Doing sampling"""
        
        try:
            chroms=bed.get_chroms()
        except AttributeError:
            raise Exception("Argument 'bed' must be given")
        
        coordinates={}
        for chrom in chroms:
            try:
                standardchrom=self.standard_chroms[chrom]
            except KeyError:
                standardchrom=chrom
            coordinates[standardchrom]=[]
            ChIP=zip(bed[chrom]['start'],bed[chrom]['end'],map(lambda x,y,:(x+y)/2,bed[chrom]['start'],bed[chrom]['end']))
            howmanyChIPs=len(ChIP)
            for i in xrange(0,howmanyChIPs):
                
                # get the begining, end, and center of a peak
                beg,end,center=ChIP[i]

                ### When using the real locations of beg
                Ns=range(center,max(0,beg-1),-1*resolution)
                Ns.reverse()
                Ns+=range(center+resolution,end+1,resolution)
                
                if Ns: coordinates[standardchrom].extend(Ns)
            
            coordinates[standardchrom].sort()
        
        return coordinates
    
class WigSampler(Sampler):
    """Class WigSampler
    
    This class samples a wig file (object)
    """

    def __init__(self,name=''):
        """Constructor"""
        
        Sampler.__init__(self,name='')
        
    def sample(self, wig, resolution):
        """Sample a wig file at the given resolution. 
        
        Parameters:
        1. wig: a wig object (see inout.py)
        2. resolution: sampling resolution
        
        Return:
        sampWig: the sampled wig object
        
        """
        
        # parameter checking
            
        try:
            chroms=wig.get_chroms()
        except AttributeError:
            raise Exception("Argument 'wig' must be given")
        
        sampWig=Wig()
        for chrom in chroms:
            try:
                standardchrom=self.standard_chroms[chrom]
            except KeyError:
                standardchrom=chrom
            
            samp=[array('l',[]),array('d',[])]
            for wc,val in itertools.izip(wig[chrom][0],wig[chrom][1]):
                coordinate=(int(round(1.0*wc/resolution)))*resolution+1
                
                if len(samp[0])==0:
                    samp[0].append(coordinate)
                    samp[1].append(val)
                    continue
                if coordinate!=samp[0][-1]:
                    samp[0].append(coordinate)
                    samp[1].append(val)
            
            #added to sampWig only if there are some point(s)
            if samp[0]: sampWig.wig[standardchrom]=samp
                
        return sampWig
    
class WigSamplerFast(Sampler):
    """Class WigSamplerFast
    
    This class samples a wig file (object) more fast. However, this sampler cannot guarantee very exact sampling.
    """

    def __init__(self,name=''):
        """Constructor"""
        
        Sampler.__init__(self,name='')
        
    def sample(self, wig, resolution):
        """Sample a wig file at the given resolution. 
        
        Parameters:
        1. wig: a wig object (see inout.py)
        2. resolution: sampling resolution
        
        Return:
        sampWig: the sampled wig object
        
        """
        
        # parameter checking
            
        try:
            chroms=wig.get_chroms()
        except AttributeError:
            raise Exception("Argument 'wig' must be given")
        
        sampWig=Wig()
        for chrom in chroms:
            try:
                standardchrom=self.standard_chroms[chrom]
            except KeyError:
                standardchrom=chrom
            
            try:
                start  = wig[chrom][0][0]
                end = wig[chrom][0][-1]
            except IndexError:
                continue
            
            samp=[array('l',[]),array('d',[])]
            cor = wig[chrom][0]
            val = wig[chrom][1]
            init = 0
            prev = -1000
            for sc in xrange(start, end, resolution):
                # get the closest one to the sampled point and save
                gotya = bisect_left(cor[init:], sc)
                
                if prev == (init+gotya): continue
                else: prev = (init+gotya)
                
                try:
                    samp[0].append(cor[init+gotya])
                    samp[1].append(val[init+gotya])
                    init += gotya
                except IndexError:
                    continue
        
            #added to sampWig only if there are some point(s)
            if samp[0]: sampWig.wig[standardchrom]=samp
                
        return sampWig
  
        
def fillupwig(wig,resolution,fillupval=0):
    """Fill up wig with a given value at the resolution
    
    Parameters:
    1. wig: a Wig object (see inout.py). 
    This Wig object is expected to be already sampled by WigSampler at a regular interval.
    2. resolution: sampling resolution
    3. fillupval: fill-up value, by default 0
    
    Return:
    1. fillupWig: filled-up Wig object
    
    """
    
    fillupWig=Wig()
    for chrom in wig.get_chroms():
        # check the chromosome is empty or only single element
        if len(wig[chrom][0])== 0 or len(wig[chrom][0])== 1:
            fillupWig[chrom]=wig[chrom][:]
            break
        
        cs=wig[chrom][0]
        vs=wig[chrom][1]
        ncs=[cs[0]]
        nvs=[vs[0]]
        pc=cs[0]
        pv=vs[0]
        for c,v in itertools.izip(cs[1:],vs[1:]): 
            if c-pc > resolution:
                a=range(pc+resolution,c,resolution)
                ncs+=a
                nvs+=[0.0]*len(a)
            else:
                ncs.append(c)
                nvs.append(v)    
            pc=c
            pv=v
    
        fillupWig.wig[chrom]=[array('l',ncs), array('d',nvs)]
    
    return fillupWig


    
        
        
            
    
            
        