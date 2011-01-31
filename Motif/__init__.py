# Time-stamp: <2008-03-19 15:01:53 Tao Liu>

"""Module to read the motif scan data file which is in binary format.

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

import re
import Cistrome
from Cistrome.TabIO import FWTrackI,RangeI

import sys
from struct import unpack as upk

# ------------------------------------
# constants
# ------------------------------------
__version__ = "MR $Revision$"
__author__ = "Tao Liu <taoliu@jimmy.harvard.edu>"
__doc__ = "Calculate Relationship of Motifs"
LOG = False

GENOME_SIZE = {"mm8":2644077689L,
               "hg18":3080419480L}

# ------------------------------------
# Misc functions
# ------------------------------------

def mlen (mfhd) :
    """Return the motif length from motif matrix data file.
    
    mfhd     : the file object for motif matrix file
    """
    mfhd.seek(0)
    return len(mfhd.readlines())-1

def mconsensus (mfhd):
    """Return the motif consensus for a motif matrix data file.

    mfhd    : the file object for motif matrix file
    """
    mfhd.seek(0)
    consensus_seq=""
    headline = mfhd.readline().rstrip()
    consensus_field_num = headline.split("\t").index("Consensus")
    for l in mfhd.readlines():
        l = l.rstrip()
        consensus_seq+=l.split("\t")[consensus_field_num]
    return consensus_seq

def read_motif_total_num (motif_fhd,species):
    """Only read the header of binary file, return the total number of
    motif scan hits regardless of cutoff.

    """
    if species == "hg18":
        chromosomes_fp = {                     # store start and number of file-pos for every chromosome in bin file
            "chr1":[0,0],"chr2":[0,0],"chr3":[0,0],
            "chr4":[0,0],"chr5":[0,0],"chr6":[0,0],
            "chr7":[0,0],"chr8":[0,0],"chr9":[0,0],
            "chr10":[0,0],"chr11":[0,0],"chr12":[0,0],
            "chr13":[0,0],"chr14":[0,0],"chr15":[0,0],
            "chr16":[0,0],"chr17":[0,0],"chr18":[0,0],
            "chr19":[0,0],"chr20":[0,0],"chr21":[0,0],
            "chr22":[0,0],"chrX":[0,0],"chrY":[0,0]
            }
    elif species == "mm8":
        chromosomes_fp = {                     # store start and number of file-pos for every chromosome in bin file
            "chr1":[0,0],"chr2":[0,0],"chr3":[0,0],
            "chr4":[0,0],"chr5":[0,0],"chr6":[0,0],
            "chr7":[0,0],"chr8":[0,0],"chr9":[0,0],
            "chr10":[0,0],"chr11":[0,0],"chr12":[0,0],
            "chr13":[0,0],"chr14":[0,0],"chr15":[0,0],
            "chr16":[0,0],"chr17":[0,0],"chr18":[0,0],
            "chr19":[0,0],"chrX":[0,0],"chrY":[0,0]
            }
    else:
        raise Exception("Only hg18/mm8 supported!")
        
    chromosomes = chromosomes_fp.keys()
    motif_fhd.seek(0)
    # unpack the start pos
    for chromosome in chromosomes:
        chromosomes_fp[chromosome][0] = upk("<i",motif_fhd.read(4))[0]
        motif_fhd.seek(124,1)
    motif_fhd.seek(0,2)
    
    # calculate number of hits
    total_motif_hits = 0
    for i in range(len(chromosomes)-1):
        mh = (chromosomes_fp[chromosomes[i+1]][0]-chromosomes_fp[chromosomes[i]][0])/8
        chromosomes_fp[chromosomes[i]][1] = mh
        total_motif_hits += mh
    # last one
    mh = (motif_fhd.tell()-chromosomes_fp[chromosomes[-1]][0])/8
    chromosomes_fp[chromosomes[-1]][1]=mh
    total_motif_hits += mh
 
    return total_motif_hits


def read_motif (motif_fhd,species,cutoff=0):
    """Read motif scan result, and return a TabIO.FWTrackI object
    containing the motif locations.

    motif_fhd : a file handler for binary motif scan result
    species   : must be "mm8" for mouse or "hg18" for human
    cutoff    : cutoff for the motif scan score
    """
    motif_range_list = FWTrackI(fw=0)
    if species == "hg18":
        chromosomes_fp = {                     # store start and number of file-pos for every chromosome in bin file
            "chr1":[0,0],"chr2":[0,0],"chr3":[0,0],
            "chr4":[0,0],"chr5":[0,0],"chr6":[0,0],
            "chr7":[0,0],"chr8":[0,0],"chr9":[0,0],
            "chr10":[0,0],"chr11":[0,0],"chr12":[0,0],
            "chr13":[0,0],"chr14":[0,0],"chr15":[0,0],
            "chr16":[0,0],"chr17":[0,0],"chr18":[0,0],
            "chr19":[0,0],"chr20":[0,0],"chr21":[0,0],
            "chr22":[0,0],"chrX":[0,0],"chrY":[0,0]
            }
    elif species == "mm8":
        chromosomes_fp = {                     # store start and number of file-pos for every chromosome in bin file
            "chr1":[0,0],"chr2":[0,0],"chr3":[0,0],
            "chr4":[0,0],"chr5":[0,0],"chr6":[0,0],
            "chr7":[0,0],"chr8":[0,0],"chr9":[0,0],
            "chr10":[0,0],"chr11":[0,0],"chr12":[0,0],
            "chr13":[0,0],"chr14":[0,0],"chr15":[0,0],
            "chr16":[0,0],"chr17":[0,0],"chr18":[0,0],
            "chr19":[0,0],"chrX":[0,0],"chrY":[0,0]
            }
    else:
        raise Exception("Only hg18/mm8 supported!")
        
    chromosomes = chromosomes_fp.keys()
    motif_fhd.seek(0)
    # unpack the start pos
    for chromosome in chromosomes:
        chromosomes_fp[chromosome][0] = upk("<i",motif_fhd.read(4))[0]
        motif_fhd.seek(124,1)
    motif_fhd.seek(0,2)
    
    # calculate number of hits
    total_motif_hits = 0
    for i in range(len(chromosomes)-1):
        mh = (chromosomes_fp[chromosomes[i+1]][0]-chromosomes_fp[chromosomes[i]][0])/8
        chromosomes_fp[chromosomes[i]][1] = mh
        total_motif_hits += mh
    # last one
    mh = (motif_fhd.tell()-chromosomes_fp[chromosomes[-1]][0])/8
    chromosomes_fp[chromosomes[-1]][1]=mh
    total_motif_hits += mh

    # read and write
    read_motif_hits = 0
    portion = 0
    for chromosome in chromosomes:
        motif_fhd.seek(chromosomes_fp[chromosome][0],0)
        for i in range(chromosomes_fp[chromosome][1]):
            read_motif_hits += 1
            portion = float(read_motif_hits)/total_motif_hits
            if LOG:
                sys.stdout.write("\r%.1f%% %s" % (portion*100,"#"*int(portion*50)))
                sys.stdout.flush()
            loc = upk("<i",motif_fhd.read(4))[0]
            score = upk("<f",motif_fhd.read(4))[0]
            motif_fhd.read(4)
            if score < 0:
               strand = -1
               score = score*-1
            else:
               strand = 1
            #ofhd.write("%s\t%d\t%d\t%s_%s_%d\t%.2f\t%s\n" % (chromosome,loc-1,loc+motif_len-1,motif,chromosome,i,score,strand))
            if score > cutoff:
		#print score,cutoff
                motif_range_list.add_range(chromosome,RangeI(start=loc-1,end=loc,strand=strand))
            #print loc-1
    #sys.stdout.write("\n")
    motif_range_list.merge_overlap()
    return motif_range_list

def read_motif2 (motif_fhd,species,cutoff=0):
    """Read motif scan result, and return a TabIO.FWTrackI object
    containing the motif locations.

    * If the motif scan data file is not big, use this function to
      load the whole file into memory. It may be faster than
      read_motif().

    motif_fhd : a file handler for binary motif scan result
    species   : must be "mm8" for mouse or "hg18" for human
    cutoff    : cutoff for the motif scan score
    """
    motif_range_list = FWTrackI(fw=0)
    if species == "hg18":
        chromosomes_fp = {                     # store start and number of file-pos for every chromosome in bin file
            "chr1":[0,0],"chr2":[0,0],"chr3":[0,0],
            "chr4":[0,0],"chr5":[0,0],"chr6":[0,0],
            "chr7":[0,0],"chr8":[0,0],"chr9":[0,0],
            "chr10":[0,0],"chr11":[0,0],"chr12":[0,0],
            "chr13":[0,0],"chr14":[0,0],"chr15":[0,0],
            "chr16":[0,0],"chr17":[0,0],"chr18":[0,0],
            "chr19":[0,0],"chr20":[0,0],"chr21":[0,0],
            "chr22":[0,0],"chrX":[0,0],"chrY":[0,0]
            }
        chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6",
                       "chr7","chr8","chr9","chr10","chr11","chr12",
                       "chr13","chr14","chr15","chr16","chr17","chr18",
                       "chr19","chr20","chr21","chr22","chrX","chrY"]
    elif species == "mm8":
        chromosomes_fp = {                     # store start and number of file-pos for every chromosome in bin file
            "chr1":[0,0],"chr2":[0,0],"chr3":[0,0],
            "chr4":[0,0],"chr5":[0,0],"chr6":[0,0],
            "chr7":[0,0],"chr8":[0,0],"chr9":[0,0],
            "chr10":[0,0],"chr11":[0,0],"chr12":[0,0],
            "chr13":[0,0],"chr14":[0,0],"chr15":[0,0],
            "chr16":[0,0],"chr17":[0,0],"chr18":[0,0],
            "chr19":[0,0],"chrX":[0,0],"chrY":[0,0]
            }
        chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6",
                       "chr7","chr8","chr9","chr10","chr11","chr12",
                       "chr13","chr14","chr15","chr16","chr17","chr18",
                       "chr19","chrX","chrY"]
    else:
        raise Exception("Only hg18/mm8 supported!")
        
    motif_fhd.seek(0)
    data = motif_fhd.read()
    # unpack the start pos
    p = 0
    for chromosome in chromosomes:
        chromosomes_fp[chromosome][0] = upk("<i",data[p:p+4])[0]
        p += 128

    # calculate number of hits
    total_motif_hits = 0
    for i in range(len(chromosomes)-1):
        mh = (chromosomes_fp[chromosomes[i+1]][0]-chromosomes_fp[chromosomes[i]][0])/8
        chromosomes_fp[chromosomes[i]][1] = mh
        total_motif_hits += mh
    # last one
    mh = (len(data)-chromosomes_fp[chromosomes[-1]][0])/8
    chromosomes_fp[chromosomes[-1]][1]=mh
    total_motif_hits += mh

    # read and write
    read_motif_hits = 0
    portion = 0
    p = 0

    n=0
    for chromosome in chromosomes:
        p = chromosomes_fp[chromosome][0]
        for i in range(chromosomes_fp[chromosome][1]):
            read_motif_hits += 1
            portion = float(read_motif_hits)/total_motif_hits
            if LOG:
                sys.stdout.write("\r  %.1f%% %s" % (portion*100,"#"*int(portion*50)))
                sys.stdout.flush()
            loc = upk("<i",data[p:p+4])[0]
            score = upk("<f",data[p+4:p+8])[0]
            p += 8
            if score < 0:
               strand = -1
               score = score*-1
            else:
               strand = 1
            #ofhd.write("%s\t%d\t%d\t%s_%s_%d\t%.2f\t%s\n" % (chromosome,loc-1,loc+motif_len-1,motif,chromosome,i,score,strand))
            if score > cutoff:
		#print score,cutoff
                n+=1
                motif_range_list.add_range(chromosome,RangeI(start=loc-1,end=loc,strand=strand))
            #print loc-1
    if LOG : sys.stdout.write("\n")
    data = None
    motif_range_list.merge_overlap()
    #print n
    return motif_range_list

def motif_count (track, motif_track):
    """Count how many motif discovered in a given track.
    
    """
    return track.include(motif_track)

