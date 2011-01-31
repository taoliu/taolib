# Time-stamp: <2009-05-08 11:38:35 Tao Liu>

"""Module Description: IO module to read chip-seq tags from BED or ELAND

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

from taolib.CoreLib.FeatIO import FWTrackI
import logging

# ------------------------------------
# Classes
# ------------------------------------

class BEDParser:
    """File Parser Class for tabular File.

    """
    def __init__ (self):
        """All parameters except for format are the same as python's
        builtin file class.

        Format parameter can be "bed","gff" or other particular string
        for other format of file. For example, if your file is '\\t'
        delimited, and the first column of the file is chromosome
        name, the second column is start position, third column is end
        position, fourth column is the strand, the coordinates are
        0-indexed and the range is closed, then you should write the
        format string as "123401\\t" (six numbers and the delimiter).

        Note: Use the fifth and sixth number in format string to
        indicate whether the coordinates are 0-indexed (0) or
        1-indexed (1) and whether the range is opened (0) or closed
        (1) (i.e. whether the end position is included in the range
        or not)
        """
        self.__format__ = "123600\t"
        self.__delimiter__ = self.__format__[6:]
        self.__chr_col__ = int(self.__format__[0])-1
        self.__start_col__ = int(self.__format__[1])-1
        self.__end_col__ = int(self.__format__[2])-1
        self.__strand_col__ = int(self.__format__[3])-1
        self.__start_shift__ = int(self.__format__[4])
        self.__end_shift__ = int(self.__format__[5])

    def build_track (self,fhd):
        pass

    def build_fwtrack (self, fhd):
        """Build FWTrackI from all lines, return a FWTrackI object.

        Note: All ranges will be merged (exclude the same
        range) then sorted after the track is built.

        If both_strand is True, it will store strand information in
        FWTrackI object.

        if do_merge is False, it will not merge the same range after
        the track is built.
        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        for thisline in fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if fpos==None or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)
        if thisline.startswith("#"): return ("comment line",None,None) # comment line is skipped
        thisfields = thisline.split()
        if len(thisfields) < 6 : # default pos strand if no strand
                                 # info can be found
            return (thisfields[0],
                    int(thisfields[1]),
                    0)
        else:
            if thisfields[5] == "+":
                return (thisfields[0],
                        int(thisfields[1]),
                        0)
            elif thisfields[5] == "-":
                return (thisfields[0],
                        int(thisfields[2]),
                        1)
            else:
                raise self.StrandFormatError(thisline,strand)

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)


class ELANDResultParser:
    """File Parser Class for ELAND result File.

    Note this parser can only work for s_N_eland_results.txt format.

    Unfiltered ELAND alignment output 
    Each line of the output file contains the following fields: 
    1. Sequence name (derived from file name and line number if format is not fasta) 
    2. Sequence 
    3. Type of match codes: 
    •NM—No match found 
    •QC—No matching done: QC failure (too many Ns) 
    •RM—No matching done: repeat masked (may be seen if repeatFile.txt was 
    specified) 
    •U0—Best match found was a unique exact match 
    •U1—Best match found was a unique 1-error match 
    •U2—Best match found was a unique 2-error match 
    •R0—Multiple exact matches found 
    •R1—Multiple 1-error matches found, no exact matches 
    •R2—Multiple 2-error matches found, no exact or 1-error matches 
    4. Number of exact matches found 
    5. Number of 1-error matches found 
    6. Number of 2-error matches found 
    7. The following fields are only used if a unique best match was found: 
    8. Genome file in which match was found 
    9. Position of match (bases in file are numbered starting at 1) 
    10.Direction of match (F=forward strand, R=reverse) 
    11.How N characters in read were interpreted (“.”=not applicable, “D”=Detection, 
    “I”=Insertion) 
    The following field is only used in the case of a unique inexact match: 
    12.Position and type of first substitution error (A numeral refers to a run of matching 
    bases, an upper case base or N refers to a base in the reference that differs from the 
    read. For example, 11A: after 11 matching bases, base 12 is A in the reference but 
    not in the read) 

    """
    def __init__ (self):
        """
        """
        pass

    def build_fwtrack (self, fhd):
        """Build FWTrackI from all lines, return a FWTrackI object.

        Note only unique matches are kept.
        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        for thisline in fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        if not thisline: return (None,None,None)
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)

        if thisline[0] == "#": return ("comment line",None,None) # comment line is skipped
        thisfields = thisline.split("")
        thistaglength = len(thisfields[1])

        if thisfields[2] == "U0" or thisfields[2]=="U1" or thisfields[2]=="U2":
            strand = thisfields[8]
            if strand == "F":
                return (thisfields[6],
                        int(thisfields[7])-1,
                        0)
            elif strand == "R":
                return (thisfields[6].rstrip(".fa"),
                        int(thisfields[7])+thistaglength-1,
                        1)
            else:
                raise self.StrandFormatError(thisline,strand)
        else:
            return (None,None,None)

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)


class ELANDMultiParser:
    """File Parser Class for ELAND multi File.

    Note this parser can only work for s_N_eland_multi.txt format.

    Each line of the output file contains the following fields: 
    1. Sequence name 
    2. Sequence 
    3. Either NM, QC, RM (as described above) or the following: 
    4. x:y:z where x, y, and z are the number of exact, single-error, and 2-error matches 
    found 
    5. Blank, if no matches found or if too many matches found, or the following: 
    BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1 
    This says there are two matches to BAC_plus_vector.fa: one in the reverse direction 
    starting at position 160322 with one error, one in the forward direction starting at 
    position 170128 with two errors. There is also a single-error match to E_coli.fa.
    """
    def __init__ (self):
        """
        """
        pass

    def build_fwtrack (self, fhd):
        """Build FWTrackI from all lines, return a FWTrackI object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        for thisline in fhd:
            (chromosome,fpos,strand) = self.__fw_parse_line(thisline)
            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
            if not fpos or not chromosome:
                continue
            fwtrack.add_loc(chromosome,fpos,strand)
        return fwtrack
    
    def __fw_parse_line (self, thisline ):
        if not thisline: return (None,None,None)
        thisline = thisline.rstrip()
        if not thisline: return ("blank",None,None)

        if thisline[0] == "#": return ("comment line",None,None) # comment line is skipped
        thisfields = thisline.split()
        thistagname = thisfields[0]         # name of tag
        thistaglength = len(thisfields[1]) # length of tag

        if len(thisfields) < 4:
            return (None,None,None)
        else:
            thistaghits = sum(map(int,thisfields[2].split(':')))
            if thistaghits > 1:
                # multiple hits
                return (None,None,None)
            else:
                (name,pos) = thisfields[3].split(':')
                strand  = pos[-2]
                if strand == "F":
                    return (name.rstrip('.fa'),
                            int(pos[:-2])-1,
                            0)
                elif strand == "R":
                    return (name.rstrip('.fa'),
                            int(pos[:-2])+thistaglength-1,
                            1)
                else:
                    raise self.StrandFormatError(thisline,strand)

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)

class PairEndELANDMultiParser:
    """File Parser Class for two ELAND multi Files for Pair-End sequencing.

    Note this parser can only work for s_N_eland_multi.txt format.

    Each line of the output file contains the following fields: 
    1. Sequence name 
    2. Sequence 
    3. Either NM, QC, RM (as described above) or the following: 
    4. x:y:z where x, y, and z are the number of exact, single-error, and 2-error matches 
    found 
    5. Blank, if no matches found or if too many matches found, or the following: 
    BAC_plus_vector.fa:163022R1,170128F2,E_coli.fa:3909847R1 
    This says there are two matches to BAC_plus_vector.fa: one in the reverse direction 
    starting at position 160322 with one error, one in the forward direction starting at 
    position 170128 with two errors. There is also a single-error match to E_coli.fa.
    """
    def __init__ (self):
        """
        """
        pass

    def build_fwtrack (self, lfhd, rfhd, dist=200):
        """Build FWTrackI from all lines, return a FWTrackI object.

        lfhd: the filehandler for left tag file
        rfhd: the filehandler for right tag file
        dist: the best distance between two tags in a pair

        The score system for pairing two tags:

        score = abs(abs(rtag-ltag)-200)+error4lefttag+error4righttag

        the smaller score the better pairing. If the score for a
        pairing is bigger than 200, this pairing will be discarded.

        Note only the best pair is kept. If there are over two best
        pairings, this pair of left and right tags will be discarded.

        Note, the orders in left tag file and right tag file must
        match, i.e., the Nth left tag must has the same name as the
        Nth right tag.

        Note, remove comment lines beforehand.
        """
        fwtrack = FWTrackI()
        i = 0
        m = 0
        lnext = lfhd.next
        rnext = rfhd.next
        self.dist = dist
        try:
            while 1:
                lline = lnext()
                rline = rnext()
                (chromosome,fpos,strand) = self.__fw_parse_line(lline,rline)
                i+=1
                if i == 1000000:
                    m += 1
                    logging.info(" %d" % (m*1000000))
                    i=0
                if not fpos or not chromosome:
                    continue
                fwtrack.add_loc(chromosome,fpos,strand)
        except StopIteration:
            pass
        return fwtrack
    
    def __fw_parse_line (self, leftline, rightline ):
        # >HWI-EAS275_5:4:100:340:1199/1	GTGCTGGTGGAGAGGGCAAACCACATTGACATGCT	2:1:0	chrI.fa:15061365F0,15068562F0,chrIV.fa:4783988R1
        # >HWI-EAS275_5:4:100:340:1199/2	GGTGGTGTGTCCCCCTCTCCACCAGCACTGCGGCT	3:0:0	chrI.fa:15061451R0,15068648R0,15071742R0

        leftfields = leftline.split()
        lefttaglength = len(leftfields[1]) # length of tag
        rightfields = rightline.split()
        righttaglength = len(rightfields[1]) # length of tag

        if len(rightfields) < 4 or len(leftfields) < 4:
            # one of the tag cann't be mapped to genome
            return (None,None,None)
        else:
            lefthits = self.__parse_line_to_dict(leftfields[3])
            righthits = self.__parse_line_to_dict(rightfields[3])            
            parings = []

            for seqname in lefthits.keys():
                if not righthits.has_key(seqname):
                    continue
                else:
                    leftpses = lefthits[seqname] # pse=position+strand+error
                    rightpses = righthits[seqname]
                    for (lp,ls,le) in leftpses:
                        for (rp,rs,re) in rightpses:
                            # try to pair them
                            if ls == 'F':
                                if rs == 'R':
                                    score = abs(abs(rp-lp)-self.dist)+le+re
                                    if score < 200:
                                        parings.append((score,seqname,int((lp+rp)/2),0) )
                                else:
                                    # strands don't match
                                    continue
                            else:
                                if rs == 'F':
                                    score = abs(abs(rp-lp)-self.dist)+le+re
                                    if score < 200:
                                        parings.append((score,seqname,int((lp+rp)/2),1) )
                                else:
                                    # strands don't match
                                    continue
            if not parings:
                return (None,None,None)
            parings.sort()
            if len(parings)>1 and parings[0][0] == parings[1][0]:
                # >2 best paring, reject!
                return (None,None,None)
            else:
                return parings[0][1:]                                
                    
    def __parse_line_to_dict ( self, linestr ):
        items = linestr.split(',')
        hits = {}
        for item in items:
            if item.find(':') != -1:
                # a seqname section
                (n,pse) = item.split(":") # pse=position+strand+error
                hits[n]=[]
                try:
                    sindex = pse.rindex('F')
                except ValueError:
                    sindex = pse.rindex('R')
                p = int(pse[:sindex])
                s = pse[sindex]
                e = int(pse[sindex+1:])
                hits[n].append((p,s,e))
            else:
                # only pse section
                try:
                    sindex = pse.rindex('F')
                except ValueError:
                    sindex = pse.rindex('R')
                p = int(pse[:sindex])
                s = pse[sindex]
                e = int(pse[sindex+1:])
                hits[n].append((p,s,e))
        return hits

    class StrandFormatError(Exception):
        def __init__ (self, string, strand):
            self.message = "Strand information can not be recognized in this line: \"%s\",\"%s\"" % (string,strand)

        def __str__ (self):
            return repr(self.message)
