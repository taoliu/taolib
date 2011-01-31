# Time-stamp: <2008-07-29 18:20:54 zhenhua>

"""Module Description

Copyright (c) 2008 Zhenhua Wu <jeremywu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Zhenhua Wu
@contact: jeremywu@jimmy.harvard.edu
"""


# ------------------------------------
# python modules
# ------------------------------------

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Class
# ------------------------------------


ambiguous_nucleotide_code_codec = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "N": "ACGT",
    "AC": "M",
    "AG": "R",
    "AT": "W",
    "CG": "S",
    "CT": "Y",
    "GT": "K",
    "ACG": "V",
    "ACT": "H",
    "AGT": "D",
    "CGT": "B",
    "ACGT": "N"
    }


class pwMatrix(object):
    nucleotide_column = {"A":0, "C":1, "G":2, "T":3, "a":0, "c":1, "g":2, "t":3, 0:"A", 1:"C", 2:"G", 3:"T"}

    def __init__(self):
        self.pwm = []


    def __str__(self):
        lines = []

        lines.append("\tA\tC\tG\tT\tCon.\tDeg.")

        for index in range(len(self)):
            lines.append("%d\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s" % \
                         (index + 1,\
                          self.getProbabilityWeight(index, "A"),\
                          self.getProbabilityWeight(index, "C"),\
                          self.getProbabilityWeight(index, "G"),\
                          self.getProbabilityWeight(index, "T"),\
                          self.consensus(index), \
                          self.degenerate(index) \
                          )\
                         )

        return "\n".join(lines)

    def toDBstring(self):
        lines = []
#        lines.append("\tA\tC\tG\tT\tCon.\tDeg.")
        for index in range(len(self)):
            lines.append("%.4f\t%.4f\t%.4f\t%.4f" % \
                         (self.getProbabilityWeight(index, "A"),\
                          self.getProbabilityWeight(index, "C"),\
                          self.getProbabilityWeight(index, "G"),\
                          self.getProbabilityWeight(index, "T"),\
                          )\
                         )

        return "\n".join(lines)

    def append(self, nucleotide_weight_dict):
        if len(nucleotide_weight_dict) != 4:
            raise Exception, "The parameter should have length of 4!"
        
        row = [0, 0, 0, 0]
        
        for nucleotide in nucleotide_weight_dict.keys():
            row[self.__getColumn__(nucleotide)] = nucleotide_weight_dict[nucleotide]
            
        self.pwm.append(row)

    def __getColumn__(self, nucleotide):
        return pwMatrix.nucleotide_column[nucleotide]

    def getProbabilityWeight(self, position, nucleotide):
        return self.pwm[position][ self.__getColumn__(nucleotide) ]

    def __len__(self):
        return len(self.pwm)

    def clear(self):
        self.pwm = []

    def consensus(self, position):
        
        weightRow = self.pwm[position]
        
        maxP = max(weightRow)
        baseWithMaxP = ""
        
        for index in range(4):
            if weightRow[index] == maxP:
                baseWithMaxP += pwMatrix.nucleotide_column[index]

        return ambiguous_nucleotide_code_codec[baseWithMaxP]

    def degenerate(self, position):
        weightRow = self.pwm[position]
        
        degenerateBases = ""
        
        for index in range(4):
            if weightRow[index] >= 0.25:
                degenerateBases += pwMatrix.nucleotide_column[index]

        return ambiguous_nucleotide_code_codec[degenerateBases]

    def getICat(self, position):
        ic = 2.0
        letters = "ACGT"

        for letter in letters:
            p = self.getProbabilityWeight(position, letter)
            if p != 0:
                ic += p * math.log(p, 2)

        return ic

    def getIC(self):
        total_ic = 0
        for position in range(len(self)):
            total_ic += self.getICat(position)

        return total_ic
