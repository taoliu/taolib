# Time-stamp: <2008-07-29 18:19:22 zhenhua>

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

"""
Each motif should be stored as:

Position  A C G T
1        Probability matrix
2
3

Source, literature, concensus (the highest number of each row) (sequence logo), if total information containted is less than 0.5, don't consider)

source: Transfac, JA, TR, BU, Li, LiuLab
gene name, (gene that codes the transcription factor)
symbol, (like a ID,)

domain, (which domain of the transcription factor binds this motif)
structure categrory, (structure of the protein)
organism (search online)
"""

import math
import random
import os


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
            

class transfacRecord(object):

    def __init__(self, factorName = None, positionWeightMatrix = None, species = None, accessionNumber = None, PMID = None):
        self.pwmatrix = positionWeightMatrix
        self.factorName = factorName
        self.species = species
        self.accessionNumber = accessionNumber
        self.PMID = PMID

    def save(self):
        if self.accessionNumber == None:
            raise Exception, "invalid accessionNumber"

        print self.accessionNumber
        
        outputName = self.accessionNumber + ".motif"
        
        outputfile = open(outputName, 'w')
        
        outputfile.write("#Position Probability Weight Matrix:\n")
        outputfile.write("%s\n" % str(self.pwmatrix))
        
        outputfile.write("#Source:\nTransfac\n")
        
        outputfile.write("#Factor Name:\n%s\n" % ", ".join(self.factorName))

        if self.species == None:
            outputfile.write("#Species:\nN/A\n")
        else:
            outputfile.write("#Species:\n%s\n" % ", ".join(self.species) )

        if self.PMID == None:
            outputfile.write("#PMID:\nN/A\n")
        else:
            outputfile.write("#PMID:\n%s\n" % self.PMID )

        outputfile.write("#Domain:\nN/A\n")
        outputfile.write("#Structure Category:\nN/A\n")

        outputfile.close()


def transfacRecordParser(lines):

    if len(lines) == 0:
        return None

    record = transfacRecord()
    
    for line in lines:

        if line.startswith("AC"):
            record.accessionNumber = line.split()[1]
            
        elif line.startswith("NA"):
            line = line.strip()
            line = line[4:]

            if not "C/EB" in line:
                line = line.replace("/", ",")

            if not "Su(H)" in line:
                line = line.replace("(", ",")
                line = line.replace(")", ",")
            
            line = line.replace(":", ",")
            line = line.replace("+", ",")
            
            items =  line.split(",")
            record.factorName = []
            for item in items:
                item = item.strip()
                if item != "":
                    record.factorName.append(item)
                    
        elif line.startswith("BF"):
            triplets = line.strip().split(";")

            if len(triplets) != 3:
                print triplets
                raise Exception, "factors format exception: %s" % line
                

            if record.species == None:
                record.species = []

            speciesTemp = triplets[2].replace(".", " ").split(":")[1].strip().split(",")

            for temp in speciesTemp:
                temp = temp.strip()
                if not temp in record.species:
                    record.species.append(temp)

        elif line.startswith("P0"):
            items = line.split()
            nucleotideList = items[1:5]
            record.pwmatrix = pwMatrix()
            
            
        elif line[0:2].isdigit():
            items = line.split()
            lineIndex = int(items[0])

            if lineIndex != len(record.pwmatrix) + 1:
                raise Exception, "wrong matrix"

            totalCount = float(items[1]) + float(items[2]) + float(items[3]) + float(items[4])
                    
            record.pwmatrix.append( { nucleotideList[0]:float(items[1]) / totalCount, \
                                      nucleotideList[1]:float(items[2]) / totalCount, \
                                      nucleotideList[2]:float(items[3]) / totalCount, \
                                      nucleotideList[3]:float(items[4]) / totalCount} )
        elif line.startswith("RX"):
            line = line.replace(".", " ")
            record.PMID = line.split()[2].strip()


    return record

def transfacMatrixTableParser(transfacMT):

    lines = []
    while True:

        line = transfacMT.readline()

        if line == '':
            break;

        if line.startswith("AC"):

            if len(lines) != 0:
                record = transfacRecordParser(lines)
                record.save()
                lines = []
                
            lines.append(line)
            
        elif len(lines) != 0:
            lines.append(line)



def transfacRawFileParser(rawFile):

    species = []
    factorName = []

    while True:

        line = rawFile.readline()

        if line.startswith("Accession Number"):
            line = rawFile.readline()
            accessionNumber = line.strip()
            
        if line.startswith("Factor name"):
            line = rawFile.readline()
            
            line = line.strip()

            line = line.replace("/", ",")
            line = line.replace(":", ",")
            line = line.replace("(", ",")
            line = line.replace(")", ",")
            line = line.replace("+", ",")
#            line = line.replace(" ", ",")

            
            items =  line.split(",")
            for item in items:
                item = item.strip()
                item = item.replace("&", "/")
                if item != "":
                    factorName.append(item)

        if line.startswith("Factors"):
            line = rawFile.readline().strip()

            factors = line.split(". T")
#            print factors

            for factorItem in factors:

                if factorItem == '':
                    continue
                
                triplets = factorItem.strip().split(";")

                if len(triplets) != 3:
                    print triplets
                    raise Exception, "factors format exception: %s" % line

                if triplets[1].strip() in factorName :
                    speciesTemp = triplets[2].split(":")[1].strip().split(",")

                    for temp in speciesTemp:
                        temp = temp.strip()
                        if not temp in species:
                            species.append(temp)

        if line.startswith("Binding Matrix"):
            
            lineIndex = 0

            pwm = pwMatrix()
    
            while True:
                line = rawFile.readline()

                items = line.split()

                if len(items) != 5:
                    break

                if lineIndex == 0:
                    nucleotideList = items[0:4]
                else:
                    totalCount = float(items[0]) + float(items[1]) + float(items[2]) + float(items[3])
                    
                    pwm.append( { nucleotideList[0]:float(items[0]) / totalCount, \
                                  nucleotideList[1]:float(items[1]) / totalCount, \
                                  nucleotideList[2]:float(items[2]) / totalCount, \
                                  nucleotideList[3]:float(items[3]) / totalCount} )
                    
                lineIndex = lineIndex + 1

        if not line:
            break

    return transfacRecord(factorName, pwm, species, accessionNumber)



def transfac_matrix_convert(transfacMatrixFile):

    lineIndex = 0

    pwm = pwMatrix()
    
    for line in transfacMatrixFile:
        items = line.split()

        if lineIndex == 0:
            assert(len(items) == 5)
            nucleotideList = items[0:4]
        else:
            pwm.append( { nucleotideList[0]:float(items[1]) / 100.0, \
                          nucleotideList[1]:float(items[2]) / 100.0, \
                          nucleotideList[2]:float(items[3]) / 100.0, \
                          nucleotideList[3]:float(items[4]) / 100.0} )

        lineIndex = lineIndex + 1

    return pwm


def transfac_rawdata_convert(transfactRawFileName):


    for rawFileName in transfactRawFileName:
        print rawFileName.strip()
        rawFile = open(rawFileName.strip(), 'r')
        
        
        transfacRecord = transfacRawFileParser(rawFile)
        
        outputName = transfacRecord.accessionNumber + "." + ":".join([x.replace("/", "_slash_") for x in transfacRecord.factorName]) + ".motif"
        
        outputfile = open(outputName, 'w')
        
        outputfile.write("#Position Probability Weight Matrix:\n")
        outputfile.write("%s\n" % str(transfacRecord.pwmatrix))
        
        outputfile.write("#Source:\nTransfac\n")
        
        outputfile.write("#Factor Name:\n%s\n" % ", ".join(transfacRecord.factorName))

        if len(transfacRecord.species) == 0:
            outputfile.write("#Species:\nN/A\n")
        else:
            outputfile.write("#Species:\n%s\n" % ", ".join(transfacRecord.species) )

        outputfile.write("#Domain:\nN/A\n")
        outputfile.write("#Structure Category:\nN/A\n")

        rawFile.close()
        outputfile.close()


####################################################################
#
# Jaspar motif parser
#
#
####################################################################


class jasparRecord(object):
    
    def __init__(self, factorName = None, positionWeightMatrix = None, species = None, PMID = None):
        self.pwmatrix = positionWeightMatrix
        self.factorName = factorName
        self.species = species
        self.PMID = PMID


def jasparRawFileParser(rawFile):

    dataMatrix = []

    while True:
        line = rawFile.readline()

        if line == '':
            break
        
        items = line.split()

        dataMatrix.append([int(item) for item in items])

    assert( len(dataMatrix) == 4)

    motifLength = len(dataMatrix[0])

    assert( motifLength == len(dataMatrix[1]) )
    assert( motifLength == len(dataMatrix[2]) )
    assert( motifLength == len(dataMatrix[3]) )

    pwm = pwMatrix()
    
    nucleotideList = ['A', 'C', 'G', 'T']
    
    for position in range(motifLength):
        
        totalCount = dataMatrix[0][position] + dataMatrix[1][position] + dataMatrix[2][position] + dataMatrix[3][position]

        pwm.append( { nucleotideList[0]:float( dataMatrix[0][position] ) / totalCount, \
                      nucleotideList[1]:float( dataMatrix[1][position] ) / totalCount, \
                      nucleotideList[2]:float( dataMatrix[2][position] ) / totalCount, \
                      nucleotideList[3]:float( dataMatrix[3][position] ) / totalCount} )
                    
    return jasparRecord(positionWeightMatrix = pwm)


def jaspar_rawdata_convert(annotationFile, jasparRawFileName):


    annotationDict = {}

    for line in annotationFile:
        items = line.split('\t')
        
        motifID = items[0]
        
        if not annotationDict.has_key(motifID):
            annotationDict[motifID] = {}

        annotationDict[motifID][items[1]] = items[2].strip()



    for rawFileName in jasparRawFileName:
        print rawFileName.strip()

        motifID = rawFileName.split('.')[0]

        rawFile = open(rawFileName.strip(), 'r')
        
        jasparRecord = jasparRawFileParser(rawFile)

        if annotationDict[motifID].has_key("name"):
            jasparRecord.factorName = annotationDict[motifID]["name"]
            
        if annotationDict[motifID].has_key("species"):
            jasparRecord.species = annotationDict[motifID]["species"]

        if annotationDict[motifID].has_key("medline"):
            jasparRecord.PMID = annotationDict[motifID]["medline"]
        
        
        outputName = motifID + ".motif"
        
        outputfile = open(outputName, 'w')
        
        outputfile.write("#Position Probability Weight Matrix:\n")
        outputfile.write("%s\n" % str(jasparRecord.pwmatrix))
        
        outputfile.write("#Source:\nJaspar\n")
        
        outputfile.write("#Factor Name:\n%s\n" % jasparRecord.factorName)

        if len(jasparRecord.species) == 0 or jasparRecord.species == None or jasparRecord.species == '-':
            outputfile.write("#Species:\nN/A\n")
        else:
            outputfile.write("#Species:\n%s\n" % jasparRecord.species )

        if len(jasparRecord.PMID) == 0 or jasparRecord.PMID == None or jasparRecord.PMID == '-':
            outputfile.write("#PMID:\nN/A\n")
        else:
            outputfile.write("#PMID:\n%s\n" % jasparRecord.PMID )

        outputfile.write("#Domain:\nN/A\n")
        outputfile.write("#Structure Category:\nN/A\n")
        
        rawFile.close()
        outputfile.close()


####################################################################
#
# DFCI  motif parser
#
#
####################################################################

class dfciRecord(object):
    
    def __init__(self, factorName = None, positionWeightMatrix = None):
        self.pwmatrix = positionWeightMatrix
        self.factorName = factorName

def dfciRawFileParser(rawFile):

    lineIndex = 0
    
    pwm = pwMatrix()
    
    while True:
        line = rawFile.readline()

        if line == '':
            break

        items = line.split()
        
        if lineIndex == 0:
            nucleotideList = items[0:4]
        else:
            totalCount = float(items[1]) + float(items[2]) + float(items[3]) + float(items[4])
            
            pwm.append( { nucleotideList[0]:float(items[1]) / totalCount, \
                          nucleotideList[1]:float(items[2]) / totalCount, \
                          nucleotideList[2]:float(items[3]) / totalCount, \
                          nucleotideList[3]:float(items[4]) / totalCount} )
                    
        lineIndex = lineIndex + 1

    return dfciRecord( positionWeightMatrix = pwm )

def dfci_rawdata_convert(dfciRawFileName):
    
    for rawFileName in dfciRawFileName:
        print rawFileName.strip()

        factorName = rawFileName.split('.')[1]

        rawFile = open(rawFileName.strip(), 'r')
        
        dfciRecord = dfciRawFileParser(rawFile)

        dfciRecord.factorName = factorName
        
        outputName = "DFCI" + "." + dfciRecord.factorName + ".motif"
        
        outputfile = open(outputName, 'w')
        
        outputfile.write("#Position Probability Weight Matrix:\n")
        outputfile.write("%s\n" % str(dfciRecord.pwmatrix))
        
        outputfile.write("#Source:\nDFCI\n")
        
        outputfile.write("#Factor Name:\n%s\n" % dfciRecord.factorName)

        outputfile.write("#Species:\nN/A\n")

        outputfile.write("#PMID:\nN/A\n")

        outputfile.write("#Domain:\nN/A\n")

        outputfile.write("#Structure Category:\nN/A\n")
        
        rawFile.close()
        outputfile.close()


####################################################################
#
# Martha Bulyk motif parser
#
#
####################################################################


class marthaBulykRecord(object):

    def __init__(self, factorName = None, positionWeightMatrix = None, species = None, accessionNumber = None, PMID = None):
        self.pwmatrix = positionWeightMatrix
        self.factorName = factorName
        self.species = species
        self.accessionNumber = accessionNumber
        self.PMID = PMID

    def save(self):
        if self.accessionNumber == None:
            raise Exception, "invalid accessionNumber"

        print self.accessionNumber
        
        outputName = self.accessionNumber + ".motif"
        
        outputfile = open(outputName, 'w')
        
        outputfile.write("#Position Probability Weight Matrix:\n")
        outputfile.write("%s\n" % str(self.pwmatrix))
        
        outputfile.write("#Source:\nMartha Bulyk\n")
        
        outputfile.write("#Factor Name:\n%s\n" % self.factorName)

        if self.species == None:
            outputfile.write("#Species:\nN/A\n")
        else:
            outputfile.write("#Species:\n%s\n" % ", ".join(self.species) )

        if self.PMID == None:
            outputfile.write("#PMID:\nN/A\n")
        else:
            outputfile.write("#PMID:\n%s\n" % self.PMID )

        outputfile.write("#Domain:\nN/A\n")
        outputfile.write("#Structure Category:\nN/A\n")

        outputfile.close()


def marthaBulykRecordParser(lines):

    if len(lines) == 0:
        return None

    record = marthaBulykRecord()

    pwmData = {}
    
    for line in lines:

        if '_' in line:
            record.factorName = line.split('_')[0].strip()
            
        elif line.startswith("A:"):
            items = line.split('\t')
            
            pwmData['A'] = [ float(item) for item in items[1:] ]
            
        elif line.startswith("C:"):
            items = line.split('\t')
            
            pwmData['C'] = [ float(item) for item in items[1:] ]
            
        elif line.startswith("G:"):
            items = line.split('\t')
            
            pwmData['G'] = [ float(item) for item in items[1:] ]
            
        elif line.startswith("T:"):
            items = line.split('\t')
            
            pwmData['T'] = [ float(item) for item in items[1:] ]


    motifLen = len(pwmData['A'])
    assert(motifLen == len(pwmData['C']))
    assert(motifLen == len(pwmData['G']))
    assert(motifLen == len(pwmData['T']))

    record.pwmatrix = pwMatrix()

    for lineIndex in range(motifLen):

        totalCount = pwmData['A'][lineIndex] + pwmData['C'][lineIndex] + pwmData['G'][lineIndex] + pwmData['T'][lineIndex]
        
        record.pwmatrix.append( { 'A': pwmData['A'][lineIndex] / totalCount, \
                                  'C': pwmData['C'][lineIndex] / totalCount, \
                                  'G': pwmData['G'][lineIndex] / totalCount, \
                                  'T': pwmData['T'][lineIndex] / totalCount} )

    return record

def marthaBulykMatrixTableParser(mbMT):

    lines = []

    motifIndex = 1
    
    while True:

        line = mbMT.readline()

        if line == '':
            if len(lines) != 0:
                record = marthaBulykRecordParser(lines)

                motifID = "MB" + ("%4d" % motifIndex).replace(' ', '0')
                motifIndex += 1
                record.accessionNumber = motifID
                
                record.save()

            break;
        elif line == "\n":
            continue

        if '_' in line:

            if len(lines) != 0:
                record = marthaBulykRecordParser(lines)

                motifID = "MB" + ("%4d" % motifIndex).replace(' ', '0')
                motifIndex += 1
                record.accessionNumber = motifID
                
                record.save()
                lines = []
                
            lines.append(line)
            
        elif len(lines) != 0:
            lines.append(line)




####################################################################
#
# Functions for cistrom motif database
#
#
####################################################################


class cistromMotifRecord(object):

    def __init__(self, \
                 factorName = None, \
                 positionWeightMatrix = None, \
                 species = None, \
                 source = None, \
                 PMID = None, \
                 domain = None, \
                 structureCategory = None\
                 ):
        
        self.pwmatrix = positionWeightMatrix
        self.factorName = factorName
        self.species = species
        self.source = source
        self.PMID = PMID
        self.domain = domain
        self.structureCategory = structureCategory

    def parser(self, mFile):
        
        self.species = []
        self.factorName = []
        
        while True:
            
            line = mFile.readline()
            
            if line.startswith("#Position Probability Weight Matrix:"):
                
                lineIndex = 0

                self.pwmatrix = pwMatrix()

                while True:
                    line = mFile.readline()

                    if line[0] == '#':
                        break
                    
                    items = line.split()

                    if lineIndex == 0:
                        nucleotideList = items[0:4]
                    else:
                        self.pwmatrix.append( { nucleotideList[0]:float(items[1]), \
                                                nucleotideList[1]:float(items[2]), \
                                                nucleotideList[2]:float(items[3]), \
                                                nucleotideList[3]:float(items[4]) } )
                    
                    lineIndex = lineIndex + 1

            if line.startswith("#Source"):
                line = mFile.readline()
                self.source = line.strip()
                
            if line.startswith("#Factor Name"):
                line = mFile.readline()
                
                line = line.strip()
                
                items =  line.split(",")
                for item in items:
                    item = item.strip()
                    if item != "":
                        self.factorName.append(item)

            if line.startswith("#Species:"):
                line = mFile.readline()
                
                line = line.strip()
                
                items =  line.split(",")
                for item in items:
                    item = item.strip()
                    if item != "":
                        self.species.append(item)


            if line.startswith("#Domain:"):
                line = mFile.readline()
                
                self.domain = line.strip()

            if line.startswith("#Structure Category:"):
                line = mFile.readline()
                
                self.structureCategory = line.strip()

            if not line:
                break



def motifSeqGenerator(pwmatrix, number):

    """
    Genertor to generate DNA sequances to match the position weight matrix
    """

    baseNum = []

    amplifiedFactor = number
    
    for position in range( len(pwmatrix) ):
        baseNum.append([ int(pwmatrix.getProbabilityWeight(position, 'A')*amplifiedFactor), \
                         int(pwmatrix.getProbabilityWeight(position, 'C')*amplifiedFactor), \
                         int(pwmatrix.getProbabilityWeight(position, 'G')*amplifiedFactor), \
                         int(pwmatrix.getProbabilityWeight(position, 'T')*amplifiedFactor)
                         ])

    for position in range( len(baseNum) ):
        total = sum(baseNum[position])
        if total != amplifiedFactor:
            delta = amplifiedFactor - total
            maxBase = max(baseNum[position])
            baseNum[position][baseNum[position].index(maxBase)] += delta

    for index in range(number):

        seq = ''

        for position in range( len(baseNum) ):
            
            if baseNum[position][0] > 0:
                seq += 'A'
                baseNum[position][0] -= 1
            elif baseNum[position][1] > 0:
                seq += 'C'
                baseNum[position][1] -= 1
            elif baseNum[position][2] > 0:
                seq += 'G'
                baseNum[position][2] -= 1
            elif baseNum[position][3] > 0:
                seq += 'T'
                baseNum[position][3] -= 1
            else:
                raise Exception
            
        yield(seq)


def makeMotifSeqFasta(mFile, seqNumber, output):

    """
    Make sequence fasta file that satisfies the position weight matrix of the motif
    This fasta file can be used to make weblogo of the motif
    """
    
    record = cistromMotifRecord()

    record.parser(mFile)

#    print len(record.pwmatrix)

#    print record.pwmatrix

    index = 0
    for seq in motifSeqGenerator(record.pwmatrix, seqNumber):

        output.write("> %d\n" % index)
        output.write("%s\n" % seq)

        index += 1

def makeMotifWebLogo(mFileNameList, seqNumber):

    for mFileName in mFileNameList:
        mFileName = mFileName.strip()
        print mFileName
        mFile = open(mFileName, 'r')
        output = open("seq.fa", 'w')

        makeMotifSeqFasta(mFile, 10000, output)

        output.flush()
        output.close()

        comStr = ("seqlogo -n -Y -p -c -a -B 2 -F PNG -f seq.fa -o %s" % mFileName.split('.')[0] )

        os.system(comStr)
        
