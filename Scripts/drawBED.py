#!/usr/bin/env python

import sys
import urllib2
import tempfile
import gzip
import os
import sys
from reportlab.lib.colors import *

def get_chrom_length ( dbname ):
    # first try to find a local file called dbname
    try:
        fhd = open(dbname,"r")
        chrom_len = []
        for l in fhd:
            fs = l.split()
            try:
                chrom_len.append( (fs[0],int(fs[1])) )
            except:
                pass
        fhd.close()
    except:
        # get chromosome length from UCSC db download page
        f = urllib2.urlopen(UCSC_chrom_URL % (dbname))
        tmpfname = tempfile.mkstemp(prefix="drawBED")[1]
        tmpf = open(tmpfname,'w')
        tmpf.write(f.read())                # write file content to temp file
        tmpf.close()
        f.close
        # read it
        fhd = gzip.open(tmpfname,'r')
        chrom_len = []
        for l in fhd:
            fs = l.split()
            try:
                chrom_len.append( (fs[0],int(fs[1])) )
            except:
                pass
        fhd.close()
        os.unlink(tmpfname)
    return chrom_len

def getCol(colorname):
    try:
        col = toColor(colorname)
    except ValueError:
        col = toColor('black')

def main():
    from taolib.CoreLib.Parser.BedIO import parse_BED
    from Bio.Graphics import BasicChromosome

    if len(sys.argv) < 3:
        sys.stderr.write("Draw Chromosome Figure\nneed 2 paras: %s <dbname> <color> <bed file>\n" % sys.argv[0])
        sys.exit(1)

    try:
        entries = get_chrom_length(sys.argv[1])
    except:
        error("Error!")
        sys.exit(1)

    col = getCol(sys.argv[2])
    
    bpeaks = parse_BED(open(sys.argv[3]))
    pdffile = sys.argv[3]+".pdf"

    max_length = max([x[1] for x in entries])
    chr_diagram = BasicChromosome.Organism()
    for name, length in entries:
        cur_chromosome = BasicChromosome.Chromosome(name)
        cur_chromosome.title_size = 0.5
        #Set the length, adding and extra 20 percent for the tolomeres:
        cur_chromosome.scale_num = max_length * 1.1
        # Add an opening telomere
        start = BasicChromosome.TelomereSegment()
        start.scale = 0.05 * max_length
        start.fill_color=gray
        cur_chromosome.add(start)
        #Add a body - using bp as the scale length here.
        try:
            cpeaks = bpeaks.peaks[name]
        except:
            cpeaks = []
        body_regions = []
        last_pos = 0
        for p in cpeaks:
            body_regions.append( (p[0]-last_pos,white) ) # outside regions
            body_regions.append( (p[1]-p[0],col) ) # enriched regions
            last_pos = p[1]
            assert p[1] < length
        body_regions.append( (length-last_pos,white) ) # last part

        for b,c in body_regions:
            body = BasicChromosome.ChromosomeSegment()
            body.fill_color= c
            body.scale = b
            cur_chromosome.add(body)
        
        #Add a closing telomere
        end = BasicChromosome.TelomereSegment(inverted=True)
        end.scale = 0.05 * max_length
        end.fill_color=gray
        cur_chromosome.add(end)
        #This chromosome is done
        chr_diagram.add(cur_chromosome)
        
    chr_diagram.draw(pdffile, "Highlight regions" )

main()
