#!/usr/bin/env python
# Time-stamp: <2010-09-08 02:38:38 Tao Liu>

"""Module Description

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

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

import os
import sys
import re
import csv
import logging
from optparse import OptionParser
import reportlab
import Bio

from taolib.CoreLib.FeatIO import WigTrackI
from taolib.CoreLib.BasicStat.Func import mean,median,std

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )
error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info
# ------------------------------------
# Misc functions
# ------------------------------------

def andfilter ( cvsfile, write_func, *args ):
    """
    
    """
    argv = args[0]
    if len(argv) < 2:
        sys.stderr.write("Need two extra arguments for 'organize', e.g. command: <1,2,3> <4,5,6> means the first 1,2,3 will be used as dependent variables/response, and 4,5,6 will be used as independent variables/terms/predictors.\n")
        sys.exit()

    responses_num = map(int,argv[0].split(","))
    predictors_num = map(int,argv[1].split(","))
    
    fields = cvsfile.fieldnames
    responses_label = map(lambda x:"res."+fields[x],responses_num)
    predictors_label = map(lambda x:"pre."+fields[x],predictors_num)

    responses_name = map(lambda x:fields[x],responses_num)
    predictors_name = map(lambda x:fields[x],predictors_num)
    
    #write_func( "#%s\t%s\n" \
    #            % ( ",".join(map( lambda x:str(x[0])+":"+str(x[1]) , zip(responses_num,responses_name) )),
    #                ",".join(map( lambda x:str(x[0])+":"+str(x[1]) , zip(predictors_num,predictors_name) ))) )
    write_func( "%s\t%s\n" \
                % ( ",".join(map( lambda x:str(x) , responses_name )),
                    ",".join(map( lambda x:str(x) , predictors_name ))) )

    for l in cvsfile:
        # for responses
        t_str_list = []
        for t in responses_name:
            t_str_list.append(l.setdefault(t,"NA"))
        # for predictors
        v_str_list = []
        for v in predictors_name:
            v_str_list.append(l.setdefault(v,"NA"))

        write_func( "\t".join( (",".join(t_str_list),",".join(v_str_list)) ) )
        write_func( "\n" )


def combcall2draw ( cvsfile, write_func, *args ):
    """User specifies several columns to consider, this tool will call
    regions where either of the column is above its threshold.
    
    """
    argv = args[0]
    if len(argv) < 6:
        sys.stderr.write("Need 6 extra arguments for 'combcall2draw', options <loc column> <score column1[,score column2,...]> <cutoff1[,cutoff2,cutoff3]> <min length> <max gap> <pdf filename>\ne.g. command: <0> <1,2,3> <0.5,0.6,0.7> <10000> <2000> <a.pdf>, means to use the first column as genome coordinations to call enriched regions from the combinition of #1, #2 and #3, the thresholds to call enriched region are 0.5 for column 1, 0.6 for column 2 and 0.7 for column 3, the minimum length of region is 10k, and the maximum gap to link two nearby regions is 2k. Then the figure will be saved in a.pdf.\n")
        sys.exit()        
    cor_column = cvsfile.fieldnames[int(argv[0])]
    var_columns = map(lambda x:cvsfile.fieldnames[int(x)],argv[1].split(","))
    cutoffs = map(float,argv[2].split(","))
    
    min_len = int(argv[3])
    max_gap = int(argv[4])
    wtrack = WigTrackI()                 # combined track containing 1 if either of track is above cutoff
    add_func = wtrack.add_loc
    
    for l in cvsfile:
        cor = l.setdefault(cor_column,None)
        if not cor or cor =="NA":
            continue
        
        for i in range(len(var_columns)):
            var_column = var_columns[i]
            cutoff = cutoffs[i]
            var = l.setdefault(var_column,None)
            if var and var != "NA" and float(var) > cutoff:
                (chrom,start,end) = cor.split(".")
                add_func(chrom,int(start),1.1)
                break
            
    wtrack.span = int(end)-int(start)
    bpeaks = wtrack.call_peaks(cutoff=1.0,min_length=min_len,max_gap=max_gap)
    #f = argv[5]
    fhd = open(argv[5].replace("pdf","bed"),"w")
    fhd.write(bpeaks.tobed())
    
    from Bio.Graphics import BasicChromosome
    from reportlab.lib.colors import gray, black, white
    entries = [("chrI", 15072419),
               ("chrII", 15279316),
               ("chrIII", 13783681),
               ("chrIV", 17493784),
               ("chrV", 20919398),
               ("chrX", 17718852)]
    max_length = max([x[1] for x in entries])
    chr_diagram = BasicChromosome.Organism()
    for name, length in entries:
        cur_chromosome = BasicChromosome.Chromosome(name)
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
            body_regions.append( (p[1]-p[0],black) ) # enriched regions
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
        
    chr_diagram.draw(argv[5], "Highlight regions in Caenorhabditis elegans" )


def call1draw ( cvsfile, write_func, *args ):
    """Call regions, then plot it in chromosome figure.

    A combination of drawchrom and call1
    
    """
    argv = args[0]
    if len(argv) < 6:
        sys.stderr.write("Need 6 extra arguments for 'call1draw', options <loc column> <score column> <cutoff> <min length> <max gap> <pdf filename>\ne.g. command: <0> <1> <0.5> <10000> <2000> <a.pdf>, means to use the first column as genome coordinations to call enriched regions from the second column, the threshold to call enriched region is 0.5, the minimum length of region is 10k, and the maximum gap to link two nearby regions is 2k. Then the figure will be saved in a.pdf.\n")
        sys.exit()
    cor_column = cvsfile.fieldnames[int(argv[0])]
    var_column = cvsfile.fieldnames[int(argv[1])]
    cutoff = float(argv[2])
    min_len = int(argv[3])
    max_gap = int(argv[4])
    wtrack = WigTrackI()
    add_func = wtrack.add_loc
    for l in cvsfile:
        cor = l.setdefault(cor_column,None)
        var = l.setdefault(var_column,None)        
        if cor and var and cor != "NA" and var != "NA":
            (chrom,start,end) = cor.split(".")
            add_func(chrom,int(start),float(var))
    wtrack.span = int(end)-int(start)
    bpeaks = wtrack.call_peaks(cutoff=cutoff,min_length=min_len,max_gap=max_gap)
    fhd = open(argv[5].replace("pdf","bed"),"w")
    fhd.write(bpeaks.tobed())
    
    from Bio.Graphics import BasicChromosome
    from reportlab.lib.colors import gray, black, white
    entries = [("chrI", 15072419),
               ("chrII", 15279316),
               ("chrIII", 13783681),
               ("chrIV", 17493784),
               ("chrV", 20919398),
               ("chrX", 17718852)]
    max_length = max([x[1] for x in entries])
    chr_diagram = BasicChromosome.Organism()
    for name, length in entries:
        cur_chromosome = BasicChromosome.Chromosome(name)
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
            body_regions.append( (p[1]-p[0],black) ) # enriched regions
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
        
    chr_diagram.draw(argv[5], "%s regions in Caenorhabditis elegans" % (var_column) )

        
def drawchrom ( cvsfile, write_func, *args ):
    """Draw CE chromosome tool.

    Doesn't need any parameters.
    """
    from Bio.Graphics import BasicChromosome
    from reportlab.lib.colors import gray, black
    entries = [("chrI", 15072419),
               ("chrII", 15279316),
               ("chrIII", 13783681),
               ("chrIV", 17493784),
               ("chrV", 20919398),
               ("chrX", 17718852)]
    max_length = max([x[1] for x in entries])
    chr_diagram = BasicChromosome.Organism()
    for name, length in entries:
        cur_chromosome = BasicChromosome.Chromosome(name)
        #Set the length, adding and extra 20 percent for the tolomeres:
        cur_chromosome.scale_num = max_length * 1.1
        # Add an opening telomere
        start = BasicChromosome.TelomereSegment()
        start.scale = 0.05 * max_length
        start.fill_color=black
        cur_chromosome.add(start)
        #Add a body - using bp as the scale length here.
        body = BasicChromosome.ChromosomeSegment()
        body.fill_color=gray
        body.scale = length
        cur_chromosome.add(body)
        #Add a closing telomere
        end = BasicChromosome.TelomereSegment(inverted=True)
        end.scale = 0.05 * max_length
        end.fill_color=black
        cur_chromosome.add(end)
        #This chromosome is done
        chr_diagram.add(cur_chromosome)

    chr_diagram.draw("simple_chrom.pdf", "Caenorhabditis elegans" )

def summary ( cvsfile, write_func, *args ):
    """Show the column names.
    """
    fsnames = cvsfile.fieldnames
    data_dict = {}
    for f in fsnames:
        data_dict[f]=[]
    #print "\n".join(map( lambda x:":".join(map(str,x)) ,enumerate(fsnames))  )
    for l in cvsfile:
        for f in fsnames:
            v = l.setdefault(f,None)
            if v and v!="NA":
                data_dict[f].append(v)
    write_func( "colnum:colname\tsum,mean,median,std,cutoff\n" )
    for (i,f) in enumerate(fsnames):
        try:
            v_array = map(float,data_dict[f])
            v_sum = "%.2f" % sum(v_array)
            v_mean = "%.2f" % mean(v_array)
            v_median = "%.2f" % median(v_array)
            v_std = "%.2f" % std(v_array, float(v_mean))
            v_cutoff = "%.2f" % (float(v_mean)+float(v_std))
        except ValueError:
            (v_sum,v_mean,v_median,v_std,v_cutoff)=["NA"]*5
        write_func( "%d:%s\t%s,%s,%s,%s,%s\n" % (i,f,v_sum,v_mean,v_median,v_std,v_cutoff ))

def organize ( cvsfile, write_func, *args ):
    """Re-organize the columns for data-mining.
    """
    argv = args[0]
    if len(argv) < 2:
        sys.stderr.write("Need two extra arguments for 'organize', e.g. command: <1,2,3> <4,5,6> means the first 1,2,3 will be used as dependent variables/response, and 4,5,6 will be used as independent variables/terms/predictors.\n")
        sys.exit()

    responses_num = map(int,argv[0].split(","))
    predictors_num = map(int,argv[1].split(","))
    
    fields = cvsfile.fieldnames
    responses_label = map(lambda x:"res."+fields[x],responses_num)
    predictors_label = map(lambda x:"pre."+fields[x],predictors_num)

    responses_name = map(lambda x:fields[x],responses_num)
    predictors_name = map(lambda x:fields[x],predictors_num)
    
    #write_func( "#%s\t%s\n" \
    #            % ( ",".join(map( lambda x:str(x[0])+":"+str(x[1]) , zip(responses_num,responses_name) )),
    #                ",".join(map( lambda x:str(x[0])+":"+str(x[1]) , zip(predictors_num,predictors_name) ))) )
    write_func( "%s\t%s\n" \
                % ( ",".join(map( lambda x:str(x) , responses_name )),
                    ",".join(map( lambda x:str(x) , predictors_name ))) )

    for l in cvsfile:
        # for responses
        t_str_list = []
        for t in responses_name:
            t_str_list.append(l.setdefault(t,"NA"))
        # for predictors
        v_str_list = []
        for v in predictors_name:
            v_str_list.append(l.setdefault(v,"NA"))

        write_func( "\t".join( (",".join(t_str_list),",".join(v_str_list)) ) )
        write_func( "\n" )
    

def call1 (cvsfile, write_func, *args ):
    """Call enrich regions from certain column
    """
    argv = args[0]
    if len(argv) < 5:
        sys.stderr.write("Need 5 extra arguments for 'call', options <loc column> <score column> <cutoff> <min length> <max gap>\ne.g. command: <0> <1> <0.5> <10000> <2000>, means to use the first column as genome coordinations to call enriched regions from the second column, the threshold to call enriched region is 0.5, the minimum length of region is 10k, and the maximum gap to link two nearby regions is 2k.\n")
        sys.exit()
    cor_column = cvsfile.fieldnames[int(argv[0])]
    var_column = cvsfile.fieldnames[int(argv[1])]
    cutoff = float(argv[2])
    min_len = int(argv[3])
    max_gap = int(argv[4])
    wtrack = WigTrackI()
    add_func = wtrack.add_loc
    for l in cvsfile:
        cor = l.setdefault(cor_column,None)
        var = l.setdefault(var_column,None)        
        if cor and var and cor != "NA" and var != "NA":
            (chrom,start,end) = cor.split(".")
            add_func(chrom,int(start),float(var))
    wtrack.span = int(end)-int(start)
    write_func( "# regions called from %s:%s\n" % (argv[1],var_column) )
    
    bpeaks = wtrack.call_peaks(cutoff=cutoff,min_length=min_len,max_gap=max_gap)
    write_func( bpeaks.tobed() )

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog [options]"
    description = "Script to analyze C. elegans histone marks data matrix."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-i","--ifile",dest="ifile",type="string",
                         help="input file")
    optparser.add_option("-o","--ofile",dest="ofile",
                         help="output file, default: stdout") 

    (options,args) = optparser.parse_args()

    command_list = {"summary":summary,
                    "organize":organize,
                    "call1":call1,
                    "drawchrom":drawchrom,
                    "call1draw":call1draw,
                    "combcall2draw":combcall2draw,
                    }
    command_des = {"summary":"Show the column names.",
                   "organize":"Re-organize the file for data-mining.",
                   "call1":"Call enriched regions for certain column.",
                   "drawchrom":"Draw ce chromosomes.",
                   "call1draw":"Call enriched regions and then draw chromosome figures.",
                   "combcall2draw":"Call enriched regions where any of the tracks is above threshold and draw them on chromosome figures.",
                   }

    if not options.ifile or not args:
        optparser.print_help()
        sys.exit()

    if options.ofile:
        write_func = open(options.ofile,"w").write
    else:
        write_func = sys.stdout.write

    if command_list.has_key(args[0]):
        com = command_list[args[0]]
        com_args = args[1:]
    else:
        optparser.print_help()
        sys.stderr.write("Avialable Commands:\n\n")
        for c in command_list.keys():
            sys.stderr.write(c+": "+command_des[c]+"\n")
        sys.exit()

    cvsfilereader = csv.DictReader(open(options.ifile,"r"),delimiter="\t")
    # run commands
    com(cvsfilereader,write_func,com_args)
    
    

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
