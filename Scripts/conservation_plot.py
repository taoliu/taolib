#!/usr/bin/env python
# Time-stamp: <2009-11-18 17:48:58 Tao Liu>

"""Description: Draw correlation plot for many wiggle files.

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
import logging
from optparse import OptionParser
import urllib2
import subprocess
import tempfile
import gzip
from taolib.CoreLib.Parser import WiggleIO, BedIO
from taolib.CoreLib.BasicStat.Func import * 

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )


hgWiggle = 'hgWiggle'                   # hgWiggle executable in the system

# ------------------------------------
# Misc functions
# ------------------------------------

error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog <-d path> [options] <bed files> ..."
    description = "Draw conservation plot for many bed files."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option('-H','--height', dest='height',type='int',default=10, help="height of plot")
    optparser.add_option('-W','--width',dest='width',type='int',default=10, help="width of plot")
    optparser.add_option('-w',dest='w',type='int',default=1000, help="window width centered at middle of bed regions,default: 1000")    
    optparser.add_option('-t','--title',dest='title',help="title of the figure. Default: 'Average Phastcons around the Center of Sites'",default= 'Average Phastcons around the Center of Sites')
    optparser.add_option('-d','--phasdb',dest='phasdb',help= 'The directory to store phastcons scores in the server')
    optparser.add_option("-l","--bed-label",dest="bedlabel",type="string",action="append",
                         help="the BED file labels in the figure. No space is allowed. This option should be used same times as -w option, and please input them in the same order as BED files. default: will use the BED file filename as labels.")
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    
    (options,bedfiles) = optparser.parse_args()

    bedfiles = map(os.path.abspath,bedfiles)
    bedfilenames = map(os.path.basename,bedfiles)

    bedfilenum = len(bedfiles)

    if bedfilenum < 1 or not options.phasdb:
        optparser.print_help()
        sys.exit(1)

    if options.bedlabel and len(options.bedlabel) == bedfilenum:
        bedlabel = options.bedlabel
    else:                               # or use the filename
        bedlabel = map(lambda x:os.path.basename(x),bedfiles)

    if options.height < 10:
        error("Height can not be lower than 10!")
        sys.exit(1)
    if options.width < 10:
        error("Width can not be smaller than 10!")
        sys.exit(1)

    # check the files
    for f in bedfiles:
        if not os.path.isfile(f):
            error("%s is not valid!" % f)
            sys.exit(1)

    # check phastcons db
    if not os.path.isdir(options.phasdb):
        error("%s is not valid!" % options.phasdb)
        sys.exit(1)

    # change wd to phastcons db path
    olddir = os.path.abspath('.')
    os.chdir(options.phasdb)

    phas_chrnames = []
    
    files_phasdb = os.listdir('.')
    for file_phasdb in files_phasdb:
        if file_phasdb.endswith('.wib'):
            name = file_phasdb.rstrip('.wib')
            phas_chrnames.append(name)

    if not phas_chrnames:
        error("%s has no valid phastcons db wib&wig files!" % options.phasdb)
        sys.exit(1)
        
    info("number of bed files: %d" % bedfilenum)

    avgValues = []

    # for each bed file
    for f in bedfiles:
        info("extract phastcons scores using %s" % f)
        scores = extract_phastcons(f,phas_chrnames, options.w)
        avgValues.append(scores)
    makeBmpFile(avgValues,olddir,options.height,options.width,options.w,options.title,bedlabel)

def extract_phastcons ( bedfile, phas_chrnames, width ):
    """Extract phastcons scores from a bed file.

    Return the average scores
    """
    info("read bed file...")
    bfhd = open(bedfile)
    bed = BedIO.parse_BED(bfhd)

    # calculate the middle point of bed regions then extend left and right by 1/2 width
    bchrs = bed.peaks.keys()
    bchrs.sort()

    chrs = []

    for c in phas_chrnames:
        if c in bchrs:
            chrs.append(c)

    sumscores = [0]*width
    n = 0
    tmpfname = tempfile.mkstemp(prefix="consplotscore")[1]
    tmpbedfname = tempfile.mkstemp(prefix="consplotbed")[1]
    # fix regions in bed file
    for chrom in chrs:
        pchrom = bed.peaks[chrom]
        for i in range(len(pchrom)):
            mid = int((pchrom[i][0]+pchrom[i][1])/2)
            left = int(mid - width/2)
            right = int(mid + width/2)

            if left < 0:
                pchrom[i] = (0,width,1,1,1,1,1,1)
            else:
                pchrom[i] = (left,right,1,1,1,1,1,1)

    bedfhd = open (tmpbedfname,'w')
    bedfhd.write(bed.tobed())
    bedfhd.close()

    for chrom in chrs:
        tmpf = open(tmpfname,'w')
        info ("extract chromosome %s" % (chrom))
        p = subprocess.Popen([hgWiggle,'-bedFile=%s' % tmpbedfname,chrom],stdout=tmpf)
        p.communicate()
        tmpf.close()

        wio = WiggleIO.WiggleIO(open(tmpfname))
        wtrack = wio.build_wigtrack()
        wtrack.sort()

        scores = bed.extract_wiggle_pv(wtrack)
        add_scores ( sumscores, scores )
        n +=  len(scores)
    os.unlink(tmpfname)
    os.unlink(tmpbedfname)
    # calculate average score
    return map(lambda x:float(x)/n, sumscores)
        
    
def add_scores ( sumscores, scores ):
    """Calculate the avg scores for each positions in a width window.
    
    """
    n = len(scores)
    for score_win in scores:
        for (p,s) in score_win:
            sumscores[p-1] += s
    return True
    
def makeBmpFile(avgValues, wd, h,w, width, title, bedlabel):
    
    #creating R file in which to write the rscript which defines the correlation plot
    #create and save the file in the current working directory

    fileName = os.path.join(wd, 'tmp')
    rFile = open(fileName+'.R','w')
    bmpname = fileName+'.bmp'
    rscript = 'sink(file=file("/dev/null", "w"), type="message")\n'
    rscript += 'sink(file=file("/dev/null", "w"), type="output")\n'    
    rscript += 'bitmap("%s",height=%d,width=%d)\n' %(bmpname,h,w)
    xInfo = range(int(-width/2),int(width/2))
    rscript += 'x<-c('+','.join(map(str,xInfo[:-1]))+')\n' # throw the last point which may be buggy
    for i in range(len(avgValues)):
        avgscores = avgValues[i]
        tmpname = 'y'+str(i)
        rscript += tmpname+'<-c('+','.join(map(str,avgscores[:-1]))+')\n' # throw the last point which may be buggy

    tmplist = []
    for i in range(len(avgValues)):
        tmplist.append( "y%d" % i )
    
    rscript += "ymax <- max("+ ",".join(tmplist) +")\n"
    rscript += "ymin <- min("+ ",".join(tmplist) +")\n"    
    rscript += "yquart <- (ymax-ymin)/4\n"

    rscript += 'plot(x,y0,type="l",col=rainbow(%d)[1],main=\"%s\",xlab="Relative Distance from the Center (bp)",ylab="Average Phastcons",ylim=c(ymin-yquart,ymax+yquart))\n' % (len(avgValues),title)
    for i in range(1,len(avgValues)):
        rscript += 'lines(x,y'+str(i)+',col=rainbow(%d)[%d])\n' % (len(avgValues),i+1)
    rscript += 'abline(v=0)\n'
    legend_list = map(lambda x:"'"+x+"'", bedlabel)
    rscript += 'legend("topright",c(%s),col=rainbow(%d),lty=c(%s))\n' % (','.join(legend_list),len(avgValues),','.join(['1']*len(avgValues)))
        
    rscript += 'dev.off()\n'
    rFile.write(rscript)
    rFile.close()
    #executing the R file and forming the pdf file
    data = subprocess.call(['Rscript',fileName+'.R'])

    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
