#!/usr/bin/env python
# Time-stamp: <2011-01-31 17:17:00 Tao Liu>

# simple normalization script for a bedGraph file from wig2bedGraphBins.py script

import re
import sys
from taolib.CoreLib.BasicStat.Func import mean,std

if len(sys.argv) < 2:
    sys.stderr.write("normalize bg file from wig2bedGraphBins.py, output zscores.\nneed 1 para: %s <bedGraph with fixed bins>\n" % sys.argv[0])
    sys.exit(1)

fhd = open(sys.argv[1])

d = []
v = []
for i in fhd:
    fs = i.rstrip().split()
    d.append( ["\t".join(fs[:3]),fs[3]] )
    if fs[3] != "NA":
        v.append(float(fs[3]))

m = mean(v)
s = std(v,n_mean=m)

for (i,j) in d:
    if j!="NA":
        print "%s\t%.6f" % (i,(float(j)-m)/s)
    else:
        print i+"\t"+j
