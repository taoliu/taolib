# Time-stamp: <2009-09-22 17:24:16 Tao Liu>

"""Module Description: IO for BED

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

import sys
import re
from taolib.CoreLib.FeatIO import PeakIO

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
def parse_BED (fhd):
    """Parse a tab-delimited bed file

    Return a PeakIO object containing peak regions.
    """
    peaks = PeakIO()
    for thisline in fhd:
        thisline = thisline.rstrip()
        if not thisline: continue #return ("blank",None,None)
        if thisline.startswith("#"): continue #return ("comment line",None,None) # comment line is skipped
        if thisline.startswith("track"): continue
        if thisline.startswith("browser"): continue        
        thisfields = thisline.split()
        startpos = max(0,int(thisfields[1]))

        peaks.add(thisfields[0],startpos,int(thisfields[2]),1,1,1,1,1,1)
    return peaks

# ------------------------------------
# Classes
# ------------------------------------
