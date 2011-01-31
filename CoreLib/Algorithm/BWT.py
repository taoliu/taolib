# Time-stamp: <2009-05-08 12:18:16 Tao Liu>

"""Module Description: Burrows-Wheeler transform

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

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
def bwt(s):
    s = s + '\0'
    return ''.join([x[-1] for x in                             # Take last character
                    sorted([s[i:] + s[:i] for i in range(len(s))])])    # Of each sorted row

def ibwt(s):
    L = [''] * len(s)                                          # Make empty table
    for i in range(len(s)):
        L = sorted([s[i] + L[i] for i in range(len(s))])         # Sort rows after adding s as first column
    return [x for x in L if x.endswith('\0')][0][:-1]          # Get row ending with null, [:-1] removes the null
