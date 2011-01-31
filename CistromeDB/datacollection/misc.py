# Time-stamp: <2009-01-22 17:44:11 Tao Liu>

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

import sys
import os
import mimetypes
import gzip
import bz2
# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
def combine (sdir,tfile,suffix=""):
    fs = os.walk(sdir)
    tfhd = open(tfile,"wb")
    for path,dirs,files in fs:
        for ff in files:
            f = os.path.join(path,ff)
            if f.endswith(suffix):
                if mimetypes.guess_type(f)[1] == 'gzip':
                    sfhd = gzip.open(f,'rb')
                    cached = sfhd.read(100000)
                    while cached:
                        tfhd.write(cached)
                        cached = sfhd.read(100000)
                    sfhd.close()
                elif mimetypes.guess_type(f)[1] == 'bzip2':
                    sfhd = bz2.BZ2File(f,'rb')
                    cached = sfhd.read(100000)
                    while cached:
                        tfhd.write(cached)
                        cached = sfhd.read(100000)
                    sfhd.close()
                elif not mimetypes.guess_type(f)[0] and not mimetypes.guess_type(f)[1]:
                    sfhd = file(f,'rb')
                    cached = sfhd.read(100000)
                    while cached:
                        tfhd.write(cached)
                        cached = sfhd.read(100000)                        
                    sfhd.close()
                else:
                    raise Exception("Invalid Format!!\n")
    tfhd.close()


# ------------------------------------
# Classes
# ------------------------------------
