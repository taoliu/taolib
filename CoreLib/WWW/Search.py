# Time-stamp: <2009-05-08 11:39:35 Tao Liu>

"""Module Description: Search Engine Modules

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
from taolib.CoreLib.DB import dbRawChIP, dbMotif

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def cgi_search_DBRawChIP (dbfile, query):
    db = dbRawChIP.DBRawChIP(file=dbfile)
    where_clause = ""
    if query:
        terms = parse_query_text(query)
        where_clause = "where "+" and ".join( map(lambda x:str(x[0])+" like \"%%"+str(x[1])+"%%\"", terms) )
        
    c = db.execute('select id,factor,organism,cell,condition,platform,lab,pmid from chips %s order by id' % where_clause)
    return c

def parse_query_text ( text ):
    """Return a dictionary with keywords and values paires.
    
    """
    keywds = re.findall('\[\w+\]',text)
    fields = re.split('\[\w+\]',text)
    keywds.append('[any]')
    return zip(keywds,fields)
    


# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def test():
    #test function
    pass

if __name__ == '__main__':
    test()
