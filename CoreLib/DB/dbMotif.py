# Time-stamp: <2008-08-05 17:19:55 Tao Liu>

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

import sys
import re
import sqlite3

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Class
# ------------------------------------

class DBMotif(object):

    def __init__(self, file = None):
        """initialize the connection.
        
        All the data will be saved in @parameter file.

        If file is ':memory:', all data will be saved in memory.
        """
        self.connection = sqlite3.connect(file)

    def create_table (self, tablename = 'motif'):
        """Create table for motif data.

        Table name: motif (default)
        
        columns: matrix, source, factorName, species, pmid, domain, structureCategory

        The pmid field is integer type, whereas others are text type.
        """

        c = self.connection.cursor()

        # create

        c.execute('''create table ?  (matrix text, source text, factorName text, species
        text, pmid integer, domain text, structureCategory text )''', (tablename))

        self.connection.commit()

        c.close()

    def add (self, matrix, source, factorName, species, pmid, domain, structureCategory):
        c = self.connection.cursor()
        t = (matrix, source, factorName, species, pmid, domain, structureCategory)
        c.execute ('insert into chips values (?,?,?,?,?,?,?)', t)
        self.connection.commit()
        c.close()

    def add_from_motif_text (self, fhd ):
        """Add motif entry from dotmotif file.

        fhd : file handler
        """
        # read the whole file
        text = fhd.read()
        parts = text.split('#')

    def execute (self, string):
        c = self.connection.cursor()
        c.execute(string)
        return c

