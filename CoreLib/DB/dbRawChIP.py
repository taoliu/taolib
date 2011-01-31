# Time-stamp: <2009-05-08 11:29:39 Tao Liu>

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

import sqlite3
import csv
from taolib.CoreLib.DB.NDD import *
# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# functions
# ------------------------------------

# ------------------------------------
# Class
# ------------------------------------
# class DBRawChIP:
    
#     def __init__ (self,file=None,url=None):
#         """initialize the connection.
        
#         Note, only 'file' parameter is available. Once it is set, all the
#         data will be saved in that file.

#         If file is ':memory:', all data will be saved in memory.
#         """
#         self.conn = sqlite3.connect(file)

#     def create_table (self, tablename='chips'):
#         """Create table for RawChIP data.

#         Table name: chips (default)
#         columns: id, factor, organism, cell, condition,
#         platform, lab, pmid, rawdata, output. and any

#         The id and pmid field is integer type, whereas others are text
#         type.
#         """
#         c = self.conn.cursor()
#         # create
#         c.execute("create table %s (id integer, factor text, organism text, cell text, condition text, platform text, lab text, pmid integer, rawdata text, output text, any text)" % tablename)
#         self.conn.commit()
#         c.close()

#     def add_chip (self,id,factor, organism, cell,condition,platform, lab, pmid, rawdata, output):
#         """Add a chIP entry by giving all its information including:

#         1. id       : id number
#         2. factor   : TF
#         3. organism : the assembly like hg18 or mm9
#         4. cell     : cell type
#         5. condition:
#         6. platform : platform for array, like affy WG version 2 or Nimblegen Encode
#         7. lab      :
#         8. pmid     : pubmed id
#         9. rawdata  : directory pathname of the raw chIP data
#         10. output  : directory pathname for the peak data (bed and wiggle)
#         """
#         c = self.conn.cursor()
#         t = [id,factor,organism,cell,condition,platform,lab,pmid,rawdata,output]
#         t2 = " ".join(map(str,t))
#         t.append(t2)
#         c.execute ('insert into chips values (?,?,?,?,?,?,?,?,?,?,?)', t)
#         self.conn.commit()
#         c.close()

#     def add_from_table (self, fhd):
#         """Batchly add entries from excel-tab-delimited file.

#         parameter:
#         fhd : the file handler

#         Note: file must have a header line for every columns defined
#         in add_chip() function.
#         """
        
#         reader = csv.DictReader(fhd,dialect='excel-tab')
#         c = self.conn.cursor()

#         for row in reader:
#             t = [row['id'],row['factor'], row['organism'], row['cell'], row['condition'], row['platform'], row['lab'], row['pmid'], row['rawdata'], row['output']]
#             t2 = " ".join(map(str,t))
#             t.append(t2)
#             c.execute ('insert into chips values (?,?,?,?,?,?,?,?,?,?,?)',t)
#         self.conn.commit()
#         c.close()

#     def execute ( self,string ):
#         """To execute a SQL command.
        
#         """
#         c = self.conn.cursor()
#         c.execute(string)
#         return c
    

#     def id2rawpath (self,rawid):
#         """id to rawdata pathname matcher.
        
#         """
#         c = self.execute('select rawdata from chips where id==%d' % (rawid))
#         for r in c:
#             return r[0]

#     def id2processedpath (self,rawid):
#         """id to processed data filename matcher.
        
#         """
#         c = self.execute('select output from chips where id==%d' % (rawid))
#         for r in c:
#             return r[0]


# # another database

class DBRawData:
    
    def __init__ (self,file=None,url=None):
        """initialize the connection.
        
        Note, only 'file' parameter is available. Once it is set, all the
        data will be saved in that file.

        If file is ':memory:', all data will be saved in memory.
        """
        self.conn = sqlite3.connect(file)

    def init_tables (self):
        """Create table for RawChIP data. Run at the first time,
        otherwise, the old database file will be overwritten.

        12 tables needed.
        """
        c = self.conn.cursor()
        # create 11 reference tables
        c.execute("create table Factor (id text NOT NULL PRIMARY KEY, comment blob)")
        c.execute("create table Organism (id text NOT NULL PRIMARY KEY, comment blob);")
        c.execute("create table Cell (id text NOT NULL PRIMARY KEY, comment blob);")
        c.execute("create table Platform (id text NOT NULL PRIMARY KEY, comment blob);")
        c.execute("create table Status (id text NOT NULL PRIMARY KEY, comment blob);")
        c.execute("create table Lab (id text NOT NULL PRIMARY KEY, comment blob);")
        c.execute("create table Paper (id text NOT NULL PRIMARY KEY, comment blob);")
        c.execute("create table Condition (id text NOT NULL PRIMARY KEY, comment blob);")
        c.execute("create table Type (id text NOT NULL PRIMARY KEY, comment blob);")
        c.execute("create table Submit (id text NOT NULL PRIMARY KEY, comment blob);")
        c.execute("create table Assembly (id text NOT NULL PRIMARY KEY, comment blob, Organismid text NOT NULL REFERENCES Organism);")
        c.execute("""
create table dataset
(id text NOT NULL PRIMARY KEY,
Factorid text NOT NULL REFERENCES Factor,
Organismid text NOT NULL REFERENCES Organism,
Assemblyid text NOT NULL REFERENCES Assembly,
Cellid text NOT NULL REFERENCES Cell,
Platformid text NOT NULL REFERENCES Platform,
Statusid text NOT NULL REFERENCES Status,
Labid text NOT NULL REFERENCES Lab,
Paperid text NOT NULL REFERENCES Paper,
Conditionid text NOT NULL REFERENCES Condition,
Typeid text NOT NULL REFERENCES Type,
Submitid text NOT NULL REFERENCES Submit,
Raw blob,
XLS blob,
Wiggle blob,
GEO blob,
PubMED integer(20),
Comment blob,
MD5 text,
Timestamp text)
""")
        c.close()

    def add_NDDdata (self,ndddata):
        """Add a chIP entry from a NDDdata object.
        """
        assert isinstance(ndddata, NDDdata)
        c = self.conn.cursor()
        
        factor = ndddata.get('Factor')
        if not c.execute('select id from Factor where id="%s"' % (factor)).fetchall():
            c.execute('insert into Factor values ("%s","")' % (factor))

        organism = ndddata.get('Organism')
        if not c.execute('select id from Organism where id="%s"' % (organism)).fetchall():
            c.execute('insert into Organism values ("%s","")' % (organism))

        cell = ndddata.get('Cell')
        if not c.execute('select id from Cell where id="%s"' % (cell)).fetchall():
            c.execute('insert into Cell values ("%s","")' % (cell))

        platform = ndddata.get('Platform')
        if not c.execute('select id from Platform where id="%s"' % (platform)).fetchall():
            c.execute('insert into Platform values ("%s","")' % (platform))

        status = ndddata.get('Status')
        if not c.execute('select id from Status where id="%s"' % (status)).fetchall():
            c.execute('insert into Status values ("%s","")' % (status))

        lab = ndddata.get('Lab')
        if not c.execute('select id from Lab where id="%s"' % (lab)).fetchall():
            c.execute('insert into Lab values ("%s","")' % (lab))

        paper = ndddata.get('Paper')
        if not c.execute('select id from Paper where id="%s"' % (paper)).fetchall():
            c.execute('insert into Paper values ("%s","")' % (paper))

        condition = ndddata.get('Condition')
        if not c.execute('select id from Condition where id="%s"' % (condition)).fetchall():
            c.execute('insert into Condition values ("%s","")' % (condition))

        ty = ndddata.get('Type')
        if not c.execute('select id from Type where id="%s"' % (ty)).fetchall():
            c.execute('insert into Type values ("%s","")' % (ty))

        submit = ndddata.get('Submit')
        if not c.execute('select id from Submit where id="%s"' % (submit)).fetchall():
            c.execute('insert into Submit values ("%s","")' % (submit))

        assem = ndddata.get('Assembly')
        if not c.execute('select id from Assembly where id="%s"' % (assem)).fetchall():
            c.execute('insert into Assembly values ("%s","","%s")' % (assem,organism))

        ids = c.execute('select id from dataset').fetchall()
        lastnum = max([int(x) for (x,) in ids])
        
        c.execute('insert into dataset values (%d,"%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s")' %
                  (lastnum+1,factor,organism,assem,cell,platform,status,lab,paper,
                   condition,ty,submit,
                   ndddata.get('Raw'),
                   ndddata.get('XLS'),
                   ndddata.get('Wiggle'),
                   ndddata.get('GEO'),
                   ndddata.get('PubMED'),
                   ndddata.get('Comment'),
                   ndddata.md5key,
                   ndddata.timestamp,
                   )
                  )
        self.conn.commit()
        c.close()
        return lastnum+1

    def update_NDDdata (self, iid, ndddata):
        """Add a chIP entry from a NDDdata object.
        """
        assert isinstance(ndddata, NDDdata)
        c = self.conn.cursor()
        
        factor = ndddata.get('Factor')
        organism = ndddata.get('Organism')
        cell = ndddata.get('Cell')
        platform = ndddata.get('Platform')
        status = ndddata.get('Status')
        lab = ndddata.get('Lab')
        paper = ndddata.get('Paper')
        condition = ndddata.get('Condition')
        ty = ndddata.get('Type')
        submit = ndddata.get('Submit')
        assem = ndddata.get('Assembly')
        c.execute('update dataset set Factorid = "%s", Organismid = "%s", Assemblyid = "%s", Cellid = "%s", Platformid = "%s", Statusid = "%s", Labid = "%s", Paperid = "%s", Conditionid = "%s", Typeid = "%s", Submitid = "%s",Raw = "%s",XLS = "%s", Wiggle = "%s", GEO = "%s", PubMED = "%s", Comment = "%s" where id=%s' %
                  (factor,organism,assem,cell,platform,status,lab,paper,
                   condition,ty,submit,
                   ndddata.get('Raw'),
                   ndddata.get('XLS'),
                   ndddata.get('Wiggle'),
                   ndddata.get('GEO'),
                   ndddata.get('PubMED'),
                   ndddata.get('Comment'),
                   str(iid) ))
        
        self.conn.commit()
        c.close()
        return True

    def add_NDDfile (self, filename):
        fhd = open(filename,'r')
        data = parse(fhd)
        return self.add_NDDdata(data)

    def execute ( self,string ):
        """To execute a SQL command.
        
        """
        c = self.conn.cursor()
        c.execute(string)
        return c

    def delete (self, i):
        try:
            self.execute ('delete from dataset where id=\"%s\"' % ((str(i))))
            self.conn.commit()
        except:
            return False
        else:
            return True

    def select (self, i):
        d = self.execute ('select * from dataset where id=%s' % (str(i))).fetchall()[0]
        ndddata = NDDdata()
        ndddata.set('Index',d[0])
        ndddata.set('Factor',d[1])
        ndddata.set('Organism',d[2])
        ndddata.set('Assembly',d[3])
        ndddata.set('Cell',d[4])
        ndddata.set('Platform',d[5])
        ndddata.set('Status',d[6])
        ndddata.set('Lab',d[7])
        ndddata.set('Paper',d[8])
        ndddata.set('Condition',d[9])
        ndddata.set('Type',d[10])
        ndddata.set('Submit',d[11])
        ndddata.set('Raw',d[12])
        ndddata.set('XLS',d[13])
        ndddata.set('Wiggle',d[14])
        ndddata.set('GEO',d[15])
        ndddata.set('PubMED',d[16])
        ndddata.set('Comment',d[17])
        ndddata.md5key = d[18]
        ndddata.timestamp = d[19]
        #print ndddata
        return ndddata
        
