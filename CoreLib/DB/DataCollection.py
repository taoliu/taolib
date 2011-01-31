# Time-stamp: <2009-05-08 11:28:58 Tao Liu>

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
import MySQLdb
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
class DBDataCollection:
    """Data Collection Database. Interact with NDD objects.

    Member variables:
    1. self.conn, a connection object for DB

    Provide high-end and low-end operations:
    1. init_tables, to initialze the DB object
    2. 
    
    """
    def __init__ (self,host="localhost",user=None,passwd=None,db="DataCollection",port=3306):
        """initialize MySQL connection.
        
        Note, only 'file' parameter is available. Once it is set, all the
        data will be saved in that file.

        """
        if user and passwd:
            self.conn = MySQLdb.connect(host=host,user=user,
                                        passwd=passwd,db=db,
                                        port=port)
        elif user:
            self.conn = MySQLdb.connect(host=host,user=user,
                                        db=db,
                                        port=port)

    def init_tables (self):
        """Create table for RawChIP data. Run at the first time,
        otherwise, the old database file will be overwritten.

        12 tables needed.
        """
        c = self.conn.cursor()
        # create 11 reference tables
        c.execute("DROP TABLE dataset")
        c.execute("""
CREATE TABLE dataset
(id INTEGER(20) NOT NULL PRIMARY KEY,
Factor TEXT NOT NULL,
Organism TEXT NOT NULL,
Assembly TEXT NOT NULL,
Cell TEXT NOT NULL Cell,
Platform TEXT NOT NULL,
Status TEXT NOT NULL,
Lab TEXT NOT NULL,
Paper TEXT NOT NULL,
Condition TEXT NOT NULL,
Type TEXT NOT NULL,
Submit TEXT NOT NULL,
Raw BLOB,
XLS BLOB,
Wiggle BLOB,
GEO BLOB,
PubMED INTEGER(20),
Comment BLOB,
MD5 TEXT,
Timestamp TIMESTAMP)
""")
        c.close()

    def add_NDDdata (self,ndddata):
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
        ids = c.execute('select id from dataset').fetchall()
        lastnum = max([int(x) for (x,) in ids])
        
        c.execute('INSERT INTO dataset VALUES (%d,"%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s")' %
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
        """Update a chIP entry for a NDDdata object.
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
        c.execute('UPDATE dataset SET Factor = "%s", Organism = "%s", Assembly = "%s", Cell = "%s", Platform = "%s", Status = "%s", Lab = "%s", Paper = "%s", Condition = "%s", Type = "%s", Submit = "%s",Raw = "%s",XLS = "%s", Wiggle = "%s", GEO = "%s", PubMED = "%s", Comment = "%s" where id=%s' %
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

        Return a MySQLdb.Cursor object.
        """
        c = self.conn.cursor()
        c.execute(string)
        return c

    def delete (self, i):
        """Delete an entry.
        """
        try:
            self.execute ('DELETE FROM dataset WHERE id=\"%s\"' % ((str(i))))
            self.conn.commit()
        except:
            return False
        else:
            return True

    def select (self, i):
        """Select an entry by id. Return NDD object.
        
        """
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
        
