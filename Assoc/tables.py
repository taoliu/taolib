
"""Module Description

Copyright (c) 2008 H. Gene Shin <shin@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  H. Gene Shin
@contact: shin@jimmy.harvard.edu
"""

# ------------------------------------
# Python modules
# ------------------------------------
import sys,os,re,operator,copy,sqlite3,warnings,time

# ------------------------------------
# My own Python modules
# ------------------------------------
from Cistrome.Assoc.inout import *

#-------------------------------------
# classes
#-------------------------------------  
     
class Table:
    """Class Table
    
    This class is a super-class to AnnotTable, Summary, and P"""
    
    def __init__(self,name=''):
        """Constructor"""
        
        self.table={}
        self.columns=['chrom']
        self.set_name(name)
        
    def init_table(self,chrom):
        """Initialize the table of a chromosome"""
        pass
    
    def readfile(self,filepath=''):
        """Read a table from a file
        
        Parameters:
        1. filepath: the path to a file to read
    
        """
        
        #open a file
        fn=os.path.join(filepath,self.get_name()+'.xls')
        try:
            f=open(fn,'r')
        except IOError:
            raise Exception('No table, %s exists in %s' %(self.get_name(),filepath))
        
        # reset the current table
        if not self.isempty(): 
            warnings.warn("Resetting the current table, %s") %(self.get_name())
            self.table={}
         
        chrom=''
        for line in f.xreadlines():
            # obtaining the parameters
            line=line.strip()
            if line=='' or line[0]=="#": continue
                
            # get chromosomes
            elements=self.parse_line(line)
            newchrom=elements[0]
            # when new chromosome is met
            if chrom!=newchrom:
                chrom=newchrom
                self.init_table(chrom)
            self.add_row(chrom,elements[1:])
        
        f.close()

        
    def readdb(self,Db=''):
        """Read a table from a sqlite3 db file
        
        Parameters:
        1. Db: a sqlite3 db file name
        """
        
        # connect to the local database file
        dbconnect=sqlite3.connect(Db)
        cursor=dbconnect.cursor()
        tablename=self.get_name()
        
        # get the summary table
        sqlcmd="""select """+', '.join(self.columns)+""" from """ +tablename
        cursor.execute(sqlcmd)
        data=cursor.fetchall()
        
        if data: # if there is some data
            if not self.isempty(): 
                warnings.warn("Resetting the current table, %s") %(tablename)
                self.table={}
            chrom=''
            for l in data:
                newchrom=str(l[0])
                if newchrom!=chrom:
                    chrom=newchrom
                    self.init_table(chrom)
                elements=self.parse_db_line(l)
                self.add_row(chrom,elements)
 
    
    def savefile(self,filepath=''):
        """Store the table as a xls file (tab-delimited)
        
        Parameters:
        1. filepath: the path to the xls file to make
        """
        
        if self.isempty(): raise Exception('Empty table')
        #open a file
        fn=os.path.join(filepath,self.get_name()+'.xls')
        f=open(fn,'w')
        
        f.write(self.to_xls())
        
        #close the file
        f.close()
    
    def add_row(self,chrom,elements):
        """Add a row to the table"""
        pass
    
    def get_row(self,chrom):
        """Return the row of the chromosome"""
        pass
    
    def parse_line(self,line):
        """Parse a line read from a file"""
        pass
    
    def parse_db_line(self,dbline):
        """Parse a line read from a db"""
        pass
    
    def fit_to_db(self,chrom,i):
        """Fit the ith line into a db"""
        pass
    
    def to_xls(self):
        """Turn the table into a string representation"""
        pass
    
    def __getitem__(self,key):
        """Emulate __getitem__ method"""
        
        return self.table[key]
                
    def __setitem__(self,key,value):
        """Emulate __setitem__ method"""
        
        self.table[key]=value
        
    def __delitem__(self,key):
        """Emulate __delitem__ method"""

        del self.table[key]
        
    def get_column_names(self):
        """Return the columns of the table"""
        
        return self.columns      
    
    def get_chroms(self):
        
        chroms=self.table.keys()
        chroms.sort()
        return chroms
    
    def get_name(self):
        """Get the table name"""
        
        return self.__name__
    
    def set_name(self,name):
        """Set the table name"""
        
        self.__name__=name
        
    def add_column_names(self,col_names):
        """Add column names to the table"""
        
        self.columns+=col_names
        
    def isempty(self):
        """Check if the table is empty"""
        
        if not self.table: return True
        else: False
        
    def get_column_num(self):
        """Return the number of columns"""
        
        return len(self.columns)
    
    def has_chrom(self,chrom):
        """Return whether this table has the given chromosome"""
        
        return self.table.has_key(chrom)
                    
    def __str__(self):
        """Emulate __str__ method"""
        
        return self.to_xls()
    
    
class AnnotTable(Table):
    """Class AnnotTable
    
    This class is a container of annotation. This class inherits Table"""
    
    def __init__(self,name=''):
        """Constructor"""
        Table.__init__(self,name)
        self.add_column_names(['coordinate','promoter','bipromoter','downstream',\
                               'gene','rel_loc','rel_loc_cds','roi'])

        
    def init_table(self,chrom):
        """Initialize the table of a chromosome"""
        
        self[chrom]={'coordinate':[],'promoter':[], 'bipromoter':[], 'downstream':[],\
                     'gene':[],'rel_loc':[],'rel_loc_cds':[],'roi':[]}
    
    def readdb(self,Db=None,chrom=''):
        """Read an AnnotTable from a sqlite3 db file. If chrom is given, only the data of the chromosome are read.
        
        1. Db: the name of a db file to read
        2. chrom: a chromosome to read. If not given, the whole table is read.
        
        """
        # connect to the local database file
        dbconnect=sqlite3.connect(Db)
        cursor=dbconnect.cursor()
        tablename=self.get_name()
        
        # get the summary table
        if chrom:
            sqlcmd="""select """+', '.join(self.columns)+""" from """+tablename+""" where chrom='%s'""" %(chrom)
        else:  
            sqlcmd="""select """+', '.join(self.columns)+""" from """ +tablename
            
        # fetch the data 
        cursor.execute(sqlcmd)
        data=cursor.fetchall()
        
        if data:    # if data is NOT empty
            if not self.isempty(): 
                warnings.warn("Resetting the current table, %s") %(tablename)
                self.table={}
            chrom=''
            for l in data:
                newchrom=str(l[0])
                if newchrom!=chrom:
                    chrom=newchrom
                    self.init_table(chrom)
                elements=self.parse_db_line(l)
                self.add_row(chrom,elements)
        
             
    def savefile(self,filepath='',resolution=None,prom=None,biprom=None,down=None,gene_div=None):
        """Save the genome annotation table as a xls file. This method overrides the method of Table.
        
        Parameters
        1. filepath: the path to the file to open.
        2. resolution: annotation resolution.
        3. prom: promoter range
        4. biprom: bidirectional promoter range
        5. down: downstream range
        6. gene_div: the number of divisions for relative location within a gene or CDS
        """
        
        if self.isempty(): raise Exception('Empty table')
        #open a file
        fn=os.path.join(filepath,self.get_name()+'.xls')
        f=open(fn,'w')
        
        f.write(self.to_xls(resolution=resolution,prom=prom,biprom=biprom,down=down,gene_div=gene_div))
        
        #close the file
        f.close()
    
    def savedb(self,Db=None,overwrite=False):
        """Store the table in a db
    
        Parameters:
        1. Db: the sqlite3 db file to open
        2. overwrite: a boolean switch. 
            False: do not overwrite if there is a table with the same name already
            True: overwrite if there is a table with the same name already
        """
        
        # if table is empty, raise an exception
        if self.isempty(): raise Exception('Empty table')
        
        # connect to the local database file
        dbconnect=sqlite3.connect(Db)
        cursor=dbconnect.cursor()
        tablename=self.get_name()
        
        # create a table
        sqlcmd="""create table %s (""" %(tablename)+','.join(self.columns)+""")"""
        try:
            cursor.execute(sqlcmd)
        except sqlite3.OperationalError,e:
            if overwrite:
                warnings.warn("Overwritting the existing table")
                cursor.execute("""drop table %s""" %(tablename))
                cursor.execute(sqlcmd)
            else:
                raise sqlite3.OperationalError(e)
         
        # write the annotation in the table
        num_columns=len(self.columns)
        chroms=self.get_chroms()
        sqlcmd="""insert into %s (""" %(tablename)+', '.join(self.columns)+""") values (""" +', '.join(['?' for ix in range(0,num_columns)])+""")""" 
        for chrom in chroms:
            length=len(self.table[chrom][self.columns[1]])
            for i in xrange(0,length):
                export_data=[chrom]
                export_data.extend(self.fit_to_db(chrom,i))
                export_data=tuple(export_data)
                cursor.execute(sqlcmd,export_data)
                
        # save the change
        dbconnect.commit()
        #close the database file
        cursor.close()
        dbconnect.close()
        
    def append(self,Db=None):
        """Append more info to the existing db
        
        Note that the currrent version of this method does not have any way of checking if the table where to be appended 
        already exists. This method runs assuming that the table already exists and it only appends more data to the table. 
        Here, the table referes to the database table. (not self.table)
        
        """
        
         # if table is empty, raise an exception
        if self.isempty(): raise Exception('Empty table')
        
        # connect to the local database file
        dbconnect=sqlite3.connect(Db)
        cursor=dbconnect.cursor()
        tablename=self.get_name()
    
        # write the annotation in the table
        num_columns=len(self.columns)
        chroms=self.get_chroms()
        sqlcmd="""insert into %s (""" %(tablename)+', '.join(self.columns)+""") values (""" +', '.join(['?' for ix in range(0,num_columns)])+""")""" 
        for chrom in chroms:
            length=len(self.table[chrom][self.columns[1]])
            for i in xrange(0,length):
                export_data=[chrom]
                export_data.extend(self.fit_to_db(chrom,i))
                export_data=tuple(export_data)
                try:
                    cursor.execute(sqlcmd,export_data)
                except sqlite3.OperationalError,e:
                    a=str(e)
                    if a.startswith('no such table:'):
                        sqlcmd2="""create table %s (""" %(tablename)+','.join(self.columns)+""")"""
                        cursor.execute(sqlcmd2)
                        cursor.execute(sqlcmd,export_data)
                    else:
                        raise sqlite3.OperationalError(e)
                        
        # save the change
        dbconnect.commit()
        #close the database file
        cursor.close()
        dbconnect.close()
               
    def parse_line(self,line):
        """Parse a line from a file"""
        
        l=line.strip().split()
        elements=[l[0]]+[int(i) for i in l[1:]]
        
        return elements
    
    def parse_db_line(self,l):
        """Parse a line read from a db"""
        
        elements=l[1:6]+tuple([int(p) for p in l[6].rstrip(',').split(',')]+[int(p) for p in l[7].rstrip(',').split(',')]+[l[-1]])
        
        return elements
    
    def fit_to_db(self,chrom,i):
        """Fit the ith row of the table to the db"""
        
        export_data=[]
        for column in self.columns[1:]:
            element=self.table[chrom][column][i]
            if type(element).__name__=='list': element=str(element)[1:-1]+','
            try:
                export_data.append(abs(element))
            except TypeError:
                export_data.append(element)
                
        return export_data
        
    def add_row(self,chrom,elements):
        """Add a row to the table"""
        
        try:
            self[chrom]['coordinate'].append(elements[0])
            self[chrom]['promoter'].append(elements[1])
            self[chrom]['bipromoter'].append(elements[2])
            self[chrom]['downstream'].append(elements[3])
            self[chrom]['gene'].append(elements[4])
            self[chrom]['rel_loc'].append([elements[5],elements[6]])
            self[chrom]['rel_loc_cds'].append([elements[7],elements[8]])
            self[chrom]['roi'].append(elements[-1])
        except IndexError:
            raise Exception("Invalid table format: %s" %('\t'.join([chrom]+[str(e) for e in elements])))  
    
    def to_xls(self,resolution=None,prom=None,biprom=None,down=None,gene_div=None):
        """Turn the genome annotation table into a string representation
        
        Parameters
        1. filepath: the path to the file to open.
        2. resolution: annotation resolution.
        3. prom: promoter range
        4. biprom: bidirectional promoter range
        5. down: downstream range
        6. gene_div: the number of divisions for relative location within a gene or CDS 
        """
        
        columns=self.columns
        chroms=self.get_chroms()
        # put headers
        xlsout="# This file was generated by CEAS at %s.\n" %time.strftime("%H:%M:%S %a, %d %b %Y", time.localtime())
        xlsout+="# name: %s\n" %(self.get_name())   
        if resolution: xlsout+="# resolution: %d\n" %(resolution)
        else:xlsout+="# resolution: \n"
        if prom: xlsout+="# promoter: %d bp\n" %(prom)
        else: xlsout+="# promoter: \n"
        if biprom: xlsout+="# bipromoter: %d bp\n" %(biprom)
        else: xlsout+="# bipromoter: \n"
        if down: xlsout+="# downstream: %d bp\n" %(down)
        else: xlsout+="# downstream: \n"
        xlsout+="# gene: 3'UTR,5'UTR,exon,intron\n"
        if gene_div: xlsout+="# relative location within gene or CDS: %d,%d divisions\n" %gene_div
        else: xlsout+="# relative location within gene or CDS: \n"
        
        xlsout+="# chr\tcoordinate\tpromoter\tbipromoter\tdownstream\tgene\trel loc gene\trel loc gene CDS\tregion of interest\n"
        format_str="%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        for chrom in chroms:
            length=len(self[chrom][self.columns[1]])
            for i in xrange(0,length):
                export_data=[chrom]
                for column in columns[1:]:
                    element=self[chrom][column][i]
                    if type(element).__name__=='list': export_data.extend(element)
                    else: export_data.append(abs(element))
                export_data=tuple(export_data)
                xlsout+= format_str %export_data
        
        return xlsout.strip()
    
    def size(self,chrom=''):
        """Return a tuple of the table size"""
        
        numcol=len(self.columns)
        numrow=0
        if not chrom:
            chrs=self.get_chroms()
            try:
                for chr in chrs:
                    numrow+=len(self.table[chr][self.columns[1]])
            except IndexError:
               numrow=0
        else:
            try:
                numrow=len(self.table[chrom][self.columns[1]])
            except IndexError:
                pass
        
        return (numrow,numcol)
    
         
class Summary(Table):
    """Class Summary
    
    This class is a container of an annotation summary. This class inherits Table"""
    
    def __init__(self,name=''):
        """Constructor"""
        Table.__init__(self,name)
        self.add_column_names(['promoter','bipromoter','downstream','gene','rel_loc','rel_loc_cds','roi','Ns'])
        
    def init_table(self,chrom):
        """Return the table of a chromosome
        """
        
        self[chrom]={'promoter':[0,0,0],'bipromoter':[0,0],'downstream':[0,0,0],'gene':[0,0,0,0,0],\
                     'rel_loc':[[0,0,0],[0,0,0,0,0]],'rel_loc_cds':[[0,0,0],[0,0,0,0,0]],'roi':0,'Ns':0}
        
        
    def parse_line(self,line):
        """Parse a line from a file"""
        
        l=line.strip().split()
        elements=[l[0]]+[int(i) for i in l[1:]]
        
        return elements
    
    def parse_db_line(self,l):
        """Parse a line read from a db"""
        
        elements=[]
        for e in l[1:]:
            type_e=type(e)
            if type_e==str or type_e==unicode:
                elements.append([int(i) for i in e.rstrip(',').split(',')])
            elif type_e==int:
                elements.append(e)
                
        # re-parse rel_loc and rel_loc_cds
        i=self.columns.index('rel_loc')-1
        elements[i]=[elements[i][:3],elements[i][3:]]
        i=self.columns.index('rel_loc_cds')-1
        elements[i]=[elements[i][:3],elements[i][3:]]
        
        return elements
    
    def fit_to_db(self,chrom):
        """Fit the ith row to the db"""

        export_data=[chrom]
        for column in self.columns[1:]:
            element=self.table[chrom][column]
            if column=='rel_loc' or column=='rel_loc_cds': element=element[0]+element[1]
            if type(element)==list: element=str(element)[1:-1]+','
            try:
                export_data.append(abs(element))
            except TypeError:
                export_data.append(element)
                
        return export_data
    
    def savedb(self,Db=None,overwrite=False):
        """Store the table in a db
        
        Parameters:
        1. Db: the sqlite3 db file to open
        2. overwrite: a boolean switch. 
            False: do not overwrite if there is a table with the same name already
            True: overwrite if there is a table with the same name already
        """
        
        # if table is empty, raise an exception
        if self.isempty(): raise Exception('Empty table')
        
        # connect to the local database file
        dbconnect=sqlite3.connect(Db)
        cursor=dbconnect.cursor()
        tablename=self.get_name()
        
        # create a table
        sqlcmd="""create table %s (""" %(tablename)+','.join(self.columns)+""")"""
        try:
            cursor.execute(sqlcmd)
        except sqlite3.OperationalError,e:
            if overwrite:
                warnings.warn("Overwritting the existing table")
                cursor.execute("""drop table %s""" %(tablename))
                cursor.execute(sqlcmd)
            else:
                raise sqlite3.OperationalError(e)
            
        # write the annotation in the table
        num_columns=len(self.columns)
        chroms=self.get_chroms()
        sqlcmd="""insert into %s (""" %(tablename)+', '.join(self.columns)+""") values (""" +', '.join(['?' for ix in range(0,num_columns)])+""")""" 
        for chrom in chroms:
            export_data=self.fit_to_db(chrom)
            cursor.execute(sqlcmd,tuple(export_data))
                
        # save the change
        dbconnect.commit()
        #close the database file
        cursor.close()
        dbconnect.close()
        
    def add_row(self,chrom,elements):
        """Add a row to the table"""
        
        self.init_table(chrom)
        
        for column,el in itertools.izip(self.columns[1:],elements):
            self[chrom][column]=el
          
    def get_row(self,chrom):
        """Return the row of a chromosome"""
        
        elements=[]
        for column in self.columns[1:]:
            elements.append(self[chrom][column])
            
        return elements

        
    def summarize(self):
        """Summarize the table in 'whole' row"""
        
        try:
            chroms=self.get_chroms()
            if not chroms: raise ValueError
        except (AttributeError,ValueError):
            return
                    
        # obtain a summary statistics
        self.init_table('whole')
        array_adder=lambda x,y: [x[i]+y[i] for i in range(0,len(x))]
        for chrom in chroms:
            if chrom!='whole':
                self['whole']['promoter']=array_adder(self['whole']['promoter'],self[chrom]['promoter'])
                self['whole']['bipromoter']=array_adder(self['whole']['bipromoter'],self[chrom]['bipromoter'])
                self['whole']['downstream']=array_adder(self['whole']['downstream'],self[chrom]['downstream'])
                self['whole']['gene']=array_adder(self['whole']['gene'],self[chrom]['gene'])
                self['whole']['rel_loc'][0]=array_adder(self['whole']['rel_loc'][0],self[chrom]['rel_loc'][0])
                self['whole']['rel_loc'][1]=array_adder(self['whole']['rel_loc'][1],self[chrom]['rel_loc'][1])
                self['whole']['rel_loc_cds'][0]=array_adder(self['whole']['rel_loc_cds'][0],self[chrom]['rel_loc_cds'][0])
                self['whole']['rel_loc_cds'][1]=array_adder(self['whole']['rel_loc_cds'][1],self[chrom]['rel_loc_cds'][1])
                self['whole']['roi']+=self[chrom]['roi']
                self['whole']['Ns']+=self[chrom]['Ns']
                
    def get_p(self):
        """Return a P object of probabilites"""
        
        p=P()
        try:
            chroms=self.get_chroms()
            if not chroms: raise ValueError
        except (AttributeError,ValueError):
            return p
        
        for chrom in chroms:
            total=self[chrom]['Ns']
            p.init_table(chrom)
            # check if the denominator is zero
            try:
                p[chrom]['promoter']=map(lambda x: 1.0*x/total,self[chrom]['promoter'])
            except ZeroDivisionError:
                total=1
                p[chrom]['promoter']=map(lambda x: 1.0*x/total,self[chrom]['promoter'])
            p[chrom]['bipromoter']=map(lambda x: 1.0*x/total,self[chrom]['bipromoter'])
            p[chrom]['downstream']=map(lambda x: 1.0*x/total,self[chrom]['downstream'])
            p[chrom]['gene']=map(lambda x: 1.0*x/total,self[chrom]['gene'])
            
            #relative locations
            total_rel_loc=sum(self[chrom]['rel_loc'][0])
            total_rel_loc_cds=sum(self[chrom]['rel_loc_cds'][0])
        # check if the denominator is zero
            try:
                p_rel_loc0=map(lambda x: 1.0*x/total_rel_loc,self[chrom]['rel_loc'][0])
            except ZeroDivisionError:
                total_rel_loc=1
                p_rel_loc0=map(lambda x: 1.0*x/total_rel_loc,self[chrom]['rel_loc'][0])
            p_rel_loc1=map(lambda x: 1.0*x/total_rel_loc,self[chrom]['rel_loc'][1])
            
            try:
                p_rel_loc_cds0=map(lambda x: 1.0*x/total_rel_loc_cds,self[chrom]['rel_loc_cds'][0])
            except ZeroDivisionError:
                total_rel_loc_cds=1
                p_rel_loc_cds0=map(lambda x: 1.0*x/total_rel_loc_cds,self[chrom]['rel_loc_cds'][0])
            p_rel_loc_cds1=map(lambda x: 1.0*x/total_rel_loc_cds,self[chrom]['rel_loc_cds'][1])
           
            p[chrom]['rel_loc']=[p_rel_loc0,p_rel_loc1]
            p[chrom]['rel_loc_cds']=[p_rel_loc_cds0,p_rel_loc_cds1]
        
        
            try:
                p[chrom]['roi']=1.0*self[chrom]['roi']/total
            except ZeroDivisionError:
                p[chrom]['roi']=1.0*self[chrom]['roi']
            except KeyError:
                pass
        
            
            try:
                p[chrom]['chroms']=1.0*self[chrom]['Ns']/self['whole']['Ns']
            except ZeroDivisionError:
                p[chrom]['chroms']=total
                
        return p

            
    
    def to_xls(self):
        """Turn the summary table into a xls representation"""
        
        columns=self.columns
        chroms=self.get_chroms()
        
        xlsout="# This file was generated by CEAS at %s.\n" %time.strftime("%H:%M:%S %a, %d %b %Y", time.localtime())
        xlsout+="# name: %s\n" %(self.get_name())
        xlsout+="# "+"\t".join(self.columns)+"\n"
        for chrom in chroms:
            elements=[chrom]
            for column in columns[1:5]: elements+=self[chrom][column]
            elements+=self[chrom]['rel_loc'][0]+self[chrom]['rel_loc'][1]+self[chrom]['rel_loc_cds'][0]+self[chrom]['rel_loc_cds'][1]+[self[chrom]['roi']]+[self[chrom]['Ns']]
            fm="%s\t" +"\t".join(["%d"]*(len(elements)-1))+"\n"
            xlsout+=fm %tuple(elements)
        
        return xlsout.strip()
      
class P(Table):
    """Class P
    
    This class is a container of annotation probabilities. This class inherits Table"""
    
    def __init__(self,name=''):
        """Constructor"""
        Table.__init__(self,name)
        self.add_column_names(['promoter','bipromoter','downstream','gene','rel_loc','rel_loc_cds','roi','chroms'])
        
    def init_table(self,chrom):
        """Return the table of a chromosome
        
        Parameters:
        1. chrom: chromosome to create
        2. bnprom: the number of bins for promoter
        3. bnbiprom: the number of bins for bipromoter
        4. bndown; the number of bins for downstream
        
        """
        
        self[chrom]={'promoter':[0.0,0.0,0.0],'bipromoter':[0.0,0.0],'downstream':[0.0,0.0,0.0],'gene':[0.0,0.0,0.0,0.0,0.0],\
                     'rel_loc':[[0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0]],'rel_loc_cds':[[0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0]],'roi':0.0,'chroms':0.0}
                    
    def parse_line(self,line):
        """Parse a line from a file"""
        
        l=line.strip().split()
        elements=[l[0]]+[float(i) for i in l[1:]]
        
        return elements
    
    def parse_db_line(self,l):
        """Parse a line read from a db"""
        
        elements=[]
        for e in l[1:]:
            type_e=type(e)
            if type_e==str or type_e==unicode:
                elements.append([float(i) for i in e.rstrip(',').split(',')])
            elif type_e==float:
                elements.append(e)
        
        # re-parse rel_loc and rel_loc_cds
        i=self.columns.index('rel_loc')-1
        elements[i]=[elements[i][:3],elements[i][3:]]
        i=self.columns.index('rel_loc_cds')-1
        elements[i]=[elements[i][:3],elements[i][3:]]
        
        return elements
    
    def fit_to_db(self,chrom):
        """Fit the ith row to the db"""

        export_data=[chrom]
        for column in self.columns[1:]:
            element=self.table[chrom][column]
            if column=='rel_loc' or column=='rel_loc_cds': element=element[0]+element[1]
            if type(element)==list: element=str(element)[1:-1]+','
            try:
                export_data.append(abs(element))
            except TypeError:
                export_data.append(element)
                
        return export_data
    
    def savedb(self,Db=None,overwrite=False):
        """Store the table in a db
        
        Parameters:
        1. Db: the sqlite3 db file to open
        2. overwrite: a boolean switch. 
            False: do not overwrite if there is a table with the same name already
            True: overwrite if there is a table with the same name already
        """
        
        # if table is empty, raise an exception
        if self.isempty(): raise Exception('Empty table')
        
        # connect to the local database file
        dbconnect=sqlite3.connect(Db)
        cursor=dbconnect.cursor()
        tablename=self.get_name()
        
        # create a table
        sqlcmd="""create table %s (""" %(tablename)+','.join(self.columns)+""")"""
        try:
            cursor.execute(sqlcmd)
        except sqlite3.OperationalError,e:
            if overwrite:
                warnings.warn("Overwritting the existing table")
                cursor.execute("""drop table %s""" %(tablename))
                cursor.execute(sqlcmd)
            else:
                raise sqlite3.OperationalError(e)
            
        # write the annotation in the table
        num_columns=len(self.columns)
        chroms=self.get_chroms()
        sqlcmd="""insert into %s (""" %(tablename)+', '.join(self.columns)+""") values (""" +', '.join(['?' for ix in range(0,num_columns)])+""")""" 
        for chrom in chroms:
            export_data=self.fit_to_db(chrom)
            cursor.execute(sqlcmd,tuple(export_data))
                
        # save the change
        dbconnect.commit()
        #close the database file
        cursor.close()
        dbconnect.close()
         
    def add_row(self,chrom,elements):
        """Add a row to the table"""
        
        for column,el in itertools.izip(self.columns[1:],elements):
            self[chrom][column]=el
            
  
        
    def get_row(self,chrom):
        """Return the row of a chromosome"""
        
        elements=[]
        
        for column in self.columns[1:]:
            elements.append(self[chrom][column])
        
        return elements
            
    
    def to_xls(self):
        """Turn the p table into a xls representation"""
        
        columns=self.columns
        chroms=self.get_chroms()
        
        xlsout="# This file was generated by CEAS at %s.\n" %time.strftime("%H:%M:%S %a, %d %b %Y", time.localtime())
        xlsout+="# name: %s\n" %(self.get_name())
        xlsout+="# "+"\t".join(self.columns)+"\n"
        for chrom in chroms:
            elements=[chrom]
            for column in columns[1:5]: elements+=self[chrom][column]
            elements+=self[chrom]['rel_loc'][0]+self[chrom]['rel_loc'][1]+self[chrom]['rel_loc_cds'][0]+self[chrom]['rel_loc_cds'][1]+[self[chrom]['roi']]+[self[chrom]['chroms']]
            fm="%s\t" +"\t".join(["%1.5e"]*(len(elements)-1))+"\n"
            xlsout+=fm  %tuple(elements)            
        return xlsout.strip()
    

class SummaryGBG(Summary):

    def __init__(self,name='',numprom=3,numbiprom=2,numdown=3):
        """Constructor"""
        
        Summary.__init__(self,name=name)
        self.numprom=numprom
        self.numbiprom=numbiprom
        self.numdown=numdown
        
    def set_dim(self,numprom,numbiprom,numdown):
        """Set the dimenstions of promoter, bipromoter and downstream"""
        
        self.numprom=numprom
        self.numbiprom=numbiprom
        self.numdown=numdown
        
    def init_table(self,chrom):
        """Return the table of a chromosome
            
        """
        
        self[chrom]={'promoter':[0]*self.numprom,'bipromoter':[0]*self.numbiprom,'downstream':[0]*self.numdown,'gene':[0,0,0,0,0],\
                     'rel_loc':[[0,0,0],[0,0,0,0,0]],'rel_loc_cds':[[0,0,0],[0,0,0,0,0]],'roi':0,'Ns':0}


    def get_p(self):
        """Return a P object of probabilites"""
        
        p=PGBG(numprom=self.numprom,numbiprom=self.numbiprom,numdown=self.numdown)
        
        try:
            chroms=self.get_chroms()
            if not chroms: raise ValueError
        except (AttributeError,ValueError):
            return p
        
        for chrom in chroms:
            total=self[chrom]['Ns']
            p.init_table(chrom)
            # check if the denominator is zero
            try:
                p[chrom]['promoter']=map(lambda x: 1.0*x/total,self[chrom]['promoter'])
            except ZeroDivisionError:
                total=1
                p[chrom]['promoter']=map(lambda x: 1.0*x/total,self[chrom]['promoter'])
            p[chrom]['bipromoter']=map(lambda x: 1.0*x/total,self[chrom]['bipromoter'])
            p[chrom]['downstream']=map(lambda x: 1.0*x/total,self[chrom]['downstream'])
            p[chrom]['gene']=map(lambda x: 1.0*x/total,self[chrom]['gene'])
            
            #relative locations
            total_rel_loc=sum(self[chrom]['rel_loc'][0])
            total_rel_loc_cds=sum(self[chrom]['rel_loc_cds'][0])
        # check if the denominator is zero
            try:
                p_rel_loc0=map(lambda x: 1.0*x/total_rel_loc,self[chrom]['rel_loc'][0])
            except ZeroDivisionError:
                total_rel_loc=1
                p_rel_loc0=map(lambda x: 1.0*x/total_rel_loc,self[chrom]['rel_loc'][0])
            p_rel_loc1=map(lambda x: 1.0*x/total_rel_loc,self[chrom]['rel_loc'][1])
            
            try:
                p_rel_loc_cds0=map(lambda x: 1.0*x/total_rel_loc_cds,self[chrom]['rel_loc_cds'][0])
            except ZeroDivisionError:
                total_rel_loc_cds=1
                p_rel_loc_cds0=map(lambda x: 1.0*x/total_rel_loc_cds,self[chrom]['rel_loc_cds'][0])
            p_rel_loc_cds1=map(lambda x: 1.0*x/total_rel_loc_cds,self[chrom]['rel_loc_cds'][1])
           
            p[chrom]['rel_loc']=[p_rel_loc0,p_rel_loc1]
            p[chrom]['rel_loc_cds']=[p_rel_loc_cds0,p_rel_loc_cds1]
        
        
            try:
                p[chrom]['roi']=1.0*self[chrom]['roi']/total
            except ZeroDivisionError:
                p[chrom]['roi']=1.0*self[chrom]['roi']
            except KeyError:
                pass
        
            
            try:
                p[chrom]['chroms']=1.0*self[chrom]['Ns']/self['whole']['Ns']
            except ZeroDivisionError:
                p[chrom]['chroms']=total
                
        return p


class PGBG(P):
    
    def __init__(self,name='',numprom=3,numbiprom=2,numdown=3):
        """Constructor"""
        
        P.__init__(self,name=name)
        self.numprom=numprom
        self.numbiprom=numbiprom
        self.numdown=numdown
        
    def set_dim(self,numprom,numbiprom,numdown):
        """Set the dimenstions of promoter, bipromoter and downstream"""
        
        self.numprom=numprom
        self.numbiprom=numbiprom
        self.numdown=numdown
        
    def init_table(self,chrom):
        """Return the table of a chromosome
            
        """
        
        self[chrom]={'promoter':[0.0]*self.numprom,'bipromoter':[0.0]*self.numbiprom,'downstream':[0.0]*self.numdown,'gene':[0,0,0,0,0],\
                     'rel_loc':[[0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0]],'rel_loc_cds':[[0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0]],'roi':0.0,'chroms':0.0}

