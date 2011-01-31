

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
# python modules
# ------------------------------------
import sys,time,re,operator,sqlite3,warnings
from array import *
import itertools
import random
import copy
from bisect import insort,bisect_left,bisect_right

# ------------------------------------
# My own Python modules
# ------------------------------------
import Cistrome.Assoc.R as R
import Cistrome.Assoc.graphics as graphics
from Cistrome.Assoc.corelib import *


#-------------------------------------
# exception and warning class
#------------------------------------- 
class CEASError(Exception):
    """CEASError class"""
    def __init__(self,message='',value=None):
        self.message=message
        self.value=value
        
    def __str__(self):
        return self.message
    
class CEASWarning(UserWarning):
    """CEASWarning class"""
    
    def __init__(self,message=''):
        self.message=message
        
    def __str__(self):
        return self.message
        
# ------------------------------------
# necessary python packages
# ------------------------------------
MYSQL=True
try:
    import MySQLdb
except ImportError:
    MYSQL=False
    warnings.warn("sqlite3 is used instead of MySQLdb because MySQLdb is not installed")

# ------------------------------------
# some literals
# ------------------------------------
standard_chroms={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
# to determine the byte size
if array('H',[1]).itemsize == 2:
    BYTE2 = 'H'
else:
    raise Exception("BYTE2 type cannot be determined!")

if array('I',[1]).itemsize == 4:
    BYTE4 = 'I'
elif array('L',[1]).itemsize == 4:
    BYTE4 = 'L'
else:
    raise Exception("BYTE4 type cannot be determined!")

if array('f',[1]).itemsize == 4:
    FBYTE4 = 'f'
elif array('d',[1]).itemsize == 4:
    FBYTE4 = 'd'
else:
    raise Exception("BYTE4 type cannot be determined!")

#-------------------------------------
# classes
#-------------------------------------  
class Bed(object):
    """class Bed
    
    This class reads a bed file and stores."""
    
    def __init__(self):
        """Constructor"""
        
        #------ public
        self.bed=None
        self.__name__=None
        
        #------ private
        self.__columns=['chr','start','end','name','score','strand','thickStart','thickEnd',\
                        'itemRgb','blockCount','blockSizes','blockStarts']
        self.__types=['str','int','int','str','float','str','int','int','list','int','list','list']
        self.__format_cols=['%s','%d','%d','%s','%f','%s','%d','%d','%s','%d','%s','%s']
        self.__list_i=[i for i in range(0,len(self.__types)) if self.__types[i]=='list']
        self.__chroms=None
        
        
    def read(self,fn=''):
        """Read a bed file
        
        Parameters:
        1. fn: file name
        """
        
        #open the bed file
        f=open(fn,'r')
        self.__name__=fn
        
        standard_chroms={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
        types=self.__types
        type_coversion={'int':'l','float':'d'}
        self.bed={}
        chrom=''
        for line in f.xreadlines():
            if re.search(r'track',line) or line=='' or line=='\n' or line[0]=='#':continue
            l=line.strip().split()
            
            # when new chromosome met, initialize the dictionary 
#            try:
#                newchrom=standard_chroms[l[0]]
#            except KeyError:
#                newchrom=l[0]
#            if chrom!=newchrom:

            try:
                chrom=standard_chroms[l[0]]
            except KeyError:
                chrom=l[0]
            if not self.bed.has_key(chrom):
                self.bed[chrom]={}
                for i in range(1,len(l)):
                    if types[i]=='int':        # when column is integer
                        self.bed[chrom][self.__columns[i]]=array('l',[])
                    elif types[i]=='float':    # when column is float
                        self.bed[chrom][self.__columns[i]]=array('d',[])
                    elif types[i]=='str'or types[i]=='list':       # when string or blob type (comma separated string)
                        self.bed[chrom][self.__columns[i]]=[]
                        
            for i in range(1,len(l)):
                if types[i]=='int':
                    element=int(l[i])
                elif types[i]=='float':
                    element=float(l[i])
                elif types[i]=='str':
                    element=l[i]
                elif types[i]=='list':
                    element=[int(e) for e in l[i].rstrip(',').split(',')]
                self.bed[chrom][self.__columns[i]].append(element)
        
        f.close()

    
    def write(self,fn=None):
        """Write the object to a bed file
        
        Parameters:
        1. fn: file name
        """
        
        # open the bed file to write
        f=open(fn,'w')        
        f.write(self.__to_bed())
        f.close()  
        
    def sort(self,key='start'):
        """Sort the bed according to the key.
        
        Parameters:
        1. key: the key according to which the bed is sorted (default='start')"""

        # iterate through the choromosomes
        chroms=self.get_chroms()
        
        for chrom in chroms:
            try:
                key_column=self.bed[chrom][key]
            except KeyError, e:
                raise Exception('%s is not a valid column name' %(e))
            
            # if key_column is not list, TypeError is raised. Then, convert key_column into a list and do sorting again.
            try:
                sort_order=self.__argsort(key_column)
            except TypeError:
                key_column=list(key_column)
                sort_order=self.__argsort(key_column)
            
            # sort the other columns according to sort_order
            try:
                columns=self.bed[chrom].keys()
            except AttributeError:
                raise Exception('Empty bed in %s' %(chrom))
            
            for column in columns:
                if column=='chrom': continue
                if column=='start' or column=='end' or column=='thickStart' or column=='thickEnd' or column=='blockCount':
                    self.bed[chrom][column]=array('l',[self.bed[chrom][column][so] for so in sort_order])
                elif column=='score':
                    self.bed[chrom][column]=array('d',[self.bed[chrom][column][so] for so in sort_order])
                else:
                    self.bed[chrom][column]=[self.bed[chrom][column][so] for so in sort_order]
    
            
    def __argsort(self,L):
        """house-made argsort function not to import NumPy.
        
        This method is equivalent to 'argsort' function in corelib.py
        """
    
        index_element=[(ix,L[ix]) for ix in range(0,len(L))]
        index_element.sort(key=operator.itemgetter(1))
        return [ix for ix,el in index_element]
        
        
    def __getitem__(self,chrom):
        """Emulate the built-in __getitem__ method"""
                
        return self.bed[chrom]
    
    def keys(self):
        """Emulate the built-in keys() method"""
        
        return self.bed.keys()
    
    def has_key(self,key):
        """Emulate the built-in has_key() method"""
        
        return key in self.bed.keys()
    
    def get_chroms(self):
        """Return the chromosomes"""
        try:
            chroms=self.bed.keys()
        except AttributeError:
            raise Exception('Empty bed')
        
        return sort_chroms(chroms)
    
    def __to_bed(self):
        """Turn the bed object into a bed-format string"""
        
        chroms=self.get_chroms()

        bedout=''
        for chrom in chroms:
            columns=self.bed[chrom].keys()
            len_columns=len(columns)
            formated='\t'.join(self.__format_cols[:len_columns+1])
            for i in xrange(0,len(self.bed[chrom][self.__columns[1]])):
                out=[chrom]+[self.bed[chrom][column][i] for column in self.__columns[1:len_columns+1]]
                try:
                    for j in self.__list_i: out[j]=str(out[j])[1:-1]+','
                except IndexError:
                    pass
                bedout+= formated %tuple(out)+'\n'
        
        return bedout
    
    def __str__(self):
        """Emulate __str__() method"""
        
        return self.__to_bed()
    
    def add_line(self,chrom,l):
        """Add a line"""
        
        try:
            for i in range(0,len(l)):
                if self.__types[i+1]=='int':
                    element=int(l[i])
                elif self.__types[i+1]=='float':
                    element=float(l[i])
                elif self.__types[i+1]=='str':
                    element=l[i]
                elif self.__types[i+1]=='list':
                    element=[int(e) for e in l[i].rstrip(',').split(',')]
                self.bed[chrom][self.__columns[i+1]].append(element)
        except TypeError:    # the first shot
            self.bed={}
            self.bed[chrom]={}
            for i in range(0,len(l)):
                if self.__types[i+1]=='int':        # when column is integer
                    self.bed[chrom][self.__columns[i+1]]=array('l',[int(l[i])])
                elif self.__types[i+1]=='float':    # when column is float
                    self.bed[chrom][self.__columns[i+1]]=array('d',[float(l[i])])
                elif self.__types[i+1]=='str'or self.__types[i+1]=='list':       # when string or blob type (comma separated string)
                    self.bed[chrom][self.__columns[i+1]]=[l[i]]
                elif self.__types[i+1]=='list':
                    self.bed[chrom][self.__columns[i+1]]=[[int(e) for e in l[i].rstrip(',').split(',')]]
        except KeyError:     # if new chromosome is met
            self.bed[chrom]={}
            for i in range(0,len(l)):
                if self.__types[i+1]=='int':        # when column is integer
                    self.bed[chrom][self.__columns[i+1]]=array('l',[int(l[i])])
                elif self.__types[i+1]=='float':    # when column is float
                    self.bed[chrom][self.__columns[i+1]]=array('d',[float(l[i])])
                elif self.__types[i+1]=='str'or self.__types[i+1]=='list':       # when string or blob type (comma separated string)
                    self.bed[chrom][self.__columns[i+1]]=[l[i]]
                elif self.__types[i+1]=='list':
                    self.bed[chrom][self.__columns[i+1]]=[[int(e) for e in l[i].rstrip(',').split(',')]]
                    
    def filter_by_score(self,lowercutoff,uppercutoff=None):
        """Filter elements by their scores.
        
        Parameters:
        1. lowercutoff: lower cut-off for the scores
        2. uppercutoff: upper cut-off for the socres. If this parameter is not given, only lowercutoff will be applied."""
        
        # get the look-up table to know the types of columns
        lookup={}
        for i,x in enumerate(self.__columns): lookup[x]=i
        
        # the new filtered Bed object
        filtered=Bed()
        chroms=self.get_chroms()
        ix=[]
        for chrom in chroms:
            for i in xrange(len(self.bed[chrom]['score'])):
                if uppercutoff:
                    if lowercutoff <=self.bed[chrom]['score'][i] >= uppercutoff: ix.append(i)
                else:
                    if lowercutoff <=self.bed[chrom]['score'][i]: ix.append(i)
            
            filtered.bed[chrom]={}
            columns=self.bed[chrom].keys()
            for column in columns:
                if column=='start' or column=='end' or column=='thickStart' or column=='thickEnd' or column=='blockCount':
                    filtered.bed[chrom][column]=array('l',[self.bed[chrom][column][j] for j in ix])
                elif column=='score':
                    filtered.bed[chrom][column]=array('d',[self.bed[chrom][column][j] for j in ix])
                else:
                    filtered.bed[chrom][column]=[self.bed[chrom][column][j] for j in ix]
                    
        return filtered
    
    
    def size(self, chrom=''):
        """Return how many genomic regions in this object
        
        Parameters:
        1. chrom: chromosome whose size will be obtained. If '', return a dictionary of all the chromosomes
        """
        
        # a specific chromosome
        if chrom:
            try:
                n = len(self[chrom]['start'])
            except KeyError:
                n = 0
        else:
            chroms = self.get_chroms()
            n = {}
            a = 0
            for chrom in chroms:
                try:
                    n[chrom] = len(self[chrom]['start'])
                except KeyError:
                    n[chrom] = 0
                a += n[chrom]
                
            n['whole']= a 

        return n
                    
          
class GeneTable(object):
    """Class GeneTable
    
    This class reads and store a gene annotation table from UCSC or a local sqlite3 db file
    """
     
    def __init__(self,name=''):
        """Constructor"""
        
        #------- public attributes
        self.table=None
        self.__name__=name
        
        #------- private attributes
        self.columns=None
        
    def reset(self):
        """Reset the GeneTable"""
        
        self.table = {}
        self.__name__=''
        self.columns = ()
        
        
    def read(self, Host=None, User=None, Db=None, annotation=None, columns=('name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts', 'exonEnds'), which='', where=()): 
        """read a gene annotation table from a database (UCSC or local sqlite3)
        
        Parameters:
        When ucsc genome database is used:
        1. Host: "genome-mysql.cse.ucsc.edu"
        2. User: "genome"
        3. Db: organism name (e.g., "ce4" (C elegans))
        4. annotation: the name of the table to read ("refGene" (refSeq genes))
        
        When a local sqlite3 db file is used:
        1. Db: the local sqlite3 db file name
        2. annotation: the name of the table to read (e.g., "refGene") 
        
        Other parameters:
        1. columns: the table column names to read
        2. which: the column name that is used to read certain specific genes from the table
        3. where: a list of values for the column specified by 'which'.
        """     
        
        # connect and read the database table
        if Host and User:
            try:
                dbconnect=MySQLdb.connect(host=Host,user=User)
                cursor=dbconnect.cursor()
                cursor.execute('use '+Db)
            except NameError:
                dbconnect=sqlite3.connect(Db)
                cursor.execute()
        else:
            dbconnect=sqlite3.connect(Db)
            cursor=dbconnect.cursor()
        
        if which and (which not in columns): columns=columns+(which,)
        
        # store the columns to read as an attribute
        self.columns=columns
        
        # get the schema and their indicies
        schema={}
        for i,x in enumerate(columns): schema[x]=i
        if which and where:
            num_where=len(where)
            if num_where>500:
                table=[]
                sqlcmd="SELECT "+ ", ".join(columns)+" FROM "+annotation+" WHERE "
                for j in range(0,num_where/500+1):
                    sqlcmd="SELECT "+ ", ".join(columns)+" FROM "+annotation+" WHERE "
                    it=iter(where[j*500:(j+1)*500])
                    for i in it:
                        if type(i).__name__=='str': 
                            sqlcmd+= which + "='" + i+"'" + " or "
                        elif type(i).__name__=='list' or type(i).__name__=='tuple' or type(i).__name__=='array':
                            sqlcmd+=which + "=" + str(i)[1:-1]+ " or " 
                        else: sqlcmd+=which + "=" + str(i)+ " or " 
                    sqlcmd=sqlcmd[:-4]
                    cursor.execute(sqlcmd)
                    table+=cursor.fetchall()
            else:
                sqlcmd="SELECT "+ ", ".join(columns)+" FROM "+annotation+" WHERE "
                it=iter(where)
                for i in it:
                    if type(i).__name__=='str': 
                        sqlcmd+= which + "='" + i+"'" + " or "
                    elif type(i).__name__=='list' or type(i).__name__=='tuple' or type(i).__name__=='array':
                        sqlcmd+=which + "=" + str(i)[1:-1]+ " or " 
                    else: sqlcmd+=which + "=" + str(i)+ " or " 
                sqlcmd=sqlcmd[:-4]
                cursor.execute(sqlcmd)
                table=cursor.fetchall()
        else:  
            sqlcmd="SELECT "+ ", ".join(columns)+" FROM "+annotation
            cursor.execute(sqlcmd)
            table=cursor.fetchall()
    
        # parse the database table.
        # lambda expression for parsing exonStarts and exonEnds
        parser_commasep_array=lambda strarray: array('l',[int(k) for k in strarray.tostring().rstrip(',').split(',')])
        parser_commasep_str=lambda strlist: array('l',[int(k) for k in strlist])
        
        # get where 'chrom' is. 'chrom' must be included
        try:
            i=schema['chrom']
        except KeyError:
            raise Exception("'chrom' must be included in 'columns'")
        
        # extract the 'chrom' column and get the unique chromosome names
        chroms=set(map(operator.itemgetter(i),table))
        chroms=list(chroms)
        chroms.sort()
        
        # set the table up with the chromosome names
        self.table={}
        for chrom in chroms: self.table[str(chrom)]={}
        
        # run through the table
        for tb in table:
            for column in columns:
                if column=='chrom': continue
                value=tb[schema[column]]
                ctype=type(tb[schema[column]]).__name__
                
                # if blob (array in the case of MySQLdb, a comma separate string in the case of sqlite3), parse them
                if ctype=='array':
                    value=parser_commasep_array(value)
                elif ctype=='str' or ctype=='unicode': 
                    try:                    # try block for checking if the string (value) is empty
                        if value[-1]==',':    # this is comma-separated blob
                            try:
                                temp=value.rstrip(',').split(',')
                                value=parser_commasep_str(temp)
                            except ValueError:
                                pass
                        else: value=str(value)
                    except IndexError:     
                        pass
                # when no array or list was formed, make. Otherwise, just append the new data to the existing array or list
                try:
                    self.table[tb[schema['chrom']]][column].append(value)
                except KeyError:
                    if ctype=='int' or ctype=='long': self.table[tb[schema['chrom']]][column]=array('l',[value])
                    elif ctype=='float': self.table[tb[schema['chrom']]][column]=array('d',[value])
                    else: self.table[tb[schema['chrom']]][column]=[value]
                
        #close cursor and connection
        cursor.close()
        dbconnect.close()
        
    def __getitem__(self,chrom):
        """Emulate the built-in __getitem__() method"""
        
        return self.table[chrom]
    
    def __setitem__(self,chrom,value):
        
        self.table[chrom]=value

    def get_genelengths(self):
        """Return the minimum, average, and maximum gene lengths with std"""
        
        chroms=self.get_chroms()
        
        lengths={}
        total_min_gene=0
        total_max_gene=0
        total_mean_gene=0
        total_min_exon=0
        total_max_exon=0
        total_mean_exon=0
        total_min_intron=0
        total_max_intron=0
        total_mean_intron=0
        total_gene_length=0
        total_exon_length=0
        total_intron_length=0
        total_gene_num=0
        for chrom in chroms:
            gene_length=[end-start for start,end in zip(self[chrom]['txStart'],self[chrom]['txEnd'])]
            exon_length=[]
            intron_length=[]
            current_gene_num=len(self[chrom]['txStart'])
            for i in xrange(0,current_gene_num):
                exons=zip(self[chrom]['exonStarts'][i],self[chrom]['exonEnds'][i])
                try:
                    exon_length.append(sum([end-start for start,end in exons]))
                    intron_length.append(gene_length[i]-exon_length[i])
                except:
                    pass
            current_gene_length=sum(gene_length)
            current_exon_length=sum(exon_length)
            current_intron_length=current_gene_length-current_exon_length
            current_min_gene=min(gene_length)
            current_max_gene=min(gene_length)
            current_min_exon=min(exon_length)
            current_max_exon=max(exon_length)
            current_min_intron=min(intron_length)
            current_max_intron=max(intron_length)
            current_mean_gene=1.0*current_gene_length/current_gene_num
            current_mean_exon=1.0*current_exon_length/current_gene_num
            current_mean_intron=1.0*current_intron_length/current_gene_num
            
            lengths[chrom]={'gene':[current_min_gene,current_mean_gene,current_max_gene], 'exon':[current_min_exon,current_mean_exon,current_max_exon],\
                            'intron':[current_min_intron,current_mean_intron,current_max_intron]}
            
            if total_min_gene==0: total_min_gene=current_min_gene
            else: total_min_gene=min(total_min_gene,current_min_gene)
            if total_min_exon==0: total_min_exon=current_min_exon
            else: total_min_exon=min(total_min_exon,current_min_exon)
            if total_min_intron==0: total_min_intron=current_min_intron
            else: total_min_intron=min(total_min_intron,current_min_intron)
            
            total_max_gene=max(total_max_gene,current_max_gene)
            total_max_exon=max(total_max_exon,current_max_exon)
            total_max_intron=max(total_max_intron,current_max_intron)
            
            total_gene_length+=current_gene_length
            total_exon_length+=current_exon_length
            total_intron_length+=current_intron_length
            total_gene_num+=current_gene_num
            
            
        total_mean_gene=1.0*total_gene_length/total_gene_num
        total_mean_exon=1.0*total_exon_length/total_gene_num
        total_mean_intron=1.0*total_intron_length/total_gene_num
        
        lengths['whole']={'gene':[total_min_gene,total_mean_gene,total_max_gene], 'exon':[total_min_exon,total_mean_exon,total_max_exon],\
                          'intron':[total_min_intron,total_mean_intron,total_max_intron]}
        
        return lengths
    
    def remove(self,chrom):
        """Remove the gene annotation table of the given chromosome.
        
        Parameters:
        1. chrom: the chromosome to remove
        """
        
        try:
            del self.table[chrom]
        except KeyError:
            pass
            
            
    def get_gene_lens(self, chrom=''):
        """Return the gene lengths"""
        
        # get chroms
        chroms=self.get_chroms()
        
        lengths=[]
        if chrom != '' and chrom in chroms: # a specific chromomosome is given
            lengths = map(lambda start, end: end - start, self[chrom]['txStart'], self[chrom]['txEnd'])
            return lengths
        elif chrom == '':   # all the chromosomes are considered
            for chrom in chroms:
                lengths+= map(lambda start, end: end - start, self[chrom]['txStart'], self[chrom]['txEnd'])
            return lengths
        else:
            return lengths
    
    
    def turn2dict(self):
        """Return a dictionary whose keys are gene names. Values such as txStart and txEnd are accessed by key.
        
        Genes (or transcripts) with the same name are kept only once
        """
        
        # if empty table, just return an empty dictionary.
        if not self.table: return {}
        
        chroms=self.get_chroms()
        
        dict={}
        for chrom in chroms:
            for i in xrange(len(self[chrom]['name'])):
                # only if such a name does not exist
                if not dict.has_key(self[chrom]['name'][i]):
                    dict[self[chrom]['name'][i]]=[chrom]
                    for column in self.columns[2:]:
                        dict[self[chrom]['name'][i]].append(self[chrom][column][i])
        
        return dict
                    
    def savedb(self,Db=None,annotation='refGene',overwrite=False):
        """write a database file and create a table in it.
        
        If a database file already exists, just create a table"""
        
        # connect to the local database
        dbconnect=sqlite3.connect(Db)
        cursor=dbconnect.cursor()

        # create a db table
        other_columns=[column for column in self.columns if column!='chrom']
        columns=['chrom']
        columns.extend(other_columns)
        sqlcmd="""create table %s (""" %(annotation)+','.join(columns)+""")"""
        try:
            cursor.execute(sqlcmd)
        except sqlite3.OperationalError,e:
            if overwrite:
                warnings.warn("Overwritting the existing table with the current one")
                cursor.execute("""drop table %s""" %(annotation))
                cursor.execute(sqlcmd)
            else:
                raise sqlite3.OperationalError(e)
        
        # write the table into the db table
        num_columns=len(columns)
        chroms=self.table.keys()
        chroms.sort()
        for chrom in chroms:
            for i in xrange(0,len(self.table[chrom][other_columns[0]])):
                sqlcmd="""insert into %s values (""" %(annotation)+', '.join(['?' for ix in range(0,num_columns)])+""")""" 
                export_data=[chrom]    
                for column in other_columns:
                    element=self.table[chrom][column][i]
                    if type(element).__name__=='array': element=str(element.tolist())[1:-1]+','
                    export_data.append(element)
                export_data=tuple(export_data)
                cursor.execute(sqlcmd,export_data)
                
        # save the change
        dbconnect.commit()
        #close the database file
        cursor.close()
        dbconnect.close()   
    
    def droptable(self,Db=None,table=None):
        """Drop a table from the database"""
        
        dbconnect=sqlite3.connect(Db)
        cursor=dbconnect.cursor()
        
        sqlcmd="""drop table %s""" %(table)
        cursor.execute(sqlcmd)
        dbconnect.commit()
        cursor.close()
        dbconnect.close()
    
    def sort(self,key='txStart'):
        """Sort the gene table in asscending order"""

        # iterate through the choromosomes
        chroms=self.table.keys()
        for chrom in chroms:
            try:
                key_column=self.table[chrom][key]
            except KeyError, e:
                raise Exception('%s is not a valid column name' %(e))
            
            # if key_column is not list, TypeError is raised. Then, convert key_column into a list and do sorting again.
            try:
                sort_order=self.__argsort(key_column)
            except TypeError:
                key_column=list(key_column)
                sort_order=self.__argsort(key_column)
            
            # sort the other columns according to sort_order
            for column in self.columns:
                if column=='chrom': continue
                if type(self.table[chrom][column]).__name__=='array':
                    self.table[chrom][column]=array(self.table[chrom][column].typecode, [self.table[chrom][column][so] for so in sort_order])
                else:
                    self.table[chrom][column]=[self.table[chrom][column][so] for so in sort_order]
    
            
    def __argsort(self,L):
        """house-made argsort function not to import NumPy"""
    
        index_element=[(ix,L[ix]) for ix in range(0,len(L))]
        index_element.sort(key=operator.itemgetter(1))
        return [ix for ix,el in index_element]
    
    
    def get_chroms(self):
        """Return the chromsomes of the gene annotation table"""
        
        try:
            chroms=self.table.keys()
        except AttributeError:
            raise Exception('Empty gene annotation table')
        
        return sort_chroms(chroms)
    
    
    def get_gene_num(self):
        """Return the number of genes"""
        
        gene_num={}
        chroms=self.get_chroms()
        for chrom in chroms:
            gene_num[chrom]=len(self[chrom]['txStart'])
        
        gene_num['total']=sum([gene_num[chrom] for chrom in chroms])
        
        return gene_num
    
    def get_exon_intron_lens(self, chrom=''):
        """Return exon and introns lengths from a given GeneTable.
    
        Parameters:
        1. chrom: chromosome to investigate. If not given, do about all the chromosomes.
        """
    
        # if exonStarts and exonEnds are not in the GeneTable, just return empty lists
        if 'exonStarts' not in self.columns or 'exonEnds' not in self.columns:
            return [],[]
    
        if chrom =='': 
            elengths=[]
            ilengths=[]
            for chrom in self.get_chroms():
                for es,ee in itertools.izip(self[chrom]['exonStarts'], self[chrom]['exonEnds']):
                    elengths.extend(array_adder(ee,map(lambda x: -1*x, es))[0])
                    ilengths.extend(array_adder(es[1:],map(lambda x: -1*x,ee[:-1]))[0])
        else:
            # if the given chromosome is not in the GeneTable, just return empty lists.
            if chrom not in self.get_chroms(): 
                return [],[]
            elengths=[]
            ilengths=[]
            for es,ee in itertools.izip(self[chrom]['exonStarts'], self[chrom]['exonEnds']):
                elengths.extend(array_adder(ee,map(lambda x: -1*x, es))[0])
                ilengths.extend(array_adder(es[1:],map(lambda x: -1*x,ee[:-1]))[0])
            
        elengths.sort()
        ilengths.sort()
        
        return elengths, ilengths
    
    def get_cat_exon_intron_lens(self, chrom=''):
        """Return concatenated exons and introns lenghts of this GeneTable
        
        Parameters:
        1. chrom: chromosome to take count in. If not given, all of the chromsomes will be considered.
        """
        
        if 'exonStarts' not in self.columns or 'exonEnds' not in self.columns:
            return [], []
        
        if chrom == '':
            cate = []
            cati = []
            
            chroms = self.get_chroms()
            for chrom in chroms:
                 cate += map(lambda estarts, eends: sum(map(lambda s, e: e-s, estarts, eends)), self[chrom]['exonStarts'], self[chrom]['exonEnds'])
                 cati += map(lambda estarts, eends: sum(map(lambda s, e: s-e, estarts[1:], eends[:-1])), self[chrom]['exonStarts'], self[chrom]['exonEnds'])
            
            return cate, cati
        else:
            chroms = self.get_chroms()
            
            if chrom in chroms:
                cate = map(lambda estarts, eends: sum(map(lambda s, e: e-s, estarts, eends)), self[chrom]['exonStarts'], self[chrom]['exonEnds'])
                cati = map(lambda estarts, eends: sum(map(lambda s, e: s-e, estarts[1:], eends[:-1])), self[chrom]['exonStarts'], self[chrom]['exonEnds'])
            else:
                cate = []
                cati = []
            
            return cate, cati
                
            
class Wig:
    """Class Wig
    
    This class reads and stores a wig file as in a dictionary. This class only can read a variableStep wig file.
    For other types of wig files, use Tao's Cistrome
    """
    
    def __init__(self):
        """Constructor"""
        
        self.wig={}
        self.description=''
        
    def read(self,fnwig=''):
        """Read a wig file.
        
        Parameters:
        1. fnwig: wig file name
        """
        
        # for C elegans
        standard_chroms={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
        
        fwig=open(fnwig,'r')
        chrom=''
        wig={}
        for line in fwig.xreadlines():
            if re.search(r'track',line): 
                try:
                    description=re.search(r'description="([A-Za-z0-9]+)"\s',line).group(1)
                    if not self.description: self.description=description
                except AttributeError:
                    pass
                continue
            if re.search(r'chrom=([A-Za-z0-9]+)\s',line):
                chrom=re.search(r'chrom=([A-Za-z0-9]+)\s',line).group(1)
                try:
                    chrom=standard_chroms[chrom]
                except KeyError:
                    pass
                
                wig[chrom]=[array('l',[]),array('d',[])]
                continue
            # get the nearest loction and compare it with the last element of locs. If not the same, add it
            l=line.strip().split()
            wig[chrom][0].append(int(l[0]))
            wig[chrom][1].append(float(l[1]))
                
        fwig.close()
        self.wig=wig
        
    def write(self,fnwig,span=1):
        """Write the wig in a variable step wiggle format.
        
        Parameters:
        1. fnwig: wig file name
        2. span: span for each bar (default=1bp)
        """
        
        fwig=open(fnwig,'w')
        chroms=self.get_chroms()
        
        for chrom in chroms:
            fwig.write('track type=wiggle_0\nvariableStep chrom=%s span=%d\n' %(chrom,span))
            for i in xrange(0,len(self.wig[chrom][0])):
                fwig.write("%d\t%f\n" %(self.wig[chrom][0][i],self.wig[chrom][1][i]))
        
        fwig.close()
                
    def __getitem__(self,key):
        """Emulate __getitem__ method"""
        
        try:
            return self.wig[key]
        except TypeError:
            raise Exception('Empty wig')
        except KeyError:
            chroms=self.wig.keys()
            return_wig={}
            for chrom in chroms:
                if key=='coordinate':
                    return_wig[chrom]=self.wig[chrom][0]
                elif key=='value':
                    return_wig[chrom]=self.wig[chrom][1]
                else:
                    raise KeyError,key
            return return_wig
        
    def __setitem__(self,key,value):
        """Emulate __setitem__ method"""
        
        if len(value)!=2 or type(value[0])!=array or type(value[1])!=array:
            raise Exception('Value must be a list of two arrays')
        
        self.wig[key]=value 
    
    def get_chroms(self):
        """Return the chromosomes of the wig"""
        
        try:
            chroms=self.wig.keys()
        except AttributeError:
            raise Exception('Empty wig')
        
        return sort_chroms(chroms)
    
    def get_name(self):
        """Return the name of the gene table"""
        
        return self.__name__
    
    def set_name(self,name):
        """Set the gene table name"""
        
        self.__name__=name
        
    def add_line(self,chrom,l):
        """Add a line"""
        
        try:
            self.wig[chrom][0].append(int(l[0]))
            self.wig[chrom][1].append(float(l[1]))
        except TypeError:    # the first shot
            self.wig={}
            self.wig[chrom]=[array('l',[]),array('d',[])]
            self.wig[chrom][0].append(int(l[0]))
            self.wig[chrom][1].append(float(l[1]))
        except KeyError:     # if new chromosome is met
            self.wig[chrom]=[array('l',[]),array('d',[])]
            self.wig[chrom][0].append(int(l[0]))
            self.wig[chrom][1].append(float(l[1]))
            
    def length(self,chrom=''):
        """Return the length of the wig file (similar to line number)"""
        
        if chrom:    
            try:
                return len(self.wig[chrom][0])
            except KeyError:
                return 0
        
        length=0
        for chrom in self.get_chroms():
            length+=len(self.wig[chrom][0])
            
        return length
    
    def chrom_length(self,chrom):
        """Return a very simply calculated chromosome legnth of this wig file.
        
        Simple chromosome length=the last wig point - the first wig point
        
        If there are only one or even zero point in the chromosome, just returns 0"""

        try:
            if len(self.wig[chrom][0])==0 or len(self.wig[chrom][0])==1:
                return 0
            else:
                return self.wig[chrom][0][-1]-self.wig[chrom][0][0]
        except KeyError:
            return 0
        
    def __add__(self,other):
        """Emulate __add__ method"""
        
        
        selfchroms=self.get_chroms()
        otherchroms=other.get_chroms()
        chroms=list(set(selfchroms).union(otherchroms))
        chroms.sort()
        addwig=Wig()
        
        for chrom in chroms:
        
            # if one of them do not have this chromosome, just simply copy the other's
            if chrom in selfchroms and chrom not in otherchroms:
                addwig.wig[chrom]=self.wig[chrom][:]
                continue
            elif chrom not in selfchroms and chrom in otherchroms:
                addwig.wig[chrom]=other.wig[chrom][:]
                continue
            
            END1=False
            END2=False
            # p--index of f, m--index of f2
            p=0
            m=0
            # signal processing to reduce the effect of unbalanced tag distribution between + and - strands.
            selflen=len(self.wig[chrom][0])
            otherlen=len(other.wig[chrom][0])
            if selflen==0 and otherlen==0:
                END1==True
                END2==True
            
            addwig.wig[chrom]=[array('l',[]), array('d',[])]
            # signal processing loop.
            while END1!=True or END2!=True:      # keep running when the end of either strand is reached.
                if END1!=True and END2!=True:    # not the ends of P and M
                    if self.wig[chrom][0][p]==other.wig[chrom][0][m]:
                        addwig.add_line(chrom,[self.wig[chrom][0][p],self.wig[chrom][1][p]+other.wig[chrom][1][m]])
                        p+=1
                        m+=1
                    elif self.wig[chrom][0][p]<other.wig[chrom][0][m]:     
                        addwig.add_line(chrom,[self.wig[chrom][0][p],self.wig[chrom][1][p]])    
                        p+=1
                    else:
                        addwig.add_line(chrom,[other.wig[chrom][0][m],other.wig[chrom][1][m]])
                        m+=1
                elif END1!=True and END2==True:    # The end of M but not P
                    addwig.add_line(chrom,[self.wig[chrom][0][p],self.wig[chrom][1][p]])
                    p+=1
                elif END1==True and END2!=True:    # the end of P but not M
                    addwig.add_line(chrom,[other.wig[chrom][0][m],other.wig[chrom][1][m]])
                    m+=1
                
                #check the ends.
                if p==selflen: END1=True
                if m==otherlen: END2=True
                
        return addwig
    
    def max(self,chrom=''):
        """Return the max of wig
        
        Parameters:
        1. chrom: if '', return the max of the whole wig; otherwise, return the maxi of the chrom
        """
        
        if chrom:
            try: 
                return max(self[chrom][1])
            except (KeyError, IndexError):
                return None
        else:
            mx=None
            for chrom in self.get_chroms():
                if mx==None:
                    mx=max(self[chrom][1])
                else:
                    mx=max(mx, max(self[chrom][1]))
            return mx
                
        
    def min(self,chrom=''):
        """Return the min of wig
        
        Parameters:
        1. chrom: if '', return the max of the whole wig; otherwise, return the maxi of the chrom
        """
        
        if chrom:
            try: 
                return min(self[chrom][1])
            except (KeyError, IndexError):
                return None
        else:
            mn=None
            for chrom in self.get_chroms():
                if mn==None:
                    mn=min(self[chrom][1])
                else:
                    mn=min(mn, min(self[chrom][1]))
            return mn
        
    def pseudomedian(self,chrom=''):
        """Return the pseudo median. Here Pseudo median means the median of medians of chromosomes. 
        
        Parameters:
        1. chrom: if '', return the max of the whole wig; otherwise, return the maxi of the chrom
        
        Note that this method is under construction. The method using sorting might take a long time and a large amount of memory.
        """
        
        if chrom:
            meds = []
            for chr in self.get_chroms():
                meds.append(median(self[chrom][1]))
            return median(meds)
        else:
            try:
                return median(self[chrom][1])
            except KeyError:
                return None
        
    def mean(self,chrom=''):
        """Return the mean of wig
        
        Parameters:
        1. chrom: if '', return the max of the whole wig; otherwise, return the maxi of the chrom
        """
        
        if chrom:
            n=self.length(chrom=chrom)
            try:
                return 1.0*sum(self[chrom][1])/n
            except (KeyError, IndexError, ZeroDivisionError):
                return None
        else:
            cumsum=0
            cumcnt=0
            for chrom in self.get_chroms():
                try:
                    cumsum+=sum(self[chrom][1])
                    cumcnt+=self.length(chrom=chrom)
                except (KeyError, IndexError):
                    pass
            
            try: 
                return 1.0*cumsum/cumcnt
            except ZeroDivisionError:
                return None
            
    def histogram(self, breaks, chrom=''):
        """Return a histogram as a dictionary of {'bin':[], 'count':[]}
        
        Parameters:
        1. breaks: this can be a list (bins) or an integer (number of bins). If integer, it will be split the range from min to max to 'breaks' bins.
        2. chrom: if NULL string, return the histogram of the whole WIG.
        
        WARNING: this method might take too much memory if the Wig object is large.
        """
        
        if type(breaks)==list:      # if bins are given
            histogram = {'bin':breaks[:-1]}
            if chrom:
                try:
                    histogram['count']=bin(self[chrom][1], breaks)
                    return histogram
                except KeyError:
                    return None
            else:
                w = []
                for chrom in self.get_chroms():
                    w += self[chrom][1]
                histogram['count'] = bin(w, breaks)
                
                return histogram
                
        elif type(breaks)==int: # only number of bins are given
            histogram = {}
            if chrom:
                min_this_chr = self.min(chrom)
                max_this_chr = self.max(chrom)
                
                bins = linspace(min_this_chr, max_this_chr, breaks+1)
                histogra['bin'] = bins[:-1]
                try:
                    histogram['count'] = bin(self[chrom][1], bins)
                    return histogram
                except KeyError:
                    return None   
            else:
                w = []
                for chrom in self.get_chroms():
                    w += self[chrom][1]
                min_all = min(w)
                max_all = max(w)
                
                bins = linspace(min_this_chr, max_this_chr, breaks+1)
                histogram['bin'] = bins
                histogram['count'] = bin(w, bins)
                return histogram
        else:
            return None
        
        
    def get_middle_percent(self, percent):
        """Return the left and right values of the middle %
        
        Parameters:
        1. percent: the middle percent
        """
        
        w = []
        for chrom in self.get_chroms():
            w += self[chrom][1]
        w.sort()
        wlen= len(w)
        left = int(((100-percent)/200.0)*(wlen-1))
        return w[left], w[-1*left]
    
class XLS:
    """Class XLS
    
    This class is created to handle an xls file IO for CEAS.
    
    Note that the first column must be chromosome
    """
    def __init__(self,name=''):
        """Constructor"""

        self.__name__=name
        self.columns=[]
        self.types=[]
        self.xls={}
        
    def read(self,fn,columns,types):
        """Read an xls file"""
    
        self.set_columns(columns)
        self.set_types(types)
        
        for line in open(fn,'r').xreadlines():
            if line.startswith('#') or line.startswith('Chr') or line.startswith('chr') or line=='': continue
            
            ls=line.strip().split()
            
            if not self.xls.has_key(ls[0]):
                self.xls[ls[0]]={}
                for column, type, l in itertools.izip(columns[1:],types[1:],ls[1:]):
                    self.xls[ls[0]][column]=[self.convert(l,type)]
            else:
                for column, type, l in itertools.izip(columns[1:],types[1:],ls[1:]):
                    self.xls[ls[0]][column].append(self.convert(l,type))
        
    def write(self,fn=''):
        """Write as an xls file. 
        
        Note that to write to an xls file, self.columns and self.types must be set beforehand
        """
        
        if fn:
            f=open(fn,'w')
        else:
            f=open(self.get_name()+'.xls', 'w')
            
        f.write(self.toxls())
        f.close()
        
    def add(self,line):
        """Add a line"""
        
        if type(line)==str or type(line)==unicode:
            ls=line.strip().split()
        else:
            ls=line
            
        if self.xls.has_key(ls[0]):
            for column,tp,l in itertools.izip(self.columns[1:],self.types[1:],ls[1:]):
                self.xls[ls[0]][column].append(self.convert(l,tp))
        else:
            self.xls[ls[0]]={}
            for column,tp,l in itertools.izip(self.columns[1:],self.types[1:],ls[1:]):
                self.xls[ls[0]][column]=[self.convert(l,tp)]

                    
    def get_chroms(self):
        """Get the chromosomes"""
        
        chroms=self.xls.keys()
        
        return sort_chroms(chroms)
    
    def get_name(self):
        """Get the name"""
        
        return self.__name__
    
    def set_name(self,name):
        """Set the name"""

        self.__name__=name
        
    def toxls(self):
        """To a xls format"""
        
        script='# %s\n' %self.get_name()
        script+='\t'.join(self.columns)+'\n'
        for chrom in self.get_chroms():

            for i in xrange(len(self.xls[chrom][self.columns[1]])):
                line=[str(self.xls[chrom][column][i]) for column in self.columns[1:]]
                script+= chrom+ '\t' + '\t'.join(line) + '\n'

        return script
        
    def set_columns(self,columns):
        """Set the columns"""
        
        self.columns=columns

    def set_types(self,types):
        """Set the formats"""

        self.types=types
    
    def importfromdict(self,dict):
        """Import data from a dictionary
        
        Note that the existing data in self.xls will be erased and self.columns and self.types will be erased too.
        Note that the data in dict will be copied into XLS. Thus, it is necessary to consider memory usage.
        
        """
        self.xls=copy.deepcopy(dict)
        self.columns=[]
        self.types=[]
        
            
    def exporttodict(self):
        """Export the data in self.xls to a dictionary.
        
        Like self.importfromdict, a new dictionary is copied from XLS. Thus, it is necessary to consider memory usage in case that
        XLS takes a large amount of memory.
        
        """
        
        return copy.deepcopy(self.xls)
            
    def convert(self,i,type):
        """do conversion to 'int' or 'float' from 'str'"""
        if type=='int':
            return int(i)
        elif type=='float':
            return float(i)
        else:
            return i
        
    def __str__(self):
        """Emulate str method"""
        
        return self.toxls()
    
    def __getitem__(self,chrom):
        """A simple emulation of __getitem__ method"""
        
        return self.xls[chrom]
    
    
class DataFrame:
    """Class DataFrame
    
    This class works similarly to R's dataframe. Data is added in a columnwise-manner.
    This class is used by the gene-centered annotation (GeneAnnotator in annotator.py
    """
    
    def __init__(self,name=''):
        """Constructor"""
        
        self.rownames=[]
        self.colnames=[]
        self.table={}
        self.__name__=name    # table name
        
        
    def set_name(self, name):
        """Set the data frame name"""
        
        self.__name__=name
        
    def append_column(self,column,colname=''):
        """Append a column to the table"""
        
        if not self.table: self.rownames=['%d' %i for i in range(len(column))]
        
        if len(self.rownames)!=len(column):
            raise Exception('differing number of rows, %d rows in the table, %d in the new column' %(len(self.rownames),len(column)))
        
        if colname:
            self.table[colname]=column
            self.colnames.append(colname)
        else:
            ncols=len(self.table.keys())
            self.table['C%d' %ncols]=column
            self.colnames.append('C%d' %ncols)
    
    def remove_column(self,colname):
        """Remove a column from the table"""
        
        try:
            del self.table[colname]
            self.colnames=[col for col in self.colnames if col!=colname]
        except KeyError:
            pass
        
    def pop_column(self):
        """Pop the last column"""
        
        try:
            lastcolname=self.colnames[-1]
            del self.table[lastcolname]
            self.colnames.pop()
        except IndexError:
            pass
            
    def append_row(self,row,rowname=''):
        """Append a row to the end of the table"""
        
        ncols=len(self.colnames)
        lasti=0
        for i in range(min(ncols,len(row))):
            try:
                self.table[self.colnames[i]].append(row[i])
            except KeyError: # when the table is empty
                for colname in self.colnames:
                    self.table[colname]=[]
                self.table[self.colnames[i]].append(row[i])
            
            lasti=i
        
        if lasti<ncols-1:    # when the given row is shorter than the row length of the table, add 0s
            for j in range(i,ncols):
                self.table[self.colnames[j]].append(0)
        
        if rowname:
            self.rownames.append(rowname)
        else:
            self.rownames.append('%d' %len(self.rownames))
    
    def pop_row(self):
        """Remove the last row from the table"""
        
        ncols=len(self.colnames)
        for i in range(ncols):
            self.table[self.colnames[i]].pop()
        
        self.rownames.pop()
            
    def toxls(self,wrownames=False):
        """To xls file format (tab-delimited)"""
        
        ncols=len(self.colnames)
        nrows=len(self.rownames)
        
        # if table is not NULL
        if self.table:
            header='#'+'\t'.join(self.colnames)+'\n'
            rows=''
            for i in xrange(nrows):
                if wrownames:
                    rows+=self.rownames[i]+'\t'+'\t'.join([str(self.table[colname][i]) for colname in self.colnames])+'\n'
                else:
                    rows+='\t'.join([str(self.table[colname][i]) for colname in self.colnames])+'\n'
            return header+rows
        else:
            return 'NULL data frame with %d rows and %d columns' %(nrows,ncols)
                    
    def set_rownames(self,rownames):
        """Set the rownames"""
        
        if not self.table:
            self.rownames=rownames
        else:
            if len(rownames)!=len(self.table[self.colnames[0]]):
                raise Exception('Differing number of rows')
            self.rownames=rownames
                
    def get_colnames(self):
        """Return the column names"""
        
        return self.colnames
    
    def get_rownames(self):
        """Return the row names"""
        
        return self.rownames
    
    def __str__(self):
        """Emulate __str__ method"""
        
        return self.toxls()
    
    def __getitem__(self,key):
        """Emulate __getitem__ method"""
        
        try:
            return self.table[key]
        except KeyError,e:
            if type(key)==tuple:
                if len(key)==2:
                    try:
                        rowix=self.rownames.index(key[0])
                        return self.table[key[1]][rowix]
                    except:
                        raise Exception("Reference error")
                elif len(key)==1:
                    try:
                        rowix=self.rownames.index(key[0])
                    except ValueError:
                        raise Exception("No such a row name exists,%s" %key[0])
                    return [self.table[col][rowix] for col in self.colnames]
            else:
                raise KeyError,e
            
    def __setitem__(self,key,value):
        """Emulate __setitem__ method
        probably it is not that fast."""
        
        try:
            self.table[key]=value
        except KeyError,e:
            if type(key)==tuple:
                if len(key)==2:
                    try:
                        rowix=self.rownames.index(key[0])
                        self.table[key[1]][rowix]=value
                    except:
                        raise Exception("Reference error")
                elif len(key)==1:
                    try:
                        rowix=self.rownames.index(key[0])
                    except ValueError:
                        raise Exception("No such a row or column name exists, %s" %key[0])
                    
                    try:
                        for col,val in itertools.izip(self.colnames,value): self.table[col][rowix]=val
                    except TypeError:
                        raise Exception("Value must be iterable and have the same length as each row has")
                else:
                    raise KeyError, e
                    
                    
    def fill_a_row(self,rowname):
        """Fill a row"""    
        pass
    
    def size(self):
        """Return the table dimension as a tuple (rownum,colnum)"""
        
        return (len(self.rownames),len(self.colnames))
    
    def write(self, header=''):
        """Write the DataFrame as a xls file
        """
        
        fn = self.__name__ + '.xls'
        f = open(fn, 'w')
        f.write(header)
        f.write(self.toxls())
        f.close()
    
      
#-------------------------------------
# function
#-------------------------------------  
def filter_genes(gene_table,limits=[2900,3100]):
    """Leave genes whose lengthsa re within the lmit"""
    
    chroms=gene_table.get_chroms()
    new_gt=GeneTable()
    new_gt.table={}
    new_gt.columns=gene_table.columns[:]
    for chrom in chroms:
        new_gt[chrom]={}
        columns=gene_table[chrom].keys()
        for col in columns:
            new_gt[chrom][col]=[]
        for i in xrange(0,len(gene_table[chrom]['txStart'])):
            length=gene_table[chrom]['txEnd'][i]-gene_table[chrom]['txStart'][i]  
            if length>=limits[0] and length<=limits[1]:
                for col in columns:
                    new_gt[chrom][col].append(gene_table[chrom][col][i])
    
    return new_gt

def read_gene_subsets(fns):
    """Read gene names from a single txt file or multiple files"""
    
    names=[]
    for fn in fns:
        name=[]
        for line in open(fn,'r').xreadlines(): 
            try:
                name.append(line.strip().split()[0])
            except IndexError:    # sometimes null line was parsed and null array is returned, in this case, just pass
                pass
        names.append(name)
    
    return names

def get_genome_coordinates(fnwig=None,resolution=100):
    """Get genome locations from a wig file to annotate"""
    
    #for Celegans
    standard_chroms={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
    
    fwig=open(fnwig,'r')
    chrom=''
    coordinates={}
    for line in fwig.xreadlines():
        if re.search(r'track',line): continue
        if re.search(r'chrom=([A-Za-z0-9]+)\s',line):
            chrom=re.search(r'chrom=([A-Za-z0-9]+)\s',line).group(1)
            try:
                chrom=standard_chroms[chrom]
            except KeyError:
                pass
                
            coordinates[chrom]=array('l',[])
            continue
        # get the nearest loction and compare it with the last element of locs. If not the same, add it
        coordinate=(int(round(1.0*int(line.strip().split()[0])/resolution)))*resolution+1
        if not coordinates[chrom] or coordinate!=coordinates[chrom][-1]:
            coordinates[chrom].append(coordinate)
    
    fwig.close()
    return coordinates

def read_genome_coordinates(Db=None,table=None,resolution=100):
    
    standard_chroms={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
    
    dbconnect=sqlite3.connect(Db)
    cursor=dbconnect.cursor()
    
    sqlcmd="""select chrom coordinate from ?"""
    cursor.execute(sqlcmd,(table,))
    data=cursor.fetchall()
    chroms=set(map(str(operator.itemgetter(0)),data))
    chroms=list(chroms)
    temp=[]
    for chrom in chroms:
        try:
            temp.append(standard_chroms[str(chrom)])
        except KeyError:
            temp.append(str(chrom))
    chroms=temp
    chroms.sort()
    
    coordinates={}
    # set the table up with the chromosome names
    for chrom in chroms: coordinates[chrom]=array('l',[])
    for line in data:
        chrom=line[0]
        coordinate=line[1]
        coordinates[chrom].append(coordinate)
    
    return coordinates

def parse_wig(line):
    """Parse a wig line"""
    
    standard_chroms={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','V':'chrV','M':'chrM','X':'chrX'}
    
    if re.search(r'track',line): 
        try:
            description=re.search(r'description="([A-Za-z0-9]+)"\s',line).group(1)
            if not self.description: self.description=description
        except AttributeError:
            pass
    elif re.search(r'chrom=([A-Za-z0-9]+)\s',line):
        chrom=re.search(r'chrom=([A-Za-z0-9]+)\s',line).group(1)
        try:
            chrom=standard_chroms[chrom]
        except KeyError:
            pass

def estimate_wig_interval(wig):
    """When a Wig object (variable step) is given, estimate the median interval between points
        Just take a few random samples, and estimate the median interval.
    """

    chroms = wig.get_chroms()
    if not chroms: return None
    
    n_random_positions = 10
    intervals = []
    for chrom in chroms:
        len_this_chr = len(wig[chrom][0])
        a = randints(0, len_this_chr, 2 * n_random_positions)
        a.sort()
        starts = [a[i] for i in xrange(len(a)) if a%2 == 0]
        ends = [a[i] for i in xrange(len(a)) if a%2 == 1]

        for start, end in itertools.izip(starts, ends):
             intervals.append(median(diff(wig[chrom][start:end])))
             
    return median(intervals)
        
        
def select_ChIP_by_score_cutoff(ChIP, scorecutoff = 0, descend = True):
    """Select ChIP regions according to a score cutoff
    
    Parameters:
    1. ChIP: Bed object
    2. scorecutoff: the cut-off for scores
    3. descend: If True, select ChIP regions with higher scores first.
    """
    
    # if no chroms, just return None
    chroms = ChIP.get_chroms()
    if not chroms: return None
    
    # if any of chroms does not have 'score', just return None
    YESSCORE = True
    for chrom in chroms:
        YESSCORE = YESSCORE and ChIP[chrom].has_key('score')
    if not YESSCORE: return None
    
    newChIP = Bed()
    for chrom in chroms:
        if descend:
            b = map(lambda x: x>= scorecutoff, ChIP[chrom]['score'])
        else:
            b = map(lambda x: x<= scorecutoff, ChIP[chrom]['score'])
        newChIP.bed = {}
        newChIP.bed[chrom]={}
        for key in ChIP[chrom].keys():
            if key == 'start' or key =='end':
                newChIP.bed[chrom][key] = array('l', array_extractor(ChIP[chrom]['start'], b))
            elif key == 'score':
                newChIP.bed[chrom][key] = array('d', array_extractor(ChIP[chrom]['score'], b))
            else:
                newChIP.bed[chrom][key] = array_extractor(ChIP[chrom][key])
                
    return newChIP
        
        
def determine_score_cutoff_by_n_peaks(ChIP, n_peaks= 1000, descend = True):
    """Determine the score cutoff when want to only include a certain number of peaks
    
    Parameters:
    1. ChIP: a Bed object
    2. n_peaks: the number of peaks that want to be kept
    3. descend: if True, count in the highest scores first
    """
    
    # if no chroms, just return None
    chroms = ChIP.get_chroms()
    if not chroms: return None
    
    # if any of chroms does not have 'score', just return None
    YESSCORE = True
    for chrom in chroms:
        YESSCORE = YESSCORE and ChIP[chrom].has_key('score')
    if not YESSCORE: return None
    
    # concatenate the scores and sort
    scores = []
    for chrom in chroms:
        scores += ChIP[chrom]['score']
    scores.sort()
    
    # if n_peaks is more than the actual scores, just return
    if n_peaks >= len(scores): 
        try:
            if descend:
                return scores[-1]
            else:
                return scores[0]
        except IndexError:
            return None
    else:
        if descend:
            return scores[-1*n_peaks]
        else:
            return scores[n_pekas]


def select_ChIP_by_n_peaks(ChIP, n_peaks = 1000, descend=True):
    """Select n_peaks ChIP regions by score
    
    Parameters:
    1. ChIP: a Bed object
    2. n_peaks: the number of peaks to be selected
    3. descend: if True, from the top; otherwise from the bottom of scores
    """
    
     # if no chroms, just return None
    chroms = ChIP.get_chroms()
    if not chroms: return None
    
    # if any of chroms does not have 'score', just return None
    YESSCORE = True
    for chrom in chroms:
        YESSCORE = YESSCORE and ChIP[chrom].has_key('score')
    if not YESSCORE: return None
    
    # concatenate the scores and sort
    scores = []
    for chrom in chroms:
        scores += ChIP[chrom]['score']
    scores.sort()
    
    # if more peaks will be kept than we have, just return ChIP
    if n_peaks >= len(scores):
        return ChIP
    
    # get the cutoff 
    if descend:
        cutoff = scores[-1*n_peaks]
    else:
        cutoff = scores[n_peaks]
    
    # select ChIP regions        
    newChIP = Bed()
    for chrom in chroms:
        if descend:
            b = map(lambda x: x>= cutoff, ChIP[chrom]['score'])
        else:
            b = map(lambda x: x<= cutoff, ChIP[chrom]['score'])
        newChIP.bed = {}
        newChIP.bed[chrom]={}
        for key in ChIP[chrom].keys():
            if key == 'start' or key =='end':
                newChIP.bed[chrom][key] = array('l', array_extractor(ChIP[chrom]['start'], b))
            elif key == 'score':
                newChIP.bed[chrom][key] = array('d', array_extractor(ChIP[chrom]['score'], b))
            else:
                newChIP.bed[chrom][key] = array_extractor(ChIP[chrom][key], b)
                
    return newChIP 


def get_percent_from_BED(bed, percentage, twoside=False):
    """Get the lower and upper bondaries of scores given a percentage
    
    Parameters:
    1. bed: a bed object
    2. percentage: integer from 0 to 100
    3. twosided: If True, get the middle part of scores within the percentage; otherwise, started with the lowest score.
    
    WARNING: bed must have 'score'; otherwise, just [] will be returned.
    WARNING: there must be an enough number of elements
    """
    
    # if no chroms, just return None
    chroms = bed.get_chroms()
    if not chroms: return []
    
    # if any of chroms does not have 'score', just return None
    YESSCORE = True
    for chrom in chroms:
        YESSCORE = YESSCORE and bed[chrom].has_key('score')
    if not YESSCORE: return []
    
    # put percentage into the right range
    percentage = min(100, percentage)
    percentage = max(0, percentage)
    
    # concatenate the scores and sort
    scores = []
    for chrom in chroms:
        scores += bed[chrom]['score']
    scores.sort()
    
    l = len(scores)
    if twoside:
        left = int(((100-percentage)/200.0)*(l-1))
        right = l - left
    else:
        left = 0
        right = int((percentage/100.0)*(l-1))
    
    return scores[left], scores[right]
    

def get_histogram_BED_scores(bed, breaks, chrom=''):
    """Get the histogram of bed scores
    
    Parameters:
    1. bed: a Bed object
    2. breaks: bins for the histogram
    """
    
    # if no chroms, just return None
    chroms = bed.get_chroms()
    if not chroms: return {}
    
    # if any of chroms does not have 'score', just return None
    YESSCORE = True
    for chrom in chroms:
        YESSCORE = YESSCORE and bed[chrom].has_key('score')
    if not YESSCORE: return {}
    
    
    if type(breaks)==list:      # if bins are given
        histogram = {'bin':breaks[:-1]}
        if chrom:
            try:
                histogram['count']=bin(bed[chrom]['score'], breaks)
            except KeyError:
                histogram = {}
        else:
            w = []
            for chrom in chroms:
                try:
                    w += bed[chrom]['score']
                except KeyError:
                    pass
            try:
                histogram['count'] = bin(w, breaks)
            except:
                histogram = {}
                
        return histogram
                
    elif type(breaks)==int: # only number of bins are given
        histogram = {}
        if chrom:
            try:
                min_this_chr = min(bed[chrom]['score'])
                max_this_chr = max(bed[chrom]['score'])
            except KeyError:
                return {}
                
            try:
                bins = linspace(min_this_chr, max_this_chr, breaks+1)
                histogra['bin'] = bins[:-1]
                histogram['count'] = bin(bed[chrom]['score'], bins)
                return histogram
            except:
                return {}
               
        else:
            w = []
            for chrom in chroms:
                try:
                    w += bed[chrom]['score']
                except KeyError:
                    pass
        
            try:   
                min_all = min(w)
                max_all = max(w)
                bins = linspace(min_this_chr, max_this_chr, breaks+1)
                histogram['bin'] = bins
                histogram['count'] = bin(w, bins)
                return histogram
            except:
                return {}
    else:
        return {}


def check_if_yes_score(bed):
    """Check if the input bed has 'score's
    
    This function returns True only if all the chromosomes have scores
    
    Parameters:
    1. bed: a Bed object
    """
    
    try:
        chroms = bed.get_chroms()
    except AttributeError:
        return False
    
    if not chroms:
        return False
    
    YESSCORE = True
    for chrom in chroms:
        YESSCORE = YESSCORE and bed[chrom].has_key('score')
    
    return YESSCORE


def fill_up_scores_w_val(bed, val = 1):
    """Fill up scores with a value
    
    Parameters:
    1. bed: a Bed object
    2. val: a value to fill with
    """
    #if empty bed, just return None
    try:
        chroms = bed.get_chroms()
    except AttributeError:
        return
    
    if not chroms: return
    
    for chrom in chroms:
        
        if not bed[chrom].has_key('name'):
            bed[chrom]['name'] = [0] * len(bed[chrom]['start'])
        
        bed[chrom]['score'] = array('d', [1.0]*len(bed[chrom]['start']))

   
def downsample_ChIP(ChIP, interval = 5000):
    """Down-sample the ChIP regions according to the given interval
    
    Parameters:
    1. ChIP: a Bed object
    2. interval: an interval in bp. within this interval, only one ChIP region remains. Default = 5kb
    
    NOTE that the input bed object must be sorted by 'start'
    """
    
    chroms = ChIP.get_chroms()
    if not chroms: return None
    
    sampledChIP = Bed()
    sampledChIP.bed = {}
    for chrom in chroms:
        sampledChIP.bed[chrom] = {}
        
        # initialize
        keys = ChIP[chrom].keys()
        for key in keys:
            if key == 'start' or key == 'end':
                sampledChIP[chrom][key] = array('l', [])
            elif key == 'score':
                sampledChIP[chrom][key] = array('d', [])
            else:
                sampledChIP[chrom][key] = []
                
        # if there is only one ChIP region, just copy and go
        len_this_chrom = len(ChIP[chrom]['start'])
        if len_this_chrom == 1:
            sampledChIP.bed[chrom] = copy.deepcopy(ChIP.bed[chrom])
            continue
        
        # downsampling
        prev = ChIP[chrom]['start'][0]
        # update the first element for each column
        for key in keys: sampledChIP[chrom][key].append(ChIP[chrom][key][0])
        
        # iterate through this chromosome
        for i, start in itertools.izip(xrange(1,len_this_chrom), ChIP[chrom]['start'][1:]):

            # keep only if the new bed region is far enough from the previously sampled point
            if start - prev > interval:
                for key in keys: sampledChIP[chrom][key].append(ChIP[chrom][key][i])
                prev = start
    
    return sampledChIP
                
#-------------------------------------
# functions related to R script out
#------------------------------------- 
def _pval_str(pvals):
    
    pval_str=[]
    for pval in pvals:
        if pval==4.94e-324: pval_str.append('(<=%.1e)' %(pval))
        elif pval<1e-3: pval_str.append('(%.1e)' %(pval))
        else: pval_str.append('(%.3f)' %(pval))

    return pval_str

def _get_percentage(fraction):
    
    return [100*f for f in fraction]

def _percent_str(percentages):
    """Convert percentages values into string representations"""
    
    pstr = []
    for percent in percentages:
        if percent >= 0.1:
            pstr.append('%.1f' %percent+' %')
        elif percent >= 0.01:
            pstr.append('%.2f' %percent+' %')
        else:
            pstr.append('%.3f' %percent+' %')
    
    return pstr

def _cat_percent_pval_strs(percent_str,pval_str,sep='\n'):
    
    return [sep.join(a) for a in zip(percent_str,pval_str)]

         
def _draw_barplot_w_val(genome,ChIP,pval,names,main=None,xlab=None,ylab=None,horiz=False,upsidedown=False, cex_names = 1, cex_text = 1):
    """Draw a barplot with p values
    
    Parameters:
    1. genome: a list of genome background annotation results.
    2. ChIP: a list of ChiP annotation results.
    3. pval: a list of p values.
    4. name: bar plot name
    5. main: the mail title
    6. xlab: x axis label
    7. ylab: y axis label
    8. horiz: horizontal or vertical
    9. upsidedown: reverse the order of bars
    """
    
    # draw bar plot
    # get margins
    maxp=max(max(genome),max(ChIP))
    marg=maxp/1.2 # margin for writing percentage and p value. Thus, 1/4 of the height of the plot is taken by texts                      

    if horiz:
        xlim=[0,maxp+marg]
        ylim=None
        if upsidedown:
            genome.reverse()
            ChIP.reverse()
            pval.reverse()
            names.reverse()
    else:
        xlim=None
        ylim=[0,maxp+marg]
    rscript=R.barplot([genome,ChIP],names=names,beside=True,horiz=horiz,col=['lightblue','mistyrose'],\
                    main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim, cex_names=cex_names)
    
    # write p values
    pct_genome=_percent_str(genome)
    pct_ChIP=_percent_str(ChIP)
    pval_ChIP=_pval_str(pval)
    if horiz:
        pct_pval_ChIP=_cat_percent_pval_strs(pct_ChIP,pval_ChIP,sep=' ')
        x1=genome
        x2=ChIP
        y1='mp[1,]'
        y2='mp[2,]'
        pos=4
    else:
        pct_pval_ChIP=_cat_percent_pval_strs(pct_ChIP,pval_ChIP,sep='\n')
        x1='mp[1,]'
        x2='mp[2,]'
        y1=genome
        y2=ChIP
        pos=3
    
    # write percentages and pvalues   
    if cex_text == 1:
        rscript+=R.text(x1,y1,pct_genome,pos=pos,offset=0.2)
        rscript+=R.text(x2,y2,pct_pval_ChIP,pos=pos,offset=0.2)
    else:
        rscript+=R.text(x1,y1,pct_genome,pos=pos,offset=0.2,cex=cex_text)
        rscript+=R.text(x2,y2,pct_pval_ChIP,pos=pos,offset=0.2,cex=cex_text)

    return rscript

def _draw_pie(x,names,main,pval=None,border=False, cols=None):
    """Draw a pie chart """
    
    percentages=_percent_str([i*100.0 for i in x])
    #cols=["#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0"]
    if not pval:
        labels=['%s: %s' %(n,pct) for n,pct in zip(names,percentages)]
    else:
        pval_str=_pval_str(pval)
        labels=['%s: %s %s' %(n,pct,pv) for n,pct,pv in zip(names,percentages,pval_str)]
    
    if type(cols) == list or type(cols) ==tuple:
        rscript=R.pie(x,labels=labels,main=main,col=cols[:len(x)],clockwise=True,border=border)
    elif type(cols) == str:
        rscript=R.pie(x,labels=labels,main=main,col=cols,clockwise=True,border=border)
    else:
        rscript=R.pie(x,labels=labels,main=main,clockwise=True,border=border)
        
    return rscript

def draw_rel_loc_pie(genome_p,ChIP_p,pval,gene_div):
    
    label=[str(i)+'/'+str(gene_div[0]) for i in range(1,gene_div[0]+1)]
    rscript=_draw_pie(genome_p['whole']['rel_loc'][0],label,r'Relative Location Gene (%d)\nGenome' %gene_div[0])
    rscript+=_draw_pie(ChIP_p['whole']['rel_loc'][0],label,r'\nChIP (p-value)',pval=pval['whole']['rel_loc'][0])
    label=[str(i)+'/'+str(gene_div[1]) for i in range(1,gene_div[1]+1)]
    rscript+=_draw_pie(genome_p['whole']['rel_loc'][1],label,r'Relative Location Gene (%d)\nGenome' %gene_div[1])
    rscript+=_draw_pie(ChIP_p['whole']['rel_loc'][1],label,r'\nChIP (p-value)',pval=pval['whole']['rel_loc'][1])
    
    return rscript

def draw_rel_loc_cds_pie(genome_p,ChIP_p,pval,gene_div):
    
    label=[str(i)+'/'+str(gene_div[0]) for i in range(1,gene_div[0]+1)]
    rscript=_draw_pie(genome_p['whole']['rel_loc_cds'][0],label,r'Relative Location CDS (%d)\n(Genome)' %gene_div[0])
    rscript+=_draw_pie(ChIP_p['whole']['rel_loc_cds'][0],label,r'\nChIP (p-value)',pval=pval['whole']['rel_loc_cds'][0])
    label=[str(i)+'/'+str(gene_div[1]) for i in range(1,gene_div[1]+1)]
    rscript+=_draw_pie(genome_p['whole']['rel_loc_cds'][1],label,r'Relative Location CDS (%d)\nGenome' %gene_div[1])
    rscript+=_draw_pie(ChIP_p['whole']['rel_loc_cds'][1],label,r'\nChIP (p-value)',pval=pval['whole']['rel_loc_cds'][1])

    return rscript


def draw_pie_distribution_of_elements(ChIP_EP, prom=3000, down=3000):
    """Draw the pie charts of the overall distributions
    """
    
    # get the labels (legend) for the pie chart
    names = ['Promoter (=<%d bp)' %(prom/3), 'Promoter (%d-%d bp)' %(prom/3, 2*prom/3), 'Promoter (%d-%d bp)' %(2*prom/3, prom),\
             'Downstream (<=%d bp)' %(prom/3), 'Downstream (%d-%d bp)' %(prom/3, 2*prom/3), 'Downstream (%d-%d bp)' %(2*prom/3, prom),\
             "5'UTR","3'UTR", "Coding exon", "Intron",\
             'Distal intergenic']
     
    
    # get the proportions to draw
    x = ChIP_EP['whole']['promoter'] + ChIP_EP['whole']['downstream'] + ChIP_EP['whole']['gene'] + [ChIP_EP['whole']['enhancer']]
    x_percent = _percent_str(map(lambda a: 100.0*a, x))
    names_w_percent = map(lambda x, y: x + ': ' + y, names, x_percent)
    # make x values less than .1% .5% because they are too small to see in the pie chart. But x_percent does not change
    x = map(max, x, [0.01]*len(x))
    
    #
    # producing R script return
    #
    
    # put header
    rscript = '\n'
    rscript += R.comment('')
    rscript += R.comment('Distribution of ChIP regions over cis-regulatory element')
    rscript += R.comment('Note that the x may be modified for better graphics in case a value is too small')
    rscript += R.comment('Thus, look at the labels of the pie chart to get the real percentage values' )
    rscript += R.comment('')
    rscript += '\n'
    
    # some graphical parameters
    s = 0.7
    v = 0.8
    init_angle = 45
    main = 'Distribution of ChIP Regions'
    colval = 'piecol'
    #newcolval = 'newpiecol'
    mar = mar=[4,4,5,3.8]
    oma=[4,2,4,2]
    mfrow = [2, 1]
    rscript += R.par(mar=mar, oma=oma, mfrow=mfrow)
    
    # R script
    rscript += R.rainbow(len(x), s = s, v = v, rval=colval)
    #rscript += '%s <- piecol[c(%s)]\n' %(newcolval, str(map(lambda a: a+1, newix))[1:-1])
    rscript += R.pie(x, labels="NA", main=main, col=colval,clockwise=True, border=colval, radius=0.9,init_angle=init_angle, cex=0.8)
    rscript += R.plot([0,1],[0,1], tp="n", axes=False, xlab="", ylab="", main="", frame=False)
    rscript += R.legend(x='top', legend=names_w_percent, pch=15, col=colval, bty="n")
    
    return rscript


def draw_chrom_barplot(genome_p,ChIP_p,pval,horiz=True,upsidedown=True):
    
    # get chromosome
    chroms=set([chrom for chrom in genome_p.get_chroms() if chrom!='whole']).intersection(set([chrom for chrom in ChIP_p.get_chroms() if chrom != 'whole']))
    chroms = sort_chroms(list(chroms))
    n_chroms = len(chroms)
    
    # if no chromosomes shaped between genome and ChIP, return null rscript
    if not chroms: return ''
    
    # get percentages of the chromsomes
    genome_chrom_p=[]
    ChIP_chrom_p=[]
    pval_chrom=[]
    for chrom in chroms:
        genome_chrom_p.append(genome_p[chrom]['chroms'])
        ChIP_chrom_p.append(ChIP_p[chrom]['chroms'])
        pval_chrom.append(pval[chrom]['chroms'])
    
    # turn into percentages    
    gp=_get_percentage(genome_chrom_p)
    cp=_get_percentage(ChIP_chrom_p)
    chroms = [chrom.rsplit('chr',1)[-1] for chrom in chroms]    # get rid of 'chr' from chromosome name
    
    # adjust the font size for percentages and p values
    cex_text = lininterpol([6, 0.9], [24, 0.7], n_chroms)
    cex_text = max(min(0.9, cex_text), 0.7)
    
    # adjust the font size for bar names
    cex_names = lininterpol([6, 1], [24, 0.8], n_chroms)
    cex_names = max(min(1, cex_names), 0.8)
    
    # draw barplot
    rscript=_draw_barplot_w_val(gp,cp,pval_chrom,chroms,'Chromosomal Distribution of ChIP Regions',xlab='Percentage %',ylab='Chromosome',horiz=horiz,upsidedown=upsidedown, cex_names = cex_names, cex_text=cex_text)
    rscript+=R.legend(x='right',legend=['Genome','ChIP (p-value)'],pch=15,col=['lightblue','mistyrose'],bty='n')
    return rscript

def draw_promoter_barplot(genome_p,ChIP_p,pval,prom):
    
    names=['<=%d bp' %(prom/3),'<=%d bp' %(2*prom/3),'<=%d bp' %prom]
    gp=_get_percentage(genome_p['whole']['promoter'])
    cp=_get_percentage(ChIP_p['whole']['promoter'])
    rscript=_draw_barplot_w_val(gp,cp,pval['whole']['promoter'],names,'Promoter',ylab='Percentage %')
    rscript+=R.legend(x='topleft',legend=['Genome','ChIP (p-value)'],pch=15,col=['lightblue','mistyrose'],bty='n')
    return rscript

def draw_bipromoter_barplot(genome_p,ChIP_p,pval,biprom):
    
    names=['<=%d bp' %(biprom/2),'<=%d bp' %biprom]
    gp=_get_percentage(genome_p['whole']['bipromoter'])
    cp=_get_percentage(ChIP_p['whole']['bipromoter'])
    rscript=_draw_barplot_w_val(gp,cp,pval['whole']['bipromoter'],names,'Bidirectional Promoter',ylab='Percentage %')
    rscript+=R.legend(x='topleft',legend=['Genome','ChIP (p-value)'],pch=15,col=['lightblue','mistyrose'],bty='n')
    return rscript

def draw_downstream_barplot(genome_p,ChIP_p,pval,down):
    
    names=['<=%d bp' %(down/3),'<=%d bp' %(2*down/3),'<=%d bp' %down]
    gp=_get_percentage(genome_p['whole']['downstream'])
    cp=_get_percentage(ChIP_p['whole']['downstream'])
    rscript=_draw_barplot_w_val(gp,cp,pval['whole']['downstream'],names,'Downstream',ylab='Percentage %')
    rscript+=R.legend(x='topleft',legend=['Genome','ChIP (p-value)'],pch=15,col=['lightblue','mistyrose'],bty='n')
    return rscript

def draw_gene_barplot(genome_p,ChIP_p,pval):
    
    names=["5'UTR","3'UTR","Coding Exon","Intron","All"]
    gp=_get_percentage(genome_p['whole']['gene'])
    cp=_get_percentage(ChIP_p['whole']['gene'])
    rscript=_draw_barplot_w_val(gp,cp,pval['whole']['gene'],names,'Gene',ylab='Percentage %')
    rscript+=R.legend(x='topleft',legend=['Genome','ChIP (p-value)'],pch=15,col=['lightblue','mistyrose'],bty='n')
    return rscript

def draw_roi_barplot(genome_p,ChIP_p,pval):
    
    names=['']
    gp=_get_percentage([genome_p['whole']['roi']])
    cp=_get_percentage([ChIP_p['whole']['roi']])
    rscript=_draw_barplot_w_val(gp,cp,[pval['whole']['roi']],names,'Regions of Interest',ylab='Percentage %')
    rscript+=R.legend(x='topleft',legend=['Genome','ChIP (p-value)'],pch=15,col=['lightblue','mistyrose'],bty='n')
    return rscript


def draw_colorbar(lim, n_cols = 128, cmap = 'cmap', vertical=False, hist=None, linecol = 'cyan'):
    """Draw a colorbar"""
    
    # if you want to draw a vertical color bar
    if vertical:
        left=0
        right=1
        mid=sum(lim)/2.0
        rscript=R.plot([left,right],lim,ylim=lim,xlim=[left,right],tp="n",frame=False, xaxt="n", xlab="", ylab="", main="")
        rscript+=R.polygon([left,right,right,left],[lim[0],lim[0],lim[1],lim[1]])
        vals = linspace(lim[0], lim[1], max(2, n_cols))
        rscript+=R.heatmap_bar(vals, vals, lim, lim, bot=left, top=right, vertical=vertical, cmap=cmap)
    # if you want to draw a horizontal color bar
    else:
        top=1
        bot=0
        mid=sum(lim)/2.0
        rscript=R.plot(lim,[bot,top], xlim=lim, ylim=[bot,top], tp="n", frame=False, yaxt="n", xlab="", ylab="",main='')
        rscript+=R.polygon([lim[0],lim[1],lim[1],lim[0]],[bot,bot,top,top])
        vals = linspace(lim[0], lim[1], max(2, n_cols))  
        rscript+=R.heatmap_bar(vals, vals, lim, lim, bot=bot, top=top, cmap = cmap)
        # if histogram is given
        if hist:
            N = 0.8                         # ratio for normalization
            max_count=max(hist['count'])    # get the maximum count and normalize the other counts with it
            probability = [N * h / max_count for h in hist['count']]    # normalize
            within_lim = arglim(hist['bin'],lim)  # get the index of elements within lim
            bins2draw = array_extractor(hist['bin'], within_lim)    # extract the corresponding bins
            prob2draw = array_extractor(probability, within_lim)    # extract the corresponding probabilities
            rscript += R.lines(bins2draw, prob2draw, col=[linecol], xlim=lim, ylim=[bot, top])  # draw using R.lines
    
    return rscript


def draw_histogram_bar(histogram, lim, main ='', xlab = '', linecol='cyan', bgcol='black', cex=1):
    """Draw a histogram bar
    
    Parameters:
    1. histogram: a dictionary of histogram. This dictionary must have 'bin' and 'count', which must be the same length
    2. lim: lower and upper boundary for drowing histogram
    3. linecol: line color. Default = cyan
    4. bgcol: background color. Default = black
    """

    top=1
    bot=0
    mid=sum(lim)/2.0
    if cex == 1:
        rscript=R.plot(lim,[bot,top], xlim=lim, ylim=[bot,top], tp="n", frame=False, yaxt="n", xlab=xlab, ylab="",main=main)
    else:
        rscript=R.plot(lim,[bot,top], xlim=lim, ylim=[bot,top], tp="n", frame=False, yaxt="n", xlab=xlab, ylab="",main=main, cex=cex)
    rscript+=R.polygon([lim[0],lim[1],lim[1],lim[0]],[bot,bot,top,top], col = [bgcol])
    
    # draw a histogram
    N = 0.8                         # ratio for normalization
    max_count=max(histogram['count'])    # get the maximum count and normalize the other counts with it
    probability = [N * h / max_count for h in histogram['count']]    # normalize
    within_lim = arglim(histogram['bin'],lim)  # get the index of elements within lim
    bins2draw = array_extractor(histogram['bin'], within_lim)    # extract the corresponding bins
    prob2draw = array_extractor(probability, within_lim)    # extract the corresponding probabilities
    rscript += R.lines(bins2draw, prob2draw, col=[linecol], xlim=lim, ylim=[bot, top], lwd =2)  # draw using R.lines
    
    return rscript
      
def draw_nullbar():
    """Draw a null bar"""
    lim=[-3, 3] # some dummy value
    top=1
    bot=0
    rscript=R.plot(lim,[bot,top], xlim=lim, ylim=[bot,top], tp="n", axes=False, xlab="", ylab="", main='')
    
    return rscript


def draw_chrom_frame(chrominfo, chrom, maxchrom, centromere):
    """Draw the chromosome frame.
    
    Parameters:
    1. chrominfo: a dictionary of a chromosome (BED[chrom])
    2. centromere: """
    bot=-1
    top=1
    mid=(bot+top)/2
    rscript=R.plot(chrominfo,[bot,top], xlim=maxchrom, ylim=[bot,top], tp="n", frame=False, yaxt="n", xlab="", ylab="", main=chrom)
    if centromere:
        p_start=centromere[0][0]
        p_end=centromere[1][0]
        q_start=centromere[0][1]
        q_end=centromere[1][1]
        rscript+=R.polygon([chrominfo[0], p_start, p_end, p_start,chrominfo[0],chrominfo[0]], [bot, bot, mid, top, top, top])
        rscript+=R.polygon([q_start, q_end, chrominfo[1], chrominfo[1], q_end, q_start], [mid, bot, bot, top, top, mid])
        rscript+=R.polygon([p_start, p_end, p_start], [bot, mid, top], col="gray40")
        rscript+=R.polygon([q_start, q_end, q_end],[mid,bot,top],col="gray40")
    else:
        rscript+=R.polygon([chrominfo[0], chrominfo[1], chrominfo[1], chrominfo[0]], [bot, bot, top, top])
    
    return rscript
        
    

def draw_a_chrom_heatmap(wig, chrom, ylim, chrominfo, maxchrom, centromere=None, bot = -1, top = 1, line=True, linecol='cyan'):
    """Draw a chromosome heatmap"""
        
    # draw the frame
    rscript = ''
    #rscript=draw_chrom_frame(chrominfo, chrom, maxchrom, centromere)
        
    if centromere:
        # get the parts of the chromosome before and after centromere
        wig_bef_cent=[[], []]
        wig_aft_cent=[[], []]
        i=0
        between=[]
        for c,v in intertools.izip(wig[0],wig[1]):
            if c < centromeres[0][0]:
               wig_bef_cent[0].append(c)
               wig_bef_cent[1].append(v)
            elif c > centromeres[1][1]:
                wig_aft_cent[0].append(c)
                wig_aft_cent[1].append(v)
            else:
                between=[c,v]
        
        wig_bef_cent[0].append(centromeres[0][0])
        if between:
            wig_aft_cent[0].insert(0,centromeres[1][1])
            wig_aft_cent[1].insert(0,between[1])
        
        rscript+=R.heatmap_bar(wig_bef_cent[0],wig_bef_cent[1],[chrominfo[0],centromere[0][0]],ylim,bot=bot,top=top)
        rscript+=R.heatmap_bar(wig_aft_cent[0],wig_aft_cent[1],[centromere[1][1],chrominfo[1]],ylim,bot=bot,top=top)
            
    else:
        rscript += R.heatmap_bar(wig[0], wig[1], chrominfo, ylim, bot=bot, top=top)
        if line: 
            med = (top + bot)/2         # med position to draw the line
            hdelta = (top - bot)/2      # get the height from the mid to the top
            scaling = 0.8               # scaling factor to shrink the line
            scaled_wig =  map(lambda x: hdelta* scaling * ((max(min(ylim[1],x), ylim[0]) - ylim[1] - ylim[0])/(ylim[1] - ylim[0])) + med, wig[1])
            rscript += R.lines(wig[0], scaled_wig, col=[linecol])
        
    return rscript
  
    
def draw_single_profile(breaks,prfl,col=[],main='',xlab='',ylab='',ylim=[],v=None):
    
    if not col: col=['red']
    if ylim: rscript=R.plot(breaks,prfl,col=col,main=main,xlab=xlab,ylab=ylab,ylim=ylim,lwd=2)
    else: rscript=R.plot(breaks,prfl,col=col,main=main,xlab=xlab,ylab=ylab,lwd=2)
    
    if v!=None: rscript+=R.abline(v=v,lty=2,col=['black'])
    
    return rscript


def draw_multiple_profiles(breaks, prfls, cols=[], main='', xlab='', ylab='', ylim=[], v=None, legends=[]):
    
    n_prfls = len(prfls)
    if not ylim: ylim=[min(min_col_by_col(prfls)),max(max_col_by_col(prfls))]
    
    if cols:
        rscript +=R.plot(breaks,prfls[0],col=cols[0],main=main,xlab=xlab,ylab=ylab,ylim=ylim,lwd=2) 
        for i in range(1, n_prfls):
            rscript+=R.lines(breaks,prfls[i],col=[cols[i]],lwd=2)
        rscript+=R.legend(x='topleft', legend=legends, pch=15, col=col[:len(prfls)], bty='o')  
            
    else:
        rval = 'linecols'
        rscript = R.rainbow(n_prfls-1, rval='%s' %rval)
        rscript += '%s <- c(%s, "black")\n' %(rval, rval)
        rscript +=R.plot(breaks,prfls[0],col='%s[1]' %rval, main=main, xlab=xlab, ylab=ylab, ylim=ylim, lwd=2) 
        for i in range(1, n_prfls):
            rscript += R.lines(breaks, prfls[i], col='%s[%d]' %(rval, i+1), lwd=2)
        if not legends:
            legends=['Group %d' %i for i in range(1,len(prfls)+1)]
        rscript+=R.legend(x='topleft', legend=legends, pch=15, col='%s' %(rval), bty='o')
        
    if v!=None: rscript+=R.abline(v=v,lty=2,col=['black'])
    
    return rscript


def draw_profile_plot(breaks,avg_upstream,avg_downstream,metagene_breaks,avg_metagene,metacatexon_breaks,avg_metacatexon,metacatintron_breaks,avg_metacatintron,metagene_breaks_lim=[-1000,1000]):
    
    #comment
    R.comment('')
    R.comment('Draw wig profiles')
    R.comment('')
    
    rscript=R.par(mar=[4,4,5,3.8],oma=[4,2,4,2])
    mat=[[1,2],[3,3],[4,5]]
    rscript+=R.layout(mat,widths=[1,1],heights=[1,1,1])
    
    # for metagene plot
    step = breaks[1]-breaks[0]
    start_up,end_up=where2(max(metagene_breaks_lim[0],breaks[0]), 0, breaks, step)
    start_down,end_down=where2(step, min(metagene_breaks_lim[1]+1,breaks[-1]+1), breaks, step)
    
    # KEEP THIS IN MIND: the first and last points of metagene are not used in metagene plots
    comb_breaks=breaks[start_up:end_up] + metagene_breaks[1:-1] + map(lambda b: b + metagene_breaks[-1], breaks[start_down:end_down])
    comb_metagene=avg_upstream[start_up:end_up] + avg_metagene[1:-1] + avg_downstream[start_down:end_down]  
    
    # draw
    #ylim=[min(avg_upstream)-abs(min(avg_upstream))*0.1,max(avg_upstream)+abs(max(avg_upstream))*0.1]
    rscript+=draw_single_profile(breaks,avg_upstream,main='Average Profile near TSS',xlab='Relative Distance to TSS (bp)',ylab='Average Profile',v=0)#,ylim=ylim)
    
    #ylim=[min(avg_downstream)-abs(min(avg_downstream))*0.1,max(avg_downstream)+abs(max(avg_downstream))*0.1]
    rscript+=draw_single_profile(breaks,avg_downstream,main='Average Profile near TTS',xlab='Relative Distance to TTS (bp)',ylab='Average Profile',v=0)#,ylim=ylim)
    
    #ylim=[min(comb_metagene)-abs(min(comb_metagene))*0.1,max(comb_metagene)+abs(max(comb_metagene))*0.1]
    rscript+=draw_single_profile(comb_breaks,comb_metagene,main='Average Gene Profile',xlab='Upstream (bp), %d bp of Meta-gene, Downstream (bp)' %int(metagene_breaks[-1]),ylab='Average Profile')#,ylim=ylim)
    rscript+=R.abline(v=metagene_breaks[0],lty=2,col=['black'])
    rscript+=R.abline(v=metagene_breaks[-1],lty=2,col=['black'])
    
    ylim=[min(min_w_nan(avg_metacatexon),min_w_nan(avg_metacatintron)),max(max_w_nan(avg_metacatexon),max_w_nan(avg_metacatintron))]
    rscript+=draw_single_profile(metacatexon_breaks,avg_metacatexon,main='Average Concatenated Exon Profile',xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim)
    rscript+=draw_single_profile(metacatintron_breaks,avg_metacatintron,main='Average Concatenated Intron Profile',xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim)
    
    return rscript


def draw_exon_intron_profile_plot(metaexon_breaks,avg_metaexon,metaintron_breaks,avg_metaintron,elowers,euppers,ilowers,iuppers):    
    
    n_rows = 3
    # exon and intron profiles
    rscript = ''
    
    temp_metaexon = filter(lambda x: x, avg_metaexon)
    temp_metaintron = filter(lambda x: x, avg_metaintron)
    
    # ge the global minimum and global maximum
    ylim = [0, 0]
    if temp_metaexon:
        mins_avg_metaexon = map(min_w_nan, temp_metaexon)
        maxs_avg_metaexon = map(max_w_nan, temp_metaexon)
    else:
        mins_avg_metaexon = []
        maxs_avg_metaexon = []
    
    if temp_metaintron:
        mins_avg_metaintron = map(min_w_nan, temp_metaintron)
        maxs_avg_metaintron = map(max_w_nan, temp_metaintron)
    else:
        mins_avg_metaintron = []
        maxs_avg_metaintron = []
    
    try:
        ylim[0] = min(mins_avg_metaexon + mins_avg_metaintron)
    except ValueError:
        ylim[0] = 0
        
    try:
        ylim[1] = max(maxs_avg_metaexon + maxs_avg_metaintron)
    except ValueError:
        ylim[1] = 0
        
    ylim = [ylim[0], ylim[1]+0.2*ylim[1]]
    row = 0
    #for metaexon_break, metaintron_break, ame, ami, elower, eupper, ilower, iupper in itertools.izip(metaexon_breaks, metaintron_breaks, avg_metaexon, avg_metaintron, elowers, euppers, ilowers, iuppers):
    for metaexon_break, metaintron_break, ame, ami, elower, eupper, ilower, iupper in map(None, metaexon_breaks, metaintron_breaks, avg_metaexon, avg_metaintron, elowers, euppers, ilowers, iuppers):
        if row % n_rows == 0: rscript+=R.par(mar=[4,4,5,3.8],oma=[4,2,4,2],mfrow=[n_rows,2])
        if ame:
            rscript+=draw_single_profile(metaexon_break, ame, main='Average Exon Profile\n(%d <= length < %d bp)' %(elower, eupper),xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim)
        else:   # when nothing, just draw null plot
            rscript+=draw_nullbar()
            
        if ami:
            rscript+=draw_single_profile(metaintron_break, ami,main='Average Intron Profile\n(%d <= length < %d bp)' %(ilower, iupper),xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim)
        else:
            rscript+=draw_nullbar()
        row += 1
        
    return rscript


def draw_profile_plots(breaks,avg_upstreams,avg_downstreams,metagene_breaks,avg_metagenes,metaexon_breaks,avg_metaexons,metaintron_breaks,avg_metaintrons,metagene_breaks_lim=[-1000,1000],legends=[]):

    # comments
    R.comment('')
    R.comment('Draw wig profiles')
    R.comment('')
    
    # layout
    rscript=R.par(mar=[4,4,5,3.8],oma=[4,2,4,2])
    mat=[[1,2],[3,3],[4,5]]
    rscript+=R.layout(mat,widths=[1,1],heights=[1,1,1])
    n_profiles=len(avg_upstreams)
    
    # legends
    if legends:
        legends+=['All']
    else:
        legends=['Group %d' %i for i in range(1,n_profiles)]+['All']
    
    # colors
    #default_cols=['red','blue','hotpink','cyan','green','bown']
    #cols=default_cols[:n_profiles-1]+['black']
    
    # for metagene plot
    start_up,end_up=where(max(metagene_breaks_lim[0],breaks[0]),0,breaks)
    start_down,end_down=where(1,min(metagene_breaks_lim[1]+1,breaks[-1]+1),breaks)
    comb_breaks=breaks[start_up:end_up]+metagene_breaks+[b+metagene_breaks[-1] for b in breaks[start_down:end_down]]
    comb_metagenes=[up[start_up:end_up]+mg+down[start_down:end_down] for up,mg,down in zip(avg_upstreams,avg_metagenes,avg_downstreams)]
    
    ylim =[0, 0]
    no_empty_upstreams = filter_out_empty_rows(avg_upstreams)
    if no_empty_upstreams:
        ylim=[min(min_row_by_row(no_empty_upstreams,True))-abs(min(min_row_by_row(no_empty_upstreams,True)))*0.1,max(max_row_by_row(no_empty_upstreams,True))+abs(max(max_row_by_row(no_empty_upstreams,True)))*0.1]
    rscript+=draw_multiple_profiles(breaks,avg_upstreams,main='Average Profiles near TSS',xlab='Relative Distance (bp)',ylab='Average Profile',v=0,legends=legends,cols=[], ylim=ylim)
    
    ylim = [0, 0]
    no_empty_downstreams = filter_out_empty_rows(avg_downstreams)
    if no_empty_downstreams:
        ylim=[min(min_row_by_row(no_empty_downstreams,True))-abs(min(min_row_by_row(no_empty_downstreams,True)))*0.1,max(max_row_by_row(no_empty_downstreams,True))+abs(max(max_row_by_row(no_empty_downstreams,True)))*0.1]
    rscript+=draw_multiple_profiles(breaks,avg_downstreams,main='Average Profiles near TTS',xlab='Relative Distance (bp)',ylab='Average Profile',v=0,legends=legends, cols=[], ylim=ylim)
    
    ylim = [0, 0]
    no_empty_metagenes = filter_out_empty_rows(comb_metagenes)
    if no_empty_metagenes:
        ylim=[min(min_row_by_row(no_empty_metagenes,True))-abs(min(min_row_by_row(no_empty_metagenes,True)))*0.1,max(max_row_by_row(no_empty_metagenes,True))+abs(max(max_row_by_row(no_empty_metagenes,True)))*0.1]
    rscript+=draw_multiple_profiles(comb_breaks,comb_metagenes,main='Average Gene Profiles',xlab='Upstream (bp), %d bp of Meta-gene, Downstream (bp)' %int(metagene_breaks[-1]),ylab='Average Profile',legends=legends, cols=[], ylim=ylim)
    rscript+=R.abline(v=metagene_breaks[0],lty=2,col=['black'])
    rscript+=R.abline(v=metagene_breaks[-1],lty=2,col=['black'])
    
    ylim = [0, 0]
    no_empty_metaexons = filter_out_empty_rows(avg_metaexons)
    no_empty_metaintrons = filter_out_empty_rows(avg_metaintrons)
    if no_empty_metaexons and no_empty_metaintrons:
        ylim=[min(min(min_row_by_row(no_empty_metaexons,True)),min(min_row_by_row(no_empty_metaintrons,True))),max(max(max_row_by_row(no_empty_metaexons,True)),max(max_row_by_row(no_empty_metaintrons,True)))+abs(max(max(max_row_by_row(no_empty_metaexons,True)),max(max_row_by_row(no_empty_metaintrons,True))))*0.2]
    elif not no_empty_metaexons and no_empty_metaintrons:
        ylim = [min(min_row_by_row(no_empty_metaintrons,True))-abs(min(min_row_by_row(no_empty_metaintrons,True)))*0.1,max(max_row_by_row(no_empty_metaintrons,True))+abs(max(max_row_by_row(no_empty_metaintrons,True)))*0.1]
    elif no_empty_metaexons and not no_empty_metaintrons:
        ylim = [min(min_row_by_row(no_empty_metaexons,True))-abs(min(min_row_by_row(no_empty_metaexons,True)))*0.1,max(max_row_by_row(no_empty_metaexons,True))+abs(max(max_row_by_row(no_empty_metaexons,True)))*0.1]
    rscript+=draw_multiple_profiles(metaexon_breaks,avg_metaexons,main='Average Concatenated Exon Profiles',xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim,cols=[], legends=legends)
    rscript+=draw_multiple_profiles(metaintron_breaks,avg_metaintrons,main='Average Concatenated Intron Profiles',xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim,cols=[], legends=legends)
    
    return rscript


def draw_exon_intron_profile_plots(metaexon_breaks,avg_metaexons,metaintron_breaks,avg_metaintrons,elowers, euppers, ilowers, iuppers, legends=[]):
    
    n_rows = 3  # number of rows per page
    rscript = ''

    # ge the global minimum and global maximum
    ylim = [0, 0]
    temp_metaexons = map(lambda x: filter(lambda y: y, x), avg_metaexons)
    temp_metaexons = filter(lambda x: x,temp_metaexons)
    if temp_metaexons:            
        mins_avg_metaexons = map(min, map(min_row_by_row, temp_metaexons,[True]*len(temp_metaexons)))
        maxs_avg_metaexons = map(max, map(max_row_by_row, temp_metaexons,[True]*len(temp_metaexons)))
    else:
        mins_avg_metaexons = []
        maxs_avg_metaexons = []
    
    temp_metaintrons = map(lambda x: filter(lambda y: y, x), avg_metaintrons)
    temp_metaintrons = filter(lambda x: x,temp_metaintrons)
    if temp_metaintrons:
        mins_avg_metaintrons = map(min, map(min_row_by_row, temp_metaintrons,[True]*len(temp_metaintrons)))
        maxs_avg_metaintrons = map(max, map(max_row_by_row, temp_metaintrons,[True]*len(temp_metaintrons))) 
    else:
        mins_avg_metaintrons = []
        maxs_avg_metaintrons = []
    
    try:
        ylim[0] = min(mins_avg_metaexons + mins_avg_metaintrons)
    except ValueError:
        ylim[0] = 0
        
    try:
        ylim[1] = max(maxs_avg_metaexons + maxs_avg_metaintrons)
    except ValueError:
        ylim[1] = 0
        
    #ylim = [min(map(min, map(min_row_by_row, avg_metaexons,[True]*n_ranges)) + map(min, map(min_row_by_row, avg_metaintrons,[True]*n_ranges))), max(map(max, map(max_row_by_row, avg_metaexons,[True]*n_ranges)) + map(max, map(max_row_by_row, avg_metaintrons,[True]*n_ranges)))]
    ylim = [ylim[0], ylim[1]+0.2*ylim[1]]
    rows = 0
    #for metaexon_break, metaintron_break, ame, ami, elower, eupper, ilower, iupper in itertools.izip(metaexon_breaks, metaintron_breaks, avg_metaexons, avg_metaintrons, elowers, euppers, ilowers, iuppers):
    for metaexon_break, metaintron_break, ame, ami, elower, eupper, ilower, iupper in map(None, metaexon_breaks, metaintron_breaks, avg_metaexons, avg_metaintrons, elowers, euppers, ilowers, iuppers): 
        if rows % n_rows == 0: rscript+=R.par(mar=[4,4,5,3.8],oma=[4,2,4,2],mfrow=[n_rows,2])
        if ame:
            rscript+=draw_multiple_profiles(metaexon_break,ame,main='Average Exon Profiles\n(%d <= length < %d bp)' %(elower, eupper),xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim,cols=[],legends=legends)
        else:
            rscript+=draw_nullbar()
            
        if ami:
            rscript+=draw_multiple_profiles(metaintron_break,ami,main='Average Intron Profiles\n(%d <= length < %d bp)' %(ilower, iupper),xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim,cols=[],legends=legends)
        else:
            rscript+=draw_nullbar()
        rows += 1
        
    return rscript


def draw_metaexon_metaintron_profile_plot(breaks,avg_upstream,avg_downstream,metaexon_breaks,avg_metaexon,metaintron_breaks,avg_metaintron):

    # dimension: 3 by 2 the first two rows have only single column and the 3rd row has two columns
    rscript=R.par(mar=[4,4,5,3.8],oma=[4,2,4,2])
    mat=[[1,1],[2,2],[3,4]]
    rscript+=R.layout(mat,widths=[1,1],heights=[1,1,1])
    
    rscript+=draw_single_profile(breaks,avg_upstream,main='Average Profile near TSS',xlab='Relative Distance to TSS (bp)',ylab='Average Profile',v=0)
    rscript+=draw_single_profile(breaks,avg_downstream,main='Average Profile near TTS',xlab='Relative Distance to TTS (bp)',ylab='Average Profile',v=0)
    ylim=[min(min(avg_metaexon),min(avg_metaintron)),max(max(avg_metaexon),max(avg_metaintron))]
    rscript+=draw_single_profile(metaexon_breaks,avg_metaexon,main='Average Exon Profile',xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim,v=0)
    rscript+=draw_single_profile(metaintron_breaks,avg_metaintron,main='Average Intron Profile',xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim,v=0)

    return rscript


def draw_metaexon_metaintron_profile_plots(breaks,avg_upstreams,avg_downstreams,metaexon_breaks,avg_metaexons,metaintron_breaks,avg_metagenes,legends=[]):

    rscript=R.par(mar=[4,4,5,3.8],oma=[4,2,4,2])
    mat=[[1,1],[2,2],[3,4]]
    rscript+=R.layout(mat,widths=[1,1],heights=[1,1,1])
    rscript+=draw_multiple_profiles(breaks,avg_upstreams,main='Average Profiles near TSS',xlab='Relative Distance (bp)',ylab='Average Profile',v=0,cols=[],legends=legends)
    rscript+=draw_multiple_profiles(breaks,avg_upstreams,main='Average Profilesnear TTS',xlab='Relative Distance (bp)',ylab='Average Profile',v=0,cols=[],legends=legends)
    
    ylim=[min(min(min_col_by_col(avg_metaexon)),min(min_col_by_col(avg_metaintron))),max(max(max_col_by_col(avg_metaexon)),max(max_col_by_col(avg_metaintron)))]
    rscript+=draw_multiple_profiles(breaks,avg_upstreams,main='Average Exon Profiles',xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim,v=0,cols=[],legends=legends)
    rscript+=draw_multiple_profiles(breaks,avg_upstreams,main='Average Intron Profiles',xlab='Relative Location (%)',ylab='Average Profile',ylim=ylim,v=0,cols=[],legends=legends)
    
    return rscript


def draw_single_profile_heatmap(profl,col='greenred(100)',main='',xlab='',ylab=''):
    
    #rscript=R.par(mar=[4,3.8,5,4],oma=[4,2,4,2])
    rscript=''
    # upstreams
    rscript+=heatmap_2(x=profl,Rowv=False,Colv=False,col=col,main='Heatmap near TSS',xlab='Relative Location (%)', ylab='Genes',key=True,density_info='density',keysize=0.7,labRow='')
    
    return rscript
   

def draw_CEAS(genome_p,ChIP_p,pval,bg_res=100,chip_res=600,prom=3000,biprom=5000,down=3000,gene_div=(3,5)):
    
    # got the time stamp
    timestamp=time.strftime("%H:%M:%S %a, %d %b %Y", time.localtime())
    
    # parameters
    rscript = '\n'
    rscript+=R.comment('%s' %timestamp)
    rscript+=R.comment('')
    rscript+=R.comment('ChIP annotation')
    rscript+=R.comment('')
    rscript+='\n'
    
    # PG1
    rscript+='\n'
    rscript+=R.comment('')
    rscript+=R.comment('Chromosomal Distribution')
    rscript+=R.comment('')
    rscript+= '\n'
    # set margin
    rscript+=R.par(mar=[4,4,5,3.8],oma=[4,2,4,2])

    #chromosome distribution
    rscript+=draw_chrom_barplot(genome_p,ChIP_p,pval)
    # put timestamp
    #rscript+=mtext(text=timestamp,side=1,line=2,font=3,outer=True,cex=0.8)
    
    # PG2
    # set the margin
    rscript+='\n'
    rscript+=R.comment('')
    rscript+=R.comment('Promoter,Bipromoter,Downstream, Gene and Regions of interest')
    rscript+=R.comment('')
    rscript+='\n'
    roi=True
    if genome_p['whole']['roi']==0.0: roi=False
    
    if roi: rscript+=R.par(mfrow=[5,1],mar=[4,4,5,3.8],oma=[4,2,4,2])
    else: rscript+=R.par(mfrow=[4,1],mar=[4,4,5,3.8],oma=[4,2,4,2])

    #promoter,bipromoter,downstream,gene
    rscript+=draw_promoter_barplot(genome_p,ChIP_p,pval,prom)
    rscript+=draw_bipromoter_barplot(genome_p,ChIP_p,pval,biprom)
    rscript+=draw_downstream_barplot(genome_p,ChIP_p,pval,down)
    rscript+=draw_gene_barplot(genome_p,ChIP_p,pval)
    if roi: rscript+=draw_roi_barplot(genome_p,ChIP_p,pval)
    # mark time stamp 
    #rscript+=mtext(text=timestamp,side=1,line=2,font=3,outer=True,cex=0.8)
        
    return rscript


def draw_wig_heatmap_overview(swig, ylim=None, chrominfo=None, colormap='GBR256'):
    """Draw wig heatmap overview"""
    
    # header print out in the R script
    rscript = '\n'
    rscript += R.comment('')
    rscript += R.comment('Wig heatmap overview')
    rscript += R.comment('')
    rscript += '\n'
    
    #    
    # get chromosome information
    #
    chroms=swig.get_chroms()
    n_chroms=len(chroms)
    if n_chroms==0: return rscript
    
    # parmeters
    # the height from the mid point (vertical) of the middle of the chromosome
    # when n_chroms is 6, 0.2 and n_chroms is 24, 0.4.
    hdelta = (0.4 - 0.2)/(24 - 6) * (n_chroms - 6) + 0.2 
    delta = 100 # delta added to the last point of wig
    widths = [1, 1]
    heights=[1, 5.5]
    mat= [[1, 0], [2, 2]]
    
    # get the largest chrom length and set the lower and upper length
    mins=[]
    maxs=[]
    if chrominfo:
        mins=[chrominfo[chrom]['start'][0] for chrom in chrominfo.get_chroms()]
        maxs=[chrominfo[chrom]['end'][-1] for chrom in chrominfo.get_chroms()] 
    else:
        mins=[swig[chrom][0][0] for chrom in swig.get_chroms()]
        maxs=[swig[chrom][0][-1] for chrom in swig.get_chroms()] 
    mins.sort()
    maxs.sort()
    maxchrom=[mins[0], maxs[-1]]
    
    #
    # get the lower and upper limits of the wig signal
    #   
    #if ylim is not given, decide the ylim as [-max of wig, max of wig]
    if not ylim: 
        percentage = 98
        left, right = swig.get_middle_percent(percentage)
        ylim=[-1*max(abs(left),abs(right)), max(abs(left),abs(right))]
        maxmax = max(ylim)
    
    # get the histogram of wig
    n_bins = 127
    histogram = swig.histogram(breaks=linspace(ylim[0], ylim[1], n_bins+1))
    
    #
    # draw the color bar and heatmaps
    #
        
    # put color code
    cmd = 'cm = graphics.%s' %colormap
    exec(cmd) 
    cmap = 'cmap'
    rscript += R.vector(cm, varname=cmap)
    
    # draw color bar 
    rscript += R.par(mar=[4,4,5,3.8],oma=[4, 2, 4, 2])
    rscript += R.layout(mat=mat,widths=widths,heights=heights)
    rscript += draw_colorbar(ylim, vertical=False, hist=histogram, linecol='cyan')
    rscript += R.plot(maxchrom, [1 - hdelta, n_chroms + hdelta], tp="n", col="", ylab="Chromosome", xlim=maxchrom, ylim=[1 - hdelta, n_chroms + hdelta], frame=False, xaxt="s", yaxt="n", main="Enrichment over Chromosomes", xlab="Chromosome Size (bp)")
    counter = n_chroms
    for chrom in chroms:
        bot = counter - hdelta
        top = counter + hdelta
        rscript += draw_a_chrom_heatmap(swig[chrom], chrom, ylim, [swig[chrom][0][0],swig[chrom][0][-1] + delta], maxchrom, centromere=None, bot=bot, top=top, line=True, linecol = 'cyan')
        rscript += R.mtext(chrom.rsplit('chr', 1)[-1], side=2, at=(top+bot)/2)
        counter -= 1
        
    return rscript


def draw_ChIP_over_genome_var_colors(ChIP, wig_size, n_peaks = 10000):
    """Draw heatmap bars of ChIP regions over the genome
    
    Parameters:
    ChIP: a Bed object of ChIP regions with scores
    wig_size: a dictionary contains the start and end positions of chromosomes
    """

    #
    # parameters for the graphics and mathematical process
    #
    
    # page dimension
    width = 8.5         # inch
    height = 11.5       # inch
    pixel = 1.0/72      # line width
    mar=[4,4,5,3.8] # marbin
    oma=[4, 2, 4, 2]    # outer margin
    rough_plot_width = 6.5
    
    # page segmentation
    widths = [1, 1]
    heights = [1, 5.5]
    mat = [[1, 0], [2, 2]]
     
    # percentage of scores to consider
    percentage = 98    
    
    # number of colors in the colormap
    n_colors = 16
    
    #
    # get info from chromosomes and ChIP regions
    #
    
    # get the chromosomes
    if wig_size != None:
        chroms = wig_size.keys()
        chroms = sort_chroms(chroms)
    else:
        chroms = ChIP.get_chroms()
    n_chroms = len(chroms)
    if n_chroms == 0: return rscript
    
    # bar height: adjusted according to the number of chromosomes
    hdelta = (0.4 - 0.2)/(24 - 6) * (n_chroms - 6) + 0.2 
    
    # get the maximum chromosome size
    if wig_size != None:
        minx = min([wig_size[chrom][0] for chrom in chroms])
        maxx = max([wig_size[chrom][1] for chrom in chroms]) 
    else:
        minx = min([ChIP[chrom]['start'][0] for chrom in chroms])
        maxx = max([ChIP[chrom]['end'][-1] for chrom in chroms])
    maxchrom = [minx, maxx]
    
    # if there are more than n_peaks, just select n_peaks by score
    reducedChIP = select_ChIP_by_n_peaks(ChIP, n_peaks = n_peaks, descend = True)
    # because the resulting pdf does not have infinite resolution, we need to get rid of redundant signals in the pdf.
    interval = (maxchrom[1] - maxchrom[0]) / (rough_plot_width * int(1/pixel))   # numerator: maximum chromosome length, denominator: how many bars can be drawn in rough_plot_width (6.5 in by default), interval is the bar width in bp
    reducedChIP = downsample_ChIP(reducedChIP, interval = interval)
    
    # get the histogram of bed
    left, right = get_percent_from_BED(reducedChIP, percentage, twoside=False)
    ylim=[left, right]
    maxmax = max([abs(ylim[0]), ylim[1]])       # get the largest boundary
    
    n_bins = 16
    if ylim[0] < 0: # if negative score, make symmetric
        ylim = [-1*maxmax, maxmax]
        breaks = linspace(-1*maxmax, maxmax, n_bins+1)
    else: # asymmetric
        ylim = [0, ylim[1]] 
        breaks = linspace(0, ylim[1]+0.2, n_bins+1) # add a small value to the maximum to include the max value
    histogram = get_histogram_BED_scores(reducedChIP, breaks)
    
    #
    # header print out in the R script
    #
    rscript = '\n'
    rscript+=R.comment('')
    rscript+=R.comment('ChIP regions over the genome')
    rscript+=R.comment('')
    rscript+= '\n'
    
    #
    # draw the color bar and heatmaps
    #
        
    # set the margins
    rscript += R.par(mar=mar,oma=oma)
    rscript += R.layout(mat=mat,widths=widths,heights=heights)
    
    # select color scheme
    if ylim[0] < 0:
        colormap = "GBR"
    else:
        colormap = "BR"
    cmd = 'cm = graphics.%s' %colormap
    exec(cmd) 
    
    # generate color map using R's colorRampPalette function
    rscript += R.colorRampPalette(cm, rval=colormap)
    cmap = 'cmap'   # color map name
    rscript += '%s <- %s(%s)\n' %(cmap, colormap, str(n_colors))
    rscript += draw_colorbar(ylim, vertical=False, cmap = cmap, hist=histogram, linecol='cyan')
    
    # iterate through chromosomes
    rscript += R.plot(maxchrom, [1 - hdelta, n_chroms + hdelta], tp="n", col="", ylab="Chromosome", xlim=maxchrom, ylim=[1 - hdelta, n_chroms + hdelta], frame=False, xaxt="s", yaxt="n", main="ChIP Regions over Chromosomes", xlab="Chromosome Size (bp)")
    count = n_chroms
    for chrom in chroms:
        bot = count - hdelta
        top = count + hdelta
        try:
            rscript += R.rectangles_with_heights_and_colors(reducedChIP[chrom]['start'], reducedChIP[chrom]['end'], reducedChIP[chrom]['score'], ylim, bot=bot, top=top, xaxt = "n", cmap = cmap)
            rscript += R.mtext(chrom.rsplit('chr', 1)[-1], side=2, at=(top+bot)/2)
        except KeyError:    # when ChIP does not have this chromosome
            rscript += R.rectangles_with_heights_and_colors(None, None, None, [minx, maxx], ylim, bot=bot, top=top, xaxt = "n")
            rscript += R.mtext(chrom.rsplit('chr', 1)[-1], side=2, at=(top+bot)/2)
        count -= 1
        
    return rscript
     
     
def draw_ChIP_over_genome_mono_col(ChIP, wig_size, n_peaks = 10000, barcol='#CC0000'):
    """Draw heatmap bars of ChIP regions over the genome
    
    Parameters:
    ChIP: a Bed object of ChIP regions with scores
    wig_size: a dictionary contains the start and end positions of chromosomes
    """

    #
    # parameters for the graphics and mathematical process
    #
    
    # page dimension
    width = 8.5         # inch
    height = 11.5       # inch
    pixel = 1.0/72      # line width
    mar=[4,4,5,3.8]  # margin
    oma=[4, 2, 4, 2]    # outer margin
    rough_plot_width = 7    # rough plot width - need improving
    
    # page segmentation
    widths = [1, 1]
    heights = [1, 5]
    mat = [[1, 0], [2, 2]]
     
    # percentage of scores to consider
    percentage = 98    
    
    #
    # get info from chromosomes and ChIP regions
    #
    
    # get the chromosomes
    if wig_size != None:
        chroms = wig_size.keys()
        chroms = sort_chroms(chroms)
    else:
        chroms = ChIP.get_chroms()
    n_chroms = len(chroms)
    if n_chroms == 0: return rscript
    
    # bar height: adjusted according to the number of chromosomes
    hdelta = (0.4 - 0.2)/(24 - 6) * (n_chroms - 6) + 0.2 
    
    # get the maximum chromosome size
    if wig_size != None:
        minx = min([wig_size[chrom][0] for chrom in chroms])
        maxx = max([wig_size[chrom][1] for chrom in chroms]) 
    else:
        minx = min([ChIP[chrom]['start'][0] for chrom in chroms])
        maxx = max([ChIP[chrom]['end'][-1] for chrom in chroms])
    maxchrom = [minx, maxx]
    
    # if there are more than n_peaks, just select n_peaks by score
    reducedChIP = select_ChIP_by_n_peaks(ChIP, n_peaks = n_peaks, descend = True)
    
    # because the resulting pdf does not have infinite resolution, we need to get rid of redundant signals in the pdf.
    interval = (maxchrom[1] - maxchrom[0]) / (rough_plot_width * int(1/pixel))   # numerator: maximum chromosome length, denominator: how many bars can be drawn in rough_plot_width (6.5 in by default), interval is the bar width in bp
    reducedChIP = downsample_ChIP(reducedChIP, interval = interval)
    
    # get the histogram of bed
    left, right = get_percent_from_BED(reducedChIP, percentage, twoside=False)
    ylim=[left, right]
    maxmax = max([abs(ylim[0]), ylim[1]])       # get the largest boundary
    
    n_bins = 16
    if ylim[0] < 0: # if negative score, make symmetric
        ylim = [-1*maxmax, maxmax]
        breaks = linspace(-1*maxmax-0.2, maxmax+0.2, n_bins+1)
    else: # asymmetric
        ylim = [0, ylim[1]] 
        breaks = linspace(0, ylim[1]+0.2, n_bins+1)
    histogram = get_histogram_BED_scores(reducedChIP, breaks)
    
    #
    # header print out in the R script
    #
    rscript = '\n'
    rscript+=R.comment('')
    rscript+=R.comment('ChIP regions over the genome')
    rscript+=R.comment('')
    rscript+= '\n'
    
    #
    # draw the color bar and heatmaps
    #
        
    # set the margins
    rscript += R.par(mar=mar,oma=oma)
    rscript += R.layout(mat=mat,widths=widths,heights=heights)
        
    # draw histogram
    rscript += draw_histogram_bar(histogram, ylim, main = "Distribution of Peak Heights", linecol='cyan', bgcol='black', cex=0.9)
    
    # iterate through chromosomes
    rscript += R.plot(maxchrom, [1 - hdelta, n_chroms + hdelta], tp="n", col="", ylab="Chromosome", xlim=maxchrom, ylim=[1 - hdelta, n_chroms + hdelta], frame=False, xaxt="s", yaxt="n", main="ChIP Regions (Peaks) over Chromosomes", xlab="Chromosome Size (bp)")
    count = n_chroms
    for chrom in chroms:
        bot = count - hdelta
        top = count + hdelta
        try:
            rscript += R.rectangles_with_heights(reducedChIP[chrom]['start'], reducedChIP[chrom]['end'], reducedChIP[chrom]['score'], ylim, bot=bot, top=top, xaxt = "n", col = [barcol])
            rscript += R.mtext(chrom.rsplit('chr', 1)[-1], side=2, at=(top+bot)/2)
        except KeyError:    # when ChIP does not have this chromosome
            rscript += R.rectangles_with_heights(None, None, None, ylim, bot=bot, top=top, xaxt = "n")
            rscript += R.mtext(chrom.rsplit('chr', 1)[-1], side=2, at=(top+bot)/2)
        count -= 1
        
    return rscript
   

def draw_wig_overview(wig,plots_per_pg=4):
    """Draw wig overview"""
    
    # parameters
    rscript=R.comment('')
    rscript+=R.comment('Wig Overview over the Genome')
    rscript+=R.comment('')
    
    # count how many chromosomes and draw four chromosomes per page by default
    chroms=wig.get_chroms()
    if len(chroms)==0: return rscript
    
    chromcnt=1
    for chrom in chroms:
        if chromcnt%plots_per_pg==1: rscript+=R.par(mar=[4,4,5,3.8],oma=[4,2,4,2])
        rscript+=R.barplot(wig[chrom][1],names=wig[chrom][0],col=["red"],border=False,space=0,main=chrom,ylab="Profile",xlab="Chromosome (nt)")
        chromcnt+=1
           
    return rscript
    
    
    

   


   