# Time-stamp: <2008-10-20 15:52:52 Tao Liu>

"""Module Description: Parser for NDD file

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
import re
# ------------------------------------
# constants
# ------------------------------------
NDDFIELDS = ["Index","Factor","Organism","Assembly",
             "Cell","Platform","Status","Raw","XLS",
             "Wiggle","Lab","Paper","Condition","GEO",
             "PubMED","Type","Comment","Submit"]
    
# ------------------------------------
# Misc functions
# ------------------------------------
def parse (fhd):
    """Parse a NDD file object.

    fhd must be a file object.
    """
    #assert type(fhd)==file
    data = NDDdata()
    for l in fhd:
        if l.startswith("#"):
            continue
        #print l
        if l.find('=') == -1:
            continue
        (name,value) = l.strip().split('=')
        name = name.strip()
        value = value.strip()
        #print name,value
        data.set(name,value)
    return data
    
# ------------------------------------
# Classes
# ------------------------------------
class NDDdata:
    """Store NDD data.
    
    """
    def __init__ (self):
        self.__d = ['NA']*19
        self.timestamp = ''
        self.md5key = ''

    def get ( self, name):
        try:
            index = NDDFIELDS.index(name)
        except ValueError:
            return False
        return self.__d[index]

    def __td_input_id ( self, name):
        try:
            index = NDDFIELDS.index(name)
        except ValueError:
            return False
        t =u"""
    <tr>
        <td bgcolor=\"d0d0d0\" colspan=2>
	  %s:
	  <input name=\"%s\" type=\"text\" value=\"%s\" size=%d /><br />
	</td>
    <tr/>
    """ % (name,name,self.__d[index],100-len(name))
        return t

    def __td_id ( self, name):
        try:
            index = NDDFIELDS.index(name)
        except ValueError:
            return False
        t =u"""
    <tr>
        <td bgcolor=\"d0d0d0\" colspan=2>
	  <b>%s</b>:
	  <input name=\"%s\" type=\"text\" value=\"%s\" readonly size=%d /><br />
	</td>
    <tr/>
    """ % (name,name,self.__d[index],100-len(name))
        return t


    def __str__ (self):
        #t = "\n".join([x for x in self.__d if x])
        #return t+"\n"
        return ">>"+str(self.__d)+"<<"
        
    def set ( self, name, value):
        try:
            index = NDDFIELDS.index(name)
        except ValueError:
            return False
        self.__d[index]=value
        return True

    def htmlize (self):
        t = u"""
        <b>id</b> : %s <a href=\"/cgi-bin/delete_by_id.py?id=%s\">[delete]</a>  <a href=\"/cgi-bin/update_by_id.py?id=%s\">[update]</a><br />
        <b> Submitted</b> : %s <b>by</b> <a href=\"/cgi-bin/search_by_submit.py?submit=%s\">%s</a> <br />
        <b>Paper Title</b> : <a href=\"/cgi-bin/search_by_paper.py?paper=%s\">%s</a><br />
        <b>Organism</b> : <a href=\"/cgi-bin/search_by_organism.py?organism=%s\">%s</a><br />
        <b>Factor</b> : %s<br />        
        <b>Cell</b> : <a href=\"/cgi-bin/search_by_cell.py?cell='%s'\">%s</a><br />
        <b>Platform</b> : %s<br />
        <a href=\"/cgi-bin/detail_by_id.py?id=%s\">[detail]</a></b>
        """ % (self.get('Index'),self.get('Index'),self.get('Index'),
               self.timestamp,self.get('Submit'),self.get('Submit'),
               self.get('Paper'),self.get('Paper'),
               self.get('Organism'),self.get('Organism'),
               self.__Factor2link(),
               self.get('Cell'),self.get('Cell'),
               self.get('Platform'),
               self.get('Index')
               )
        t= t.encode('utf-8')
        return t

    def htmlize_update (self):
        t =u"""
<form name=\"arcform\" action=\"/cgi-bin/update_by_id_confirmed.py\"
method=\"query\">
<table border=1>
"""
        end = u"""
  <tr>
        <td bgcolor=\"d0d0d0\">
	  Admin passcode:
	  <input name=\"passcode\" type=\"password\" /><br />
	</td>
	<td bgcolor=\"d0d0d0\">
	  MD5:
	  <input name=\"md5key\" type=\"text\" /><br />
	</td>
  </tr>
  <tr>
	<td bgcolor=\"d0d0d0\" colspan=2>
	  <input name=\"process\" type=\"submit\" value=\"update\" />
	</td>
  </tr>
</table>
</form>
"""
        t += self.__td_id("Index")
        t += u"""
    <tr>
        <td bgcolor=\"d0d0d0\" colspan=2>
	  <b>Time</b>:
	  <input name=\"Time\" type=\"text\" value=\"%s\" readonly size=%d /><br />
	</td>
    <tr/>
    """ % (self.timestamp,100-len("Time"))
        t += self.__td_id("Submit")
        t += self.__td_id("Paper")
        t += self.__td_input_id("Organism")
        t += self.__td_input_id("Assembly")
        t += self.__td_input_id("Factor")
        t += self.__td_input_id("Cell")
        t += self.__td_input_id("Platform")
        t += self.__td_input_id("Lab")
        t += self.__td_input_id("Condition")
        t += self.__td_input_id("Type")
        t += self.__td_input_id("Status")
        t += self.__td_input_id("PubMED")
        t += self.__td_input_id("GEO")
        t += self.__td_input_id("Raw")
        t += self.__td_input_id("XLS")
        t += self.__td_input_id("Wiggle")
        t += self.__td_input_id("Comment")
        t += end
        t= t.encode('utf-8')
        return t


    def htmlize_detail (self):
        t = u"""
        <b>id</b> : %s <a href=\"/cgi-bin/delete_by_id.py?id=%s\">[delete]</a>  <a href=\"/cgi-bin/update_by_id.py?id=%s\">[update]</a><br />
        <b> Submitted</b> : %s <b>by</b> <a href=\"/cgi-bin/search_by_submit.py?submit=%s\">%s</a> <br />
        <b>Paper Title</b> : <a href=\"/cgi-bin/search_by_paper.py?paper=%s\">%s</a><br />
        <b>PubMED</b> : <a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&cmd=search&term=%s[pmid]\">%s</a><br />        
        <b>Organism</b> : <a href=\"/cgi-bin/search_by_organism.py?organism=%s\">%s</a><br />
        <b>Assembly</b> : %s<br />
        <b>Factor</b> : %s<br />
        <b>Cell</b> : <a href=\"/cgi-bin/search_by_cell.py?cell='%s'\">%s</a><br />
        <b>Platform</b> : %s<br />
        <b>Lab</b> : %s<br />
        <b>Condition</b> : %s<br />
        <b>Type</b> : %s<br />
        <b>Status</b> : %s<br />
        <b>GEO</b> : %s<br />
        <b>Raw data</b> : %s<br />
        <b>Peak XLS</b> : %s<br />
        <b>Score Wiggle</b> : %s<br />
        
        """ % (self.get('Index'),self.get('Index'),self.get('Index'),
               self.timestamp,self.get('Submit'),self.get('Submit'),
               self.get('Paper'),self.get('Paper'),
               self.get('PubMED'),self.get('PubMED'),               
               self.get('Organism'),self.get('Organism'),
               self.get('Assembly'),
               self.__Factor2link(),
               self.get('Cell'),self.get('Cell'),
               self.get('Platform'),
               self.get('Lab'),
               self.get('Condition'),
               self.get('Type'),
               self.get('Status'),
               self.__GEO2link(),
               self.get('Raw'),
               self.get('XLS'),
               self.get('Wiggle')
               )
        t= t.encode('utf-8')
        return t

    def __GEO2link ( self ):
        s = self.get("GEO")
        return re.sub(r"(G\w{2}\d+)",r'<a href="http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=\1">\1</a>',s)

    def __Factor2link ( self ):
        f = self.get("Factor")
        fs = ",".join(map(lambda x:"<a href=\"/cgi-bin/search_by_factor.py?factor="+x+"\">"+x+"</a>",f.split(",")))
        return fs
