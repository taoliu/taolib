#!/Library/Frameworks/Python.framework/Versions/Current/bin/python2.6
# Time-stamp: <2008-10-27 15:51:34 Tao Liu>

import re
import os
import cgi
import cgitb; cgitb.enable()
import md5
import time
from Cistrome.CoreLib.DB import *

def update2NDD( nd,form ):
    l = ["Organism","Assembly","Factor","Cell","Platform",
         "Lab","Condition","Type","Status","PubMED","GEO",
         "Raw","XLS","Wiggle","Comment"]
    for i in l:
        v = form.getvalue(i)
        if v:
            nd.set(i,v)
    return nd


print "content-type: text/html\n"     # HTML is following

HTML_HEAD = '''<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <title>Cistrome - Main page</title>
  <link rel="author" href="http://liulab.dfci.harvard.edu/~taoliu/web/AboutMe.html" />
  <link rel="stylesheet" type="text/css" charset="utf-8" media="all" href="/cistrome/css/common.css" />
<link rel="stylesheet" type="text/css" charset="utf-8" media="screen" href="/cistrome/css/screen.css" />
<link rel="stylesheet" type="text/css" charset="utf-8" media="print" href="/cistrome/css/print.css" />
  <meta name="description" content="Cistrome website by Liu lab and
  Brown lab" />
  <meta name="keywords" content="chip-seq, chip-chip, solexa, solid, biology, computer" />
  <meta name="generator" content="Python" />
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
</head>

<body>
<div class="menu">
  <div class="menuitem">
    <a title="About this site" href="/cistrome/About.html">About</a>
  </div>
  <div class="menuitem">
    <a title="The tutorial" href="/cistrome/Tutorial.html">Tutorial</a>
  </div>
  <div class="menuitem">
    <a title="The Services" href="/cistrome/Services.html">Services</a>
  </div>
  <div class="menuitem">
    <a title="The Frequently Asked Questions" href="/cistrome/FAQ.html">FAQ</a>
  </div>
  <div class="menuitem">
    <a title="The Changelog" href="/cistrome/ChangeLog.html">ChangeLog</a>
  </div>
  <div class="menuitem">
    <a title="The TODO list" href="/cistrome/Todo.html">TODOs</a>
  </div>
</div><!-- menu ends here -->

<h1 id="top">
<a href="index.html"
style="text-decoration: none">
<img src="/cistrome/pic/Cistrome_logo.png" class="logo"
alt="Cistrome logo" />
</a>::
Search NDD
</h1>

<hr />
'''
HTML_END = '''
    <div class="navfoot">
      <hr />
      <table width="100%" border="0" summary="Footer navigation">
        <col width="90%" /><col width="5%" /><col width="5%" />
        <tr>
          <td align="left">
            [ <a href="/cistrome/index.html">Cistrome</a>| <a href="http://liulab.dfci.harvard.edu/"> Liu Lab</a> ]
          </td>
          <td>
            <a href="http://validator.w3.org/check?uri=referer"
               style="text-decoration: none;">
              <img src="http://www.w3.org/Icons/valid-xhtml10"
                   alt="Valid XHTML 1.0!" height="31" width="88" />
            </a>
          </td>
          <td>
            <a href="http://jigsaw.w3.org/css-validator/check/referer"
               style="text-decoration: none;">
              <img style="border:0;width:88px;height:31px"
                   src="http://jigsaw.w3.org/css-validator/images/vcss"
                   alt="Valid CSS!" />
            </a>
          </td>
        </tr>
      </table>
    </div>
  </body>
</html>'''

form = cgi.FieldStorage()

#dbname = form.getvalue('db')

#HTML_HEAD = re.sub("value=\"%s\"" % (dbname),"value=\"%s\" selected=\"selected\"" % (dbname),HTML_HEAD)

print HTML_HEAD

queryid = form.getvalue("Index")
passcode = form.getvalue('passcode')
md5key = form.getvalue('md5key')

db = DataCollection.DBDataCollection(host="localhost",user="root",db="dataCollection")

nd = db.select(queryid)

if passcode == 'cistrome' or md5key == nd.md5key:
    if db.update_NDDdata(queryid,update2NDD(nd,form)):
        print "Success!"
        print "<a href=\"/cgi-bin/search_by_id.py?id=%s\">go back</a>" % queryid
    else:
        print "Failed!"
else:
    print "<b>Passcode or md5key incorrect! You must be the owner or admin to delete an entry from database!</b><br />"

print HTML_END

