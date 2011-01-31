#!/Library/Frameworks/Python.framework/Versions/Current/bin/python2.6
# Time-stamp: <2008-10-27 15:47:11 Tao Liu>

import re
import os
import cgi
import cgitb; cgitb.enable()
import md5
import time
from Cistrome.CoreLib.DB import *

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

#HTML_HEAD = re.sub("value=\"%s\"" % (dbname),"value=\"%s\" selected=\"selected\"" % (dbname),HTML_HEAD)

print HTML_HEAD

db = DataCollection.DBDataCollection(host="localhost",user="root",db="dataCollection")

papers = db.execute('select Paperid from dataset group by Paperid').fetchall()

print 'Here are all the papers in Cistrome database: (total <b>%s</b>)<br />' % (str(len(papers)))
print '[number of entries] title of paper'

print '<table>'

for p in papers:
    paper = p[0]
    print "<tr><td>"
    n = db.execute('select COUNT(*) from dataset where Paperid="%s"' % (paper)).fetchall()[0][0]
    print "[%s]" % (str(n))
    print "<a href='/cgi-bin/search_by_paper.py?paper=%s'>%s</a>" % (paper,paper)
    print "</td></tr>"

#sqlsentence = 'select id from dataset where Paperid="%s"' % (paper)

#ids = [x[0] for x in db.execute(sqlsentence)]##

#for i in ids:
#    print "<hr />"
#    nd = db.select(i)
#    print nd.htmlize()

print HTML_END
