#!/Library/Frameworks/Python.framework/Versions/Current/bin/python2.6
# Time-stamp: <2008-10-27 16:50:17 Tao Liu>

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

form = cgi.FieldStorage()

#dbname = form.getvalue('db')

#HTML_HEAD = re.sub("value=\"%s\"" % (dbname),"value=\"%s\" selected=\"selected\"" % (dbname),HTML_HEAD)

print HTML_HEAD

queryid = form.getvalue('id')

db = DataCollection.DBDataCollection(host="localhost",user="root",db="dataCollection")

try:
    nd = db.select(queryid)
except:
    print "No such entry with id %s!<br />" % (str(queryid))    
else:
    print "following entry is about to be removed!<hr />"
    print nd.htmlize()
    print "<hr />"
    
    print "Please input md5sum key for this entry or administrator passcode<br />"
    print "<hr />"
    print """
<form name=\"arcform\" action=\"/cgi-bin/delete_by_id_confirmed.py\"
method=\"query\">
<table border>
  <tr>
        <td bgcolor=\"d0d0d0\">
	  Admin passcode:
	  <input name=\"passcode\" type=\"password\" /><br />
	</td>
	<td bgcolor=\"d0d0d0\">
	  MD5:
	  <input name=\"md5key\" type=\"text\" /><br />
	</td>
	<td bgcolor=\"d0d0d0\">
	  <input name=\"id\" type=\"hidden\" value=%s />
	</td>
	<td bgcolor=\"d0d0d0\">
	  <input name=\"process\" type=\"submit\" value=\"delete\"/>
	</td>
  </tr>
</table>
</form>
""" % (queryid)

print HTML_END
