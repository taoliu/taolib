#!/Library/Frameworks/Python.framework/Versions/Current/bin/python2.6
# Time-stamp: <2008-10-27 15:39:40 Tao Liu>

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
Submit NDD
</h1>

<form name="arcform" action="/cgi-bin/submit.py"
method="post" enctype="multipart/form-data">
<table border>
  <tr>
        <td bgcolor="d0d0d0">
	  Submit in
          <select name="db">
            <option value="RawChIP">Raw ChIP DB</option>
            <option value="Motif">Motif DB</option>
            <option value="GEO">GEO DB</option>
          </select>
        </td>
	<td bgcolor="d0d0d0">
	  Submit your NDD file:
	  <input name="ifile" type="file" /><br />
	</td>
	<td bgcolor="d0d0d0">
	  <input name="process" type="submit" />
          <a href="/cistrome/HelpSubmit.html">Help</a>
	</td>
  </tr>
</table>
</form>

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

dbname = form.getvalue('db')

HTML_HEAD = re.sub("value=\"%s\"" % (dbname),"value=\"%s\" selected=\"selected\"" % (dbname),HTML_HEAD)

print HTML_HEAD

nddfile = form['ifile']

ndddata = NDD.parse(nddfile.file)

db = DataCollection.DBDataCollection(host="localhost",user="root",db="dataCollection")

# time stamp for this data
tstamp = time.asctime() + ' '+ time.tzname[0]
ndddata.timestamp = tstamp

# md5 key for this data
md5key = md5.new(tstamp+str(ndddata)).hexdigest()
ndddata.md5key = md5key

# index number for this data
theindex = db.add_NDDdata(ndddata)

print "Thank you %s!<br />" % (ndddata.get('Submit'))
print "Your submission is completed!<br />"
print "The md5 key for this entry is '<font color=\"red\" ><i>%s</i></font>'. The key will be required when you need to modify or delete your entry!<br />" % (md5key)

yoursubmitteddata = db.execute('select COUNT(*) from dataset where Submitid="%s"' % (ndddata.get('Submit'))).fetchall()[0][0]

yoursubmittedpub = db.execute('select COUNT(DISTINCT PubMED) from dataset where Submitid="%s" group by PubMED' % (ndddata.get('Submit'))).fetchall()[0][0]

# link to your submission
link2submit = "/cgi-bin/search_by_submit.py?submit=%s" % (ndddata.get('Submit'))
link2submit_pubmed = "/cgi-bin/search_by_submit.py?submit=%s" % (ndddata.get('Submit'))
link2index = "/cgi-bin/search_by_id.py?id=%s" % (theindex)

print "All the datasets you submitted: %d [<a href=\"%s\">link</a>]<br />" % (yoursubmitteddata,link2submit)
print "All Publications involved: %d [<a href=\"%s\">link</a>]<br />" % (yoursubmittedpub,link2submit_pubmed)
print "<hr />"
print "The entry you submitted is #%d [<a href=\"%s\">link</a>]:<br /><br />" % (theindex,link2index)

text = db.execute ('select id,Organismid,Cellid,Factorid,Paperid,Statusid,Timestamp from dataset where id=%d' % (theindex)).fetchall()[0]
print "<b>Time submitted: %s</b><br /><br />" % (text[6])
print "<b>id</b>: %s<br />" % (text[0])
print "<b>Organism</b>: %s<br />" % (text[1])
print "<b>Cell</b>: %s<br />" % (text[2])
print "<b>Factor</b>: %s<br />" % (text[3])
print "<b>Paper</b>: %s<br />" % (text[4])
print "<b>Status</b>: %s<br />" % (text[5])
#nddseldata = db.select(theindex)

print HTML_END
