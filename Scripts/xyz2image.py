#!/usr/bin/env python
# Time-stamp: <2009-04-08 16:54:22 Tao Liu>

import os
import sys
import re
from PIL import Image, ImageDraw

# ------------------------------------
# Main function
# ------------------------------------

def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Restore the microarray image.\n need 2 parameter: %s <X-Y-Value-file> <x_ext> <y_ext>\n X-Y-Value-file is a tab-delimited file with 1st column as X-corrdinates, 2nd column as Y, and 3rd column as log-ratio value\n" % sys.argv[0])
        sys.exit(1)

    fhd = open (sys.argv[1])
    x_ext = int(sys.argv[2])
    y_ext = int(sys.argv[3])

    a = Image.new("RGB",(1100,4400),"white")
    
    d = ImageDraw.Draw(a)

    for i in fhd:
        i.strip()
        if not re.search("^\d+",i):
            continue
        (x,y,z) = i.split()
        x=int(x)
        y=int(y)
        z=int(float(z)*10)
        if z>=0:
            c = "hsl(%d%%,100%%,90%%)" % (max(0,100-z))
        else:
            c = "hsl(%d%%,100%%,90%%)" % (max(0,100+z))            
        d.rectangle([(int(x*x_ext),int(y*y_ext)),(int((x+1)*x_ext),int((y+1)*y_ext))],outline=c,fill=c)

    a.save(sys.argv[1]+".tiff")
    print "check %s!" % (sys.argv[1]+".tiff")


if __name__ == '__main__':
    main()

