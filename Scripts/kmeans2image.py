#!/usr/bin/env python
# Time-stamp: <2009-04-14 14:07:21 Tao Liu>

import os
import sys
import re
from PIL import Image, ImageDraw

# ------------------------------------
# Main function
# ------------------------------------

help_message = """
Draw the K-means clustering result.
need 6 parameter: %s <kmeans_file> <lim> <x_points> <y_points> <x_ext> <y_ext>

kmeans_file : tab-delimited plain text file. First column is cluster number by k-means, and following columns are data columns.
lim         : data value limit
x_points    : number of data value columns
y_points    : number of rows
x_ext       : pixels extended in x-axis
y_ext       : pixels extended in y-axis
""" % sys.argv[0]

def main():
    if len(sys.argv) < 7:
        sys.stderr.write(help_message)
        sys.exit(1)

    fhd = open (sys.argv[1])
    lim = int(sys.argv[2])
    x_points = int(sys.argv[3])
    y_points = int(sys.argv[4])
    x_ext = int(sys.argv[5])
    y_ext = int(sys.argv[6])

    a = Image.new("RGB",(x_points*x_ext,y_points*y_ext),"white")
    
    d = ImageDraw.Draw(a)

    y = 0
    for i in fhd:
        y += 1
        i.strip()
        if not re.search("^\d+",i):
            continue
        values = map(float,i.split())
        x = 0
        cl = values[0]
        for v in values[1:]:
            x += 1
            c = "hsl(%d,100%%,%d%%)" % (cl*70,min(1,v/lim)*90.0)
            d.rectangle([(int(x*x_ext),int(y*y_ext)),(int((x+1)*x_ext),int((y+1)*y_ext))],outline=c,fill=c)

    a.save(sys.argv[1]+".png")
    print "check %s!" % (sys.argv[1]+".png")


if __name__ == '__main__':
    main()

