#!/usr/bin/env python
# Time-stamp: <2011-02-04 16:01:52 Tao Liu>
# Split the wig files into individual ones for each chromosome


import os
import sys

from taolib.CoreLib.Parser import WiggleIO
from taolib.CoreLib.FeatIO import WigTrackI

# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Extract data for all chromosomes from a wiggle file.\n")
        sys.stderr.write("need 2 paras: %s <wig> <output_prefix>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)

    wigfhd = open(sys.argv[1])
    prefix = sys.argv[2]
    
    wigtrack = WiggleIO.WiggleIO(wigfhd).build_wigtrack()
    wigtrack.sort()
    
    for chrom in wigtrack.get_chr_names():
        wig_chr = wigtrack.get_data_by_chr(chrom)
        newwigtrack = WigTrackI()
        newwigtrack.span = wigtrack.span
	
        (wig_chr_p,wig_chr_s) = wig_chr

        for i in range(len(wig_chr_p)):
            newwigtrack.add_loc(chrom,wig_chr_p[i],wig_chr_s[i])

        newwigfile = sys.argv[2]+"."+chrom+".wig"
        newwigfhd = open(newwigfile,"w")
        newwigtrack.write_wig(newwigfhd,name="for chromosome %s" % chrom)
        newwigfhd.close()
        

if __name__ == '__main__':
    main()
