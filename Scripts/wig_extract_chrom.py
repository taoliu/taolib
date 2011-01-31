#!/usr/bin/env python
# Time-stamp: <2009-08-18 23:13:45 Tao Liu>

import os
import sys

from taolib.CoreLib.Parser import WiggleIO
from taolib.CoreLib.FeatIO import WigTrackI
# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 4:
        sys.stderr.write("Extract data for certain chromosome from a wiggle file.\n")
        sys.stderr.write("need 3 paras: %s <chr> <wig> <newwig>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)

    chrom = sys.argv[1]
    wigfhd = open(sys.argv[2])
    wigtrack = WiggleIO.WiggleIO(wigfhd).build_wigtrack()
    wig_chr = wigtrack.get_data_by_chr(chrom)
    if not wig_chr:
        sys.stderr.write("No data for chromosome %s!\n" % chrom)
        sys.exit(1)

    newwigtrack = WigTrackI()
    newwigtrack.span = wigtrack.span

    (wig_chr_p,wig_chr_s) = wig_chr

    for i in range(len(wig_chr_p)):
        newwigtrack.add_loc(chrom,wig_chr_p[i],wig_chr_s[i])

    newwigfile = sys.argv[3]
    newwigfhd = open(newwigfile,"w")
    newwigtrack.write_wig(newwigfhd,name="for chromosome %s" % chrom)
    newwigfhd.close()
        

if __name__ == '__main__':
    main()
