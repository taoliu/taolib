#!/usr/bin/env python
# Time-stamp: <2009-03-03 16:58:17 Tao Liu>

import os
import sys
# ------------------------------------
# Main function
# ------------------------------------
def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Convert fastq file to fasta file\nneed 2 paras: %s <in.fq> <out.fa>\n" % sys.argv[0])
        sys.exit(1)
    infhd  = open(sys.argv[1],"r")
    outfhd = open(sys.argv[2],"w")
    readname = None
    for l in infhd:
        if readname:
            outfhd.write(">"+readname+l)
            readname = None
        else:
            if l.startswith("@"):
                readname = l[1:]
    infhd.close()
    outfhd.close()

if __name__ == '__main__':
    main()
