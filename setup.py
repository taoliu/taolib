#!/usr/bin/env python

import os
import sys
from distutils.core import setup, Extension

def main():
    if not float(sys.version[:3])>=2.4:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4! python 2.6.2 is recommended!\n")
        sys.exit(1)
    setup(name="taolib",
          version="1.0",
          description="Tao's libraries",
          author='Tao (Foo) Liu',
          author_email='vladimir.liu@gmail.com',
          url='http://vladimirliu.com/~taoliu/',
          package_dir={'taolib' : '.'},
          packages=['taolib','taolib.CoreLib',
                    'taolib.CoreLib.DB','taolib.CoreLib.FeatIO',
                    'taolib.CoreLib.BasicStat','taolib.CoreLib.WWW',
                    'taolib.CoreLib.Parser','taolib.CoreLib.SeqIO',
                    'taolib.CoreLib.BinKeeper','taolib.CoreLib.Algorithm',
                    'taolib.Assoc',
                    'taolib.ExtApp',
                    'taolib.Motif',
#                     'taolib.IntegrativeBioinformatics',
#                     'taolib.IntegrativeBioinformatics.elements',
#                     'taolib.IntegrativeBioinformatics.networks',
#                     'taolib.IntegrativeBioinformatics.algos',
#                     'taolib.IntegrativeBioinformatics.features',
#                     'taolib.IntegrativeBioinformatics.links',
#                     'taolib.IntegrativeBioinformatics.apache',
                    ],

          scripts=['Scripts/motif_enrich.py',
                   'Scripts/qc_chIP_peak.py',
                   'Scripts/qc_chIP_whole.py',
                   'Scripts/count_probes_in_peaks.py',
                   'Scripts/count_probes_in_ranges.py',
                   'Scripts/xyz2image.py',
                   'Scripts/refine_peak.py',
                   'Scripts/fq2fa.py',
                   'Scripts/wiggle_reformat.py',
                   'Scripts/wig_correlation.py',
                   'Scripts/wig_correlation_in_bed_file.py',
                   'Scripts/conservation_plot.py',
                   'Scripts/wig_extract_chrom.py',
                   'Scripts/bed_extract_chrom.py',
                   'Scripts/wig_call_peaks.py',                   
                   'Scripts/wig_call_peaks2.py',
                   'Scripts/naive_call_peaks.py',
                   'Scripts/wig2bedGraphBins.py',
                   'Scripts/bed_correlation.py',                   
                   'Scripts/ce_histone_matrix.py',
                   'Scripts/rand_pos',
                   'Scripts/draw_BED.py',
                   'Scripts/norm.py',
#                   'Scripts/hmm_conception.py',
                   ],

          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: Database',
            ],
          requires=['MySQL_python','PIL']
          )

if __name__ == '__main__':
    main()

