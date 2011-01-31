#!/usr/bin/env python

import os
import sys
from distutils.core import setup, Extension
#from Cython.Distutils import build_ext

def main():
    if not float(sys.version[:3])>=2.4:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4! python 2.6.2 is recommended!\n")
        sys.exit(1)
#    ext_modules = [Extension("taolib.CoreLib.FeatIO.WigTrackII", ["CoreLib/FeatIO/WigTrackII.pyx"])]
    setup(name="taolib",
          version="1.0",
          description="Tao's libraries",
          author='Tao (Foo) Liu',
          author_email='taoliu@jimmy.harvard.edu',
          url='http://cistrome.dfci.harvard.edu/~taoliu/',
#          cmdclass = {'build_ext': build_ext},
#          ext_modules = ext_modules,
          package_dir={'taolib' : '.'},
          packages=['taolib','taolib.CoreLib',
                    'taolib.CistromeDB','taolib.CistromeDB.datacollection',
                    'taolib.CoreLib.DB','taolib.CoreLib.FeatIO',
                    'taolib.CoreLib.BasicStat','taolib.CoreLib.WWW',
                    'taolib.CoreLib.Parser','taolib.CoreLib.SeqIO',
                    'taolib.CoreLib.BinKeeper','taolib.CoreLib.Algorithm',
                    'taolib.Assoc',
                    'taolib.ExtApp','taolib.Motif',
#                     'taolib.IntegrativeBioinformatics',
#                     'taolib.IntegrativeBioinformatics.elements',
#                     'taolib.IntegrativeBioinformatics.networks',
#                     'taolib.IntegrativeBioinformatics.algos',
#                     'taolib.IntegrativeBioinformatics.features',
#                     'taolib.IntegrativeBioinformatics.links',
#                     'taolib.IntegrativeBioinformatics.apache',
                    ],

          scripts=['Scripts/motif_enrich.py','Scripts/qc_chIP_peak.py','Scripts/qc_chIP_whole.py',
                   'Scripts/count_probes_in_peaks.py','Scripts/count_probes_in_ranges.py',
                   'Scripts/xyz2image.py','Scripts/refinePeak.py','Scripts/fq2fa.py',
                   'Scripts/wiggle_reformat.py','Scripts/wig_correlation.py',
                   'Scripts/wig_correlation_in_bed_file.py','Scripts/conservation_plot.py',
                   'Scripts/wig_extract_chrom.py','Scripts/bed_extract_chrom.py','Scripts/wig_call_peaks.py',
                   'Scripts/naive_call_peaks.py',
                   'Scripts/wig2bedGraphBins.py',
                   'Scripts/bed_correlation.py',                   
                   'Scripts/bidirection_promoter.py',
                   'Scripts/heatmap.py',
                   'Scripts/wig_call_peaks2',
                   'Scripts/ceHistoneMatrix.py',
                   'Scripts/randPos',
                   'Scripts/drawBED.py',
                   'Scripts/sitepro2heatmap.py',
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

