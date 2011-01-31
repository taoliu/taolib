from taolib.CoreLib.FeatIO import PeakIO
from csv import DictReader

def parse_pMA2C_xls (fhd):
    """Parse a tab-delimited xls file from python version MA2C.

    fhd is a file handler for python version of MA2C xls result.

    Return a taolib.CoreLib.FeatIO.PeakIO object.
    """
    reader = DictReader(Filter(fhd),dialect='excel-tab')
    peaks = PeakIO()
    for row in reader:
        peaks.add(row['Chr'],int(row['Start']),int(row['End']),
                  int(row['Summit']),float(row['MA2Cscore']),
                  0,float(row['-10*log10(PValue)']),0,float(row['FDR']))
    return peaks


def parse_MAT_xls (fhd):
    """Parse a tab-delimited xls file from MAT.

    fhd is a file handler for MAT xls result.

    Return a taolib.CoreLib.FeatIO.PeakIO object.
    """
    peaks = PeakIO()
    for row in fhd:
        if row.startswith("Chr") or row.startswith("#") or row.startswith("track"):
            continue
        else:
            f = row.split()
            peaks.add(f[0],int(f[1]),int(f[2]),
                      int(f[8]),float(f[5]),
                      0,float(f[4]),0,float(f[7]))
            break
    
    for row in fhd:
        f = row.split()
        peaks.add(f[0],int(f[1]),int(f[2]),
                  int(f[8]),float(f[5]),
                  0,float(f[4]),0,float(f[7]))

    return peaks

def parse_MACS_xls (fhd):
    """Parse a tab-delimited xls file from MACS.

    fhd is a file handler for MACS xls result.

    Return a taolib.CoreLib.FeatIO.PeakIO object.
    """
    reader = DictReader(Filter(fhd),dialect='excel-tab')
    peaks = PeakIO()
    for row in reader:
        peaks.add(row['chr'],int(row['start']),int(row['end']),
                  int(row['summit']),0,
                  int(row['tags']),float(row['-10*log10(pvalue)']),
                  float(row['fold_enrichment']),float(row['FDR(%)']))
    return peaks


class Filter:
    def __init__(self, fhd):
        self.f = fhd
    
    def _nextLine(self):
        
        result = None
        
        while not result:
            result = self.f.next().strip()
            commentPos = result.find('#')
            if commentPos != - 1:
                result = result[:commentPos]

        return result
    
    def __iter__(self):
        return self
    
    def next(self):
        return self._nextLine()
