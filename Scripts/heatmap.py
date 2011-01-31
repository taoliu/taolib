#!/usr/bin/env python

"""Draw heatmaps from aligned WIG profiles."""

import sys
import optparse

import numpy
import scipy
import scipy.cluster
import pylab



class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class WIGError(Error):
    """Exception raised for errors in WIG file format."""
    pass

def build_binKeeper (file):
    """Use this function to return a dictionary of BinKeeper
    objects.
    NOTE: this method was ripped directly from taolib.CoreLib.Parser.WiggleIO
    
    """
    fhd = None
    if type(file) == str:
        fhd = open(file, "r")
    elif type(file) == file:
        fhd = file
    else:
        raise Exception ("f must be a filename or a file handler")
    
    data = {}
    chrom = "Unknown"
    for i in fhd:
        if i.startswith("track"):
            continue
        elif i.startswith("browse"):
            continue
        elif i.startswith("variableStep"):
            ci = i.rfind("chrom=")  # where the 'chrom=' is
            si = i.rfind("span=")   # where the 'span=' is
            if ci != -1:
                chrom = i[i.rfind("chrom=")+6:].strip().split()[0]
            else:
                chrom = "Unknown"
            if si != -1:
                span = int(i[i.rfind("span=")+5:].strip().split()[0])
            else:
                span = 0
            data[chrom]=BinKeeperI(bin=8,chromosomesize=250000000)
            add = data[chrom].add
        else:
            (pos,value) = i.split()
            add(int(pos),float(value))
    fhd.seek(0)
    return data
    
def findSpan(file):
        """ Takes in a specific wig file and returns the span of the file in integer form. """
        f=open(file,'r')
        spanTmp=''
        text = f.read()
        spanPos = text.find('span')
        i=5
        while text[spanPos+i]!='\n':
            spanTmp=spanTmp+text[spanPos+i]
            i+=1
        span=int(spanTmp)
        if span < 10:
            span = 10                   # minimum span is 10
        return span


def percentile_score(data, percentile):
    """ takes in data (a list containing floats or integers) and returns the value of a certain percentile in float form.  For example if the percentile is .5 this function will return the median of the list. """
    data.sort()
    length = len(data)
    if (length+1)*percentile%100 == 0:
        value = data[(length+1)*percentile/100-1]
    else: value = (data[int(length+1)*percentile/100-1]+data[int(length+1)*percentile/100])/2.0
    return value
            

class Heatmap:
    def __init__(self, data=None, xlim=None, groups=None, scores=None):
        """Initialize a heatmap.
        
        If given, initialize data with a list of lists or a 2-D NumPy array.
        If given, initialize xmin and xmax from xlim = (xmin, xmax).
        If given, initialize groups as a set of row indices.
        If given, initialize scores as a numeric array.
        
        xlim, groups, and scores will be initialized only if data is initialized."""
        if data is None:
            self.data = None
            self.xmin = None
            self.xmax = None
            self.groups = None
            self.scores = None
        else:
            try:
                self.data = numpy.array(data, float)
            except ValueError:
                raise ValueError, "ERROR Heatmap can only be initialized with a list of lists or a 2-D NumPy array"
            if self.data.ndim != 2:
                raise ValueError, "ERROR Heatmap can only be initialized with a list of lists or a 2-D NumPy array"
            if xlim is None:
                self.xmin = 0.5
                self.xmax = 0.5 + float(self.data.shape[1])
            else:
                try:
                    self.xmin = float(xlim[0])
                    self.xmax = float(xlim[1])
                except IndexError:
                    raise IndexError, "ERROR xlim can only be initialized with a pair of numbers: (xmin, xmax)"
                except ValueError:
                    raise ValueError, "ERROR xlim can only be initialized with a pair of numbers: (xmin, xmax)"
            if groups is None:
                self.groups = set()
            else:
                try:
                    self.groups = set(groups)
                except ValueError:
                    raise ValueError, "ERROR groups must be a set of row indices from the heatmap data"
                for elem in self.groups:
                    if elem is not int or elem < 0 or self.data.shape[1] <= elem:
                        raise ValueError, "ERROR groups must be a set of row indices from the heatmap data"
                self.groups.discard(0)  # remove 0 from group indices; 0 is implicitly understood to be the start of the first group
            if scores is None:
                self.scores = numpy.zeros(self.data.shape[0], float)
            else:
                try:
                    self.scores = numpy.array(scores, float)
                except ValueError:
                    raise ValueError, "ERROR scores must be a numeric list or a 1-D NumPy array"
                if self.scores.ndim != 1:
                    raise ValueError, "ERROR scores must be a numeric list or a 1-D NumPy array"
                if len(self.scores) != self.data.shape[0]:
                    raise ValueError, "ERROR scores must have length equal to the number of rows in the heatmap data"
    
    def align_wig(self, wigfilename, bedfilename, upstream, downstream, step, offset, span=None):
        # store parameters
        self.xmin = int(round(float(-upstream - offset) / step)) * step + offset
        self.xmax = int(round(float(downstream - offset) / step)) * step + offset


        # for output bed file
        self.locations = []
        # read alignment locations and scores from file
        locations = []  # list to store alignment locations
        bedfile = open(bedfilename, 'r')
        linenum = 0
        for line in bedfile:
            linenum += 1
            # split line by tabs
            row = line.strip().split('\t')
            # retrieve chromosome name
            try:
                chromosome = row[0]
            except IndexError:
                print >> sys.stderr, "WARNING Cannot retrieve chromosome from line %d of '%s'" % (linenum, bedfilename)
                print >> sys.stderr, "        This line will be skipped:"
                print >> sys.stderr, "        '%s'" % line[:-1]
                continue
            # retrieve chromosome strand
            try:
                strand = row[5]
            except IndexError:
                strand = '+'
            if not (strand == '+' or strand == '-'):
                strand = '+'
            # retrieve chromosome position
            try:
                if strand == '+':
                    position = int(row[1])
                else: # strand == '-'
                    position = int(row[2])
            except IndexError:
                print >> sys.stderr, "WARNING Cannot retrieve position from line %d of '%s'" % (linenum, bedfilename)
                print >> sys.stderr, "        This line will be skipped:"
                print >> sys.stderr, "        '%s'" % line[:-1]
                continue
            except ValueError:
                print >> sys.stderr, "WARNING Cannot convert position to an integer in line %d of '%s'" % (linenum, bedfilename)
                print >> sys.stderr, "        This line will be skipped:"
                print >> sys.stderr, "        '%s'" % line[:-1]
                continue
            # retrieve score (if provided)
            try:
                score = float(row[3])
            except IndexError:
                score = None
            except ValueError:
                score = None
            # store alignment location
            self.locations.append( (chromosome,row[1],row[2],strand) )
            locations.append((chromosome, position, strand, score))
        bedfile.close()
        locations.sort()
        self.locations.sort()
        self.locations = numpy.array(self.locations)
        
        # calculate parameters for generation of WIG alignment
        num_locations = len(locations)                # number of genome locations to align
        num_backsteps = -self.xmin / step             # number of steps upstream from each genome location to include in each alignment
qq

num_samples = (self.xmax - self.xmin) / step  # total number of steps to include in each alignment
        
        # generate WIG alignment data
        self.data = numpy.empty((num_locations, num_samples), float)
        self.groups = set()  # row indices of the first rows of each new group of rows
        wig_chromosome = None

        bk = build_binKeeper(wigfilename) #building bins with point data from each chromosome
        if not span:
            span = findSpan(wigfilename)
        
        for i, (chromosome, position, strand, score) in enumerate(locations):
            # retrieve WIG data if alignment location is in a new chromosome
            if chromosome != wig_chromosome:
###               wig_chromosome = chromosome
                self.groups.add(i)
            # round chromosome position to the nearest sampled location
            position = int(round(float(position - offset) / step)) * step + offset
            # set start position, end position, and increment for sampling
            if strand == '+':
                startpos = position - (step * num_backsteps)
                endpos = startpos + (step * num_samples)
                increment = step
            else:  # strand == '-'
                startpos = position + (step * num_backsteps)
                endpos = startpos - (step * num_samples)
                increment = -step
            # retrieve and store aligned data points from WIG data
            for j, position in enumerate(xrange(startpos, endpos, increment)):
                if chromosome in bk:
                    
                    values = bk[chromosome].pp2v(position-span,position)
                    if len(values)>0:
                        if increment > 0: # + strand
                            self.data[i,j] = values[-1]
                        else:           # - strand
                            self.data[i,j] = values[0]
                    else: self.data[i,j] = 0.0
                
        self.groups.discard(0)  # remove 0 from group indices; 0 is implicitly understood to be the start of the first group
        
        # store scores (if provided for every location)
        scores = [score for (chromosome, position, strand, score) in locations]
        num_missing_scores = sum(score is None for score in scores)
        if num_missing_scores == 0:
            self.scores = numpy.array(scores, float)
        else:
            self.scores = numpy.zeros(num_locations, float)
        # store 5' locations
        #self.locations = [chromosome+":"+str(position) for (chromosome,position,strand,score) in locations]
    
    def summarize_wigs(self, wigfilenames, bedfilename, step, offset, summary='mean', span=None):
        # set function for summarizing WIG data
        if summary == 'mean':
            _summarize = numpy.mean
        elif summary == 'median':
            _summarize = numpy.median
        elif summary == 'max':
            _summarize = numpy.max
        elif summary == 'min':
            _summarize = numpy.min
        else:
            raise ValueError, "'summary' parameter must be 'mean', median', 'max', or 'min'"
        
        # read genome regions from BED file
        regions = []  # list to store genome regions
        bedfile = open(bedfilename, 'r')
        linenum = 0
        for line in bedfile:
            linenum += 1
            # split line by tabs
            row = line.strip().split('\t')
            # retrieve chromosome name
            try:
                chromosome = row[0]
            except IndexError:
                print >> sys.stderr, "WARNING Cannot retrieve chromosome from line %d of '%s'" % (linenum, bedfilename)
                print >> sys.stderr, "        This line will be skipped:"
                print >> sys.stderr, "        '%s'" % line[:-1]
                continue
            # retrieve start position
            try:
                start = int(row[1])
            except IndexError:
                print >> sys.stderr, "WARNING Cannot retrieve start position from line %d of '%s'" % (linenum, bedfilename)
                print >> sys.stderr, "        This line will be skipped:"
                print >> sys.stderr, "        '%s'" % line[:-1]
                continue
            except ValueError:
                print >> sys.stderr, "WARNING Cannot convert start position to an integer in line %d of '%s'" % (linenum, bedfilename)
                print >> sys.stderr, "        This line will be skipped:"
                print >> sys.stderr, "        '%s'" % line[:-1]
                continue
            # retrieve end position
            try:
                end = int(row[2])
            except IndexError:
                print >> sys.stderr, "WARNING Cannot retrieve end position from line %d of '%s'" % (linenum, bedfilename)
                print >> sys.stderr, "        This line will be skipped:"
                print >> sys.stderr, "        '%s'" % line[:-1]
                continue
            except ValueError:
                print >> sys.stderr, "WARNING Cannot convert end position to an integer in line %d of '%s'" % (linenum, bedfilename)
                print >> sys.stderr, "        This line will be skipped:"
                print >> sys.stderr, "        '%s'" % line[:-1]
                continue
            # store alignment location
            regions.append((chromosome, start, end))
        bedfile.close()
        regions.sort()
        
        # calculate parameters for generation of WIG summary
        num_regions = len(regions)       # number of genome regions to summarize
        num_samples = len(wigfilenames)  # total number of steps to include in each alignment
        
        # generate WIG summary data
        self.data = numpy.empty((num_regions, num_samples), float)
        self.groups = set()  # row indices of the first rows of each new group of rows
        wig_chromosome = None

        bks = []
        spans = []
        for wigfilename in wigfilenames:
            bk = build_binKeeper(wigfilename)
            if not span:
                span = findSpan(wigfilename)
            bks.append(bk)
            spans.append(span)
        
        for i, (chromosome, start, end) in enumerate(regions):
            # retrieve WIG data if genome region is in a new chromosome
            if chromosome != wig_chromosome:
                wig_chromosome = chromosome
                self.groups.add(i)
            # round start and end positions to the nearest sampled location
            start = int(round(float(start - offset) / step)) * step + offset
            end = int(round(float(end - offset) / step)) * step + offset + step  # add extra step to make include end position
            # retrieve data points from WIG data and store summarized data
            num_wigvalues = (end - start) / step
            for j, bk in enumerate(bks):
                wigvalues = numpy.empty(num_wigvalues, float)
                for k, position in enumerate(xrange(start, end, step)):
                    if chromosome in bk:
                        values = bk[chromosome].pp2v(position-span,position)
                        if len(values)>0:
                            wigvalues[k] = values[-1]
                        else:
                            wigvalues[k] = 0.0
                self.data[i,j] = _summarize(wigvalues)
        self.groups.discard(0)  # remove 0 from group indices; 0 is implicitly understood to be the start of the first group
        
        # set default for other values
        self.xmin = 0.5
        self.xmax = 0.5 + float(self.data.shape[1])
        self.scores = numpy.zeros(num_regions, float)
    
    def sort(self, scores=None, group_threshold=None):
        """Sort rows of data based on scores.
        
        Heatmap scores are used as the default scores if scores is not specified as a list.
        
        Group threshold is the minimum difference between scores that specify a new group.
        
        For example, suppose that scores, after sorting, is [1, 3, 4, 5, 8, 9].
        If the group threshold is given as 2, there will be three groups:
        1. [1]
        2. [3, 4, 5]
        3. [8, 9]
        
        By default, group_threshold is None, meaning no groups will be specified."""
        if self.data is None:
            print >> sys.stderr, "WARNING Data has not been initilized; nothing to sort"
            return
        if scores is None:
            scores = self.scores
        else:
            scores = numpy.array(scores)
        I = range(self.data.shape[0])      # I = new sorted order of old unsorted indices
        I.sort(key=lambda i: scores[i])
        scores = scores[I]
        self.data = self.data[I]
        self.scores = self.scores[I]
        self.locations = self.locations[I]
        self.groups = set()
        if type(group_threshold) in (int, long, float):
            self.groups.update((numpy.diff(scores) >= group_threshold).nonzero()[0] + 1)
    
    def kmeans(self, k):
        """Cluster rows of data based on k-means."""
        if self.data is None:
            print >> sys.stderr, "WARNING Data has not been initilized; nothing to cluster"
            return
        # sort rows of data based on the Euclidean norm of each row of data
        data_norms = [numpy.linalg.norm(row) for row in self.data]
        self.sort(data_norms)
        # normalize standard deviation of each column of data to unit variance
        normalized = scipy.cluster.vq.whiten(self.data)
        # cluster using k-means algorithm
        cluster_centers, cluster_indices = scipy.cluster.vq.kmeans2(normalized, k, 72, minit='points')
        # sort clusters based on the Euclidean norm of each cluster center
        cluster_norms = [numpy.linalg.norm(center) for center in cluster_centers]
        I = range(len(cluster_norms))             # I = new sorted order of old unsorted cluster indices
        I.sort(key=lambda i: cluster_norms[i])
        J = [I.index(i) for i in xrange(len(I))]  # J = old unsorted order of new sorted cluster indices
        J = numpy.array(J)
        cluster_centers = cluster_centers[I]
        cluster_indices = J[cluster_indices]
        # sort rows of data based on cluster indices
        #print len(I)
        #print len(J)
        #print len(cluster_centers)
        #print len(cluster_indices)
        self.cluster_indices = sorted(cluster_indices)
        self.sort(scores=cluster_indices, group_threshold=1)

        # return data
        return cluster_centers, cluster_indices
    
    def draw(self, heatmapfilename=None, size=(16,12), zmin=None, zmax=None, xlabel=None, ylabel=None, title=None, axhline=False, axvline=False, grid=False, colormap='default', colorbar=False):
        """Draw heatmap of data."""
        if self.data is None:
            print >> sys.stderr, "WARNING Data has not been initilized; nothing to draw"
            return
        print colormap
        if colormap in pylab.cm.datad:
            colormap = pylab.get_cmap(colormap)
        else:
            colormap = None
        if zmin is None or zmax is None:
            flattened_data = self.data.flatten()
            (Q1, Q2, Q3) = [percentile_score(flattened_data, percentile) for percentile in (25, 50, 75)]
            IQR = Q3 - Q1
        if zmin is None:
            zmin = max(self.data.min(), Q2 - 1.5 * IQR)
        if zmax is None:
            zmax = min(self.data.max(), Q2 + 1.5 * IQR)
        ymin =  0.5 + self.data.shape[0]
        ymax = 0.5
        h = pylab.figure(figsize=size)
        pylab.rc('font', size=(2 * size[1]))
        left = 0.1
        bottom = 0.1
        width = 0.8
        height = 0.8
        pylab.axes([left, bottom, width, height])
        #print colormap
        #print self.data
        pylab.imshow(self.data, cmap=colormap, aspect='auto', interpolation='nearest', vmin=zmin, vmax=zmax, extent=(self.xmin, self.xmax, ymin, ymax))
        #print self.data,self.xmin,self.xmax,ymin,ymax,zmin,zmax
        if type(xlabel) is str:
            pylab.xlabel(xlabel, fontsize='x-small')
#         if type(ylabel) is str:
#             pylab.ylabel(ylabel, fontsize='x-small')
        if type(title) is str:
            pylab.title(title, fontsize='small')
        if axhline is True:
            for group_index in self.groups:
                pylab.axhline(0.5 + group_index, color='black', linewidth=2)
            pylab.ylim(ymin, ymax)
        if axvline is True:
            pylab.axvline(color='black')
            pylab.xlim(self.xmin, self.xmax)
        if grid is True:
            pylab.grid()
        pylab.xticks(fontsize='xx-small')
        pylab.yticks(fontsize='xx-small')
        if colorbar is True:
            cb = pylab.colorbar()
            pylab.axes(cb.ax)
            pylab.yticks(fontsize='xx-small')
        if heatmapfilename is None:
            pylab.show()
        else:
            pylab.savefig(str(heatmapfilename))
            pylab.close(h)

class Heatmaps(list):
    def __init__(self, iterable=[]):
        list.__init__(self, iterable)
        for elem in self:
            if not isinstance(elem, Heatmap):
                raise ValueError, "not an Heatmap object: %s" % repr(elem)
    
    def __setitem__(self, index, value):
        list.__setitem__(self, index, value)
        if not isinstance(value, Heatmap):
            raise ValueError, "not an Heatmap object: %s" % repr(elem)
    
    def __setslice__(self, index_start, index_end, iterable):
        list.__setslice__(self, index_start, index_end, iterable)
        for elem in self:
            if not isinstance(elem, Heatmap):
                raise ValueError, "not an Heatmap object: %s" % repr(elem)
    
    def append(self, value):
        list.append(self, value)
        if not isinstance(value, Heatmap):
            raise ValueError, "not an Heatmap object: %s" % repr(elem)
    
    def extend(self, iterable):
        list.extend(self, iterable)
        for elem in self:
            if not isinstance(elem, Heatmap):
                raise ValueError, "not an Heatmap object: %s" % repr(elem)
    
    def insert(self, index, value):
        list.append(self, index, value)
        if not isinstance(value, Heatmap):
            raise ValueError, "not an Heatmap object: %s" % repr(elem)
    
    def kmeans(self, k):
        """Cluster rows of data across all heatmaps based on k-means.
        
        All heatmaps must have the same number of rows of data."""
        # check that all Heatmap objects have been initialized with data
        for heatmap in self:
            if heatmap is None:
                print >> sys.stderr, "WARNING not all Heatmap objects have been initialized with data; cannot cluster"
                return
        # check that all Heatmap objects have the same number of rows
        m = None
        for heatmap in self:
            if m is None:
                m = heatmap.data.shape[0]
            elif heatmap.data.shape[0] != m:
                raise ValueError, "Heatmap objects have unequal numbers of rows"
        # horizontally concatenate all Heatmap objects
        data = numpy.hstack([heatmap.data for heatmap in self])
        # sort rows of data based on the Euclidean norm of each row of data
        data_norms = [numpy.linalg.norm(row) for row in data]
        for heatmap in self:
            heatmap.sort(data_norms)
        # horizontally concatenate all Heatmap objects again
        data = numpy.hstack([heatmap.data for heatmap in self])
        # normalize standard deviation of each column of data to unit variance
        normalized = scipy.cluster.vq.whiten(data)
        # cluster using k-means algorithm
        cluster_centers, cluster_indices = scipy.cluster.vq.kmeans2(normalized, k, 72, minit='points')
        # sort clusters based on the Euclidean norm of each cluster center
        cluster_norms = [numpy.linalg.norm(center) for center in cluster_centers]
        I = range(len(cluster_norms))             # I = new sorted order of old unsorted cluster indices
        I.sort(key=lambda i: cluster_norms[i])
        J = [I.index(i) for i in xrange(len(I))]  # J = old unsorted order of new sorted cluster indices
        J = numpy.array(J)
        cluster_centers = cluster_centers[I]
        cluster_indices = J[cluster_indices]
        self.cluster_indices = sorted(cluster_indices)
        # sort rows of data based on cluster indices
        for heatmap in self:
            heatmap.sort(scores=cluster_indices, group_threshold=1)
        return cluster_centers, cluster_indices            
    
    def draw(self, heatmapfilename=None, size=(16,12), zmin=None, zmax=None, xlabel=None, ylabel=None, title=None, subtitle=None, axhline=False, axvline=False, grid=False, colormap='default', colorbar=False):
        """Draw heatmap of data."""
        # check that all Heatmap objects have been initialized with data
        for heatmap in self:
            if heatmap is None:
                print >> sys.stderr, "WARNING not all Heatmap objects have been initialized with data; cannot draw"
                return
        # set color map for mapping values to colors
        if colormap in pylab.cm.datad:
            colormap = pylab.get_cmap(colormap)
        else:
            colormap = None
        # set minimum and maximum values for color map
        if zmin is None or zmax is None:
            flattened_data = numpy.hstack([heatmap.data.flatten() for heatmap in self])

            (Q1, Q2, Q3) = [percentile_score(flattened_data, percentile) for percentile in (25, 50, 75)]
            IQR = Q3 - Q1
        if zmin is None:
            zmin = max(flattened_data.min(), Q2 - 1.5 * IQR)
        if zmax is None:
            zmax = min(flattened_data.max(), Q2 + 1.5 * IQR)
        # set range of values on y-axis
        ymin = 0.5 + self[0].data.shape[0]
        ymax = 0.5
        # set x-axis label for each heatmap
        if type(xlabel) not in (tuple, list):
            xlabel = len(self) * [xlabel]
        # set y-axes label for the first heatmap
        if type(ylabel) not in (tuple, list):
            ylabel = [ylabel]
        # set subtitles for each heatmap
        if type(subtitle) not in (tuple, list):
            subtitle = len(self) * [subtitle]
        # create figure
        h = pylab.figure(figsize=size)
        # scale font sizes to the size of figure
        pylab.rc('font', size=(2 * size[1]))
        # set parameters for determining the position and size of each heatmap
        n = len(self)  # number of heatmaps
        p = 0.10       # padding along edge of figure as a percentage of figure height/width (same for all sides)
        s = 0.10       # spacing between heatmaps as a percentage of heatmap width
        if colorbar is True:
            c = 0.25   # spacing for color bar as a percentage of heatmap width
        else:
            c = 0.00   # spacing for color bar as a percentage of heatmap width
        width = (1 - 2*p) / (n + (n-1)*s + c)  # heatmap width
        height = 1 - 2*p                       # heatmap height
        bottom = p                             # bottom boundary of heatmap
        # create and draw heatmaps
        for j, heatmap in enumerate(self):
            # set range of values on x-axis
            # set another parameter for determining the position the heatmap
            left = p + j * (1 + s) * width  # left boundary of heatmap
            # create axes to hold heatmap (use correct position and size)
            if j == 0:
                ax = pylab.axes([left, bottom, width, height])
                ax1 = ax
            elif j == (n - 1):
                ax = pylab.axes([left, bottom, (1 + c) * width, height], sharey=ax1)
            else:
                ax = pylab.axes([left, bottom, width, height], sharey=ax1)
            # draw heatmap
            pylab.imshow(heatmap.data, cmap=colormap, aspect='auto', interpolation='nearest', vmin=zmin, vmax=zmax, extent=(heatmap.xmin, heatmap.xmax, ymin, ymax))
            # label x-axis
            if j < len(xlabel) and type(xlabel[j]) is str:
                pylab.xlabel(xlabel[j], fontsize='x-small')
            # label y-axis
            # only for the first y-axis
            if j == 0 and type(ylabel[j]) is str:
                pylab.ylabel(ylabel[j], fontsize='x-small')
            # write subtitle
            if j < len(subtitle) and type(subtitle[j]) is str:
                pylab.title(subtitle[j], fontsize='small')
            # draw horizontal lines to divide rows into groups
            if axhline is True:
                for group_index in heatmap.groups:
                    pylab.axhline(0.5 + group_index, color='black', linewidth=2)
                pylab.ylim(ymin, ymax)
            # draw vertical line at x = 0
            if axvline is True:
                pylab.axvline(color='black')
                pylab.xlim(heatmap.xmin, heatmap.xmax)
            # draw grid
            if grid is True:
                pylab.grid()
            # specify font size for numerical labels along the x and y axdes
            pylab.xticks(fontsize='xx-small')
            if j == 0:
                pylab.yticks(fontsize='xx-small')
            else:
                pylab.yticks([])
    
            # draw color bar
            if j == (n - 1) and colorbar is True:
                cb = pylab.colorbar()
                pylab.axes(cb.ax)
                pylab.yticks(fontsize='xx-small')
        # write title
        #if type(title) is str:
        #    pylab.title(title)
        # display or save heatmap
        if heatmapfilename is None:
            pylab.show()
        else:
            pylab.savefig(str(heatmapfilename))
            pylab.close(h)

if __name__ == '__main__':
    # parse arguments
    sep = ','  # separator character for delimiting a list of WIG files, BED files, or other options
    usage = """Usage: wigheatmap.py [options] WIGFILENAME(S) BEDFILENAME(S)

For multiple WIG/BED files, provide a list of file names separated by '%s'.
User must specify the only one BED file.""" % sep
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-o', '--output-file', dest='heatmapfilename', default='output.png')
    parser.add_option('-u', '--upstream', dest='upstream', default='500')
    parser.add_option('-d', '--downstream', dest='downstream', default='500')
    parser.add_option('-p', '--step', dest='step', default='10')
    parser.add_option('-f', '--offset', dest='offset', default='0')
    parser.add_option('-s', '--sort', action='store_true', dest='sort', default=False)
    parser.add_option('-k', '--k-means', type='int', dest='kmeans', default=0)
    parser.add_option('-W', '--width', type='int', dest='width', default=1600)
    parser.add_option('-H', '--height', type='int', dest='height', default=1200)
    parser.add_option('-m', '--minimum-scale', type='float', dest='zmin', default=None)
    parser.add_option('-M', '--maximum-scale', type='float', dest='zmax', default=None)
    parser.add_option('-x', '--x-label', dest='xlabel', default=None)
    parser.add_option('-y', '--y-label', dest='ylabel', default=None)
    parser.add_option('-t', '--title', dest='title', default=None)
    parser.add_option('-l', '--subtitle', dest='subtitle', default=None)
    parser.add_option('-z', '--horizontal-line', action='store_true', dest='axhline', default=False)
    parser.add_option('-v', '--vertical-line', action='store_true', dest='axvline', default=False)
    parser.add_option('-g', '--grid', action='store_true', dest='grid', default=False)
    parser.add_option('-c', '--color-map', dest='colormap', default='default')
    parser.add_option('-b', '--color-bar', action='store_true', dest='colorbar', default=False)
    parser.add_option('-S', '--span', dest='span',type='int',default=None)
    (opts, args) = parser.parse_args(sys.argv)
    # check for 2 required arguments
    if len(args) != 3:
        parser.print_help()
        sys.exit(1)
        parser.error("must provide 2 required arguments specifying WIGFILENAME(S) and BEDFILENAME(S)")
    # decide whether user has specified a single WIG file or multiple WIG files
    if sep not in args[1]:
        wigfilename = args[1]
        bedfilename = args[2]
        # parse numeric options
        def _parse_int_option(opt_name, opt_value):
            try:
                parsed_opt = int(opt_value)
            except ValueError:
                parser.error("option %s: invalid integer value: %s" % (opt_name, opt_value))
            return parsed_opt
        upstream = _parse_int_option('--upstream', opts.upstream)
        downstream = _parse_int_option('--downstream', opts.downstream)
        step = _parse_int_option('--step', opts.step)
        offset = _parse_int_option('--offset', opts.offset)
        figsize = (opts.width / 100.0, opts.height / 100.0)
        # align WIG file and draw heatmap
        A = Heatmap()
        A.align_wig(wigfilename, bedfilename, upstream, downstream, step, offset, span=opts.span)
        if opts.sort is True:
            A.sort()
        if opts.kmeans > 1:
            (tmp1,tmp2) = A.kmeans(opts.kmeans)
            #bfhd = open (bedfilename)
            outfhd = open (opts.heatmapfilename+".bed","w")
            #i = 0
            for i in xrange(len(A.locations)):
                l = A.locations[i]
                c = A.cluster_indices[i]
                outfhd.write("%s\t%s\t%s\tNA\t%d\t%s\n" % (l[0],l[1],l[2],c,l[3]))
            #for l in bfhd:
            #    if l.startswith('#') or l.startswith('browse') or l.startswith('track'):
            #        continue
            #    else:
            #        ls = l.split()
            #        if len(ls) == 3:
            #            outfhd.write( "\t".join(map(str,(ls[0],ls[1],ls[2],'.',tmp2[i])))+"\n" )
            #        elif len(ls) == 4:
            #            ls.append(tmp2[i])
            #            outfhd.write( "\t".join(map(str,ls))+"\n" )
            #        else:
            #            ls[4] = tmp2[i]
            #            outfhd.write( "\t".join(map(str,ls))+"\n" )
            #        i+=1
            #bfhd.close()
            outfhd.close()
        A.draw(opts.heatmapfilename, figsize, opts.zmin, opts.zmax, opts.xlabel, opts.ylabel, opts.title, opts.axhline, opts.axvline, opts.grid, opts.colormap, opts.colorbar)
    else:
        # multiple wiggle files
        wigfilenames = args[1].rstrip(sep).split(sep)
        bedfilenames = len(wigfilenames) * [args[2]]

        # parse numeric options that could have multiple inputs separated by the 'sep' delimiter
        # if multiple inputs are provided, make sure there is one for each WIG file

        def _parse_multiple_int_option(opt_name, opt_value):
            if sep in opt_value:
                try:
                    parsed_opt = [int(elem) for elem in opt_value.split(sep)]
                except ValueError:
                    parser.error("option %s: invalid integer value(s): %s" % (opt_name, opt_value))
                if len(parsed_opt) != len(wigfilenames):
                    parser.error("option %s: if a list, must provide the same number of values as the number of WIG files" % opt_name)
            else:
                try:
                    parsed_opt = len(wigfilenames) * [int(opt_value)]
                except ValueError:
                    parser.error("option %s: invalid integer value(s): %s" % (opt_name, opt_value))
            return parsed_opt
        upstreams = _parse_multiple_int_option('--upstream', opts.upstream)
        downstreams = _parse_multiple_int_option('--downstream', opts.downstream)
        steps = _parse_multiple_int_option('--step', opts.step)
        offsets = _parse_multiple_int_option('--offset', opts.offset)
        figsize = (opts.width / 100.0, opts.height / 100.0)
        # parse string options that could have multiple inputs separated by the 'sep' delimiter
        def _parse_multiple_str_option(opt_name, opt_value):
            if opt_value is not None and sep in opt_value:
                parsed_opt = opt_value.split(sep)
            else:
                parsed_opt = opt_value
            return parsed_opt
        xlabel = _parse_multiple_str_option('--x-label', opts.xlabel)
        ylabel = _parse_multiple_str_option('--y-label', opts.ylabel)
        subtitle = _parse_multiple_str_option('--subtitle', opts.subtitle)
        # align WIG files and draw heatmaps
        A = []
        for i in xrange(len(wigfilenames)):
            a = Heatmap()
            a.align_wig(wigfilenames[i], bedfilenames[i], upstreams[i], downstreams[i], steps[i], offsets[i],span=opts.span)
            A.append(a)
        if opts.sort is True:
            for a in A:
                a.sort()
        B = Heatmaps(A)
        if opts.kmeans > 1:
            (tmp1,tmp2) = B.kmeans(opts.kmeans)
            #bfhd = open (bedfilenames[i])
            outfhd = open (opts.heatmapfilename+".bed","w")
            for i in xrange(len(B[0].locations)):
                l = B[0].locations[i]
                c = B.cluster_indices[i]
                outfhd.write("%s\t%s\t%s\tNA\t%d\t%s\n" % (l[0],l[1],l[2],c,l[3]))
            #i = 0
            #for l in bfhd:
            #    if l.startswith('#') or l.startswith('browse') or l.startswith('track'):
            #        continue
            #    else:
            #        ls = l.split()
            #        if len(ls) == 3:
            #            outfhd.write( "\t".join(map(str,(ls[0],ls[1],ls[2],'.',tmp2[i])))+"\n")
            #        elif len(ls) == 4:
            #            ls.append(tmp2[i])
            #            outfhd.write( "\t".join(map(str,ls))+"\n")
            #        else:
            #            ls[4] = tmp2[i]
            #            outfhd.write( "\t".join(map(str,ls))+"\n")
            #        i+=1
            #bfhd.close()
            outfhd.close()

                    
        B.draw(opts.heatmapfilename, figsize, opts.zmin, opts.zmax, xlabel, ylabel, opts.title, subtitle, opts.axhline, opts.axvline, opts.grid, opts.colormap, opts.colorbar)
