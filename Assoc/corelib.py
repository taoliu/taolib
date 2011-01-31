
"""Module Description

Copyright (c) 2008 H. Gene Shin <shin@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  H. Gene Shin
@contact: shin@jimmy.harvard.edu

This module contains functions that produce R scripts corresponding to certain R functions. 
For example, R.plot(x,y) produces a R script of plot(x,y).

For the details on the parameters of the below functions, refer to R documents.
"""

# ------------------------------------
# Python modules
# ------------------------------------
import sys,operator,copy
import itertools
from array import *
import random
import re
import bisect

# ------------------------------------
# functions
# ------------------------------------
def cumsum(arr):
    """Do cummulative summation"""
    
    cs=[arr[0]]
    for a in arr[1:]: cs.append(cs[-1]+a)
    
    return cs

def sum2(arr, NaN=False):
    """Return the sum of an array.
    If NaN==True, NaNs are ingnored and return the sum of the other elements"""
    
    if NaN:
        s=0
        for a in arr:
            if a==a:
                s+=a
    else:
        s=sum(arr)
    
    return s

def mean(arr, NaN=False):
    """Take a mean of array. 
    If NaN==True, NaNs are excluded"""
    
    if NaN:
        c = sum(map(lambda x: x==x, arr))
        s = sum(map(lambda x: (x==x) and x or 0, arr))
        
        try:
            return 1.0*s/c
        except ZeroDivisionError:
            return float('nan')
    
    else:
        try:
            return 1.0*sum(arr)/len(arr)
        except ZeroDivisionError:
            return 0.0
    
def weight_mean(arr, wei, NaN=False):
    """Perform the weight mean.
    
    Parameters:
    1. arr: the array of numbers
    2. wei: the array of weights.
    
    Each number in arr is weighted by its corresponding weight in wei and added. Then divided by the total weight (sum of weight)
    
    If NaN=True, NaNs in arr are rejected in calculating the weighted mean. 
    BE CAUTIONS! if NaN=False, but there are NaNs in arr, an wanted result may come.
    """
    
    if NaN:
        sw = sum(map(lambda x, y: (x==x) and y, arr, wei))
        swa = sum(map(lambda x, y: (x==x) and x*y or 0, arr, wei))
        
        try:
            return 1.0*swa/sw
        except ZeroDivisionError:
            return float('nan')
    
    else:   # when NaNs are not considered. This is faster, but make sure that there are no NaNs in arr and wei.
        sw = sum(wei)
        swa = sum(map(lambda x, y: x * y, arr, wei))
        
        try:
            return 1.0*swa/sw
        except ZeroDivisionError:
            return 0.0
    
def array_mutex(x,y):
    """Mutual exclusion for x and y"""
    
    return [(not x and y) or (x and not y) for x,y in itertools.izip(x,y)]
    
def array_union(x,y):
    """Union x and y"""
    
    return [x or y for x,y in zip(x,y)]

def array_intersection(x,y):
    """Intersection x and y"""

    return [x and y for x,y in zip(x,y)]

def array_extractor(x,y):
    """Extract array elements indicated by y. 
    
    y could be a numeric array or boolean array"""
    
    if type(y[0]).__name__=='bool':
        return [x[i] for i in xrange(0,len(y)) if y[i]==True]  
    else:
        return [x[i] for i in y]
    
def array_rearranger(x,y):
    """Rearrange x according to y"""
    
    return [x[i] for i in y]

def array_negator(x):
    """Negate the given array
    
    Parameters:
    1. x: an array of boolean values or any values. If x is an array of non-boolean values, just check if each element is NULL or not
    """
    
    return map(lambda b: not b, b)

def is_nan(n):
    """Return True if n is NaN"""
    
    try:
        return map(lambda x: not (x==x), n)
    
    except TypeError:   # probably just integer
        
        return not (n == n)
        
    
def is_inf(n):

    ntype=type(n).__name__
    if ntype=='list' or ntype=='tuple' or ntype=='array':
        return [not i==float('inf') for i in n]
    else:
        return n==n
    
def array_adder(x,y):
    """Add two arrays together
    
    if either of elements is nan, have non-nan value. If both of then are nan, put nan"""
    
    add=[]
    xnan=is_nan(x)
    ynan=is_nan(y)
    and_nan=array_intersection(xnan,ynan)
    
    x2=[]
    for a,o in zip(x,xnan):
        if o: x2.append(0)
        else: x2.append(a)
    y2=[]
    for b,o in zip(y,ynan):
        if o: y2.append(0)
        else: y2.append(b)
         
    add=[a+b for a,b in zip(x2,y2)]
    for i in xrange(0,len(add)):
        if and_nan[i]: add[i]=float("nan")

    count=[1*(not xn)+1*(not yn) for xn,yn in zip(xnan,ynan)]
    
    return add,count

    
def sum_col_by_col(arrays):
    """Perform col by col summation"""
    
    ncol=len(arrays[0])
    nrow=len(arrays)
    colsum=[]
    for j in xrange(0,ncol):
        cs=0
        for i in xrange(0,nrow):
            if arrays[i][j]==arrays[i][j]:
                cs+=arrays[i][j]
        colsum.append(cs)
    
    return colsum
    
        
def mean_col_by_col(arrays, counts=True, NaN=True):
    """The second version of mean_col_by_col"""
    
    ncols= map(lambda a: len(a), arrays)
    ncols.sort()
    try:
        ncol = ncols[-1]
    except IndexError:
        if counts:
            return [],[]
        else:
            return []
    
    colmean = [] 
    colcount = []
    try:
        for i in xrange(ncol):
            col = map(operator.itemgetter(i), arrays)
            colmean.append(mean(col, NaN=NaN))
            colcount.append(sum(map(lambda x: x==x, col)))
    except IndexError:  # when some [] is in the array
        temp = [a for a in arrays if a]
        for i in xrange(ncol):
            col = map(operator.itemgetter(i), temp)
            colmean.append(mean(col, NaN=NaN))
            colcount.append(sum(map(lambda x: x==x, col)))
            
    if counts:
        return colmean, colcount
    else:
        return colmean
        
        
def weight_mean_col_by_col(arrays, weights, counts=True, NaN=False):
    """The second version of mean_col_by_col"""
    
    ncols= map(lambda a: len(a), arrays)
    ncols.sort()
    try:
        ncol = ncols[-1]
    except IndexError:
        if counts:
            return [],[]
        else:
            return []
    
    colmean = [] 
    colcount = []
    
    try:
        for i in xrange(ncol):
            col = map(operator.itemgetter(i), arrays)
            wei = map(operator.itemgetter(i), weights)
            colmean.append(weight_mean(col, wei, NaN=NaN))
            colcount.append(sum(filter(lambda x: x==x, wei)))
    except IndexError:  # when some [] is in the array
        temp = [a for a, w in itertools.izip(arrays, weights) if a and w]
        tempweight = [w for a, w in itertools.izip(arrays, weights) if a and w]
        for i in xrange(ncol):
            col = map(operator.itemgetter(i), temp)
            wei = map(operator.itemgetter(i), tempweight)
            colmean.append(weight_mean(col, wei, NaN=NaN))
            colcount.append(sum(filter(lambda x: x==x, wei)))
            
    if counts:
        return colmean, colcount
    else:
        return colmean
    

def mean_col_by_col_multi_bins(mularrays,counts=False):
    """Weighted colum by column mean for multiple bins"""
    
    nbins=len(mularrays)
    mulcolmean=[[] for i in xrange(nbins)]
    if counts:
        mulcolcount=[[] for i in xrange(nbins)]
    for i in xrange(nbins):
        if counts:
            mulcolmean[i],mulcolcount[i]=mean_col_by_col(mularrays[i],counts=counts)
        else:
            mulcolmean[i]=mean_col_by_col(mularrays[i],counts=counts)
    
    if counts:
        return mulcolmean,mulcolcount
    else:
        return mulcolmean
    
    
def weight_mean_col_by_col_multi_bins(mularrays,mulweights,counts=False):
    """Weighted colum by column mean for multiple bins"""
    
    nbins=len(mularrays)
    mulcolmean=[[] for i in xrange(nbins)]
    if counts:
        mulcolcount=[[] for i in xrange(nbins)]
    for i in xrange(nbins):
        if counts:
            mulcolmean[i],mulcolcount[i]=weight_mean_col_by_col(mularrays[i],mulweights[i],counts=counts)
        else:
            mulcolmean[i]=weight_mean_col_by_col(mularrays[i],mulweights[i],counts=counts)
    
    if counts:
        return mulcolmean,mulcolcount
    else:
        return mulcolmean
    
    
def mean_binwise(profiles):
    """Assuming 'profiles' has n bins, perform binwise averaging"""
    
    ms,cs=[],[]
    for p in profiles:
        s=extend_list_series(p)
        m,c=mean_col_by_col(s,counts=True)
        ms.append(m)
        cs.append(c)
        
    return ms,cs


def sum_row_by_row(arrays):
    """Perform row by row summation
    
    arrays is a multi-dimensional list. The rows must have the same length"""
    
    rowsum=[]
    for arr in arrays: rowsum.append(sum2(arr, NaN=True))
    
    return rowsum


def mean_row_by_row(arrays):
    """Row-by-row mean
    
    Probably, needs to be upgraded in the future"""
    
    rowmean=[]
    for arr in arrays:
        rowmean.append(mean(arr,NaN=False))
    
    return rowmean
    
    
def min_col_by_col(arrays):
    """Perform column by column min"""
    
    colmin=[]
    for arr in arrays:
        if not colmin:
            colmin=arr[:]
        else:
            if NaN:
                colmin=[min_w_nan([colmin[i],arr[i]]) for i in xrange(0,len(arr))]
            else:
                colmin=[min(colmin[i],arr[i]) for i in xrange(0,len(arr))]
    return colmin


def min_row_by_row(arrays,NaN=False):
    """Perform row by row min"""
    
    if NaN:
        rowmin=[min_w_nan(arr) for arr in arrays]
    else:
        rowmin=[min(arr) for arr in arrays]
    
    return rowmin

def max_col_by_col(arrays,NaN=False):
    """Perform column by column max"""
    
    colmax=[]
    for arr in arrays:
        if not colmax:
            colmax=arr[:]
        else:
            if NaN:
                colmax=[max_w_nan([colmax[i],arr[i]]) for i in xrange(0,len(arr))]
            else:
                colmax=[max(colmax[i],arr[i]) for i in xrange(0,len(arr))]
    
    return colmax

def max_row_by_row(arrays,NaN=False):
    """Perform row by row min"""
    
    if NaN:
        rowmax=[max_w_nan(arr) for arr in arrays]
    else:
        rowmax=[max(arr) for arr in arrays]
    
    return rowmax


def filter_out_empty_rows(matrix):
    """Filter out empty rows in a matrix"""
    
    try:
        return filter(lambda x: x, matrix)
    except:
        return []

"""argmax and argmin lambda expression
    Reference to http://www.daniel-lemire.com/blog/archives/2004/11/25/computing-argmax-fast-in-python/"""
argmax=lambda array: max(itertools.izip(array, xrange(len(array))))[1]
argmin=lambda array: min(itertools.izip(array, xrange(len(array))))[1]
arglim=lambda array,lim: [a for a in array if a >= lim[0] and a <= lim[1]]

def arglim(array, lim):
    """Return the index of element between a lower-limit and upper-limit
    """
    
    index = []
    for a,i in itertools.izip(array,xrange(len(array))):
        if a >= lim[0] and a <= lim[1]: 
            index.append(i)
    
    return index

def max_w_nan(array):
    """Return max ignoring nan"""
    
    return max([a for a in array if a==a])

def min_w_nan(array):
    """Return min ignoring nan"""
    
    return min([a for a in array if a==a])

def max_diff_w_nan(array):
    """Return the maximum different between element in an array"""
    
    arr_no_nan=[a for a in array if a==a]
    arr_no_nan.sort()
    
    return arr_no_nan[-1]-arr_no_nan[0]

def median(array):
    """Return the median of the array
    
    """
    sa = sorted(array)
    salen=len(sa)
    
    return (salen%2==0)*((sa[salen/2-1]+sa[salen/2])/2)+(salen%2!=0)*sa[salen/2]


def diff(array):
    """Return the differences between adjacent elements in an array"""
    
    if len(array)==0 or len(array)==1: return array
    d=[]
    for i in xrange(1,len(array)):
        d.append(array[i]-array[i-1])
    
    return d
    
def intersect_chroms(chroms1,chroms2):
    """Return the intersection of two chromosome sets"""
    
    chroms=set(chroms1).intersection(set(chroms2))
    chroms=list(chroms)
    chroms.sort()
    
    return chroms

def find_nearest(array,point,search_start):
    """Find the nearest value to the point in the array"""
    
    diff=abs(array[search_start]-point)
    length=len(array)
    for ix in xrange(search_start+1,length):
        newdiff=abs(array[ix]-point)
        if newdiff>=diff: return ix-1
        else: diff=newdiff
        
    return length-1

def where(value1,value2,ls):
    """Find where the value1 and value2 are in ls
    
    This function returns the slice of ls between value1 and value2. 
    Suppose that where returns (s1, s2) given value1, value2, and ls.
    
    The use can get the slice between value1 and value2 by ls[s1:s2].
    """
    
    length=len(ls)
    if length == 0: return 0, 0
    if value1 > ls[-1]: return length,length
    
    start = bisect.bisect_left(ls, value1)
    if start == length: 
        end = start
        return start, end
    
    end = bisect.bisect_right(ls[start:], value2) + start
#    for i in xrange(0, length):
#        if value1 <= ls[i]: 
#            break
#    start = i
#    
#    if value2 >= ls[-1]: return start, length
#    
#    for i in xrange(start, length):
#        if value2 < ls[i]:
#            break
#    end = i
        
    return start, end


def where2(value1, value2, ls, interval):
    """Find where the value1 and value2 are located in ls, where the interval between neighboring elements are given.
    
    This function may be faster than where, but might be slower in case that the interval is irregular.
    
    This function returns the slice of ls between value1 and value2. 
    Suppose that where returns (s1, s2) given value1, value2, and ls, and value1 <= value2
    
    The use can get the slice between value1 and value2 by ls[s1:s2].
    """
    
    length = len(ls)
    start = 0
    end = 0
    
    # if empty array, just return 0, 0
    if length ==0: return start, end
    
    diff1 = int(value1 - ls[0])
    if diff1 >= 0: 
        start = min(diff1/interval, length-1)
    
    # adjust 'start'    
    while start > 0 and value1 < ls[start]:
        start-=1
        
    while start < length and value1 > ls[start]:
        start+=1
    
    diff2 = int(value2 - value1)
    if diff2 >= 0:
        try:
            end = min(start + diff2/interval, length-1)
        except ZeroDivisionError:
            interval = 1
            end = min(start + diff2/interval, length-1)
    else:   # if value1 > value2, just return start, start
        return start, start
    
    # adjust 'end'
    while end > 0 and value2 < ls[end]:
        end-=1   
    
    while end < length and value2 >= ls[end]:
        end+=1
     
    end = min(length, end)
    
    return start, end   
    
def findbin(value,bins):
    """Return the bin where value belongs.
    
    Parameters:
    1. value: a value whose bin to know
    2. bins: an array of bins. bins[i] represents a bin of bins[i]<= <bins[i+1]. 
    
    Return
    i: if i==-1, value is less than bins[0]. If i==len(bins), value >=bins[-1].
    
    """
    
    if value < bins[0]: return -1
    elif value >= bins[-1]: return len(bins)
    else:
        i=0
        while bins[i]<value:
            i+=1
        return i-1
    

def argsort(L):
    """house-made argsort function not to import NumPy
    Return the sorted list as the first output argument and sorted indices as the second
    """
    
    index_element=[(ix,L[ix]) for ix in range(0,len(L))]
    index_element.sort(key=operator.itemgetter(1))
    sorted=[]
    sorted_ix=[]
    for ix,el in index_element:
        sorted.append(el)
        sorted_ix.append(ix)
        
    return sorted,sorted_ix

def linspace(s,e,n):
    """Get n points that are equi-spaced between s and e"""
    
    l = e-s
    return map(lambda i: s + 1.0 * l * i / (n-1), range(0,n))

def extend_list_series(nestlist):
    """Extend nested lists in lists"""
    
    series=[]
    for n in nestlist:
        series+=n

    return series


def bin(x,bins):
    """Do binning for x.
    
    Parameters:
    1. x: an array of data to do binning for
    2. bins: an array of bins. b[i-1]<= <b[i]
    
    Return:
    1. binned: an array of binned data. Each array element is the count within each bin
    """
    
#    binlen=len(bins)
#    binned=[0]*(binlen+1)
#    sx=sorted(x)
#    xlen=len(sx)
#    j=0
#    
#    #x < bins[0]
#    while j<xlen and sx[j]<bins[0]:
#        binned[0]+=1
#        j+=1
#    
#    # bins[i] <= x < bins[i+1]
#    for i in xrange(0,binlen-1):
#        while j<xlen and sx[j]>=bins[i] and sx[j]<bins[i+1]:
#            binned[i+1]+=1
#            j+=1
#    
#    # x>=bins[-1]
#    while j < xlen and sx[j]>=bins[-1]:
#        binned[-1]+=1
#        j+=1
#                 
#    return binned 

    binlen = len(bins)
    binned = [0]* (binlen-1)
    sx = sorted(x)
    xlen = len(sx)
    j=0
    
    # prefilter
    while j < xlen and sx[j]<bins[0]:
        j+=1
        
    for i in xrange(0, binlen-1):
        while j < xlen and sx[j] >= bins[i] and sx[j] < bins[i+1]:
            binned[i] += 1
            j+=1
    
    return binned


def binxy(bins, x, y, binfunc='mean', NaN=False, right=False):
    """Do binning on y in order of x.
        
    In general, x represents time or genomic coordinates and y is corresponding signal.
    Note that x is assumed to be sorted in ascending order; thus, this function is more 
    specific than 'bin' function.
        
    Parameters:
    1. bins: bins
    2. x: x (time or genomic corrdinates)
    3. y: y (a list of signal values at x points)
    4. binfunc: function that select the representative for each bin (eg, first: the first element in a bin, last: last element in a bin etc)
    4. NaN: True - Put float('nan') as a representative value if no values are in a bin.
    """
                
    length_bins=len(bins)
    if length_bins < 2: raise Exception('bins must have more than two elements')
    if len(x) != len(y): raise Exception('x and y must have the same length')
    if len(x) == 0: return []
    
    binned=[[] for i in bins[:-1]]
    xlength=len(x)
    if xlength == 0: return []
    
    j=0
    # do binning, x must be sorted from the smallest to largest
    if right:
        for i in xrange(0,length_bins-1):
            while j<xlength and x[j]>bins[i] and x[j]<=bins[i+1]:
                binned[i].append(y[j])
                j+=1
    else:
        for i in xrange(0,length_bins-1):
            while j<xlength and x[j]>=bins[i] and x[j]<bins[i+1]:
                binned[i].append(y[j])
                j+=1
    
    if NaN:
        if binfunc =='mean':
            binfunc = lambda xs: xs and 1.0*sum(xs)/len(xs) or float('nan')
        elif binfunc =='first':
            binfunc = lambda xs: xs and xs[0] or float('nan')
        elif binfunc =='last':
            binfunc = lambda xs: xs and xs[-1] or float('nan')
        elif binfunc == 'middle':
            binfunc = lambda xs: xs and xs[len(xs)/2] or float('nan')
        elif binfunc =='median':
            binfunc = lambda xs: xs and median(xs) or float('nan')
        elif binfunc =='min':
            binfunc = lambda xs: xs and min(xs) or float('nan')
        elif binfunc == 'max':
            binfunc = lambda xs: xs and max(xs) or float('nan')
        elif binfunc == 'raw':
            binfunc = lambda xs: xs
    else:
        if binfunc =='mean':
            binfunc = lambda xs: xs and 1.0*sum(xs)/len(xs) or 0.0
        elif binfunc =='first':
            binfunc = lambda xs: xs and xs[0] or 0.0
        elif binfunc == 'middle':
            binfunc = lambda xs: xs and xs[len(xs)/2] or 0.0
        elif binfunc =='last':
            binfunc = lambda xs: xs and xs[-1] or 0.0
        elif binfunc =='median':
            binfunc = lambda xs: xs and median(xs) or 0.0
        elif binfunc =='min':
            binfunc = lambda xs: xs and min(xs) or 0.0
        elif binfunc == 'max':
            binfunc = lambda xs: xs and max(xs) or 0.0
        elif binfunc == 'raw':
            binfunc = lambda xs: xs
     
    binned = map(binfunc, binned)
       
    return binned 

def binxy_equibin(bins, x, y, binfunc='mean', NaN=False):
    """Specialized bin function for equi-spaced bins"""
    
    
    length_bins = len(bins)    
    if length_bins < 2: raise Exception('bins must have more than two elements')
    if len(x) != len(y): raise Exception('x and y must have the same length')
    if len(x) == 0: return []
    
    binned = [[] for i in bins[:-1]]
    first = bins[0]
    last = bins[-1]
    interval = bins[1]-bins[0]
    for i, j in itertools.izip(x,y):
        if i< first: continue
        elif i >= last: break
        
        try:
            binned[int((i-first)/interval)].append(j)
        except ZeroDivisionError:
            return []
        
    if NaN:
        if binfunc =='mean':
            binfunc = lambda xs: xs and 1.0*sum(xs)/len(xs) or float('nan')
        elif binfunc =='first':
            binfunc = lambda xs: xs and xs[0] or float('nan')
        elif binfunc =='last':
            binfunc = lambda xs: xs and xs[-1] or float('nan')
        elif binfunc == 'middle':
            binfunc = lambda xs: xs and xs[len(xs)/2] or float('nan')
        elif binfunc =='median':
            binfunc = lambda xs: xs and median(xs) or float('nan')
        elif binfunc =='min':
            binfunc = lambda xs: xs and min(xs) or float('nan')
        elif binfunc == 'max':
            binfunc = lambda xs: xs and max(xs) or float('nan')
        elif binfunc == 'raw':
            binfunc = lambda xs: xs
    else:
        if binfunc =='mean':
            binfunc = lambda xs: xs and 1.0*sum(xs)/len(xs) or 0.0
        elif binfunc =='first':
            binfunc = lambda xs: xs and xs[0] or 0.0
        elif binfunc == 'middle':
            binfunc = lambda xs: xs and xs[len(xs)/2] or 0.0
        elif binfunc =='last':
            binfunc = lambda xs: xs and xs[-1] or 0.0
        elif binfunc =='median':
            binfunc = lambda xs: xs and median(xs) or 0.0
        elif binfunc =='min':
            binfunc = lambda xs: xs and min(xs) or 0.0
        elif binfunc == 'max':
            binfunc = lambda xs: xs and max(xs) or 0.0
        elif binfunc == 'raw':
            binfunc = lambda xs: xs
     
    binned = map(binfunc, binned)
       
    return binned 
        
def seq(fr,to,by):
    """An analogous function to 'seq' in R
    
    Parameters:
    1. fr: from
    2. to: to
    3. by: by (interval)
    """
    if fr<to:
        return range(fr,to+abs(by),abs(by))
    elif fr>to:
        if by>0:
            aseq = range(fr,to-by,-1*by)
        else:
            aseq = range(fr,to+by,by)
    else:
        aseq = [fr]
    
    if aseq[-1]>to: return aseq[:-1]
    else: return aseq
    
    
def lininterpol(first,second,x):
    """Perform linear interpolation for x between first and second.
    
    Parameters:
    1. first: [x1, y1], where x1 is the first x value and y1 is its y value
    2. second: [x2, y2], where x2 is the second x value and y2 is its y value
    3. x: the x value whose y value will be obtained.
    
    Return:
    y
    
    """
    
    x1,y1=first
    x2,y2=second
    
    a=x-x1
    b=x2-x
    
    try:
        y=1.0*(b*y1+a*y2)/(a+b)
    except ZeroDivisionError:
        y=y1
    
    return y


def ordinate(values,maxrange,levels):
    """Ordinate values given a maximum data range and number of levels
    
    Parameters:
    1. values: an array of continuous values to ordinate
    2. maxrange: the maximum data range. Values larger than this will be saturated.
    3. levels: the number of levels at which values are ordinated
    
    """
    
    quantizer=lambda dist,maxrange,levels: int(1.0*max(1,dist-1)*levels/maxrange)+1
    if type(values)==list or type(values)==tuple or type(values)==array:
        ordinated=[]
        for v in values:
            if v==0:
                ordinated.append(v)
            else:
                ordinated.append(quantizer(v,maxrange,levels))
        return ordinated
    else:
        if values==0:
            return values
        else:
            return quantizer(values,maxrange,levels)


def get_certain_part(things, percentlimit=[50,100]):
    """Get a certain portion (or part) from a set of numbers
    
    Parameters:
    1. things: a list of numbers to consider
    2. percentlimit: a list of a lower and upper limits
    
    """
    # if neither list nor array nor tuple, just return None
    if type(things)!=list and type(things)!=array and type(things)!=tuple: return None
    
    n_things = len(things)
    # return None if empty array
    if n_things == 0: return None
    
    sthings = sorted(list(things))
    
    lower = int(round(1.0*percentlimit[0]*n_things/100))-1
    upper = int(round(1.0*percentlimit[1]*n_things/100))
    
    return sthings[lower:upper]


def get_boundary_median(things, lower, upper):
    """Return the boundary values and median of a given percentage range
    
    Parameters:
    1. things: a list of numbers
    2. lower: lower percentage limit
    3. upper: upper percentage limit
    
    Return:
    lower, median, upper
    """
    
    # if neither list nor array nor tuple, just return None
    if type(things)!=list and type(things)!=array and type(things)!=tuple: return None, None, None
    
    n_things = len(things)
    
    if n_things == 0: return None, None, None
    
    sthings = sorted(list(things))
    
    l = int(round(1.0*lower*n_things/100))-1
    r = int(round(1.0*upper*n_things/100))
    
    return sthings[l], median(sthings[max(0, l):min(n_things, r+1)]), sthings[r]
    

def get_boundaries_medians(things, lowers=[], uppers=[]):
    """Return the boundaries and medians of given percentage ranges.
    
    Parameters:
    1. things: a list of numbers
    2. lowers: lower percentage limits
    3. uppers: upper percentage limits
    
    Returns:
    lower, median, upper
    """
    
    # if neither list nor array nor tuple, just return None
    if type(things)!=list and type(things)!=array and type(things)!=tuple: return [], [], []
    
    n_things = len(things)
    
    if n_things == 0: return [], [], []
    
    sthings = sorted(list(things))
    
    l = map(lambda x: int(round(1.0*x*n_things/100))-1, lowers)
    r = map(lambda x: int(round(1.0*x*n_things/100)), uppers)
            
    return map(lambda x: sthings[x], l), map(lambda x, y: median(sthings[max(0, x):min(n_things, y+1)]), l, r), map(lambda y: sthings[y], r)
    

def simply_fill_NaNs(arr):
    """fill NaNs in an array by interpolation
    
    Parameters:
    1. arr: input array to fill
    
    WARNING: if one NaN is enclosed by a NaN or NaNs, the value will be filled with NaN again.
    """
    
    filled = []
    l=len(arr)
    for i in xrange(l):
        if arr[i]!=arr[i]:
            if i==0: 
                filled.append(arr[1])
            elif i==l:
                filled.append(arr[-2])
            else:
                filled.append(sum([arr[i-1],arr[i+1]])/2)
        else:
            filled.append(arr[i])
                
    return filled

def scale_array(arr, s):
    """Scale an array by s
    
    Parameters:
    1. array: a numeric array or list
    2. s: scaling factor, real number
    """
    
    return [a*s for a in arr]

def randints(a, b, N=1):
    """Generate N random numbers between a and b
    
    Parameters:
    1. a: lower limit
    2. b: upper limit
    3. N: the number of random numbers
    """
    
    return [random.randint(a, b) for i in xrange(N)] 


def find_nearest_multiple(num, n):
    """Find the nearest multiple of n to num
    
    1. num: a number, num must be larger than n to have a meaningful result.
    2. n : a number to be multipled
    """
    
    return int(round(num / n)) * n
    
def sort_chroms(chroms):
    """Sort chromosomes by their true order, not by simple string order.
    
    Note that chroms must have a form of 'chrN', where N is alpha-numeric characters (eg, I, 1, 2L, or 1_random etc).
    In case that 'I', 'II', these must be changed into a standard format.
    
    """
        
    # split the chromosome names into several parts ('chr', number(s), 'alphabet(s)', '_random')
    startwchr = []
    startwochr = []
    for chrom in chroms:
        try:
            startwchr.append(re.search(r'(chr)([0-9]*)([A-Za-z]*)(_[_A-Za-z0-9]+)?', chrom).groups())
        except AttributeError: # when the chrom name does not start with 'chr'
            #startwochr.append(re.search(r'([A-Za-z]*)([0-9]*)([A-Za-z]*)(_[_A-Za-z0-9]+)?', chrom).groups())
            startwochr.append(chrom)
    
    # handle normal chromosomes
    if startwchr:
        chr_num_letter_tail = [s for s in startwchr if s[-1] and s[1]]
        chr_letter_tail = [s for s in startwchr if s[-1] and not s[1] and s[2]]
        chr_num_letter = [s for s in startwchr if not s[-1] and s[1]]
        chr_letter = [s for s in startwchr if not s[-1] and not s[1] and s[2]]
    
        chr_num_letter_tail = sorted(chr_num_letter_tail, key=operator.itemgetter(3))   # sort by tail
        chr_num_letter_tail = sorted(chr_num_letter_tail, key=operator.itemgetter(2))   # sort by alphabet (X, Y or L R etc)
        chr_num_letter_tail = sorted(map(lambda x: (x[0], int(x[1]), x[2], x[3]), chr_num_letter_tail), key=operator.itemgetter(1)) # sort by number
        chr_num_letter_tail = map(lambda x: [x[0], str(x[1]), x[2], x[3]], chr_num_letter_tail)
        chr_num_letter_tail = map(lambda x: ''.join(x), chr_num_letter_tail) # join the elements back to a single string
    
        chr_letter_tail = sorted(chr_letter_tail, key=operator.itemgetter(3))   # sort by tail
        chr_letter_tail = sorted(chr_letter_tail, key=operator.itemgetter(2))   # sort by alphabet (X, Y or L R etc)
        chr_letter_tail = map(lambda x: ''.join(x), chr_letter_tail)
    
        chr_num_letter = sorted(chr_num_letter, key=operator.itemgetter(2))   # sort by alphabet (X, Y or L R etc)
        chr_num_letter = sorted(map(lambda x: (x[0], int(x[1]), x[2], x[3]), chr_num_letter), key=operator.itemgetter(1)) # sort by number
        chr_num_letter = map(lambda x: [x[0], str(x[1]), x[2], x[3]], chr_num_letter)
        chr_num_letter = map(lambda x: ''.join(x[:3]), chr_num_letter) # join the elements back to a single string
    
        chr_letter = sorted(chr_letter, key=operator.itemgetter(3))   # sort by tail
        chr_letter = sorted(chr_letter, key=operator.itemgetter(2))   # sort by alphabet (X, Y or L R etc)
        chr_letter = map(lambda x: ''.join(x[:3]), chr_letter)
    
        startwchr = chr_num_letter + chr_letter + chr_num_letter_tail + chr_letter_tail
    
    # handle not normal chromosomes
    if startwochr:

        startwochr.sort()
    
    return startwchr + startwochr 
            
    
    
    

    
    
    
    
    

        
    
    
    