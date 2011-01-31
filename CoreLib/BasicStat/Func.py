# Time-stamp: <2010-12-14 21:06:20 Tao Liu>

"""Module Description

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

from math import sqrt
# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
def mean (nums):
    """Calculate Mean.

    Parameters:
    nums:  list of numbers
    Return Value:
    mean value
    """
    return float(sum(nums))/len(nums)


def median (nums):
    """Calculate Median.

    Parameters:
    nums:  list of numbers
    Return Value:
    median value
    """
    p = sorted(nums)
    l = len(p)
    if l%2 == 0:
        return (p[l/2]+p[l/2-1])/2
    else:
        return p[l/2]


def std (nums,n_mean=None):
    """Calculate unbias standard deviation.

    Parameters:
    nums:  list of numbers
    Return Value:
    std value
    """
    if not n_mean:
        n_mean = mean(nums)
    n = len(nums)
    if n == 1:
        return 0.0
    variance = 0.0
    for i in xrange(n):
        tmp = (nums[i]-n_mean)
        variance += (tmp*tmp)
        
    variance /= n-1
    return sqrt(variance)

def normalize (nums):
    """Normalize given numbers.

    Parameters:
    nums:  list of numbers
    Return Value:
    normalized numbers
    """
    n_mean = mean(nums)
    n_std  = std(nums,n_mean=n_mean)
    for i in xrange(len(nums)):
        nums[i] = float(nums[i]-n_mean)/n_std
    return nums

def centering (nums):
    """center given numbers such that the mean be 0.

    Parameters:
    nums:  list of numbers
    Return Value:
    centered numbers
    """
    n_mean = mean(nums)
    for i in xrange(len(nums)):
        nums[i] = float(nums[i]-n_mean)
    return nums

def cor (list1,list2):
    p1 = []
    p2 = []
    p1.extend(list1)
    p2.extend(list2)
    p1 = centering(p1)
    p2 = centering(p2)

def ciw_95_normal (list1):
    """Return the upper/lower width of 95% CI for a normal distributed data set.
    
    """
    sd = std(list1)
    return 1.96*sd/sqrt(len(list1))


# ------------------------------------
# Classes
# ------------------------------------
