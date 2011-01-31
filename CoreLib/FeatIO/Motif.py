# Time-stamp: <2008-07-24 02:58:09 Tao Liu>

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

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

class Motif:
    """Motif class.
    
    """
    def __init__ (self,name="",pwm=None,filehandle=None):
        self.name = name
        self.pwd = pwd
        if self.pwd:
            self.__parse_pwm__()
        pass

    def __init_from_text__ (self, text=""):
        pass

    def __init_from_filehandle__ (self, fhd = None ):
        pass

    def __parse_pwm__ (self):
        self.A = None
        self.T = None
        self.C = None
        self.G = None

# ------------------------------------
# Main function
# ------------------------------------
def test():
    #test function
    pass

if __name__ == '__main__':
    test()
