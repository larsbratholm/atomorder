"""
atomorder/utils.py

Utility functions that are not object specific

"""

from __future__ import print_function
import __builtin__
import sys
import argparse
from . import settings
import numpy as np

def oprint(level, string):
    """
    Helper for printing to stdout

    Parameters
    ---------
    level: integer
        print level of the string (0: None, 1: results, 3: +progress, 4+: excess)
    string: string
        string to be printed

    """

    if settings.print_level >= level:
        __builtin__.print(string)

def eprint(level, string):
    """ 
    Helper for printing to stderr

    Parameters
    ---------
    level: integer
        print level of the string (0: None, 1: errors, 2: +warnings, 4+: excess)
    string: string
        string to be printed

    """

    if settings.print_level >= level:
        __builtin__.print(string, file = sys.stderr)

def get_distance(x, y, axis):
    return np.sum((x-y)**2, axis=axis)**0.5

# http://stackoverflow.com/a/13849249
def vector_angle(v1,v2):
    """
    Angle between two normed vectors
    """
    return np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))

