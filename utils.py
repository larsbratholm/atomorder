from __future__ import print_function
import __builtin__
import sys

def oprint(level, string)
    """
    Helper for printing to stdout

    Parameters
    ---------
    level: integer
        print level of the string (0: None, 1: results, 3: +progress, 4+: excess)
    string: string
        string to be printed

    """

    if args.print_level >= level:
        __builtin__.print(string)
        #__builtin__.print(*args, **kwargs, file = sys.stderr)

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

    if args.print_level >= level:
        __builtin__.print(string, file = sys.stderr)
