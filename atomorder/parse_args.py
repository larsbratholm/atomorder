"""
atomorder/parse_args.py

Parses command line arguments and overwrites setting defaults

"""
from . import settings
import argparse
import sys

description = ""
epilog = ""

parser = argparse.ArgumentParser(
        description = description,
        formatter_class = argparse.RawDescriptionHelpFormatter,
        epilog = epilog)

parser = argparse.ArgumentParser(description='Fit probability density functions to data-files')
parser.add_argument('-r', '--reactants', help='Reactant structures in a coordinate file format.', action='store', type=str, nargs='+')
parser.add_argument('-p', '--products', help='Product structures in a coordinate file format.', action='store', type=str, nargs='+')
parser.add_argument('--print-level', help='Print-level -  0: quiet, 1: results and errors, 2: +warnings, 3: +progress, 4: excess, 5: EXTREME',
                                     action='store', choices = range(0,6), default=1, type=int)
parser.add_argument('-f', '--format', help='File format', type=str, action='store', default='guess', choices=["guess","xyz","pdb"])
parser.add_argument('-m', '--method', help='Method to use.\n \
                                            rotate: Ignore bond order, align a single reactant and product molecule and match all atoms\n \
                                            no-bond: Atom matching by rotation and atomic similarity\n \
                                            full: Atom matching by rotation and bond similarity\n \
                                            info: Information about molecule sybyl atom types, bond types and conjugated sub systems',
                                            choices = ['rotate', 'full', 'info', 'no-bond'], action='store', default='full')
parser.add_argument('-o', '--output', help='Given a filename, output the reordered product in xyz format instead of printing to stdout', action='store', type=str, default=sys.stdout)
parser.add_argument('--atomic-sybyl-weight', action='store', default=1, type=float)
parser.add_argument('--bond-weight', action='store', default=1, type=float)
# TODO output to folder
# TODO output atom mapping oneline, save reordered products
# TODO allow possibility to give pickle with reaction object
# TODO output sybyl
# TODO batch reactions
# TODO output aromatic/conjugated subgroups

args = parser.parse_args()

# override setting defaults
settings.update(args)
