from __future__ import print_function
import __builtin__
import sys
import argparse
from . import settings

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

def get_csd_atom_type(atom):
    """
    CSD non-matched (3d) deterministic sybyl atom type matching
    from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
    Comments refer to the definitions given from the above site

    Parameters
    ----------
    atom: object
        atom object

    Returns
    -------
    sybyl: string
        sybyl atom type

    """

    # TODO name->element_symbol

    # 2.2 If element_symbol is D then atom_type is H
    if atom.element_symbol == "D":
        sybyl = "H"
    # 2.3 If element_symbol is P then atom_type is P.3
    elif atom.element_symbol == "P":
        sybyl = "P.3"
    # 2.5 If element_symbol is C then
    elif atom.element_symbol == "C":
        sybyl = get_csd_atom_type_C(atom)
    # 2.6 If element_symbol is O then
    elif atom.element_symbol == "O":
        sybyl = get_csd_atom_type_O(atom)
    # 2.7 If element_symbol is N then
    elif atom.element_symbol == "N":
        sybyl = get_csd_atom_type_N(atom)
    # 2.8 If element_symbol is S then
    elif atom.element_symbol == "S":
        sybyl = get_csd_atom_type_S(atom)
    # 2.10 If element_symbol is none of the above then atom_type is element_symbol
    else:
        sybyl = atom.element_symbol

    return sybyl

def get_csd_atom_type_C(atom):
    """
    CSD non-matched (3d) deterministic sybyl atom type matching
    from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

    Parameters
    ----------
    atom: object
        atom object
    """

    # 2.5.1 If num_bond .ge. 4 then atom_type is C.3
    if atom.num_bond >= 4: return "C.3"
    # 2.5.2 If num_bond .eq. 1 then calculate bond_distance
    if atom.num_bond == 1:
        # TODO calc bond distance
        bond_distance = None

        # 2.5.2.1 If bond_distance .gt. 1.41A then atom_type is C.3
        if bond_distance > 1.41: return "C.3"
        # 2.5.2.2 If bond_distance .le. 1.22A then atom_type is C.1
        if bond_distance <= 1.22: return "C.1"
        # 2.5.2.3 If bond_distance is none of the above then atom_type is C.2
        return "C.2"
    # 2.5.3 If element_symbol is C and none of the above then calculate average_angle about C
    # TODO calc average angle
    average_angle = None
    # 2.5.3.1 If average_angle .le. 115 deg then atom_type is C.3
    if average_angle <= 115: return "C.3"
    # 2.5.3.2 If average_angle .gt. 160 deg then atom_type is C.1
    if average_angle > 160: return "C.1"
    # 2.5.3.3 If average_angle is none of the above then atom_type is C.2
    return "C.2"

def get_csd_atom_type_O(atom):
    """
    CSD non-matched (3d) deterministic sybyl atom type matching
    from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

    Parameters
    ----------
    atom: object
        atom object

    """

    # 2.6.1 If num_nonmet .eq. 1 then
    if atom.num_bond == 1:
        # TODO get bonding partner this

        # 2.6.1.1 If bond is to carbon .AND. carbon forms a total of 3
        #         bonds, 2 of which are to an oxygen forming only 1
        #         non-metal bond then atom_type is O.co2
        if this.element_symbol == "C" and this.num_bonds == 3:
            for this_neighbor in this.bonded: # TODO ?
                # TODO ?
                if this_neighbor.element_symbol == "O" and this_neighbor.id != this.id \
                        and this_neighbor.num_bond == 1:
                            return "O.co2"
        # 2.6.1.2 If bond is to phosphorus .AND. phosphorous forms at 
        #         least 2 bonds to an oxygen forming only 1 non-metal 
        #         bond then atom_type is O.co2
        if this.element_symbol == "P":
            for this_neighbor in this.bonded:
                if this_neighbor.element_symbol == "O" and this_neighbor.id != this.id \
                        and this_neighbor.num_bond == 1:
                            return "O.co2"
    # 2.6.3 If num_bond .ge. 2 then atom_type is O.3
    if atom.num_bond >= 2: return "O.3"

    # 2.6.4 If element_symbol is O and none of the above then atom_type is O.2
    return "O.2"

def get_csd_atom_type_N(atom):
    """
    CSD non-matched (3d) deterministic sybyl atom type matching
    from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

    Parameters
    ----------
    atom: object
        atom object

    """

    # 2.7.1 If num_nonmet .eq. 4 then atom_type is N.4
    if atom.num_bond == 4: return "N.4"

    # 2.7.2 If num_nonmet .eq. 1 then calculate bond_distance
    if atom.num_bond == 1:
        # TODO calc bond distance
        bond_distance = None

        # 2.7.2.1 If bond_distance .le. 1.2A then atom_type is N.1
        if bond_distance >= 1.2: return "N.1"

        # 2.7.2.2 If bond_distance .gt. 1.2A then atom_type is N.3
        return "N.3"

    # 2.7.3 If num_nonmet .eq. 3 .AND. one bond is to C--O or C--S then atom_type is N.am
    if atom.num_bond == 3:
        for this in atom.bonded:
            if this.element_symbol == "C":
                for this_neighbor in this.bonded:
                    # At this stage it's not known if a bond is double or single
                    # so just check that O and S only bonds to the parent carbon
                    if this_neighbor.element_symbol in ("O","S") and \
                            this_neighbor.num_bond == 1:
                                return "N.am"
            # 2.7.4 If num_nonmet .eq. 3 otherwise then calculate sum_of_angles around N
            average_angle = None # TODO

            # 2.7.4.1 If sum_of_angles .ge. 350 deg then atom_type is N.pl3
            if average_angle >= 350: return "N.pl3"
            # 2.7.4.2 If sum_of_angles .lt. 350 deg then atom_type is N.3
            return "N.3"

    # 2.7.5 If element_symbol is N and none of the above then calculate average_angle about N
    # check if average_angles is already calculated
    if average_angles not in locals():
        average_angle = None # TODO

    # 2.7.5.1 If average_angle .gt. 160 deg then atom_type is N.1
    if average_angle > 160: return "N.1"
    # 2.7.5.2 If average_angle .le. 160 deg then atom_type is N.2
    return "N.2"

def get_csd_atom_type_S(atom):
    """
    CSD non-matched (3d) deterministic sybyl atom type matching
    from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

    Parameters
    ----------
    atom: object
        atom object

    """

    # 2.8.1 If num_nonmet .eq. 3 .AND. 1 bond is to an oxygen with only one
    #       non-metal bond then atom_type is S.o
    if atom.num_bond == 3:
        for this in atom.bonds:
            if this.element_symbol == "O" and this.num_bond == 1: return "S.o"
    # 2.8.2  If num_nonmet .eq. 4 .AND. 2 bonds are to an oxygen with only 
    #        one non-metal bond then atom_type is S.o2
    if atom.num_bond == 4:
        count = 0
        for this in atom.bonds:
            if this.element_symbol == "O" and this.num_bond == 1: count += 1
            if count == 2: return "S.o2"
    # 2.8.3 If num_bond .ge. 2 then atom_type is S.3
    if atom.num_bond >= 2: return "S.3"
    # 2.8.4 If element_symbol is S and none of the above then atom_type is S.2
    return "S.2"

def parse_args():
    """
    Parses command line arguments and overwrites setting defaults

    """

    description = ""
    epilog = ""

    parser = argparse.ArgumentParser(
            description = description,
            formatter_class = argparse.RawDescriptionHelpFormatter,
            epilog = epilog)

    parser = argparse.ArgumentParser(description='Fit probability density functions to data-files')
    parser.add_argument('-r', '--reactants', help='Reactant structures in a coordinate file format.', action='store', type=str, nargs='+')
    parser.add_argument('-p', '--products', help='Product structures in a coordinate file format.', action='store', type=str, nargs='+')
    parser.add_argument('--print-level', help='Print-level.', type=int, action='store', default=1) # 0: quiet, 1: results and errors, 2: +warnings, 3: +progress, 4+: excess
    parser.add_argument('-f', '--format', help='File format', type=str, action='store', default='guess', choices=["guess","xyz","pdb"])
    # TODO output atom mapping oneline, save reordered products
    # TODO parameter object

    args = parser.parse_args()

    # override setting defaults
    settings.update(args)

