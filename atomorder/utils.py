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
    parser.add_argument('--print-level', help='Print-level -  0: quiet, 1: results and errors, 2: +warnings, 3: +progress, 4+: excess', action='store', default=1, type=int) 
    parser.add_argument('-f', '--format', help='File format', type=str, action='store', default='guess', choices=["guess","xyz","pdb"])
    parser.add_argument('--simple-rotate', help='Ignore bond order, align a single reactant and product molecule and match all atoms', action='store_true', default=False)
    # TODO output atom mapping oneline, save reordered products
    # TODO parameter object
    # TODO allow possibility to give pickle with reaction object
    # TODO output sybyl

    args = parser.parse_args()

    # override setting defaults
    settings.update(args)

def get_distance(x, y, axis):
    return np.sum((x-y)**2, axis=axis)**0.5

def read_coordinates(filename, format_ = "guess"):
    """
    Parses coordinates file

    Parameters
    ---------
    filename: string
        file name
    format_: string
        file format

    """
    # Guess from extension if file format is not given.
    if format_ == "guess":
        format_ = filename.split(".")[-1]

    if format_ == "xyz":
        return read_coordinates_xyz(filename)
    elif self._format == "pdb":
        return read_coordinates_pdb(filename)
    else:
        quit("Error: Unknown file format %s" % format_)

def read_coordinates_xyz(filename):
    """
    Read and parse xyz coordinate file

    Parameters
    ----------
    filename: string
        Coordinates filename

    Returns
    -------
    atoms: array-like
        N-sized array of element symbols, e.g. ['H','C']
    coordinates: array-like
        Nx3-sized array of euclidian coordinates

    """
    with open(filename) as f:
        lines = f.readlines()
        try:
            n_atoms = int(lines[0])
        except ValueError:
            quit("Error reading XYZ-file. Expected no. of atoms in first line \
                  of filefile. Got %s" % lines[0])
        atoms = np.empty(n_atoms, dtype = np.str)
        coordinates = np.zeros((n_atoms, 3), dtype = float)
        for l, line in enumerate(lines[2:]):
            try:
                tokens = line.split()
                atoms[l] = tokens[0]
                coordinates[l] = np.asarray(tokens[1:4], dtype=float)
            except IndexError, ValueError:
                quit("Error reading line %d in inputfile %s: \n line %s" % (l+3, filename, line))
    return atoms, coordinates

def get_coordinates_pdb(filename):
    """
    get_coordinates_pdb(filename)

    Get coordinates from the first chain in a pdb file

    Parameters
    ----------
    filename: string
        PDB3 filename

    Returns
    -------
    atoms: array-like
        N-sized array of element_symbols, e.g. ['H','C']
    coordinates: array-like
        Nx3-sized array of euclidian coordinates

    """
    # PDB files tend to be a bit of a mess. The x, y and z coordinates
    # are supposed to be in column 31-38, 39-46 and 47-54, but this is not always the case.
    # Because of this the three first columns containing a decimal is used.
    # Since the format doesn't require a space between columns, we use the above
    # column indices as a fallback.
    x_column = None
    coordinates = []
    # Same with atoms and atom naming. The most robust way to do this is probably
    # to assume that the element symbol can be inferred from the element symbol given in column 3.
    atoms = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("TER") or line.startswith("END"):
                break
            if line.startswith("ATOM"):
                tokens = line.split()
                # Try to get the element symbol
                try:
                    atom = tokens[2][0]
                    if atom in ["H", "C", "N", "O", "S", "P"]:
                        atoms.append(atom)
                    else:
                        atom = tokens[2][1]
                        if atom in ["H", "C", "N", "O", "S", "P"]:
                            atoms.append(atom)
                        else:
                            exit("Error parsing element symbol for the following line: \n%s" % line)
                except IndexError:
                        exit("Error parsing element symbol for the following line: \n%s" % line)

                if x_column == None:
                    try:
                        # look for x column
                        for i, x in enumerate(tokens):
                            if "." in x and "." in tokens[i+1] and "." in tokens[i+2]:
                                x_column = i
                                break
                    except IndexError:
                        exit("Error parsing coordinates for the following line: \n%s" % line)
                # Try to read the coordinates
                try:
                    coordinates.append(np.asarray(tokens[x_column:x_column+3],dtype=float))
                except:
                    # If that doesn't work, use hardcoded indices
                    try:
                        x = line[30:38]
                        y = line[38:46]
                        z = line[46:54]
                        coordinates.append(np.asarray([x,y,z],dtype=float))
                    except IndexError, TypeError:
                        exit("Error parsing input for the following line: \n%s" % line)

    coordinates = np.asarray(V)
    atoms = np.asarray(atoms)
    if coordinates.shape[0] != atoms.size:
        error("Mismatch in number of parsed element symbols (%d) and number of parsed coordinates (%d) from PDB file: %s" \
                % (coordinates.shape[0], atoms.size, filename))
    return atoms, coordinates

def read_coordinates(filename, format_ = "guess"):
    """
    Parses coordinates file

    Parameters
    ---------
    filename: string
        file name
    format_: string
        file format

    """
    # Guess from extension if file format is not given.
    if format_ == "guess":
        format_ = filename.split(".")[-1]

    if format_ == "xyz":
        return read_coordinates_xyz(filename)
    elif self._format == "pdb":
        return read_coordinates_pdb(filename)
    else:
        quit("Error: Unknown file format %s" % format_)

def read_coordinates_xyz(filename):
    """
    Read and parse xyz coordinate file

    Parameters
    ----------
    filename: string
        Coordinates filename

    Returns
    -------
    atoms: array-like
        N-sized array of element symbols, e.g. ['H','C']
    coordinates: array-like
        Nx3-sized array of euclidian coordinates

    """
    with open(filename) as f:
        lines = f.readlines()
        try:
            n_atoms = int(lines[0])
        except ValueError:
            quit("Error reading XYZ-file. Expected no. of atoms in first line \
                  of filefile. Got %s" % lines[0])
        atoms = np.empty(n_atoms, dtype = np.str)
        coordinates = np.zeros((n_atoms, 3), dtype = float)
        for l, line in enumerate(lines[2:]):
            try:
                tokens = line.split()
                atoms[l] = tokens[0]
                coordinates[l] = np.asarray(tokens[1:4], dtype=float)
            except IndexError, ValueError:
                quit("Error reading line %d in inputfile %s: \n line %s" % (l+3, filename, line))
    return atoms, coordinates

def get_coordinates_pdb(filename):
    """
    get_coordinates_pdb(filename)

    Get coordinates from the first chain in a pdb file

    Parameters
    ----------
    filename: string
        PDB3 filename

    Returns
    -------
    atoms: array-like
        N-sized array of element_symbols, e.g. ['H','C']
    coordinates: array-like
        Nx3-sized array of euclidian coordinates

    """
    # PDB files tend to be a bit of a mess. The x, y and z coordinates
    # are supposed to be in column 31-38, 39-46 and 47-54, but this is not always the case.
    # Because of this the three first columns containing a decimal is used.
    # Since the format doesn't require a space between columns, we use the above
    # column indices as a fallback.
    x_column = None
    coordinates = []
    # Same with atoms and atom naming. The most robust way to do this is probably
    # to assume that the element symbol can be inferred from the element symbol given in column 3.
    atoms = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("TER") or line.startswith("END"):
                break
            if line.startswith("ATOM"):
                tokens = line.split()
                # Try to get the element symbol
                try:
                    atom = tokens[2][0]
                    if atom in ["H", "C", "N", "O", "S", "P"]:
                        atoms.append(atom)
                    else:
                        atom = tokens[2][1]
                        if atom in ["H", "C", "N", "O", "S", "P"]:
                            atoms.append(atom)
                        else:
                            exit("Error parsing element symbol for the following line: \n%s" % line)
                except IndexError:
                        exit("Error parsing element symbol for the following line: \n%s" % line)

                if x_column == None:
                    try:
                        # look for x column
                        for i, x in enumerate(tokens):
                            if "." in x and "." in tokens[i+1] and "." in tokens[i+2]:
                                x_column = i
                                break
                    except IndexError:
                        exit("Error parsing coordinates for the following line: \n%s" % line)
                # Try to read the coordinates
                try:
                    coordinates.append(np.asarray(tokens[x_column:x_column+3],dtype=float))
                except:
                    # If that doesn't work, use hardcoded indices
                    try:
                        x = line[30:38]
                        y = line[38:46]
                        z = line[46:54]
                        coordinates.append(np.asarray([x,y,z],dtype=float))
                    except IndexError, TypeError:
                        exit("Error parsing input for the following line: \n%s" % line)

    coordinates = np.asarray(V)
    atoms = np.asarray(atoms)
    if coordinates.shape[0] != atoms.size:
        error("Mismatch in number of parsed element symbols (%d) and number of parsed coordinates (%d) from PDB file: %s" \
                % (coordinates.shape[0], atoms.size, filename))
    return atoms, coordinates
