#!/usr/bin/env python2

import . utils
from .utils import oprint, eprint
import numpy as np
import sys
import .definitions
#import re
#import itertools
#import scipy.sparse
#import scipy.sparse.linalg
#import scipy.linalg

def parse_args():
    """
    Parses command line arguments.

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

    return parser.parse_args()

class reaction(object):
    """
    reaction(args)

    Everything involving both reactants and products


    Attributes
    ----------
    reactants: object
        Mixture object.
    products: object
        Mixture object.

    """

    # TODO construct the objective

    def __init__(self):
        if args.print_level == 3:
            print "Creating reactants"
        self.reactants = mixture(args.reactants)
        if args.print_level == 3:
            print "Creating products"
        self.products = mixture(args.products)

class mixture(object):
     """
    mixture(filenames)

    Collection of one or several molecules.

    Parameters
    ----------
    filenames: list
        Coordinate filenames.

    Attributes
    ----------
    filenames: list
        Coordinate filenames.
    molecules: list
        Molecule objects.
    """

    def __init__(self, filenames):
        self.filenames = filenames
        self.molecules = self.make_molecules()

    # A molecule object will contain everything in a single coordinate file
    def make_molecules(self):
        """
        make_molecules()

        Creates list of all given molecules.

        """

        molecules = []
        for filename in self.filenames:
            if args.print_level == 3:
                print "Creating molecule from %s" % filename
            mol = molecule(filename)
            quit()
            # If there are disconnected atoms in the molecule
            # split the molecules to separate objects
            mols = molecule.split()
            molecules.extend(mols)
        return molecules

class molecule(object):
    """
    molecule(filename)

    All atoms from a single coordinate file.

    Parameters
    ----------
    filename: string
        Coordinate filename

    Atributes
    ---------
    filename: string
        Coordinate filename
    atomnames: array_like
        N-size array of names of the molecule atoms
        e.g. ['H', 'C']
    coordinates: array_like
        Nx3-size array of euclidian coordinates of the molecule atoms
    size: integer
        Number of atoms in the molecule
    distance_matrix: array_like
        NxN-size array with inter-atom distances.
    monovalent: array_like
        Indices of monovalent atoms

    """

    def __init__(self, filename):
        self.filename = filename
        self.atomnames, self.coordinates = read_coordinates(self.filename)
        self.size = self.atomnames.size
        self.distance_matrix = self.get_distance(self.coordinates[None,:,:],self.coordinates[:,None,:])
        self.monovalent = np.in1d(self.atomnames,definitions.monovalent,invert=False)
        self.atoms = self.make_atoms()
        self.find_bonds()
        quit()
        #self.get_bond_details()
        #self.make_graph()
        # TODO make bond2

    def split(self):
        """
        split()

        Splits the molecule into several molecule objects if two or more
        disconnected regions are present (e.g. ligands or several reactants
        in the same file)

        Returns
        -------
        molecules: tuple
            one or more molecule objects
        """

        # TODO

    def get_distance(self, x, y, axis = -1):
        return np.sum((x-y)**2, axis=axis)**0.5

    def make_atoms(self):
        """
        make_atoms()

        Make array of atom objects

        Returns
        -------
        atoms: array-like
            atom objects

        """
        atoms = np.empty(self.size, dtype=object)
        for i in range(self.size):
            atoms[i] = atom(i,self)
        return atoms

    def find_bonds(self):
        """
        find_bonds()

        Determine which atoms form bonds.
        Credible bond lengths, bond angles and out of plane bending
        are taken from the CCDC 2016 database.

        """

        # To avoid redundancy, only analyze unique pairs
        pair_indices = np.triu_indices(self.size, 1)
        # ignore monovalent pairs
        multivalent = ~(self.monovalent[pair_indices[0]] & self.monovalent[pair_indices[1]])
        pair_indices = pair_indices[0][multivalent], pair_indices[1][multivalent]

        pair_distances = self.distance_matrix[pair_indices]

        # create arrays of typical bond length limits (size 4 x n_pairs)
        limits = np.asarray([definitions.bond_length_limits[(self.atomnames[i],self.atomnames[j])] for i,j in zip(*pair_indices)]
                            , dtype=float).T

        # bonds within the 'usual' range
        credible_bonds = (pair_distances >= limits[1]) & (pair_distances <= limits[2])

        # bonds outside the usual range but within a looser restricted range
        less_credible_bonds = (~credible_bonds) & (pair_distances >= limits[0]) & (pair_distances <= limits[3])

        credible_bond_indices = np.asarray(zip(*pair_indices))[credible_bonds]
        less_credible_bond_indices = np.asarray(zip(*pair_indices))[less_credible_bonds]

        # Associate the bonds with the involved atoms
        for index1, index2 in credible_bond_indices:
            self.atoms[index1].add_bond(index2)
            self.atoms[index2].add_bond(index1)

        # Check that all atoms have a realistic set of bonds
        self.check_bonds(less_credible_bond_indices)
        quit()

    def check_bonds(bond_indices, subset = None):
        if subset == None:
            subset = np.arange(self.size)

        # check all atoms, or a subset of atoms
        for atom in self.atoms[subset]:
            n_bonds = len(atom.bonded_atom_indices)
            status = atom.check_planarity(n_bonds) #TODO both too many and too few?
            status = False # TODO remove this
            if status:
                # flag the atom as validated
                atom.validate()
                continue
            else:
                # check if any bonds can be formed outside the usual
                # range, that can make the number of bonds and bond
                # angles consistent

                # go through any possible bond that the current atom
                # can be involved
                for index1, index2 in bond_indices:
                    # only need to check index1, as index1 is always less than index2
                    if index1 == atom.index and self.atoms[index2].validated == False:

class atom(object):
    """
    atom(index, molecule)

    Parameters
    ----------
    index: integer
        index of atom in the coordinate file
    molecule: object
        molecule object

    Attributes
    ----------
    atype: string
        atomtype, e.g. 'H' or 'Cl'
    molecule_index: integer
        index of atom in the molecule
    mixture_index: integer
        index of atom in the mixture
    coordinates: array-like
        3-sized array of euclidian coordinates
    distances: array-like
        N-sized array of euclidean distances between the atom and all
        atoms in the molecule.
    bonded_atom_indices: list
        Indices of atoms in the molecule that possibly form a covalent bond
        with the atom. The list is pruned to be consistent if an excess
        number of atoms are within the bond length cut-off.
    validated: bool
        Wether the atom has passed consistency checks on bond lengths and angles

    """

    def __init__(self, index, molecule):
        self.atomtype = molecule.atomnames[index]
        self.molecule_index = index
        self.mixture_index = None
        self.coordinates = molecule.coordinates[index]
        self.distances = molecule.distance_matrix[index]
        self.bonded_atom_indices = []
        self.validated = False

    def set_mixture_index(self, index):
        self.mixture_index = index

    def set_molecule_index(self, index):
        self.molecule_index = index

    def add_bond(self, atom):
        self.bonded_atom_indices.append(atom)

    def validate(self):
        self.validated = True

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
        N-sized array of atomtypes, e.g. ['H','C']
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
        N-sized array of atomtypes, e.g. ['H','C']
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
    # to assume that the atomtype can be inferred from the atomname given in column 3.
    atoms = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("TER") or line.startswith("END"):
                break
            if line.startswith("ATOM"):
                tokens = line.split()
                # Try to get the atomtype
                try:
                    atom = tokens[2][0]
                    if atom in ["H", "C", "N", "O", "S", "P"]:
                        atoms.append(atom)
                    else:
                        atom = tokens[2][1]
                        if atom in ["H", "C", "N", "O", "S", "P"]:
                            atoms.append(atom)
                        else:
                            exit("Error parsing atomtype for the following line: \n%s" % line)
                except IndexError:
                        exit("Error parsing atomtype for the following line: \n%s" % line)

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
        error("Mismatch in number of parsed atomtypes (%d) and number of parsed coordinates (%d) from PDB file: %s" \
                % (coordinates.shape[0], atoms.size, filename))
    return atoms, coordinates

if __name__ == "__main__":

    import argparse

    # global
    args = parse_args()
    quit()
    r = reaction()
    # initialize optimization
    # begin optimization
    # post processing
