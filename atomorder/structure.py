import numpy as np
import itertools
import warnings
from . import settings, constants, utils
from .utils import oprint, eprint

class Reaction(object):
    """
    Reaction()

    Everything involving both reactants and products


    Attributes
    ----------
    reactants: object
        Mixture object.
    products: object
        Mixture object.

    """

    def __init__(self):
        oprint(3, "Creating reactants")
        self.reactants = Mixture(settings.reactant_filenames)
        oprint(3, "Creating products")
        self.products = Mixture(settings.product_filenames)

class Mixture(object):
    """
    Mixture(filenames)

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
            if settings.print_level == 3:
                oprint(3, "Creating molecule from %s" % filename)
            mol = Molecule(filename)
            quit()
            # If there are disconnected atoms in the molecule
            # split the molecules to separate objects
            mols = mol.split()
            molecules.extend(mols)
        return molecules

class Molecule(object):
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
    element_symbols: array_like
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
        self.element_symbols, self.coordinates = read_coordinates(self.filename)
        self.size = self.element_symbols.size
        self.distance_matrix = self.get_distance_matrix()
        self.monovalent = self.get_monovalent_indices()
        self.atoms = self.make_atoms()
        self.find_bonds()
        self.get_sybyl_types()
        quit()
        #TODO sybyl types
        #TODO single/double/triple bond? - iterative from sybyl
        # TODO make bond2
        # TODO split
        #self.make_graph() - move?

    def get_sybyl_types(self):
        """
        Set sybyl atom types.

        """

        for atom in self.atoms:
            atom.get_sybyl_type()


    def get_distance_matrix(self):
        return utils.get_distance(self.coordinates[None,:,:],self.coordinates[:,None,:], axis = -1)

    def get_monovalent_indices(self):
        return np.in1d(self.element_symbols, constants.monovalent, invert = False)

    # TODO
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



    def make_atoms(self):
        """
        Make array of atom objects

        Returns
        -------
        atoms: array-like
            atom objects

        """
        atoms = np.empty(self.size, dtype=object)
        for i in xrange(self.size):
            atoms[i] = Atom(i,self)
        return atoms

    def find_bonds(self):
        """
        Determine which atoms form bonds.
        Credible bond lengths, bond angles and out of plane bending
        are taken from the CCDC 2016 database.

        """

        oprint(3, "Locating bonded atom pairs")

        # To avoid redundancy, only analyze unique pairs
        pair_indices = np.triu_indices(self.size, 1)
        # ignore monovalent pairs
        multivalent = ~(self.monovalent[pair_indices[0]] & self.monovalent[pair_indices[1]])
        pair_indices = pair_indices[0][multivalent], pair_indices[1][multivalent]

        pair_distances = self.distance_matrix[pair_indices]

        # create arrays of typical bond length limits (size 4 x n_pairs)
        limits = np.asarray([constants.bond_length_limits[(self.element_symbols[i],self.element_symbols[j])] for i,j in zip(*pair_indices)]
                            , dtype=float).T

        # bonds within the 'usual' range
        credible_bonds = (pair_distances >= limits[1]) & (pair_distances <= limits[2])

        # NOTE: could probably do without the less_credible_bonds range, but shouldn't be a bottleneck of any sorts
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

    def check_bonds(self, bond_indices, subset = None):
        """
        Checks to assert that the bonds are not abnormal.
        If any bonds are missing for a physically correct structure,
        then check within the looser credible range interval.

        Parameters:
        -----------
        bond_indices: tuple
            Indices outside the credible bond length range
            but within the less credible (but still realistic) range.

        """

        # Indices to be checked in the next iteration
        next_subset = []

        if subset == None:
            subset = np.arange(self.size)

        # TODO redo when better criteria is added
        # check all atoms, or a subset of atoms
        for atom in self.atoms[subset]:
            n_bonds = len(atom.bonded_atom_indices)

            # check if the number of bonds is a valid value
            num_bonds_check = atom.check_bond_count(n_bonds)
            # OK
            if num_bonds_check == 0:
                pass
            # missing bond(s)
            elif num_bonds_check < 0:

                # check if any bonds can be formed outside the usual
                # range, that can make the number of bonds and bond
                # angles consistent

                # go through any possible bond that the current atom
                # can be involved in
                for index1, index2 in bond_indices:
                    # only need to check index1, as index1 is always less than index2
                    # TODO: double check that the second equality is not too restrictive
                    if index1 < atom.index:
                        continue
                    elif index1 == atom.index and self.atoms[index2].validated == False:
                        atom.bonded_atom_indices.append(index2)
                        self.atoms[index2].bonded_atom_indices.append(index1)

                        next_subset.extend([index1,index2])
                        break
                    else:
                        eprint(2, "WARNING: Unusual number of bonds to atom %s" % atom.index)
                        break
                continue
            # extra bond(s)
            else:
                # TODO implement solution
                eprint(2, "WARNING: Unusual number of bonds to atom %s" % atom.index)

            #status = atom.check_planarity(n_bonds) #TODO both too many and too few?

            # flag the atom as validated
            # TODO update if extra conditions are added
            atom.validate()

        # reiterate with non-validated atoms if needed:
        if len(next_subset) > 0:
            self.check_bonds(bond_indices, subset = next_subset)

        else:
            # warning if not all atoms are validated
            for atom in self.atoms:
                if atom.validated == False:
                    eprint(2, "WARNING: something's off with the geometry of atom %s" % atom.index)



class Atom(object):
    """
    atom(index, molecule)

    Parameters
    ----------
    index: integer
        index of atom in the coordinate file
    molecule: object
        Molecule object

    Attributes
    ----------
    element_symbol: string
        atomtype, e.g. 'H' or 'Cl'
    molecule_index: integer
        index of atom in the molecule
    mixture_index: integer
        index of atom in the mixture
    coordinates: array-like
        3-sized array of euclidian coordinates
    relative_normed_coordinates: array-like
        Nx3 sized array of normed direction vectors
    distances: array-like
        N-sized array of euclidean distances between the atom and all
        atoms in the molecule.
    bonded_atom_indices: list
        Indices of atoms in the molecule that possibly form a covalent bond
        with the atom. The list is pruned to be consistent if an excess
        number of atoms are within the bond length cut-off.
    validated: bool
        Whether the atom has passed consistency checks on bond lengths and angles
    sybyl: string
        Sybyl atom type
    num_bonds: int
        number of bonded atoms

    """

    def __init__(self, index, molecule):
        self.element_symbol = molecule.element_symbols[index]
        self.molecule_index = index
        self.mixture_index = None
        self.coordinates = molecule.coordinates[index]
        self.distances = molecule.distance_matrix[index]
        self.relative_normed_coordinates = self.get_relative_normed_coordinates(molecule)
        self.bonded_atom_indices = []
        self.validated = False
        self.sybyl = self.element_symbol
        self.num_bond = 0

    def set_sybyl_type(self, s):
        self.sybyl = s

    def get_sybyl_type(

    def get_relative_normed_coordinates(self, molecule):
        # ignore divide by zero warning
        #with warnings.catch_warnings():
        #    warnings.filterwarnings("ignore", r'invalid value encountered in divide')
        with np.errstate(divide="ignore", invalid="ignore"):
            return (molecule.coordinates - self.coordinates[None,:]) / self.distances[:,None]

    def set_mixture_index(self, index):
        self.mixture_index = index

    # TODO: Less ugly returns
    def check_bond_count(self, n_bonds):
        if self.element_symbol not in constants.number_bonds:
            eprint(2, "WARNING: element %s not completely implemented, but ordering should still work")
            return True
        possible_num_bonds = constants.number_bonds[self.element_symbol]
        # return 0 if passed check. 1 if theres too many bonds. -1 if theres too few.
        if n_bonds in possible_num_bonds: return 0
        if n_bonds < possible_num_bonds.min(): return possible_num_bonds.min() - n_bonds
        if n_bonds > possible_num_bonds.max(): return n_bonds - possible_num_bonds.max()


    def set_molecule_index(self, index):
        self.molecule_index = index

    def add_bond(self, atom):
        self.bonded_atom_indices.append(atom)
        self.num_bond += 1

    def validate(self):
        self.validated = True

    def get_bond_distances(self):
        return [self.distances[i] for i in self.bonded_atom_indices]

    def get_average_bond_angle(self):
        bond_angles = self.get_bond_angles()
        return np.mean(bond_angles)

    def get_bond_angles(self):
        # http://stackoverflow.com/a/13849249
        def angle_v(v1,v2):
            return np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))

        angles = []
        pair_iterator = itertools.combinations(self.bonded_atom_indices,2)
        for i,j in pair_iterator:
            v1 = self.relative_normed_coordinates[i]
            v2 = self.relative_normed_coordinates[j]
            angle = angle_v(v1,v2) * (180/np.pi)
            # TODO remove at release
            assert(angle >= 0)
            assert(angle <= 180)

            angles.append(angle)

        return angles




    # p1 central point
    def angle_p(p):
        p1,p2,p3 = p
        d12 = sum((p1-p2)**2)
        d13 = sum((p1-p3)**2)
        d23 = sum((p2-p3)**2)
        max_idx = np.argmax([d12,d13,d23])
        if max_idx == 0:
            v1 = p1-p3
            v2 = p2-p3
        elif max_idx == 1:
            v1 = p1 - p2
            v2 = p3 - p2
        else:
            v1 = p2 - p1
            v2 = p3 - p1
        return angle_v(v1,v2)*(180/np.pi)

    #def get_distance(self, x, axis = -1):
    #    """
    #    Returns distance between the atom and a set of coordinates
    #    or an atom index, depending on type

    #    Parameters:
    #    -----------
    #    x: array or index
    #        coordinates or coordinate index

    #    """
    #    try:
    #        return self.distances
    #    except TypeError

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
