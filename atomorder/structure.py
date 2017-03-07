import numpy as np
import itertools
import warnings
from . import settings, constants, utils
from .utils import oprint, eprint, read_coordinates

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
    atoms: array_like
        atoms in the molecule
    bonds: array_like
        bonds in the molecule

    """

    def __init__(self, filename):
        self.filename = filename
        self.element_symbols, self.coordinates = read_coordinates(self.filename)
        self.size = self.element_symbols.size
        self.distance_matrix = self.get_distance_matrix()
        self.monovalent = self.get_monovalent_indices()
        self.atoms = self.make_atoms()
        self.bonds = self.find_bonds()
        quit()
        self.get_sybyl_types()
        self.get_bond_order()
        # TODO add matched sybyl type check
        quit()
        # TODO make bond2
        # TODO split
        #self.make_graph() - move?

    def get_bond_order(self):
        """
        Find which bonds are single/double/triple/aromatic etc. self consistently.
        This is done by assigning electrons not forming sigma bonds to pi bonds
        or lone pairs (charges).

        """

        # Generate list of all non-sigma electrons and
        # list of 



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
    #def split(self):
    #    """
    #    split()

    #    Splits the molecule into several molecule objects if two or more
    #    disconnected regions are present (e.g. ligands or several reactants
    #    in the same file)

    #    Returns
    #    -------
    #    molecules: tuple
    #        one or more molecule objects
    #    """

    def make_atoms(self):
        """
        Make array of atom objects

        Returns
        -------
        atoms: array-like
            Atom objects

        """

        oprint(3, "Creating atoms")

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

        oprint(3, "Locating covalent bonds")

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
            self.add_bond(index1, index2)
        quit()

        # Check that all atoms have a realistic set of bonds
        # If not, try to add bonds between atoms from less_credible_bond_indices
        # until everything looks fine
        self.check_bonds(less_credible_bond_indices)

    def add_bond(self, index1, index2):
        bond = Bond(index1, index2)
        self.bonds.append(bond)
        atom1 = self.atoms[index1]
        atom2 = self.atoms[index2]
        atom1.add_bond(atom2)
        atom2.add_bond(atom1)

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
                    if index1 < atom.molecule_index:
                        continue
                    elif index1 == atom.molecule_index and self.atoms[index2].validated == False:
                        atom.bonded_atom_indices.append(index2)
                        self.atoms[index2].bonded_atom_indices.append(index1)

                        next_subset.extend([index1,index2])
                        break
                    else:
                        eprint(2, "WARNING: Unusual number of bonds to atom %s" % atom.molecule_index)
                        break
                continue
            # extra bond(s)
            elif num_bonds_check > 0:
                # TODO implement solution
                eprint(2, "WARNING: Unusual number of bonds to atom %s" % atom.molecule_index)
                continue

            # not implemented element / type
            else:
                continue


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
                    eprint(2, "WARNING: something's off with the geometry of atom %s" % atom.molecule_index)



class Atom(object):
    """
    Atom(index, molecule)

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
    bonded_atoms: list
        Atoms with indices given by bonded_atom_indices
    validated: bool
        Whether the atom has passed consistency checks on bond lengths and angles
    sybyl: string
        Sybyl atom type
    num_bonds: int
        number of bonded atoms
    charge: float
        charge in fractions of plus/minus 1

    """

    def __init__(self, index, molecule):
        self.element_symbol = molecule.element_symbols[index]
        self.molecule_index = index
        self.mixture_index = None
        self.coordinates = molecule.coordinates[index]
        self.distances = molecule.distance_matrix[index]
        self.relative_normed_coordinates = self.get_relative_normed_coordinates(molecule)
        self.bonded_atom_indices = []
        self.bonded_atoms = []
        self.validated = False
        self.sybyl = self.element_symbol
        self.num_bond = 0
        self.charge = 0

    def set_bonded_atoms(self, molecule):
        return [molecule.atoms[i] for i in self.bonded_atom_indices]

    def set_sybyl_type(self, s):
        self.sybyl = s

    def get_sybyl_type(self):
        """
        CSD non-matched (3d) deterministic sybyl atom type matching
        from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
        Comments refer to the definitions given from the above site

        """

        # 2.2 If element_symbol is D then atom_type is H
        if self.element_symbol == "D":
            sybyl = "H"
        # 2.3 If element_symbol is P then atom_type is P.3
        elif self.element_symbol == "P":
            sybyl = "P.3"
        # 2.5 If element_symbol is C then
        elif self.element_symbol == "C":
            sybyl = self.get_sybyl_type_C()
        # 2.6 If element_symbol is O then
        elif self.element_symbol == "O":
            sybyl = self.get_sybyl_type_O()
        # 2.7 If element_symbol is N then
        elif self.element_symbol == "N":
            sybyl = self.get_sybyl_type_N()
        # 2.8 If element_symbol is S then
        elif self.element_symbol == "S":
            sybyl = self.get_sybyl_type_S()
        # 2.10 If element_symbol is none of the above then atom_type is element_symbol
        else:
            sybyl = self.element_symbol

        self.set_sybyl_type(sybyl)

    def get_sybyl_type_C(self):
        """
        CSD non-matched (3d) deterministic sybyl atom type matching
        from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

        """

        # 2.5.1 If num_bond .ge. 4 then atom_type is C.3
        if self.num_bond >= 4: return "C.3"
        # 2.5.2 If num_bond .eq. 1 then calculate bond_distance
        if self.num_bond == 1:
            bond_distance = self.get_bond_distances()[0]

            # 2.5.2.1 If bond_distance .gt. 1.41A then atom_type is C.3
            if bond_distance > 1.41: return "C.3"
            # 2.5.2.2 If bond_distance .le. 1.22A then atom_type is C.1
            if bond_distance <= 1.22: return "C.1"
            # 2.5.2.3 If bond_distance is none of the above then atom_type is C.2
            return "C.2"
        # 2.5.3 If element_symbol is C and none of the above then calculate average_angle about C
        average_angle = self.get_average_bond_angle()
        # 2.5.3.1 If average_angle .le. 115 deg then atom_type is C.3
        if average_angle <= 115: return "C.3"
        # 2.5.3.2 If average_angle .gt. 160 deg then atom_type is C.1
        if average_angle > 160: return "C.1"
        # 2.5.3.3 If average_angle is none of the above then atom_type is C.2
        return "C.2"

    def get_sybyl_type_O(self):
        """
        CSD non-matched (3d) deterministic sybyl atom type matching
        from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

        """

        # 2.6.1 If num_nonmet .eq. 1 then
        if self.num_bond == 1:
            # 2.6.1.1 If bond is to carbon .AND. carbon forms a total of 3
            #         bonds, 2 of which are to an oxygen forming only 1
            #         non-metal bond then atom_type is O.co2
            neighbour = self.bonded_atoms[0]
            if neighbour.element_symbol == "C" and neighbour.num_bond == 3:
                for neighbours_neighbour in neighbour.bonded_atoms:
                    if neighbours_neighbour.element_symbol == "O" and neighbours_neighbour.molecule_index != self.molecule_index \
                            and this_neighbor.num_bond == 1:
                        return "O.co2"
            # 2.6.1.2 If bond is to phosphorus .AND. phosphorous forms at
            #         least 2 bonds to an oxygen forming only 1 non-metal
            #         bond then atom_type is O.co2
            if neighbour.element_symbol == "P":
                for neighbours_neighbour in neighbour.bonded_atoms:
                    if neighbours_neighbour.element_symbol == "O" and neighbours_neighbour.molecule_index != self.molecule_index \
                            and neighbours_neighbour.num_bond == 1:
                        return "O.co2"
        # 2.6.3 If num_bond .ge. 2 then atom_type is O.3
        if self.num_bond >= 2: return "O.3"

        # 2.6.4 If element_symbol is O and none of the above then atom_type is O.2
        return "O.2"

    def get_sybyl_type_N(self):
        """
        CSD non-matched (3d) deterministic sybyl atom type matching
        from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

        """

        # 2.7.1 If num_nonmet .eq. 4 then atom_type is N.4
        if self.num_bond == 4: return "N.4"

        # 2.7.2 If num_nonmet .eq. 1 then calculate bond_distance
        if self.num_bond == 1:
            bond_distance = self.get_bond_distances()[0]

            # 2.7.2.1 If bond_distance .le. 1.2A then atom_type is N.1
            if bond_distance <= 1.2: return "N.1"

            # 2.7.2.2 If bond_distance .gt. 1.2A then atom_type is N.3
            return "N.3"

        # 2.7.3 If num_nonmet .eq. 3 .AND. one bond is to C--O or C--S then atom_type is N.am
        if self.num_bond == 3:
            for neighbour in self.bonded_atoms:
                if neighbour.element_symbol == "C":
                    for neighbours_neighbour in neighbour.bonded_atoms:
                        # At this stage it's not known if a bond is double or single
                        # so just check that O and S only bonds to the parent carbon
                        if neighbours_neighbor.element_symbol in ("O","S") and \
                                this_neighbor.num_bond == 1:
                                    return "N.am"
            # 2.7.4 If num_nonmet .eq. 3 otherwise then calculate sum_of_angles around N
            angle_sum = self.get_angle_sum()

            # 2.7.4.1 If sum_of_angles .ge. 350 deg then atom_type is N.pl3
            if angle_sum >= 350: return "N.pl3"
            # 2.7.4.2 If sum_of_angles .lt. 350 deg then atom_type is N.3
            return "N.3"

        # 2.7.5 If element_symbol is N and none of the above then calculate average_angle about N
        average_angle = self.get_average_angle()

        # 2.7.5.1 If average_angle .gt. 160 deg then atom_type is N.1
        if average_angle > 160: return "N.1"
        # 2.7.5.2 If average_angle .le. 160 deg then atom_type is N.2
        return "N.2"

    def get_sybyl_type_S(self):
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
        if self.num_bond == 3:
            for neighbour in self.bonded_atoms:
                if neighbour.element_symbol == "O" and neighbour.num_bond == 1: return "S.o"
        # 2.8.2  If num_nonmet .eq. 4 .AND. 2 bonds are to an oxygen with only
        #        one non-metal bond then atom_type is S.o2
        if self.num_bond == 4:
            count = 0
            for neighbour in self.bonded_atoms:
                if neighbour.element_symbol == "O" and neighbour.num_bond == 1: count += 1
                if count == 2: return "S.o2"
        # 2.8.3 If num_bond .ge. 2 then atom_type is S.3
        if self.num_bond >= 2: return "S.3"
        # 2.8.4 If element_symbol is S and none of the above then atom_type is S.2
        return "S.2"

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

    def add_bond(self, bonded_atom):
        bonded_atom_index = bonded_atom.molecule_index
        self.bonded_atom_indices.append(bonded_atom_index)
        self.bonded_atoms.append(bonded_atom)
        self.num_bond += 1

    def validate(self):
        self.validated = True

    def get_bond_distances(self):
        return [self.distances[i] for i in self.bonded_atom_indices]

    def get_bond_angle_sum(self):
        bond_angles = self.get_bond_angles()
        return np.sum(bond_angles)

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

class Bond(object):
    """
    Bond(atom1, atom2)

    Parameters
    ----------
    atom1: object
        Atom object
    atom2: object
        Atom object

    Attributes
    ----------
    is_type: string
        What type of bond this is
    can_be_single: bool
        If the bond can be of type single
    can_be_double: bool
        If the bond can be of type double
    can_be_triple: bool
        If the bond can be of type triple
    can_be_conjugated: bool
        If the bond can be of type conjugated
    atom1: object
        Atom object
    atom2: object
        Atom object
    molecule_index: int
        index of bond in molecule

    """
    def __init__(self, atom1, atom2):
        self.is_type = "Unknown"
        self.can_be_single = True
        self.can_be_double = True
        self.can_be_triple = True
        self.can_be_conjugated = True
        self.atom1 = atom1
        self.atom2 = atom2
        self.molecule_index = None
