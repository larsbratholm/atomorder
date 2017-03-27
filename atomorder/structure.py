import numpy as np
import itertools
import warnings
from . import settings, constants, utils
from .utils import oprint, eprint, vector_angle
from .parsers import read_coordinates

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
    # TODO charge/atom conservation

    def __init__(self):
        oprint(3, "Creating reactants")
        self.reactants = Mixture(settings.reactant_filenames)
        self.num_reactant_atoms = self.reactants.get_molecule_sizes()
        oprint(3, "Creating products")
        self.products = Mixture(settings.product_filenames)
        self.num_product_atoms = self.products.get_molecule_sizes()


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
    num_mol: int
        number of molecules
    molecules: list
        Molecule objects.
    atoms: array
        array of atoms in all the molecules
    num_atoms: int
        number of atoms in all the molecules
    coordinates: array
        coordinates of atoms in all the molecules
    element_symbols: array
        element symbols of atoms in all the molecules
    sybyl_atom_types:
        sybyl types of atoms in all the molecules
    bond_matrix: array
        connectivity between atoms

    """
    # TODO total charge

    def __init__(self, filenames):
        self.filenames = filenames
        self.num_mol = len(filenames)
        self.molecules = self.make_molecules()
        self.atoms = self.get_atoms()
        self.num_atoms = self.atoms.size
        self.set_atom_mixture_indices()
        self.coordinates = self.get_coordinates()
        self.element_symbols = self.get_element_symbols()
        if settings.create_atoms:
            self.sybyl_atom_types = self.get_sybyl_atom_types()
            self.bond_matrix = self.get_bond_matrix()

    def get_molecule_sizes(self):
        return [molecule.size for molecule in self.molecules]

    def get_coordinates(self):
        return np.concatenate([molecule.coordinates for molecule in self.molecules])

    def get_sybyl_atom_types(self):
        return np.concatenate([molecule.sybyl_atom_types for molecule in self.molecules])

    def get_element_symbols(self):
        return np.concatenate([molecule.element_symbols for molecule in self.molecules])

    # A molecule object will contain everything in a single coordinate file
    def make_molecules(self):
        """
        make_molecules()

        Creates list of all given molecules.

        """

        molecules = np.empty(self.num_mol, dtype=object)
        for index, filename in enumerate(self.filenames):
            oprint(3, "Creating molecule from %s" % filename)
            molecules[index] = Molecule(index, filename)
        return molecules

    def get_atoms(self):
        """
        Construct array of all atoms in the mixture

        """
        return  np.concatenate([mol.atoms for mol in self.molecules])


    def set_atom_mixture_indices(self):
        """
        In case of multiple molecules, set unique identifier
        for each atom

        """
        if self.num_mol > 1:
            for i, atom in enumerate(self.atoms):
                atom.set_mixture_index(i)

    def get_bond_matrix(self):
        """
        Create binary matrix with the bond network

        """

        connectivity = np.zeros((self.atoms.size, self.atoms.size), dtype=bool)

        for i, atom_i in enumerate(self.atoms):
            for atom_j in atom_i.bonded_atoms:
                connectivity[atom_i.molecule_index, atom_j.molecule_index] = 1

        return connectivity


class Molecule(object):
    """
    molecule(index, filename)

    All atoms from a single coordinate file.

    Parameters
    ----------
    index: int
        index of the molecule in the mixture
    filename: string
        Coordinate filename

    Atributes
    ---------
    mixture_index: int
        index of the molecule in the mixture
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
    num_atom: int
        number of atoms
    bonds: array_like
        bonds in the molecule
    num_bonds: int
        number of bonds
    num_atoms: int
        number of atoms
    sybyl_atom_types: array
        sybyl type of the atoms

    """
    #TODO total charge

    def __init__(self, index, filename):
        self.mixture_index = index
        self.filename = filename
        self.element_symbols, self.coordinates = read_coordinates(self.filename)
        self.monovalent = self.get_monovalent_indices()
        self.size = self.element_symbols.size
        self.distance_matrix = None
        self.atoms = []
        self.bonds = []
        self.num_bond = 0
        self.num_atoms = self.element_symbols.size

        # stop here if there's no need for bond information
        if settings.create_atoms == False:
            return

        self.create_atoms_and_bonds()
        self.sybyl_atom_types = self.get_sybyl_atom_types()


    def get_sybyl_atom_types(self):
        return np.asarray([atom.sybyl for atom in self.atoms])

    def create_atoms_and_bonds(self):
        self.distance_matrix = self.get_distance_matrix()
        self.atoms = self.make_atoms()
        self.get_bonds()
        self.num_bonds = len(self.bonds)
        self.get_sybyl_types()
        return
        self.get_bond_order()
        # TODO add matched sybyl type check
        quit()
        # TODO make bond2
        # TODO split
        #self.make_graph() - move?

    def get_bond_order(self):
        #"""
        #Find which bonds are single/double/triple/aromatic etc. self consistently.
        #This is done by a breadth-first search to determine possible conjugation
        #in the molecule, followed by construction of resonance structures

        #"""

        # keep track of which atoms have all bonds types determined
        completed_atoms = np.zeros(self.size, dtype=bool)

        # no need to check monovalent atoms to the queue
        completed_atoms[self.monovalent] = True


        # keep track of atoms which have had bonds to them assigned
        atom_queue = np.zeros(self.size, dtype=bool)

        # First assign all bonds that can only be single
        for i, bond in enumerate(self.bonds):
            was_single = bond.check_if_single()
            # add involved atoms to queue
            if was_single:
                index1 = bond.atom1.molecule_index
                index2 = bond.atom2.molecule_index

                # only add atom to queue if its bonds are not completely determined
                if completed_atoms[index1] == False:
                    atom_queue[bond.atom1.molecule_index] = True
                if completed_atoms[index2] == False:
                    atom_queue[bond.atom2.molecule_index] = True

        print [atom.sybyl for atom in self.atoms]

        # Flag to keep track if there's any changes in the iteration
        converged = False
        # max 1000 iterations to avoid getting stuck in while loop
        for it in range(1000):
            print "iteration", it
            if converged == True:
                print "converged"
                break
            converged = True

            # go through queue
            for atom_index in np.where(atom_queue == True)[0]:
                print "index", atom_index
                atom = self.atoms[atom_index]

                # check if any bonds involving the current atom can be determined
                neighbour_indices = atom.determine_bond_type()
                print neighbour_indices

                # neighbour_index is None if no new bonds had their types determined.
                if neighbour_indices == None:
                    continue
                # Something changed, so we're not converged
                converged = False
                # remove atom from queue since bonds were determined
                atom_queue[atom_index] = False
                # register that we're forever done with this atom
                completed_atoms[atom_index] = True
                # if bonds involving the neighbour_index is still not fully determined, add it to the queue
                for neighbour_index in neighbour_indices:
                    print neighbour_index
                    if completed_atoms[neighbour_index] == False:
                        atom_queue[neighbour_index] = True
        print completed_atoms


        ## find unassigned pi electrons and lone pairs
        #electrons = []
        #lone_pairs = []
        #for atom in self.atoms:
        #    n_electrons = atom.sybyl_properties[0]
        #    n_lone_pairs = atom.sybyl_properties[1]

        #    for _ in xrange(n_electrons):
        #        electrons.append(atom)
        #    for _ in xrange(n_lone_pairs):
        #        lone_pairs.append(atom)

        ## find possible pi bonds
        #bonds = []
        #for bond in self.bonds:
        #    if bond.atom1 in electrons and bond.atom2 in electrons:
        #        bonds.append(bond)

        ## construct assignment matrix
        #m = np.zeros( (len(electrons), len(lone_pairs) + len(bonds)) )

        ## initialize by naive distribution



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

    def get_bonds(self):
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
        limits = np.asarray([constants.get_bond_length_limits(self.element_symbols[i],self.element_symbols[j]) for i,j in zip(*pair_indices)]
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

        # Check that all atoms have a realistic set of bonds
        # If not, try to add bonds between atoms from less_credible_bond_indices
        # until everything looks fine
        self.check_bonds(less_credible_bond_indices)

    def add_bond(self, index1, index2):
        bond = Bond(self, index1, index2)
        self.bonds.append(bond)

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
            # check if the number of bonds is a valid value
            num_bonds_check = atom.check_bond_count()

            # OK
            if num_bonds_check == 0:
                pass
            # missing bond(s)
            elif num_bonds_check > 0:

                # check if any bonds can be formed outside the usual
                # range, that can make the number of bonds and bond
                # angles consistent

                # go through any possible bond that the current atom
                # can be involved in
                for index1, index2 in bond_indices:
                    # only need to check index2, as index1 is always less than index2
                    # this should make sure that excess bonds is not added to an atom
                    if index2 > atom.molecule_index:
                        continue
                    # TODO add consistency check, because this might fail
                    if index2 == atom.molecule_index:# and self.atoms[index1].validated == False:
                        self.add_bond(index1, index2)
                        next_subset.extend([index1,index2])
                        break
                    if index1 >= atom.molecule_index:
                        eprint(2, "WARNING: Unusual number of bonds to atom %s" % atom.molecule_index)
                        break
                continue
            # extra bond(s)
            elif num_bonds_check < 0:
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
            next_subset = np.unique(next_subset)
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
        Index of atom in the coordinate file
    molecule: object
        Molecule object

    Attributes
    ----------
    molecule: object
        Parent Molecule object
    element_symbol: string
        Atomtype, e.g. 'H' or 'Cl'
    molecule_index: integer
        Index of atom in the molecule
    mixture_index: integer
        Index of atom in the mixture
    coordinates: array-like
        3-sized array of euclidian coordinates
    relative_normed_coordinates: array-like
        Nx3 sized array of normed direction vectors
    distances: array-like
        N-sized array of euclidean distances between the atom and all
        atoms in the molecule.
    bonded_atoms: list
        atoms in the molecule that possibly form a covalent bond
        with the atom. The list is pruned to be consistent if an excess
        number of atoms are within the bond length cut-off.
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
        self.mixture_index = index
        self.coordinates = molecule.coordinates[index]
        self.distances = molecule.distance_matrix[index]
        self.relative_normed_coordinates = self.get_relative_normed_coordinates(molecule)
        self.bonded_atoms = []
        self.bonds = []
        self.validated = False
        self.sybyl = self.element_symbol
        self.num_bonds = 0
        self.charge = None
        self.sybyl_properties = None


    def set_sybyl_type(self, s):
        """
        Sets sybyl type and various properties

        """
        self.sybyl = s
        # TODO make this not fail if weird structures are tried
        self.sybyl_properties = constants.sybyl_bonds[self.sybyl][self.num_bonds]
        self.charge = self.sybyl_properties[2]

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

        # 2.5.1 If num_bonds .ge. 4 then atom_type is C.3
        if self.num_bonds >= 4: return "C.3"
        # 2.5.2 If num_bonds .eq. 1 then calculate bond_distance
        if self.num_bonds == 1:
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
        if self.num_bonds == 1:
            # 2.6.1.1 If bond is to carbon .AND. carbon forms a total of 3
            #         bonds, 2 of which are to an oxygen forming only 1
            #         non-metal bond then atom_type is O.co2
            neighbour = self.bonded_atoms[0]
            if neighbour.element_symbol == "C" and neighbour.num_bonds == 3:
                for neighbours_neighbour in neighbour.bonded_atoms:
                    if neighbours_neighbour.element_symbol == "O" and neighbours_neighbour.molecule_index != self.molecule_index \
                            and neighbours_neighbour.num_bonds == 1:
                        return "O.co2"
            # 2.6.1.2 If bond is to phosphorus .AND. phosphorous forms at
            #         least 2 bonds to an oxygen forming only 1 non-metal
            #         bond then atom_type is O.co2
            if neighbour.element_symbol == "P":
                for neighbours_neighbour in neighbour.bonded_atoms:
                    if neighbours_neighbour.element_symbol == "O" and neighbours_neighbour.molecule_index != self.molecule_index \
                            and neighbours_neighbour.num_bonds == 1:
                        return "O.co2"
        # 2.6.3 If num_bonds .ge. 2 then atom_type is O.3
        if self.num_bonds >= 2: return "O.3"

        # 2.6.4 If element_symbol is O and none of the above then atom_type is O.2
        return "O.2"

    def get_sybyl_type_N(self):
        """
        CSD non-matched (3d) deterministic sybyl atom type matching
        from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

        """

        # 2.7.1 If num_nonmet .eq. 4 then atom_type is N.4
        if self.num_bonds == 4: return "N.4"

        # 2.7.2 If num_nonmet .eq. 1 then calculate bond_distance
        if self.num_bonds == 1:
            bond_distance = self.get_bond_distances()[0]

            # 2.7.2.1 If bond_distance .le. 1.2A then atom_type is N.1
            if bond_distance <= 1.2: return "N.1"

            # 2.7.2.2 If bond_distance .gt. 1.2A then atom_type is N.3
            return "N.3"

        # 2.7.3 If num_nonmet .eq. 3 .AND. one bond is to C--O or C--S then atom_type is N.am
        if self.num_bonds == 3:
            for neighbour in self.bonded_atoms:
                if neighbour.element_symbol == "C":
                    for neighbours_neighbour in neighbour.bonded_atoms:
                        # At this stage it's not known if a bond is double or single
                        # so just check that O and S only bonds to the parent carbon
                        if neighbours_neighbour.element_symbol in ("O","S") and \
                                neighbours_neighbour.num_bonds == 1:
                            return "N.am"
            # 2.7.4 If num_nonmet .eq. 3 otherwise then calculate sum_of_angles around N
            angle_sum = self.get_angle_sum()

            # 2.7.4.1 If sum_of_angles .ge. 350 deg then atom_type is N.pl3
            if angle_sum >= 350: return "N.pl3"
            # 2.7.4.2 If sum_of_angles .lt. 350 deg then atom_type is N.3
            return "N.3"

        # 2.7.5 If element_symbol is N and none of the above then calculate average_angle about N
        average_angle = self.get_average_bond_angle()

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
        if self.num_bonds == 3:
            for neighbour in self.bonded_atoms:
                if neighbour.element_symbol == "O" and neighbour.num_bonds == 1: return "S.o"
        # 2.8.2  If num_nonmet .eq. 4 .AND. 2 bonds are to an oxygen with only
        #        one non-metal bond then atom_type is S.o2
        if self.num_bonds == 4:
            count = 0
            for neighbour in self.bonded_atoms:
                if neighbour.element_symbol == "O" and neighbour.num_bonds == 1: count += 1
                if count == 2: return "S.o2"
        # 2.8.3 If num_bonds .ge. 2 then atom_type is S.3
        if self.num_bonds >= 2: return "S.3"
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
    def check_bond_count(self):
        # return 0 if passed check. 1 if theres too many bonds. -1 if theres too few.
        if self.element_symbol not in constants.number_bonds:
            eprint(2, "WARNING: element %s not completely implemented, but ordering should still work")
            return True
        possible_num_bonds = constants.number_bonds[self.element_symbol]
        if self.num_bonds in possible_num_bonds: return 0
        if self.num_bonds < possible_num_bonds.min(): return possible_num_bonds.min() - self.num_bonds
        if self.num_bonds > possible_num_bonds.max(): return self.num_bonds - possible_num_bonds.max()

    def set_molecule_index(self, index):
        self.molecule_index = index

    def add_bond(self, bond):
        if bond.atom1 == self:
            bonded_atom = bond.atom2
        else:
            bonded_atom = bond.atom1

        self.bonds.append(bond)
        self.bonded_atoms.append(bonded_atom)
        self.num_bonds += 1

    def validate(self):
        self.validated = True

    def get_bond_distances(self):
        return [atom.distances[self.molecule_index] for atom in self.bonded_atoms]


    def get_angle_sum(self):
        bond_angles = self.get_bond_angles()
        return np.sum(bond_angles)

    def get_average_bond_angle(self):
        bond_angles = self.get_bond_angles()
        return np.mean(bond_angles)

    def get_bond_angles(self):
        angles = []
        pair_iterator = itertools.combinations(self.bonded_atoms,2)
        for atom_i,atom_j in pair_iterator:
            v1 = self.relative_normed_coordinates[atom_i.molecule_index]
            v2 = self.relative_normed_coordinates[atom_j.molecule_index]
            angle = vector_angle(v1,v2) * (180/np.pi)
            # TODO remove at release
            assert(angle >= 0)
            assert(angle <= 180)

            angles.append(angle)

        return angles

    def determine_bond_type(self):
        """
        Tries to determine type of any bonds involving self.

        If bonds were determined, return indices
        of the other atoms. Else return None.

        """
        # TODO shorten this when confirmed working

        n_unassigned_electrons = self.sybyl_properties[0]
        determined_bonds = np.zeros(self.num_bonds, dtype=bool)
        lone_pairs = self.sybyl_properties[1]

        for i, bond in enumerate(self.bonds):
            if bond.is_type == None:
                continue
            determined_bonds[i] = True
            if bond.is_type == "double":
                n_unassigned_electrons -= 1
            elif bond.is_type == "triple":
                n_unassigned_electrons -= 2

        if n_unassigned_electrons == 2 and (self.num_bonds - determined_bonds.sum()) == 1 and lone_pairs == 0:
            bond_index = np.where(determined_bonds == False)[0][0]
            self.bonds[bond_index].set_type("triple")
            return [self.bonded_atoms[bond_index].molecule_index]

        if n_unassigned_electrons == 0 and (self.num_bonds - determined_bonds.sum()) == 0 and lone_pairs == 0:
            # all bonds were already determined
            # just return index of an arbitrary bond
            return [self.bonded_atoms[0].molecule_index]

        if n_unassigned_electrons == 1 and (self.num_bonds - determined_bonds.sum()) == 0 and lone_pairs == 1:
            # there's a left over electron on the atom
            self.charge -= 1
            # all bonds were already determined
            # just return index of an arbitrary bond
            return [self.bonded_atoms[0].molecule_index]

        if n_unassigned_electrons == 0 and (self.num_bonds - determined_bonds.sum()) == 0 and lone_pairs == 1:
            # all bonds were already determined
            # just return index of an arbitrary bond
            return [self.bonded_atoms[0].molecule_index]

        if n_unassigned_electrons == 1 and (self.num_bonds - determined_bonds.sum()) == 2 and lone_pairs == 0:
            return None

        if n_unassigned_electrons == 1 and (self.num_bonds - determined_bonds.sum()) == 1 and lone_pairs == 0:
            bond_index = np.where(determined_bonds == False)[0][0]
            self.bonds[bond_index].set_type("double")
            return [self.bonded_atoms[bond_index].molecule_index]

        if n_unassigned_electrons == 0 and (self.num_bonds - determined_bonds.sum()) == 1 and lone_pairs == 0:
            bond_index = np.where(determined_bonds == False)[0][0]
            self.bonds[bond_index].set_type("single")
            return [self.bonded_atoms[bond_index].molecule_index]

        if n_unassigned_electrons == 0 and (self.num_bonds - determined_bonds.sum()) == 2 and lone_pairs == 0:
            bond_index1, bond_index2 = np.where(determined_bonds == False)[0]
            self.bonds[bond_index1].set_type("single")
            self.bonds[bond_index2].set_type("single")
            return [self.bonded_atoms[bond_index1].molecule_index,self.bonded_atoms[bond_index2].molecule_index]

        if n_unassigned_electrons == 1 and (self.num_bonds - determined_bonds.sum()) == 1 and lone_pairs == 1:
            # e.g. N.am
            return None

        if n_unassigned_electrons == 1 and (self.num_bonds - determined_bonds.sum()) == 2 and lone_pairs == 1:
            return None

        if n_unassigned_electrons == 0 and (self.num_bonds - determined_bonds.sum()) == 1 and lone_pairs == 1:
            bond_index = np.where(determined_bonds == False)[0][0]
            self.bonds[bond_index].set_type("single")
            return [self.bonded_atoms[bond_index].molecule_index]

        print n_unassigned_electrons, determined_bonds, self.sybyl_properties[1]
        print self.sybyl, self.num_bonds, self.sybyl_properties
        quit("Add to determine_bond_type")

class Bond(object):
    """
    Bond(molecule, index1, index2)

    Parameters
    ----------
    molecule: object
        parent molecule
    index1: int
        molecule_index of first atom
    index2: int
        molecule_index of second atom

    Attributes
    ----------
    is_type: string / Nonetype
        What type of bond this is
    atom1: object
        Atom object
    atom2: object
        Atom object
    molecule_index: int
        index of bond in molecule



    """
    def __init__(self, molecule, index1, index2):
        self.is_type = None
        self.molecule = molecule
        self.atom1 = molecule.atoms[index1]
        self.atom2 = molecule.atoms[index2]
        self.molecule_index = len(molecule.bonds)
        self.add_bond_to_atoms()

    def add_bond_to_atoms(self):
        self.atom1.add_bond(self)
        self.atom2.add_bond(self)

    def check_if_single(self):
        """
        Checks sybyl types of both atoms to infer if the bond
        can only be single.

        Returns:
        --------
        was_single: bool
            If the bond was single

        """

        # if an atom does not have any free pi electrons to contribute to the bond
        # the bond must be single
        n_pi1 = self.atom1.sybyl_properties[0]
        n_pi2 = self.atom2.sybyl_properties[0]

        min_pi = min(n_pi1, n_pi2)
        if min_pi == 0:
            self.set_type("single")
            was_single = True
        else:
            was_single = False

        return was_single

    def set_type(self, s):
        self.is_type = s.lower()



    # TODO set_pi_electrons
