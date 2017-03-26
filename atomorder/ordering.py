"""
atomorder/ordering.py

The actual ordering

"""

import numpy as np
from . import settings, objectives
from .utils import oprint, eprint

# TODO when the atoms much match, and there's the same number of atoms
# the match matrix must be symmetric, and we only have to calculate the upper half
class Ordering(object):
    """
    Ordering(reaction)

    The actual assignment / ordering.

    Parameters:
    -----------
    reaction: object
        Reaction object

    Attributes:
    -----------
    reaction: object
        Reaction object

    #TODO

    """

    def __init__(self, reaction):
        # TODO allow for matching between different elements
        oprint(3, "Initializing ordering")
        self.reaction = reaction
        self.num_reactants = self.reaction.reactants.num_mol
        self.num_products = self.reaction.products.num_mol
        self.reactants_elements = self.reaction.reactants.element_symbols
        self.products_elements = self.reaction.products.element_symbols
        self.reactants_coordinates = self.reaction.reactants.coordinates
        self.products_coordinates = self.reaction.products.coordinates
        if settings.create_atoms:
            self.reactants_sybyl_types = self.reaction.reactants.sybyl_atom_types
            self.products_sybyl_types = self.reaction.products.sybyl_atom_types
        self.check_size_consistency()
        self.element_type_subset_indices = self.get_element_type_subset_indices()
        self.reactant_subset_indices, self.product_subset_indices = self.get_molecule_subset_indices()
        self.match_matrix = self.initialize_match_matrix()
        self.old_match_matrix = self.match_matrix.copy()
        self.score = self.initialize_objectives()
        self.inverse_temperature = 0
        self.Q = np.zeros(self.match_matrix.shape,dtype=float)
        self.initial_inverse_temperature, self.final_inverse_temperature, self.inverse_temperature_increment, \
                self.max_annealing_iterations, self.max_relaxation_iterations, self.max_softassign_iterations, \
                self.annealing_convergence_threshold, self.relaxation_convergence_threshold, self.softassign_convergence_threshold = self.set_parameters()
        self.it = 0

        # do the annealing
        self.converged = self.annealing()
        self.finalize()

    def initialize_objectives(self):
        """
        Create the scoring function used in the ordering

        Returns:
        --------
        objectives: function
            Function that returns the energy of the current match matrix

        """

        self.obj = []

        if settings.rotation_objective:
            self.obj.append(objectives.Rotation(self))

        if settings.atomic_objective:
            self.obj.append(objectives.Atomic(self))

        if settings.bond_objective:
            self.obj.append(objectives.Bond(self))


        # return function that returns the sum of all objective scoring functions
        return (lambda x, obj=self.obj: sum((fun.score(x) for fun in obj)))


    def get_element_type_subset_indices(self):
        """
        It is currently required that the element of two matching atoms is the same.
        This constructs indices to e.g. the carbon-carbon submatrix.

        """
        # TODO: this is redundant if the elements does not have to match
        unique_elements = np.unique(self.reactants_elements)
        subset_indices = np.empty(unique_elements.size, dtype=object)
        for i, element in enumerate(unique_elements):
            rows = np.where(self.reactants_elements == element)[0]
            cols = np.where(self.products_elements == element)[0]
            subset_indices[i] = np.ix_(rows,cols)

        return subset_indices

    def get_molecule_subset_indices(self):
        """
        Get indices for each reactant and each product molecule

        """
        reactant_subset_indices = np.empty(self.num_reactants, dtype=object)
        product_subset_indices = np.empty(self.num_products, dtype=object)
        n_count, m_count = 0, 0
        for i, n in enumerate(self.reaction.num_reactant_atoms):
            reactant_subset_indices[i] = xrange(n_count, n_count +n)
            n_count += n
        for j, m in enumerate(self.reaction.num_product_atoms):
            product_subset_indices[j] = xrange(m_count, m_count + m)
            m_count += m

        # TODO remove on release
        assert(n_count == self.reactants_elements.size)
        assert(m_count == self.products_elements.size)

        return reactant_subset_indices, product_subset_indices

    def check_size_consistency(self):
        """
        Check that there's the same elements in reactants and products

        """
        # TODO add slack and remove the need for this

        if np.array_equal(np.sort(self.reactants_elements), np.sort(self.products_elements)) == False:
            eprint(1, "ERROR: It is currently required that there's the same elements in \
                       reactants and products")
            quit()

    def initialize_match_matrix(self):
        """
        Construct the initial match matrix.

        Returns:
        --------
        match_matrix: array
            The match matrix

        """
        # TODO add possibility for slack

        match_matrix = np.zeros((self.reactants_elements.size, self.products_elements.size))

        # set sub blocks of the match matrix to one plus random pertubation
        # followed by column normalization
        for indices in self.element_type_subset_indices:
            match_matrix[indices] = 1 + 1e-6 * np.random.random(match_matrix[indices].shape)
            match_matrix[indices] /= match_matrix[indices].sum(0)

        #match_matrix = np.eye(match_matrix.shape[0])
        #for i,j in [(0,0),(0,4),(1,1),(1,3),(4,0),(4,4),(3,3),(3,1),(7,7),(7,11),(8,8),(8,10),(20,20),(20,24),(21,21),(21,23),(11,7),(11,11),(10,8),(10,10),(24,20),(24,24),(23,23),(23,21)]:
        #    match_matrix[i,j] = 0.5

        return match_matrix

    def row_dominance(self, return_indices = False):
        """
        Checks if assigning the largest value in a row leads to a unique assignment

        Parameters:
        -----------
        M: ndarray
            Assignment probability matrix.
        return_indices: bool
            Whether to return the indices of the largest values in each row (True),
            else return whether there's row dominance.

        """

        M = self.match_matrix

        # get indices of the largest values in each row
        largest_row_indices = M.argmax(1)
        if return_indices:
            return largest_row_indices
        # TODO add possibility for slack
        unique_indices, counts = np.unique(largest_row_indices, return_counts=True)
        if (counts == 1).all() == False:
            return False
        return True

    def softassign(self):
        """
        Run the softassign algorithm until convergence.

        """
        # TODO add possibility of slack

        for i, indices in enumerate(self.element_type_subset_indices):
            M = self.match_matrix[indices]
            for it in xrange(self.max_softassign_iterations):
                # normalize across rows (except slack)
                M /= np.sum(M,axis=1)[:,None]
                # normalize across columns (except slack)
                M /= np.sum(M,axis=0)

                max_row_normalization_error = np.max(abs(np.sum(M, axis = 1)-1))
                # break if converged
                if max_row_normalization_error < self.softassign_convergence_threshold:
                    oprint(5, "Softassign algorithm for subset %d converged in iteration %d" % (i, it+1))
                    break

                if it == (self.max_softassign_iterations - 1):
                    eprint(3, "WARNING: Softassign algorithm for subset %d did not converge to %.2g (reached %.2g) in %d iterations" % (i, self.softassign_convergence_threshold, max_row_normalization_error, self.max_softassign_iterations))

            # M is NOT a view, but a copy
            self.match_matrix[indices] = M


    def set_parameters(self):
        """
        Returns:
        --------
        initial_inverse_temperature: float
            start inverse temperature
        final_inverse_temperature: float
            end inverse temperature
        inverse_temperature_increment: float
            annealing rate
        max_annealing_iterations: float
            max iterations in annealing loop
        max_relaxation_iterations: float
            max iterations in relaxation loop
        max_softassign_iterations: float
            max iterations in softassign loop
        annealing_convergence_threshold: float
            convergence criteria for the match_matrix l2-norm in the annealing loop
        relaxation_convergence_threshold: float
            convergence criteria ?
        softassign_convergence_threshold: float
            convergence criteria for the maximum deviation from 1 in the match_matrix l1-norm in the softassign loop

        """

        initial_inverse_temperature = settings.initial_inverse_temperature
        final_inverse_temperature = settings.final_inverse_temperature
        max_annealing_iterations = settings.max_annealing_iterations
        max_relaxation_iterations = settings.max_relaxation_iterations
        max_softassign_iterations = settings.max_softassign_iterations
        annealing_convergence_threshold = settings.annealing_convergence_threshold
        relaxation_convergence_threshold = settings.relaxation_convergence_threshold
        softassign_convergence_threshold = settings.softassign_convergence_threshold

        if settings.annealing_method == "multiplication":
            inverse_temperature_increment = (final_inverse_temperature / float(initial_inverse_temperature)) ** (1.0 / max_annealing_iterations)
        elif settings.annealing_method == "addition":
            inverse_temperature_increment = (final_inverse_temperature - initial_inverse_temperature) / float(max_annealing_iterations)
        else:
            eprint(1, "ERROR: Annealing method %s not recognized" % settings.annealing_method)
            quit()

        return (initial_inverse_temperature, final_inverse_temperature, inverse_temperature_increment,
                max_annealing_iterations, max_relaxation_iterations, max_softassign_iterations,
                annealing_convergence_threshold, relaxation_convergence_threshold, softassign_convergence_threshold)

    def annealing(self):
        """
        Begin the actual ordering by deterministic annealing.

        """
        oprint(3, "Beginning deterministic annealing")

        # TODO annealing b_n = b_0 * r**n (t_n = t_0 * r**-n) vs b_n = b_0 + n*r (t_n = t_0/(1 + t_0*n*r))
        # TODO allow soft assignment of atom types
        # TODO allow outliers
        # TODO add fix for symmetry when there's no rotation term / use rotation only to break symmetry
        # TODO add convergence criteria for change in match matrix
        # sum(Mai) term - anti sparseness ?
        # sum(Mai(1-Mai)) term - binaryness ?

        N, M = self.match_matrix.shape
        # initialize inverse temperature
        self.inverse_temperature = self.initial_inverse_temperature

        # begin Deterministic annealing
        for it in xrange(self.max_annealing_iterations):
            self.it = it
            if ((it+1)*100) % self.max_annealing_iterations == 0 and it > 0:
                oprint(4, "On iteration %d of %d in deterministic annealing with convergence at %.2g of %.2g" % (it+1, self.max_annealing_iterations, row_l2_deviation, self.annealing_convergence_threshold))
                print np.sum(self.match_matrix.diagonal())


            self.relaxation()

            row_l2_deviation = abs(1 - np.sum(self.match_matrix**2)/N)

            if row_l2_deviation < self.annealing_convergence_threshold and self.row_dominance():
                oprint(4, "Annealing algorithm converged in iteration %d" % (it+1))
                return True

            self.update_inverse_temperature()

        # The annealing terminated without a definitive assignment
        eprint(3, "WARNING: Annealing algorithm did not converge to %.2g (reached %.2g) in %d iterations and row dominance was %s" % (self.annealing_convergence_threshold, row_l2_deviation, self.max_annealing_iterations, str(self.row_dominance())))
        return False

    def update_inverse_temperature(self):
        """
        Updates inverse temperature dependant on selected annealing scheme

        """

        if settings.annealing_method == "multiplication":
            self.inverse_temperature *= self.inverse_temperature_increment
        elif settings.annealing_method == "addition":
            self.inverse_temperature += self.inverse_temperature_increment
        else:
            eprint(1, "ERROR: Annealing method %s not recognized" % settings.annealing_method)
            quit()

    def relaxation(self):
        """
        Converge the match matrix at the current temperature

        """

        for it in xrange(self.max_relaxation_iterations):

            Q = self.score(self.match_matrix)
            # minimum in Q will be 0, so there's no risk of overflow in the exponential
            # OR IS THERE?
            #scores = [fun.score(self.match_matrix) for fun in self.obj]
            #real_scores = [fun.score(np.eye(16)) for fun in self.obj]
            #for i in range(3):
            #    print i
            #    for j in range(16):
            #        print j
            #        print (self.match_matrix*scores[i])[j]
            #        print (np.eye(16)*real_scores[i])[j]
            #print ""

            #np.exp(-self.inverse_temperature*(Q - min(max(Q.min(0)),max(Q.min(1)))) - 100 , out=self.match_matrix)
            np.exp(-self.inverse_temperature*Q, out=self.match_matrix)
            #print self.match_matrix.sum(0)
            #print self.match_matrix.sum(1)
            #print ""
            if np.isnan(self.match_matrix).any():
                for i in xrange(self.match_matrix.shape[0]):
                    for j in xrange(self.match_matrix.shape[1]):
                        print i,j, self.match_matrix[i,j]
                quit("lol")
            # TODO check if this changes anything
            #self.match_matrix /= self.match_matrix.sum()**0.5

            self.softassign()

            # break if converged
            #mean_squared_difference = np.sum((self.old_match_matrix-self.match_matrix)**2)/self.match_matrix.shape[0]
            mean_squared_difference = np.max(abs(self.old_match_matrix-self.match_matrix))
            if mean_squared_difference < self.relaxation_convergence_threshold:
                oprint(5, "Relaxation algorithm converged in iteration %d" % (it+1))
                break

            if it == self.max_relaxation_iterations-1:
                eprint(3, "WARNING: Relaxation algorithm did not converge to %.2g (reached %.2g) in %d iterations" % (self.relaxation_convergence_threshold, mean_squared_difference, self.max_relaxation_iterations))

            self.backup_match_matrix()

    def backup_match_matrix(self):
        np.copyto(self.old_match_matrix,self.match_matrix)

    def finalize(self):
        """
        Saves output files, prints information etc.

        """

        if self.converged:
            oprint(4, "Match matrix diagonal: " + str(self.match_matrix.diagonal()))
            oprint(4, "Match matrix diagonal sum: " + str(self.match_matrix.diagonal().sum()))
            N = self.match_matrix.shape[0]
            oprint(4, "%d out of %d matched in diagonal" % (sum(self.match_matrix.diagonal() > 0.5), N))
            rd = self.row_dominance(True)
            oprint(4, "Order: " + str(rd))
            M1 = np.eye(N, dtype=bool)
            M2 = np.zeros((N,N), dtype=bool)
            for i in range(N):
                M2[i,rd[i]] = 1

            scores1 = [fun.score(M1) for fun in self.obj]
            scores2 = [fun.score(M2) for fun in self.obj]
            #for i in range(N):
            #    for j in range(M):
            #        if M1[i,j] + M2[i,j] > 0.9:
            #            print i, j, [(scores1[k][i,j],scores2[k][i,j]) for k in range(3)]
            oprint(4, "Scores of diagonal ordering compared to found ordering: " + str([(np.sum(M1*scores1[k]), np.sum(M2*scores2[k])) for k in range(3)]))
        else:
            oprint(4, "Row dominance: " + str(self.row_dominance))

