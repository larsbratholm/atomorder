"""
atomorder/objectives.py

all scoring functions (objectives)

"""

import numpy as np
import scipy.optimize
from . import settings
import itertools
from .utils import eprint, oprint
from .parsers import write_xyz, write_mol2

# TODO options ignore monovalent / ignore hydrogens

class Rotation(object):
    """
    Rotation(M)

    Objective based on how well reactants and products can be aligned to each other.

    Based on Gold et al. 1998 doi:10.1016/S0031-3203(98)80010-1 and
    Walker et al. 1991 doi:10.1016/1049-9660(91)90036-O

    Parameters:
    -----------
    M: object
        Ordering object

    Attributes:
    -----------
    M: object
        Ordering object
    X: array
        N*3 sized array with reactant coordinates
    Y: array
        N*3 sized array with product coordinates
    reactant_hydrogen_mask: array
        bool array where all reactant non-hydrogen atoms are True
    product_hydrogen_mask: array
        bool array where all product non-hydrogen atoms are True
    W: array
        Page 363 Walker91
    Q: array
        Page 363 Walker91
    Qt_dot_W: array
        Page 363 Walker91
    W_minus_Q: array
        Page 363 Walker91
    q: array
        Each rotation and transformation is defined by a dual quaternion.
        q is the set of all quaternions, where q[0,0] = s1, q[0,1] = r1, q[1,0] = s2 etc.
    # TODO update


    """

    def __init__(self, M):
        oprint(3, "Initializing rotation objective")
        self.M = M
        self.X = self.M.reactants_coordinates
        self.Y = self.M.products_coordinates
        self.reactant_hydrogen_mask = (self.M.reactants_elements == "H") | (self.M.reactants_elements == "D")
        self.product_hydrogen_mask = (self.M.products_elements == "H") | (self.M.products_elements == "D")
        self.score = self.set_solver()
        self.q = self.initialize_quaternions()
        self.W, self.Q, self.Qt_dot_W, self.W_minus_Q = self.interaction_arrays()
        # only used in numeric part
        self.squared_distances = np.zeros(self.M.match_matrix.shape, dtype=float)

    def initialize_quaternions(self):
        """
        Each transformation (rotation+translation) is described by a dual quaternion pair (s,r).
        self.q holds all dual quaternions needed for all the transformations, such that
        self.q = np.concatenate([s1,r1,s2,r2,...]).
        For N reactants and M products, (N+M-1) transformations is required.

        """
        # initialize as s = (0,0,0,0) and r = (0,0,0,1)
        N, M = self.M.num_reactants, self.M.num_products
        # only use three elements and enforce constraints in the fourth
        q = np.zeros((N-1, 2, 3))
        #q[:,1,3] = 1
        return q

    def set_solver(self):
        """
        Sets analytical, iterative or numerical solver
        based on the number of reactants and products

        """
        # TODO: solve general case analytically or iterative solver
        if self.M.num_reactants == 1 or self.M.num_products == 1:
            return self.analytical_solver
        else:
            return self.semi_analytical_solver

    def analytical_solver(self, m):
        """
        See page 363 Walker91

        Calculate the optimal translation and rotation between
        all reactants and products, given the current
        match matrix.
        Will only give the correct result when there's either only 1 reactant or 1 product

        Parameters:
        -----------
        match: ndarray
            array of atom order probabilities

        Returns:
        --------
        E: float
            Energy

        """

        #if self.M.it < 0:
        #    return np.zeros(match.shape)
        E = 0
        match = m.copy()

        # Pick first reactant as reference state
        ref_indices = self.M.reactant_subset_indices[0]
        Y0 = self.Y[ref_indices]

        squared_distances = np.zeros(self.M.match_matrix.shape, dtype=float)

        match[np.ix_(self.reactant_hydrogen_mask, self.product_hydrogen_mask)] *= settings.hydrogen_rotation_weight
        all_X  = []
        for i, reactant_indices in enumerate(self.M.reactant_subset_indices):
            X = self.X[reactant_indices]
            for j, product_indices in enumerate(self.M.product_subset_indices):
                Y = self.Y[product_indices]
                sub_matrix_indices = np.ix_(reactant_indices, product_indices)
                match_sub_matrix = match[sub_matrix_indices]
                Qt_dot_W = self.Qt_dot_W[sub_matrix_indices]
                W_minus_Q = self.W_minus_Q[sub_matrix_indices]

                C1 = -2*np.sum(match_sub_matrix[:,:,None,None]*Qt_dot_W,axis=(0,1))
                C2 = match_sub_matrix.sum()
                C3 = 2*np.sum(match_sub_matrix[:,:,None,None]*W_minus_Q,axis=(0,1))
                A = 0.5*(C3.T.dot(C3)/(2.0*C2)-C1-C1.T)
                # TODO remove
                assert(np.allclose(A,A.T))
                eigen = np.linalg.eigh(A)
                r = eigen[1][:,-1]
                # TODO remove
                assert(np.allclose(A.dot(r), eigen[0][-1]*r))
                s = -C3.dot(r)/(2.0*C2)
                rot, trans = self.transform(r,s)

                squared_distances[sub_matrix_indices] = np.sum((Y[None,:,:] - trans[None,None,:] - rot.dot(X.T).T[:,None,:])**2, axis=2)
            all_X.append(trans + rot.dot(X.T).T)

        squared_distances[np.ix_(self.reactant_hydrogen_mask, self.product_hydrogen_mask)] *= settings.hydrogen_rotation_weight
        write_mol2(np.concatenate(all_X), self.M.reaction.reactants.atoms, "reactant%d.mol2" % self.M.it,self.M.reaction.num_reactant_atoms)
        write_mol2(self.Y+np.asarray([10,0,0]), self.M.reaction.products.atoms, "product%d.mol2" % self.M.it, self.M.reaction.num_product_atoms, match)

        return squared_distances

    def numerical_solver(self, match):
        """
        See page 363 Walker91

        Numerical version for multiple reactants and products

        Parameters:
        -----------
        match: ndarray
            array of atom order probabilities

        Returns:
        --------
        E: float
            Energy
        q: ndarray
            array of fitted quaternions

        """

        def objective(flat_q, match, self):#, jac):
            # TODO jacobian could be written in matrix form.
            # TODO add penalty for clashing in the numerical version

            # size
            N, M = self.M.num_reactants, self.M.num_products
            # energy
            E = 0

            #if jac:
            #    J = np.zeros(8*(N+M-1))

            q = flat_q.reshape(N+M-1,2,4)
            # There's N+M-1 pairs of r,s.
            # construct rotation matrices and translation vectors
            trans = np.zeros((N+M-1, 3))
            rot = np.zeros((N+M-1, 3, 3))
            for j in xrange(M):
                rot[j], trans[j] = self.transform(q[j, 1], q[j, 0])
            for i in xrange(N-1):
                rot[M+i], trans[M+i] = self.transform(q[M+i,1], q[M+i, 0])

            # use product 1 as reference
            ref_indices = self.M.product_subset_indices[0]
            Y0 = self.Y[ref_indices]
            squared_distances = np.zeros(match.shape)
            all_X = []
            all_Y = []
            for j, reactant_indices in enumerate(self.M.reactant_subset_indices):
                # contribution between reactants and product 1
                X = self.X[reactant_indices]
                #match_sub_matrix = match[np.ix_(ref_indices, product_indices)]
                squared_distances[np.ix_(reactant_indices, ref_indices)] = np.sum((Y0[None,:,:] - trans[j][None,None,:] - rot[j].dot(X.T).T[:,None,:])**2, axis=2)


                # contributions between products and remaining reactants
                for i, product_indices in enumerate(self.M.product_subset_indices[1:]):
                    Y = self.Y[product_indices]
                    #match_sub_matrix = match[np.ix_(reactant_indices, product_indices)]
                    squared_distances[np.ix_(reactant_indices, product_indices)] = np.sum((trans[M+i][None,None,:] + (rot[M+i].dot(Y.T).T)[None,:,:] - trans[j][None,None,:] - (rot[j].dot(X.T).T)[:,None,:])**2, axis=2)

            for j, reactant_indices in enumerate(self.M.reactant_subset_indices):
                all_X.append(trans[j] + rot[j].dot(self.X[reactant_indices].T).T)
            all_Y.append(Y0)
            for i, reactant_indices in enumerate(self.M.reactant_subset_indices[1:]):
                all_Y.append(trans[M+i] + rot[M+i].dot(self.Y[product_indices].T).T)
            write_xyz(np.concatenate(all_X), self.M.reactants_elements, "reactant%d.xyz" % self.M.it, str(np.sum(match*squared_distances)))
            write_xyz(np.concatenate(all_Y), self.M.products_elements, "product%d.xyz" % self.M.it, str(np.sum(match*squared_distances)))
            E = np.sum(match*squared_distances)
            self.squared_distances = squared_distances.copy()

            return E


        def rr_constraint_jacobian(x, j, N, M):
            # jacobian for use with SLSQP
            jac = np.zeros(8*(M+N-1))
            jac[8*j+4:8*j+8] = 2*x[8*j+4:8*j+8]
            return jac

        def rs_constraint_jacobian(x, j, N, M):
            # jacobian for use with SLSQP
            jac = np.zeros(8*(M+N-1))
            jac[8*j:8*j+4] = x[8*j+4:8*j+8]
            jac[8*j+4:8*j+8] = x[8*j:8*j+4]
            return jac


        N, M = self.M.num_reactants, self.M.num_products
        # create constraints
        cons = np.empty(2*(N+M-1), dtype=dict)

        for j in range(M+N-1):
            cons[2*j]   = {"type": "eq", "fun": lambda x: np.sum(x[8*j+4:8*j+8]*x[8*j+4:8*j+8]) -1} # r.T*r = 1
            cons[2*j+1] = {"type": "eq", "fun": lambda x: np.sum(x[8*j+4:8*j+8]*x[8*j:8*j+4])} # r.T*s = 0

            # jacobians for SLSQP
            cons[2*j]["jac"] = lambda x: rr_constraint_jacobian(x, j, N, M)
            cons[2*j+1]["jac"] = lambda x: rs_constraint_jacobian(x, j, N, M)


        bounds = []
        for j in range(M+N-1):
            bounds.extend([(None, None)]*4)
            bounds.extend([(-1, 1)]*4)

        q = np.zeros((N+M-1, 2, 4))
        q[:,1,3] = 1

        objective(q, match, self)
        opt = scipy.optimize.minimize(objective, q, constraints=cons, method="SLSQP", options={"maxiter": 500, "disp": 0, "ftol": 1e-4}, args=(match, self), bounds=bounds)
        squared_distances = self.squared_distances.copy()
        #self.update_quaternions(opt.x.reshape((N+M-1),2,4))
        return self.squared_distances

    def semi_analytical_solver(self, match):
        """
        See page 363 Walker91

        Semi-numerical version for multiple reactants and products

        Parameters:
        -----------
        match: ndarray
            array of atom order probabilities

        """

        def objective(flat_q, m, self):
            match = m.copy()

            # size
            N, M = self.M.num_reactants, self.M.num_products
            squared_distances = np.zeros(match.shape)
            # energy
            E = 0

            # reshape
            q = flat_q.reshape(N-1,2,3)
            # add restraints as in the fourth element in the quaternions
            q = np.concatenate([q, np.zeros((N-1,2,1))], axis=2)

            # There's N-1 pairs of r,s that has to be solved numerically
            # construct rotation matrices and translation vectors
            trans = np.zeros((N-1, 3))
            rot = np.zeros((N-1, 3, 3))
            for j in xrange(N-1):
                # r restraints
                r4_sq = 1 - q[j,1,0]**2 - q[j,1,1]**2 - q[j,1,2]**2
                if r4_sq < 0:
                    return np.inf
                q[j,1,3] = np.sqrt(r4_sq)
                # s restraints
                s4 = - (q[j,0,0]*q[j,1,0] + q[j,0,1]*q[j,1,1] + q[j,0,2]*q[j,1,2])/q[j,1,3]
                rot[j], trans[j] = self.transform(q[j, 1], q[j, 0])

            ## use product 1 as reference
            #match[np.ix_(self.reactant_hydrogen_mask, self.product_hydrogen_mask)] *= settings.hydrogen_rotation_weight
            #for i, reactant_indices in enumerate(self.M.reactant_subset_indices):
            #    X = self.X[reactant_indices]
            #    for j, product_indices in enumerate(self.M.product_subset_indices):
            #        # solve everything other than internal product rotations analytically
            #        if j == 0:
            #            Y = self.Y[product_indices]
            #        else:
            #            Y = trans[j-1] + rot[j-1].dot(self.Y[product_indices].T).T

            #        sub_matrix_indices = np.ix_(reactant_indices, product_indices)
            #        match_sub_matrix = match[sub_matrix_indices]
            #        #Qt_dot_W = self.Qt_dot_W[sub_matrix_indices]
            #        #W_minus_Q = self.W_minus_Q[sub_matrix_indices]
            #        Q = np.asarray([self.makeQ(*y) for y in Y])
            #        Qt_dot_W = np.asarray([[np.dot(q.T,w) for q in Q] for w in self.W[reactant_indices]])
            #        W_minus_Q = np.asarray([[w - q for q in Q] for w in self.W[reactant_indices]])

            #        C1 = -2*np.sum(match_sub_matrix[:,:,None,None]*Qt_dot_W,axis=(0,1))
            #        C2 = match_sub_matrix.sum()
            #        C3 = 2*np.sum(match_sub_matrix[:,:,None,None]*W_minus_Q,axis=(0,1))
            #        A = 0.5*(C3.T.dot(C3)/(2.0*C2)-C1-C1.T)
            #        # TODO remove
            #        assert(np.allclose(A,A.T))
            #        eigen = np.linalg.eigh(A)
            #        r = eigen[1][:,-1]
            #        # TODO remove
            #        assert(np.allclose(A.dot(r), eigen[0][-1]*r))
            #        s = -C3.dot(r)/(2.0*C2)
            #        xrot, xtrans = self.transform(r,s)

            #        #all_X.append(xtrans + xrot.dot(X.T).T)
            #        #all_Y.append(Y)

            #        squared_distances[sub_matrix_indices] = np.sum((Y[None,:,:] - xtrans[None,None,:] - xrot.dot(X.T).T[:,None,:])**2, axis=2)

            #squared_distances[np.ix_(self.reactant_hydrogen_mask, self.product_hydrogen_mask)] *= settings.hydrogen_rotation_weight
            #E += np.sum(m*squared_distances)
            #self.squared_distances = squared_distances.copy()

            Y = []
            for j,  product_indices in enumerate(self.M.product_subset_indices):
                if j == 0:
                    Y.append(self.Y[product_indices])
                else:
                    Y.append(trans[j-1] + rot[j-1].dot(self.Y[product_indices].T).T)
            Y = np.concatenate(Y)
            all_X = []
            for i, reactant_indices in enumerate(self.M.reactant_subset_indices):
                X = self.X[reactant_indices]
                Q = np.asarray([self.makeQ(*y) for y in Y])
                Qt_dot_W = np.asarray([[np.dot(q.T,w) for q in Q] for w in self.W[reactant_indices]])
                W_minus_Q = np.asarray([[w - q for q in Q] for w in self.W[reactant_indices]])
                match_sub_matrix = match[reactant_indices,:]

                C1 = -2*np.sum(match_sub_matrix[:,:,None,None]*Qt_dot_W,axis=(0,1))
                C2 = match_sub_matrix.sum()
                C3 = 2*np.sum(match_sub_matrix[:,:,None,None]*W_minus_Q,axis=(0,1))
                # 1e-9 for stability
                A = 0.5*(C3.T.dot(C3)/(2.0*C2+1e-9)-C1-C1.T)
                # TODO remove
                assert(np.allclose(A,A.T))
                eigen = np.linalg.eigh(A)
                r = eigen[1][:,-1]
                # TODO remove
                assert(np.allclose(A.dot(r), eigen[0][-1]*r))
                s = -C3.dot(r)/(2.0*C2)
                xrot, xtrans = self.transform(r,s)

                #all_X.append(xtrans + xrot.dot(X.T).T)
                #all_Y.append(Y)
                all_X.append(xtrans + xrot.dot(X.T).T)

                squared_distances[reactant_indices,:] = np.sum((Y[None,:,:] - xtrans[None,None,:] - xrot.dot(X.T).T[:,None,:])**2, axis=2)
            for i, v in enumerate(all_X):
                all_X[i] += np.asarray([1*i,0,0])
            #write_xyz(np.concatenate(all_X), self.M.reactants_elements, "reactant%d.xyz" % self.M.it, str(np.sum(match*squared_distances)))
            #write_xyz(Y, self.M.products_elements, "product%d.xyz" % self.M.it, str(np.sum(match*squared_distances)))
            write_mol2(np.concatenate(all_X), self.M.reaction.reactants.atoms, "reactant%d.mol2" % self.M.it,self.M.reaction.num_reactant_atoms)
            write_mol2(Y+np.asarray([10,0,0]), self.M.reaction.products.atoms, "product%d.mol2" % self.M.it, self.M.reaction.num_product_atoms, match)

            squared_distances[np.ix_(self.reactant_hydrogen_mask, self.product_hydrogen_mask)] *= settings.hydrogen_rotation_weight
            E += np.sum(m*squared_distances)
            self.squared_distances = squared_distances.copy()
            #for i, xi in enumerate(all_X):
            #    for xj in all_X[i+1:]:
            #        clash = (np.sum((xi[:,None,:]-xj[None,:,:])**2, axis=2)**0.5 < 1.2)
            #        if clash.any():
            #            d = (np.sum((xi[:,None,:]-xj[None,:,:])**2, axis=2)**0.5)[clash].ravel()
            #            E += np.sum((d-1.2)**2)*10000

            #        #clash = (np.sum((X[~self.reactant_hydrogen_mask[reactant_indices]][:,None,:] - XK[~self.reactant_hydrogen_mask[K]][None,:,:])**2, axis=2)**0.5 < 1.8)
            #        #if clash.any():
            #        #    quit()
            #        #    E += clash.sum()*1000
            #print E
            return E

        N, M = self.M.num_reactants, self.M.num_products
        bounds = []
        for j in range(N-1):
            bounds.extend([(None, None)]*3)
            bounds.extend([(-1, 1)]*3)
        q = self.q.copy().ravel()
        #q = np.zeros(q.shape)
        opt = scipy.optimize.minimize(objective, q, method="l-bfgs-b", options={"maxiter": 500, "disp": 0, "ftol": 1e-6}, args=(match, self), bounds=bounds)
        #opt = scipy.optimize.minimize(objective, q, method="nelder-mead", options={"maxiter": 500, "disp": 0, "ftol": 1e-6}, args=(match, self))
        #assert(np.allclose(self.squared_distances, self.analytical_solver(match)))
        self.update_quaternions(opt.x.reshape((N-1),2,3))
        return self.squared_distances 

    def update_quaternions(self, q):
        self.q = q.copy()

    def interaction_arrays(self):
        """
        Generate W, Q, Qt_dot_W, W_minus_Q arrays from page 363 Walker91
        disregarding weights (for now)

        """
        # TODO: add ignore hydrogens for calculating the rotations
        W = np.asarray([self.makeW(*x) for x in self.X])
        Q = np.asarray([self.makeQ(*y) for y in self.Y])
        Qt_dot_W = np.asarray([[np.dot(q.T,w) for q in Q] for w in W])
        W_minus_Q = np.asarray([[w - q for q in Q] for w in W])
        return W, Q, Qt_dot_W, W_minus_Q

    def makeW(self, r1,r2,r3,r4=0):
        # eqn 16 Walker91
        W = np.asarray([
                 [r4, r3, -r2, r1],
                 [-r3, r4, r1, r2],
                 [r2, -r1, r4, r3],
                 [-r1, -r2, -r3, r4] ])
        return W

    def makeQ(self, r1,r2,r3,r4=0):
        # eqn 15 Walker91
        Q = np.asarray([
                 [r4, -r3, r2, r1],
                 [r3, r4, -r1, r2],
                 [-r2, r1, r4, r3],
                 [-r1, -r2, -r3, r4] ])
        return Q

    def transform(self, r, s):
        Wt_r = self.makeW(*r).T
        Q_r = self.makeQ(*r)
        rot = Wt_r.dot(Q_r)[:3,:3]
        trans = Wt_r.dot(s)[:3]
        return rot, trans

class Atomic(object):
    """
    Atomic(M)

    Atomic matching contribution

    """

    def __init__(self, M):
        oprint(3, "Initializing atomic objective")
        self.M = M
        self.score_matrix = self.get_score_matrix()

    def get_score_matrix(self):
        # TODO expand on this

        # the score matrix is created such that a perfect match is 0
        # and imperfect matches are positive
        score_matrix = np.zeros(self.M.match_matrix.shape)

        # sybyl match matrix
        score_matrix += settings.atomic_sybyl_weight * (self.M.reactants_sybyl_types[:, None] != self.M.products_sybyl_types[None,:])

        # enforce that elements only match other elements of same type

        score_matrix += 1e6 * (self.M.reactants_elements[:, None] != self.M.products_elements[None,:])

        return score_matrix

    def score(self, match):
        """
        Return the one bond scoring matrix

        Parameters:
        -----------
        match: array
            match_matrix

        Returns:
        --------
        Score matrix

        """

        # match matrix is not used for atomic terms, but is included as a parameter for convenience

        return self.score_matrix

class Bond(object):
    """
    Bond(M)

    Bond matching contribution

    """

    def __init__(self, M):
        oprint(3, "Initializing bond objective")
        self.M = M
        self.matrix_function = self.get_matrix_function()
        self.evaluate_score_matrix_vect = np.vectorize(self.evaluate_score_matrix, otypes=[np.float], excluded=['match'])

    def get_matrix_function(self):

        reactant_elements = self.M.reaction.reactants.element_symbols
        product_elements = self.M.reaction.products.element_symbols
        reactant_bond_matrix = self.M.reaction.reactants.bond_matrix
        product_bond_matrix = self.M.reaction.products.bond_matrix

        # Construct a matrix function for easy evaluation of the objective
        C = np.empty(self.M.match_matrix.shape, dtype = object)
        # Fill it with functions that returns zero
        for a,i in itertools.product(xrange(reactant_elements.size),xrange(product_elements.size)):
            C[a,i] = lambda x: 0

        # temporary helper dict to construct the matrix
        pairs = {}

        # fill pairs with all possible bond matches
        for a,b in zip(*np.where(reactant_bond_matrix)):
            element_a, element_b = reactant_elements[[a,b]]
            for i,j in zip(*np.where(product_bond_matrix)):
                element_i, element_j = reactant_elements[[i,j]]
                if (element_a == element_i and element_b == element_j) or (element_a == element_j and element_b == element_i):
                    if (a,i) not in pairs: pairs[(a,i)] = []
                    pairs[a,i].append((b,j))

        # construct the matrix function
        for (a,i), container in pairs.items():
            # y=container makes sure that the values of b,j can be used outside the local scope
            C[a,i] = lambda x, y=container: np.sum((x[b,j] for (b,j) in y))

        return C

    def evaluate_score_matrix(self, a, i, match):
        return self.matrix_function[a,i](match)

    def score(self, match):
        """
        Return the one bond scoring matrix

        Parameters:
        -----------
        match: array
            match_matrix

        Returns:
        --------
        score_matrix

        """

        score_matrix = -settings.bond_weight * np.fromfunction(self.evaluate_score_matrix_vect, match.shape, match=match, dtype=int)

        return score_matrix

