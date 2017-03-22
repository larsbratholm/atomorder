def quaternion_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD
    """
    P_rot = quaternion_rotate(P, Q)
    return rmsd(P_rot, Q)


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix

    Returns:
    U -- Rotation matrix

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U

def centroid(X):
    """
    Calculate the centroid from a vectorset X
    """
    C = sum(X)/len(X)
    return C

def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)

def write_coordinates(atoms, V):
    """
    Print coordinates V
    """
    N, D = V.shape

    print str(N)
    print

    for i in xrange(N):
        line = "{0:2s} {1:15.8f} {2:15.8f} {3:15.8f}".format(atoms[i], V[i, 0], V[i, 1], V[i, 2])
        print line



# get normal vector to the plane spanned by r1,r2,r3,
# then determine if r4 is above or below the plane
def get_cross_dot(r1,r2,r3,r4):
    v1 = r2 - r1
    v2 = r3 - r1
    v3 = r4 - r1
    n = np.cross(v1,v2)
    return (n.dot(v3) > 0)

# get bond angle with r1 being the central atom
# http://stackoverflow.com/a/13849249
def get_angle(r1,r2,r3):
    v1 = r2 - r1
    v2 = r3 - r1
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1,v2),-1.0,1.0))*180/np.pi



# matrix involved in quaternion rotation

# generate matrices needed in the calculation
def quaternion_rotate(X, Y):
    N = X.shape[0]
    W = np.asarray([makeW(*Y[k]) for k in range(N)])
    P = np.asarray([makeP(*X[k]) for k in range(N)])
    Pt_dot_W = np.asarray([np.dot(P[k].T,W[k]) for k in range(N)])
    W_minus_P = np.asarray([W[k] - P[k] for k in range(N)])
    C1 = -np.sum(Pt_dot_W,axis=0)
    C2 = 0.5*N
    C3 = np.sum(W_minus_P,axis=0)
    A = np.dot(C3.T,C3)*C2-C1
    eigen = np.linalg.eigh(A)
    r = eigen[1][:,eigen[0].argmax()]
    rot = quaternion_transform(r)
    return np.dot(X,rot)



# from http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
def rand_rotation_matrix(deflection=1.0, randnums=None):
    """
    Creates a random rotation matrix.

    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
    rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c

    if randnums is None:
        randnums = np.random.uniform(size=(3,))

    theta, phi, z = randnums

    theta = theta * 2.0*deflection*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0*deflection  # For magnitude of pole deflection.

    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.

    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
        )

    st = np.sin(theta)
    ct = np.cos(theta)

    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))

    # Construct the rotation matrix  ( V Transpose(V) - I ) R.

    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M


def closest_nodes(nodes):
    nodes = np.asarray(nodes)
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)

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
#if cut:
#    # weights is an attempt to make eps3 an averaged limit
#    weights = np.zeros(M.shape)
#    for a in range(M.shape[0]):
#        for i in range(M.shape[1]):
#            weights[a,i] = max(np.sum(~mask[a]), np.sum(~mask[:,i]))
#    weps = cut_eps / weights
#    # NOTE: do this?
#    w = np.where((M < weps) != mask)
#    if w[0].size > 0:
#        # update the idx_lookup to remove these indices
#        # since they are set to zero anyways
#        # first pop ai
#        for key in zip(*w):
#            idx_lookup.pop(key)
#        # then pop bj's
#        while True:
#            flag = True
#            # iterate through all items
#            for key, value in idx_lookup.items():
#                # find any bj's that should be removed
#                w = np.where((v in idx_lookup.keys() for v in value))[0]
#                # if all bj's for a give ai is removed, pop the ai key
#                if w.size == len(value):
#                    idx_lookup.pop(key)
#                    flag = False
#                # else just pop the individual bj
#                elif w.size > 0:
#                    for i in w[::-1]:
#                        value.pop(i)
#            # if no keys were removed, break
#            if flag:
#                break

#        mask[np.where(M < weps)] = True
#def angle_from_points(p):
#    """
#    Angle between three points with the first point being the central one.
#
#    Parameters:
#    -----------
#    p: array-like
#        tuple or array of three points in carteesian coordinates
#
#    """
#
#    p1,p2,p3 = p
#    d12 = sum((p1-p2)**2)
#    d13 = sum((p1-p3)**2)
#    d23 = sum((p2-p3)**2)
#    max_idx = np.argmax([d12,d13,d23])
#    if max_idx == 0:
#        v1 = p1-p3
#        v2 = p2-p3
#    elif max_idx == 1:
#        v1 = p1 - p2
#        v2 = p3 - p2
#    else:
#        v1 = p2 - p1
#        v2 = p3 - p1
#    return angle_v(v1,v2)*(180/np.pi)
