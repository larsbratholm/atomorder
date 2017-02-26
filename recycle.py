#def quaternion_rmsd(P, Q):
#    """
#    Rotate matrix P unto Q and calculate the RMSD
#    """
#    P_rot = quaternion_rotate(P, Q)
#    return rmsd(P_rot, Q)
#
#def rotate(X, Y, M):
#
#    # algorithm from http://dx.doi.org/10.1016/S0031-3203(98)80010-1
#    def quaternion_transform(r, s):
#        Wt_r = makeW(*r).T
#        P_r = makeP(*r)
#        rot = Wt_r.dot(P_r)[:3,:3]
#        trans = Wt_r.dot(s)[:3]
#        return rot, trans
#
#    X -= centroid(X)
#    Y -= centroid(Y)
#
#    Nx = X.shape[0]
#    Ny = Y.shape[0]
#    # generate W, P, Pt_dot_W, W_minus_P arrays
#    W = np.asarray([makeW(*X[j]) for j in range(Nx)])
#    P = np.asarray([makeP(*Y[k]) for k in range(Ny)])
#    Pt_dot_W = np.asarray([[np.dot(P[k].T,W[j]) for k in xrange(Ny)] for j in xrange(Nx)])
#    W_minus_P = np.asarray([[W[j] - P[k] for k in xrange(Ny)] for j in xrange(Nx)])
#
#    # create quaternion representation of P and Q by adding a zero column
#    x = np.concatenate([X,np.zeros((Nx,1))],axis=1)
#    y = np.concatenate([Y,np.zeros((Ny,1))],axis=1)
#
#
#    C1 = -np.sum(M[:,:,None,None]*Pt_dot_W,axis=(0,1))
#    msum = M.sum()
#    C3 = np.sum(M[:,:,None,None]*W_minus_P,axis=(0,1))
#    # find optimal rotation and translation
#    A = 0.5*C3.T.dot(C3)/msum-C1
#    eigen = np.linalg.eigh(A)
#    r = eigen[1][:,-1]
#    s = -C3.dot(r)/msum
#    rot, trans = quaternion_transform(r,s)
#    distances = np.sum((Y[None,:,:] - trans[None,None,:] - rot.dot(X.T).T[:,None,:])**2,axis=2)
#    return distances
#
#def kabsch(P, Q):
#    """
#    The optimal rotation matrix U is calculated and then used to rotate matrix
#    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
#    calculated.
#
#    Using the Kabsch algorithm with two sets of paired point P and Q,
#    centered around the center-of-mass.
#    Each vector set is represented as an NxD matrix, where D is the
#    the dimension of the space.
#
#    The algorithm works in three steps:
#    - a translation of P and Q
#    - the computation of a covariance matrix C
#    - computation of the optimal rotation matrix U
#
#    http://en.wikipedia.org/wiki/Kabsch_algorithm
#
#    Parameters:
#    P -- (N, number of points)x(D, dimension) matrix
#    Q -- (N, number of points)x(D, dimension) matrix
#
#    Returns:
#    U -- Rotation matrix
#
#    """
#
#    # Computation of the covariance matrix
#    C = np.dot(np.transpose(P), Q)
#
#    # Computation of the optimal rotation matrix
#    # This can be done using singular value decomposition (SVD)
#    # Getting the sign of the det(V)*(W) to decide
#    # whether we need to correct our rotation matrix to ensure a
#    # right-handed coordinate system.
#    # And finally calculating the optimal rotation matrix U
#    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
#    V, S, W = np.linalg.svd(C)
#    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
#
#    if d:
#        S[-1] = -S[-1]
#        V[:, -1] = -V[:, -1]
#
#    # Create Rotation matrix U
#    U = np.dot(V, W)
#
#    return U
#
#def centroid(X):
#    """
#    Calculate the centroid from a vectorset X
#    """
#    C = sum(X)/len(X)
#    return C
#
#def rmsd(V, W):
#    """
#    Calculate Root-mean-square deviation from two sets of vectors V and W.
#    """
#    D = len(V[0])
#    N = len(V)
#    rmsd = 0.0
#    for v, w in zip(V, W):
#        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
#    return np.sqrt(rmsd/N)
#
#def write_coordinates(atoms, V):
#    """
#    Print coordinates V
#    """
#    N, D = V.shape
#
#    print str(N)
#    print
#
#    for i in xrange(N):
#        line = "{0:2s} {1:15.8f} {2:15.8f} {3:15.8f}".format(atoms[i], V[i, 0], V[i, 1], V[i, 2])
#        print line
#
#
## bi: start inverse temperature
## bf: end inverse temperature
## br: annealing rate
## IB: max iterations in B loop
## ID: max iterations in D loop
#def test(atomsP, P, atomsQ, Q, bi = 1e-2, bf = 1e1, br = 1.01, IB = 1000, ID = 1000, epsD=1e-6, eps_gamma = 0.01):#, epsC=1e-6,
#         #cut_eps = 1e-6, cut=0, eps_gamma = 0.01, M_init = 1,
#         #hard_atom = True, has_scipy=1, use_mask=0, offset=0, normalize=0, simple_convergence=0, add_gamma=1, alt_annealing=False, restrict = "CR"):
#    # http://stackoverflow.com/a/13849249
#    def norm(v):
#        return v/np.linalg.norm(v)
#
#    def angle_v(v1,v2):
#        v1_u = norm(v1)
#        v2_u = norm(v2)
#        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
#
#    # p1 central point
#    def angle_p(p):
#        p1,p2,p3 = p
#        d12 = sum((p1-p2)**2)
#        d13 = sum((p1-p3)**2)
#        d23 = sum((p2-p3)**2)
#        max_idx = np.argmax([d12,d13,d23])
#        if max_idx == 0:
#            v1 = p1-p3
#            v2 = p2-p3
#        elif max_idx == 1:
#            v1 = p1 - p2
#            v2 = p3 - p2
#        else:
#            v1 = p2 - p1
#            v2 = p3 - p1
#        return angle_v(v1,v2)*(180/np.pi)
#    # check that assigning the largest value in a row doesn't break
#    def row_dominance(m, output = False):
#        largest_row_indices = m.argmax(1)
#        # remove the outlier indices
#        #largest_row_indices = largest_row_indices[largest_row_indices != m.shape[1]-1]
#        unique_indices, counts = np.unique(largest_row_indices, return_counts=True)
#        if output:
#            return largest_row_indices
#        if (counts == 1).all():
#            return True
#        else:
#            return False
#
#
#
#    ## possibilities
#    # soft/hard atom-atom requirement
#    # arbitrary delta or calced delta
#    # sum beta or multiply beta
#    # soft outliers or hard outliers
#    # analytisk loesning / sparse numerisk
#    # cutoff mask
#    # large molecule support
#    # scipy
#    # score dot product or null/1/-2
#    # sum(Mai) term - anti sparseness
#    # sum(Mai(1-Mai)) term - binaryness
#
#    molP = molecule(atomsP,P)
#    molQ = molecule(atomsQ,Q)
#
#    unique_atypes_P, P_counts = np.unique(molP.atypes, return_counts=1)
#    unique_atypes_Q, Q_counts = np.unique(molQ.atypes, return_counts=1)
#
#    # TODO expand for not same atoms
#    assert((unique_atypes_P == unique_atypes_Q).all())
#    assert((P_counts == P_counts).all())
#
#    # mask is True where atomtypes differ between the two molecules
#    # NOTE: size N+1 instead of N due to dummy for outliers (not implemented)
#    mask = np.zeros((molP.size,molQ.size), dtype=bool)
#    mask[np.where(atomsP[:,None] != atomsQ[None,:])] = True
#    #mask[-1,-1] = True
#    ## NOTE: temporary
#    #mask[-1,:] = True
#    #mask[:,-1] = True
#
#    # softassign definitions
#    Q = np.zeros((molP.size,molQ.size),dtype=float)
#    # Set N to the larger of two molecule sizes plus 1 for outliers
#    # N = max(molP.size, molQ.size)
#    N = molP.size
#    # assignment matrix
#    #if M_init == 0:
#    #    M = 1 + np.random.random((molP.size,molQ.size))*1e-6
#    #else:
#    f = (molP.atypes[:,None] == molQ.atypes[None,:])
#    #M = N**-1 + np.random.random((N,N))*1e-3
#    M = f*(1+np.random.random(f.shape)*1e-3/N)
#
#    # find non zero indices for scoring
#    w1 = np.where(molP.bonds)
#    w2 = np.where(molQ.bonds)
#    # weird behavior of itertools.izip, so using zip
#    bonds1 = zip(*w1)
#    bonds2 = zip(*w2)
#
#    idx_lookup = {}
#    for a,b in bonds1:
#        for i,j in bonds2:
#            if molP.atypes[a] == molQ.atypes[i] and molP.atypes[b] == molQ.atypes[j]:
#                if (a,i) not in idx_lookup: idx_lookup[(a,i)] = []
#                idx_lookup[(a,i)].append((b,j))
#
#    #if has_scipy:
#    C = scipy.sparse.dok_matrix((N**2,N**2), dtype=float)
#    #else:
#    #    C = np.zeros((N**2, N**2), dtype=float)
#
#    for (a,i), v in idx_lookup.items():
#        for (b,j) in v:
#            C[(a*N+i, b*N+j)] = 1
#    #if has_scipy:
#    C = C.tocsr()
#    assert((C!=C.T).nnz == 0)
#    #else:
#    #    assert(np.allclose(C,C.T))
#
#    #if restrict == "CR":
#    #r1 = np.eye(N)-N**(-1) * np.ones((N,N))
#    #r2 = r1
#    #elif restrict == "C":
#    #    r1 = np.eye(N)-N**(-1) * np.ones((N,N))
#    #    r2 = np.eye(N)
#    #else:
#    r1 = np.eye(N)
#    r2 = r1
#
#    #if has_scipy:
#    R = scipy.sparse.kron(r1,r2) # only with column restriction
#    R = ((C.dot(R).T).dot(R.T)).T
#    l = scipy.sparse.linalg.eigsh(R, k=1, return_eigenvectors=False, which="SA", tol=1e-4)[0]
#    #else:
#    #    R = np.kron(r1,r2)
#    #    R = ((C.dot(R).T).dot(R.T)).T
#    #    l = np.linalg.eigvalsh(R)[0]
#
#
#    #if add_gamma:
#    #gamma = -l + eps_gamma
#    #else:
#    gamma = 0
#
#    #con = 0
#    #for j in range(N):
#    #    maxj = 0
#    #    for a,b in itertools.product(*([xrange(N)]*2)):
#    #        #diff = np.max(C[a*N+i, :] - C[b*N+i,:])
#    #        diff = np.max(molP.bonds[a,b]*molQ.bonds*f[a,:,None]*f[b,None,:])
#    #        if diff > maxj:
#    #            maxj = diff
#    #    con += maxj
#
#    #con2 = 0
#    #for j in range(molQ.size):
#    #    con2 += max(molQ.bonds[:,j].max()*max([molP.bonds[:,c].max() - molP.bonds[:,c].min() for c in range(molP.size)]),
#    #                molQ.bonds[:,j].min()*min([molP.bonds[:,c].min() - molP.bonds[:,c].max() for c in range(molP.size)]))
#    #assert(con == con2)
#
#    con3 = np.sum(np.max(molQ.bonds[None,None,None,:,:]*f[None,None,:,None,:]*(molP.bonds[:,None,:,None,None]*f[:,None,None,:,None]-molP.bonds[None,:,:,None,None]*f[None,:,None,:,None]), axis=(0,1,2,3)))
#    con = con3
#    #assert(con3 == con)
#    #con4 = 0
#    #v = (molP.bonds[:,None,:,None]*f[:,None,None,:]-molP.bonds[None,:,:,None]*f[None,:,None,:])[:,:,:,:,None]
#    #for j in range(molQ.size):
#    #    con4 += np.max(molQ.bonds[None,None,None,:,j]*f[None,None,:,None,j]*v)
#
#    #assert(con4 == con)
#
#    # inverse temperature
#    b = bi
#    Mlast = M.copy()
#    flag = False
#    # begin A: Deterministic annealing
#    for itA in range(10000000): # avoid while loop
#        if b > bf: # annealing stop
#            print M
#            print 1 - np.sum(M**2)/N < 0.1
#            print row_dominance(M)
#            print M.max(0)
#            break
#        if itA > 0 and 1 - np.sum(M**2)/N < 0.1:
#            if row_dominance(M):
#                print "pewpew"
#                print M.diagonal()
#                rd = row_dominance(M,True)
#
#                M1 = np.eye(N)
#                M2 = np.zeros((N,N))
#                for i in range(N):
#                    M2[i,rd[i]] = 1
#                #np.random.shuffle(ASD)
#                E1M1 = 0
#                E2M1 = 0
#                E4M1 = 0
#                E1M2 = 0
#                E2M2 = 0
#                E4M2 = 0
#                for (a,i),v in idx_lookup.items():
#                    E1M1 += -0.5 * M1[a,i]*(sum((M1[b,j] for b,j in v)))
#                    E2M1 += -0.5 * M1[a,i]*(sum((M1[b,j] for b,j in v)))
#                    E4M1 += -0.5 * M1[a,i]*(sum((M1[b,j] for b,j in v)))
#                    E1M2 += -0.5 * M2[a,i]*(sum((M2[b,j] for b,j in v)))
#                    E2M2 += -0.5 * M2[a,i]*(sum((M2[b,j] for b,j in v)))
#                    E4M2 += -0.5 * M2[a,i]*(sum((M2[b,j] for b,j in v)))
#
#                print E1M1, E1M2,E1#, dE1
#                print E2M1, E2M2,E2#, dE2
#                print E4M1, E4M2,E4#, dE4
#
#
#
#                break
#            else: # check for symmetry
#                # check if there's values that are close to 0.5, which will be the case when there's symmetry
#                w2 = np.where((M - epsD**0.5 < 0.5) & (M + epsD**0.5 > 0.5))
#                pairs = np.asarray(zip(*w2))
#                # check if removing the rows and colums where these are in, we have row dominance
#                mask2 = np.ones((N,N), dtype=bool)
#                mask2[np.unique(w2[0]),:] = False
#                mask2[:,np.unique(w2[1])] = False
#                M2 = M[mask2]
#                M2 = M2.reshape((int(M2.size**0.5),int(M2.size**0.5)))
#                if row_dominance(M2):
#                    distances = rotate(molP.coords, molQ.coords, M)
#                    closest = distances[w2].argmin()
#                    best_pair = pairs[closest]
#                    # get pairs that
#                    idx2 = pairs[np.where(((pairs[:,0] == best_pair[0]) | (pairs[:,1] == best_pair[1])) & (~(pairs == best_pair).all(1)))]
#                    # update the idx_lookup to remove these indices
#                    # since they are set to zero anyways
#                    # first pop ai
#                    #for key in idx2:
#                    #    if tuple(key) in idx_lookup: idx_lookup.pop(tuple(key))
#                    ## then pop bj's
#                    #while True:
#                    #    flag = True
#                    #    # iterate through all items
#                    #    for key, value in idx_lookup.items():
#                    #        # find any bj's that should be removed
#                    #        w = np.where((v in idx_lookup.keys() for v in value))[0]
#                    #        # if all bj's for a give ai is removed, pop the ai key
#                    #        if w.size == len(value):
#                    #            idx_lookup.pop(key)
#                    #            flag = False
#                    #        # else just pop the individual bj
#                    #        elif w.size > 0:
#                    #            for i in w[::-1]:
#                    #                value.pop(i)
#                    #    # if no keys were removed, break
#                    #    if flag:
#                    #        break
#                    for pair in idx2:
#                        print pair
#                        mask[pair[0],pair[1]] = True
#                        M[pair[0],pair[1]] = 0
#                    M[best_pair[0],best_pair[1]] = 1
#                    #print molP.bonds[1]
#                    #print molQ.bonds[1]
#                    #quit()
#
#
#        # begin B: Relaxation
#        for itB in range(IB):
#            #Mold = M.copy()
#            #if simple_convergence:
#            #    if np.sum(abs(Mold-M)) < epsB:
#            #        break
#            #elif itB > 0:
#            if itB > 0:
#                #delta = 2 * np.sqrt(epsD*(con + np.log((N-1+epsD)/(1-epsD))/b)/(eps_gamma*N))
#                #print delta, np.sum(M**2)/N
#                #assert(delta < 0.05)
#                delta = 1e-2
#                if np.sqrt(np.sum((Mlast-M)**2)/N) < delta:
#                    break
#
#
#            for (a,i),v in idx_lookup.items():
#                Q[a,i] = sum((M[b,j] for b,j in v))
#            Q += np.eye(N) * M * gamma
#
#            # begin softassign
#            # more than ~exp(100) overflows, so move the range.
#            # VERY small values are unimportant anyways
#            #if offset:
#            #np.exp(b*Q-np.max(b*Q)+100, out=M)
#            #else:
#            np.exp(b*Q, out=M)
#            # apply mask (enforce atomtype equality and already determined assignments
#            #if use_mask:
#            M[mask] = 0
#            #if flag:
#            #    print M.sum(0)
#            #    print M.sum(1)
#            #    quit()
#            #if normalize:
#            #    M *= N/(M[:-1].sum())
#
#            # begin C
#            for itD in range(ID):
#                #Mold = M.copy()
#                # normalize across rows (except slack)
#                # TODO slack
#                M /= (np.sum(M,axis=1))[:,None]
#                # normalize across columns (except slack)
#                M /= (np.sum(M,axis=0))
#                # break if converged
#                #if simple_convergence:
#                #    if np.sum(abs(Mold-M)) < epsC:
#                #        break
#                #else:
#                if (abs(np.sum(M, axis = 1) - 1) < epsD).all():
#                    break
#                if itD == ID-1:
#                    print np.max(abs(np.sum(M, axis = 1) - 1)), epsD
#
#            # end C
#            # end Softassign
#
#            #if cut:
#            #    # weights is an attempt to make eps3 an averaged limit
#            #    weights = np.zeros(M.shape)
#            #    for a in range(M.shape[0]):
#            #        for i in range(M.shape[1]):
#            #            weights[a,i] = max(np.sum(~mask[a]), np.sum(~mask[:,i]))
#            #    weps = cut_eps / weights
#            #    # NOTE: do this?
#            #    w = np.where((M < weps) != mask)
#            #    if w[0].size > 0:
#            #        # update the idx_lookup to remove these indices
#            #        # since they are set to zero anyways
#            #        # first pop ai
#            #        for key in zip(*w):
#            #            idx_lookup.pop(key)
#            #        # then pop bj's
#            #        while True:
#            #            flag = True
#            #            # iterate through all items
#            #            for key, value in idx_lookup.items():
#            #                # find any bj's that should be removed
#            #                w = np.where((v in idx_lookup.keys() for v in value))[0]
#            #                # if all bj's for a give ai is removed, pop the ai key
#            #                if w.size == len(value):
#            #                    idx_lookup.pop(key)
#            #                    flag = False
#            #                # else just pop the individual bj
#            #                elif w.size > 0:
#            #                    for i in w[::-1]:
#            #                        value.pop(i)
#            #            # if no keys were removed, break
#            #            if flag:
#            #                break
#
#            #        mask[np.where(M < weps)] = True
#
#            E1 = 0
#            E2 = 0
#            E4 = 0
#            for (a,i),v in idx_lookup.items():
#                E1 += -0.5 * M[a,i]*(sum((M[b,j] for b,j in v)))
#                E2 += -0.5 * M[a,i]*(sum((M[b,j] for b,j in v)))
#                E4 += -0.5 * M[a,i]*(sum((M[b,j] for b,j in v)))
#                #E4 += M[a,i] * np.log(M[a,i]) * b**-1
#
#            E2 += -0.5*np.trace(M) * gamma
#            E4 += -0.5*np.trace(M) * gamma
#
#            if itB > 0:
#                dE1 = E1old - E1
#                dE2 = E2old - E2
#                dE4 = E4old - E4
#                #print dE1, dE2, dE4, M
#                #assert(dE1 > 0)
#                #assert(dE2 > 0)
#                #assert(dE4 > 0)
#                print dE1, dE2, dE4
#            E1old = E1
#            E2old = E2
#            E4old = E4
#            np.copyto(Mlast,M)
#
#        # end B
#        #if alt_annealing:
#        #    b += br
#        #else:
#        b *= br
#    # end A
#        #print b, m.diagonal()
#
#
#
## get normal vector to the plane spanned by r1,r2,r3,
## then determine if r4 is above or below the plane
#def get_cross_dot(r1,r2,r3,r4):
#    v1 = r2 - r1
#    v2 = r3 - r1
#    v3 = r4 - r1
#    n = np.cross(v1,v2)
#    return (n.dot(v3) > 0)
#
## get bond angle with r1 being the central atom
## http://stackoverflow.com/a/13849249
#def get_angle(r1,r2,r3):
#    v1 = r2 - r1
#    v2 = r3 - r1
#    v1 /= np.linalg.norm(v1)
#    v2 /= np.linalg.norm(v2)
#    return np.arccos(np.clip(np.dot(v1,v2),-1.0,1.0))*180/np.pi
#
#
## matrix involved in quaternion rotation
#def makeW(r1,r2,r3,r4=0):
#    W = np.asarray([
#             [r4, r3, -r2, r1],
#             [-r3, r4, r1, r2],
#             [r2, -r1, r4, r3],
#             [-r1, -r2, -r3, r4] ])
#    return W
#
## matrix involved in quaternion rotation
#def makeP(r1,r2,r3,r4=0):
#    P = np.asarray([
#             [r4, -r3, r2, r1],
#             [r3, r4, -r1, r2],
#             [-r2, r1, r4, r3],
#             [-r1, -r2, -r3, r4] ])
#    return P
#
## generate matrices needed in the calculation
#def quaternion_rotate(X, Y):
#    N = X.shape[0]
#    W = np.asarray([makeW(*Y[k]) for k in range(N)])
#    P = np.asarray([makeP(*X[k]) for k in range(N)])
#    Pt_dot_W = np.asarray([np.dot(P[k].T,W[k]) for k in range(N)])
#    W_minus_P = np.asarray([W[k] - P[k] for k in range(N)])
#    C1 = -np.sum(Pt_dot_W,axis=0)
#    C2 = 0.5*N
#    C3 = np.sum(W_minus_P,axis=0)
#    A = np.dot(C3.T,C3)*C2-C1
#    eigen = np.linalg.eigh(A)
#    r = eigen[1][:,eigen[0].argmax()]
#    rot = quaternion_transform(r)
#    return np.dot(X,rot)
#
#
#
## from http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
#def rand_rotation_matrix(deflection=1.0, randnums=None):
#    """
#    Creates a random rotation matrix.
#
#    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
#    rotation. Small deflection => small perturbation.
#    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
#    """
#    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
#
#    if randnums is None:
#        randnums = np.random.uniform(size=(3,))
#
#    theta, phi, z = randnums
#
#    theta = theta * 2.0*deflection*np.pi  # Rotation about the pole (Z).
#    phi = phi * 2.0*np.pi  # For direction of pole deflection.
#    z = z * 2.0*deflection  # For magnitude of pole deflection.
#
#    # Compute a vector V used for distributing points over the sphere
#    # via the reflection I - V Transpose(V).  This formulation of V
#    # will guarantee that if x[1] and x[2] are uniformly distributed,
#    # the reflected points will be uniform on the sphere.  Note that V
#    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
#
#    r = np.sqrt(z)
#    Vx, Vy, Vz = V = (
#        np.sin(phi) * r,
#        np.cos(phi) * r,
#        np.sqrt(2.0 - z)
#        )
#
#    st = np.sin(theta)
#    ct = np.cos(theta)
#
#    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
#
#    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
#
#    M = (np.outer(V, V) - np.eye(3)).dot(R)
#    return M
#
#
#def closest_nodes(nodes):
#    nodes = np.asarray(nodes)
#    deltas = nodes - node
#    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
#    return np.argmin(dist_2)
