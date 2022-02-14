
from sage.all import *
from itertools import product,combinations

def border_apolarity_cycl_inv(T,reps,r,verbose=True):
    # Check that the data of reps actually stabilize $T$ and compute the simple roots
    simple_roots = check_T_is_stabilized(T,reps)
    # Find the representation data for $T(C^*)^\perp$ and report the matrix of
    # the embedding into $A^*\ot B^*$
    reptp,em = Tperp_rep_and_embedding(T,reps) 
    # Find the weight diagram of $T(C^*)^\perp$
    wtd = weight_diagram(reptp,simple_roots)

    a = len(T)
    b = T[0].nrows()
    cand110 = []

    # Enumerate Borel fixed subspace families of $T(C^*)^\perp$. These are given
    # by matrices Stperp over a polynomial ring and a corresponding ideal I the
    # coefficients of Stperp are considered modulo, i.e., the family takes
    # parameters in $V(I)$.
    for Stperp,I in borel_fixed_subspaces(wtd, em.nrows()-r, verbose):
        # Embed into $A^*\ot B^*$
        Sab = em*Stperp

        # Restrict to the closed subset of $V(I)$ passing the $(210)$ test
        S210 = matrix_11_to_21(Sab,a)
        I = minors_ideal(S210, S210.nrows() - r + 1, I)

        # Restrict to the closed subset of $V(I)$ passing the $(120)$ test
        S120 = matrix_11_to_21(transpose_tensor(Sab,a),b)
        I = minors_ideal(S120, S120.nrows() - r + 1, I)

        # Check if the family is nonempty
        if QQ(1) not in I:
            # If there is a single solution in the family add it to cand110,
            # otherwise raise an error. This is sufficient for $M_{\langle 3\rangle}$ and $\operatorname{det}_3$, 
            # and more generally one can go on to apply the $(111)$ test to
            # parameterized families of candidate triples
            cand110.append(Sab.apply_map(lambda e: QQ(I.reduce(e)), QQ))

    cand111 = []
    for xs in product(*map(enumerate,[cand110]*3)):
        ixs = tuple(i for i,x in xs)
        if ixs != min([ixs,ixs[1:]+ixs[:1],ixs[2:]+ixs[:2]]):
            # Skip triples equal to others we check modulo cyclic permutation
            continue
        if verbose:
            print(ixs,end=' ')
            sys.stdout.flush()

        W = matrix_to_111(*[x for _,x in xs])
        if W.rank() <= W.nrows() - r:
            cand111.append(W)
    return cand111,cand110

def check_T_is_stabilized(T,reps):
    a = len(T)
    b,c = T[0].dimensions()

    # There is one weight for each basis vector
    assert all(len(wts) == d for d,(wts,_) in zip((a,b,c),reps))
    # The given representations of $\mathfrak{n}$ operate on the correct dimensional space
    assert all(x.nrows() == d and x.ncols() == d 
            for d,(_,xs) in zip((a,b,c),reps) for x in xs)

    # compute roots of the simple roots vectors xs
    simple_roots = []
    for xi in range(len(reps[0][1])):
        for wts,xs in reps:
            if not xs[xi].is_zero():
                i,j = xs[xi].nonzero_positions()[0]
                simple_roots.append(wts[i] - wts[j])
                break

    # Check xs are actally root vectors for the corresponding roots and that 
    # the simple roots agree for each rep in reps
    for wts,xs in reps:
        assert len(simple_roots) == len(xs)
        for root,x in zip(simple_roots,xs):
            for i,j in x.nonzero_positions():
                assert wts[i]-wts[j] == root

    # Check T is weight zero in $A\ot B\ot C$
    for i,m in enumerate(T):
        for j,k in m.nonzero_positions():
            assert (reps[0][0][i] + reps[1][0][j] + reps[2][0][k]).is_zero()
    
    # Check T is closed under $\mathfrak{n}$
    (_,xsA), (_,xsB), (_,xsC) = reps
    for xa,xb,xc in zip(xsA, xsB, xsC):
        xaT = [sum(xa[i,j] * T[j] for j in range(a)) for i in range(a)]
        xbT = [xb*m for m in T]
        xcT = [m*xc.T for m in T]
        # we should have xaT + xbT + xcT = 0
        assert all((m1+m2+m3).is_zero() for m1,m2,m3 in zip(xaT, xbT, xcT))
    
    return simple_roots

def borel_fixed_subspaces(wtd,D,verbose=False):
    for i,dlambda in enumerate(possible_dim_assignments(wtd,D)):
        if verbose:
            print(i,end=' ')
        for S in borel_fixed_subspaces_dlambda(wtd,dlambda,verbose):
            yield S

def borel_fixed_subspaces_dlambda(wtd,dlambda,verbose=False):

    dlist = [(wt,d,len(wtd.get_vertex(wt))) for wt, d in dlambda.items()]

    if verbose: 
        missing = [wt for wt,d in dlambda.items() if d < len(wtd.get_vertex(wt))]
        missing.sort(key=lambda wt: (sum(wt),wt))
        print('missing',missing)


    if all(m == f or m == 0 for wt,m,f in dlist):
        # $S$ is full in full in every weight space for which it is nonzero.
        # Such a set satisfies the weight diagram arrow conditions due to the
        # necessary conditions on $d_{\lambda}$. No need to check again
        yield (identity_matrix(QQ, sum(f for _,_,f in dlist))[:,
                [i for wt,m,f in dlist if m == f for i in wtd.get_vertex(wt)]],QQ.ideal(0))
        return

    # Otherwise, we need parameters. Here we check all combinations of choices
    # of pivot columns in each grassmannian

    for nzsi,nzs in enumerate(product(*[combinations(range(f),m) for wt,m,f in dlist])):

        if verbose:
            print(nzsi,end=' ')
            sys.stdout.flush()

        nvars = sum( (nzi+1)*(j-i-1) for (wt,m,f),nz in zip(dlist,nzs)
            for nzi,(i,j) in enumerate(zip(nz,nz[1:]+(f,))))
        R = PolynomialRing(QQ,'t',nvars,implementation='singular') if nvars > 0 else QQ

        S = {}
        xi = 0
        for (wt,m,f),nz in zip(dlist,nzs):
            t = matrix(R,f,m,sparse=True)
            t[nz,:] = identity_matrix(R,m,sparse=True)
            for j,ks in enumerate(zip(nz,nz[1:]+(f,))):
                inc = (ks[1]-ks[0]-1)*(j+1)
                t[ks[0]+1:ks[1],:j+1] = matrix(R,ks[1]-ks[0]-1,j+1,R.gens()[xi:xi+inc])
                xi += inc 

            S[wt] = (t, nz) # remember the pivots for equations below

        # Now restrict the parameters appearing so that $S$ is closed under the
        # arrows of the weight diagram
        eqs = []
        for wta,wtb,x in wtd.edges():
            S1,_ = S[wta]
            S2,S2pvts = S[wtb]
            # Reduce the columns of xS1 modulo S2 using the pivots of S2.
            xS1 = (x*S1).T
            S2 = S2.T
            for i,j in enumerate(S2pvts):
                for k in xS1.nonzero_positions_in_column(j):
                    xS1[k] -= xS1[k,j]*S2[i]
            # xS1 now is in a normal form modulo S2, so the equations of
            # containment are the conditions that xS1 == 0
            eqs.extend(xS1.dict().values())

        I = R.ideal(eqs)
        if R.one() in I:
            continue

        Smat = block_diagonal_matrix([S[wt][0] for wt,_,_ in dlist],sparse=True)

        ix = [i for wt,_,_ in dlist for i in wtd.get_vertex(wt)]
        ixi = [None]*len(ix)
        for j,i in enumerate(ix):
            ixi[i] = j
        Smat = Smat[ixi,:]

        yield (Smat,I)

    if verbose: 
        print()

def possible_dim_assignments(wtd,dim):
    lp = MixedIntegerLinearProgram()
    for wt in wtd:
        lp.set_min(lp[wt],0)
        lp.set_max(lp[wt],len(wtd.get_vertex(wt)))

    lp.add_constraint(lp.sum(lp[wt] for wt in wtd) == dim)

    for wt in wtd:
        # Condition $\eqref{primal}$
        for k in range(1,len(wtd.outgoing_edges(wt))+1):
            for es in combinations(wtd.outgoing_edges(wt),k):
                kdim = block_matrix([[xr] for _,_,xr in es]).right_kernel().dimension()
                lp.add_constraint(lp[wt] <= kdim + lp.sum(lp[wtr] for _,wtr,_ in es))

        # Condition $\eqref{dual}$
        # This when $k = 1$ occurs also for condition $\eqref{primal}$ as well, so we skip it.
        for k in range(2,len(wtd.incoming_edges(wt))+1):
            for es in combinations(wtd.incoming_edges(wt),k):
                cokdim = block_matrix([[-xr.transpose()] 
                    for _,_,xr in es]).right_kernel().dimension()
                lp.add_constraint(len(wtd.get_vertex(wt)) - lp[wt] <= cokdim + lp.sum(
                    len(wtd.get_vertex(wtl)) - lp[wtl] for wtl,_,xr in es))

    # Condition $\eqref{composition}$
    for wtstart in wtd:
        for _,wtmid,x1 in wtd.outgoing_edges(wtstart):
            for _,wtlast,x2 in wtd.outgoing_edges(wtmid):
                kdim = (x2*x1).right_kernel().dimension()
                lp.add_constraint(lp[wtstart] <= kdim + lp[wtlast])

    # We have computed the polytope, now enumerate the integer points

    # Fix an ordering of the weights $\lambda_i$
    wts = wtd.topological_sort()

    # Given assignments of $d_{\lambda_j}$, $j\le i-1$, set $d_{\lambda_i}$ to the values for which the 
    # polytope still intersects the coordinate hyperplane corresponding to the
    # current partial assignment and recursively apply this procedure.
    from sage.numerical.mip import MIPSolverException
    def dfs(i):
        if i == len(wts):
            # We have set all coordinates to integer values and remain in the
            # polytope, report this as a solution
            yield {wt : int(lp.get_min(lp[wt])) for wt in wts}
            return
        wt = wts[i]
        mult = int(lp.get_max(lp[wt]))
        for dcur in range(0, mult+1):
            lp.set_min(lp[wt],dcur)
            lp.set_max(lp[wt],dcur)
            try:
                # The polytope has nonzero intersection with set corresponding
                # to fixing the first $i$ coordinates as we have. Recursively try 
                # to set the remaining coordinates and report all solutions.
                lp.solve()
                for pt in dfs(i+1):
                    yield pt
            except MIPSolverException:
                # There are no points in the polytope with the first $i$
                # coordinates fixed to the values we have used. Don't search
                # further.
                pass
        lp.set_min(lp[wt],0)
        lp.set_max(lp[wt],mult)

    return dfs(0)

def matrix_11_to_21(B,a):
    b = B.nrows() // a
    S = B.ncols()
    W = {}
    # i,k < a
    # j < b
    # s < S
    for I,s in B.nonzero_positions():
        i,j = I // b, I % b
        v = B[I,s]
        for k in range(a):
            mi,ma = min(i,k),max(i,k)
            ix = (( binomial(ma+1,2)+mi )*b+j,s*a+k)
            W[ix] = W.get(ix,0) + v
    W = matrix(B.base_ring(),binomial(a+1,2)*b,S*a,W)
    return W

def matrix_to_111(A,B,C):
    a = int(sqrt(B.nrows()*C.nrows()/A.nrows()))
    b = C.nrows() // a
    c = B.nrows() // a
    W = {}
    for i in range(a):
        for I,l in A.nonzero_positions():
            j,k = I // c, I % c
            W[((i*b+j)*c+k, i*A.ncols()+l)] = A[j*c+k,l]
    for j in range(b):
        for I,l in B.nonzero_positions():
            k,i = I // a, I % a
            W[((i*b+j)*c+k, a*A.ncols() + j*B.ncols()+l)] = B[k*a+i,l]
    for k in range(c):
        for I,l in C.nonzero_positions():
            i,j = I // b, I % b
            W[((i*b+j)*c+k, a*A.ncols() + b*B.ncols() +\
                    k*C.ncols()+l)] = C[i*b+j,l]
    W = matrix(A.base_ring(),a*b*c,a*A.ncols()+b*B.ncols()+c*C.ncols(),W)
    return W

def transpose_tensor(B,a):
    b = B.nrows() // a
    S = B.ncols()
    Bp = {}
    for I,k in B.nonzero_positions():
        i,j = I // b, I % b
        Bp[(j*a+i,k)] = B[I,k]
    return matrix(B.base_ring(),B.nrows(),B.ncols(),Bp)

def Tperp_rep_and_embedding(T,reps,missing=2):
    reps = reps[missing+1:] + reps[:missing]
    repAB = module_product(*[module_dual(rep) for rep in reps])

    Tcycl = T
    for i in range(missing):
        Tcycl = tensor_cycl(Tcycl)
    Tflattening = matrix(QQ,[m.list() for m in Tcycl],sparse=True)
    em = Tflattening.right_kernel_matrix().transpose().sparse_matrix()

    reptp = submodule_from_basis(repAB,em)
    return reptp,em

def module_product(repa,repb):
    wtsa, xsa = repa
    wtsb, xsb = repb
    wtsab = [wta + wtb for wta,wtb in product(wtsa,wtsb)]
    xsab = []
    for xa,xb in zip(xsa,xsb):
        xab = xa.tensor_product(identity_matrix(QQ,xb.nrows()))
        xab += identity_matrix(QQ,xa.nrows()).tensor_product(xb)
        xsab.append(xab)
    return (wtsab,xsab)

def module_dual(rep):
    wts, xs = rep
    wtsd = [-wt for wt in wts]
    xsd = [-x.transpose() for x in xs]
    return (wtsd,xsd)

def submodule_from_basis(rep,em):
    wts, xs = rep

    # check em describes a basis
    assert em.rank() == em.ncols()

    wtsw = []
    # check the columns of em are weight vectors and find their weights
    for v in em.columns():
        wt = wts[v.nonzero_positions()[0]]
        assert all(wt == wts[j] for j in v.nonzero_positions())
        wtsw.append(wt)

    xsw = []
    # check $W$ is fixed under the simple root vectors and compute the matrix of
    # their action
    for x in xs:
        # find y solving x*em == em*y, and raises an error if there is none
        y = em.solve_right(x*em) 
        # Since em is basis, the eqution x*em == em*y is sufficient to
        # guarantee the claims
        xsw.append(y)

    return (wtsw,xsw)

def weight_diagram(rep,simple_roots):
    # tot = simultaneous_eigenspace(a[2])
    wts, xs = rep
    dim = len(wts)

    assert all(x.nrows() == dim and x.ncols() == dim for x in xs)

    wtd = DiGraph()
    for bi,wt in enumerate(wts):
        wt.set_immutable()
        if wt not in wtd:
            wtd.add_vertex(wt)
            wtd.set_vertex(wt,[])
        wtd.get_vertex(wt).append(bi)

    for wt in wtd:
        wt.set_immutable()
        for root,x in zip(simple_roots,xs):
            wt_raised = wt + root
            wt_raised.set_immutable()
            if wt_raised in wtd:
                wtd.add_edge(wt,wt_raised,
                        x[wtd.get_vertex(wt_raised),wtd.get_vertex(wt)])
    
    return wtd

def minors_ideal(m,r,I=None):
    from collections import Counter

    S = m.base_ring()
    if I is None: I = S.ideal()

    if S is QQ:
        return QQ.ideal(1) if QQ(1) in I or m.rank() >= r else QQ.ideal(0)

    R = PolynomialRing(S.base_ring(),S.gens()+('d',))
    dv = R.gens()[-1]

    def rec(m,I,d,r):
        if m.is_zero():
            return I.elimination_ideal(dv).change_ring(S)
        if r == 1:
            return (I + m.coefficients()).elimination_ideal(dv).change_ring(S)

        # select the most common nonzero element of 
        f = Counter(m.dict().values()).most_common(1)[0][0]

        # Case 1: pass to the closed set where f == 0
        Icur = I + f
        mcur,lo = elimination_by_units(m,Icur)
        I1 = rec(mcur[lo:,lo:],Icur,d,r-lo) if lo < r else S.ideal(1)

        # Case 2: pass to the open set where f not identically zero, e.g., localize at f
        mcur = m.subs({dv : dv*f})
        Icur = I.elimination_ideal(dv) + [dv*d*f-1]
        mcur,lo = elimination_by_units(mcur,Icur)
        I2 = rec(mcur[lo:,lo:],Icur,d*f,r-lo) if lo < r else S.ideal(1)

        # Combine the result of both cases: $I_1\cap I_2$ corresponds to the
        # union of the corresponding parameter sets
        return I1.intersection(I2)

    return rec(m.change_ring(R), I.change_ring(R), R.one(), r)

def elimination_by_units(m,I):
    m = m.apply_map(lambda e: I.reduce(e))

    from sage.libs.singular.function import singular_function
    lift = singular_function('lift')

    # computes the inverse mod $I$ if it exists, otherwise None
    def try_inverse(e):
        try:
            return lift(I+e,1).list()[-1]
        except RuntimeError:
            return None

    r = 0
    while True:
        if r == min(m.nrows(),m.ncols()):
            break

        einv = None
        for (i,j),e in m[r:,r:].dict().items():
            einv = try_inverse(e)
            if einv is not None:
                break
        if einv is None:
            break
        i += r
        j += r

        m.swap_rows(r,i)
        m.swap_columns(r,j)
        m[r,:] *= einv
        for k in m.nonzero_positions_in_row(r):
            m[r,k] = I.reduce(m[r,k])
        assert m[r,r] == 1
        for i in m.column(r)[r+1:].nonzero_positions():
            i += r+1
            m.add_multiple_of_row(i,r,-m[i,r])
            for k in m.nonzero_positions_in_row(i):
                m[i,k] = I.reduce(m[i,k])
        r += 1

    return m,r

def tensor_cycl(T):
    m = len(T)
    n,s = T[0].dimensions()
    S = [zero_matrix(QQ,s,m) for i in range(n)]
    for i,j,k in product(range(m),range(n),range(s)):
        S[j][k,i] = T[i][j,k]
    return S
