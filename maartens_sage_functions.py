from sage.all import *

def counts(list):
    s=set(list)
    return [(i,list.count(i)) for i in sorted(s)]

def tate_normal_form(E,p):
    assert p!=E([0,1,0])
    #assert p*2!=E([0,1,0])
    #assert p*3!=E([0,1,0])
    E=E.change_weierstrass_model(1,p[0],0,p[1])
    E=E.change_weierstrass_model(1,0,E.a4()/E.a3(),0)
    u=E.a3()/E.a2()
    a1,a2,a3,a4,a6=E.ainvs()
    ainvs=[a1/u,a2/u**2,a3/u**3,a4/u**4,a6/u**6]
    return ainvs



def ambient_integral_structure_matrix(M,verbose=False):
    """
    In certatain cases (weight two any level) this might return the integral structure of a
    an ambient modular symbol space. I wrote this because for high level this is very slow
    in sage because there is no good sparse hermite normal form code in sage for huge matrices.
    """
    if verbose: tm = cputime(); mem = get_memory_usage(); print "Int struct start"
    #This code is the same as the firs part of M.integral_structure
    G = set([i for i, _ in M._mod2term])
    G = list(G)
    G.sort()
    #if there is a two term relation between two manin symbols we only need one of the two
    #so that's why we only use elements from G instead of all manin symbols.
    if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "G"
    B = M._manin_gens_to_basis.matrix_from_rows(list(G)).sparse_matrix()
    if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "B"
    #The collums of B now span M.integral_structure as ZZ-module
    B, d = B._clear_denom()
    if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "Clear denom"
    if d == 1:
        #for explanation see d == 2
        assert len(set([B.nonzero_positions_in_row(i)[0] for i in xrange(B.nrows()) if len(B.nonzero_positions_in_row(i)) == 1 and B[i, B.nonzero_positions_in_row(i)[0]] == 1])) == B.ncols(), "B doesn't contain the Identity"
        if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "Check Id"
        ZZbasis = MatrixSpace(QQ, B.ncols(), sparse=True)(1)
        if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "ZZbasis"
    elif d == 2:
        #in this case the matrix B will contain 2*Id as a minor this allows us to compute the hermite normal form of B in a very efficient way. This will give us the integral basis ZZbasis.
        #if it turns out to be nessecarry this can be generalized to the case d%4==2 if we don't mind to only get the right structure localized at 2
        assert len(set([B.nonzero_positions_in_row(i)[0] for i in xrange(B.nrows()) if len(B.nonzero_positions_in_row(i)) == 1 and B[i, B.nonzero_positions_in_row(i)[0]] == 2])) == B.ncols(), "B doesn't contain 2*Identity"    
        if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "Check 2*Id"
        E = matrix_modp(B,sparse=True)
        if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "matmodp"
        E = E.echelon_form()
        if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "echelon"
        ZZbasis = MatrixSpace(QQ, B.ncols(), sparse=True)(1)
        for (pivot_row, pivot_col) in zip(E.pivot_rows(), E.pivots()):
            for j in E.nonzero_positions_in_row(pivot_row):
                ZZbasis[pivot_col, j] = QQ(1) / 2 
        if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "ZZbasis"
    else:
        raise NotImplementedError
    return ZZbasis



def cuspidal_integral_structure_matrix(M,verbose=False):
    """
    Computes the cuspidal integral structure matrix of a modular symbols space in a runningtime hopefully faster then that of sage
    """
    #now we compute the integral kernel of the boundary map with respect to the integral basis. This will give us the integral cuspidal submodule.
    ZZbasis = ambient_integral_structure_matrix(M,verbose=verbose)              
    boundary_matrix = M.boundary_map().matrix()
    if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "Boundary matrix"
    ZZboundary_matrix=(ZZbasis*boundary_matrix).change_ring(ZZ)
    if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "ZZBoundary matrix"
    left_kernel_matrix=ZZboundary_matrix.transpose().dense_matrix()._right_kernel_matrix(algorithm='pari')
    if type(left_kernel_matrix)==tuple:
        left_kernel_matrix=left_kernel_matrix[1]
    if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "kernel matrix"
    ZZcuspidal_basis=left_kernel_matrix*ZZbasis
    if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "ZZkernel matrix"
    S=M.cuspidal_subspace()
    assert ZZcuspidal_basis.change_ring(QQ).echelon_form()==S.basis_matrix() , "the calculated integral basis does not span the right QQ vector space" # a little sanity check. This shows that the colums of ZZcuspidal_basis really span the right QQ vectorspace
    if verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "finnished"
    return ZZcuspidal_basis



def galois_action(self, t, N):
    from sage.modular.modsym.p1list import lift_to_sl2z
    if self.is_infinity():
        return self
    if not isinstance(t, Integer): 
        t = Integer(t)

    a = self._Cusp__a
    b = self._Cusp__b * t.inverse_mod(N)
    if b.gcd(a) != ZZ(1):
        _,_,a,b = lift_to_sl2z(a,b,N)
        a = Integer(a); b = Integer(b)

    # Now that we've computed the Galois action, we efficiently
    # construct the corresponding cusp as a Cusp object.
    return Cusp(a,b,check=False)

def galois_orbit(cusp,G):
    N=G.level()
    orbit=set([])
    for i in xrange(1,N):
        if gcd(i,N)==1:
            orbit.add(G.reduce_cusp(galois_action(cusp,i,N)))
    return tuple(sorted(orbit))
    
def galois_orbits(G):
    return set([galois_orbit(c,G) for c in G.cusps()])


def cuspidal_integral_to_rational_basis(M):
    """
    Returns the base change matrix for the cuspidal subspace of M
    """
    S=M.cuspidal_subspace()
    return Matrix([S.coordinate_vector(i) for i in cuspidal_integral_structure_matrix(M)])

def cuspidal_rational_to_integral_basis(M):
    return cuspidal_integral_to_rational_basis(M)**(-1)

def period_mapping(M):
    N=M.level()
    S=M.cuspidal_subspace()
    q=exists(primes(N+1,N**3+3), lambda x: x%N==1)[1] 
    #q is the smallest prime equal to 1 mod N
    Tq=M.hecke_matrix(q)
    dq=M.diamond_bracket_operator(q).matrix()
    to_cusp_space=(Tq-dq*q-1) #maybe use smaller primes in the future*(Tq-dq-q) 
    #to_cusp_space is a map whose image (over QQ) is the cuspidal subspace and is invertible when 
    #restricted to the cuspidal subspace so this means that we can compute the intergral period 
    #mapping (projection onto the cusps) by composing to_cusp_space by its inverse on the cuspidal
    #subspace
    codomain_restricted = to_cusp_space.restrict_codomain(S)
    restricted = codomain_restricted.restrict_domain(S)
    return codomain_restricted*restricted**(-1)

@cached_function
def integral_period_mapping(M):
    return period_mapping(M)*cuspidal_rational_to_integral_basis(M)   
    
def modular_unit_lattice(G,return_ambient=True):
    M=G.modular_symbols()
    S=M.cuspidal_subspace()
    period_mapping=integral_period_mapping(M)
    H1QQ=QQ**period_mapping.ncols()
    cusps=G.cusps()
    assert cusps[-1]==Cusp(infinity)
    n=len(cusps)
    D=ZZ**n
    D0=D.submodule_with_basis([D.gen(i)-D.gen(n-1) for i in xrange(n-1)])
    period_images=[M.coordinate_vector(M([c,infinity]))*period_mapping for c in cusps if not c.is_infinity()]
    m=Matrix(period_images).transpose().augment(H1QQ.basis_matrix()).transpose()
    mint=(m*m.denominator()).change_ring(ZZ)
    kernel=mint.kernel()
    F0=kernel.basis_matrix().change_ring(ZZ).transpose()[:n-1].transpose()
    F0inD0=D0.submodule([D0.linear_combination_of_basis(i) for i in F0])
    if return_ambient:
        return F0inD0,D0
    return F0inD0
    
def rational_modular_unit_lattice(G,return_ambient=True):
    cusps=G.cusps()
    n=len(cusps)
    D=ZZ**n
    F0inD0,D0=modular_unit_lattice(G,return_ambient=True)
    Drational=D.submodule([sum([D.gen(cusps.index(i)) for i in j]) for j in galois_orbits(G)])
    F0rational=F0inD0.intersection(Drational)
    if return_ambient:
        return F0rational,D0.intersection(Drational)
    return F0rational

  
def rational_cuspidal_classgroup(G):
    functions,divisors=rational_modular_unit_lattice(G,return_ambient=True)
    return divisors.quotient(functions)

