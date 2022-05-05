"""
First import everything::

    sage: from mdsage import *
    sage: counts("113")
    [('1', 2), ('3', 1)]
"""

import sage
from sage.all import (AbelianGroup,
                      cached_function,
                      ceil,
                      cputime,
                      Cusp,
                      exists,
                      Gamma0,
                      GammaH,
                      gcd,
                      get_memory_usage,
                      infinity,
                      Integer,
                      Integers,
                      Matrix,
                      MatrixSpace,
                      ModularSymbols,
                      pari,
                      primes,
                      prod,
                      QuadraticForm,
                      QQ,
                      vector,
                      zero_matrix,
                      ZZ)

from .kamiennys_criterion import matrix_modp
from copy import copy

def gonality_lower_bound(G,lambda1 = 0.238):
    """
    Return the gonality bound of the modular curve X(G)
    as given by Dan Abramovich in "A linear lower bound on the
    gonality of modular curves"
    
    lambda1 is conjectured to 1/4th, but has only been proven to
    be > 0.238 [1, p. 3], bigger lambda1 means a better lowerbound.
    
    EXAMPLES::

        sage: from mdsage import *
        sage: gonality_lower_bound(Gamma1(171))
        129

    [1] On the torsion of elliptic curves over quartic number fields - Daeyol Jeon, Chang Heon Kim and Euising Park
    """
    return ceil(G.projective_index()*lambda1/24)


def counts(list):
    """
    On input of some iterable l, this function returns a list of pairs (s,c) where s runs trough
    the distinct elements of l and c indicates how often each element occurs.
    
    EXAMPLES::

        sage: from mdsage import *
        sage: l = "A string is also iterable!"
        sage: counts(l)
        [(' ', 4), ('!', 1), ('A', 1), ('a', 2), ('b', 1), ('e', 2), ('g', 1), ('i', 3), ('l', 2), ('n', 1), ('o', 1), ('r', 2), ('s', 3), ('t', 2)]
    """
    s=set(list)
    return [(i,list.count(i)) for i in sorted(s)]

def positive_part(v):
    """
    Returns the vector consisting only of the postive entries of the given vector
    
    EXAMPLES::

        sage: from mdsage import *
        sage: v = vector([1,-1,2,-2,3,-3])
        sage: positive_part(v)
        (1, 0, 2, 0, 3, 0)
    """
    return v.parent().ambient_module()([i if i > 0 else 0 for i in v])

def tate_normal_form(E,p):
    """
    Writes a pair (E,p) of elliptic curve with a point in tate normal form.
    
    On input a pair (E,p) where E is an elliptic curve and p a point of order > 3 on E, this function
    returns the a invariants a1,a2,a3,a4,a6 after writing the pair (E,p) in the tate normal form.
    I.e. the form in which p=(0,0), a4=a6=0 and a2=a3 .
   
    EXAMPLES:

    We show in this example how to calculate the diamond operators
    on X_1(5) which is isomorphic to P1:

        sage: from mdsage import *
        sage: R.<b> = QQ[]

    This is a model of the universal elliptic curve over X_1(5):

        sage: E_univ = EllipticCurve(R.fraction_field(),[b+1,b,b,0,0])
        sage: E_univ([0,0])*5
        (0 : 1 : 0)
        sage: tate_normal_form(E_univ,2*(E_univ(0,0)))
        [(-b + 1)/-b, -1/b, 1/-b, 0, 0]
        
    This shows that <2> works on the above model of X_1(5) by sending b to -1/b.
    """
    assert p!=E([0,1,0])
    #assert p*2!=E([0,1,0])
    #assert p*3!=E([0,1,0])
    E=E.change_weierstrass_model(1,p[0],0,p[1])
    E=E.change_weierstrass_model(1,0,E.a4()/E.a3(),0)
    u=E.a3()/E.a2()
    a1,a2,a3,a4,a6=E.ainvs()
    ainvs=[a1/u,a2/u**2,a3/u**3,0,0]
    return ainvs

def diamond_operator(E,d):
    """
    Returns the tate normal form of the pair (E,d*(0,0))
    
    INPUT:
    
    - E - an elliptic curve in tate normal form
    - d - an integer
        
    OUTPUT:
    
    - a1,a2,a3,a4,a6 - the a invatiants of E after writing  (E,d*(0,0)) in tate normal form.
        
    EXAMPLES::

        sage: from mdsage import *
        sage: E = EllipticCurve([11,10,10,0,0])
        sage: diamond_operator(E,2)
        [9/10, -1/10, -1/10, 0, 0]
    """
    return tate_normal_form(E,E([0,0])*d)

def diamond_orbit(E,N=None):
    """
    """
    if N==None:
        N=E([0,0]).order()
    for d in xrange(1,N):
        if gcd(d,N)==1:
            yield diamond_operator(E,d)

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
    
    Input:
    
    - M - a modular symbols space
    
    Output:
    
    - a matrix whose rows give a ZZ basis for the cuspidals supspace with respect to the standard basis of M
    
    Tests::

        sage: from mdsage import *
        sage: M = ModularSymbols(Gamma1(15))
        sage: cuspidal_integral_structure_matrix(M)
        [ 0  0  0  0  0  0  0  0  0  0  1  0  0  0  0 -1  0]
        [ 0  0  0  0  0  1  0  0  1  0  0 -1  0  0  0  0  0]
    
    """
    #now we compute the integral kernel of the boundary map with respect to the integral basis. This will give us the integral cuspidal submodule.
    if verbose: tm = cputime(); mem = get_memory_usage(); 
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
    
def modular_unit_lattice(G,return_ambient=True,ambient_degree_zero=True):
    M=G.modular_symbols()
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
    if return_ambient and ambient_degree_zero:
        return F0inD0,D0
    if return_ambient:
        return F0inD0,D
    return F0inD0
    
def rational_modular_unit_lattice(G,return_ambient=True,ambient_degree_zero=True):
    cusps=G.cusps()
    n=len(cusps)
    D=ZZ**n
    F0inD0,D0=modular_unit_lattice(G,return_ambient=True,ambient_degree_zero=ambient_degree_zero)
    Drational=D.submodule([sum([D.gen(cusps.index(i)) for i in j]) for j in galois_orbits(G)])
    F0rational=F0inD0.intersection(Drational)
    if return_ambient:
        return F0rational,D0.intersection(Drational)
    return F0rational

  
def rational_cuspidal_classgroup(G,degree_zero=True):
    """
    On input a congruence subgroup G this function returns the lattice of sums of galois invariant orbits of cusps
    modulo divisors of modular units.
    """
    functions,divisors=rational_modular_unit_lattice(G,return_ambient=True,ambient_degree_zero=degree_zero)
    return divisors.quotient(functions)
    
def generators_of_subgroups_of_unit_group(R):
    """
    INPUT:
    
    - R - a commutative ring whose unit group is finite
    
    OUTPUT:
    
    - An iterator which yields a set of generators for each subgroup of R^*
    
    EXAMPLES::

        sage: from mdsage import *
        sage: list(generators_of_subgroups_of_unit_group(Integers(28)))
        [[15, 17], [11], [3], [13, 15], [15], [27], [17], [9], [13], []]
    
    """
    gens = R.unit_gens()
    invariants = [g.multiplicative_order() for g in gens]
    assert all(i!=0 for i in invariants)
    A=AbelianGroup(invariants)
    for G in A.subgroups():
        yield [prod(f**e for f,e in zip(gens,g.exponents())) for g in G.gens()]

def congruence_groups_between_gamma0_and_gamma1(N):
    """
    INPUT:
        
    - N - an integer
    
    OUTPUT:
        
    - The set of all congruence subgroups contained in Gamma0(N) that contain Gamma1(N)
    
    EXAMPLES::

        sage: from mdsage import *
        sage: congruence_groups_between_gamma0_and_gamma1(1)
        {Modular Group SL(2,Z)}
    
        sage: congruence_groups_between_gamma0_and_gamma1(15)
        {Congruence Subgroup Gamma_H(15) with H generated by [14],
         Congruence Subgroup Gamma_H(15) with H generated by [4, 11, 14],
         Congruence Subgroup Gamma0(15)}
    """
    if N == 1:
        return set([Gamma0(1)])
    level_N_modular_groups = set([])
    for gens in generators_of_subgroups_of_unit_group(Integers(N)):
        gens = gens+[-1]
        H = sage.modular.arithgroup.congroup_gammaH._list_subgroup(N,gens)
        G = GammaH(N,H)
        level_N_modular_groups.add(G)
    return level_N_modular_groups

def count_points_J_H(G,p):
    """
    Returns the number of points on J_H(GF(p))
    
    INPUT:
        
    - G - a congruence subgroup of type GammaH
    - p - a prime
    
    OUTPUT:
        
    - the number of F_p points on the jacobian of the modular curve X_H
    
    EXAMPLES::

        sage: from mdsage import *
        sage: count_points_J_H(Gamma0(29),19)
        196
    """
    M=ModularSymbols(G,sign=1)
    S=M.cuspidal_subspace()
    dq=S.diamond_bracket_matrix(p)
    #print "computing tq"
    Tq=S.hecke_matrix(p)
    #assert (Tq-dq-p).det() == (Tq-p*dq-1).det()
    return (dq+p-Tq).det()


def Gamma11(m,n):
    """
    sage doesn't have  the congruence subgroup
    for the curve X_1(m,mn) but the
    group G defined here is conjugate to the group
    defining the modular curve X_1(m,mn)the one we are
    interested in so give isomorphic modular curves
    where the isomorphism even has a moduli interpretation
    In particular I use the surjection X_1(m^2n) -> X_1(m,mn)
    Which sends (E,P) to (E/<mnP>, Q mod <mnP>, P mod <mnP>) where Q in E
    is such that <Q,mnP> = zeta_m
    note that this surjection is defined over Q <=> m=2.
    
    If d is a unit in Z/m^2nZ then sometimes (E,P) and (E,dP) are
    mapped to the same point on X_1(m,mn)
    This happens precicely if d is +- 1 modulo mn if j(E) is not 0 or 1728
    """
    return GammaH(m**2*n,[-1]+[d*m*n+1 for d in range(1,m)])

      
def QuadraticForm_from_quadric(Q):
    """
    On input a homogeneous polynomial of degree 2 return the Quadratic Form corresponding to it.
    """
    R = Q.parent()
    assert all([sum(e)==2 for e in Q.exponents()])
    M = copy(zero_matrix(R.ngens(),R.ngens()))
    for i in range(R.ngens()):
        for j in range(R.ngens()):
            if i==j:
                M[i,j]=2*Q.coefficient(R.gen(i)*R.gen(j))
            else:
                M[i,j]=Q.coefficient(R.gen(i)*R.gen(j))
   
    return QuadraticForm(M)



def has_modular_unit_of_degree(G,deg,rational = True, verbose = False,qfminim_bound = 10**5,l2_step=0):
    """
    Returns True,v if the modular curve X(G) has a modular unit v of degree equal to deg, and false,None otherwise.
    
    INPUT:
        
    - ``G`` - a congruence subgroup
    - ``deg`` - int, the degree of modular unit to search for
    - ``rational`` - bool, true means modular unit should be defined over QQ
    - ``verbose`` - bool (default = false), wether or not to print progress
    - ``qfminim_bound`` - int (default - 10^5), given to pari's qfminim command, and is an upper bound on
                          how many vectors of short l2 norm are returned by pari
                          this function will raise an error if pari finds more short
                          vectors then it returns
    - ``l2_step`` - int (default = 0) If l2_step>0 this function first searches the modular units with l2 norm equal to l2_step
                                      then 2*l2_step, 3*l2_step, e.t.c. instead of searching all vectors with l2 norm 2*deg^2.
                                      The l2 norm of a modular unit with divisor n1*c1 + ... + nk*ck is the l2 norm of (n1,...nk). 
    """
    if rational:
        L,D=rational_modular_unit_lattice(G)
    else:
        L,D=modular_unit_lattice(G)
    
    M = L.basis_matrix().change_ring(ZZ).LLL()
    for v in M:
        if v.norm(1)/2 == deg:
            return True,L(v)
       
    GS_matrix=M*M.transpose()
    pari_gs=pari(GS_matrix)
    
    
    #just to speed up positive results
    if l2_step > 0:
        for l2 in range(l2_step,deg**2*2-l2_step+1,l2_step):
            short_vectors=pari_gs.qfminim(l2,qfminim_bound)
            if verbose:
                print short_vectors[:2]
            
            count = 0
            for i in short_vectors[2]:
                count+=1
                if verbose and count%10000==0:
                    print count
                v=vector(QQ,i.list())*M
                if v.norm(1)/2 == deg:
                    return True,L(v)
    
    
    short_vectors=pari_gs.qfminim(deg**2*2,qfminim_bound)
    
    if verbose:
        print short_vectors[:2]
    
    count = 0
    for i in short_vectors[2]:
        count+=1
        if verbose and count%10000==0:
            print count
        v=vector(QQ,i.list())*M
        if v.norm(1)/2 == deg:
            return True,L(v)
    assert short_vectors[0].sage() < 2*qfminim_bound
    return False,None

def modular_units_of_degree(G,deg,rational = True, verbose = False,qfminim_bound = 10**5):
    """
    Returns an iterator over all modular units on the curve X(G).
    
    INPUT:
        
    - ``G`` - a congruence subgroup
    - ``deg`` - int, the degree of modular unit to search for
    - ``rational`` - bool, true means modular unit should be defined over QQ
    - ``verbose`` - bool (default = false), wether or not to print progress
    - ``qfminim_bound`` - int (default - 10^5), given to pari's qfminim command, and is an upper bound on
                          how many vectors of short l2 norm are returned by pari
                          this function will raise an error if pari finds more short
                          vectors then it returns

    """
    if rational:
        L,D=rational_modular_unit_lattice(G)
    else:
        L,D=modular_unit_lattice(G)
    
    M = L.basis_matrix().change_ring(ZZ).LLL()
    
       
    GS_matrix=M*M.transpose()
    pari_gs=pari(GS_matrix)
    
    
    short_vectors=pari_gs.qfminim(deg**2*2,qfminim_bound)
    
    if verbose:
        print short_vectors[:2]
    
    count = 0
    for i in short_vectors[2]:
        count+=1
        if verbose and count%10000==0:
            print count
        v=vector(QQ,i.list())*M
        if v.norm(1)/2 == deg:
            yield L(v)
    assert short_vectors[0].sage() < 2*qfminim_bound
