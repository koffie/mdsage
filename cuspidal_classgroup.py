from maartens_sage_functions import *

def intersection(self,other):
    """Returns the intersection of two quotient modules self = V1/W and other V2/W
    V1 and V2 should be submodulus of the same ambient module.
    
    EXAMPLE::
        
        sage: W = (ZZ^3).scale(100)
        sage: V1 = (ZZ^3).span([[5,5,0],[0,10,20],[0,0,50]]) + W
        sage: V2 = (ZZ^3).span([[5,0,0],[0,25,50]]) + W
        sage: V1/W;V2/W
        Finitely generated module V/W over Integer Ring with invariants (2, 10, 20)
        Finitely generated module V/W over Integer Ring with invariants (4, 20)
        sage: intersection(V1/W,V2/W)
        Finitely generated module V/W over Integer Ring with invariants (2, 4)
    
    """
    assert self.W()==other.W()
    return self.V().intersection(other.V()).quotient(self.W())


def cuspidal_rational_subgroup_mod_rational_cuspidal_subgroup(G):
    """On input a congruence subgroup G that contains Gamma1(N) checks whether there are cuspsums that are not equivalent to a cuspsum defined over Q but that are equivalent to a divisor defined over Q.
    
    INPUT:
    
    - G - a congruence subgroup
    
    OUTPUT:
    
    - the invariants of the cuspidal rational subgroup divided out by the rational cuspidal subgroup 
    
    EXAMPLES::
        
        sage: cuspidal_rational_subgroup_mod_rational_cuspidal_subgroup(Gamma1(26))
        ()
    """
    if G.genus() == 0:
        return tuple()
        
    N=G.level()
    ZZcusp=ZZ**G.ncusps()
    unit_gens=Integers(N).unit_gens()
    L,D=modular_unit_lattice(G)
    Lrat,Drat=rational_modular_unit_lattice(G)
    DmodL=D.quotient(L)
    DmodLrat=Drat.quotient(Lrat)
    kernels=[]
    for g in unit_gens:
        m=Matrix([ZZcusp.gen(G.cusps().index(G.reduce_cusp(galois_action(c,g,N)))) for c in G.cusps() ])
        f = DmodL.hom([DmodL(i.lift()*(m-1)) for i in  DmodL.gens()])
        kernels.append(f.kernel())
        #print kernels[-1]
    rat_cusp_tors=reduce(intersection,kernels)
    return rat_cusp_tors.V().quotient(Drat+L).invariants()
