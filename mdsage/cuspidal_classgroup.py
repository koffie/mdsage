"""
A module for checking wether the torsion in a modular jacobian is
generated by cusps. It returns a multiplicative upperbound on the index of
the cuspidal rational group in the rational group

EXAMPLES::

    sage: from mdsage import *
    sage: d = rational_cuspidal_classgroup(Gamma1(23)).cardinality()
    sage: upper_bound_index_cusps_in_JG_torsion(Gamma1(23),d)
    1

"""

from sage.all import (Cusp,
                      gcd,
                      infinity,
                      Integers,
                      prime_range,
                      oo,
                      prod,
                      QQ,
                      Matrix,
                      ModularSymbols,
                      ZZ)

from .maartens_sage_functions import (
    integral_period_mapping,
    cuspidal_integral_structure_matrix,
    galois_action,
    galois_orbits,
    modular_unit_lattice,
    rational_modular_unit_lattice
)
from functools import reduce

def intersection(self, other):
    """Returns the intersection of two quotient modules self = V1/W and other V2/W
    V1 and V2 should be submodulus of the same ambient module.

    EXAMPLES::

        sage: from mdsage import *
        sage: W = (ZZ^3).scale(100)
        sage: V1 = (ZZ^3).span([[5, 5, 0], [0, 10, 20], [0, 0, 50]]) + W
        sage: V2 = (ZZ^3).span([[5, 0, 0], [0, 25, 50]]) + W
        sage: V1/W;V2/W
        Finitely generated module V/W over Integer Ring with invariants (2, 10, 20)
        Finitely generated module V/W over Integer Ring with invariants (4, 20)
        sage: intersection(V1/W,V2/W)
        Finitely generated module V/W over Integer Ring with invariants (2, 4)

    """
    assert self.W() == other.W()
    return self.V().intersection(other.V()).quotient(self.W())


def cuspidal_rational_subgroup_mod_rational_cuspidal_subgroup(G):
    """On input a congruence subgroup G that contains Gamma1(N) checks whether
    there are cuspsums that are not equivalent to a cuspsum defined over `\\QQ`
    but that are equivalent to a divisor defined over `\\QQ`.

    INPUT:

    - G - a congruence subgroup

    OUTPUT:

    - the invariants of the cuspidal rational subgroup divided out by the rational cuspidal subgroup

    EXAMPLES::

        sage: from mdsage import *
        sage: cuspidal_rational_subgroup_mod_rational_cuspidal_subgroup(Gamma1(26))
        ()
    """
    if G.genus() == 0:
        return tuple()

    N = G.level()
    ZZcusp = ZZ ** G.ncusps()
    unit_gens = Integers(N).unit_gens()
    L, D = modular_unit_lattice(G)
    Lrat, Drat = rational_modular_unit_lattice(G)
    DmodL = D.quotient(L)
    #DmodLrat = Drat.quotient(Lrat)
    kernels = []
    for g in unit_gens:
        m = Matrix([ZZcusp.gen(G.cusps().index(G.reduce_cusp(galois_action(c, g, N))))
                    for c in G.cusps()])
        f = DmodL.hom([DmodL(i.lift() * (m-1)) for i in  DmodL.gens()])
        kernels.append(f.kernel())
        #print kernels[-1]
    rat_cusp_tors = reduce(intersection, kernels)
    return rat_cusp_tors.V().quotient(Drat + L).invariants()

def upper_bound_index_cusps_in_JG_torsion(G, d, bound=60):
    """
    INPUT:

    - G - a congruence subgroup
    - d - integer, the size of the rational cuspidal subgroup
    - bound (optional, default = 60) - the bound for the primes p up to which to use
      the hecke matrix `T_p - <p> - p` for bounding the torsion subgroup

    OUTPUT:

    - an integer `i` such that `(\\#J_G(\\QQ)_{tors})/d` is a divisor of `i`.

    EXAMPLES::

        sage: from mdsage import *
        sage: d = rational_cuspidal_classgroup(Gamma1(23)).cardinality()
        sage: upper_bound_index_cusps_in_JG_torsion(Gamma1(23),d)
        1

    """
    N = G.level()
    M = ModularSymbols(G)
    Sint = cuspidal_integral_structure(M)
    kill_mat = M.star_involution().matrix().restrict(Sint) - 1
    kill = kill_mat.transpose().change_ring(ZZ).row_module()

    for p in prime_range(3, bound):
        if not N % p == 0:
            kill += kill_torsion_coprime_to_q(p, M).restrict(Sint).change_ring(ZZ).transpose().row_module()
        if kill.matrix().is_square() and kill.matrix().determinant() == d:
            break
    kill_mat = kill.matrix().transpose()
    #print N,"index of torsion in stuff killed",kill.matrix().determinant()/d

    if kill.matrix().determinant() == d:
        return 1
    pm = integral_period_mapping(M)
    period_images1 = [sum([M.coordinate_vector(M([c, infinity])) for c in cusps]) * pm
                      for cusps in galois_orbits(G)]

    m = (Matrix(period_images1) * kill_mat).stack(kill_mat)
    diag = m.change_ring(ZZ).echelon_form().diagonal()
    #print diag,prod(diag)
    assert prod(diag) == kill.matrix().determinant() / d

    period_images2 = [M.coordinate_vector(M([c, infinity]))*pm
                      for c in G.cusps() if c != Cusp(oo)]
    m = (Matrix(period_images2) * kill_mat).stack(kill_mat)
    m, denom = m._clear_denom()
    diag = (m.change_ring(ZZ).echelon_form()/denom).diagonal()
    #print diag
    #print prod(i.numerator() for i in diag),"if this is 1 then :)"
    return prod(i.numerator() for i in diag)


def JG_torsion_upperbound(G, bound=60):
    """
    INPUT:

    - G - a congruence subgroup
    - bound (optional, default = 60) - the bound for the primes p up to which to use
      the hecke matrix `T_p - <p> - p` for bounding the torsion subgroup

    OUTPUT:

    - A subgroup of `(S_2(G) \\otimes \\QQ) / S_2(G)` that is guaranteed to contain
      the rational torison subgroup, together with a subgroup generated by the
      rational cusps.
      The subgroup is given as a subgroup of `S_2(G)/NS_2(G)` for a suitable integer N

    EXAMPLES::

        sage: from mdsage import *
        sage: d = rational_cuspidal_classgroup(Gamma1(23)).cardinality()
        sage: upper_bound_index_cusps_in_JG_torsion(Gamma1(23),d)
        1

    """
    N = G.level()
    M = ModularSymbols(G)
    Sint = cuspidal_integral_structure(M)
    kill_mat = (M.star_involution().matrix().restrict(Sint) - 1)
    kill = kill_mat.transpose().change_ring(ZZ).row_module()

    for p in prime_range(3, bound):
        if not N % p == 0:
            kill += kill_torsion_coprime_to_q(p, M).restrict(Sint).change_ring(ZZ).transpose().row_module()
        #if kill.matrix().is_square() and kill.matrix().determinant() == d:
        #    #print p
        #    break
    kill_mat = kill.matrix().transpose()
    #print N,"index of torsion in stuff killed",kill.matrix().determinant()/d
    #if kill.matrix().determinant()==d:
    #    return 1
    d = prod(kill_mat.smith_form()[0].diagonal())
    pm = integral_period_mapping(M)
    #period_images1 = [sum([M.coordinate_vector(M([c, infinity])) for c in cusps]) * pm
    #                  for cusps in galois_orbits(G)]
    period_images2 = [M.coordinate_vector(M([c, infinity])) * pm
                      for c in G.cusps() if c != Cusp(oo)]

    m = (Matrix(period_images2) * kill_mat).stack(kill_mat)
    m, d2 = m._clear_denom()
    d = gcd(d, d2)


    #diag=(m.change_ring(ZZ).echelon_form()/denom).diagonal()
    #print diag
    #print prod(i.numerator() for i in diag),"if this is 1 then :)"
    #return prod(i.numerator() for i in diag)



def kill_torsion_coprime_to_q(q, S):
    """This is with respect to the Xmu model"""
    #print "computing diamond",q
    dq = S.diamond_bracket_matrix(q)
    #print "computing tq"
    Tq = S.hecke_matrix(q)
    return Tq - dq - q

def cuspidal_integral_structure(M):
    """
    Returns the integral structure of a modular symbols space."""
    B = cuspidal_integral_structure_matrix(M)
    try:
        B_row_ambient =  B.row_ambient_module(QQ)
    except AttributeError:
        B_row_ambient = B._row_ambient_module(QQ)
    return B_row_ambient.submodule_with_basis(B).change_ring(ZZ)
