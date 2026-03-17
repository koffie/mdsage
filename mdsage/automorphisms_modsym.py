"""
Provides functions for computing automorphisms of modular curves
using modular symbols. The main focus is on automorphisms that are
not yet implemented in Sage. This includes automorphisms of X_0(N)
beyond Atkin-Lehner operators, and the extra automorphism V_5 on
X_0(25M)/w_25.
"""

from sage.all import Matrix, matrix, ZZ, Integers, ModularSymbols, PolynomialRing


def apply(m, g, ambient=False):
    """
    Applies the automorphism g to the modular symbol m.

    The automorphism g should be given as a Möbius transformation, represented
    either as a 2x2 matrix or as a tuple of four integers. The modular symbol m
    should be an element of a modular symbols space. The function returns the
    result of applying g to m, as an element of the same modular symbols space
    (or its ambient space if ambient=True).

    Note: This function does not check whether g is actually an automorphism of
    the modular curve associated to m, so it is the caller's responsibility to
    ensure that g is a valid automorphism.

    EXAMPLES::

        sage: from mdsage.automorphisms_modsym import *
        sage: M = ModularSymbols(37)
        sage: g = matrix(2,2,[0,-1,37,0])
        sage: all([apply(m,g) == M.atkin_lehner_operator(37)(m) for m in M.basis()])
        True
        sage: all([apply(m,[0,-1,37,0]) == M.atkin_lehner_operator(37)(m) for m in M.basis()])
        True

    """
    try:
        g = g.list()
    except AttributeError:
        g = list(g)
    M = m.parent()
    if ambient:
        M = M.ambient()
    gm = 0
    for a, c in m.modular_symbol_rep():
        gm += a * c.apply(g)
    return M(gm)


def upsilon(N):
    """
    Helper function that creates a matrix used in the definition of the automorphism V_5 on X_0(25M)/w_25.
    """
    return Matrix(ZZ, 2, [N, 0, 0, 1])


def b_j(N):
    """
    Helper function that creates a matrix used in the definition of the automorphism V_5 on X_0(25M)/w_25.
    """
    assert N.valuation(5) == 2
    N1 = N // 25
    j = (2 * Integers(5)(N // 25) ** (-1)).lift()
    return Matrix(ZZ, 2, [N // 25 * j + 1, -j, -N // 25, 1])


def V5_operator(S, check=True):
    """
    Constructs the automorphism V_5 on X_0(25M)/w_25, where S is the modular symbols space associated to X_0(25M)/w_25.

    EXAMPLES::

        sage: from mdsage import *
        sage: M = ModularSymbols(5**2*7)
        sage: S = M.cuspidal_subspace()
        sage: S_25 = invariant_subspace(S, [S.atkin_lehner_operator(25)])
        sage: V5 = V5_operator(S_25)
        sage: V5.matrix().charpoly().factor()
        (x - 1)^4 * (x^2 + x + 1)^6

    Note: The eigenvalues of V_5 are 1 and the primitive cube roots of unity, which is consistent with V_5 being an automorphism of order 3.
    The eigenspace for eigenvalue 1 has dimension 4, which is twice the genus of (X_0(25*7)/w_25).V_5, as expected.
    """
    N = S.level()
    M = S.ambient()
    assert N.valuation(5) == 2
    assert S.atkin_lehner_operator(25).matrix() == 1

    g0 = upsilon(5) ** (-1) * b_j(N) * upsilon(5)
    g = (g0 * g0.denominator()).list()
    rows = [apply(M(m), g) for m in S.basis()]
    V5 = S.hom([S(M.atkin_lehner_operator(25)(r) + r) / 2 for r in rows])
    return V5


def invariant_subspace(S, autormorphisms):
    """
    Computes the subspace of S that is invariant under the given automorphisms.

    EXAMPLES::

        sage: from mdsage import *
        sage: M = ModularSymbols(5**2*7)
        sage: S = M.cuspidal_subspace()
        sage: S_25 = invariant_subspace(S, [S.atkin_lehner_operator(25)])
        sage: S_25.dimension()
        16
        sage: V5 = V5_operator(S_25)
        sage: invariant_subspace(S_25, [V5]).dimension()
        4
    """
    S_inv = S
    for aut in autormorphisms:
        S_inv = S_inv.intersection((aut - 1).kernel())
    return S_inv


def V5_invariant_modular_symbols(N):
    """
    Constructs the cuspidal modular symbols space associated to the quotient of X_0(25N)/w_25 by the automorphism V_5.

    EXAMPLES::

        sage: from mdsage import *
        sage: S = V5_invariant_modular_symbols(7)
        sage: S.dimension()
        4
    """
    assert N % 5 != 0, "N = {N} must not be divisible by 5"
    M = ModularSymbols(25 * N)
    S = M.cuspidal_subspace()
    S_25 = (S.atkin_lehner_operator(25) - 1).kernel()
    V5 = V5_operator(S_25)
    return (V5 - 1).kernel()


def modsym_frobenius_polynomial(S, p):
    """
    Computes the Frobenius polynomial as acting on the modular symbols space S.

    EXAMPLES::

        sage: from mdsage import *
        sage: S = V5_invariant_modular_symbols(7)
        sage: f = modsym_frobenius_polynomial(S, 2); f
        x^4 - x^3 + 3*x^2 - 2*x + 4
        sage: f.is_weil_polynomial()
        True
    """
    assert (
        S.weight() == 2
    ), "Currently only implemented for weight 2 modular symbols spaces"
    assert (
        S.character().order() == 1
    ), "Currently only implemented for trivial character"
    N = S.level()
    assert N % p != 0, "p = {p} must not divide N = {N}"
    A = S.abelian_variety()
    # f = A.frobenius_polynomial(p)
    # the above has a bug! Instead, we compute the characteristic polynomial of Frobenius manually in terms of the hecke operator at p.
    ZZpoly = PolynomialRing(ZZ, "x")
    x = ZZpoly.gens()[0]
    Tp = A.hecke_polynomial(p)
    f = ZZpoly(x ** Tp.degree() * Tp(x + p / x))
    assert f.degree() == S.dimension()
    return f


def modsym_point_count(S, p, n=1):
    """
    This function assumes that S is the weight 2 modular symbols space associated to a (quotient) modular curve X.
    It counts the number of points on X of this quotient over the finite field with p^n elements.
    The answer is computed using the Frobenius polynomial as acting on S, which is the characteristic polynomial of Frobenius on the Jacobian of X, and applying the Weil conjectures.
    If S is not corresponding to a modular curve, the result may not be meaningful and potentially negative.

    EXAMPLES::

        sage: from mdsage import *
        sage: S = V5_invariant_modular_symbols(6)
        sage: modsym_point_count(S, 7)
        6
        sage: modsym_point_count(S, 7, n=2)
        80
    """
    f = modsym_frobenius_polynomial(S, p)
    return p**n + 1 - (matrix.companion(f) ** n).trace()
