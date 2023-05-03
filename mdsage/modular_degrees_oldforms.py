from sage.all import matrix, Gamma0, Hom, ModularSymbols
from sage.modular import abvar


def product_isogeny_map(E, N):
    """
    For a modular elliptic curve E of level M give and an integer N such that M divides N,
    this computes the map $E^n \to J_0(N)$ where $n$ is the number of divisors of N/M.
    The $i$-th coordinate of this map is given by the degeneracy map $J_0(M) \to J_0(N)$
    corresponding to the $i$-th divisor of N/M.

    Examples::

        sage: from mdsage.modular_degrees_oldforms import *
        sage: GN = Gamma0(37)
        sage: SN = GN.modular_symbols().cuspidal_subspace()
        sage: SE = SN.decomposition()[0]
        sage: E = SE.abelian_variety()
        sage: xi = product_isogeny_map(E,6*37); xi
        Abelian variety morphism:
          From: Abelian subvariety of dimension 4 of J0(37) x J0(37) x J0(37) x J0(37)
          To:   Abelian variety J0(222) of dimension 35
        sage: xi.kernel()[0].invariants()
        []

    """
    M = E.conductor()
    assert N % M == 0
    GN = Gamma0(N)
    SN = GN.modular_symbols().cuspidal_subspace()
    J0N = SN.abelian_variety()

    EtoJ0M = E.ambient_morphism()

    J0M = EtoJ0M.codomain()

    divs = (N // M).divisors()
    n = len(divs)
    En = E**n  # E**n = E^n
    J0M_to_N = [J0M.degeneracy_map(N, d) for d in divs]

    degeneracy_isogenies = []

    for M_to_N in J0M_to_N:
        degeneracy_isogenies.append([(M_to_N * EtoJ0M).matrix()])

    product_isogeny_mat = matrix.block(degeneracy_isogenies, subdivide=False)
    product_isogeny = abvar.morphism.Morphism(Hom(En, J0N), product_isogeny_mat)
    return product_isogeny


def modular_symbol_elliptic_curves(N, sign=None):
    """
    INPUT:

    - N - an integer giving the conductor of the elliptic curves to be returned
    - sign (optional) - an integer that should be +/-1. Only elliptic curves where
      atkin lehner acts with this sign on the associated newform are returned.
      If no sign is given then both elliptic curves of both signs will be returned.

    OUTPUT:

    An iterator of all strong weil elliptic curves of conductor N.
    The returned curves are returned as modular abelian varieties associated to modular
    symbol spaces rather then elliptic curves given by a Weierstrass equation.


    Examples::

        sage: from mdsage.modular_degrees_oldforms import *
        sage: list(modular_symbol_elliptic_curves(37))
        [Abelian subvariety of dimension 1 of J0(37),
         Abelian subvariety of dimension 1 of J0(37)]

        sage: list(modular_symbol_elliptic_curves(37,sign=1))
        [Abelian subvariety of dimension 1 of J0(37)]
        sage: list(modular_symbol_elliptic_curves(37,sign=-1))
        [Abelian subvariety of dimension 1 of J0(37)]

    """
    G = Gamma0(N)
    M = ModularSymbols(G)
    S = M.cuspidal_subspace()
    Snew = S.new_subspace()
    for Si in Snew.decomposition():
        if Si.dimension() != 2:
            continue
        if sign and Si.atkin_lehner_operator().matrix() != sign:
            continue
        yield Si.abelian_variety()


def modular_symbol_elliptic_curves_range(N, sign=None):
    """
    Does the same as `modular_symbol_elliptic_curves` but then all ellipic curves
    of conductor < N instead of equal to N.

    Examples::

        sage: from mdsage.modular_degrees_oldforms import *
        sage: list(modular_symbol_elliptic_curves_range(37))
         [Abelian variety J0(11) of dimension 1,
          Abelian variety J0(14) of dimension 1,
          Abelian variety J0(15) of dimension 1,
          Abelian variety J0(17) of dimension 1,
          Abelian variety J0(19) of dimension 1,
          Abelian variety J0(20) of dimension 1,
          Abelian variety J0(21) of dimension 1,
          Abelian variety J0(24) of dimension 1,
          Abelian subvariety of dimension 1 of J0(26),
          Abelian subvariety of dimension 1 of J0(26),
          Abelian variety J0(27) of dimension 1,
          Abelian subvariety of dimension 1 of J0(30),
          Abelian variety J0(32) of dimension 1,
          Abelian subvariety of dimension 1 of J0(33),
          Abelian subvariety of dimension 1 of J0(34),
          Abelian subvariety of dimension 1 of J0(35),
          Abelian variety J0(36) of dimension 1]

    """
    for M in range(1, N):
        yield from modular_symbol_elliptic_curves(M, sign=sign)
