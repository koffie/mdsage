from collections.abc import Sized

from sage.all import CuspForms, Matrix, ZZ, prime_divisors


def subhyperelliptic_primes(G, N=None):
    r"""
    INPUT:

        - G - A congruence subgroup of $\SL_2(\ZZ)$ or a list of q-expansions of one forms on a curve over $\ZZ$
        - N (optional) - If G is a list of q-expansions then N should be a multiple of all primes of bad reduction.

    OUTPUT:

        - [0,] if the curve X(G) is subhyperelliptic in characteristic 0
          of the list of primes not dividing the level of G such that the X(G) has subhyperelliptic reduction.

    EXAMPLES::

        sage: from mdsage import *
        sage: [subhyperelliptic_primes(Gamma0(N)) for N in [13, 34, 37, 48, 64, 72]]
        [[0], [], [0], [0], [], []]

    The following is the only nonhyperelliptic modular curve of shimura type that becomes hyperelliptic mod a prime::

        sage: subhyperelliptic_primes(GammaH(37,[4]))
        [2]


    """
    if isinstance(G, Sized):
        S_basis = G
        g = len(S_basis)
        prec = 4 * g - 1
        assert not N is None
    else:
        S = CuspForms(G)
        g = S.dimension()
        prec = 4 * g - 1
        S_basis = [f.qexp(prec) for f in S.integral_basis()]
        N = G.level()

    if g <= 2:
        return [0]

    qv = []
    for i in range(g):
        for j in range(i, g):
            qv.append(S_basis[i] * S_basis[j])

    M = Matrix(ZZ, [f.padded_list(prec) for f in qv])
    invariant = M.elementary_divisors()[2 * g - 1]
    invariant = invariant.prime_to_m_part(N) if invariant else 0

    return prime_divisors(invariant) if invariant else [0]
