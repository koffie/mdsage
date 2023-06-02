from collections.abc import Sized

from sage.all import CuspForms, Matrix, RR, ZZ, prime_divisors, PolynomialRing

from .utilities import CPUTimeLogger


def kernel_matrix(M):
    """
    Computing the kernel of a matrix over ZZ can sometimes be slow in sage.
    Especially since sage by default return the basis of the kernel in hermite
    normal form.

    This function can sometimes be faster in computing the kernel of an integer
    matrix. The returned matrix forms a basis of the kernel, however the output
    matrix is not put into hermite normal form.

    EXAMPLES::

        sage: from mdsage import *
        sage: M = Matrix(6, 4, range(24))
        sage: K = kernel_matrix(M); K
        [-1  0  3 -2  0  0]
        [ 0 -1  2 -1  0  0]
        [ 0  0 -2  3  0 -1]
        [ 0  0 -1  2 -1  0]
        sage: V = M.kernel(); V
        Free module of degree 6 and rank 4 over Integer Ring
        Echelon basis matrix:
        [ 1  0  0  0 -5  4]
        [ 0  1  0  0 -4  3]
        [ 0  0  1  0 -3  2]
        [ 0  0  0  1 -2  1]
        sage: K.row_module() == V
        True


    """
    logger = CPUTimeLogger(f"{__name__}.kernel_matrix")

    S, U, V = M.smith_form()
    logger.log("smith form")

    K = S.transpose().right_kernel_matrix(basis="computed")
    logger.log(f"kernel {RR(len(K.nonzero_positions()))/K.ncols()/K.nrows()}")

    K = K.sparse_matrix()
    U = U.sparse_matrix()
    logger.log("to sparse")

    prod = K * U
    logger.log("product")

    return prod


def vanishing_quadratic_forms(G):
    r"""
    On input of a congruence subgroup G of genus g compute the space of all
    quadratic forms that vanish on the image of the canonical embedding
    $X(G) \\to \\P^{g-1}$.

    EXAMPLES::

        sage: from mdsage import *
        sage: vanishing_quadratic_forms(Gamma0(47))
        [x2^2 - x1*x3 - 2*x2*x3 + 3*x3^2,
         x1*x2 - 2*x2^2 - x0*x3 + 3*x2*x3 - 4*x3^2,
         x1^2 - x0*x2 - 4*x2^2 + 8*x2*x3 - 9*x3^2]


    """
    logger = CPUTimeLogger(f"{__name__}.vanishing_quadratic_forms")

    g = G.genus()
    S = CuspForms(G, 2)
    logger.log("cuspforms")

    prec = 4 * g - 1
    S_basis = [f.qexp(prec) for f in S.integral_basis()]
    logger.log("cuspforms basis")

    base_ring = ZZ
    R = PolynomialRing(base_ring, g, names="x")

    qv2 = []
    monomials2 = []
    for i in range(g):
        for j in range(i, g):
            qv2.append(S_basis[i] * S_basis[j])
            monomials2.append(R.gen(i) * R.gen(j))
    M2 = Matrix(ZZ, [f.padded_list(prec) for f in qv2])

    logger.log("evaulate monomials")

    # The kernel of M2 doesn't change if we remove columns in such a way that
    # the rank doesn't change. It does speed up later computations so we remove
    # some superfluous columns.
    rank = M2.rank()
    rank_new = 0
    cols = rank + 2
    while rank != rank_new:
        M2_small = M2.matrix_from_columns(range(2, cols))
        rank_new = M2_small.rank()
        cols += rank - rank_new
    M2 = M2_small
    logger.log("smaller_matrix")

    quadrics_matrix = kernel_matrix(M2.change_ring(base_ring))
    logger.log("quadrics matrix")
    quadrics = []
    for coeffs in quadrics_matrix.rows():
        quadrics.append(sum(c * x for c, x in zip(coeffs, monomials2)))
    return quadrics


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


def trigonal_or_plane_quintic_primes(G):
    """
    On input of a congruence subgroup G compute the primes p such that
    the reduction associated modular curve $X(G)_{F_p}$ geometrically
    has trigonal or plane quintic reduction. If the output is [0] then
    $X(G)_{Q}$ is trigonal or plane quintic

    WARNING::

    This function assumes that neither $X(G)_{Q}$ or $X(G)_{F_p}$ ore
    sub hyperelliptic. This should be checked beforehand using
    :meth:`mdsage.canonical_rings.subhyperelliptic_primes`.

    EXAMPLES::

        sage: from mdsage import *

    A curve that is not trigonal in any characteristic:

        sage: trigonal_or_plane_quintic_primes(Gamma0(42))
        []

    A trigonal curve over Q:

        sage: trigonal_or_plane_quintic_primes(Gamma0(43))
        [0]

    The following is the only nontrigonal modular curve of shimura type that becomes trigonal mod a prime::

        sage: trigonal_or_plane_quintic_primes(Gamma0(73))
        [2]

    """
    logger = CPUTimeLogger(f"{__name__}.trigonal_or_plane_quintic_primes")
    g = G.genus()
    if g <= 4:
        return [0]
    quadrics = vanishing_quadratic_forms(G)
    R = quadrics[0].parent()
    logger.log("vanishing_quadratic_forms")

    derivatives = [[f(x0=1).derivative(x) for x in R.gens()[1:]] for f in quadrics]
    logger.log("derivatives")

    tangent_matrix = Matrix(
        [[f.constant_coefficient() for f in row] for row in derivatives]
    )
    logger.log("tangent matrix")

    invariants = tangent_matrix.elementary_divisors()
    logger.log("elementary divisors")

    divisors = prime_divisors(invariants[g - 3])
    logger.log("prime divisors")
    return divisors
