from collections import defaultdict

from sage.all import kronecker_symbol, gcd, prime_to_m_part, prime_divisors, valuation

from .class_numbers import (
    small_class_number_discriminants,
    _small_class_number_cache,
    class_number,
)


def _c1(p, d, N):
    # assert is_prime(p)
    # assert N%d == 0
    # assert gcd(d,N//d) == 1
    # assert (N//d)%p==0
    # assert d%p != 0 #consequence of previous assertions actually

    if p > 2:
        return _c2(p, d)

    if d % 4 == 1:
        if N % 4 == 0:
            return 0
        return 1

    if N % 8 == 0:
        return 3 * (1 + kronecker_symbol(-d, p))

    if N % 4 == 0:
        return 3 + kronecker_symbol(-d, p)

    return 2


def _c2(p, d):
    # assert d >= 4
    # assert d%4 == 3
    # assert is_prime(p)
    # assert d%p != 0
    # assert p>2 or d%4==3

    if d % 4 == 3:
        return 1 + kronecker_symbol(-d, p)
    return 1 + kronecker_symbol(-4 * d, p)


def atkin_lehner_ramification_degree(N, d):
    """
    Returns the ramification degree of X_0(N) -> X_0(N)/w_d. See:
    Hyperelliptic Quotients of Modular Curves X_0(N)
    by Masahiro FURUMOTO and Yuji HASEGAWA
    for correctnes. If gcd(d, N/d) != 1 it replaces d by the largest
    divisor of N that has still the same prime divisors as d.

    INPUT:

        - N - an positive integer
        - d - an positive integer dividing N

    OUTPUT:

        The ramification degree of X_0(N) -> X_0(N)/w_d

    EXAMPLES::

        sage: from mdsage import *
        sage: [(d, atkin_lehner_ramification_degree(105, d)) for d in 105.divisors()]
        [(1, 0), (3, 0), (5, 8), (7, 0), (15, 0), (21, 8), (35, 16), (105, 8)]
        sage: G = Gamma0(105)
        sage: M = ModularSymbols(G,sign=1)
        sage: S = M.cuspidal_submodule()
        sage: S.dimension()
        13
        sage: [(d, (S.atkin_lehner_operator(d) - 1).kernel().dimension()) for d in 105.divisors()]
        [(1, 13), (3, 7), (5, 5), (7, 7), (15, 7), (21, 5), (35, 3), (105, 5)]

    """
    assert d > 0 and N > 0
    if d == 1:
        return 0

    assert N % d == 0
    if not gcd(d, N // d) == 1:
        d1 = prime_to_m_part(N, d)
        d = N // d1

    assert gcd(d, N // d) == 1

    divs_N_d = prime_divisors(N // d)

    c1_prod = 1

    for p in divs_N_d:
        c1_prod *= _c1(p, d, N)

    v_d = c1_prod * class_number(-4 * d)

    c2_prod = 1

    if d == 2:
        for p in divs_N_d:
            c2_prod *= 1 + kronecker_symbol(-4, p)
        v_d += c2_prod

    elif d == 3:
        for p in divs_N_d:
            c2_prod *= 1 + kronecker_symbol(-3, p)
        v_d += c2_prod

    elif d == 4:
        for p in divs_N_d:
            v = valuation(N, p)
            c2_prod *= p ** (v // 2) + p ** ((v - 1) // 2)
        v_d += c2_prod

    elif d > 4 and d % 4 == 3:
        for p in divs_N_d:
            c2_prod *= _c2(p, d)

        v_d += c2_prod * class_number(-d)

    return v_d


def small_ramification(bound):
    """
    returns a dictionary of all N for such that the ramification degree of $X_0(N) \to X_0(N)^+$ is small

    INPUT:

        - bound (optional) - an integer up to which ramification degree to compute the levels N

    OUTPUT:

        - a dictionary whose keys are ramification degrees and whose values are lists of levels N
          such that $X_0(N) \to X_0(N)^+$  has that ramification degree.

    EXAMPLES::

        sage: from mdsage import *
        sage: small_ramification(2)
        {0: [1],
         2: [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 18, 22, 25, 28, 37, 58]}
        sage: # if the following fails the data in our article is also incorrect
        sage: {(k,len(v),max(v)) for k,v in small_ramification(100).items()}
        {(0, 1, 1),
         (2, 18, 58),
         (4, 48, 253),
         (6, 32, 652),
         (8, 128, 1012),
         (10, 39, 1318),
         (12, 173, 2608),
         (14, 35, 2293),
         (16, 329, 4048),
         (18, 62, 5692),
         (20, 225, 5377),
         (22, 40, 6637),
         (24, 576, 10432),
         (26, 63, 11302),
         (28, 257, 13297),
         (30, 88, 14422),
         (32, 790, 18748),
         (34, 56, 18397),
         (36, 482, 22768),
         (38, 74, 30493),
         (40, 785, 30178),
         (42, 130, 29437),
         (44, 375, 34318),
         (46, 78, 47338),
         (48, 1618, 41728),
         (50, 77, 43717),
         (52, 389, 50317),
         (54, 125, 48742),
         (56, 992, 62302),
         (58, 77, 48778),
         (60, 817, 83218),
         (62, 96, 85402),
         (64, 1857, 106177),
         (66, 175, 92698),
         (68, 493, 102958),
         (70, 127, 94378),
         (72, 1963, 134773),
         (74, 104, 91228),
         (76, 506, 121972),
         (78, 170, 92458),
         (80, 2309, 120712),
         (82, 118, 151237),
         (84, 1019, 166798),
         (86, 106, 137197),
         (88, 1413, 150382),
         (90, 238, 149053),
         (92, 577, 189352),
         (94, 112, 184438),
         (96, 4289, 198958),
         (98, 132, 161302),
         (100, 844, 200722)}

    """
    discriminant_to_cls_nr = _small_class_number_cache(bound)
    small_class_nr_discriminants = small_class_number_discriminants(bound)
    _small_ramification = defaultdict(list)
    todo = sum(small_class_nr_discriminants.values(), [])
    todo = [-(N // 4) for N in todo if N % 4 == 0]
    for N in todo:
        ramification = discriminant_to_cls_nr[-4 * N]
        if N % 4 == 3:
            ramification += discriminant_to_cls_nr[-N]
        if N <= 4:
            ramification = 2
        if N == 1:
            ramification = 0
        if ramification <= bound:
            _small_ramification[ramification].append(N)
    return {k: sorted(v) for k, v in _small_ramification.items()}
