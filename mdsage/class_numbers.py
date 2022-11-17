from collections import defaultdict

from sage.all import (
    BinaryQF_reduced_representatives,
    prime_divisors,
    kronecker_symbol,
    valuation,
    cached_function,
    prod,
    euler_phi,
    prime_range
)

from .utilities import load_json_data


def class_number(D):
    """
    INPUT:

    - D - an negative integer that is 0 or 1 mod 4

    OUTPUT:

    - An integer that is the classnumber of the quadratic order of discriminant D

    EXAMPLES::

        sage: from mdsage import *
        sage: class_number(-4)
        1
        sage: class_number(-47)
        5
        sage: class_number(-7*31**2)
        32
        sage: class_number(-7*29**2)
        28

    see http://www.math.toronto.edu/~ila/Cox-Primes_of_the_form_x2+ny2.pdf
    could maybe use theorem 7.24 to speed this up
    """
    assert D < 0
    qfs = BinaryQF_reduced_representatives(D, primitive_only=True, proper=True)
    return len(qfs)


def class_numbers(D, fs, clnr_D=None):
    """
    computes the class number h(D*f**2) for a bunch of f's at the same time
    using theorem 7.24 from http://www.math.toronto.edu/~ila/Cox-Primes_of_the_form_x2+ny2.pdf

    INPUT:

    - D - an negative integer that 0 or 1 mod 4
    - fs - a list of integers
    - clnr_D (optional) - an integer that equals class_number(D) to speed up the computation

    OUTPUT:

    - a list of tuples of integers equalling [(f,class_number(D*f**2)) for f in fs]

    EXAMPLES::

        sage: from mdsage import *
        sage: list(class_numbers(-7, [1, 29, 31]))
        [(1, 1), (29, 28), (31, 32)]

    """
    assert D < 0
    if not clnr_D:
        clnr_D = class_number(D)

    for f in fs:
        clnr_f = clnr_D
        for p in prime_divisors(f):
            s = kronecker_symbol(D, p)
            clnr_f *= (p - s) * p ** (valuation(f, p) - 1)
        if D == -3 and f != 1:
            clnr_f = clnr_f // 3
        if D == -4 and f != 1:
            clnr_f = clnr_f // 2
        yield f, clnr_f


@cached_function
def _small_euler_phis():
    B = prod(prime_range(12))
    assert euler_phi(B) > 300
    small_euler_phi_list = [i for i in range(1, B) if euler_phi(i) <= 100]

    small_euler_phis = defaultdict(list)
    for i in small_euler_phi_list:
        small_euler_phis[euler_phi(i)].append(i)
    return small_euler_phis


def cm_orders2(bound):
    """
    yields all triples (h(D*f**2),D,f) such that h(D*f**2) <= bound

    INPUT:

        - bound - an integer up to which to generate the cm_orders

    OUTPUT:

        an iterator over all triples (h(D*f**2),D,f) such that h(D*f**2) <= bound

    EXAMPLES::

        sage: from mdsage import *
        sage: list(cm_orders2(1))
        [(1, -3, 1),
         (1, -3, 2),
         (1, -3, 3),
         (1, -4, 1),
         (1, -4, 2),
         (1, -7, 1),
         (1, -7, 2),
         (1, -8, 1),
         (1, -11, 1),
         (1, -19, 1),
         (1, -43, 1),
         (1, -67, 1),
         (1, -163, 1)]

    """
    small_euler_phis = _small_euler_phis()
    small_class_nr_fundamental_discriminants = load_json_data("small_class_number_fundamental_discriminant.json")

    if bound > 100 or str(bound) not in small_class_nr_fundamental_discriminants:
        raise NotImplementedError("only implemented for bound <= 100")

    for h in range(1, bound + 1):
        if str(h) not in small_class_nr_fundamental_discriminants:
            raise NotImplementedError("only implemented for bound <= 100")

        fundamental_discriminants = small_class_nr_fundamental_discriminants[str(h)]

        for D in fundamental_discriminants:
            extra = 1
            if D == -3:
                extra = 3
            if D == -4:
                extra = 2
            for i in range(1, (extra * bound) // h + 1):
                phis = small_euler_phis[i]
                for f, cls_nr in class_numbers(D, phis, h):
                    if cls_nr <= bound:
                        yield cls_nr, D, f
