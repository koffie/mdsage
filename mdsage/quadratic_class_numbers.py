from collections import defaultdict

from sage.all import (
    BinaryQF_reduced_representatives,
    prime_divisors,
    kronecker_symbol,
    valuation,
    cached_function,
    prod,
    euler_phi,
    prime_range,
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


def small_class_number_fundamental_discriminants(bound=100):
    """
    returns a dictionary of all imaginary quadratic fundamental discriminants of small class number

    INPUT:

        - bound (optional) - an integer speciefiel

    OUTPUT:

        - a dictionary whose keys are class numbers and whose values are lists of discriminants having that class number

    EXAMPLES::

        sage: from mdsage import *
        sage: small_class_number_fundamental_discriminants(2)
        {1: [-3, -4, -7, -8, -11, -19, -43, -67, -163],
         2: [-15,
          -20,
          -24,
          -35,
          -40,
          -51,
          -52,
          -88,
          -91,
          -115,
          -123,
          -148,
          -187,
          -232,
          -235,
          -267,
          -403,
          -427]}
        sage: # we verify that our data matched with that of watkins
        sage: {(h,len(discs),min(discs)) for h,discs in small_class_number_fundamental_discriminants(100).items()}
        {(1, 9, -163),
         (2, 18, -427),
         (3, 16, -907),
         (4, 54, -1555),
         (5, 25, -2683),
         (6, 51, -3763),
         (7, 31, -5923),
         (8, 131, -6307),
         (9, 34, -10627),
         (10, 87, -13843),
         (11, 41, -15667),
         (12, 206, -17803),
         (13, 37, -20563),
         (14, 95, -30067),
         (15, 68, -34483),
         (16, 322, -31243),
         (17, 45, -37123),
         (18, 150, -48427),
         (19, 47, -38707),
         (20, 350, -58507),
         (21, 85, -61483),
         (22, 139, -85507),
         (23, 68, -90787),
         (24, 511, -111763),
         (25, 95, -93307),
         (26, 190, -103027),
         (27, 93, -103387),
         (28, 457, -126043),
         (29, 83, -166147),
         (30, 255, -134467),
         (31, 73, -133387),
         (32, 708, -164803),
         (33, 101, -222643),
         (34, 219, -189883),
         (35, 103, -210907),
         (36, 668, -217627),
         (37, 85, -158923),
         (38, 237, -289963),
         (39, 115, -253507),
         (40, 912, -260947),
         (41, 109, -296587),
         (42, 339, -280267),
         (43, 106, -300787),
         (44, 691, -319867),
         (45, 154, -308323),
         (46, 268, -462883),
         (47, 107, -375523),
         (48, 1365, -335203),
         (49, 132, -393187),
         (50, 345, -389467),
         (51, 159, -546067),
         (52, 770, -439147),
         (53, 114, -425107),
         (54, 427, -532123),
         (55, 163, -452083),
         (56, 1205, -494323),
         (57, 179, -615883),
         (58, 291, -586987),
         (59, 128, -474307),
         (60, 1302, -662803),
         (61, 132, -606643),
         (62, 323, -647707),
         (63, 216, -991027),
         (64, 1672, -693067),
         (65, 164, -703123),
         (66, 530, -958483),
         (67, 120, -652723),
         (68, 976, -819163),
         (69, 209, -888427),
         (70, 560, -811507),
         (71, 150, -909547),
         (72, 1930, -947923),
         (73, 119, -886867),
         (74, 407, -951043),
         (75, 237, -916507),
         (76, 1075, -1086187),
         (77, 216, -1242763),
         (78, 561, -1004347),
         (79, 175, -1333963),
         (80, 2277, -1165483),
         (81, 228, -1030723),
         (82, 402, -1446547),
         (83, 150, -1074907),
         (84, 1715, -1225387),
         (85, 221, -1285747),
         (86, 472, -1534723),
         (87, 222, -1261747),
         (88, 1905, -1265587),
         (89, 192, -1429387),
         (90, 801, -1548523),
         (91, 214, -1391083),
         (92, 1248, -1452067),
         (93, 262, -1475203),
         (94, 509, -1587763),
         (95, 241, -1659067),
         (96, 3283, -1684027),
         (97, 185, -1842523),
         (98, 580, -2383747),
         (99, 289, -1480627),
         (100, 1736, -1856563)}


    """
    if bound > 100:
        raise NotImplementedError("only implemented for bound <= 100")
    data = load_json_data("small_class_number_fundamental_discriminant.json")
    return {
        int(h): discriminants for h, discriminants in data.items() if int(h) <= bound
    }


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
    small_class_nr_fundamental_discriminants = (
        small_class_number_fundamental_discriminants(bound)
    )

    if bound > 100 or bound not in small_class_nr_fundamental_discriminants:
        raise NotImplementedError("only implemented for bound <= 100")

    for h in range(1, bound + 1):
        if h not in small_class_nr_fundamental_discriminants:
            raise NotImplementedError("only implemented for bound <= 100")

        fundamental_discriminants = small_class_nr_fundamental_discriminants[h]

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


def small_class_number_discriminants(bound=100):
    """
    returns a dictionary of all imaginary quadratic fundamental discriminants of small class number

    INPUT:

        - bound (optional) - an integer up to which class number to compute the fundamental discriminants

    OUTPUT:

        - a dictionary whose keys are class numbers and whose values are lists of discriminants having that class number

    EXAMPLES::

        sage: from mdsage import *
        sage: small_class_number_discriminants(2)
        {1: [-3, -12, -27, -4, -16, -7, -28, -8, -11, -19, -43, -67, -163],
         2: [-48,
          -75,
          -147,
          -36,
          -64,
          -100,
          -112,
          -32,
          -72,
          -99,
          -15,
          -60,
          -20,
          -24,
          -35,
          -40,
          -51,
          -52,
          -88,
          -91,
          -115,
          -123,
          -148,
          -187,
          -232,
          -235,
          -267,
          -403,
          -427]}
        sage: # if the following fails the data in our article is also incorrect
        sage: {(h,len(discs),min(discs)) for h,discs in small_class_number_discriminants(100).items()}
        {(1, 13, -163),
         (2, 29, -427),
         (3, 25, -907),
         (4, 84, -1555),
         (5, 29, -2683),
         (6, 101, -4075),
         (7, 38, -5923),
         (8, 208, -7987),
         (9, 55, -10627),
         (10, 123, -13843),
         (11, 46, -15667),
         (12, 379, -19723),
         (13, 43, -20563),
         (14, 134, -30067),
         (15, 95, -34483),
         (16, 531, -35275),
         (17, 50, -37123),
         (18, 291, -48427),
         (19, 59, -38707),
         (20, 502, -58843),
         (21, 118, -61483),
         (22, 184, -85507),
         (23, 78, -90787),
         (24, 1042, -111763),
         (25, 101, -93307),
         (26, 227, -103027),
         (27, 136, -103387),
         (28, 623, -126043),
         (29, 94, -166147),
         (30, 473, -137083),
         (31, 83, -133387),
         (32, 1231, -164803),
         (33, 158, -222643),
         (34, 261, -189883),
         (35, 111, -210907),
         (36, 1303, -217627),
         (37, 96, -158923),
         (38, 283, -289963),
         (39, 162, -253507),
         (40, 1418, -274003),
         (41, 125, -296587),
         (42, 595, -301387),
         (43, 123, -300787),
         (44, 909, -319867),
         (45, 231, -308323),
         (46, 328, -462883),
         (47, 117, -375523),
         (48, 2893, -335203),
         (49, 146, -393187),
         (50, 443, -389467),
         (51, 217, -546067),
         (52, 1003, -457867),
         (53, 130, -425107),
         (54, 806, -532123),
         (55, 177, -452083),
         (56, 1809, -494323),
         (57, 237, -615883),
         (58, 360, -586987),
         (59, 144, -474307),
         (60, 2352, -662803),
         (61, 149, -606643),
         (62, 386, -647707),
         (63, 311, -991027),
         (64, 2915, -693067),
         (65, 192, -703123),
         (66, 856, -958483),
         (67, 145, -652723),
         (68, 1227, -819163),
         (69, 292, -888427),
         (70, 702, -821683),
         (71, 176, -909547),
         (72, 4046, -947923),
         (73, 137, -886867),
         (74, 472, -951043),
         (75, 353, -916507),
         (76, 1381, -1086187),
         (77, 236, -1242763),
         (78, 921, -1004347),
         (79, 200, -1333963),
         (80, 3851, -1165483),
         (81, 338, -1030723),
         (82, 486, -1446547),
         (83, 174, -1074907),
         (84, 2990, -1225387),
         (85, 246, -1285747),
         (86, 553, -1534723),
         (87, 313, -1261747),
         (88, 2769, -1265587),
         (89, 206, -1429387),
         (90, 1508, -1548523),
         (91, 249, -1391083),
         (92, 1590, -1452067),
         (93, 354, -1475203),
         (94, 598, -1587763),
         (95, 273, -1659067),
         (96, 7265, -1684027),
         (97, 208, -1842523),
         (98, 707, -2383747),
         (99, 396, -1480627),
         (100, 2310, -1856563)}
    """
    discriminants = defaultdict(list)
    for cls_nr, D, f in cm_orders2(bound):
        discriminants[cls_nr].append(D * f**2)

    return dict(discriminants)


@cached_function
def _small_class_number_cache(bound):
    """
    Examples::

        sage: from mdsage.quadratic_class_numbers import _small_class_number_cache
        sage: _small_class_number_cache(1)
        {-163: 1,
         -67: 1,
         -43: 1,
         -28: 1,
         -27: 1,
         -19: 1,
         -16: 1,
         -12: 1,
         -11: 1,
         -8: 1,
         -7: 1,
         -4: 1,
         -3: 1}

    """
    cls_nrs = {}
    for cls_nr, D, f in cm_orders2(bound):
        cls_nrs[D * f**2] = cls_nr

    return cls_nrs
