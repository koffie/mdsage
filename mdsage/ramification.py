from collections import defaultdict

from .class_numbers import small_class_number_discriminants, _small_class_number_cache

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
    todo = sum(small_class_nr_discriminants.values(),[])
    todo = [-(N//4) for N in todo if N%4==0]
    for N in todo:
        ramification = discriminant_to_cls_nr[-4*N]
        if N%4==3:
            ramification += discriminant_to_cls_nr[-N]
        if N<=4:
            ramification=2
        if N == 1:
            ramification=0
        if ramification <= bound:
            _small_ramification[ramification].append(N)
    return {k:sorted(v) for k,v in _small_ramification.items()}