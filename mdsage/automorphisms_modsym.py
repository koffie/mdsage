"""
Provides functions for computing automorphisms of modular curves
using modular symbols. The main focus is on automorphisms that are
not yet implemented in Sage. This includes automorphisms of X_0(N)
beyond Atkin-Lehner operators, and the extra automorphism V_5 on
X_0(25M)/w_25.
"""


def apply(m,g,ambient=False):
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
  
    """
    if not isinstance(g, list):
        g = g.list()
    g = list(g)
    M = m.parent()
    if ambient:
        M = M.ambient()
    gm = 0
    for a,c in m.modular_symbol_rep():
        gm += a*c.apply(g)
    return M(gm)

