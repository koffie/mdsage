"""
This is a sage version of the original maple code by Mark van Hoeij. This maple
code can be found at:

    http://www.math.fsu.edu/~hoeij/files/X1N/cusp_divisors_program

The code is for computing divisors of the modular unit F_k as defined in:
http://arxiv.org/abs/1307.5719
"""
def Phi(d):
    """
    Function Phi in Mark his code
    """
    d = ZZ(d)
    return euler_phi(d) / d


def inverse_gcd(i,N):
    """
    Function IG in Mark his code
    """
    i = ZZ(i); N = ZZ(N)
    if N == 2*i and (i==1 or i==2):
        return 1
    return N/gcd(i,N)
    
def degree_cusp(i,N):
    """
    Function DegreeCusp in Mark his code 
    
    returns the degree over Q of the cusp $q^(i/n)\zeta_n^j$ on X_1(N)
    """
    i = ZZ(i); N = ZZ(N)
    d = euler_phi(gcd(i,N))
    if i == 0 or 2*i == N:
        return ceil(d/2)
    return d
    
def min_formula(N,t):
    """
    Function MinFormula in Mark his code 
    """
    N = ZZ(N); t = QQ(t)
    if N < 2:
        raise ValueError
    if N == 2:
        return 4 * t -1
    if N == 3:
        return 9 * min(t, ZZ(1)/3) - 8*t
    return sum(N * Phi(gcd(i,N)) * (min(t, i/N) - 4*(i/N)*(1-i/N)*t) for i in srange(1,(N-1)//2+1))
    
def divisor_F_bc(N,k):
    """
    Function Divisor_F_bc in Mark his code, return the divisor of F_k as function on X_1(N) 
    """
    N = ZZ(N); k = ZZ(k)
    return vector([min_formula(k, i/N) * inverse_gcd(i,N) for i in srange(0,N//2+1)])

def LB_c(N):
    N = ZZ(N);
    return [divisor_F_bc(N,k) for k in srange(2,N//2+2)]

def conjectural_cuspidal_classgroup(N):
    N = ZZ(N);
    cups_div = ZZ**N
    cusp_div_0 = cusp_div.span([cusp_div.gen(i)  - degree_cusp(i,N) * cusp_div.gen(1) for i in range(N//2+1)])
    return cusp_div_0.quotient(cusp_div_0.span(LB_c(N)))
    
def cusp_number_from_signature((v2, v3), N):
    """
    The cusp signature of a cusp on X_1(N) is a pair (v2,v3). 
    Where v2 is the valuation of F_2 and v3 the valuation of F_3 
    """
    v2 = ZZ(v2); v3 = ZZ(v3); N = ZZ(N)
    if v3-v2 > 0:
        return v3*N/(4*v3-v2)
    if v3-v2 < 0:
        return (3*v2+v3)*N/(8*v2+4*v3)
    return N/3

##########################
#
# From this point on these are extra functions not in Mark his original file
#
# Written by Maarten Derickx
#
##########################


def cusp_sums_from_signature_to_standard_basis(cuspsums,signatures,level):
    """
    Converts a list of cuspsums where each entry is a list of multiplicities with respect
    to the input signature order. To a list of cuspsums where each entry is a list of multiplicities
    with respect to the standard order of the cusps.
    """
    P = Permutation([cusp_number_from_signature(s,level)+1 for s in signatures]).inverse()
    return [permutation_action(P,s) for s in cuspsums]

def diamond_permutation(list,d,level):
    N = len(list)
    return [list[(i*d)%level if (i*d)%level < N else (-i*d)%level] for i in range(N)]   
    
def diamond_orbits_cuspsum(list,level):
    return [diamond_permutation(list,n,level) for n in range(1,level//2+1) if gcd(n,level) == 1]
    
def split_into_diamond_orbits(cuspsums,level):
    """
    Input a list of cuspsums w.r.t. to the standard basis. Output a list of orbits of these cuspsums.  
    """
    csp_grp = conjectural_cuspidal_classgroup(level)
    diamonds = counts([tuple(sorted(set(csp_grp(vector(i)) for i in diamond_orbits_cuspsum(c,level)))) for c in cuspsums])
    cusp_grp_to_cusp_sum = dict((csp_grp(vector(c)),c) for c in cuspsums)
    return [[cusp_grp_to_cusp_sum[k] for k in d[0]] for d in diamonds]
