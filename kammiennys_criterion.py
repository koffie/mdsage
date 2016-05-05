"""
Kamienny's Criterion

AUTHOR:

    - William Stein, February 2010
    - Maarten Derickx, 2010-2012
"""

Kamienny_Version = "1.4"

import os

from sage.all import *
from sqlalchemy import Column,Table,create_engine,MetaData
from sqlalchemy import Integer as SQLInteger
from sqlalchemy import String as SQLString
from sqlalchemy import Boolean as SQLBoolean
from sqlalchemy import select

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.structure.dynamic_class import dynamic_class


def verify_result(i,result_dir="kamienny_run",verbose=True):
    if i['result_type']=='single' and i['satisfied']:
        torsion,congruence=i['torsion_order'],i['congruence_type']
        C = KamiennyCriterion(i['torsion_order'],i['congruence_type'])
        v=None
        if i['use_rand_vec']:
            result_filename=os.path.join(result_dir,"result_%s_%s"%(torsion,congurence))
            with open(result_filename) as file:
                v=loads(file.read())
        satisfied,message,dependencies = C.verify_criterion(i['degree'],n=i['n'],
                q=i['q'],p=i['p'],v=v,use_rand_vec=i['use_rand_vec'])
        if verbose and not satisfied:
            print "The following result is corrupt: %s"%(i)
        return satisfied
    return True


def get_results_from_dir(prime,congurence_type,result_dir="kamienny_run",soft_fail=True):
    result_filename=os.path.join(result_dir,"result_%s_%s"%(prime,congurence_type))
    return get_results_from_file(result_filename)

def get_results_from_file(filename,soft_fail=True):
    try:
        with open(filename,"r") as file:
            result_str="".join(file)
    except IOError as error:
        if soft_fail:
            return []
        else:
            raise error
    exec("result="+result_str)
    return result

def get_all_results_in_range(begin,end,result_dir="kamienny_run"):
    primes=prime_range(begin,end)
    l=[]
    for prime in primes:
        for congruence_type in (0,1):
            l.extend(get_results_from_dir(prime,congruence_type,result_dir=result_dir))
    return l

def get_missing_primes_in_range_for_deg(prime_range,d,result_list,filter=lambda i:True):
    prime_range=set(prime_range)
    verified_primes=[i['torsion_order'] for i in result_list if i['satisfied'] and i['degree']==d and filter(i)]
    return prime_range.difference(verified_primes)



def load_result_list(first_prime,last_prime,result_dir="kamienny_run"):
    l=[]
    for i in prime_range(first_prime,last_prime):
        l.append(get_results_from_dir(i,0,result_dir=result_dir))
        l.append(get_results_from_dir(i,1,result_dir=result_dir))
    l=sum(l,[])
    return l



def open_database(location, verbose=False):
    """
    Creates an sqlite database at location, ready for storing data computed by this criterion, or if it already exists it opens it. 
    """
    db = create_engine('sqlite:///' + location)

    db.echo = verbose

    metadata = MetaData(db)

    results = Table('results', metadata,
        Column('run_id', SQLInteger, primary_key=True),
        Column('torsion_order', SQLInteger, index=True),
        Column('extension_degree', SQLInteger, index=True),
        Column('congruence_type', SQLBoolean, index=True),
        Column('t1n', SQLInteger),
        Column('t1mod', SQLInteger),
        Column('t2q', SQLInteger),
        Column('result', SQLBoolean, index=True),
        Column('kamienny_version', SQLString(10)),
        Column('sage_version', SQLString(80))
    )
    """
    t1data = Table('t1data', metadata,
        Column('run_id', SQLInteger, primary_key=True),
        Column('torsion_order', SQLInteger, index=True),
        Column('congruence_type', SQLBoolean, index=True),
        Column('t1n', SQLInteger,index=True),
        Column('t1mod', SQLInteger),
        Column('pol_degree', SQLInteger),
        Column('largest_d', SQLInteger),   
    )
    
    t2data = Table('t2data', metadata,
        Column('run_id', SQLInteger, primary_key=True),
        Column('torsion_order', SQLInteger,index=True),
        Column('congruence_type', SQLBoolean, index=True),
        Column('t2q', SQLInteger,index=True),
        Column('largest_d', SQLInteger),
    )
    """
    if not db.has_table('results'):
        results.create()
    """    
    if not db.has_table('t1data'):
        t1data.create()
    if not db.has_table('t2data'):
        t2data.create()
    """
        
    return db,results #,t1data,t2data


def run_criterion(result_table,iterator,congruence_type,algorithm="default",verbose=False):
    """
    The iterator should spit out things of the form (torsion_order,degree,t1n,t1mod,t2q)
    """
    insert_point=result_table.insert()
    old_torsion=None
    for torsion_order,degree,t1n,t1mod,t2q in iterator:
        if not old_torsion==torsion_order:
            old_torsion=torsion_order
            C=KamiennyCriterion(torsion_order,congruence_type=congruence_type,algorithm=algorithm,verbose=verbose)
        verified,message=C.verify_criterion(degree,n=t1n,p=t1mod,q=t2q,verbose=verbose)
        insert_point.execute(torsion_order=int(torsion_order), extension_degree=int(degree), 
                             t1n=int(t1n), t1mod=int(t1mod), t2q=int(t2q),
                             congruence_type=int(congruence_type), result=verified, 
                             kamienny_version=Kamienny_Version,sage_version=version())


def run_criterion2(torsion_order,iterator,congruence_type,use_rand_vec=False,algorithm="default",verbose=False,dump_dir="kamienny_run"):
    results=[]
    C=KamiennyCriterion(torsion_order,congruence_type=congruence_type,
                        algorithm=algorithm,verbose=verbose, dump_dir=dump_dir)

    for degree,t1n,t1mod,t2q in iterator:
        satisfied,message=C.verify_criterion(degree,n=t1n,p=t1mod,q=t2q,use_rand_vec=use_rand_vec,verbose=verbose)
        results.append({"torsion_order":torsion_order,"congruence_type":congruence_type,
                        "algorithm":algorithm,"degree":degree,"n":t1n,"p":t1mod,
                        "q":t2q,"satisfied":satisfied,"message":message,"use_rand_vec":use_rand_vec})

    if dump_dir:
        output_file = os.path.join(dump_dir,"result_%s_%s"%(torsion_order,congruence_type))
        with open(output_file,"w") as file:
            file.write(str(results))
    return results

def run_criterion3(torsion_order,degrees,n_min,n_max,q_min,q_max,congruence_type,
                   t1mod=65521,use_rand_vec=False,algorithm="custom",verbose=False,
                   dump_dir="kamienny_run",stop_if_satisfied=True):
    results=[]
    C=KamiennyCriterion(torsion_order,congruence_type=congruence_type,
                        algorithm=algorithm,verbose=verbose, dump_dir=dump_dir)

    for degree in degrees:
        result,dependancies=C.verify_criterion_range(degree,n_min,n_max,q_min,q_max,
                                t1mod,use_rand_vec=use_rand_vec, verbose=verbose,stop_if_satisfied=stop_if_satisfied)
        results.extend(result)
    if dump_dir:
        output_file = os.path.join(dump_dir,"result_%s_%s"%(torsion_order,congruence_type))
        with open(output_file,"w") as file:
            file.write(str(results))
    return results


def criterion_iterator(torsion_order,degrees=range(3,8),t1n_bound=12,t1mod=65521,t2q_bound=14):
    iterator=[]
    for deg in degrees:
        for t1n in range(2,t1n_bound):
            for q in prime_range(3,t2q_bound):
                if q!=torsion_order:
                    iterator.append((deg,t1n,t1mod,q))
    return iterator

def unverified_primes_in_range(result_table,d,start,stop):
    p=set(prime_range(start,stop))
    s=select([result_table.c.torsion_order], (result_table.c.extension_degree == d) & (result_table.c.result== True))
    l=list(p.difference(i[0] for i in s.execute))
    l.sort
    return l

def get_verify_input_for_prime(result_table,prime):
    pass

class KamiennyCriterion:
    """
    This class if for verification of Kamienny's criterion for quartic fields
    using ell=2 for a given p.

    EXAMPLES::

        sage: C = KamiennyCriterion(31); C
        Kamienny's Criterion for p=31
        sage: C.dbd(2)
        26 x 26 dense matrix over Finite Field of size 2
        sage: C.T(2)
        26 x 26 dense matrix over Finite Field of size 2
        sage: C.t1()
        26 x 26 dense matrix over Finite Field of size 2
        sage: C.t2()
        26 x 26 dense matrix over Finite Field of size 2
        sage: C.t()
        26 x 26 dense matrix over Finite Field of size 2
    """
    def __init__(self, p, congruence_type=1, algorithm="custom", verbose=False, dump_dir=None):
        """
        Create a Kamienny criterion object.

        INPUT:

            - `p` -- prime -- verify that there is no order p torsion
              over a degree `d` field
            - ``algorithm`` -- "default" or "custom" whether to use a custom (faster)
              integral structure algorithm or to use the sage builtin algortihm
            - ``verbose`` -- bool; whether to print extra stuff while
              running.

        EXAMPLES::

            sage: C = KamiennyCriterion(29, algorithm="custom", verbose=False); C
            Kamienny's Criterion for p=29
            sage: C.use_custom_algorithm
            True
            sage: C.p
            29
            sage: C.verbose
            False        
        """
        self.verbose = verbose
        self.dump_dir = dump_dir
        if self.verbose: tm = cputime(); mem = get_memory_usage(); print "init"
        if not is_prime(p):
            raise ValueError, "p must be prime"
        self.p = p
        self.congruence_type=congruence_type
        self.algorithm=algorithm
        if congruence_type==0:
            self.congruence_group=Gamma0(p)
        elif congruence_type==1:
            self.congruence_group=Gamma1(p)
        else:
            raise TypeError("congruence_type=%s but should be 0 or 1"%congruence_type)
        self.M = ModularSymbols(self.congruence_group, sign=1)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "modsym"
        self.S = self.M.cuspidal_submodule()
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "cuspsub"
        self.use_custom_algorithm = False
        if algorithm=="custom":
            self.use_custom_algorithm = True
        if self.use_custom_algorithm:
            int_struct = self.integral_cuspidal_subspace()
            if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "custom int_struct"
        else:    
            int_struct = self.S.integral_structure()
            if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "sage int_struct"
        self.S_integral = int_struct
        v = VectorSpace(GF(2), self.S.dimension()).random_element()
        self.v=v
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "rand_vect"
        if dump_dir:
            v.dump(dump_dir+"/vector%s_%s" % (p,congruence_type))
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "dump"



    def __repr__(self):
        """
        Return string representation.

        EXAMPLES::

            sage: KamiennyCriterion(29).__repr__()
            "Kamienny's Criterion for p=29"
        """
        return "Kamienny's Criterion for p=%s" % self.p

    def dbd(self, d):
        """
        Return matrix of <d>.

        INPUT:

            - `d` -- integer

        OUTPUT:

            - a matrix modulo 2

        EXAMPLES::

            sage: C = KamiennyCriterion(29)
            sage: C.dbd(2)
            22 x 22 dense matrix over Finite Field of size 2
            sage: C.dbd(2)[0]
            (0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        if self.verbose: tm = cputime(); mem = get_memory_usage(); print "dbd start"
        try: return self._dbd[d % self.p]
        except AttributeError: pass
        # Find a generator of the integers modulo p:
        z = primitive_root(self.p)
        # Compute corresponding <z> operator on integral cuspidal modular symbols
        
        X = self.M.diamond_bracket_operator(z).matrix()
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "create d"
        X = X.restrict(self.S_integral, check=False)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "restrict d"
        
        X = matrix_modp(X)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "mod d"
        # Take powers to make list self._dbd of all dbd's such that
        # self._dbd[d] = <d>
        v = [None] * self.p
        v[z] = X
        a = Mod(z, self.p)
        Y = X
        for i in range(self.p - 2): 
            Y *= X
            a *= z
            v[a] = Y
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "mul"

        assert v.count(None) == 1
        self._dbd = v
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "bdb finnished"
        return v[d % self.p]

    @cached_method
    def T(self, n):
        """
        Return matrix mod 2 of the n-th Hecke operator on the +1
        quotient of cuspidal modular symbols.

        INPUT:

            - `n` -- integer

        OUTPUT:

            matrix modulo 2

        EXAMPLES::

            sage: C = KamiennyCriterion(29)
            sage: C.T(2)
            22 x 22 dense matrix over Finite Field of size 2
            sage: C.T(2)[0]
            (0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
        """
        if self.verbose: tm = cputime(); mem = get_memory_usage(); print "T(%s) start" % (n)
        
        T = self.M.hecke_matrix(n).restrict(self.S_integral, check=False)
    
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "T created"
        if self.verbose: print "sparsity", len(T.nonzero_positions()) / RR(T.nrows()**2), T.nrows()
        T = matrix_modp(T)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "T reduced"
        #self.M._hecke_matrices={}
        #self.S._hecke_matrices={}
        #if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "T freed"
        return matrix_modp(T)

    def t1(self, n=5):
        """
        Return choice of element t1 of the Hecke algebra mod 2,
        computed using the Hecke operator $T_n$, where n is self.n
        
        INPUT:
        
            - `n` -- integer (optional default=5)
            
        OUTPUT:

            - a mod 2 matrix

        EXAMPLES::

            sage: C = KamiennyCriterion(29)
            sage: C.t1()
            22 x 22 dense matrix over Finite Field of size 2
            sage: C.t1() == 1
            True
            sage: C = KamiennyCriterion(37)
            sage: C.t1()[0]
            (0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1)
        """
        T = self.S.hecke_matrix(n)
        f = T.charpoly()
        F = f.factor()
        if prod(i[1] for i in F) != 1:
            raise ValueError("T_%s needs to be a generator of the hecke algebra"%n)

        # Compute the iterators of T acting on the winding element.
        e = self.M([0, oo]).element().dense_vector()
        t = self.M.hecke_matrix(n).dense_matrix()
        g = t.charpoly()
        Z = t.iterates(e, t.nrows(), rows=True)
        # We find all factors F[i][0] for f such that
        # (g/F[i][0])(t) * e = 0.
        # We do this by computing the polynomial
        #       h = g/F[i][0],
        # turning it into a vector v, and computing
        # the matrix product v * Z.  If the product
        # is 0, then e is killed by h(t).
        J = []
        for i in range(len(F)):
            h, r = g.quo_rem(F[i][0] ** F[i][1])
            assert r == 0
            v = vector(QQ, h.padded_list(t.nrows()))
            if v * Z == 0:
                J.append(i)


        if self.verbose: print "J =", J
        if len(J) == 0:
            # The annihilator of e is the 0 ideal.
            return matrix_modp(identity_matrix(T.nrows()))
            
        # Finally compute t1.  I'm concerned about how
        # long this will take, so we reduce T mod 2 first.

        # It is important to call "self.T(2)" to get the mod-2
        # reduction of T2 with respect to the right basis (e.g., the
        # integral basis in case use_integral_structure is true.
        Tmod2 = self.T(n) 
        g = prod(F[i][0].change_ring(GF(2)) ** F[i][1] for i in J)
        t1 = g(Tmod2)
        return t1

    @cached_method
    def hecke_polynomial(self,n):
        return self.S.hecke_matrix(n).charpoly()

    @cached_method
    def t1_prime(self, n=5, p=65521):
        """
        Return a multiple of element t1 of the Hecke algebra mod 2,
        computed using the Hecke operator $T_n$, where n is self.n.
        To make computation faster we only check if ...==0 mod p.
        Hence J will contain more elements, hence we get a multiple.
        
        INPUT:
        
            - `n` -- integer (optional default=5)
            - `p` -- prime (optional default=65521)

        OUTPUT:

            - a mod 2 matrix

        EXAMPLES::

            sage: C = KamiennyCriterion(29)
            sage: C.t1_prime()
            22 x 22 dense matrix over Finite Field of size 2
            sage: C.t1_prime() == 1
            True
            sage: C = KamiennyCriterion(37)
            sage: C.t1_prime()[0]
            (0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1)
        """
        if self.verbose: tm = cputime(); mem = get_memory_usage(); print "t1 start"
        T = self.S.hecke_matrix(n)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "hecke 1"
        f = self.hecke_polynomial(n) # this is the same as T.charpoly()
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "char 1"
        Fint = f.factor()
        if all(i[1]!=1 for i in Fint):
            return matrix_modp(zero_matrix(T.nrows()))
        #    raise ValueError("T_%s needs to be a generator of the hecke algebra"%n)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "factor 1, Fint = %s"%(Fint)
        R = f.parent().change_ring(GF(p))
        F = Fint.base_change(R)
        # Compute the iterators of T acting on the winding element.
        e = self.M([0, oo]).element().dense_vector().change_ring(GF(p))
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "wind"
        t = matrix_modp(self.M.hecke_matrix(n).dense_matrix(), p)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "hecke 2"
        g = t.charpoly()
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "char 2"
        Z = t.iterates(e, t.nrows(), rows=True)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "iter"
        # We find all factors F[i][0] for f such that
        # (g/F[i][0])(t) * e = 0.
        # We do this by computing the polynomial
        #       h = g/F[i][0],
        # turning it into a vector v, and computing
        # the matrix product v * Z.  If the product
        # is 0, then e is killed by h(t).
        J = []
        for i in range(len(F)):
            if F[i][1]!=1:
                J.append(i)
                continue
            h, r = g.quo_rem(F[i][0] ** F[i][1])
            assert r == 0
            v = vector(GF(p), h.padded_list(t.nrows()))
            if v * Z == 0:
                J.append(i)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "zero check"

        if self.verbose: print "J =", J
        if len(J) == 0:
            # The annihilator of e is the 0 ideal.
            return matrix_modp(identity_matrix(T.nrows()))

        # Finally compute t1.  I'm concerned about how
        # long this will take, so we reduce T mod 2 first.

        # It is important to call "self.T(2)" to get the mod-2
        # reduction of T2 with respect to the right basis (e.g., the
        # integral basis in case use_integral_structure is true.
        Tmod2 = self.T(n)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "hecke mod" 
        g = prod(Fint[i][0].change_ring(GF(2)) ** Fint[i][1] for i in J)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "g has degree %s"%(g.degree())
        t1 = g(Tmod2)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "t1 finnished"
        return t1        



    def t2(self, q=3):
        """
        Return mod 2 matrix t2 computed using the current choice of
        prime q, as returned by self.q.
        
        INPUT:
        
            - `q` -- integer


        OUTPUT:

            - a mod 2 matrix

        EXAMPLES::

            sage: C = KamiennyCriterion(29)
            sage: C.t2()[0]
            (0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1)
        """
        if q==2 or q==self.p:
            raise ValueError("q must be unequal to 2 or self.p")
        if self.congruence_type==1:
            return self.T(q) -  self.dbd(q) - q
        else:
            return self.T(q) - (q  + 1)

    def t(self, n=5, p=65521, q=3):
        """
        Return mod 2 matrix t, using n, p and q as options for t1 and t2.
        This is just the product of self.t1() and self.t2().
        If p=None it uses t1 else it uses t1_prime
        
        INPUT:
        
            - `n` -- integer (optional default=5), used in computing t1
            - `p` -- prime (optional default=46337), used in computing t1
            - `q` -- prime (optional default=3), used in computing t2

        OUTPUT:
                
            - a mod 2 matrix

        EXAMPLES::

            sage: C = KamiennyCriterion(29)
            sage: C.t()[0]
            (0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1)
        """
        t2 = self.t2(q)
        if p == None:
            t1 = self.t1(n)
        else:
            t1 = self.t1_prime(n, p)
        #print "t1,t2"
        #print t1
        #print t2
        #if t1==0:
        #    print "fail"
        return  t1 * t2

    def tdbdTi(self, t, k, v):
        """
        Return a list of lists which contains precomputed values of t*dbd(d)*T(i) with 0<d<p/2 and 0<i<=k 
        
        INPUT:
            
            - `k` -- integer
        
        OUTPUT:
            
            - a list of lists containing hecke operators
            
        EXAMPLES::
            
            sage: C = KamiennyCriterion(29,verbose=False)
            sage: Id=C.t().parent().identity_matrix()
            sage: C.t()==C.tdbdTi(C.t(),2,Id)[0][0]
            True
            sage: C.t()*C.T(2)==C.tdbdTi(C.t(),2,Id)[1][0]
            True
            sage: C.t()*C.dbd(2)*C.T(2)==C.tdbdTi(C.t(),2,Id)[1][1]
            True
        """
        #make sure that _dbd gets initialized
        self.dbd(2)
        return [[t * (d * (self.T(i) * v)) for d in self._dbd[1:self.p // 2 + 1]] for i in xrange(1, k + 1)]
    
    def dependancies(self, t, k, l, v):
        """
        Return the vectorspace of dependancies between the t*dbd(d)*T(i) with 0<d<p/2 and 0<i<=l
        and t*T(l+1),...,t*T(k)
        
        INPUT:
            
            - `k` -- integer the size of the largest partion occuring in the partitions you want to check
            - `l` -- integer which is larger then or equal to the second larges partition
        
        """
        T = self.T
        span = sum(self.tdbdTi(t, l, v), [])
        span += [(t * (T(i) * v)) for i in xrange(l + 1, k + 1)]
        return matrix(GF(2), [i.list() for i in span]).kernel()


    def verify_criterion_range(self,degree,n_min,n_max,q_min,q_max,t1mod,v=None,use_rand_vec=True, verbose=False,stop_if_satisfied=True):
        torsion_order=self.p
        congruence_type=self.congruence_type
        algorithm=self.algorithm
        dependancy_spaces=[]
        results=[]
        for t1n in range(n_min,n_max):
            if self.t1_prime(t1n, t1mod) == 0:
                results.append({"torsion_order":torsion_order,"congruence_type":congruence_type,
                            "algorithm":algorithm,"degree":degree,"n":t1n,
                            "satisfied":False,"message":"","use_rand_vec":use_rand_vec,
                            "result_type":"t1n_is_zero"})
                continue

            for t2q in prime_range(q_min,q_max):
                if t2q==torsion_order:
                    continue
                satisfied,message,dependancies=self.verify_criterion(degree,n=t1n,p=t1mod,q=t2q,use_rand_vec=use_rand_vec,verbose=verbose)
                dependancy_spaces.append(dependancies)
                results.append({"torsion_order":torsion_order,"congruence_type":congruence_type,
                            "algorithm":algorithm,"degree":degree,"n":t1n,"p":t1mod,
                            "q":t2q,"satisfied":satisfied,"message":message,"use_rand_vec":use_rand_vec,
                            "result_type":"single"})
                if stop_if_satisfied and satisfied:
                    break
            if stop_if_satisfied and satisfied:
                break

        print [len(i) for i in dependancy_spaces]
        intersected_dependancy_spaces=[]
        if dependancy_spaces:
            for i in range(len(dependancy_spaces[0])):
                intersection=reduce(lambda x,y:x.intersection(y), [s[i] for s in dependancy_spaces])
                intersected_dependancy_spaces.append(intersection)
            satisfied=all([d.dimension()<13 and (d.dimension()==0 or LinearCodeFromVectorSpace(d).minimum_distance()>degree) 
                                        for d in intersected_dependancy_spaces])
            results.append({"torsion_order":torsion_order,"congruence_type":congruence_type,
                            "algorithm":algorithm,"degree":degree,"n":(n_min,n_max),"p":t1mod,
                            "q":(q_min,q_max),"satisfied":satisfied,"message":"","use_rand_vec":use_rand_vec,
                            "result_type":"range"})
        return results,dependancy_spaces


        

    def verify_criterion(self, d, t=None, n=5, p=65521, q=3, v=None, use_rand_vec=True, verbose=False):
        """
        Attempt to verify the criterion at p using the input t. If t is not
        given compute it using n, p and q

        INPUT:
            
            - `t` -- hecke operator (optional default=None)
            - `n` -- integer (optional default=5), used in computing t1
            - `p` -- prime (optional default=46337), used in computing t1
            - `q` -- prime (optional default=3), used in computing t2
            - `verbose` -- bool (default: True); if true, print to
              stdout as the computation is performed.

        OUTPUT:

            - bool -- True if criterion satisfied; otherwise, False
            - string -- message about what happened or went wrong

        EXAMPLES::

        We can't get p=29 to work no matter what I try, which nicely illustrates
        the various ways the criterion can fail::
        
            sage: C = KamiennyCriterion(29)
            sage: C.verify_criterion(4,n=5)
            J = []
            partition 4=4...
            partition 4=1+3...
            partition 4=2+2...
            partition 4=1+1+2...
            partition 4=1+1+1+1...
            (False, 'Fails with partition 4=1+1+1+1 and d1=2, d2=5, d3=12')
            sage: C = KamiennyCriterion(29)
            sage: C.verify_criterion(4,n=5,q=7)
            J = []
            partition 4=4...
            partition 4=1+3...
            (False, 'Fails with partition 4=1+3 and d=12')
            sage: C = KamiennyCriterion(29)
            sage: C.verify_criterion(4,n=5,q=23)
            J = []
            partition 4=4...
            (False, 'Fails with partition 4=4')

        With p=31 the criterion is satisfied, thus proving that 31 is
        not the order of a torsion point on any elliptic curve over
        a quartic field::

            sage: C = KamiennyCriterion(31); C
            Kamienny's Criterion for p=31
            sage: C.verify_criterion(4,n=5)
            initializing t
            J = []
            ...
            Checking dependancies(4,0)...
            ...no dependancies found
            Checking dependancies(3,1)...
            ...the smallest dependancy has 15 nonzero coefficients.
            3 1 passed
            Checking dependancies(2,2)...
            ...the smallest dependancy has 9 nonzero coefficients.
            2 2 passed
            (True, 'All conditions are satified')

        """
        if self.verbose: tm = cputime(); mem = get_memory_usage(); print "verif start"
        if self.p < (1 + 2 ** (d / 2)) ** 2 and self.verbose:
            print "WARNING: p must be at least (1+2^(d/2))^2 for the criterion to work."
        if t == None:
            t = self.t(n, p, q)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "t"
        if v==None:
            if use_rand_vec:
                v=self.v
            else:
                v=t.parent()(1)
        if self.congruence_type==0:
            verified,message,dependancies=self.verify_gamma0(d,t,v)
        else:
            verified,message,dependancies=self.verify_gamma1(d,t,v)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "total verif time"
        if not verified:
            return verified,message,dependancies
        return True, "All conditions are satified for Gamma%s d=%s p=%s. Using n=%s, modp=%s, q=%s in Kamienny Version %s and %s" % (self.congruence_type,d, self.p, n, p, q, Kamienny_Version, version()),dependancies

    def verify_gamma0(self,d,t,v):
        ker = Matrix(GF(2),[(t * self.T(n) * v).list() for n in xrange(1,d+1)]).kernel()
        #if ker.dimension()<ker.degree():
        #    print ker
        return ker.dimension()==0,"",[ker]
    
    def verify_gamma1(self,d,t,v):
        dependancies=[]
        if self.verbose: tm = cputime(); mem = get_memory_usage(); print "verif gamma1 start"
        satisfied=True
        message=""
        for i in xrange(d, (d - 1) // 2, -1):
        #We only need to start at d/2 since the dependancies(k,l) contains all neccesary
        #checks for largest partition of zise at most k and second largest at most l
            dependancy = self.dependancies(t, i, d - i, v)
            dependancies.append(dependancy)
            if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "dep (%s,%s)" % (i, d - i)
            assert dependancy.degree() == self.p // 2 * (d - i) + 2 * i - d
            assert dependancy.degree() - dependancy.dimension() <= self.S.dimension()
            if dimension(dependancy) == 0:
                if self.verbose: print "...no dependancies found"
            elif dimension(dependancy) > 12:
                satisfied=False
                print "dependancy dimension to large to search trough"
            else:
                if self.verbose: print "dependancy dimension is:", dimension(dependancy)
                min_dist = LinearCodeFromVectorSpace(dependancy).minimum_distance()
                if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "min dist"
                if self.verbose: print "...the smallest dependancy has %s nonzero coefficients." % (min_dist)
                if min_dist > d:
                    if self.verbose: print i, d - i, "passed"
                else:
                    satisfied,message= False, "There is a dependancy of weigt %s in dependancies(%s,%s)" % (min_dist, i, d - i)
        return satisfied, message,dependancies

    def is_independent(self, v):
        """
        Return True if the Hecke operators in v are independent.

        INPUT:

            - `v` -- four elements of the Hecke algebra mod 2 (represented as matrices)

        OUTPUT:

            - bool

        EXAMPLES::

            sage: C = KamiennyCriterion(29)
            sage: C.is_independent([C.T(1), C.T(2), C.T(3), C.T(4)])
            True
            sage: C.is_independent([C.T(1), C.T(2), C.T(3), C.T(1)+C.T(3)])
            False        
        """
#        X = matrix(GF(2), 4, sum([a.list() for a in v], []))
#        c = sage.matrix.matrix_modn_dense.Matrix_modn_dense(X.parent(),X.list(),False,True)
#        return c.rank() == 4

        # This crashes!  See http://trac.sagemath.org/sage_trac/ticket/8301
        return matrix(GF(2), len(v), sum([a.list() for a in v], [])).rank() == len(v)
        raise NotImplementedError
    
    def integral_cuspidal_subspace(self):
        """
        In certatain cases this might return the integral structure of the cuspidal subspace.
        This code is manly a way to compute the integral structe faster than sage does now. 
        It returns None if it cannot find the integral subspace. 
        """
        if self.verbose: tm = cputime(); mem = get_memory_usage(); print "Int struct start"
        #This code is the same as the firs part of self.M.integral_structure
        G = set([i for i, _ in self.M._mod2term])
        G = list(G)
        G.sort()
        #if there is a two term relation between two manin symbols we only need one of the two
        #so that's why we only use elements from G instead of all manin symbols.
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "G"
        B = self.M._manin_gens_to_basis.matrix_from_rows(list(G)).sparse_matrix()
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "B"
        #The collums of B now span self.M.integral_structure as ZZ-module
        B, d = B._clear_denom()
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "Clear denom"
        if d == 1:
            #for explanation see d == 2
            assert len(set([B.nonzero_positions_in_row(i)[0] for i in xrange(B.nrows()) if len(B.nonzero_positions_in_row(i)) == 1 and B[i, B.nonzero_positions_in_row(i)[0]] == 1])) == B.ncols(), "B doesn't contain the Identity"
            if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "Check Id"
            ZZbasis = MatrixSpace(QQ, B.ncols(), sparse=True)(1)
            if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "ZZbasis"
        elif d == 2:
            #in this case the matrix B will contain 2*Id as a minor this allows us to compute the hermite normal form of B in a very efficient way. This will give us the integral basis ZZbasis.
            #if it turns out to be nessecarry this can be generalized to the case d%4==2 if we don't mind to only get the right structure localized at 2
            assert len(set([B.nonzero_positions_in_row(i)[0] for i in xrange(B.nrows()) if len(B.nonzero_positions_in_row(i)) == 1 and B[i, B.nonzero_positions_in_row(i)[0]] == 2])) == B.ncols(), "B doesn't contain 2*Identity"    
            if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "Check 2*Id"
            E = matrix_modp(B,sparse=True)
            if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "matmodp"
            E = E.echelon_form()
            if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "echelon"
            ZZbasis = MatrixSpace(QQ, B.ncols(), sparse=True)(1)
            for (pivot_row, pivot_col) in zip(E.pivot_rows(), E.pivots()):
                for j in E.nonzero_positions_in_row(pivot_row):
                    ZZbasis[pivot_col, j] = QQ(1) / 2 
            if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "ZZbasis"
        else:
            return None
        #now we compute the integral kernel of the boundary map with respect to the integral basis. This will give us the integral cuspidal submodule.                 
        boundary_matrix = self.M.boundary_map().matrix()
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "Boundary matrix"
        ZZboundary_matrix=(ZZbasis*boundary_matrix).change_ring(ZZ)
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "ZZBoundary matrix"
        left_kernel_matrix=ZZboundary_matrix.transpose().dense_matrix()._right_kernel_matrix(algorithm='pari')
        if type(left_kernel_matrix)==tuple:
            left_kernel_matrix=left_kernel_matrix[1]
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "kernel matrix"
        ZZcuspidal_basis=left_kernel_matrix*ZZbasis
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "ZZkernel matrix"
        assert ZZcuspidal_basis.change_ring(QQ).echelon_form()==self.S.basis_matrix() , "the calculated integral basis does not span the right QQ vector space" # a little sanity check. This shows that the colums of ZZcuspidal_basis really span the right QQ vectorspace
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "check"
        #finally create the sub-module, we delibarately do this as a QQ vector space with custom basis, because this is faster then dooing the calculations over ZZ since sage will then use a slow hermite normal form algorithm.
        ambient_module=VectorSpace(QQ,ZZcuspidal_basis.ncols())
        int_struct = ambient_module.submodule_with_basis(ZZcuspidal_basis.rows())
        if self.verbose: print "time and mem", cputime(tm), get_memory_usage(mem), "finnished"
        return int_struct

        
        

def matrix_modp(A, p=2, sparse=False):
    """
    Reduce a matrix mod p (default 2).  We use this function, 
    since there are bugs in the mod 2 matrix reduction in sage. See
    http://trac.sagemath.org/sage_trac/ticket/6904

    INPUT:

        - `A` -- matrix
        - `p` -- prime

    OUTPUT:

        - a matrix over GF(p).

    EXAMPLES::
        sage: matrix_modp(matrix(QQ,2,[1,3,5,2/3]))
        [1 1]
        [1 0]        
    """
    return matrix(GF(p), A.nrows(), A.ncols(), A.list(),sparse=sparse)



