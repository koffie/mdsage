def intersection(self,other):
    """Returns the intersection of two quotient modules self = V1/W and other V2/W
    V1 and V2 should be submodulus of the same ambient module.
    
    EXAMPLE::
        
        sage: W = (ZZ^3).scale(100)
        sage: V1 = (ZZ^3).span([[5,5,0],[0,10,20],[0,0,50]]) + W
        sage: V2 = (ZZ^3).span([[5,0,0],[0,25,50]]) + W
        sage: V1/W;V2/W
        Finitely generated module V/W over Integer Ring with invariants (2, 10, 20)
        Finitely generated module V/W over Integer Ring with invariants (4, 20)
        sage: intersection(V1/W,V2/W)
        Finitely generated module V/W over Integer Ring with invariants (2, 4)
    
    """
    assert self.W()==other.W()
    return self.V().intersection(other.V()).quotient(self.W())
