def elliptic_curve_from_hoeij_data(line):
    """Given a line of the file "http://www.math.fsu.edu/~hoeij/files/X1N/LowDegreePlaces" 
    that is actually corresponding to an elliptic curve, this function returns the triple
    N,E,d where E is an elliptic curve over a numberfield of degree d on which (0,0) is a 
    point of order N.
    """
    Rx.<x>=QQ[]
    Rxy.<y>=Rx[]
    
    N=ZZ(line.split(",")[0].split()[-1])
    x_rel=Rx(line.split(',')[-2][2:-4])
    assert x_rel.leading_coefficient()==1
    y_rel=line.split(',')[-1][1:-5]
    K.<x>=QQ.extension(x_rel)
    y_rel=Rxy(y_rel).change_ring(K)
    y_rel=y_rel/y_rel.leading_coefficient()
    if y_rel.degree()==1:
        y = - y_rel[0]
    else:
        #print "needing an extension!!!!"
        L.<y>=K.extension(y_rel)
        K=L
    #B=L.absolute_field('s')
    #f1,f2 = B.structure()
    #x,y=f2(x),f2(y)
    r = (x^2*y-x*y+y-1)/x/(x*y-1)
    s = (x*y-y+1)/x/y
    b = r*s*(r-1)
    c = s*(r-1)
    E=EllipticCurve([1-c,-b,-b,0,0])
    return N,E,K.absolute_degree()    
