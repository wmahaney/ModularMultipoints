"""
Written in Sagemath 10.7
"""

"""
Collecction of functions for computing isogenies between elliptic curves when there are multiple distinct l-isogenies between them. Auxillary function later used in main.sage.
"""


def isogenous_curves_abscissa(j1, j2, l, mult, model = None):
    """
    Parameters
    ----------
    j1 : Element of a field k
        The j-invariant of the first elliptic curve E1.
    j2 : Element of a finite field k
        The j-invariant of the second elliptic curve E2.
    l : int
        A prime integer coprime to the characteristic of the field
    mult : int
        Multiplicity of the singularity
    model : (A, B) tuple of Elements of a finite field or None.
        Optional parameter that lets the user give E1 a specific short Weierstrass model. 

    Returns
    -------

    Raises
    ------
    ValueError
        If l is not prime or not coprime to the characteristic of the field.
        If model=(A,B) is not a valid short Weierstrass model or has the wrong j-invariant.
        If Fq has characteristic 2 or 3.
        If j1 or j2 is 0 or 1728, as this would cause division by zero later.
        
    """
    #k should be a number field, C, or a finite field. 
    k = j1.parent()
    j2 = k(j2)
    p = k.characteristic()
    
    #error handling
    if p == 2 or p == 3:
        raise ValueError(f"Parameter p={p} must be a prime greater than 3.")
    
    if j1 == k(0) or j1 == k(1728) or j2 == k(0) or j2 == k(1728):
        #if j1, j2 = 0, 1728 then our algorithm doesn't work as we would encounter division by zero later
        raise ValueError(f"Parameters j1={j1}, j2={j2} can not be 0 or 1728")
        
    if not is_prime(l):
        #this algorithm should work fine if l is replaced with a composite n coprime to p. 
        raise ValueError(f"Parameter l={l} must be prime.")

    if p > 0: #don't need to check this if field is characteristic zero
        if gcd(l, p) != 1:
            #this algorithm only works for separable isogenies.
            raise ValueError(f"Parameter l={l} must be coprime to the characteristic p={p}.")

    if model == None:
        E1 = EllipticCurve_from_j(j1)
        E1=E1.short_weierstrass_model()
    else:
        #model is expected as (A, B) such that y^2 = x^3 + Ax + B has j-invariant j1
        E1 = EllipticCurve(model[0].parent(), model)
        if E1.a1() or E1.a2() or E1.a3():
            raise ValueError('E1 must be a short Weierstrass curve')
        if E1.j_invariant() != j1:
            raise ValueError(f"Parameter model={model} does not have j-invariant {j1}.")


    
    #load the lth modular polynomial and its partial derivatives
    R.<X,Y>=k[]; #get the polynomial rink k[X,Y] and place Phi_l in that ring
    Phi = classical_modular_polynomial(l); Phi = R(Phi)
    #we compute the matrix of weight mult+1 or less partial derivatives of phil with respect to X, Y
        #evaluated at (j1, j2)
    derivs = matrix(k, mult+2) 
    #This can be sped up with a smarter strategy or using instantiated modular polynomials.
    for u in range(0, mult+2):
        for v in range(0, mult+2):
            Phideriv = derivative(Phi, X, u, Y, v) 
            derivs[u,v] = Phideriv(j1,j2)

    #Computer fiber polynomial
    A=E1.a4(); B=E1.a6() #E1:y^2=x^3+Ax+B
    j1prime = 18*B/A*j1
    S.<t> = k[]
    fiber_poly =0
    for u in range(0, mult+1):
        fiber_poly += binomial(mult,u) * l^(mult-u) * j1prime^u * derivs[u, mult-u] * t^(mult-u)
    #get the roots of the fiber polynomial and then return the associated models and abscissas
    roots = [r[0] for r in fiber_poly.roots()] 
    codomain_data = []
    for r in roots:
        Atilde = ( -l^4 * r^2)/(48*j2*(j2-1728));
        Btilde = (-l^6 * r^3)/(864*j2^2*(j2-1728));
        #get the abscissa
        #first we need the mult+1 fiber function evaluated at r=j~'. 
        fiber_function_multplus1 = 0
        for u in range(0, mult+2):
            fiber_function_multplus1 += binomial(mult+1, u) * l^(mult+1-u) * j1prime^u * derivs[u, mult+1-u] * r^(mult+1-u)
        #now we take the derivative of the fiber polynomial evaluated at our fixed model Atilde, Btilde
        fiber_poly_ddt = derivative(fiber_poly, t, 1); fiber_poly_ddt = fiber_poly_ddt(r);
        if fiber_poly_ddt == 0:
        #If we found a root corresponding to multiple models we need to abort now because our later abscissa formula
            #will involve division by 0
            codomain_data.append([[Atilde, Btilde], None])
            continue
        #now compute the difference of logarithmic derivatives given by our formula
        log_deriv_diff = mult/r * fiber_function_multplus1/( binomial(mult+1, mult-1) * fiber_poly_ddt );
        #now use Drew's formula to recover the absicca
        m1 = 18*B/A; n1 = j1prime/(1728-j1)
        m2 = r/j2; n2 = r/(1728-j2) 
        abscissa = l * (log_deriv_diff/2 + (n1 - l*n2)/4 + (l*m2 - m1)/3 )
        model = [Atilde, Btilde]; 
        codomain_data.append([model, abscissa])
    return E1, codomain_data

