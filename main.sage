import os 
import sys

main_dir = os.path.dirname(os.path.abspath("main.sage"))
root_file = os.path.join(main_dir, "PROJECT_ROOT")

with open(root_file, 'r') as f:
    PROJECT_ROOT= f.read().strip()

load(os.path.join(PROJECT_ROOT, "elkies_isogeny.sage"))

def multipoint_isogeny(j1, j2, l, mult=None, model=None):
    """
    Parameters
    ----------
    j1 : Element of a field k
        The j-invariant of the first elliptic curve E1.
    j2 : Element of a finite field k
        The j-invariant of the second elliptic curve E2.
    l : int
        A prime integer coprime to the characteristic of the field
    model : (A, B) tuple of Elements of a finite field or None.
        Optional parameter that lets the user give E1 a specific short Weierstrass model. 

    Returns
    -------
    E1 : EllipticCurve
        The first elliptic curve with j-invariant j1.
    kernel_polynomials : list of kernel polynomials for all the l-isogenies from E1 to elliptic curves with j-invariant j2

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
        if p < 2*l:
            raise ValueError('characteristic must be at least 2*degree')

    if model == None:
        E = EllipticCurve_from_j(j1)
        E = E.short_weierstrass_model()
    else:
        #model is expected as (A, B) such that y^2 = x^3 + Ax + B has j-invariant j1
        E = EllipticCurve(k, model) 
        """
        This needs to be able to handle intaking an EllipticCurve object as well 
        """
        

    #proceed into main body 

    A=E.a4()
    B=E.a6() 
    R.<X,Y>=k[]
    Rt.<t>=k[]
    Phi=classical_modular_polynomial(l)
    Phi=R(Phi)

    #parameters for abscissa formula later
    mu1 = 18*B/A
    j1prime = mu1*j1 
    nu1 = mu1*j1/(1728-j1)


    #compute derivatives evaluated at (j1, j2)
    derivative_evals = {}
    #This can be sped up instantiated modular polynomials.
    if mult:
        for u in range(0, mult+2):
            for v in range(0, mult+2):
                Phi_deriv = derivative(Phi, X, u, Y, v)
                Phi_deriv  = Phi_deriv(j1, j2)
                derivative_evals[(u,v)] = Phi_deriv

    else:
        #we have to figure out the multiplicity
        mult = 1
        deriv = derivative(Phi, X, mult, Y, 0)
        deriv_eval = deriv(j1, j2) 

        while deriv_eval == 0:
            mult +=1 
            deriv = derivative(Phi, X, mult, Y, 0)
            deriv_eval = deriv(j1, j2)
        #now compute the evaluations like before 
        for u in range(0, mult+2):
            for v in range(0, mult+2):
                Phi_deriv = derivative(Phi, X, u, Y, v)
                Phi_deriv  = Phi_deriv(j1, j2)
                derivative_evals[(u,v)] = Phi_deriv

    #Compute fiber polynomial
    fiber_poly =0
    for u in range(0, mult+1):
        fiber_poly += binomial(mult,u) * l^(mult-u) * j1prime^u * derivative_evals[(u, mult-u)] * t^(mult-u)

    roots = [r[0] for r in fiber_poly.roots()] 
    isogeny_data={}
    isogenies=[]
    for r in roots:
        datum={} 
        Atilde = ( -l^4 * r^2)/(48*j2*(j2-1728))
        Btilde = (-l^6 * r^3)/(864*j2^2*(j2-1728))
        model= [Atilde, Btilde]
        Etilde_r = EllipticCurve(k, model)
        datum['codomain_curve']=Etilde_r 

        fiber_function_multplus1 = 0
        for u in range(0, mult+2):
            fiber_function_multplus1 += binomial(mult+1, u) * l^(mult+1-u) * j1prime^u * derivative_evals[(u, mult+1-u)] * r^(mult+1-u)
        #now we take the derivative of the fiber polynomial evaluated at our fixed model Atilde, Btilde
        fiber_poly_ddt = derivative(fiber_poly, t, 1); fiber_poly_ddt = fiber_poly_ddt(r);
        if fiber_poly_ddt == 0:
        #If we found a root corresponding to multiple models we need to abort now because our later abscissa formula
            #will involve division by 0
            """
            In this case our formula for the abscissa fails but one can still use ELlipticCurveIsogeny(E, None, Etilde_r, l) to compute the isogeny. We make the abscissa None to indicate something went wrong
            """
            phi = EllipticCurveIsogeny(E, None, Etilde_r, l)
            datum['isogeny_abscissa'] = None
            datum['isogeny_kernel_polynomial'] = phi.kernel_polynomial()
            isogeny_data[r] = datum
            continue
        #now compute the difference of logarithmic derivatives given by our formula
        log_deriv_diff = mult/r * fiber_function_multplus1/( binomial(mult+1, mult-1) * fiber_poly_ddt )
        #now use Drew's formula to recover the absicca
        mu2 = r/j2; nu2 = r/(1728-j2) 
        # Abscissa formula of Sutherland from Steven Galbraith's ``Mathematics of Public Key Cryptography,''
        #Chapter 25, pg555 https://www.math.auckland.ac.nz/~sgal018/crypto-book/ch25.pdf
        abscissa = l * (log_deriv_diff/2 + (nu1 - l*nu2)/4 + (l*mu2 - mu1)/3 )

        #now use Fast Elkies to compute the kernel polynomial
        kernel_poly = fast_elkies(E, Etilde_r, l, abscissa)
        datum['abscissa'] = abscissa
        datum['kernel_polynomial'] = str(kernel_poly)
        phi = EllipticCurveIsogeny(E, kernel_poly)

        isogeny_data[r]=datum
        isogenies.append(phi)

    #we will return the full dictionary of isogeny data and also a list of isogenies
    return isogeny_data, isogenies