import os 
import sys

"""
Main function
"""

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
        If None, the algorithm will compute a short Weierstrass model for E1 itself.
    mult : int or None
        Optional parameter that lets the user give the multiplicity of the root (j1, j2) in the modular polynomial Phi_l(X,Y).
        If None, the algorithm will compute the multiplicity itself.

    Returns
    -------
    isogeny_data : dict
        A dictionary where each key is a root r of the fiber polynomial and each value is another dictionary with the following keys
        'domain_curve': The elliptic curve E1 with j-invariant j1.
        'codomain_curve': The elliptic curve E2 with j-invariant j2 corresponding to the root r.
        'sigma': The sum of the x-coordinates of the affine points in the kernel of the isogeny
        'kernel_polynomial': The kernel polynomial of the isogeny as a string.
        'isogeny': The isogeny from E1 to E2 as an EllipticCurveIsogeny object.
            If the algorithm fails to compute the isogeny for a particular root r (due to division by zero in the sigma formula), the value for 'isogeny', 'sigma', and 'kernel_polynomial' will be None.

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

    #parameters for sigma formula later
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
        datum['domain_curve'] = E 
        datum['codomain_curve']=Etilde_r 

        fiber_function_multplus1 = 0
        for u in range(0, mult+2):
            fiber_function_multplus1 += binomial(mult+1, u) * l^(mult+1-u) * j1prime^u * derivative_evals[(u, mult+1-u)] * r^(mult+1-u)
        #now we take the derivative of the fiber polynomial evaluated at our fixed model Atilde, Btilde
        fiber_poly_ddt = derivative(fiber_poly, t, 1); fiber_poly_ddt = fiber_poly_ddt(r);
        if fiber_poly_ddt == 0:
        #If we found a root corresponding to multiple models we need to abort now because our later sigma formula
            #will involve division by 0
            """
            In this case our formula for sigma fails. We return a None type to indicate algorithm failure without quitting out completely
            """
            datum['sigma'] = None
            datum['isogeny_kernel_polynomial'] = None
            isogeny_data[r] = datum
            continue
        #now compute the difference of logarithmic derivatives given by our formula
        log_deriv_diff = mult/r * fiber_function_multplus1/( binomial(mult+1, mult-1) * fiber_poly_ddt )
        #now use Drew's formula to recover the absicca
        mu2 = r/j2; nu2 = r/(1728-j2) 
        # sigma formula of Sutherland from Steven Galbraith's ``Mathematics of Public Key Cryptography,''
        #Chapter 25, pg555 https://www.math.auckland.ac.nz/~sgal018/crypto-book/ch25.pdf
        sigma = l * (log_deriv_diff/2 + (nu1 - l*nu2)/4 + (l*mu2 - mu1)/3 )

        #now use Fast Elkies to compute the kernel polynomial
        kernel_poly = fast_elkies(E, Etilde_r, l, sigma)
        datum['sigma'] = sigma
        datum['kernel_polynomial'] = str(kernel_poly)
        phi = EllipticCurveIsogeny(E, kernel_poly)
        datum['isogeny']=phi

        isogeny_data[r]=datum

    #we will return the full dictionary of isogeny data and also a list of isogenies
    return isogeny_data


"""
Fast Elkies Code
"""
def fast_elkies(E1, E2, l, sigma):
    """
    Compute the kernel polynomial of the unique normalized isogeny
    of degree ``l`` between ``E1`` and ``E2``.

    Both curves must be given in short Weierstrass form, and the
    characteristic must be either `0` or no smaller than `4l+4`.

    l must be odd.

    sigma is the sum of the x-coordinates of the affine points in the kernel. 

    ALGORITHM: [BMSS2006]_, algorithm *fastElkies*.

    EXAMPLES:: TODO

        sage: from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_bmss
        sage: E1 = EllipticCurve(GF(167), [153, 112])
        sage: E2 = EllipticCurve(GF(167), [56, 40])
        sage: compute_isogeny_bmss(E1, E2, 13)
        x^6 + 139*x^5 + 73*x^4 + 139*x^3 + 120*x^2 + 88*x
    """
    # Original author of this function: Rémy Oudompheng.
    # https://github.com/remyoudompheng/isogeny_weber/blob/64289127a337ac1bf258b711e02fea02b7df5275/isogeny_weber/isogenies.py#L272-L332
    # Released under the MIT license: https://github.com/remyoudompheng/isogeny_weber/blob/64289127a337ac1bf258b711e02fea02b7df5275/LICENSE
    # Slightly adjusted for inclusion in the Sage library.
    # Updated to run FastElkies rather than FastElkies' for odd ell. 
    if l % 2 == 0:
        raise ValueError('l must be odd')
    if E1.a1() or E1.a2() or E1.a3():
        raise ValueError('E1 must be a short Weierstrass curve')
    if E2.a1() or E2.a2() or E2.a3():
        raise ValueError('E2 must be a short Weierstrass curve')
    char = E1.base_ring().characteristic()
    if char != 0 and char < 2*l:
        raise ValueError('characteristic must be at least 2*degree')
    Rx, x = E1.base_ring()["x"].objgen()
    # Compute C = 1/(1 + Ax^4 + Bx^6) mod x^4l
    A, B = E1.a4(), E1.a6()
    C = (1 + A * x**4 + B * x**6).inverse_series_trunc(l)
    # Solve differential equation
    # The number of terms doubles at each iteration.
    # S'^2 = G(x,S) = (1 + A2 S^4 + B2 S^6) / (1 + Ax^4 + Bx^6)
    # S = x + O(x^2)
    A2, B2 = E2.a4(), E2.a6()
    S = x + (A2 - A) / 10 * x**5 + (B2 - B) / 14 * x**7
    sprec = 8
    while sprec < l + 1:
        assert sprec % 2 == 0
        sprec = min(sprec, (l + 1) //2)
        # s1 => s1 + x^k s2
        # 2 s1' s2' - dG/dS(x, s1) s2 = G(x, s1) - s1'2
        s1 = S
        ds1 = s1.derivative() 
        s1pows = [1, s1]
        while len(s1pows) < 7:
            s1pows.append(s1._mul_trunc_(s1pows[-1], 2 * sprec))
        GS = C * (1 + A2 * s1pows[4] + B2 * s1pows[6])
        dGS = C * (4 * A2 * s1pows[3] + 6 * B2 * s1pows[5])
        # s2' = (dGS / 2s1') s2 + (G(x, s1) - s1'2) / (2s1')
        denom = (2 * ds1).inverse_series_trunc(2 * sprec)
        a = dGS._mul_trunc_(denom, 2 * sprec)
        b = (GS - ds1**2)._mul_trunc_(denom, 2 * sprec)
        s2 = a.add_bigoh(2 * sprec).solve_linear_de(prec=2 * sprec, b=b, f0=0)
        S = s1 + Rx(s2)
        sprec = 2 * sprec
    # Reconstruct:
    # S = x * T(x^2)
    # Compute U = 1/T^2
    # Reconstruct N(1/x) / D(1/x) = U
    T = Rx([S[2 * i + 1] for i in range((l+1)//2)])
    U = T._mul_trunc_(T, (l+1)//2).inverse_series_trunc((l+1)//2)
    h = U.list()[1:]
    d = (l-1)//2
    q = [d, sigma / 2]
    for i in range(1,(l-1)//2):
        q.append(h[i]/(4*i+2) - (2*i-1)*A*q[i-1]/(2*i+1) - (2*i-2)*B*q[i-2]/(2*i+1))
        # q.append((h[i] - (4*i-2)*A*q[i-1] - (4*i-4)*B*q[i-2])/(4*i+2))
    g = Rx([0] + [-q[i] / i for i in range(1,d+1)])
    # g = (g.add_bigoh(d+1).exp().truncate(d+1)).reverse()
    # the above should work but seems to fail for larger degree isogenies.
    g = poly_exp(g, d+1).reverse()
    return g

def reciprocal(f, n):
    R = f.parent()
    h = R(f.coefficients()[0]^(-1))
    i = 1
    while i < n:
        i = 2*i
        h = (h*(2-f*h)).truncate(i)
    return h.truncate(n)

def poly_exp(f,n):
    R = f.parent()
    if f == 0:
        return R(1)
    g = R(1)
    i = 1
    while i < n:
        i = 2*i
        h = 1 + f - poly_log(g, i)
        g = g.multiplication_trunc(h, i)
    return g.truncate(n)


# Calculats the logarithm of f to n terms
#
# Inputs:
#     f: An element of F[x] satisfying f(0)=1, where F is the finite field of size p^2
#     n: A positive integer such that n < p
# Outputs:
#     log_f: The logarithm of f in F[[x]] truncated to n terms
#
# Runs in O(n M(log p) llog p + M(n log p))
#
# Source [2]
def poly_log(f, n):
    R = f.parent()
    x = R.gen()

    f_prime = f.derivative()
    a = reciprocal(f, n-1)
    log_f = f_prime.multiplication_trunc(a, n-1)
    log_f = log_f.integral()
    return log_f