def get_ABsigma(E, j_tilde, ell):
    """
    Computes an equation for a curve Etilde with j-invariant jtilde that is ell-isogenous to E by a normalized ell-isogeny, along with the sum of the x-coordinates of the affine points in the kernel of the isogeny. 
    
    E must be in short Weierstrass form. We also require that j(E), j_tilde are not in [0,1728] and (j,jtilde) is a nonsingular point of the ellth classical modular polynomial. """
    j = E.j_invariant()
    A = E.a_invariants()[3]
    B = E.a_invariants()[4]
    F = E.base_field()
    FXY.<X,Y> = PolynomialRing(F)

    # Compute all deriviatives and partial derivatives of the ell'th modular polynomial at (j, j_tilde)
    # Takes O(ell^2 M(log p))
    Phi_ell = classical_modular_polynomial(ell)
    Phi_ell = FXY(Phi_ell)
    Phi_X = Phi_ell.derivative(X)
    Phi_Y = Phi_ell.derivative(Y)
    Phi_XX = Phi_X.derivative(X)(j,j_tilde)
    Phi_YY = Phi_Y.derivative(Y)(j,j_tilde)
    Phi_XY = Phi_X.derivative(Y)(j,j_tilde)
    Phi_X = Phi_X(j,j_tilde)
    Phi_Y = Phi_Y(j,j_tilde)

    # Recover the coefficients of E~
    # Takes O(M(log p) llog p)
    m = 18*B/A
    j_prime = m*j
    k = j_prime / (1728-j)
    j_tilde_prime = -j_prime * Phi_X / (ell*Phi_Y)
    m_tilde = j_tilde_prime / j_tilde
    k_tilde = j_tilde_prime / (1728 - j_tilde)
    A_tilde = ell^4 * m_tilde * k_tilde / 48
    B_tilde = ell^6 * m_tilde^2 * k_tilde / 864

    # Recover the sum of the abscissas of the kernel
    # Takes O(M(log p) llog p)
    r = -(j_prime^2*Phi_XX
          + 2*ell*j_prime*j_tilde_prime*Phi_XY
          + ell^2*j_tilde_prime^2*Phi_YY) / (j_prime*Phi_X)
    p1 = ell*(r/2 + (k-ell*k_tilde)/4 + (ell*m_tilde - m)/3)

    return A_tilde, B_tilde, p1
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