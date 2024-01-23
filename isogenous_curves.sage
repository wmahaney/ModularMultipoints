#Written in SageMath 10.2, Release Date 2023-12-03
def isogenous_curves(j1,j2,ell,m, model=None):
    """Intake a singular point (j1,j2) on Y_0(ell) with multiplicity m and a model E_{A,B} for E_{j1}, return all models for E_{j2} admitting normalized isogenies from E_{A, B}"""
    K = j1.parent(); R.<x,y> = K[]
    if model == None:
        E1 = EllipticCurve_from_j(j1)
    else:
        E1 = EllipticCurve(K, [K(model[0]), K(model[1])])
    A = E1.a4(); B = E1.a6()
    j1prime = K(18)*(B/A)*j1
    DB = ClassicalModularPolynomialDatabase(); Phi = DB[ell]; Phi = R(Phi)
    #First we need to compute the fiber polynomial for (j1,j2)
    Rt.<t> = K[]; F_P = 0
    for u in range(0,m+1):
        Phideriv = derivative(Phi, x, m-u, y, u)
        coeff = K(j1prime)^(m-u)*K(ell)^u*Phideriv(j1,j2)*t^u
        if u == 0 or u == m:
            F_P = F_P + coeff
        else:
            F_P = F_P + K(m)*coeff
    solutions = [r[0] for r in F_P.roots()]
    out = []
    for r in solutions:
        Atilde = ( -K(ell^4) * r^2)/(K(48)*j2*(j2-K(1728))); Btilde = ( -K(ell^6) * r^3)/(K(864)*j2^2*(j2-K(1728)))
        out.append((Atilde,Btilde))
    return out






