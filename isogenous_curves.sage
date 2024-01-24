#Written in SageMath 10.2, Release Date 2023-12-03
def isogenous_curves(j1,j2,ell,m, model=None):
    """Intake a singular point (j1,j2) on Y_0(ell) with multiplicity m and a model E_{A,B} for E_{j1}, return all models for E_{j2} admitting normalized isogenies from E_{A, B}"""
    K = j1.parent(); R.<x,y> = K[]
    if model == None:
        E1 = EllipticCurve_from_j(j1)
    else:
        E1 = EllipticCurve(K, [K(model[0]), K(model[1])])
    A = E1.a4(); B = E1.a6()
    j1prime = 18*(B/A)*j1
    DB = ClassicalModularPolynomialDatabase(); Phi = DB[ell]; Phi = R(Phi)
    #First we need to compute the fiber polynomial for (j1,j2)
    Rt.<t> = K[]; F_P = 0
    for u in range(0,m+1):
        Phideriv = derivative(Phi, x, m-u, y, u)
        F_P += K(binomial(m,u))*K(j1prime)^(m-u)*K(ell)^u*Phideriv(j1,j2)*t^u
    solutions = [r[0] for r in F_P.roots()]
    out = []
    for r in solutions:
        Atilde = ( -K(ell^4) * r^2)/(K(48)*j2*(j2-K(1728))); Btilde = ( -K(ell^6) * r^3)/(K(864)*j2^2*(j2-K(1728)))
        out.append((Atilde,Btilde))
    return out

def example1():
    p=137; ell = 5; K.<w> = GF(p^2); j1 = 43*w+130; j2 = 120*w+100; E = EllipticCurve_from_j(j1); #double from j1 to j2
    out = isogenous_curves(j1,j2,ell,2); models = [EllipticCurve(K, [model[0], model[1]]) for model in out]; maps = [EllipticCurveIsogeny(E, None, Etilde, ell) for Etilde in models];
    for isog in maps:
        print(isog.codomain())
        print(isog.is_normalized())

def example2():
    p=137; ell = 7; K.<w> = GF(p^2); j1 = 43*w+130; j2 = 94*w+114; E = EllipticCurve_from_j(j1) #triple from j1 to j2
    out = isogenous_curves(j1,j2,ell,3); models = [EllipticCurve(K, [model[0], model[1]]) for model in out]; maps = [EllipticCurveIsogeny(E, None, Etilde, ell) for Etilde in models];
    for isog in maps:
        print(isog.codomain())
        print(isog.is_normalized())

def example3():
    p=137; ell = 7; K.<w> = GF(p^2); j1 = 120*w+100; j2 = 17*w+135; E = EllipticCurve_from_j(j1) #quadruple from j1 to j2
    out = isogenous_curves(j1,j2,ell,4); models = [EllipticCurve(K, [model[0], model[1]]) for model in out]; maps = [EllipticCurveIsogeny(E, None, Etilde, ell) for Etilde in models];
    for isog in maps:
        print(isog.codomain())
        print(isog.is_normalized())



