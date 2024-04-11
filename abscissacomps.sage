#written in sagemath 10.3
#THIS CODE IS NOT READY TO BE USED
#TEST CODE: Purpose is to extend method for computing normalized codomain models
    #to give a method to compute abscissa of isogenies.

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

#this function is under testing
def isogenous_curves_abscissa(j1,j2,ell,m,model=None):
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
    Fm1 = 0
    for u in range(0, m+2):
        Phideriv = derivative(Phi,x,m+1-u,y,u)
        Fm1 += K(binomial(m+1,u))*K(j1prime)^(m+1-u)*K(ell)^u*Phideriv(j1,j2)*t^u
    F_Pprime = derivative(F_P, t,1)
    solutions = [r[0] for r in F_P.roots()]
    out = [[], []]
    for r in solutions:
        Atilde = ( -r^2)/(K(48)*j2*(j2-K(1728))); Btilde = ( -r^3)/(K(864)*j2^2*(j2-K(1728)))
        out[0].append((ell^4*Atilde,ell^6*Btilde))
        logderivdiff = (r*K(m)/K(m+1))*(Fm1(r)/F_Pprime(r))
        abscissa = K(ell/2)*(logderivdiff + K(1/2)*(A^2/B - K(ell)*Atilde^2/Btilde) -K(2/3)*(B/A - K(ell)*Btilde/Atilde))
        out[1].append(abscissa)
    return out

#this function is under testing 
def abscissadouble(j1,j2,ell,m, model=None):
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
    Fm1 = 0
    for u in range(0, m+2):
        Phideriv = derivative(Phi,x,m+1-u,y,u)
        Fm1 += K(binomial(m+1,u))*K(j1prime)^(m+1-u)*K(ell)^u*Phideriv(j1,j2)*t^u
    solutions = [r[0] for r in F_P.roots()]
    out = [[], []]
    mm = K(18)*(B/A); kk = j1prime/(K(1728)-j1);
    for r in solutions:
        mmtilde = r/j2; kktilde = r/(K(1728)-j2)
        Atilde = ( -r^2)/(K(48)*j2*(j2-K(1728))); Btilde = ( -r^3)/(K(864)*j2^2*(j2-K(1728)))
        out[0].append((ell^4*Atilde,ell^6*Btilde))
        Phiy2 = derivative(Phi, x, 0, y, 2); Phixy = derivative(Phi,x,1,y,1)
        logderivdiff = K(-1/3)*(Fm1(r)/(K(ell)^2*r^2*Phiy2(j1,j2)+K(ell)*j1prime*r*Phixy(j1,j2)))
        p1 = K(ell)*(logderivdiff/2 + (kk-K(ell)*kktilde)/4 + (K(ell)*mmtilde - mm)/3)
        out[1].append(p1)
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
