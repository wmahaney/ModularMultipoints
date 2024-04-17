#Written in SageMath 10.2
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
        Phideriv = derivative(Phi, x, u, y, m-u)
        F_P += K(binomial(m,u))*K(j1prime)^(u)*K(ell)^(m-u)*Phideriv(j1,j2)*t^(m-u)
    solutions = [r[0] for r in F_P.roots()]
    models = []
    for r in solutions:
        Atilde = ( -K(ell^4) * r^2)/(K(48)*j2*(j2-K(1728))); Btilde = ( -K(ell^6) * r^3)/(K(864)*j2^2*(j2-K(1728)))
        out.append((Atilde,Btilde))
    return out

def abssingle(j1,j2,ell):
    #code to test abscissa formula on smooth points
    K = j1.parent(); R.<x,y> = K[]
    E1 = EllipticCurve_from_j(j1)
    A = E1.a4(); B = E1.a6()
    j1prime = 18*(B/A)*j1
    DB = ClassicalModularPolynomialDatabase(); Phi = DB[ell]; Phi = R(Phi)
    Phix = derivative(Phi,x,1,y,0); Phiy = derivative(Phi,x,0,y,1);
    j2prime = (-j1prime*Phix(j1,j2))/(K(ell)*Phiy(j1,j2))
    Atilde = (  -j2prime^2)/(K(48)*j2*(j2-K(1728))); Btilde = ( -j2prime^3)/(K(864)*j2^2*(j2-K(1728)))
    model = [K(ell)^4*Atilde, K(ell)^6*Btilde]
    Phix2 = derivative(Phi,x,2,y,0); Phixy = derivative(Phi,x,1,y,1); Phiy2 = derivative(Phi,x,0,y,2);
    #get the difference of log derivatives j''/j' - ell jt''/jt'
    logderivdiff = (j1prime^2*Phix2(j1,j2)+K(2*ell)*j1prime*j2prime*Phixy(j1,j2)*K(ell)^2*j2prime^2*Phiy2(j1,j2))/(K(ell)*j2prime*Phiy(j1,j2))
    G4 = -48*A; G6 = 864*B; G4tilde = -48*Atilde; G6tilde = 864*Btilde;
    abscissa = K(ell/2)*(logderivdiff + K(1/2)*(G4^2/G6 - K(ell/2)*G4tilde^2/G6tilde) + K(2/3)*(G6/G4 - K(ell)*G6tilde/G4tilde )  )
    return model, abscissa

def singleabstest():
    p = 137; ell = 2; K.<w> = GF(p^2); j1 = K(136); j2 = K(78); model,abscissa = abssingle(j1,j2, ell);
    E = EllipticCurve_from_j(j1); Etilde = EllipticCurve(K, [model[0], model[1]]); phi = EllipticCurveIsogeny(E, None, Etilde, ell);
    print(phi.kernel_polynomial())
    print(-K(1/2)*abscissa)


def absdouble(j1,j2,ell,m,model=None):
    """Intake a double point (j1,j2) of Phi and returns the normalized codomain models for the isogenies j1 to j2 as well as the abscissa of the isogenies."""
    #old code for the models
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
        Phideriv = derivative(Phi, x, u, y, m-u)
        F_P += K(binomial(m,u))*K(j1prime)^(u)*K(ell)^(m-u)*Phideriv(j1,j2)*t^(m-u)
    solutions = [r[0] for r in F_P.roots()]
    models = []
    for r in solutions:
        Atilde = ( -r^2)/(K(48)*j2*(j2-K(1728))); Btilde = ( -r^3)/(K(864)*j2^2*(j2-K(1728)))
        models.append((K(ell)^4*Atilde,K(ell)^6*Btilde))
    #new code
    Phiy3 = derivative(Phi,x,0,y,3); Phixy2 = derivative(Phi,x,1,y,2); Phix2y = derivative(Phi,x,2,y,1); Phix3 = derivative(Phi,x,3,y,0);
    Phiy2 = derivative(Phi,x,0,y,2); Phixy = derivative(Phi,x,1,y,1); Phix2 = derivative(Phi,x,2,y,0);
    abscissas = []
    for r in solutions:
        Atilde = ( -r^2)/(K(48)*j2*(j2-K(1728))); Btilde = ( -r^3)/(K(864)*j2^2*(j2-K(1728)))
        #compute the 3rd fiber function evaluated at j1,r
        F3 = j1prime^3*Phix3(j1,j2)+K(3)*K(ell)*j1prime^2*r*Phix2y(j1,j2)+K(3)*K(ell)^2*j1prime*r^2*Phixy2(j1,j2)+K(ell)^3*r^3*Phiy3(j1,j2)
        #compute the difference j''/j' - ell jt''/jt'
        logderivdiff = -F3/(K(3)*Phiy2(j1,j2) + K(3)*ell*j1prime*r*Phixy(j1,j2))
        #compute the abscissa
        G4 = K(-48)*A; G6 = K(864)*B; G4tilde = K(-48)*Atilde; G6tilde = K(864)*Btilde
        absc = K(ell/2)*logderivdiff+K(ell/4)*(G4^2/G6 - K(ell)*G4tilde^2/G6tilde) + K(ell/3)*(G6/G4 - K(ell)*G6tilde/G4tilde);
        abscissas.append(absc)
    return models, abscissas

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

