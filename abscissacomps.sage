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
    #requires l>2
    # based on pseudocode of Sutherland from Steven Galbraith's ``Mathematics of Public Key Cryptography,'' Chapter 25, pg555 https://www.math.auckland.ac.nz/~sgal018/crypto-book/ch25.pdf
    K = j1.parent(); R.<x,y> = K[]
    E1 = EllipticCurve_from_j(j1)
    A = E1.a4(); B = E1.a6()
    DB = ClassicalModularPolynomialDatabase(); 
    Phi = DB[ell]; Phi = R(Phi)
    Phix = derivative(Phi,x,1,y,0); Phiy = derivative(Phi,x,0,y,1);
    m = 18*B/A
    j1prime = m*j1
    k = j1prime/(1728-j1)
    j2prime = (-j1prime * Phix(j1,j2)) / (ell * Phiy(j1,j2))
    m2 = j2prime / j2
    k2 = j2prime / (1728 - j2)
    Atilde = ell^4 * m2 * k2 / 48 
    Btilde = ell^6 * m2^2 * k2 / 864
    model = [Atilde, Btilde]
    Phix2 = derivative(Phi,x,2,y,0); Phixy = derivative(Phi,x,1,y,1); Phiy2 = derivative(Phi,x,0,y,2);
    #get the difference of log derivatives j''/j' - ell jt''/jt'
    logderivdiff = -(j1prime^2 * Phix2(j1,j2) + 2 * ell * j1prime * j2prime * Phixy(j1,j2) + ell^2 * j2prime^2 * Phiy2(j1,j2)) / (j1prime * Phix(j1,j2))
    # G4 = -48*A; G6 = 864*B; G4tilde = -48*Atilde; G6tilde = 864*Btilde;
    p1 = ell * (logderivdiff/2 + (k - ell * k2)/4 + (ell * m2 - m)/3)
    return EllipticCurveIsogeny(E1,None,EllipticCurve(model),ell), p1

def absdouble(j1,j2,ell, model=False):
    #requires l>2. Assumes (j1,j2) is a double point of Phi_ell
    """Intakes a double point (j1,j2) of the ellth modular polynomial"""
    K = j1.parent(); R.<x,y> = K[]
    if model:
        E1 = EllipticCurve(K, model)
    else:
        E1 = EllipticCurve_from_j(j1)
    A = E1.a4(); B = E1.a6();
    #Ideally the next few lines could be replaced with a faster loading method
    DB = ClassicalModularPolynomialDatabase(); Phi = DB[ell]; Phi = R(Phi)
    Phix = derivative(Phi,x,1,y,0); Phiy = derivative(Phi,x,0,y,1);
    Phix2 = derivative(Phi,x,2,y,0); Phixy = derivative(Phi,x,1,y,1); Phiy2 = derivative(Phi,x,0,y,2);
    Phix3 = derivative(Phi,x,3,y,0); Phix2y = derivative(Phi,x,2,y,1); Phixy2 = derivative(Phi,x,1,y,2); Phiy3 = derivative(Phi,x,0,y,3);
    #initialize modular parameters
    m1 = 18*B/A #-G6/G4
    j1prime = m1*j1 #-j*G6/G4
    k1 = j1prime/(1728-j1) #G4^2/G6
    #We need to solve for the normalized codomain model
    Kt.<t> = K[]; F2 = j1prime^2*Phix2(j1,j2)+2*ell*j1prime*Phixy(j1,j2)*t+ell^2*Phiy2(j1,j2)*t^2;
    roots = [r[0] for r in F2.roots()]
    maps = []; abscissas = [];
    for j2prime in roots:
        m2 = j2prime/j2 #1/ell * -G6tilde/G4tilde
        k2 = j2prime/(1728-j2) #G4tilde^2/G6
        Atilde = ell^4*m2*k2/48; Btilde = ell^6*m2^2*k2/864
        F3 = j1prime^3*Phix3(j1,j2)+3*ell*j1prime^2*j2prime*Phix2y(j1,j2)+3*ell^2*j1prime*j2prime^2*Phixy2(j1,j2)+ell^3*j2prime^3*Phiy3(j1,j2)
        logderivdiff = F3/(2*ell^2*j2prime^2*Phiy2(j1,j2)+2*ell*j1prime*j2prime*Phixy(j1,j2))
        p1 = ell*(logderivdiff/2 + (k1-ell*k2)/4 + (ell*m2-m1)/3)
        Etilde = EllipticCurve(K, [Atilde, Btilde]); phi = EllipticCurveIsogeny(E1, None, Etilde, ell);
        maps.append(phi); abscissas.append(p1)
    return maps, abscissas


#Test Functions
def singleabstest():
    p = 137; ell = 5; K.<w> = GF(p^2); j1 = K(22); j2 = 43*w+130; phi, p1 = abssingle(j1,j2, ell);
    print(phi.kernel_polynomial())
    print(-p1/2)

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
