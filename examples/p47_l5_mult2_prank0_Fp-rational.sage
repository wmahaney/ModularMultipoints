"""
A pair of supersingular elliptic curves connected by 2 distinct 5-isogenies.
"""

load(os.path.join(PROJECT_ROOT, "main.sage"))

p=47
l=5
k.<w>=GF(p^2)
j1=k(9)
j2=k(10)
mult=2


E = EllipticCurve_from_j(j1) 
E=E.short_weierstrass_model()

p_rank = int(E.is_ordinary())

metadata = {
    "k": str(k),    
    "k_minimal_polynomial": str(w.minimal_polynomial()),
    "p": k.characteristic(),
    "l": l,
    "j1": str(j1),
    "j2": str(j2),
    "multiplicity": str(mult),
    "domain_curve": str(E),
     "p_rank": p_rank 
}

A=E.a4()
B=E.a6() 
mu1 = 18*B/A
nu1 = mu1*j1/(1728-j1) 
j1prime = mu1*j1 



Phi = classical_modular_polynomial(l)
RXY.<X,Y>=k[]
Rt.<t>=k[]
Phi=RXY(Phi)

#assert the singularity has multiplicity 2

PhiX = derivative(Phi, X, 1, Y, 0)
PhiX=PhiX(j1,j2)
PhiY = derivative(Phi, X, 0, Y, 1)
PhiY=PhiY(j1,j2)

assert Phi(j1 ,j2) == 0
assert PhiX == 0
assert PhiY == 0



PhiXY = derivative(Phi, X, 1, Y, 1)
PhiXY=PhiXY(j1,j2)
PhiXX = derivative(Phi, X, 2, Y, 0)
PhiXX = PhiXX(j1,j2)
PhiYY = derivative(Phi, X, 0, Y, 2)
PhiYY = PhiYY(j1,j2)


PhiXXX = derivative(Phi, X, 3, Y, 0)
PhiXXX = PhiXXX(j1, j2)
PhiXXY = derivative(Phi, X, 2, Y, 1)
PhiXXY = PhiXXY(j1, j2)
PhiXYY = derivative(Phi, X, 1, Y, 2)
PhiXYY = PhiXYY(j1, j2)
PhiYYY = derivative(Phi, X, 0, Y, 3)
PhiYYY = PhiYYY(j1, j2)


#compute derivatives evaluated at (j1, j2)
derivative_evals = {}
#This can be sped up instantiated modular polynomials.
for u in range(0, mult+2):
    for v in range(0, mult+2):
        Phi_deriv = derivative(Phi, X, u, Y, v)
        Phi_deriv  = Phi_deriv(j1, j2)
        derivative_evals[(u,v)] = Phi_deriv

#Compute fiber polynomial
fiber_poly =0
for u in range(0, mult+1):
   fiber_poly += binomial(mult,u) * l^(mult-u) * j1prime^u * derivative_evals[(u, mult-u)] * t^(mult-u)

#get the roots of the fiber polynomial and then return the associated models and abscissas
metadata['fiber_polynomial'] = str(fiber_poly)

L.<z>=fiber_poly.splitting_field()

#record splitting field data for the fiber polynomial
metadata['splitting_field'] = str(L)
metadata['splitting_field_minimal_polynomial'] = str(L.modulus())



fiber_poly_ddt =fiber_poly.derivative(t)

roots_data = [r[0] for r in fiber_poly.roots()]
isogeny_data = {}
for r in roots_data:
    datum={}
    Atilde = ( -l^4 * r^2)/(48*j2*(j2-1728))
    Btilde = (-l^6 * r^3)/(864*j2^2*(j2-1728))
    model= [Atilde, Btilde]
    datum['model']=model
    Etilde_r = EllipticCurve(k, model)
    datum['codomain_curve'] = str(Etilde_r)

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
    m1 = 18*B/A; n1 = j1prime/(1728-j1)
    m2 = r/j2; n2 = r/(1728-j2) 
    abscissa = l * (log_deriv_diff/2 + (n1 - l*n2)/4 + (l*m2 - m1)/3 )

    #now use Fast Elkies to compute the kernel polynomial
    kernel_poly = fast_elkies(E, Etilde_r, l, abscissa)
    datum['abscissa'] = abscissa
    datum['kernel_polynomial'] = str(kernel_poly)

    #assert that our algorithm worked
    phi = EllipticCurveIsogeny(E, kernel_poly)
    assert phi.codomain() == Etilde_r 

    isogeny_data[r]=datum 

metadata['isogeny_data'] = isogeny_data

output_dir = os.path.join(str(PROJECT_ROOT), "examples", "metadata")
os.makedirs(output_dir, exist_ok=True) 

# Construct filename from metadata
p_val = metadata['p']
l_val = metadata['l']
mult_val = metadata['multiplicity']
p_rank_val = metadata['p_rank']
base_filename = f"metadata_p{p_val}_l{l_val}_mult{mult_val}_prank{p_rank_val}"

# Record everything in a JSON file (sanitize all values as strings)
def stringify(obj):
    if isinstance(obj, dict):
        return {str(k): stringify(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [stringify(i) for i in obj]
    else:
        try:
            return str(obj)
        except Exception:
            return repr(obj)

json_filename = os.path.join(output_dir, f"{base_filename}.json")
with open(json_filename, "w") as jsonfile:
    json.dump(stringify(metadata), jsonfile, indent=4)




