"""
A pair of supersingular elliptic curves connected by 4 distinct 7-isogenies.
"""

import os 
import sys
import json

load(os.path.join(PROJECT_ROOT, "main.sage"))

p=37
l=11
k.<w>=GF(p^2)
j1=10*w+20
j2=27*w+23 
mult=7 


E = EllipticCurve_from_j(j1)
E=E.short_weierstrass_model()


metadata = {
    "k": str(k),    
    "k_minimal_polynomial": str(w.minimal_polynomial()),
    "p": k.characteristic(),
    "l": l,
    "j1": str(j1),
    "j2": str(j2),
    "multiplicity": str(mult),
    "domain_curve": str(E),
    "is_ordinary": E.is_ordinary()
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

#get the roots of the fiber polynomial and then return the associated models and sigma valuess
metadata['fiber_polynomial'] = str(fiber_poly)

L.<z>=fiber_poly.splitting_field()

#record splitting field data for the fiber polynomial
metadata['splitting_field'] = str(L)
metadata['splitting_field_minimal_polynomial'] = str(L.modulus())


isogeny_data = multipoint_isogeny(j1,j2,l)

metadata['isogeny_data'] = isogeny_data

output_dir = os.path.join(str(PROJECT_ROOT), "examples", "metadata")
os.makedirs(output_dir, exist_ok=True) 


# Construct filename from metadata
p_val = metadata['p']
l_val = metadata['l']
mult_val = metadata['multiplicity']
base_filename = f"metadata_p{p_val}_l{l_val}_mult{mult_val}_supersingular"

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





