"""
A pair of ordinary elliptic curves connected by 2 distinct 3-isogenies. In this example the fiber polynomial does not split over the base field of the elliptic curves.
"""

import os
import sys
import csv
import json

load(os.path.join(PROJECT_ROOT, "main.sage"))

k0.<z>=QQ.extension(x^2+5)
#we extend to the Hilbert class field
k.<w>=QQ.extension(x^2-5)
l=3
mult=2

kt.<t>=k[] 
f=kt(k0.hilbert_class_polynomial())
root_data = f.roots()
j1 = root_data[0][0]
j2 = root_data[1][0]

E = EllipticCurve_from_j(j1) 
E=E.short_weierstrass_model()
A=E.a4()
B=E.a6()
j1prime = j1*18*B/A

metadata = {
    "k": str(k),
    "k_base_field": str(k.base_field()),
    "l": l,
    "j1": str(j1),
    "j2": str(j2),
    "multiplicity": str(mult),
    "domain_curve": str(E)
}

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

#get the roots of the fiber polynomial and then return the associated models and abscissas
metadata['fiber_polynomial'] = str(fiber_poly)

L.<z>=fiber_poly.splitting_field()

#record splitting field data for the fiber polynomial
metadata['splitting_field'] = str(L)
metadata['splitting_field_minimal_polynomial'] = str(L.defining_polynomial())

isogeny_data, _ = multipoint_isogeny(j1,j2,l, mult=2)
metadata['isogeny_data']=isogeny_data 


output_dir = os.path.join(str(PROJECT_ROOT), "examples", "metadata")
os.makedirs(output_dir, exist_ok=True) 


# Construct filename from metadata
field_val = "Q-sqrtpm5"
l_val = metadata['l']
mult_val = metadata['multiplicity']
base_filename = f"metadata_{field_val}_l{l_val}_mult{mult_val}"

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