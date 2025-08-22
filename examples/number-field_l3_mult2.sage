"""
A pair of ordinary elliptic curves connected by 2 distinct 3-isogenies.
"""

import os 
import sys
import csv
import json

#main.sage 
ex_dir=os.path.dirname(os.path.abspath("number-field_l3_mult2.sage"))
root_file = os.path.join(ex_dir, "PROJECT_ROOT")

with open(root_file, 'r') as f:
    PROJECT_ROOT= f.read().strip()

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

metadata = {
    "k": str(k),
    "k_base_field": str(k.base_field()),
    "l": l,
    "j1": str(j1),
    "j2": str(j2),
    "multiplicity": str(mult),
    "domain_curve": str(E),
    "p_rank": 1
}


isogeny_data, _ = multipoint_isogeny(j1,j2,l, mult=2)
metadata['isogeny_data']=isogeny_data 


output_dir = os.path.join(str(PROJECT_ROOT), "examples", "metadata")
os.makedirs(output_dir, exist_ok=True) 


# Construct filename from metadata
field_val = "Q(sqrt(-5))(sqrt(5))"
l_val = metadata['l']
mult_val = metadata['multiplicity']
p_rank_val = metadata['p_rank']
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