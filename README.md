This repo contains code for computing $\ell$-isogenies between elliptic curves from knowledge of their $j$-invariants, even when those j-invariants are connected by multiple distinct $\ell$-isogenies.

Based on "Computing Isogenies at Singular Points of the Modular Polynomial" by William E. Mahaney and Travis Morrison available at
<https://arxiv.org/abs/2402.02038>.

# Build Instructions 
0) Make sure you have [Sagemath](https://www.sagemath.org/) 10.4 or later installed. This code is intended for Linux users but should work on Mac, Windows users need [WSL](https://learn.microsoft.com/en-us/windows/wsl/install).
    If your Sage throws errors try updating to 10.7.  
1) Clone the git repo to a location of your choice. 
2) Type ```make``` to create the PROJECT_ROOT files needed for the examples and main functions work. 
3) To run all the examples change to the ```examples``` directory. Do ```sage run_all.sage``` or open Sage and run ```load('run_all.sage')```.
4) If you move the repo after installation, you must run ```make``` again to reset PROJECT_ROOT.
5) To remove the PROJECT_ROOT files and /examples/metadata files run ```make uninstall``` from the repo home directory. 

# Primary Functions

## main.sage
- ```isogeny_data, isogenies = multipoint_isogeny(j1, j2, l, mult = None, model = None)```
  Intakes a pair of $j$-invariants $j_1, j_2 \in k$ for elliptic curves over a field $k$ which is either a finite field or a subfield of $\mathbb{C}$. Optionally intakes a positive integer ```mult``` and a short Weierstrass model $(A, B)$ for an elliptic curve $E$ with $j$-invariant $j_1$; if not given ```mult``` is computed as part of the execution and ```model``` is computed by using ```E = EllipticCurve_from_j(j1)``` and coercing to a short Weierstrass model.
  
  Outputs a dictionary ```isogeny_data``` containing an entry for each $\ell$-isogeny out of $E$ to an elliptic curve $j_1$ whose values are the codomain of the isogeny, the abscissa of the isogeny, and the kernel polynomial; ```isogenies``` is a list of the EllipticCurveIsogeny objects computed from the isogeny data. 

## elkies_isogeny.sage 
- ```f = fast_elkies(E1, E2, l, sigma)```

# Examples 
- ```run_all.sage```
  Runs all example files.
  
- ```p29_l7_mult4_prank0_spinal.sage```
  A pair of supersingular elliptic curves with $j$-invariants in $\mathbb{F}_{29}$ connected by 4 distinct 7-isogenies.
  
- ```p37_l5_mult2_prank1_spinal.sage```
  A pair of ordinary elliptic curves with $j$-invariants in $\mathbb{F}_{37}$ connected by 2 distinct 5-isogenies.
  
- ```p47_l5_mult2_prank0_spinal.sage```
  A pair of supersingular elliptic curves with $j$ invariants in $\mathbb{F}_{47}$ connected by 2 distinct 5-isogenies.
