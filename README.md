This repo contains code for computing $\ell$-isogenies between elliptic curves from knowledge of their $j$-invariants, even when those j-invariants are connected by multiple distinct $\ell$-isogenies.

Based on *Computing Isogenies at Singular Points of the Modular Polynomial* by William E. Mahaney and Travis Morrison available at  
<https://arxiv.org/abs/2402.02038>.

# Build Instructions 

1. Make sure you have [Sagemath](https://www.sagemath.org/) 10.4 or later installed.  
   - This code is intended for Linux users but should work on Mac.  
   - Windows users need [WSL](https://learn.microsoft.com/en-us/windows/wsl/install).  
   - Alternatively, you can use [Cocalc](https://cocalc.com/).  
   - If your Sage throws errors, try updating to 10.7.  
   - If you have the Kohel modular polynomial database installed or your own source of modular polynomials, you can edit `main.sage` to change how modular polynomials are obtained.  

2. Start by cloning the repo and going inside `clone https://github.com/wmahaney/ModularMultipoints.git && cd ModularMultipoints`.

3. Do `make` to assign the `PROJECT_ROOT` variable in the `.sage` files.
   - This creates `ModularMultipoints/backups/` containing copies of the `.sage` files at the time of calling make. If you edit the `.sage` files after calling `make` then calling `make restore` or `make uninstall` will **delete all your changes**, go through the sage files and delete the created code injections and backup files before pushing or save your changes in separate text files before calling `make uninstall` then make the changes and push. 

5. To run all the examples:
   - Change to the `examples` directory
   - Run `sage run_all.sage`  
   - Or open Sage and run `load('run_all.sage')` 

6. If you move the repo after installation, you must run `make uninstall` and then `make` to reset `PROJECT_ROOT`. 

7. To remove the `PROJECT_ROOT` files and `/examples/metadata` files, run `make uninstall` from the repo home directory.  

# Primary Functions

## `main.sage`

- `isogeny_data, isogenies = multipoint_isogeny(j1, j2, l, mult = None, model = None)`  
  - Inputs:  
    - A pair of $j$-invariants $j_1, j_2 \in k$ for elliptic curves over a field $k$ (either a finite field or a subfield of $\mathbb{C}$).  
    - Optional:  
      - `mult`: A positive integer. 
      - `model`: A short Weierstrass model $(A, B)$ for an elliptic curve $E$ with $j$-invariant $j_1$.  
    - If not given, `mult` is computed during execution and `model` is computed using `E = EllipticCurve_from_j(j1)` and coercing to a short Weierstrass model.  

  - Outputs:  
    - `isogeny_data`: A dictionary of dictionaries containing an entry for each $\ell$-isogeny out of $E$ to an elliptic curve with $j$-invariant $j_2$.
        - The keys of `isogeny_data` are the roots of the associated fiber polynomial. The keys of the values of `isogeny_data` are: `codomain_curve`, `isogeny_abscissa`, and `isogeny_kernel_polynomial`.   
    - `isogenies`: A list of the `EllipticCurveIsogeny` objects computed from the isogeny data.  

## `elkies_isogeny.sage`
- `f = fast_elkies(E1, E2, l, sigma)`
    - Inputs:
          - A pair of $\ell$-isogenous elliptic curves $E_1, E_2$ for $\ell$ an *odd* prime.
          - The abscissa $\sigma$ of a normalized $\ell$-isogeny from $E_1$ to $E_2$. 
    - Outputs:
          - `f`: The kernel polynomial of the normalized $\ell$-isogeny from $E_1$ to $E_2$ with abscissa $\sigma$. 

# Examples 

- `run_all.sage`  
  Runs all example files and stores the metadata to `ModularMultipoints/examples/metadata/`. 

- `p29_l7_mult4_prank0_spinal.sage`  
  A pair of supersingular elliptic curves defined over $\mathbb{F}_{29}$ connected by 4 distinct 7-isogenies. 

- `p37_l5_mult2_prank1_spinal.sage`  
  A pair of ordinary elliptic curves defined over $\mathbb{F}_{37}$ connected by 2 distinct 5-isogenies.  

- `p47_l5_mult2_prank0_spinal.sage`  
  A pair of supersingular elliptic curves defined over $\mathbb{F}_{47}$ connected by 2 distinct 5-isogenies.

- `number-field_l3_mult2.sage`
  A pair of ordinary elliptic curves defined over $\mathbb{Q}(\sqrt{5}, \sqrt{-5})$ connected by 2 distinct 3-isogenies. 
  - In this example the fiber polynomial connecting $E_1$ and $E_2$ does not split over their field of definition. 
