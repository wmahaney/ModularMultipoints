This repo contains code for computing $\ell$-isogenies between elliptic curves from knowledge of their $j$-invariants, even when those j-invariants are connected by multiple distinct $\ell$-isogenies.

Based on "Computing Isogenies at Singular Points of the Modular Polynomial" by William E. Mahaney and Travis Morrison available at
<https://arxiv.org/abs/2402.02038>.

# Build Instructions 
0) Make sure you have Sagemath 10.4 or later installed. 
    If this does not work update to 10.7 or later. 
1) Clone the git repo to a location of your choice.
2) Type ```make``` to place the PROJECT_ROOT files needed for functionality of the examples.
3) To run all the examples change to the ```examples``` directory. Do ```sage run_all.sage``` or open sage and run ```load('run_all.sage')```.
4) To reset the PROJECT_ROOT if you move the git repo in your computer run ```make``` from the repo home directory. 
5) To remove the PROJECT_ROOT files and /examples/metadata files run ```make uninstall``` from the repo home directory. 

# Functions

# main.sage

## elkies_isogeny.sage 

# Examples 