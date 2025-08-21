import os 
import sys

ex_dir = os.path.dirname(os.path.abspath("run_all.sage"))

for f in os.listdir(ex_dir):
    if f != "run_all.sage" and f.endswith(".sage"):
        print(f"Running {f}")
        load(os.path.join(ex_dir, f)) 


