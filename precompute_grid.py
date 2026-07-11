import numpy as np
import os
import sys

# Ensure tier1_model is in path
sys.path.append(os.path.abspath("tier1_model"))
import rc_model

def precompute():
    print("Precomputing OLR limit grid for surrogate model...")
    kappas = np.linspace(0.02, 0.045, 6)
    t_strats = np.linspace(180, 220, 5)
    
    olr_grid = np.zeros((len(kappas), len(t_strats)))
    
    total = len(kappas) * len(t_strats)
    count = 0
    
    for i, kappa in enumerate(kappas):
        for j, t_strat in enumerate(t_strats):
            count += 1
            print(f"Computing {count}/{total}: kappa={kappa:.4f}, T_strat={t_strat:.1f}")
            _, _, olr_limit, _ = rc_model.find_olr_limit(kappa, T_strat=t_strat)
            olr_grid[i, j] = olr_limit
            
    np.savez('rc_model_grid.npz', kappas=kappas, t_strats=t_strats, olr_grid=olr_grid)
    print("Grid saved to rc_model_grid.npz")

if __name__ == "__main__":
    precompute()
