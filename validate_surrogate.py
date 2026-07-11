import numpy as np
import sys
import os
from scipy.interpolate import RegularGridInterpolator

sys.path.append(os.path.abspath("tier1_model"))
import rc_model

def validate_surrogate():
    print("Loading surrogate grid...")
    npz = np.load('rc_model_grid.npz')
    olr_interp = RegularGridInterpolator((npz['kappas'], npz['t_strats']), npz['olr_grid'])
    
    np.random.seed(101)
    N_test = 50
    kappas_test = np.random.uniform(npz['kappas'].min(), npz['kappas'].max(), N_test)
    t_strats_test = np.random.uniform(180, 220, N_test)
    
    errors = []
    
    for i, (k, t_s) in enumerate(zip(kappas_test, t_strats_test)):
        surrogate_val = olr_interp([k, t_s])[0]
        _, _, true_val, _ = rc_model.find_olr_limit(k, T_strat=t_s)
        
        err = surrogate_val - true_val
        errors.append(err)
        print(f"Test {i+1}/50 | kappa={k:.4e}, T_strat={t_s:.1f} | Surrogate={surrogate_val:.3f}, True={true_val:.3f}, Err={err:+.3f}")
        
    errors = np.array(errors)
    print("\n--- Validation Summary ---")
    print(f"Max absolute error: {np.max(np.abs(errors)):.4f} W/m^2")
    print(f"Mean absolute error: {np.mean(np.abs(errors)):.4f} W/m^2")
    print(f"RMSE: {np.sqrt(np.mean(errors**2)):.4f} W/m^2")

if __name__ == "__main__":
    validate_surrogate()
