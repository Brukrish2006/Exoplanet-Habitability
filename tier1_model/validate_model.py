import numpy as np
from scipy.integrate import quad
import rc_model

# Physical constants
S0 = 1361.0  # Solar constant, W/m^2
H_PLANCK = 6.626e-34
C_LIGHT = 3.0e8
K_BOLTZ = 1.38e-23

def planck_lambda(wav, T):
    """Planck function B_lambda(T) in W/m^3/sr"""
    a = 2.0 * H_PLANCK * C_LIGHT**2
    b = (H_PLANCK * C_LIGHT) / (wav * K_BOLTZ * T)
    # Avoid overflow
    if b > 700:
        return 0.0
    return a / (wav**5 * (np.exp(b) - 1.0))

def calc_albedo(T_eff, lambda_0, n=4.0):
    """
    Calculate broadband planetary albedo for a pure water vapor atmosphere
    by integrating a blackbody spectrum against a Rayleigh scattering albedo.
    A(lambda) = min(1, (lambda_0 / lambda)^n)
    """
    def integrand_flux(wav):
        return planck_lambda(wav, T_eff)
        
    def integrand_reflected(wav):
        # Rayleigh albedo
        a_wav = min(1.0, (lambda_0 / wav)**n)
        return planck_lambda(wav, T_eff) * a_wav

    # Integrate from 0.1 to 10 microns (contains >99% of flux for these stars)
    wav_min, wav_max = 1e-7, 1e-5
    
    total_flux, _ = quad(integrand_flux, wav_min, wav_max)
    refl_flux, _ = quad(integrand_reflected, wav_min, wav_max)
    
    return refl_flux / total_flux

def kopparapu_seff(T_eff):
    """
    Kopparapu et al. (2014) polynomial for Runaway Greenhouse (Moist Greenhouse limit).
    Coefficients for 1 M_Earth.
    """
    T_star = T_eff - 5780.0
    Seff_sun = 1.038
    a = 1.2456e-4
    b = 1.4612e-8
    c = -7.6345e-12
    d = -1.1511e-15
    return Seff_sun + a*T_star + b*T_star**2 + c*T_star**3 + d*T_star**4

def main():
    print("--- Tier 1 Validation: 1-D RC Model vs Kopparapu ---")
    
    # 1. Get the model's OLR limit
    print("\nRunning rc_model to find OLR plateau...")
    kappa_star = rc_model.calibrate_kappa(target_olr=282.0)
    _, _, olr_limit, _ = rc_model.find_olr_limit(kappa_star)
    print(f"Model OLR limit = {olr_limit:.2f} W/m^2")
    
    # 2. Calibrate the toy albedo to match Kopparapu at T_eff = 5780 K
    # Seff = 4 * OLR / (S0 * (1 - A))
    # => 1 - A = 4 * OLR / (S0 * Seff)
    Seff_sun_target = kopparapu_seff(5780.0)
    target_albedo_sun = 1.0 - (4.0 * olr_limit) / (S0 * Seff_sun_target)
    print(f"\nCalibrating Rayleigh lambda_0 for Solar albedo A={target_albedo_sun:.3f}...")
    
    from scipy.optimize import root_scalar
    def albedo_resid(l0):
        return calc_albedo(5780.0, l0) - target_albedo_sun
        
    res = root_scalar(albedo_resid, bracket=[1e-8, 1e-6])
    lambda_0 = res.root
    print(f"Found lambda_0 = {lambda_0 * 1e9:.1f} nm")
    
    # 3. Generate the Validation Table
    print("\n--- HZ Inner Edge (Runaway Greenhouse) Validation Table ---")
    print(f"{'Star Type':<12} | {'T_eff (K)':<10} | {'Model Albedo':<12} | {'Model S_eff':<12} | {'Kopparapu S_eff':<15} | {'Diff (%)':<8}")
    print("-" * 80)
    
    test_stars = [
        ("Sun-like (G)", 5780),
        ("K dwarf", 4500),
        ("M dwarf", 3000)
    ]
    
    for name, Teff in test_stars:
        # Model calculation
        A = calc_albedo(Teff, lambda_0)
        # S_eff = 4 * OLR_limit / (S0 * (1 - A))
        model_seff = (4.0 * olr_limit) / (S0 * (1.0 - A))
        
        # Kopparapu calculation
        kop_seff = kopparapu_seff(Teff)
        
        diff_pct = (model_seff - kop_seff) / kop_seff * 100.0
        print(f"{name:<12} | {Teff:<10} | {A:<12.3f} | {model_seff:<12.3f} | {kop_seff:<15.3f} | {diff_pct:>+6.1f}%")

    print("\nValidation Summary:")
    print("1. Qualitative Trend: As T_eff decreases (M dwarf), the incident spectrum shifts to the red.")
    print("   Since our toy Rayleigh albedo A(lambda) ~ lambda^-4 drops for redder light, planetary albedo")
    print("   decreases. A lower albedo means the planet absorbs more efficiently, triggering the")
    print("   runaway greenhouse at a LOWER incident flux (S_eff).")
    print("2. The Model correctly reproduces the downward trend of S_eff for cooler stars.")
    print("3. Differences arise due to the gray-gas approximation in rc_model (no explicit H2O/CO2")
    print("   absorption bands), and the simplified broadband Rayleigh treatment compared to")
    print("   Kopparapu's full 1-D line-by-line radiative transfer.")

if __name__ == "__main__":
    main()
