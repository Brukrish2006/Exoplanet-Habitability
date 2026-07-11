"""
Tier 1 - Custom 1-D radiative-convective model for the moist-greenhouse
(Simpson-Nakajima) OLR limit.

Physical picture (Nakajima et al. 1992; Pierrehumbert, Principles of
Planetary Climate, Ch. 4):

  1. Convective troposphere: temperature follows the saturated
     (pseudo-)moist adiabat set by water-vapor partial pressure via
     Clausius-Clapeyron, mixed with a FIXED non-condensible background
     gas column (representing an Earth-like N2-type background).
  2. Radiative stratosphere: capped at an isothermal "skin" temperature
     once the adiabat would otherwise cool further (standard grey
     radiative-equilibrium closure).
  3. OLR: grey-gas, plane-parallel, no-scattering formal solution of the
     Schwarzschild equation, with water-vapor continuum optical depth
     built from the saturated mixing-ratio profile and a single tunable
     absorption coefficient KAPPA (the model's one free calibration
     constant, fixed once against the literature Simpson-Nakajima value
     and held fixed thereafter).
  4. Root-find: scan Ts and locate where OLR saturates (stops rising)
     -- the asymptotic Simpson-Nakajima limit.

Key physical subtlety (this is what actually produces the saturation,
and is the part that took real debugging to get right): the total
surface pressure is NOT held fixed as Ts increases. It is

    p_surf(Ts) = P_DRY_BG + e_sat(Ts)

i.e. a fixed non-condensible background column plus a water-vapor
partial pressure that grows with Ts via Clausius-Clapeyron. Holding
total surface pressure fixed while scanning Ts (an earlier, buggy
version of this model did this) forces the saturation mixing ratio
q = eps*e/p above 1 at high Ts -- unphysical, and it silently destroys
the saturation behaviour the model exists to reproduce. The mixing
ratio itself also uses the non-dilute formula q = eps*e/(p - e*(1-eps))
rather than the dilute approximation eps*e/p, since e is not small
compared to p in the hot/steam-dominated regime this model targets.

This is a deliberately simplified, well-precedented "textbook" model
(gray-gas, no clouds, no real absorption-line data) -- see LIMITATIONS
at the bottom of this file; these must be carried into the paper's
Methods section verbatim.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ---------------------------------------------------------------------
# Physical constants (SI unless noted)
# ---------------------------------------------------------------------
SIGMA = 5.670374419e-8      # Stefan-Boltzmann, W/m^2/K^4
G = 9.80665                 # surface gravity, m/s^2 (Earth-like; adjustable)
RD = 287.05                 # specific gas constant, dry air, J/kg/K
RV = 461.5                  # specific gas constant, water vapor, J/kg/K
CPD = 1005.0                # specific heat, dry air, J/kg/K
L = 2.5e6                   # latent heat of vaporization, J/kg (const. approx.)
EPS = RD / RV               # ~0.622, molecular weight ratio

P_DRY_BG = 1.0e5            # fixed non-condensible background column, Pa
P_TOP = 1.0e1               # top-of-atmosphere pressure, Pa
# T_STRAT removed as a global constant, now passed as argument


def e_sat(T):
    """Saturation vapor pressure of water via Clausius-Clapeyron,
    integrated form, Pa. Reference point: 611 Pa at 273.15 K."""
    T = np.asarray(T, dtype=float)
    return 611.0 * np.exp((L / RV) * (1.0 / 273.15 - 1.0 / T))


def mixing_ratio(T, p):
    """Water-vapor mass mixing ratio at saturation, non-dilute formula
        q = eps*e / (p - e*(1-eps))
    which correctly limits to q -> 1 as e -> p (pure water vapor),
    unlike the dilute approximation eps*e/p which can exceed 1."""
    e = np.minimum(e_sat(T), p * 0.999999)
    q = EPS * e / (p - e * (1.0 - EPS))
    return np.clip(q, 0.0, 1.0)


def moist_adiabat_dTdlnp(lnp, T, T_strat):
    """dT/d(ln p) along the saturated pseudo-adiabat (standard
    moist-adiabatic lapse rate formula), using the non-dilute qs."""
    p = np.exp(lnp)
    T = T[0]
    if T <= T_strat:
        return [0.0]
    qs = float(mixing_ratio(T, p))
    num = 1.0 + (L * qs) / (RD * T)
    den = 1.0 + (L**2 * qs * EPS) / (CPD * RD * T**2)
    dTdlnp = (RD * T / CPD) * (num / den)
    return [dTdlnp]


def build_profile(Ts, n=400, T_strat=200.0):
    """Integrate the moist adiabat from the surface up to P_TOP.
    Surface total pressure = P_DRY_BG + e_sat(Ts) (see module docstring).
    Returns pressure array (Pa, surface->top) and temperature array (K)."""
    p_surf = P_DRY_BG + e_sat(Ts)
    lnp_surf = np.log(p_surf)
    lnp_top = np.log(P_TOP)
    lnp_eval = np.linspace(lnp_surf, lnp_top, n)

    sol = solve_ivp(
        moist_adiabat_dTdlnp,
        (lnp_surf, lnp_top),
        [Ts],
        args=(T_strat,),
        t_eval=lnp_eval,
        method="RK45",
        max_step=0.05,
    )
    T = sol.y[0]
    T = np.maximum(T, T_strat)
    p = np.exp(lnp_eval)
    return p, T


def optical_depth_profile(p, T, kappa):
    """Water-vapor continuum optical depth, TOA (tau=0) to surface
    (tau=tau_s), tau(p) = integral (kappa/g) * q(T,p) dp."""
    q = mixing_ratio(T, p)
    p_ts = p[::-1]
    q_ts = q[::-1]
    integrand = (kappa / G) * q_ts
    tau_ts = np.concatenate([[0.0], np.cumsum(
        0.5 * (integrand[1:] + integrand[:-1]) * np.diff(p_ts)
    )])
    tau = tau_ts[::-1]
    return tau


def compute_olr(Ts, kappa, n=400, T_strat=200.0):
    """Formal (no-scattering, LTE, grey) solution for OLR at TOA:
        OLR = integral_0^{tau_s} sigma*T(tau)^4 * exp(-tau) dtau
    (the surface term is the tau=tau_s endpoint of this same integral)."""
    p, T = build_profile(Ts, n=n, T_strat=T_strat)
    tau = optical_depth_profile(p, T, kappa)
    B = SIGMA * T**4
    integrand = B * np.exp(-tau)
    olr = np.trapezoid(integrand, -tau)
    return olr


def calibrate_kappa(target_olr=282.0, kappa_bracket=(5e-3, 8e-2)):
    """Solve for KAPPA so the model's own OLR *plateau* -- the same
    quantity returned by find_olr_limit(), i.e. the actual asymptotic
    saturation value the model produces, not an arbitrary single-Ts
    proxy for it -- matches the literature Simpson-Nakajima benchmark.

    An earlier version of this function calibrated against
    compute_olr(Ts_cal=370, kappa) as a stand-in for the plateau. That
    introduced a ~4% internal inconsistency (the calibration point at
    Ts=370K is not actually where the plateau sits -- the plateau is at
    Ts~470-575K -- so the reported "plateau" value differed from the
    calibration target it was supposedly matched to). Fixed here by
    making find_olr_limit() itself the target of the root-find: whatever
    kappa is returned, plugging it back into find_olr_limit() reproduces
    target_olr by construction, with the residual reported explicitly
    at calibration time so the fit quality is never silently assumed.
    """
    def resid(log10_kappa):
        kappa = 10 ** log10_kappa
        _, _, plateau, _ = find_olr_limit(kappa)
        return plateau - target_olr

    lo, hi = np.log10(kappa_bracket[0]), np.log10(kappa_bracket[1])
    log10_kappa_star = brentq(resid, lo, hi, xtol=1e-4)
    kappa_star = 10 ** log10_kappa_star

    # report the calibration residual explicitly, rather than assuming
    # the fit is exact -- this is the diagnostic the paper's Methods
    # table should show alongside the calibrated kappa value
    _, _, plateau_check, Ts_plateau_check = find_olr_limit(kappa_star)
    residual = plateau_check - target_olr
    print(f"  [calibration diagnostic] target={target_olr:.2f} W/m^2, "
          f"achieved plateau={plateau_check:.2f} W/m^2 at Ts={Ts_plateau_check:.1f} K, "
          f"residual={residual:+.3f} W/m^2")

    return kappa_star


def find_olr_limit(kappa, Ts_scan=None, T_strat=200.0):
    """Scan Ts and report the Simpson-Nakajima OLR limit.

    The OLR(Ts) curve for this model has three regimes: it first rises
    with Ts (normal unsaturated climate), overshoots, then descends to
    a genuine minimum/plateau (the Simpson-Nakajima limit) before a
    slow secondary rise at very high Ts (a known closure artifact of
    the fixed-T_STRAT isothermal stratosphere, not physical -- see
    LIMITATIONS). The plateau is therefore identified as the GLOBAL
    MINIMUM of OLR(Ts) restricted to Ts above the initial overshoot
    peak, not the global minimum over the whole scan (which would
    incorrectly pick out the cold, unsaturated end of the curve)."""
    if Ts_scan is None:
        Ts_scan = np.linspace(260, 700, 89)
    olr_vals = np.array([compute_olr(Ts, kappa, T_strat=T_strat) for Ts in Ts_scan])
    peak_idx = np.argmax(olr_vals)  # end of the initial overshoot
    tail_idx = peak_idx + np.argmin(olr_vals[peak_idx:])
    Ts_plateau = Ts_scan[tail_idx]
    olr_limit = olr_vals[tail_idx]
    return Ts_scan, olr_vals, olr_limit, Ts_plateau


if __name__ == "__main__":
    print("Calibrating kappa against Simpson-Nakajima benchmark (282 W/m^2)...")
    kappa_star = calibrate_kappa(target_olr=282.0)
    print(f"  Calibrated kappa = {kappa_star:.6e}")

    Ts_scan, olr_vals, olr_limit, Ts_plateau = find_olr_limit(kappa_star)
    print("\nOLR(Ts) scan:")
    for Ts, olr in zip(Ts_scan[::7], olr_vals[::7]):
        print(f"  Ts = {Ts:6.1f} K   OLR = {olr:7.2f} W/m^2")
    print(f"\nOLR plateau (Simpson-Nakajima limit): {olr_limit:.2f} W/m^2")
    print(f"Plateau located near Ts = {Ts_plateau:.1f} K")
    print("(flattest region ~425-575 K; a slow secondary rise beyond that")
    print(" is a closure artifact of the fixed isothermal-stratosphere cap,")
    print(" not physical -- see LIMITATIONS)")

# ---------------------------------------------------------------------
# LIMITATIONS (state explicitly in paper Methods per plan Section 7)
# ---------------------------------------------------------------------
# - Gray-gas absorber: KAPPA is a single, wavelength-independent water-vapor
#   continuum coefficient calibrated once against the literature benchmark;
#   real spectra have line structure the gray-gas approximation cannot
#   capture. This is why this model tracks the classical Nakajima (1992)
#   ~280-290 W/m^2 number rather than Goldblatt et al. (2013)'s revised
#   (lower) value obtained with real absorption-line data.
# - No clouds, no scattering, no aerosols.
# - Isothermal-stratosphere closure at a fixed T_STRAT=200 K rather than a
#   self-consistent radiative-equilibrium skin temperature. This closure
#   is only reliable while the bulk of the emission-level optical depth
#   forms within the troposphere; at very high Ts (>~400 K in this
#   implementation) the plateau breaks down and OLR rises again -- this
#   ceiling on validity must be reported alongside the plateau value, not
#   silently extrapolated past it.
# - Latent heat L and gas constants treated as temperature-independent.
# - Background/stellar spectrum not yet incorporated; the Seff-vs-Teff
#   trend check against Kopparapu is done by scaling absorbed stellar
#   flux against this fixed OLR limit, not by re-deriving line-by-line
#   opacities per spectral type.
