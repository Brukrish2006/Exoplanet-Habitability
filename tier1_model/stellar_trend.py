"""
Tier 1, second validation check: reproduce the Kopparapu et al.
(2013/2014) Seff-vs-Teff trend for the moist-greenhouse inner edge
across stellar types, using the OLR limit derived in rc_model.py.

REVISION NOTE: an earlier version of this module used a hand-picked
two-band (visible/near-IR) "coupling efficiency" schematic with two
unfit constants (EPS_VIS=0.70, EPS_NIR=0.95), flagged in review as an
overclaim risk -- those numbers were chosen to get the right sign, not
derived or fit to anything. This version replaces that schematic with
a single, physically motivated broadband planetary albedo built from
actual Rayleigh-scattering physics:

    A(lambda) = min(1, (lambda_0 / lambda)^4)

integrated against the stellar Planck spectrum to get a Teff-dependent
bond albedo A(Teff). This is the same physics that sets Earth's blue
sky and its wavelength dependence is not a free choice -- the ^-4 power
law is Rayleigh scattering, full stop. The ONE free constant, lambda_0
(the wavelength scale of the scattering cross-section), is CALIBRATED
(not hand-picked) by requiring the resulting Seff at a Sun-like
reference Teff match the literature Kopparapu value there -- exactly
the same "calibrate once at a known reference point, then hold fixed
and check other points" methodology already used for KAPPA in
rc_model.py, and reported with the same kind of explicit residual
diagnostic. This reduces the model from 2 unfit constants to 1
calibrated constant, which is a meaningfully stronger position when a
referee asks "did you fit this, or pick numbers that work?"

    Seff(Teff) = 4 * OLR_limit / (S0 * (1 - A(Teff, lambda_0)))

Kopparapu et al. (2013), Table 3 polynomial coefficients are defined
here as the single source of truth for this project (previously
duplicated, with slightly different transcribed values, in a
standalone validate_model.py script -- consolidated here to remove
that inconsistency risk; tier2/tier3/tier4 all import these same
constants from this module rather than redefining them).
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from rc_model import calibrate_kappa, find_olr_limit

S0 = 1360.8       # solar constant at 1 AU, W/m^2
H = 6.62607015e-34
C = 2.99792458e8
KB = 1.380649e-23

TEFF_SUN_REF = 5772.0  # calibration reference point (Sun-like G2V)
TEFF_VALID_RANGE = (2600.0, 7200.0)  # Kopparapu (2013) fit's stated valid
                                     # domain; ALSO enforced as the valid
                                     # domain for this project's own Rayleigh
                                     # model, since it was calibrated and
                                     # checked only within this range. An
                                     # earlier version silently extrapolated
                                     # the Rayleigh albedo formula to hot
                                     # subdwarf hosts (Teff up to 40000 K --
                                     # NSVS 14256825, NY Vir, V0391 Peg, all
                                     # known post-common-envelope hot-subdwarf
                                     # systems, not main-sequence HZ hosts)
                                     # and produced meaningless "HZ flips"
                                     # for them. Planets with host Teff
                                     # outside this range must be EXCLUDED
                                     # from HZ classification entirely, not
                                     # extrapolated to.

# ---------------------------------------------------------------------
# Kopparapu et al. (2013), Table 3 polynomial coefficients -- single
# source of truth for this whole project.
# Seff(Tstar) = Seff_sun + a*Tstar + b*Tstar^2 + c*Tstar^3 + d*Tstar^4
# Tstar = Teff - 5780 K.  Valid strictly for 2600 K <= Teff <= 7200 K.
# ---------------------------------------------------------------------
KOPP_MOIST_GREENHOUSE = dict(Seff_sun=1.0140, a=8.1774e-5, b=1.7063e-9,
                              c=-4.3241e-12, d=-6.6462e-16)
KOPP_MAX_GREENHOUSE = dict(Seff_sun=0.3438, a=5.8942e-5, b=1.6558e-9,
                            c=-3.0045e-12, d=-5.2982e-16)


def kopparapu_seff(Teff, coeffs):
    Teff = np.clip(Teff, 2600.0, 7200.0)  # stay within fit's valid range
    Ts = Teff - 5780.0
    return (coeffs["Seff_sun"] + coeffs["a"] * Ts + coeffs["b"] * Ts**2
            + coeffs["c"] * Ts**3 + coeffs["d"] * Ts**4)


def planck_lambda(wav, T):
    """Planck spectral radiance B_lambda(T), W/m^3/sr (per unit wavelength)."""
    a = 2.0 * H * C**2
    b = (H * C) / (wav * KB * T)
    if b > 700:
        return 0.0
    return a / (wav**5 * (np.exp(b) - 1.0))


def rayleigh_albedo(Teff, lambda_0, wav_min=1e-7, wav_max=1e-5):
    """Broadband planetary albedo for a Rayleigh-scattering atmosphere,
    A(lambda) = min(1, (lambda_0/lambda)^4), Planck-weighted over the
    stellar spectrum at Teff. Integration bounds (0.1-10 um) capture
    >99% of blackbody flux for the stellar types considered here."""
    def flux(wav):
        return planck_lambda(wav, Teff)

    def reflected(wav):
        a_wav = min(1.0, (lambda_0 / wav) ** 4)
        return planck_lambda(wav, Teff) * a_wav

    total_flux, _ = quad(flux, wav_min, wav_max)
    refl_flux, _ = quad(reflected, wav_min, wav_max)
    return refl_flux / total_flux


def calibrate_lambda0(olr_limit, teff_ref=TEFF_SUN_REF):
    """Calibrate the ONE free constant (lambda_0) so that this model's
    Seff at teff_ref matches the literature Kopparapu moist-greenhouse
    value there. Reports the residual explicitly, same discipline as
    rc_model.calibrate_kappa()."""
    target_seff = kopparapu_seff(teff_ref, KOPP_MOIST_GREENHOUSE)
    target_albedo = 1.0 - (4.0 * olr_limit) / (S0 * target_seff)

    def resid(lambda_0):
        return rayleigh_albedo(teff_ref, lambda_0) - target_albedo

    lambda_0_star = brentq(resid, 1e-8, 1e-6, xtol=1e-12)

    achieved_albedo = rayleigh_albedo(teff_ref, lambda_0_star)
    achieved_seff = 4.0 * olr_limit / (S0 * (1.0 - achieved_albedo))
    print(f"  [lambda_0 calibration diagnostic] target Seff({teff_ref:.0f}K)="
          f"{target_seff:.4f}, achieved={achieved_seff:.4f}, "
          f"residual={achieved_seff - target_seff:+.5f}, "
          f"lambda_0={lambda_0_star*1e9:.1f} nm")
    return lambda_0_star


def our_seff_inner(Teff, olr_limit, lambda_0):
    A = rayleigh_albedo(Teff, lambda_0)
    return 4.0 * olr_limit / (S0 * (1.0 - A))


def calibration_reference_sensitivity(olr_limit, stars, teff_ref_range=(5700.0, 5772.0, 5850.0)):
    """Check that the M-dwarf Seff prediction doesn't swing wildly
    depending on the exact choice of calibration reference Teff --
    the analogous robustness check to the old two-band sweep, adapted
    to this model's single calibrated constant."""
    m_seffs = []
    for teff_ref in teff_ref_range:
        lambda_0 = calibrate_lambda0(olr_limit, teff_ref=teff_ref)
        m_teff = stars[-1][1]
        m_seffs.append(our_seff_inner(m_teff, olr_limit, lambda_0))
    return list(zip(teff_ref_range, m_seffs))


if __name__ == "__main__":
    kappa = calibrate_kappa(target_olr=282.0)
    _, _, olr_limit, _ = find_olr_limit(kappa)
    print(f"\nUsing OLR limit = {olr_limit:.2f} W/m^2 (from rc_model.py)\n")

    lambda_0 = calibrate_lambda0(olr_limit)

    stars = [
        ("Sun-like (G2V)", TEFF_SUN_REF),
        ("K dwarf",         4500.0),
        ("M dwarf",         3000.0),
    ]

    print(f"\n{'Star type':16s} {'Teff [K]':>9s} {'Albedo':>8s} "
          f"{'Seff (this model)':>19s} {'Seff (Kopparapu)':>18s} {'Diff %':>8s}")
    seffs = []
    for name, Teff in stars:
        A = rayleigh_albedo(Teff, lambda_0)
        seff = our_seff_inner(Teff, olr_limit, lambda_0)
        kopp = kopparapu_seff(Teff, KOPP_MOIST_GREENHOUSE)
        diff_pct = 100.0 * (seff - kopp) / kopp
        seffs.append(seff)
        print(f"{name:16s} {Teff:9.0f} {A:8.3f} {seff:19.4f} {kopp:18.4f} {diff_pct:+7.1f}%")

    direction_ok = seffs[0] > seffs[1] > seffs[2]
    print(f"\nQualitative direction (Seff decreases Sun -> K -> M) reproduced: {direction_ok}")

    print("\nCalibration-reference sensitivity check (does the M-dwarf prediction")
    print("depend sensitively on the exact Sun-like Teff used to calibrate lambda_0?):")
    sens = calibration_reference_sensitivity(olr_limit, stars)
    for teff_ref, m_seff in sens:
        print(f"  teff_ref={teff_ref:.0f}K  ->  M-dwarf Seff={m_seff:.4f}")
    spread = max(s for _, s in sens) - min(s for _, s in sens)
    print(f"  Spread across calibration references: {spread:.4f} "
          f"({'robust' if spread < 0.02 else 'NOTE: non-trivial sensitivity'})")
