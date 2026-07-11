"""
Tier 3 -- Two-layer uncertainty budget + simplified hierarchical
Bayesian occurrence-rate estimate.

Layer A (measurement uncertainty): Monte Carlo propagation of catalog
parameter errors (pl_insol, st_teff -- the two quantities the HZ
classification actually depends on here) into a per-planet posterior
probability of HZ membership, P(HZ), holding the boundary MODEL fixed.
Asymmetric NASA-Archive error bars (err1 >= 0 upper, err2 <= 0 lower)
are respected via a two-piece ("split") normal distribution rather than
collapsing them into one symmetric sigma.

Layer B (model uncertainty): the Tier 2 boundary-choice spread --
whether a planet's HZ status depends on which inner-edge model
(Kopparapu vs. this project's derived boundary) is used -- reported
as its own, separate quantity. Layers A and B are NOT combined into a
single number: the whole point of this exercise (per the research
plan) is that "measurement noise" and "which physical model you trust"
are different kinds of uncertainty and conflating them hides which one
actually dominates for a given planet.

Bayesian layer: a simplified hierarchical estimate of an eta-Earth-like
occurrence rate per spectral-type bin, with a Beta-distributed
detection-completeness nuisance parameter weakly anchored to Bryson
et al. (2021)'s eta_Earth = 0.37 (+0.48/-0.22, 95% CI) for Sun-like
stars, reused as a single generic prior across all bins -- this reuse
is an explicit, stated approximation (Bryson's number is specifically
for the Kepler G-dwarf sample; treating it as representative for K/M/F
bins too is a simplification, not a re-derivation for those bins).
Posterior sampled with a hand-rolled Metropolis-Hastings sampler
(no external probabilistic-programming dependency required, and fully
inspectable/explainable end to end, per the plan's suggested fallback).

IMPORTANT CAVEAT carried through to the paper's Discussion: the
"number of systems in the catalog per spectral bin" used below as the
denominator is a stand-in for "number of stars searched with
sufficient sensitivity to detect an HZ planet," NOT the true survey
denominator (which would require target-list and detection-efficiency
data this catalog does not provide). The resulting eta_b values are
therefore illustrative of the *framework*, not a publication-ready
occurrence rate on their own -- exactly the kind of hedge the plan's
Tier 3 section requires stating explicitly rather than omitting.
"""

import csv
import numpy as np

from rc_model import calibrate_kappa, find_olr_limit
from tier2_hz_classification import (
    KOPP_MOIST_GREENHOUSE, KOPP_MAX_GREENHOUSE, kopparapu_seff, spectral_bin,
)
from stellar_trend import (
    S0, TEFF_VALID_RANGE, calibrate_lambda0, rayleigh_albedo,
)

RNG = np.random.default_rng(20260711)  # fixed seed for reproducibility
N_MC = 3000  # Monte Carlo draws per planet, Layer A

# lambda_0 is calibrated once in main() (via stellar_trend.calibrate_lambda0)
# and then used to build a fast interpolation table for rayleigh_albedo(Teff),
# since the exact function does a scipy.integrate.quad call per Teff, far too
# slow to call ~17 million times in the Monte Carlo loop below. This is a
# speed optimization only -- it changes nothing about the physics, and the
# interpolation grid spans exactly TEFF_VALID_RANGE, matching the range
# enforced everywhere else in this project (out-of-range hosts, e.g. hot
# subdwarfs, are excluded in load_catalog_full() below, not extrapolated to).
_LAMBDA_0 = None
_TEFF_GRID = np.linspace(TEFF_VALID_RANGE[0], TEFF_VALID_RANGE[1], 200)
_ALBEDO_GRID = None


def _init_fast_lookup(olr_limit):
    global _LAMBDA_0, _ALBEDO_GRID
    _LAMBDA_0 = calibrate_lambda0(olr_limit)
    _ALBEDO_GRID = np.array([rayleigh_albedo(t, _LAMBDA_0) for t in _TEFF_GRID])


def our_seff_inner_fast(Teff, olr_limit):
    A = np.interp(Teff, _TEFF_GRID, _ALBEDO_GRID)
    return 4.0 * olr_limit / (S0 * (1.0 - A))


def load_catalog_full(path):
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    rows = []
    n_out_of_range = 0
    for row in reader:
        try:
            teff = float(row["st_teff"])
            insol = float(row["pl_insol"])
        except (ValueError, KeyError):
            continue
        if not (TEFF_VALID_RANGE[0] <= teff <= TEFF_VALID_RANGE[1]):
            n_out_of_range += 1
            continue

        def ferr(key):
            try:
                v = float(row[key])
                return v if np.isfinite(v) else 0.0
            except (ValueError, KeyError):
                return 0.0

        rows.append({
            "pl_name": row["pl_name"],
            "st_teff": teff,
            "st_teff_err1": ferr("st_tefferr1"),
            "st_teff_err2": ferr("st_tefferr2"),
            "pl_insol": insol,
            "pl_insol_err1": ferr("pl_insolerr1"),
            "pl_insol_err2": ferr("pl_insolerr2"),
        })
    if n_out_of_range:
        print(f"Excluded {n_out_of_range} planets with host Teff outside "
              f"the valid {TEFF_VALID_RANGE} K range (e.g. hot subdwarfs).")
    return rows


def split_normal_sample(value, err1, err2, n, rng):
    """Two-piece ('split') normal: draws use sigma=err1 (upper, >=0)
    above the central value and sigma=|err2| (lower, err2<=0) below it.
    Falls back to a small fixed relative sigma if both errors are zero
    or missing (avoids a degenerate delta-function posterior)."""
    sig_hi = err1 if err1 > 0 else abs(value) * 0.05
    sig_lo = abs(err2) if err2 < 0 else abs(value) * 0.05
    z = rng.standard_normal(n)
    draws = np.where(z >= 0, value + z * sig_hi, value + z * sig_lo)
    return draws


def layer_a_mc(rows, olr_limit, n_mc=N_MC, rng=RNG):
    """Per-planet Monte Carlo over measurement error -> P(HZ) under
    each boundary model, holding the boundary MODEL fixed per draw."""
    results = []
    for r in rows:
        teff_draws = split_normal_sample(
            r["st_teff"], r["st_teff_err1"], r["st_teff_err2"], n_mc, rng)
        teff_draws = np.clip(teff_draws, TEFF_VALID_RANGE[0], TEFF_VALID_RANGE[1])
        insol_draws = split_normal_sample(
            r["pl_insol"], r["pl_insol_err1"], r["pl_insol_err2"], n_mc, rng)
        insol_draws = np.clip(insol_draws, 0.0, None)

        seff_in_kopp = kopparapu_seff(teff_draws, KOPP_MOIST_GREENHOUSE)
        seff_out_kopp = kopparapu_seff(teff_draws, KOPP_MAX_GREENHOUSE)
        seff_in_ours = our_seff_inner_fast(teff_draws, olr_limit)

        hz_kopp = (insol_draws <= seff_in_kopp) & (insol_draws >= seff_out_kopp)
        hz_ours = (insol_draws <= seff_in_ours) & (insol_draws >= seff_out_kopp)

        results.append({
            "pl_name": r["pl_name"],
            "st_teff": r["st_teff"],
            "bin": spectral_bin(r["st_teff"]),
            "p_hz_kopp": hz_kopp.mean(),
            "p_hz_ours": hz_ours.mean(),
        })
    return results


def layer_b_model_spread(layer_a_results):
    """Model-choice uncertainty: for each planet, whether the
    (measurement-marginalized) HZ probability differs materially
    between the two boundary choices. Kept as its OWN quantity, never
    averaged together with Layer A into one number."""
    out = []
    for r in layer_a_results:
        spread = abs(r["p_hz_kopp"] - r["p_hz_ours"])
        out.append({**r, "model_spread": spread})
    return out


# ---------------------------------------------------------------------
# Simplified hierarchical Bayesian occurrence-rate model
# ---------------------------------------------------------------------
def beta_params_from_mean_ci(mean, ci_lo, ci_hi):
    """Method-of-moments-ish fit of a Beta(alpha,beta) to a stated mean
    and an approximate 95% CI width, used to build a weakly informative
    completeness prior from Bryson et al. (2021)'s eta_Earth =
    0.37 (+0.48/-0.22, 95% CI)."""
    half_width = 0.5 * (ci_hi - ci_lo)
    sd = half_width / 1.96
    var = sd**2
    common = mean * (1 - mean) / var - 1
    alpha = max(mean * common, 0.5)
    beta = max((1 - mean) * common, 0.5)
    return alpha, beta


COMPLETENESS_PRIOR_ALPHA, COMPLETENESS_PRIOR_BETA = beta_params_from_mean_ci(
    0.37, 0.37 - 0.22, 0.37 + 0.48)


def log_prior(eta, comp):
    if not (0.0 < eta < 1.0) or not (0.0 < comp < 1.0):
        return -np.inf
    lp = 0.0
    # eta: flat Beta(1,1) prior
    lp += 0.0
    # completeness: Beta(alpha,beta) prior anchored to Bryson et al. 2021
    a, b = COMPLETENESS_PRIOR_ALPHA, COMPLETENESS_PRIOR_BETA
    lp += (a - 1) * np.log(comp) + (b - 1) * np.log(1 - comp)
    return lp


def log_likelihood(eta, comp, x_obs, n_total):
    p = np.clip(eta * comp, 1e-9, 1 - 1e-9)
    # generalized (real-valued) binomial log-likelihood via log-gamma,
    # to accommodate the fractional "effective count" x_obs from Layer A
    from scipy.special import gammaln
    n, k = n_total, x_obs
    log_binom_coeff = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    return log_binom_coeff + k * np.log(p) + (n - k) * np.log(1 - p)


def metropolis_hastings(x_obs, n_total, n_iter=20000, step=0.03, seed=0):
    rng = np.random.default_rng(seed)
    eta, comp = 0.3, 0.4  # initial state
    chain = np.zeros((n_iter, 2))
    cur_logp = log_prior(eta, comp) + log_likelihood(eta, comp, x_obs, n_total)
    n_accept = 0
    for i in range(n_iter):
        eta_new = eta + rng.normal(0, step)
        comp_new = comp + rng.normal(0, step)
        new_logp = log_prior(eta_new, comp_new)
        if np.isfinite(new_logp):
            new_logp += log_likelihood(eta_new, comp_new, x_obs, n_total)
        if np.log(rng.uniform()) < new_logp - cur_logp:
            eta, comp, cur_logp = eta_new, comp_new, new_logp
            n_accept += 1
        chain[i] = (eta, comp)
    burn = n_iter // 4
    accept_rate = n_accept / n_iter
    return chain[burn:], accept_rate


def summarize_chain(chain):
    eta_samples = chain[:, 0]
    lo, med, hi = np.percentile(eta_samples, [2.5, 50, 97.5])
    return med, lo, hi


def main(catalog_path):
    kappa = calibrate_kappa(target_olr=282.0)
    _, _, olr_limit, _ = find_olr_limit(kappa)
    _init_fast_lookup(olr_limit)
    print(f"Tier-1 OLR limit in use: {olr_limit:.2f} W/m^2\n")

    rows = load_catalog_full(catalog_path)
    print(f"Loaded {len(rows)} planets with valid st_teff & pl_insol")
    print(f"Running Layer A Monte Carlo ({N_MC} draws/planet)...\n")

    layer_a = layer_a_mc(rows, olr_limit)
    layer_ab = layer_b_model_spread(layer_a)

    # --- Layer A + B summary by spectral bin ---
    bins = ["M (<3900K)", "K (3900-5300K)", "G (5300-6000K)", "F/A (>=6000K)"]
    print("Layer A (measurement uncertainty) + Layer B (model spread), by bin:")
    print(f"{'Bin':16s} {'N':>6s} {'E[N_HZ|Kopp]':>13s} {'E[N_HZ|ours]':>13s} {'mean|spread|':>13s}")
    bin_stats = {}
    for b in bins:
        sub = [r for r in layer_ab if r["bin"] == b]
        n = len(sub)
        e_kopp = sum(r["p_hz_kopp"] for r in sub)
        e_ours = sum(r["p_hz_ours"] for r in sub)
        mean_spread = np.mean([r["model_spread"] for r in sub]) if sub else 0.0
        bin_stats[b] = dict(n=n, e_kopp=e_kopp, e_ours=e_ours)
        print(f"{b:16s} {n:6d} {e_kopp:13.2f} {e_ours:13.2f} {mean_spread:13.4f}")

    # highest per-planet model spread (most "boundary-sensitive" planets)
    top_spread = sorted(layer_ab, key=lambda r: -r["model_spread"])[:5]
    print("\nMost boundary-sensitive planets (largest Layer-B model spread):")
    for r in top_spread:
        print(f"  {r['pl_name']:20s} Teff={r['st_teff']:7.1f}  "
              f"P(HZ|Kopp)={r['p_hz_kopp']:.3f}  P(HZ|ours)={r['p_hz_ours']:.3f}  "
              f"spread={r['model_spread']:.3f}")

    # --- Bayesian hierarchical occurrence rate, per bin, per boundary choice ---
    print("\nSimplified hierarchical Bayesian occurrence-rate estimate")
    print(f"(completeness prior: Beta(alpha={COMPLETENESS_PRIOR_ALPHA:.2f}, "
          f"beta={COMPLETENESS_PRIOR_BETA:.2f}), anchored to Bryson et al. 2021)\n")
    print(f"{'Bin':16s} {'model':8s} {'eta_med':>9s} {'95% CI':>18s} {'accept%':>8s}")
    for b in bins:
        n_total = bin_stats[b]["n"]
        if n_total == 0:
            continue
        for label, key in [("Kopparapu", "e_kopp"), ("ours", "e_ours")]:
            x_obs = bin_stats[b][key]
            chain, acc = metropolis_hastings(x_obs, n_total, seed=hash((b, label)) % (2**32))
            med, lo, hi = summarize_chain(chain)
            print(f"{b:16s} {label:8s} {med:9.3f} [{lo:.3f}, {hi:.3f}]     {100*acc:6.1f}")

    print("\nCaveat (state in Discussion): the catalog-count denominator used")
    print("here is a stand-in for true survey completeness and is NOT a")
    print("publication-ready occurrence rate on its own -- see module docstring.")


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "catalog.csv"
    main(path)
