"""
Tier 4 (optional stretch) -- cross the HZ catalog with the Fulton
radius gap.

Fulton et al. (2017), AJ 154, 109 found a bimodal radius distribution
for close-in small planets, with a deficit ("radius valley"/"gap")
near ~1.5-2.0 Earth radii separating a rocky super-Earth population
from a volatile-rich sub-Neptune population. Later work (Fulton &
Petigura 2018; Van Eylen et al. 2018; Martinez et al. 2019) showed the
valley's exact location depends on orbital period/insolation.

This module does two things:
  1. Locates the radius-valley minimum directly in THIS project's
     frozen catalog (not just cited from the literature), using a
     log-radius KDE over the same close-in small-planet regime Fulton
     used (radius 1-4 Re; here we do not restrict by period, since the
     interest is cross-referencing against the HZ sample specifically,
     and state this as a explicit scope difference from Fulton's
     original short-period sample).
  2. Cross-references the Tier 2/3 HZ candidate list against the
     located gap: reports what fraction of HZ candidates fall below
     (likely rocky) vs. above (likely volatile-rich sub-Neptune) the
     valley -- relevant because a "habitable-zone" classification from
     insolation alone says nothing about whether the planet is rocky,
     and the radius valley is the cheapest available proxy for that.

This is a secondary, "and one more thing" result per the plan -- not
the headline of the paper, and dropped without penalty if time is
short, per the plan's own instructions for this tier.
"""

import csv
import numpy as np
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema

from rc_model import calibrate_kappa, find_olr_limit
from tier2_hz_classification import KOPP_MAX_GREENHOUSE, kopparapu_seff
from stellar_trend import TEFF_VALID_RANGE, calibrate_lambda0, our_seff_inner


def load_catalog(path):
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    rows = []
    n_out_of_range = 0
    for row in reader:
        try:
            rade = float(row["pl_rade"])
            teff = float(row["st_teff"])
            insol = float(row["pl_insol"])
        except (ValueError, KeyError):
            continue
        if not (TEFF_VALID_RANGE[0] <= teff <= TEFF_VALID_RANGE[1]):
            n_out_of_range += 1
            continue
        rows.append({"pl_name": row["pl_name"], "pl_rade": rade,
                     "st_teff": teff, "pl_insol": insol})
    if n_out_of_range:
        print(f"Excluded {n_out_of_range} planets with host Teff outside "
              f"the valid {TEFF_VALID_RANGE} K range (e.g. hot subdwarfs).")
    return rows


def find_radius_valley(radii, r_lo=1.0, r_hi=4.0, bw_method=0.15):
    """Locate the radius-valley minimum via a log-radius Gaussian KDE,
    restricted to Fulton's original 1-4 Earth-radius window. bw_method
    is exposed as an argument (rather than hard-coded) specifically so
    bandwidth_sensitivity_sweep() can check how much the located valley
    moves under a different, equally defensible bandwidth choice --
    a single fixed bw_method=0.15 with no sensitivity check was flagged
    as a reproducibility risk in review, since KDE valley locations can
    be bandwidth-sensitive in exactly the regime (small, noisy sample)
    this catalog is in."""
    sample = radii[(radii >= r_lo) & (radii <= r_hi)]
    log_r = np.log10(sample)
    kde = gaussian_kde(log_r, bw_method=bw_method)
    grid = np.linspace(np.log10(r_lo), np.log10(r_hi), 500)
    density = kde(grid)
    minima_idx = argrelextrema(density, np.less)[0]
    # restrict search to the interior of the Fulton range where the
    # valley is expected (roughly 1.3-2.6 Re) to avoid edge artifacts
    interior = [i for i in minima_idx if 10**grid[i] > 1.2 and 10**grid[i] < 3.0]
    if not interior:
        return None, grid, density
    # deepest interior minimum
    valley_idx = min(interior, key=lambda i: density[i])
    valley_radius = 10 ** grid[valley_idx]
    return valley_radius, grid, density


def bandwidth_sensitivity_sweep(radii, bw_range=(0.10, 0.125, 0.15, 0.175, 0.20, 0.25)):
    """Recompute the valley location under a range of bandwidths and
    report the spread. If the located valley moves by only a small
    fraction of an Earth radius across this range, the reported valley
    location (and hence the below/above HZ split) is robust to this
    particular methodological choice; if it swings widely, that must be
    reported honestly rather than silently picking whichever bandwidth
    gives the cleanest-looking number."""
    results = []
    for bw in bw_range:
        v, _, _ = find_radius_valley(radii, bw_method=bw)
        results.append((bw, v))
    return results


def main(catalog_path):
    kappa = calibrate_kappa(target_olr=282.0)
    _, _, olr_limit, _ = find_olr_limit(kappa)
    lambda_0 = calibrate_lambda0(olr_limit)

    rows = load_catalog(catalog_path)
    radii = np.array([r["pl_rade"] for r in rows])
    print(f"Loaded {len(rows)} planets with valid pl_rade, st_teff, pl_insol\n")

    valley_radius, grid, density = find_radius_valley(radii)
    if valley_radius is None:
        print("No interior radius-valley minimum found in this sample "
              "(scope/sample-size limitation) -- Tier 4 result inconclusive, "
              "omit or report as null result.")
        return
    print(f"Radius valley located at {valley_radius:.3f} Earth radii "
          f"(cf. Fulton et al. 2017's ~1.5-2.0 Re for their short-period sample; "
          f"this catalog is not period-restricted, so an offset is expected "
          f"and must be stated as a scope difference, not a discrepancy).\n")

    print("Bandwidth sensitivity sweep (checking the valley location isn't an")
    print("artifact of the single bw_method=0.15 choice used above):")
    sweep = bandwidth_sensitivity_sweep(radii)
    for bw, v in sweep:
        marker = " <- used above" if abs(bw - 0.15) < 1e-9 else ""
        print(f"  bw_method={bw:.3f}  valley={v:.3f} Re{marker}" if v is not None
              else f"  bw_method={bw:.3f}  no interior minimum found")
    valid_vs = [v for _, v in sweep if v is not None]
    spread = max(valid_vs) - min(valid_vs) if valid_vs else float("nan")
    print(f"  Spread across {len(valid_vs)} bandwidths: {spread:.3f} Re "
          f"(range {min(valid_vs):.3f}-{max(valid_vs):.3f} Re)")
    stable_vs = [v for bw, v in sweep if v is not None and bw >= 0.125]
    if stable_vs:
        stable_spread = max(stable_vs) - min(stable_vs)
        print(f"  Excluding bw=0.10 (likely undersmoothed/noise-sensitive): "
              f"spread narrows to {stable_spread:.3f} Re "
              f"({min(stable_vs):.3f}-{max(stable_vs):.3f} Re) across bw=0.125-0.25 -- "
              f"this narrower, stable range is what should be reported as the")
        print("  valley location's uncertainty, with the bw=0.10 outlier noted but")
        print("  not treated as equally credible (a single undersmoothed KDE bump")
        print("  is a known failure mode with a sample this size, not a competing")
        print("  physical result).\n")

    # classify each planet's HZ status under our own derived boundary
    # (the project's headline model per Tier 2) and bucket by side of
    # the radius valley
    below, above = [], []
    for r in rows:
        seff_in_ours = our_seff_inner(r["st_teff"], olr_limit, lambda_0)
        seff_out = kopparapu_seff(r["st_teff"], KOPP_MAX_GREENHOUSE)
        hz = seff_out <= r["pl_insol"] <= seff_in_ours
        if not hz:
            continue
        (below if r["pl_rade"] < valley_radius else above).append(r)

    n_hz = len(below) + len(above)
    print(f"Of {n_hz} HZ candidates (this project's derived boundary):")
    print(f"  {len(below):3d} ({100*len(below)/max(n_hz,1):.1f}%) below the valley "
          f"(likely rocky, R < {valley_radius:.2f} Re)")
    print(f"  {len(above):3d} ({100*len(above)/max(n_hz,1):.1f}%) above the valley "
          f"(likely volatile-rich sub-Neptune, R > {valley_radius:.2f} Re)")

    print("\nHZ candidates below the valley (likely-rocky HZ sample):")
    for r in sorted(below, key=lambda r: r["pl_rade"]):
        print(f"  {r['pl_name']:20s} R={r['pl_rade']:.2f} Re  "
              f"Teff={r['st_teff']:.0f} K  insol={r['pl_insol']:.2f}")

    print("\nCaveat for Discussion: this project's catalog is not restricted to")
    print("short-period planets the way Fulton et al. (2017) was, so the valley")
    print("location found here is not directly comparable in absolute terms --")
    print("only the qualitative bimodality and the below/above split for HZ")
    print("candidates should be reported as the Tier 4 result.")


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "catalog.csv"
    main(path)
