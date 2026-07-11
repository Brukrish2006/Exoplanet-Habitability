"""
Tier 2 -- HZ model-choice sensitivity on the frozen catalog.

Classifies every catalog planet's HZ membership TWICE:
  (a) using the standard Kopparapu et al. (2013/2014) polynomial fits
      for both the inner (moist-greenhouse) and outer (maximum-
      greenhouse) edges -- the literature baseline everyone else uses;
  (b) using this project's own Tier-1-derived inner edge (rc_model.py +
      stellar_trend.py's coupling-efficiency scaling) in place of the
      literature moist-greenhouse boundary, while keeping Kopparapu's
      outer (maximum-greenhouse) edge unchanged, since Tier 1 only
      derived the inner edge.

Reports how many planets flip HZ membership depending on which inner
boundary is trusted, broken down by host spectral type -- the headline
Tier 2 result per the research plan.
"""

import csv
import numpy as np
from rc_model import calibrate_kappa, find_olr_limit
from stellar_trend import (
    KOPP_MOIST_GREENHOUSE, KOPP_MAX_GREENHOUSE, kopparapu_seff,
    calibrate_lambda0, TEFF_VALID_RANGE,
    our_seff_inner as _our_seff_inner_raylegh,
)

# Kopparapu coefficients and the project's own Seff model now live in
# stellar_trend.py as the single source of truth (previously duplicated
# here with risk of drifting out of sync -- consolidated after review).
# our_seff_inner here just fixes lambda_0 once (calibrated in main())
# and exposes a 2-arg wrapper matching this module's original call
# signature, so the rest of this file doesn't need to change.
_LAMBDA_0 = None


def our_seff_inner(Teff, olr_limit):
    return _our_seff_inner_raylegh(Teff, olr_limit, _LAMBDA_0)


def spectral_bin(Teff):
    if Teff >= 6000:
        return "F/A (>=6000K)"
    elif Teff >= 5300:
        return "G (5300-6000K)"
    elif Teff >= 3900:
        return "K (3900-5300K)"
    else:
        return "M (<3900K)"


def load_catalog(path):
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
            # host outside the Kopparapu fit's (and this project's own
            # Rayleigh model's) valid Teff domain -- e.g. hot subdwarfs --
            # excluded rather than extrapolated to; see stellar_trend.py
            n_out_of_range += 1
            continue
        rows.append({"pl_name": row["pl_name"], "st_teff": teff,
                     "pl_insol": insol})
    if n_out_of_range:
        print(f"Excluded {n_out_of_range} planets with host Teff outside "
              f"the valid {TEFF_VALID_RANGE} K range (e.g. hot subdwarfs).")
    return rows


def main(catalog_path):
    global _LAMBDA_0
    kappa = calibrate_kappa(target_olr=282.0)
    _, _, olr_limit, _ = find_olr_limit(kappa)
    _LAMBDA_0 = calibrate_lambda0(olr_limit)
    print(f"Using Tier-1 OLR limit = {olr_limit:.2f} W/m^2\n")

    rows = load_catalog(catalog_path)
    print(f"Loaded {len(rows)} planets with valid st_teff and pl_insol\n")

    counts = {}  # bin -> [n_total, n_flip, n_kopp_hz, n_custom_hz]
    flips = []

    for r in rows:
        Teff, insol = r["st_teff"], r["pl_insol"]

        seff_in_kopp = kopparapu_seff(Teff, KOPP_MOIST_GREENHOUSE)
        seff_out_kopp = kopparapu_seff(Teff, KOPP_MAX_GREENHOUSE)
        seff_in_ours = our_seff_inner(Teff, olr_limit)

        hz_kopp = seff_out_kopp <= insol <= seff_in_kopp
        hz_ours = seff_out_kopp <= insol <= seff_in_ours

        b = spectral_bin(Teff)
        c = counts.setdefault(b, [0, 0, 0, 0])
        c[0] += 1
        c[2] += int(hz_kopp)
        c[3] += int(hz_ours)
        if hz_kopp != hz_ours:
            c[1] += 1
            flips.append((r["pl_name"], Teff, insol, hz_kopp, hz_ours))

    print(f"{'Spectral bin':18s} {'N':>6s} {'N_flip':>7s} {'N_HZ(Kopp)':>11s} {'N_HZ(ours)':>11s}")
    tot = [0, 0, 0, 0]
    for b in ["M (<3900K)", "K (3900-5300K)", "G (5300-6000K)", "F/A (>=6000K)"]:
        if b in counts:
            n, nf, nk, no = counts[b]
            print(f"{b:18s} {n:6d} {nf:7d} {nk:11d} {no:11d}")
            for i in range(4):
                tot[i] += counts[b][i]
    print(f"{'TOTAL':18s} {tot[0]:6d} {tot[1]:7d} {tot[2]:11d} {tot[3]:11d}")

    print(f"\n{len(flips)} planets flip HZ classification depending on model choice.")
    print("First 10 examples (name, Teff, insol[Earth flux], HZ_Kopparapu, HZ_ours):")
    for f in flips[:10]:
        print(" ", f)

    return rows, counts, flips


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "catalog.csv"
    main(path)
