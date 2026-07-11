# Sensitivities of Exoplanet Habitable Zone Catalogs to Climate Model Selection

This repository contains the data, models, and manuscript for the paper "Sensitivities of Exoplanet Habitable Zone Catalogs to Climate Model Selection" (Adhikary, 2026).

## Physics Insight: What We Did
In exoplanet science, the "Habitable Zone" (HZ) boundaries are usually computed using a single theoretical standard (often the Kopparapu et al. polynomial limits). However, any climate model carries deep structural uncertainty—how we parameterize water vapor absorption, stratospheric temperatures, or scattering profiles fundamentally alters the outgoing longwave radiation limit (the "moist greenhouse" runaway boundary).

We built an independent, first-principles 1-D radiative-convective climate model to explicitly compute these boundaries. Rather than just making another discrete standard, we ran a **"perturbed-physics" Monte Carlo ensemble**. We treated our internal physical constants as continuous probability distributions and propagated them directly against real catalog data from the NASA Exoplanet Archive.

**The result:** The HZ catalog demographics are incredibly robust! Despite significant variance in core physical parameters (e.g., 10% uncertainty in the water-vapor continuum, 20-25% range in stratospheric temperatures and Rayleigh scattering), the number of candidate habitable planets fluctuates by only roughly $\pm 6\%$. This mathematically proves that mission planners for future telescopes like HWO and LIFE can safely rely on simple functional boundaries without severely underestimating systematic implementation errors in their target yields.

## Repository Structure

- `Draft_Paper.md` : The complete manuscript, formatted for submission to *Icarus*.
- `PSCompPars_*.csv` : The frozen, reproducible snapshot of the NASA Exoplanet Archive used in this study.
- `hz_candidates_comparison.csv` : The derived output catalog tracking exactly which planets fell into which theoretical bounds.
- `tier1_model/` : The Python module containing the 1-D radiative-convective model, root-finding integration, and Monte Carlo engine.
- `*.py` scripts : Reproduction scripts for generating the surrogate grid, running the sensitivity audits, and plotting the figures.
- `*.png` : The generated figures used directly in the manuscript.

## License
The code and analysis scripts in this repository are distributed under the MIT License.
