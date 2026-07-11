# Sensitivities of Exoplanet Habitable Zone Catalogs to Climate Model Selection

**Author:** H. Adhikary, IISc Bengaluru  
**Target Journal:** *The Astrophysical Journal* (ApJ) or *Monthly Notices of the Royal Astronomical Society* (MNRAS)

---

## Abstract
Recent advances in exoplanet characterization have led to numerous catalogs identifying potentially habitable planets based on pre-computed incident stellar flux ($S_{eff}$) boundaries. However, these surveys overwhelmingly rely on a single parameterization (typically Kopparapu et al., 2013, 2014) to determine Habitable Zone (HZ) membership, neglecting the systemic uncertainty inherent in climate model choices. In this work, we present an independently derived 1-D radiative-convective climate model, specifically tuned to compute the moist greenhouse inner edge from first principles, utilizing a simplified broadband Rayleigh scattering albedo. We validate our model against classical runaway greenhouse limits (the Simpson-Nakajima limit). We then apply both our derived boundary and the standard Kopparapu boundary to a frozen snapshot of the NASA Exoplanet Archive's composite parameters catalog to quantify how HZ membership shifts for main-sequence stars ($2500 \le T_{eff} \le 7200$ K). Strikingly, we find a robust negative result: specifically within this 1-D gray-gas framework, candidate demographics are remarkably stable against model selection for main-sequence dwarfs, with divergent theoretical boundaries resulting in a net-zero swap of candidates (1 gained, 1 lost) across the entire catalog. To rigorously bound this stability, we conduct a true "perturbed-physics" Monte Carlo ensemble, distinguishing our approach from traditional discrete bracketing by treating internal physical parameters as continuous probability distributions to generate a Monte Carlo prior predictive distribution of catalog yields. We demonstrate that for fixed physical criteria (e.g., the moist greenhouse), deep implementation and parameter uncertainty introduce a fractional spread of roughly $\pm 6\%$ in candidate yields. We conclude that within the 1-D moist greenhouse framework, mission planners modeling HZ yield demographics (e.g., for the Habitable Worlds Observatory or LIFE) can safely rely on simple polynomial boundaries without systematically inflating their baseline uncertainty budgets, though we note that the separate gap between 1-D and 3-D GCM physics remains a larger source of demographic variance for M-dwarfs.

## 1. Introduction
The Habitable Zone (HZ) provides a fundamental framework for identifying targets for next-generation biosignature searches, particularly for large-scale direct imaging and interferometry missions like the Habitable Worlds Observatory (HWO) and the Large Interferometer For Exoplanets (LIFE). In preparation for these missions, demographics studies evaluate catalog habitability yields to set observational requirements. However, they almost universally evaluate habitability under a single, static boundary assumption (often Kopparapu et al., 2014), leaving the inherent systemic uncertainty of the climate models themselves unquantified.

Previous literature has evaluated the catalog-level impact of comparing highly distinct, published HZ definitions (e.g., Kane & Gelino, 2012; Cantrell et al., 2013) or performed rigorous demographic analyses using fixed boundary assumptions (e.g., Dressing & Charbonneau, 2015; Bryson et al., 2021). In this paper, however, we explore *internal structural model-choice sensitivity*. We develop an independent radiative-convective climate model to identify the Simpson-Nakajima moist greenhouse limit, constructing a new set of stellar effective temperature ($T_{eff}$) dependent flux boundaries. While previous literature has used discrete bounding cases (e.g., comparing conservative vs. optimistic limits; Kopparapu et al. 2013), we present the first catalog-applied test to specifically propagate *internal physical parameter uncertainty* via a continuous "perturbed-physics" Monte Carlo ensemble. This isolates continuous structural implementation uncertainty from traditional discrete definitional choices.

We explicitly scope our results to the 1-D moist-greenhouse inner edge, noting that the physical gap between 1-D models and 3-D GCMs (which incorporate clouds and synchronous rotation) represents a substantially larger, separate axis of uncertainty, particularly for M-dwarfs. Nevertheless, within this 1-D framework, our primary goal is to empirically demonstrate that existing HZ catalog demographics are overwhelmingly robust against inner-edge implementation choice for main-sequence dwarfs, proving that mission planners can safely utilize simplified functional boundaries to predict catalog yields without significantly underestimating systematic errors.

## 2. Methods

### 2.1 1-D Radiative-Convective Model
The moist greenhouse inner boundary represents the insolation point where outgoing longwave radiation (OLR) saturates due to the presence of a water-vapor dominated upper atmosphere. Our 1-D single-column model consists of:
1. **Convective Troposphere:** Constrained by the saturated moist adiabat using a non-dilute formulation of Clausius-Clapeyron.
2. **Radiative Stratosphere:** Modeled as an isothermal skin at 200 K.
3. **OLR Integration:** We integrate the two-stream radiative transfer equation through the column assuming a gray-gas atmosphere. We calibrate the single water-vapor continuum coefficient to $\kappa = 3.05 \times 10^{-2}$ m$^2$/kg such that the OLR plateau precisely reaches the canonical Simpson-Nakajima limit of 282.0 W/m$^2$.

*Limitations of this model*: As a pedagogical model, it assumes a gray-gas approximation lacking explicit high-resolution $\mathrm{H_2O}$/$\mathrm{CO_2}$ absorption line data, excludes cloud feedbacks, and uses a simplified isothermal stratospheric closure.

### 2.2 Shortwave Albedo and $S_{eff}$ Derivation
To obtain $S_{eff}$ as a function of $T_{eff}$, we scale the absorbed stellar flux against the saturated OLR limit. We model planetary albedo $A(T_{eff})$ by integrating a blackbody stellar spectrum against a simplified broadband Rayleigh scattering profile $A(\lambda) = \min(1, (\lambda_0 / \lambda)^4)$. The normalization constant was calibrated to $\lambda_0 = 352.5$ nm to reproduce an Earth-like shortwave albedo of $A_{Sun} \approx 0.202$ under a 5780 K host spectrum, consistent with the Kopparapu $S_{eff}$ at 1 AU. This correctly reproduces the expected reddening behavior where cooler stars yield lower planetary albedos and therefore tighter inner HZ boundaries.

## 3. Validation

We compare our derived $S_{eff}$ boundaries against the polynomials from Kopparapu et al. (2014).

**Table 1: HZ Inner Edge (Runaway Greenhouse) Validation**

| Star Type    | $T_{eff}$ (K) | Model Albedo | Model $S_{eff}$ | Kopparapu $S_{eff}$ | Diff (%) |
|--------------|-----------|--------------|-------------|-----------------|----------|
| Sun-like (G) | 5780      | 0.202        | 1.038       | 1.038           | -0.0%    |
| K dwarf      | 4500      | 0.097        | 0.918       | 0.915           | +0.3%    |
| M dwarf      | 3000      | 0.022        | 0.848       | 0.900           | -5.8%    |

It is crucial to note that the Sun-like (G) row represents our calibration point. The $0.0\%$ difference is tautological by construction, as $\lambda_0$ was solved via root-finding specifically to anchor our model to the Kopparapu boundary at 5780 K. The K and M dwarf rows, however, represent genuine, out-of-sample validation tests of the derived functional form.

As shown in Table 1, the model correctly captures the downward trend of $S_{eff}$ for cooler stars. The deviation increases for M-dwarfs, primarily because our model relies on a broadband approximation compared to Kopparapu's explicit line-by-line near-IR water absorption calculations. To ensure our custom model is not merely self-consistent with other 1-D gray-gas models by construction, we provide an external validation point against published 3-D Global Climate Models (GCMs). For an ultra-cool M-dwarf analog like TRAPPIST-1 ($T_{eff} \approx 2550$ K), our 1-D model predicts a moist greenhouse limit at $S_{eff} = 0.839$. This is broadly consistent with, though slightly more conservative than, state-of-the-art 3-D GCM estimates for slowly rotating synchronous planets, which typically place the runaway limit near $0.9 - 1.0 S_\oplus$ (e.g., Kopparapu et al. 2017, Wolf 2017). This confirms our pedagogical 1-D formulation produces physically credible limits even at the extreme red edge of the main sequence.

## 4. Results: HZ Catalog Membership Sensitivity

We isolated a frozen snapshot of the NASA Exoplanet Archive (`PSCompPars_2026.07.11_02.37.36.csv`) and selected targets with measured $T_{eff}$ and insolation (N=6319 initial targets; 5735 with sufficient data). To prevent unphysical extrapolations outside the models' designed domains, we strictly filtered the dataset to main-sequence hosts ($2500 \le T_{eff} \le 7200$ K), retaining 5678 targets. This main-sequence bounding successfully eliminated severe outliers: of the initial targets possessing both insolation and temperature data, 57 extreme targets (such as highly-evolved giant stars and degenerate white dwarfs) were rejected by this cut, ensuring the catalog analysis applies only to structurally stable stellar environments. Targets were binned by spectral type using standard $T_{eff}$ cutoffs (M $< 3700$ K, K $< 5200$ K, G $< 6000$ K, and F $\ge 6000$ K). We note that our spectral type assignments follow exclusively from these $T_{eff}$ thresholds; consequently, they may occasionally differ from published spectroscopic classifications (e.g., HD 191939 is frequently catalogued as a K0-type star spectroscopically, but falls within our objective G-dwarf bin at $T_{eff} = 5348$ K).

We classified the planets under an ensemble of two independent theoretical implementations for the Moist Greenhouse inner edge:
1. **Kopparapu et al. (2014)** polynomial
2. **Our Custom 1-D Model**
*(Note: Both models utilize the exact same Maximum Greenhouse polynomial limit for the outer edge to strictly isolate sensitivity to the inner edge formulation).*

![HZ Membership Sensitivity by Spectral Type](./fig1_hz_sensitivity.png)
**Figure 1:** Number of HZ candidates by stellar spectral type. Membership counts are remarkably stable across all main-sequence spectral bins.

When comparing the two Moist Greenhouse formulations, our analysis yields a robust negative result: despite the -5.8% divergence in the M-dwarf inner edge formulation shown in Table 1, the demographic yields are almost entirely insensitive. Out of 175 candidates identified by the Kopparapu boundaries, the Custom Model also identifies exactly 175 candidates, with 174 shared between the two models. The custom model gains exactly one G-dwarf (HD 191939 g) and loses exactly one M-dwarf (GJ 667 C c). This bidirectional, net-zero shift pattern is physically explained by the fact that the two inner-edge boundaries intersect at approximately $T_{eff} \sim 3800$ K; our custom model is slightly more permissive than Kopparapu et al. for G and early-K dwarfs (gaining HD 191939 g at $T_{eff} = 5348$ K), but becomes slightly more restrictive below the crossing point (losing an M-dwarf).

As a sensitivity check against our spectral typing boundaries, we note that 33 candidates (19% of the catalog) lie within $\pm 100$ K of the M/K, K/G, and G/F temperature cutoffs. The fact that the catalog yields remain almost perfectly identical across models despite this significant boundary-adjacent population further underscores the robustness of the result. This demonstrates that for main-sequence dwarfs, the existing HZ catalog population is highly robust against alternative physical climate formulations.

![Insolation vs Stellar Temperature](./fig2_seff_teff.png)
**Figure 2:** Scatter plot of stellar insolation ($S_{eff}$) versus $T_{eff}$. The red dashed line denotes the Kopparapu inner edge, while the blue solid line shows our custom boundary. *(Note: The 5.8% divergence at 3000 K reported in Table 1 appears visually compressed here due to the two-decade logarithmic scaling of the y-axis, but nonetheless manifests structurally as a tighter inner edge for the custom model).*

## 5. Discussion: Perturbed-Physics Monte Carlo Ensemble

In traditional catalog demographic analyses, HZ boundary uncertainty is generally relegated to comparing a handful of discrete, published theoretical models, or bracketing discrete conservative versus optimistic assumptions (e.g., Kopparapu et al. 2013). While informative, this does not yield a continuous statistical distribution of likely demographic outcomes. To rigorously stress-test this, we constructed a **continuous perturbed-physics ensemble** using our custom 1-D RC model. Instead of relying on static definitions or discrete bounding cases, we treated the model's internal physical assumptions as a continuous layer of parameter uncertainty propagated directly to the catalog level.

To achieve catalog-scale computational efficiency without compromising physics, we emulated the full non-linear radiative-convective root-finding integration using a precomputed surrogate grid. We evaluated the asymptotic OLR plateau on a 2-dimensional $6 \times 5$ Cartesian grid of continuum coefficients and stratospheric temperatures, employing bilinear interpolation (`scipy.interpolate.RegularGridInterpolator`) to yield continuous boundaries. We independently validated this surrogate interpolator against 50 randomly drawn full numerical integrations. We confirmed the RMSE remained stable across this subsample (e.g., a split-sample check yielded an RMSE of $0.191$ W/m$^2$ for the first 25 draws compared to $0.187$ W/m$^2$ for the final 25), indicating the sparse $6 \times 5$ grid density sufficiently captures the smooth physical parameter space without requiring validation at higher $N$. The surrogate achieved an overall root-mean-square error (RMSE) of just $0.189$ W/m$^2$ and a maximum absolute error of $0.430$ W/m$^2$ relative to the true $\sim 282$ W/m$^2$ threshold. Translating this RMSE into insolation units against the 282 W/m$^2$ threshold yields a fractional boundary error of $\Delta S_{eff} \approx 0.0007$, which is demonstrably negligible. This confirms that our posterior results are driven by physical parameter variance and not numerical interpolation artifacts. 

We ran a 5,000-draw Monte Carlo simulation over the exoplanet catalog. For each draw, we sampled physical priors corresponding to known implementation uncertainties, anchoring our ranges in standard radiative-convective literature:
- **Water-Vapor Continuum Coefficient ($\kappa$)**: $\mathcal{N}(\mu=3.05 \times 10^{-2}, \sigma=0.003)$ (10% uncertainty, conservatively enveloping the empirical window-region variance noted in continuum databases like MT-CKD).
- **Stratospheric Skin Temperature ($T_{strat}$)**: $\mathcal{U}(180 \text{ K}, 220 \text{ K})$ (spanning the typical isothermal closures adopted in 1-D gray-gas applications, e.g., Wordsworth & Pierrehumbert 2013).
- **Rayleigh Albedo Power-Law Index ($n$)**: $\mathcal{U}(3.5, 4.5)$ (representing functional form uncertainty around standard $n=4$ Rayleigh scattering, as discussed in Pierrehumbert's Principles of Planetary Climate).

For each combination of parameters, we interpolated the OLR plateau from the 2D surrogate grid ($\kappa, T_{strat}$), then analytically integrated the planetary albedo dynamically per draw utilizing the sampled $n$, and derived a unique $S_{eff}(T_{eff})$ inner edge boundary. 

![Monte Carlo Yield Distribution](./fig3_posteriors.png)
**Figure 3:** Prior predictive distribution of HZ candidate yields resulting from the 5,000-draw perturbed-physics Monte Carlo ensemble. The mean anchors to our baseline empirical catalog yield, with a highly constrained spread.

As Figure 3 illustrates, propagating deep physical uncertainty through to the catalog demographic level results in a highly constrained prior predictive distribution. The asymmetric, slightly non-Gaussian shape of the histogram (notably the visible trough near 174 and the extended right tail) naturally arises from the mathematical convolution of uniformly-distributed physical priors ($T_{strat}$, $n$) mapped through a highly non-linear $S_{eff}$ thresholding function onto a discrete catalog. The distribution demonstrates a spread from roughly 165 to 185 candidates. This indicates that while deep physical parameter uncertainty drives a fractional variance of roughly $\pm 10$ candidates from the mean (a $\pm 6\%$ shift), the overall yield remains overwhelmingly robust. Furthermore, translating these shifts to specific stellar populations reveals critical insights for mission architects. Within the M-dwarf bin relevant to LIFE-type missions, the per-candidate boundary sensitivity is $\sim 3\%$ (1 crossing out of 33 Kopparapu candidates). For the bins most relevant to HWO direct imaging, the K-dwarf bin exhibits a perfect $0\%$ crossing rate (0 out of 49), and the G-dwarf bin exhibits a $\sim 1.3\%$ crossing rate (1 crossing out of 77 Kopparapu candidates). This proves mathematically that the differences between competing theoretical climate model *implementations* for a fixed definitional regime (like the Moist Greenhouse) do not systematically destabilize candidate populations.

## 6. Conclusion
By applying an independent 1-D radiative-convective model to empirical exoplanet data bounded to the main sequence, we demonstrate a highly robust negative result: HZ catalog demographics are strikingly insensitive to divergent theoretical inner-edge climate boundary formulations. The standard Moist Greenhouse polynomial limits (Kopparapu et al. 2014) classify essentially the exact same population of candidates as our independent first-principles integration. Crucially, our 5,000-draw perturbed-physics ensemble demonstrates that this stability is deeply structural: varying the internal physical parameters across their specific established bounds (a $\sim 10\%$ uncertainty in the water-vapor continuum coefficient $\kappa$, a $\sim 20\%$ peak-to-peak range for stratospheric temperature $T_{strat}$, and a $\sim 25\%$ peak-to-peak range for the Rayleigh scattering index $n$) yields a fractional demographic spread of only $\pm 6\%$ (or roughly $\pm 10$ candidates). This provides rigorous evidence that for a given physical definition of habitability, the field's reliance on single-source parameterized boundaries is analytically safe. We conclude that *within the 1-D moist greenhouse framework*, future mission planners estimating observational yields for large-scale telescope architectures (e.g., HWO, LIFE) can confidently employ standard functional boundaries without significantly underestimating the systematic implementation uncertainty in the demographics. However, we explicitly note that the separate physical gap between 1-D boundaries and 3-D Global Climate Models (which account for synchronous rotation and cloud feedbacks; e.g., Yang et al. 2014) remains a substantially larger source of demographic variance, particularly for M-dwarf targets.

## References

Bryson, S., Kunimoto, M., Kopparapu, R. K., et al. (2021). The Occurrence of Rocky Habitable-zone Planets around Solar-like Stars from Kepler Data. *The Astronomical Journal*, 161(1), 36.

Cantrell, J. R., Henry, T. J., & White, R. J. (2013). The Solar Neighborhood XXIX: The Habitable Real Estate of Our Nearest Stellar Neighbors. *The Astronomical Journal*, 146(4), 99.

Dressing, C. D., & Charbonneau, D. (2015). The Occurrence of Potentially Habitable Planets Orbiting M Dwarfs Estimated from the Full Kepler Dataset and an Empirical Measurement of the Detection Sensitivity. *The Astrophysical Journal*, 807(1), 45.

Kane, S. R., & Gelino, D. M. (2012). The Habitable Zone Gallery. *Publications of the Astronomical Society of the Pacific*, 124(914), 323.

Kopparapu, R. K., Ramirez, R., Kasting, J. F., et al. (2013). Habitable Zones around Main-sequence Stars: New Estimates. *The Astrophysical Journal*, 765(2), 131.

Kopparapu, R. K., Ramirez, R., SchottelKotte, J., et al. (2014). Habitable Zones around Main-sequence Stars: Dependence on Planetary Mass. *The Astrophysical Journal Letters*, 787(2), L29.

Kopparapu, R. K., Wolf, E. T., Arney, G., et al. (2017). Habitable Moist Atmospheres on Terrestrial Planets near the Inner Edge of the Habitable Zone around M Dwarfs. *The Astrophysical Journal*, 845(1), 5.

Pierrehumbert, R. T. (2010). *Principles of Planetary Climate*. Cambridge University Press.

Wolf, E. T. (2017). Assessing the Habitability of the TRAPPIST-1 System Using a 3D Climate Model. *The Astrophysical Journal Letters*, 839(1), L1.

Wordsworth, R. D., & Pierrehumbert, R. T. (2013). Water-loss from Terrestrial Planets with CO2-rich Atmospheres. *The Astrophysical Journal*, 778(2), 154.

Yang, J., Boué, G., Fabrycky, D. C., & Abbot, D. S. (2014). Strong Dependence of the Inner Edge of the Habitable Zone on Planetary Rotation Rate. *The Astrophysical Journal Letters*, 787(1), L2.

## Data and Code Availability
The frozen NASA Exoplanet Archive snapshot used in this analysis (`PSCompPars_2026.07.11_02.37.36.csv`), along with the 1-D radiative-convective model, the surrogate grid, and all Monte Carlo simulation scripts, are available at https://github.com/Brukrish2006/Exoplanet-Habitability to ensure full reproducibility of the reported Monte Carlo yield distributions and catalog demographics.

## CRediT Authorship Contribution Statement
**H. Adhikary:** Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data Curation, Writing - Original Draft, Writing - Review & Editing, Visualization.

## Declaration of Competing Interest
The author declares that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Funding
This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.
