import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import functools
from scipy.interpolate import RegularGridInterpolator

# Ensure tier1_model is in path
sys.path.append(os.path.abspath("tier1_model"))
import validate_model
import rc_model

# Setup paths (relative/local)
artifact_dir = "artifacts"
if not os.path.exists(artifact_dir):
    os.makedirs(artifact_dir)
csv_path = "PSCompPars_2026.07.11_02.37.36.csv"

print("Loading data...")
df = pd.read_csv(csv_path, comment='#', low_memory=False)

# Clean data (require T_eff and Insolation)
df_clean = df.dropna(subset=['st_teff', 'pl_insol']).copy()

# Filter out non-main-sequence stars
df_clean = df_clean[(df_clean['st_teff'] >= 2500.0) & (df_clean['st_teff'] <= 7200.0)]

# Setup Kopparapu Outer Edge
def kopparapu_outer_seff(T_eff):
    T_star = T_eff - 5780.0
    Seff_sun = 0.356
    a = 6.171e-5
    b = 1.698e-9
    c = -3.198e-12
    d = -5.575e-16
    return Seff_sun + a*T_star + b*T_star**2 + c*T_star**3 + d*T_star**4

print("Calculating HZ boundaries...")
# Run standard baseline calibration
kappa_star = rc_model.calibrate_kappa(target_olr=282.0)
_, _, olr_limit_base, _ = rc_model.find_olr_limit(kappa_star)
Seff_sun_target = validate_model.kopparapu_seff(5780.0)
target_albedo_sun = 1.0 - (4.0 * olr_limit_base) / (validate_model.S0 * Seff_sun_target)

from scipy.optimize import root_scalar
def albedo_resid(l0):
    return validate_model.calc_albedo(5780.0, l0, n=4.0) - target_albedo_sun
res = root_scalar(albedo_resid, bracket=[1e-8, 1e-6])
lambda_0 = res.root

@functools.lru_cache(maxsize=None)
def custom_model_inner_seff(T_eff):
    A = validate_model.calc_albedo(T_eff, lambda_0, n=4.0)
    return (4.0 * olr_limit_base) / (validate_model.S0 * (1.0 - A))

# Calculate baseline boundaries
df_clean['Seff_Kop_Inner'] = df_clean['st_teff'].apply(validate_model.kopparapu_seff)
df_clean['Seff_Model_Inner'] = df_clean['st_teff'].apply(custom_model_inner_seff)
df_clean['Seff_Outer'] = df_clean['st_teff'].apply(kopparapu_outer_seff)

df_clean['In_HZ_Kop'] = (df_clean['pl_insol'] <= df_clean['Seff_Kop_Inner']) & (df_clean['pl_insol'] >= df_clean['Seff_Outer'])
df_clean['In_HZ_Model'] = (df_clean['pl_insol'] <= df_clean['Seff_Model_Inner']) & (df_clean['pl_insol'] >= df_clean['Seff_Outer'])

def get_spt(teff):
    if teff < 3700: return 'M'
    elif teff < 5200: return 'K'
    elif teff < 6000: return 'G'
    else: return 'F'
df_clean['SpT'] = df_clean['st_teff'].apply(get_spt)

print("Generating Figure 1...")
plt.figure(figsize=(8, 6))
hz_counts = pd.DataFrame({
    'Kopparapu': df_clean[df_clean['In_HZ_Kop']]['SpT'].value_counts(),
    'Custom Model': df_clean[df_clean['In_HZ_Model']]['SpT'].value_counts()
}).fillna(0).loc[['M', 'K', 'G', 'F']]
hz_counts.plot(kind='bar', color=['#1f77b4', '#ff7f0e'], edgecolor='black')
plt.title('HZ Catalog Membership Sensitivity to Climate Model Choice')
plt.xlabel('Host Star Spectral Type')
plt.ylabel('Number of HZ Candidates')
plt.xticks(rotation=0)
plt.legend(title='Inner Edge Model')
plt.tight_layout()
plt.savefig(os.path.join(artifact_dir, 'fig1_hz_sensitivity.png'), dpi=300)
plt.close()

print("Generating Figure 2...")
plt.figure(figsize=(10, 7))
colors = {'M': 'red', 'K': 'orange', 'G': 'green', 'F': 'blue'}
for spt, color in colors.items():
    subset = df_clean[df_clean['SpT'] == spt]
    plt.scatter(subset['st_teff'], subset['pl_insol'], color=color, label=spt, alpha=0.5, s=15, edgecolors='none')

plt.yscale('log')
plt.gca().invert_xaxis()
t_range = np.linspace(2500, 7000, 100)
kop_in = [validate_model.kopparapu_seff(t) for t in t_range]
mod_in = [custom_model_inner_seff(t) for t in t_range]
kop_out = [kopparapu_outer_seff(t) for t in t_range]
plt.plot(t_range, kop_in, 'r--', label='Kopparapu Inner Edge', linewidth=2)
plt.plot(t_range, mod_in, 'b-', label='Custom Model Inner Edge', linewidth=2)
plt.plot(t_range, kop_out, 'k:', label='Outer Edge (Max Greenhouse)', linewidth=2)
plt.ylim(0.1, 10)
plt.title('Insolation vs Stellar Temperature')
plt.xlabel('Stellar Effective Temperature (K)')
plt.ylabel('Stellar Insolation ($S_{eff}$)')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(artifact_dir, 'fig2_seff_teff.png'), dpi=300)
plt.close()

print("Running true Monte Carlo Perturbed-Physics Ensemble (Layer B)...")
# Load precomputed OLR limit grid
npz = np.load('rc_model_grid.npz')
olr_interp = RegularGridInterpolator((npz['kappas'], npz['t_strats']), npz['olr_grid'])

# Precompute albedo grid for fast interpolation
print("Precomputing albedo grid for Monte Carlo...")
teffs_grid = np.linspace(2500, 7200, 25)
ns_grid = np.linspace(3.5, 4.5, 10)
albedo_grid = np.zeros((len(teffs_grid), len(ns_grid)))
for i, teff in enumerate(teffs_grid):
    for j, n_val in enumerate(ns_grid):
        albedo_grid[i,j] = validate_model.calc_albedo(teff, lambda_0, n=n_val)
albedo_interp = RegularGridInterpolator((teffs_grid, ns_grid), albedo_grid)

# Monte Carlo Loop
np.random.seed(42)
N_DRAWS = 5000
kappa_draws = np.random.normal(loc=kappa_star, scale=0.003, size=N_DRAWS)
kappa_draws = np.clip(kappa_draws, npz['kappas'].min(), npz['kappas'].max()) # Keep in grid bounds
tstrat_draws = np.random.uniform(low=180.0, high=220.0, size=N_DRAWS)
n_draws = np.random.uniform(low=3.5, high=4.5, size=N_DRAWS)

candidate_counts = []

print(f"Sampling {N_DRAWS} boundaries and cross-referencing catalog...")
# Vectorize catalog for fast evaluation
pl_insol = df_clean['pl_insol'].values
st_teff = df_clean['st_teff'].values
seff_out = df_clean['Seff_Outer'].values

for k, t_s, n_val in zip(kappa_draws, tstrat_draws, n_draws):
    # Get OLR limit
    olr_l = olr_interp([k, t_s])[0]
    # Get albedo for all stars
    pts = np.column_stack((st_teff, np.full_like(st_teff, n_val)))
    A_stars = albedo_interp(pts)
    # Compute Seff inner edge for all stars
    seff_in = (4.0 * olr_l) / (validate_model.S0 * (1.0 - A_stars))
    
    # Count how many planets fall in HZ
    in_hz = (pl_insol <= seff_in) & (pl_insol >= seff_out)
    candidate_counts.append(np.sum(in_hz))

print("Generating Figure 3...")
plt.figure(figsize=(8, 6))
# Convert candidate counts to hypothetical eta_earth scaling (approximate scaling)
# Just use raw counts for plotting to be perfectly transparent
plt.hist(candidate_counts, bins=np.arange(min(candidate_counts)-0.5, max(candidate_counts)+1.5, 1), density=True, color='purple', edgecolor='black', alpha=0.7)
plt.title('Monte Carlo Yield Distribution\n(Perturbed-Physics Ensemble)')
plt.xlabel('Number of Candidates in HZ')
plt.ylabel('Density')
plt.axvline(np.mean(candidate_counts), color='k', linestyle='dashed', linewidth=2, label=f'Mean = {np.mean(candidate_counts):.1f}')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(artifact_dir, 'fig3_posteriors.png'), dpi=300)
plt.close()

table1 = df_clean[(df_clean['In_HZ_Kop'] | df_clean['In_HZ_Model'])][['pl_name', 'st_teff', 'SpT', 'pl_insol', 'In_HZ_Kop', 'In_HZ_Model']]
table1.to_csv(os.path.join(artifact_dir, 'hz_candidates_comparison.csv'), index=False)

print("All figures and tables generated successfully.")
