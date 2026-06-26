# The Role of Wave-Induced Stress and Drag Coefficients in Offshore Wind Power Production

**Manuscript:** 2025JC023741RR — submitted to *Journal of Geophysical Research: Oceans*

**Author:** Nícolas de Assis Bose

---

## Overview

This repository contains the analysis code and output figures for a study investigating the impact of wind–wave interactions on offshore wind energy production along the Brazilian Equatorial Margin (BEM). We compare a coupled atmosphere–wave model (WRF-SWAN, within the COAWST v3.8 framework) with a stand-alone atmospheric model (WRF) over the northern coast of Rio Grande do Norte, Brazil, during January–March (JFM) 2024.

Key findings:
- Under low wind speeds and strong swell, wave-induced stress (τ_wave) acts upward, transferring momentum from ocean to atmosphere and locally enhancing wind speeds.
- Under moderate-to-strong wind conditions (6 < Wspd₁₀ < 20 m s⁻¹), wave growth increases sea surface roughness, raising C_D and reducing wind speeds.
- The coupled WRF-SWAN model shows a ~4.3% difference in power production compared to stand-alone WRF, underscoring the importance of wave-induced stress for offshore wind resource assessment.

---

## Repository structure

```
.
├── drag_plots_owf.py               # Drag coefficient (C_D) analysis and figures
├── wave_stress_plots.py            # Wave-induced stress (τ_wave) analysis and figures
├── wind_validation_energy_analysis.py  # Wind validation vs LiDAR + energy production analysis
├── model_configuration/            # COAWST/WRF/SWAN configuration files
│   ├── namelist.input              # WRF namelist
│   ├── coupling_porto_ilha_test.in # COAWST coupling configuration
│   ├── swan_meqeast_d01.swn        # SWAN configuration domain 1
│   ├── swan_meqeast_d02.swn        # SWAN configuration domain 2
│   ├── porto_ilha_test_wrf.h       # WRF header file
│   └── build_coawst.sh             # COAWST build script
├── drag_figure/                    # Output figures from drag_plots_owf.py
├── wave_stress_figure/             # Output figures from wave_stress_plots.py
├── energy_figures/                 # Output figures from wind_validation_energy_analysis.py
└── database/                       # Placeholder — all data files are on Zenodo (see below)
```

---

## Data

All model output and input data files are archived on Zenodo:

> **Zenodo DOI:** `[DOI will be added upon acceptance]`

The following files are available there:

| File | Description | Size |
|---|---|---|
| `wrfswan_100m.nc` | WRF-SWAN 100 m wind fields (JFM 2024) |
| `wrfstand_100m.nc` | Stand-alone WRF 100 m wind fields (JFM 2024) |
| `energy_10MW_wrfswan.nc` | WRF-SWAN 10 MW power output (JFM 2024) |
| `energy_10MW_wrfstand.nc` | Stand-alone WRF 10 MW power output (JFM 2024) |
| `drag_owf_JFM2024.nc` | Drag coefficients and wave-age at OWF grid |
| `spec_1d_wrfswan_owf.nc` | WRF-SWAN 1D wave spectra at OWF grid |
| `spec_1d_tau_with_tau_chen_owf.nc` | WRF-SWAN 1D spectra + Chen τ diagnostics at OWF grid |
| `wrf_wind_owf_JFM2024.nc` | WRF wind variables at OWF grid |

**Note:** The LiDAR observational data (Porto-Ilha, ZephIR 300) are proprietary and not publicly available. Please contact the corresponding author for access inquiries.

To use the data, download the files from Zenodo and place them in the `database/` directory.

---

## Dependencies

The scripts run in the **conda base environment**. Required packages:

```
numpy
xarray
matplotlib
pandas
scipy
seaborn
```

Install any missing packages with:

```bash
conda install numpy xarray matplotlib pandas scipy seaborn
```

---

## Running the scripts

All scripts are run from the repository root. Each writes outputs to its own folder.

### 1. Wind validation and energy analysis

Validates WRF-SWAN and stand-alone WRF 100 m wind speeds against LiDAR observations and computes offshore wind power production differences.

```bash
python wind_validation_energy_analysis.py \
  --figures-dir energy_figures \
  --output-dir energy_outputs
```

Outputs figures to `energy_figures/` and CSV summary tables to `energy_outputs/`.

### 2. Drag coefficient plots

Computes and plots the drag coefficient (C_D) as a function of 10 m wind speed, wave age, and wind–wave alignment angle.

```bash
python drag_plots_owf.py \
  --output-dir drag_figure
```

### 3. Wave-induced stress plots

Plots the wind–wave stress difference (τ_WRF-SWAN − τ_WRF) and wind–wave alignment at the OWF grid.

```bash
python wave_stress_plots.py \
  --output-dir wave_stress_figure
```

---

## Contact

Corresponding author: Nícolas A. Bose — boseoceano@gmail.com
