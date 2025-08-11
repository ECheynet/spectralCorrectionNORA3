# Spectral Correction for Extreme Value Analysis with NORA3  

This MATLAB Live Script corrects the underestimation of NORA3 extreme wind speed estimates using a spectral-slope adjustment method, following Bastine et al. (2018). The method is demonstrated with ~50 years of wind speed data from the Slåtterøy lighthouse (Norway).  

## Summary  

The Live Script `Documentation.mlx` demonstrates how to:  

1. Load mast and NORA3 wind speed data from NetCDF files.  
2. Compute and compare power spectral densities (PSDs) for measurements and model data.  
3. Identify the subrange of interest and adjust the PSD slope to match the expected −5/3 law.  
4. Apply the correction to the NORA3 time series.  
5. Perform extreme value analysis to estimate X-year return values before and after correction.  
6. Visualise time series and PSD comparisons.  

The method is based on the approach described in:  

Malekmohammadi, Shokoufeh, Etienne Cheynet, and Joachim Reuder. _Observation of Kelvin–Helmholtz billows in the marine atmospheric boundary layer by a ship-borne Doppler wind lidar_ *Scientific Reports* 15.1 (2025): 5245.

Bastine, D., Larsén, X., Witha, B., Dörenkämper, M., & Gottschall, J. (2018). Extreme winds in the new European wind atlas. *Journal of Physics: Conference Series*, 1102, 012006.  

---

## Content  

The repository contains:  

- **Documentation.mlx** – MATLAB Live Script explaining and demonstrating the workflow.  
- **NORA3_lighthouse_subset.nc** – Example NORA3 dataset (time, height, wind speed `U`).  
- **slatteroy_fyr_wind_speed.nc** – Mast wind speed dataset for Slåtterøy lighthouse.  
- **binAveraging.m** – Function to perform frequency bin averaging of spectral data.  
- **correctSpectrumSlope.m** – Function to adjust PSD slope within a specified frequency range.  
- **fitGEV.m** – Function to fit the Generalized Extreme Value (GEV) distribution and return extreme wind estimates.  
- **inpaint_nans.m** – Function to interpolate and fill missing values in a series.  
- **label.m** – Helper function to add text labels to MATLAB plots.  
- **LICENSE** – License file.  
- **README.md** – This file.  
