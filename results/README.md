# Results

This directory contains the numerical results used in the analysis of the TDS model.

## Files

- `calculated_spectra/` — final calculated spectra for the five-trap finite-rate model and the ideal-sink approximation. The files also contain the release contributions associated with the five traps.
- `model_order.csv` — fit-quality metrics for the three-, four-, five-, and six-trap finite-rate models and the five-trap ideal-sink case.
- `sensitivity.csv` — fitted parameters from the sensitivity analysis. 
- `compare_TDS.py` — generates the experimental/model TDS comparison from the stored spectra and can optionally recalculate the spectra with FESTIM.
- `plot_Kr.py` — plots the temperature-dependent surface recombination coefficient obtained for all sensitivity cases.
- `tds_comparison.png` — final comparison of the experimental and calculated TDS spectra.
- `Kr_dependencies.png` — surface recombination coefficients corresponding to all sensitivity fits.

Each CSV file for calculated spectra contains:

```text
T_K,J_exp,J_calc,J_sink,J_trap1,J_trap2,J_trap3,J_trap4,J_trap5
```

where `T_K` - temperature (K), `J_exp` - experimental desorption flux (m^-2 s^-1), `J_calc` - desorption flux calculated with the finite rate of desorption from the open surface (m^-2 s^-1), `J_sink` - desorption flux calculated with the infinite rate of desorption (ideal sink) from the open surface (m^-2 s^-1), `J_trap` - effective trap contribution to the desorption flux (m^-2 s^-1).
