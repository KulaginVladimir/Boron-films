# Results

This directory contains the numerical results used in the analysis of the TDS model.

## Files

- `calculated_spectra/` — final calculated spectra for the five-state finite-rate model and the ideal-sink case. The files also contain the release contributions associated with five effective trapping states.

- `model_order.csv` — fit-quality metrics for the three-, four-, five-, and six-state finite-rate models and the five-state ideal-sink case.

- `sensitivity.csv` — fitted parameters obtained from the sensitivity analysis in which the fixed trapping and detrapping prefactors `k0` and `p0` are varied independently by factors of 0.1, 1, and 10.

- `compare_TDS.py` — generates the experimental/model TDS comparison from the stored spectra and can optionally recalculate the spectra with FESTIM.

- `plot_surface_pivot.py` — analyses the surface recombination rates obtained from `sensitivity.csv`.

- `tds_comparison.png` — final comparison of the experimental and calculated TDS spectra.

- `surface_pivot.png` — pivot plot of fitted recombination rates.

Each CSV file for calculated spectra contains:

```text
T_K,J_exp,J_calc,J_sink,J_trap1,J_trap2,J_trap3,J_trap4,J_trap5
```

where `T_K` - temperature (K), `J_exp` - experimental desorption flux (m^-2 s^-1), `J_calc` - desorption flux calculated with the finite rate of desorption from the open surface (m^-2 s^-1), `J_sink` - desorption flux calculated with the infinite rate of desorption (ideal sink) from the open surface (m^-2 s^-1), `J_trap` - effective trap contribution to the desorption flux (m^-2 s^-1).
