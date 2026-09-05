# Calculated TDS spectra

This directory contains the calculated spectra for the final five-state model and the ideal-sink case.

Each CSV file contains:

```text
T_K,J_exp,J_calc,J_sink,J_trap1,J_trap2,J_trap3,J_trap4,J_trap5
```

where `T_K` - temperature (K), `J_exp` - experimental desorption flux (m^-2 s^-1), `J_calc` - desorption flux calculated with the finite rate of desorption from the open surface (m^-2 s^-1), `J_sink` - desorption flux calculated with the infinite rate of desorption (ideal sink) from the open surface (m^-2 s^-1), `J_trap` - effective trap contribution to the desorption flux (m^-2 s^-1).