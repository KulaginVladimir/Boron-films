# TDS simulation of D from B-D co-deposited films

## Overview

This repository contains experimental data and FESTIM 2.1 scripts used to model deuterium thermal desorption spectra from B-D co-deposited films.

Five TDS spectra measured for different film thicknesses and heating rates are fitted simultaneously using a common set of effective traps. The initial deuterium inventory for each experiment is determined independently from the corresponding experimental TDS integral.

## Repository structure

```text
experimental_data/
    Raw experimental TDS spectra and experiment metadata.

fit_TDS/
    tds_model.py                  FESTIM model
    fit.py                        Simultaneous fitting of the TDS spectra
    run_fit_parallel.sh           Slurm script for parallel fitting
    sensitivity.py                Sensitivity analysis
    run_sensitivity_parallel.sh   Slurm array script for sensitivity analysis

results/
    calculated_spectra/           Final calculated TDS spectra
    model_order.csv               Model-order and ideal-sink fit metrics
    sensitivity.csv               Results of the sensitivity analysis
    compare_TDS.py                Experimental/model TDS comparison
    plot_Kr.py                    Surface recombination coefficients from sensitivity fits
    tds_comparison.png            Final TDS comparison figure
    Kr_dependencies.png           Surface recombination coefficients for sensitivity fits
    README.md                     Description of the stored results

environment.yml
    Conda environment used for the calculations.
```

## Installation

Clone the repository:
```bash
git clone https://github.com/KulaginVladimir/Boron-films.git
cd Boron-films
```

Create and activate the Conda environment:
```bash
conda env create -f environment.yml
conda activate boron-films-env
```

## TDS fitting

The thermally activated surface model can be fitted using three to six effective trapping sites.

For a local run:
```bash
python fit_TDS/fit.py 3
python fit_TDS/fit.py 4
python fit_TDS/fit.py 5
python fit_TDS/fit.py 6
```

For a parallel run on a cluster with Slurm:
```bash
sbatch fit_TDS/run_fit_parallel.sh 3
sbatch fit_TDS/run_fit_parallel.sh 4
sbatch fit_TDS/run_fit_parallel.sh 5
sbatch fit_TDS/run_fit_parallel.sh 6
```

The ideal-sink limit is fitted for the selected five-trap model:

```bash
python fit_TDS/fit.py 5 sink
```

or similar on the cluster:

```bash
sbatch fit_TDS/run_fit_parallel.sh 5 sink
```

## Sensitivity tests

The robustness of the five-trap solution with respect to the fixed trapping and detrapping kinetic prefactors can be tested by repeating the full fit for different values of `p0`, `D0`, `Ed`.

To perform the sensitivity analysis, submit the Slurm array:
```bash
sbatch fit_TDS/run_sensitivity_parallel.sh
```

The final fitted parameters from this analysis are stored in `results/sensitivity.csv`.

## Results

To generate the TDS comparison from the spectra stored in `results/calculated_spectra/`:

```bash
python results/compare_TDS.py
```

To recalculate the final spectra with FESTIM and overwrite the stored calculated spectra:

```bash
python results/compare_TDS.py --recalculate
```

The surface-rate Arrhenius family derived from the sensitivity results can be analysed with:

```bash
python results/plot_surface_pivot.py
```

The numerical results used for the model-order, surface-boundary, and sensitivity analyses are documented in `results/README.md`.
