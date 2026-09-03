# Experimental TDS data

This directory contains the experimental deuterium TDS spectra used in the analysis.

Each experimental spectrum is stored as `<experiment_id>.csv`. The numeric data are read from the first four columns:
1. time, s;
2. temperature, K;
3. desorption-rate signal from the first measured channel, D/s;
4. desorption-rate signal from the second measured channel, D/s.

General information for all experiments, including the experiment identifier, heating rate, and boron-film thickness, is provided in ``experiments.csv``. The subset used for the simultaneous TDS fitting is listed in `experiments_5.csv`.