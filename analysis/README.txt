This folder has code the performs the analysis described in the paper and produces the figures.

1. estimate_CESM_models.R estimates the models for changes in mean temperature, temperature variability, and dew point quantiles from CESM-LE. It also outputs the supplement diagnostic figures for the gridcells corresponding to the best and worst model fit criteria described in the supplement. The estimated models are saved to be input into the simulations.

2. simulate_Temp_Dewp.R creates the observation-based simulations from GSOD data using the above estimated changes. We also make plots for the Minneapolis gridcell illustrating the simulation. Results of the simulation are saved. 

3. map_figures.R contains code for all of the figures that show maps analyzing output over CONUS.

4. helper_functions.R contain all of the functions needed to estimate the models and produce the simulations from steps 1 and 2.

5. plotting_functions.R contain the functions that make plots at individual gridcells as a part of steps 1 and 2.