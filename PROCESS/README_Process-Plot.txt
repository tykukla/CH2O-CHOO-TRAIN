# ---------------------------------------------------------- # 
# NOTES ON THE PROCESS/PLOT CODE FOR CH2O-CHOO TRAIN         #
# ---                                                        #
# T Kukla (Colostate Univ. 2022)                             #
# ---------------------------------------------------------- #

## GENERAL: 
## ... This directory has scripts for processing and plotting model results.
## ... The user should ONLY need to modify PLOT_RESULTS.R to generate plots.

## ... /PLOT_RESULTS.R
## ... ... This is the main script for the user. Running it will display plots
## ... ... and, if selected, save plots to a .../figures subdirectory in the 
## ... ... main model output directory.

## ... /parameters_plotting.R
## ... ... Code sourced by PLOT_RESULTS.R which sets some plotting dimensions/decisions.

## ... /processing_functions.R
## ... ... Functions for reading in and bringing data together, and then saving output plots. 

## ... /Calculate_Clim_Timeseries.R 
## ... ... Code for recovering climate timeseries data if it was missing from the original run.