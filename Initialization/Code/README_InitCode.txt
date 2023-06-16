# ---------------------------------------------------------- # 
# NOTES ON THE INITIALIZATION CODE FOR CH2O-CHOO TRAIN       #
# ---                                                        #
# T Kukla (Colostate Univ. 2022)                             #
# ---------------------------------------------------------- #

## GENERAL:
## ... This subdirectory contains all of the code needed to run the model
## ... We go through a brief overview of each file's contents below.

## ... /CH2O-CHOO_LogFile-WriterFXN.R
## ... ... Function to write the simulation log file

## ... /CH2O-CHOO_p_func_bistable-trackBC.R ; /CH2O-CHOO_p_func_bistable-trackBC-INSOL.R
## ... ... The functions for running the CH2O-CHOO TRAIN without insolation (first 
## ... ... listed above) or with insolation (...-INSOL.R). This is the primary set of 
## ... ... functions, including weathering functions, the carbon cycle solver, and 
## ... ... the model initialization function.

## ... /Initialize_Loop_and_Code_Fxns.R
## ... ... This script reads in data from the config.R file to build a data.table 
## ... ... that serves as a road-map for the CH2O-CHOO TRAIN simulations. When the 
## ... ... user opts to run more than one simulation at a time, a function in this 
## ... ... script defines each individual run. Other functions simply source the code
## ... ... for the individual model components. 

## ... /MEBM_constants.R
## ... ... Physical constants used in the MEBM

## ... /MEBM_hydrofun.R
## ... ... Function that reads in the results of the ODE and solves the zonal mean hydroclimate

## ... /MEBM_main_run_SENSITIVITY.R
## ... ... Function that reads in parList and uses it to initialize basic MEBM parameters

## ... /MEBM_ODEfun.R 
## ... ... The boundary value problem set of equations and solver

## ... /MEBM_solve-bistable-glwx.R
## ... ... Script to call the ODE and hydrofun to solve and return the full MEBM output

## ... /MultiStable_NAVIGATOR-v2_sensitivityVersion.R ; /MultiStable_NAVIGATOR-v2_sensitivityVersion-INSOL.R
## ... ... Functions to navigate to a stable climate state in case the model crashes or returns a NULL result.
## ... ... The first file applies to the constant insolation case, and second file to variable insolation.

## ... /ParSave_TRAIN.R
## ... ... A helper function which defines variable names that can be saved in a parList.RDS
## ... ... format for later application.


