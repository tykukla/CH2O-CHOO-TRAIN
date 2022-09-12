# ---------------------------------------------------------- # 
# NOTES ON THE DIRECTORIES FOR CH2O-CHOO TRAIN MODEL         #
# ---                                                        #
# T Kukla (Colostate Univ. 2022)                             #
# ---------------------------------------------------------- #

## GENERAL:
## ... This directory contains subdirectories for the model code, scripts to run the model,
## ... scripts to build the model inputs, generate figures, and more. Each sub directory
## ... has its own README file to walk through its contents. 

## ... /Example_Results
## ... ... Includes example output from the run produced by the defulat 'config.R' file
## ... ... as well as figures produced by the PROCESS/PLOT_RESULTS.R file.
## ... ... NOTE: Results directories are automatically created whenever the model is run,
## ... ... users don't need to create one in advance (but they can). 

## ... /Geography
## ... ... Scripts to create new zonal mean geographies that the model can read, as well 
## ... ... as some idealized and paleo examples. 

## ... /Initialization
## ... ... This directory has the model code, initial input files, and more. 
## ... ... Any model inputs must be saved in the appropriate Initialization subdir.
## ... ... Otherwise, the user shouldn't need to do much with this folder.

## ... /Input_Parameters
## ... ... Code for users to create model input files for climate and carbon cycling 
## ... ... that can be read by the CH2O CHOO TRAIN. 

## ... /Insolation
## ... ... Code to create, and examples of, insolation input files.

## ... /PROCESS
## ... ... Scripts to plot model results.

## ... /RUN
## ... ... Scripts to run the model. 