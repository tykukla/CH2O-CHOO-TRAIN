# ---------------------------------------------------------- # 
# NOTES ON THE RUN SCRIPTS FOR CH2O-CHOO TRAIN               #
# ---                                                        #
# T Kukla (Colostate Univ. 2022)                             #
# ---------------------------------------------------------- #

## GENERAL: 
## ... This directory includes the main two files the user will interact
## ... with to run the model. 

## WORKFLOW: 
## [1] config.R -- first, use this file to set the experiments for CH2O-CHOO TRAIN to run.
## [2] CH2O-CHOO-TRAIN_RUN.R -- second, run this file to execute the simulations and save results. 

## ... /config.R 
## ... ... The primary lines the user may want to modify are denoted by:
## ... ... ... '## ******************************************** # -- START MODIFY LINES'
## ... ... ... '## ******************************************** # -- END MODIFY LINES'
## ... ... Use this file to define the input files, set the perturbations, etc. 
## ... ... This file is sourced by CH2O-CHOO-TRAIN_RUN.R

## ... /CH2O-CHOO-TRAIN_RUN.R
## ... ... User can run this file fully to run the model. 
## ... ... The script first generates the forcing data.table, then saves that 
## ... ... data.table, then prints a pre-run summary for the user so the user 
## ... ... is not surprised by how many simulations are run. 
## ... ... The train.tracks() function runs the model. 


## !! NOTE ON CONFIG !! 
## ... We ran into file corruption errors when trying to save files of the same
## ... name in the same folder as a previous run (even if we deleted the folder).
## ... Didn't get to the bottom of this, but to avoid this issue we suggest 
## ... setting a new run name for each simulation. 