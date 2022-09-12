# ---------------------------------------------------------- # 
# NOTES ON THE INSOLATION CODE FOR CH2O-CHOO TRAIN           #
# ---                                                        #
# T Kukla (Colostate Univ. 2022)                             #
# ---------------------------------------------------------- #

## GENERAL:
## ... This directory includes example insolation inputs and an R script 
## ... to build your own insolation input

## EXAMPLE INSOL INPUTS:
## ... /insol_5Myr_ANN.RDS
## ... ... annual mean insolation for 5 million years of time
## ... /insol_10Myr_ANN.RDS
## ... ... annual mean insolation for 10 million years of time
## ... /insol_10Myr-2500steps_ANN.RDS
## ... ... annual mean insolation in 2500 year steps for 10 million years
## ... /insol_10Myr-2500steps_ANN-CONSTANT_T0.RDS
## ... ... same file structure as above, but the first timeslice insolation
## ... ... is repeated for 10 million years. This is because the first timeslice
## ... ... is not equivalent to the no-insol MEBM initial. So for direct comparison
## ... ... between insol off vs insol on we need to use the insol file. 

## CODE TO CREATE INSOL:
## ... /Generate_insol_inputs.R
## ... ... Change yrMin, yrMax, and yrStep if needed.
## ... ... If you want to use the output file for a run, set save.2.run = TRUE
## ... ... then the file can be accessed by config.R
