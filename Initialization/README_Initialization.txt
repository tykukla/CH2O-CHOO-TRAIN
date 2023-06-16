# ---------------------------------------------------------- # 
# NOTES ON THE INITIALIZATION DIRECTORY FOR CH2O-CHOO TRAIN  #
# ---                                                        #
# T Kukla (Colostate Univ. 2022)                             #
# ---------------------------------------------------------- #

## GENERAL:
## ... This is the directory where the config.R file looks to find any
## ... model inputs that are prescribed. Geographies, insolation,
## ... control parameters, etc. MUST have their relevant files 
## ... located appropriately here. 

## ... Most subdirectories has its own readme. Still, a quick overview of their contents is below: 
## ... (in some cases, the info below is sufficient and there is no ReadMe in the subdir). 


## C_cycle_dat:
## ... This directory includes text data files for C-cycle initialization and forcing files
## ... The seawater values through time file is the most important for CH2O-CHOO (others are 
## ... not always read in.

## Code:
## ... all of the model code is here.

## Control_Params: 
## ... RDS files that are designed to be read by CH2O-CHOO TRAIN and 
## ... including all of the parameters required for input (examples 
## ... include climate sensitivity, ice sheet albedo, diffusivity, etc).
## ... To make your own input file, see 
## ... ... CODE_DISTRIBUTE/Input_Parameters/build_input_par_files

## Geography: 
## ... RDS files that can be used as geography input files for CH2O-CHOO TRAIN.
## ... To make your own input file, see
## ... ... CODE_DISTRIBUTE/Geography/generate_geography

## Insolation:
## ... RDS files that can be used as insolation input files for CH2O-CHOO TRAIN.
## ... To make your own input file, see 
## ... ... CODE_DISTRIBUTE/Insolation
