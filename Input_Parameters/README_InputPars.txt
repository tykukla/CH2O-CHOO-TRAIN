# ---------------------------------------------------------- # 
# NOTES ON THE INPUT PARAMETERS SUBDIR FOR CH2O-CHOO TRAIN   #
# ---                                                        #
# T Kukla (Colostate Univ. 2022)                             #
# ---------------------------------------------------------- #

## GENERAL:
## ... This subdirectory contains a folder to build files of the 
## ... parList.RDS format (to be read into CH2O-CHOO TRAIN) and 
## ... a folder to save some example files. NOTE: The example files
## ... saved here WILL NOT be readable to the model! They MUST also
## ... be saved in CODE_DISTRIBUTE/Initialization/Control_Params
## ... to be called from the config file. 

## ... /build_input_par_files
## ... ... Contains a script that serves as a template for building
## ... ... your own input parameter files of the parList.RDS format

## ... /input_par_files
## ... ... Contains example parList.RDS files that can be read in 
## ... ... to the CH2O-CHOO TRAIN (if, of course, they're also saved 
## ... ... in the Initialization/Control_Params directory). Files 
## ... ... DO NOT need to be saved here to be read by the config file.
## ... ... Instead, this directory serves as workspace for testing 
## ... ... different versions before reaching a final version. 