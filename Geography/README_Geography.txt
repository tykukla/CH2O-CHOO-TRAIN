# ---------------------------------------------------------- # 
# NOTES ON THE GEOGRAPHY GENERATION CODE FOR CH2O-CHOO TRAIN #
# ---                                                        #
# T Kukla (Colostate Univ. 2022)                             #
# ---------------------------------------------------------- #

## GENERAL: 
## ... This directory contains the code to generate geography input
## ... files for the CH2O-CHOO TRAIN model as well as the input files

## /geog_inputfiles
## ... Series of .RDS files containing geography inputs that can be 
## ... read by CH2O-CHOO TRAIN. Note that files with `LIST' appended
## ... are multiple geography inputs ordered as a `list' element in R

## /generate_geography
## ... /idealized_geography
## ... ... Scripts to build various idealized geographic configurations.
## ... ... each reads in the PaleoGeo_dataframes file to use as a template.
## ... /paleogeography
## ... ... Only need to run the ``PaleoGeo_dataframes_MAKE.R'' file. The
## ... ... `GeogConstructDF_Fxn.R' file is a function to translate 2d data
## ... ... into 1D zonal mean input files. Note that it is not designed to 
## ... ... work with a wide range of input formats, so be careful... 
## ... ... The `_MAKE' file includes code to build the SpatialPolygons_5MyrSteps.R 
## ... ... file but it is commented out. It hasn't been tested in a few years. 