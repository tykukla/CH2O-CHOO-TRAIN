# ------------------------------------------------- # 
#                   ALL ABOARD!!                    # 
#                                                   #
#     Script to calculate CH2O-CHOO TRAIN climate   #
#          for timesteps that weren't saved         #  
#                                                   #
#                                                   #
# *************************************             #
# T Kukla (Colostate Univ. 2022)                    #
# KV Lau (Penn State Univ. 2022)                    #
# DE Ibarra (Brown Univ. 2022)                      #
# JKC Rugenstein (Colostate Univ. 2022)             #
# *************************************             # 
#                                                   #
# ------------------------------------------------- # 

## NOTE: This script is designed to calculate climate 
## for all timesteps for some model output. This is 
## a slightly computationally faster way to calculate 
## missing steps, but we suggest just changing the 
## config file to save the correct number of steps


library(data.table)  # preferred data format
library(rstudioapi)  # for setting working dir to file path
library(stringr)     # for finding correct filenames

# clear environment
rm(list=ls())

# set working dir to file path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# read in the functions to help processing
source('processing_functions.R')

# find results
RESULTS.FOLDER <- 'MyTestRun_Results'
dat <- collect.dat(results.path = RESULTS.FOLDER)

# Initialize MEBM (self-consistent with simulation)
prep.MEBM(results.path = RESULTS.FOLDER)

## CALCULATE MISSING CLIM
out.clim <- CLIM.RUN()

