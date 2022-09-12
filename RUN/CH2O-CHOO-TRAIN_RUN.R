# ------------------------------------------------- # 
#                   ALL ABOARD!!                    # 
#                                                   #
#      Script to run the CH2O-CHOO TRAIN model      #
#    Check the config.R file for correct settings   # 
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
library(rstudioapi)  # for setting working dir to file path

# clear environment
rm(list=ls())

# set working dir to file path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in the config file
CONFIG.FN <- 'config.R'      # name of the config file
parent.path.split <- strsplit(getwd(), '/')[[1]]
parent.path.idx <- which(parent.path.split == 'CODE_DISTRIBUTE')
parent.path <- paste0(parent.path.split[c(1:parent.path.idx)], collapse='/')
source(paste(parent.path, 'RUN', CONFIG.FN, sep='/'))


## [1] BUILD AND SAVE A FORCING TABLE (1 row per simulation) ---------------------------------------
source(code.init.fxn.loc)
FORCE.DF <- as.data.table(initialize.forcing())
pre_run.summary(FORCE.DF)     # quick check to know what the code will do
# save the forcing data.table
saveRDS(FORCE.DF, paste(saveDir, 'forcing_iters.RDS', sep='/'))


## [2] INITIALIZE THE MODEL CODE -------------------------------------------------------------------
initialize.mebm()      # moist energy balance model functions
initialize.c.cycle()   # carbon cycle and weathering functions

## [3] RUN THE CH2O-CHOO TRAIN ---------------------------------------------------------------------
# results will automatically save in the defined or auto-generated output folder
train.tracks()

