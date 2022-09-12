# ------------------------------------------------------ # 
# Script to create an input control simulation so the    #
# same values can be applied to perturbation runs and    #
# other sensitivity tests                                #
# ---                                                    #
# T. Kukla (Stanford Univ. 2021)                         #
#                                                        #
# Date created: 4 Dec 2021                               #
# Last modified: 4 Dec 2021                              #
# ------------------------------------------------------ # 
#### 
## This is an example script to show how one could change 
## the control case input file using tests on the MEBM 
## output 
#### 

rm(list=ls())

library(deSolve) 
library(rstudioapi)

# set working directory to source file directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# ------------------------------------------------------------------------- #
# LOAD MEBM files ********************************************************* #
# ------------------------------------------------------------------------- #
base.path.split <- strsplit(getwd(), '/')[[1]]
base.path.idx <- which(base.path.split == 'CODE_DISTRIBUTE')
sourcefile.Path <- paste0(base.path.split[c(1:base.path.idx)], collapse='/')
sourcefile.Dir <- paste(sourcefile.Path, 'Initialization', sep='/')
# bring in the initialization scripts
source(paste(sourcefile.Dir, 'ParSave_TRAIN.R', sep='/'))   # list of parameters to save for reproducibility
source(paste(sourcefile.Dir, 'CH2O-CHO_p_func_bistable.r', sep='/'))    # source CH2O-CHO functions
source(paste(sourcefile.Dir, 'MEBM_solve-bistable.R', sep='/'))        # solver
source(paste(sourcefile.Dir, 'MEBM_constants.R', sep='/'))    # physical constants
source(paste(sourcefile.Dir, 'MEBM_ODEfun.R', sep='/'))       # MEBM ODE to get mse
source(paste(sourcefile.Dir, 'MEBM_hydrofun.R', sep='/'))     # solve hydrologic fluxes
# source(paste(sourcefile.Dir, 'MEBM_main_run_INSOL.R', sep='/'))   # MEBM main solver - insolation-capable
source(paste(sourcefile.Dir, 'MultiStable_NAVIGATOR_sensitivityVersion.R', sep='/')) # script to navigate temp boundary conditions and other useful fxns
source(paste(sourcefile.Dir, 'CH2O-CHO_LogFile-WriterFXN.R', sep='/'))  # log file writer
# ... Special MEBM version so we can change it without affecting anything else
source(paste(sourcefile.Dir, 'MEBM_main_run_SENSITIVITY.R', sep='/'))   # MEBM main solver - insolation-capable


## ------------------------------------------------------------------------------------------------------ #
## INITIAL PARLIST TO UPDATE ---------------------------------------------------------------------------- 
# LOAD CONTROL INPUT 
listPath <- paste(sourcefile.Path, 'Input_Parameters/input_par_files', sep='/')
listFile <- 'PAR-LIST_bistable.RDS'
parList <<- readRDS(paste(listPath, listFile, sep='/'))   # MUST BE called parList for the MEBM_TRAIN.sensitivity fxn to work!!

## ------------------------------------------------------------------------------------------------------ #
## LAND INPUT ------------------------------------------------------------------------------------------- 
listPath.geog <- paste(sourcefile.Path, 'Geography/geog_inputfiles', sep='/')
LandFracFileName <- "Cat-eye_geography_2deg.RDS"
myLandFrac <- readRDS(paste(listPath.geog, LandFracFileName, sep='/'))
thisPaleogeo <- myLandFrac

## ------------------------------------------------------------------------------------------------------ #
## TEST ALTERNATIVE PARS ---------------------------------------------------------------------------------- 
newpar.vec.1  <- c(15, 20, 25)
newpar.vec.2 <- 3.35  # constant
newpar.vec.3 <- 222.5 # constant
pCO2.check <- c(200, 400, 800, 1600, 3200)  # doubling
# inits and lat area wts
LBtemp <<- -8 ; RBtemp <<- 8
e_sq <- 0.00669437999014 ; Re <- 6370000
deg2rad <- function(deg) {(deg * pi) / (180)}
lat_area_wt <- (pi*Re*cos(deg2rad(thisPaleogeo$Lat)) / (180 * sqrt(1-e_sq * sin(deg2rad(thisPaleogeo$Lat))**2)))
idx <- 1
for(i in 1:length(newpar.vec.1)){
  this.par <- newpar.vec.1[i]
  parList$MEBMpars$A0.slope <- this.par
  parList$MEBMpars$B0 <- newpar.vec.2
  parList$MEBMpars$A0.intercept <- newpar.vec.3
  parList$TRAINpars$iceAlf <- 0.75
  parList$MEBMpars$modCO2 <- 280
  for(j in 1:length(pCO2.check)){
    this.CO2 <- pCO2.check[j]
    df <- MEBM_TRAIN.sensitivity(pCO2_in = this.CO2, LandFrac_tx = thisPaleogeo, get_pars=FALSE )
    # plot(df$temp_C, type='l)
    # lines(df$temp_C, type='l')
    # get global temperature
    lat_area_wt <- (pi*Re*cos(deg2rad(df$Lat)) / (180 * sqrt(1-e_sq * sin(deg2rad(df$Lat))**2)))
    global_temp <- weighted.mean(x = df$temp_C, w = lat_area_wt)
    # save result
    if(idx == 1){
      outdf <- as.data.table(cbind(this.par, this.CO2, global_temp))
    } else{
      tempdf <- as.data.table(cbind(this.par, this.CO2, global_temp))
      outdf <- as.data.table(rbind(outdf, tempdf))
    }
    idx <- idx + 1
  }
  print(idx)
}
colnames(outdf) <- c("par", "pCO2", "tempC")
# get mean doubling
thisparval <- 3
thisCO2.val <- seq(length(pCO2.check), 2, by=-1)
outdf[par == newpar.vec.1[thisparval]][thisCO2.val] - outdf[par == newpar.vec.1[thisparval]][thisCO2.val-1]
# plot results
ggplot(outdf) + 
  coord_cartesian(expand=c(0,0), ylim=c(10, 45)) +
  geom_line(aes(x=pCO2, y=tempC, group=par, color=par)) +
  geom_line(data=outdf1, aes(x=pCO2, y=tempC, group=par, color=par), linetype='dashed') +
  scale_color_viridis_b() +
  scale_x_log10() +
  theme_few()

outdf1 <- outdf


## CHANGE FINAL RESULT ------------------------------------------------------------------------------------------- 
# [1] A0.slope
parList$MEBMpars$A0.slope <- 25
# [2] B0
parList$MEBMpars$B0 <- 3.35
# [3] A0.intercept
parList$MEBMpars$A0.intercept <- 222.5
# [4] iceAlf
parList$TRAINpars$iceAlf <- 0.75


## CLIMATE INPUT ------------------------------------------------------------------------------------------- 
# pCO2 <- 600
## ------------------------------------------------------------------------------------------------------ #
## UPDATE PARLIST STEP-BY-STEP; OR CONTINUE TO GENERATE A NEW ONE --------------------------------------- 
## --- ATMOSPHERIC CO2
#      update in MEBM and TRAIN
# parList$MEBMpars$pCO2_in <- parList$TRAINpars$pCO2.i <- 600

# ********************************************** # 
# SAVE OUTPUT ---------------------------------- # 
# saveHere <- 'CONTROL-INPUT'
saveName <- 'PAR-LIST_bistable.RDS'
# saveRDS(parList, paste(saveHere, saveName, sep='/'))
# ********************************************** # 



