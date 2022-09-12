# --------------------------------------------- # 
# CH2O-CHOO TRAIN CONFIGUATION FILE             # 
#                                               # 
# Assign initial model conditions, forcings,    #
# and input files for time-transient            #
# CH2O-CHOO TRAIN simulations                   #
# --------------------------------------------- # 

## ******************************************** # 
## NOTE: Lines bounded by asterisks (*) such as #
## these are those which many users will want   # 
## to change. They include things like:         #
## -- forcings, the time grid, geography file   #
## -- parameters to change, and more            #
## Other lines may be modified too, but may     #
## require a more advanced understandign of the # 
## code to be sure nothing breaks.              #
##                                              #
## For questions, please contact:               #
## -- Tyler Kukla (tykukla@colostate.edu)       #
## -- (check tylerkukla.com for updated email)  #
## ******************************************** #
library(deSolve) 
library(data.table)
library(palinsol) 

## [1] RUN NAME, DIRECTORY, AND NOTES --------------------------------------------------------------
## ******************************************** # -- START MODIFY LINES
RunName_BASE <<- 'ExampleRun'    # for log file and output naming [note: log file automatically adds date to file name]
thisRunDir <<- 'Example_Results'  # directory to store results (will be automatically created if doesn't exist already!)
myName <- 'Dr. Test Run'      # user's name for log file (in case multiple researchers working on same project)
myNotes <- 'This is a test run. This note will be stored in the log file. I like to use the note to jog my memory of why I am conducting this specific run'
## ******************************************** # -- END MODIFY LINES
# find parent directory 
config.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)  # defaults to this file's path
parent.path.split <- strsplit(config.dir, '/')[[1]]
parent.path.idx <- which(parent.path.split == 'CODE_DISTRIBUTE')
parent.path <- paste0(parent.path.split[c(1:parent.path.idx)], collapse='/')
# build the save directory if it doesn't exist
saveDir <- paste(parent.path, thisRunDir, sep='/')
if(!file.exists(saveDir)){dir.create(file.path(saveDir))}


## [2] SELECT PAR TO CHANGE FROM CTRL VALUE --------------------------------------------------------
## ******************************************** # -- START MODIFY LINES
# name of the variable you're testing (!!must be consistent with its name in parList, read in below!!)
SENSITIVITY.PAR <<- 'budy_om'     # default set to budyko omega (moisture recycling efficiency parameter)
# vector of parameter values to test (will run one simulation for each)
SENSITIVITY.VEC <<- c(2.0, 3.5)   # input parameters to test (code will test CTRL value even if not provided here)
# other variables commonly modified
glwx.vector <<- c(0)              # fraction of runoff that causes weathering beneath an ice sheet
pCO2.vector <<- c(500)            # [ppmv] initial partial pressure of atmospheric co2 
# scale fluxes to initial CO2?
scale.to.init.co2 <<- TRUE        # [TRUE ; FALSE] if TRUE, then wx is scaled to CO2 and temperature at the first CO2 value (weathering scaler is constant for all iterations)
## ******************************************** # -- END MODIFY LINES
# read in the rest of the control parameters
parList.path <- paste(parent.path, 'Initialization/Control_Params', sep='/')
parList.name <- 'PAR-LIST.RDS'
parList <<- readRDS(paste(parList.path, parList.name, sep='/'))
# make a control index vector
par.here <- ifelse(any(names(parList[[1]]) %in% SENSITIVITY.PAR), 1, 2) # Find if the par is in list element 1 or 2
if(parList[[par.here]][[SENSITIVITY.PAR]] %in% SENSITIVITY.VEC){ # if the control value is already in the sensitivity vector, make the control id vector
  ctrl.here <- which(SENSITIVITY.VEC == parList[[par.here]][[SENSITIVITY.PAR]]) # identify where ctrl is
  CONTROL.idx <- rep('N', length(SENSITIVITY.VEC))   # establish ctrl vec
  CONTROL.idx[ctrl.here] <- 'Y'    # over-write the ctrl value
} else{   # if the control value is not in the original vector, add it
  ctrl.value <- parList[[par.here]][[SENSITIVITY.PAR]]  # extract control value
  new.SENSITIVITY.VEC <- c(ctrl.value, SENSITIVITY.VEC)  # add it to the existing vector
  SENSITIVITY.VEC <- new.SENSITIVITY.VEC[order(new.SENSITIVITY.VEC)]  # re-order to ascending order 
  CONTROL.idx <- ifelse(SENSITIVITY.VEC==ctrl.value, 'Y', 'N')  # write the control index
}


## [3] GEOGRAPHY / SPATIAL GRID --------------------------------------------------------------------
## ******************************************** # -- START MODIFY LINES
# name of the geography file (make sure it's in Initialization/Geography)
# can be a single data.table or a list of data.tables
LandFracFileName <<- 'ModernGeo_dataframe_2deg.RDS'
## ******************************************** # -- END MODIFY LINES
geog.input <- paste(parent.path, 'Initialization/Geography', LandFracFileName, sep='/')
geog.list <- readRDS(geog.input)


## [4] TIME GRID AND CARBON EMISSION FORCING -------------------------------------------------------
## ******************************************** # -- START MODIFY LINES
# time grid
duration <- 1e6       # [years] duration of model run
dt <- 5000            # [years] timestep
t <- seq(0, duration, dt)  # [years] vector of times
# Which timesteps to save climate data
save.step <- 2       # save climate result every [X] timesteps
save.clim.here <<- t[seq(1, length(t), save.step)]
# assign custom insolation? 
custom_insol <- FALSE  # [FALSE; TRUE] false uses non-changing insolation, if True provide an insol file following example format
insol.fn <- 'insol_10Myr-2500steps_ANN.RDS'  # if above is false, this file is not read in
# perturbation
t.exc.start <- 0.25e6  # [year] start of perturbation in model time
t.exc.end <- 0.45e6 # [year] end of perturbation
m <- c(0)   # [moles of C added] carbon emissions integrated across forcing (no forcing if zero)
n <- (t.exc.end-t.exc.start)/dt   # number of timesteps
d13C.extra <<- -5   # [per mille] isotopic composition of volcanic perturbation
## ******************************************** # -- END MODIFY LINES
if(custom_insol == TRUE){
  insol.file <- paste(parent.path, 'Initialization/Insolation', insol.fn, sep='/')
  insol_df <<- readRDS(insol.file)
}


## !! [NOTE - No more lines that users are likely to modify below this point] !! 


## [5] MEBM TIMEOUT AND SOLVER CONFIG --------------------------------------------------------------
cpu.timeMax <<- timeMax <<- 10   # [seconds] max time to allow MEBM to seek a solution before quitting
elapsed.timeMax <<- Inf          # [seconds] elapsed time (can cut the cpu.timeMax short, so we leave at infinity)
Bistable.Run <<- "N"             # ['Y' or 'N']  (yes is bistable, no is not bistable; turns off and on the LB and RB temp options)
SnowballFinder <<- "on"          # ["on" or "off"]  (CAN BE ON W/O BISTABILITY)  # this is needed if the solution isn't found in 10 secs (otherwise it crashes)
ClimStateTester <<- "off"        # ["on" or "off"]  (CAN BE ON W/O BISTABILITY)
impose_loop1_BCs <<- "Y"         # ['Y' or 'N'] (yes means the temperature boundary counditions in the first loop are applied to the second (instead of just the climate state))
# set temperature boundary conditions by climate state
LBtemp.options <<- c("both poles" = -8, "one pole" = -8, "no ice" = 8)   # [degC] Left bound temp (south pole) starting conditions for three possible climate states
RBtemp.options <<- c("both poles" = 5, "one pole" = 5, "no ice" = 8)     # [degC] Right bound temp (north pole) starting conditions for three possible climate states
if(Bistable.Run == "N"){LBtemp <<- LBtemp.initial <<- LBtemp.options[['no ice']] ; RBtemp <<- RBtemp.initial <<-  RBtemp.options[['no ice']]}
# Rules for navigating climate state changes (only matters if Bistable.Run == 'Y')
# ... ctrl+f `this.nav.count.conditions` in `MultiStable_NAVIGATOR-...` for more info
CLIM.NAV.CONDITIONS <<- c("new" = 4, "old" = 0, "other" = 1e3) # ORIGINAL = c("new"=4, "old"=0, "other"=1e3)
# Clim Resolution (in number of decimal points for rounding)
CO2_res <<- 4               # [] number of decimal points for CO2 inputs
q_res <<- 4                 # [] number of decimal points for runoff estimate (for any point in space) when runoff is [m/yr]
temp_res <<- 4              # [] number of decimal points for temperature inputs


## [6] MODEL RUN FILES -----------------------------------------------------------------------------
# Functions for initializing code
code.init.fn <-  'Initialize_Loop_and_Code_Fxns.R'
code.init.fxn.loc <<- paste(parent.path, 'Initialization/Code', code.init.fn, sep='/')
# MEBM run files
sourcefile.Dir <<- paste(parent.path, 'Initialization/Code', sep='/')  # where MEBM files are located
MEBMmain.file <<- "MEBM_main_run_SENSITIVITY.R"                # MEBM main solver 
ParSave.file <<- "ParSave_TRAIN.R"                             # list of parameters to save for reproducibility
MEBMsolver.file <<- "MEBM_solve-bistable-glwx.R"               # solver
MEBMconst.file <<- "MEBM_constants.R"                          # physical constants
MEBMode.file <<- "MEBM_ODEfun.R"                               # MEBM ODE to get mse
MEBMhydro.file <<- "MEBM_hydrofun.R"                           # solve hydrologic fluxes
CH2OCHOwritelog.file <<- 'CH2O-CHO_LogFile-WriterFXN.R'   # log file writer
# files that depend on insolation
if(custom_insol == TRUE){
  CH2OCHOfxn.file <<- "CH2O-CHO_p_func_bistable-trackBC-INSOL.R"       # source CH2O-CHO functions
  MEBMmultistable_nav.file <<- "MultiStable_NAVIGATOR-v2_sensitivityVersion-INSOL.R"  # script to navigate temp boundary conditions and other useful fxns
} else{
  CH2OCHOfxn.file <<- "CH2O-CHO_p_func_bistable-trackBC.R"       # source CH2O-CHO functions
  MEBMmultistable_nav.file <<- "MultiStable_NAVIGATOR-v2_sensitivityVersion.R"  # script to navigate temp boundary conditions and other useful fxns
}
# C-cycle run files
sourcefile.Ccycle <<- paste(parent.path, 'Initialization/C_cycle_dat', sep='/')  # where c-cycle chem files are stored
