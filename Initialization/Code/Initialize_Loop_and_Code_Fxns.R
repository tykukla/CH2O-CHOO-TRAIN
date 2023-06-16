# --------------------------------------------------- #
#            INITIALIZE BISTABILITY FORCING           #
#           AND MEBM + CARBON CYCLE SCRIPTS           #
# --------------------------------------------------- #


## [1] INITIALIZE FORCING
## Function to initialize a list of parameters
## Each element is a vector of length = model.iters
## where the length depends on the number of unique 
## model cases (calculated by this fxn)
## 
## For example, to run two CO2 levels at two ice sheet wx coefficients
## we would get four total iterations. Other terms that don't vary 
## would just have their value repeated four times in the list
initialize.forcing <- function(co2.force = pCO2.vector, # atmos CO2 vector
                               glwx.force = glwx.vector, # glacier weathering vector
                               m.force = m,   # vector of mass of carbon injected (or removed)
                               t.start.force = t.exc.start, # time of volcanic anomaly start
                               t.end.force = t.exc.end,   # time of volcanic anomaly end
                               SENSITIVITY.force = SENSITIVITY.VEC, # vector of sensitivity values
                               duration.force = duration,   # [years] how long model is run
                               dt.force = dt,    # [years] timestep
                               CTRL.idx.force = CONTROL.idx,  # which sensitivity.vec values are equal to control case value
                               geog.force = geog.list   # geography input files
){
  # check whether geography is a list of many geogs, or a single geog
  if(is.null(nrow(geog.force))){  # if it doesn't have rows, then the vec is the list
    geog.vec <- names(geog.force)
  } else{   # if it does have rows, then the vec is just this file
    geog.vec <- LandFracFileName 
  }
  
  # Get iters associated with time (three start and three end times are three simulations, not nine)
  length.time.force <- max(c(length(t.start.force), length(t.end.force), length(duration.force), length(dt.force)))
  # make a fake vec to expand the grid 
  fake.vec <- rep(-9999, length.time.force)
  
  # expand grid (get every combination of relevant inputs)
  grid.expansion <- expand.grid(co2.force, glwx.force, m.force, SENSITIVITY.force, geog.vec, fake.vec)
  colnames(grid.expansion) <- c("co2.force", "glwx.force", "m.force", "SENSITIVITY.force", "geography.force", "x")
  # we put fake vec last because we know it will be repeated in row (rep(THIS, each=nrow(grid.expansion)/length(fake.vec)))
  # ... remove fake vec (was a placeholder)
  grid.expansion <- subset(grid.expansion, select=-x)
  # ... add the time force vecs
  grid.expansion$t.start.force <- rep(t.start.force, each = (nrow(grid.expansion)/length(t.start.force)) )
  grid.expansion$t.end.force <- rep(t.end.force, each = (nrow(grid.expansion)/length(t.end.force)) )
  grid.expansion$duration.force <- rep(duration.force, each = (nrow(grid.expansion)/length(duration.force)) )
  grid.expansion$dt.force <- rep(dt.force, each = (nrow(grid.expansion)/length(dt.force)) )
  
  # ... add control index
  if(length(CTRL.idx.force) < 2){
    grid.expansion$CTRL.idx <- CTRL.idx.force[1]
  } else{
    # equation relating index to value in sensitivity.vec
    sens.idx.fun <- approxfun(x=SENSITIVITY.force, y=c(1:length(SENSITIVITY.force)))
    sens.idx <- sens.idx.fun(grid.expansion$SENSITIVITY.force)
    grid.expansion$CTRL.idx <- CTRL.idx.force[sens.idx]
  }
  
  # order by control index so control iter is first
  grid.expansion <- as.data.table(grid.expansion)  # for the 'order' fxn to work
  grid.out <- grid.expansion[order(-CTRL.idx)]
  
  # return the result
  return(grid.out)
}




# -------------------------------------------------------------------- # 
# -------------------------------------------------------------------- # 
# -------------------------------------------------------------------- # 



## [2] INITIALIZE MEBM FILES 
##     terms are sourced (fxn does not need to be run to a variable)
initialize.mebm <- function(sourcefile.Dir.in = sourcefile.Dir,  # where are the files located 
                            MEBMmain.file.in = MEBMmain.file,  # MEBM main solver 
                            ParSave.file.in = ParSave.file,   # list of parameters to save for reproducibility
                            CH2OCHOfxn.file.in = CH2OCHOfxn.file,  # source CH2O-CHO functions
                            MEBMsolver.file.in = MEBMsolver.file,  # solver
                            MEBMconst.file.in = MEBMconst.file,  # physical constants
                            MEBMode.file.in = MEBMode.file,  # MEBM ODE to get mse
                            MEBMhydro.file.in = MEBMhydro.file,  # solve hydrologic fluxes
                            MEBMmultistable_nav.file.in = MEBMmultistable_nav.file,  # script to navigate temp boundary conditions and other useful fxns
                            CH2OCHOwritelog.file.in = CH2OCHOwritelog.file   # log file writer
){  
  # bring in the initialization scripts
  source(paste(sourcefile.Dir.in, MEBMmain.file.in, sep='/')) 
  source(paste(sourcefile.Dir.in, ParSave.file.in, sep='/'))   
  source(paste(sourcefile.Dir.in, CH2OCHOfxn.file.in, sep='/')) 
  source(paste(sourcefile.Dir.in, MEBMsolver.file.in, sep='/'))       
  source(paste(sourcefile.Dir.in, MEBMconst.file.in, sep='/'))    
  source(paste(sourcefile.Dir.in, MEBMode.file.in, sep='/'))    
  source(paste(sourcefile.Dir.in, MEBMhydro.file.in, sep='/'))     
  source(paste(sourcefile.Dir.in, MEBMmultistable_nav.file.in, sep='/')) 
  source(paste(sourcefile.Dir.in, CH2OCHOwritelog.file.in, sep='/'))  
}



# -------------------------------------------------------------------- # 
# -------------------------------------------------------------------- # 
# -------------------------------------------------------------------- # 



## [3] INITIALIZE CARBON CYCLE 
##     terms are sourced or just read right to the environment (fxn does not need to be run to a variable)
initialize.c.cycle <- function(sourcefile.Dir.in = sourcefile.Ccycle){
  # ------------------------------------------------------------------------- #
  # LOAD GEOCHEM INITIALIZATION FILES *************************************** #
  # ------------------------------------------------------------------------- #
  # CLiBeSO Files
  sw.file <<- "seawater values over time_interp.txt"   # these go to global env. because they're read into the log file
  Berner.volc.file <<- "Berner Volcanism Interp.txt"
  # sens.file <<- "PaleozoicSensValues.txt"
  CO2.file <<- "InterpCO2_Aug2018.txt"
  # read in
  sw.chem <<- read.table(file = paste(sourcefile.Dir.in, sw.file, sep='/'), header = TRUE)
  Berner.volc <<- read.table(file = paste(sourcefile.Dir.in, Berner.volc.file, sep='/'), header = TRUE)
  # sens.parm <<- read.table(file = paste(sourcefile.Dir, sens.file, sep='/'), header = TRUE)
  interpCO2 <<- read.table(file = paste(sourcefile.Dir.in, CO2.file, sep='/'), header = TRUE)
  
  Ca.i <<- vector()
  Mg.i <<- vector()
  k.run <<- vector()
  clim.sens <<- Clim.Sens <<- vector()
  degassing <<- vector()
  # pCO2 <<- vector()
  # temp.i <<- vector()
  Age <<- vector()
  output <<- vector(mode='list',length=length(sw.chem$Age))
  
  ## MODIFY BELOW TO INITIALIZE WITH PALEO CONDITIONS -- (Turned off for now)
  ## NOTE: some of these are not used in the calculations for certain scripts
  f=1    # set to run for one time step (use for-loop for Phanerozoic run instead)
  y=1   # time index in input data
  Age[f] <<- sw.chem$Age[y]
  Ca.i[f] <<- sw.chem$Ca[y]
  Mg.i[f] <<- sw.chem$Mg[y]
  degassing[f] <<- 1     # Berner.volc$Fvolc[y]
}



# -------------------------------------------------------------------- # 
# -------------------------------------------------------------------- # 
# -------------------------------------------------------------------- # 



## [4] PRINT a summary of the simulation
pre_run.summary <- function(df.force){
  # COLLECT THE NUMBER OF UNIQUE VALUES FOR EACH CASE 
  n.times.start <- length(unique(df.force$t.start.force))
  n.times.end <- length(unique(df.force$t.end.force))
  n.dur <- length(unique(df.force$duration.force))
  n.time <- max(n.times.start, n.times.end, n.dur)
  n.geog <- length(unique(df.force$geography.force))
  # n.sens.par <- length(unique(df.force$SENSITIVITY.force)) # we don't track sens par names here
  n.sens.val <- length(unique(df.force$SENSITIVITY.force))
  n.m <- length(unique(df.force$m.force))
  n.glwx <- length(unique(df.force$glwx.force))
  n.co2 <- length(unique(df.force$co2.force))
  
  force.values <- c("time grids" = n.time, "geographies" = n.geog, 
                    "param values" = n.sens.val, "init co2s" = n.co2, 
                    "carbon perturbations" = n.m, "glacial wx factors" = n.glwx)
  print.idx <- which(force.values > 1)
  
  # total sims
  n.sims <- nrow(df.force)
  n.summary <- paste(force.values[print.idx], names(force.values[print.idx]), collapse=' X ')
  
  # ------------------------------------- #
  # Print result
  cat(paste("=======================================================", 
            paste("TOTAL SIMULATIONS TO RUN:", n.sims),
            "**",
            paste("One simulation for each:", n.summary),
            "=======================================================", 
            sep='\n'))
  # ------------------------------------- #
}




