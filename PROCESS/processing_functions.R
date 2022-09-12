# ------------------------------------------------ # 
# Some functions to process CH2O-CHOO TRAIN output #
# ---                                              #
# T Kukla (Colostate Univ. 2022)                   #
# ------------------------------------------------ #

## [1] COLLECT AND SYNTHESIZE OUTPUT
collect.dat <- function(results.path, wd = getwd()){
  parent.path.split <- strsplit(wd, '/')[[1]]
  parent.path.idx <- which(parent.path.split == 'CODE_DISTRIBUTE')
  parent.path <- paste(paste0(parent.path.split[c(1:parent.path.idx)], collapse='/'), results.path, sep='/')
  
  ## ... Read in synthesized results
  suffixes <- list('init' = 'INITCLIM.RDS',    # MEBM results for t = 0 for each simulation
                   'full' = 'FULL.RDS',        # CH2O-CHOO TRAIN timeseries results for each simulation
                   'summary' = 'SUMMARY.RDS',
                   'clim' = 'CLIM_TS')  # global mean summary results for each simulation
  fns <- list.files(parent.path, pattern='.RDS', full.names = FALSE)  # filenames
  fn.idx <- lapply(X = suffixes, FUN = function(x) which(str_detect(fns, x)))  # list of filename indices
  # define data.tables
  df.init <- readRDS(paste(parent.path, fns[fn.idx$init], sep='/'))  # init clim file
  df.ts <- readRDS(paste(parent.path, fns[fn.idx$full], sep='/'))    # timeseries file
  df.sum <- readRDS(paste(parent.path, fns[fn.idx$summary], sep='/'))# summary file
  
  ## ... collect clim timeseries data if available
  for(i in 1:length(fn.idx$clim)){  # loop thru iters and save together
    tmp.clim <- readRDS(paste(parent.path, fns[fn.idx$clim[i]], sep='/'))
    if(i == 1){
      df.climts <- tmp.clim
    } else{
      df.climts <- as.data.table(rbind(df.climts, tmp.clim))
    }
  }
  
  ## ... return as list 
  out.list <- list(df.ts, df.init, df.climts, df.sum)
  names(out.list) <- c('timeseries', 'initial_clim', 'timeseries_clim', 'summary')
  return(out.list)
}

## [2] FUNCTION TO SAVE FIGURES
plot.save <- function(save.name, this.plot, 
                      save.path = save.folder,
                      results.path = RESULTS.FOLDER, cp.size = FALSE, # cp.size = whether the size is crossplot or not
                      wd = getwd(), save_or_not = permit.save, 
                      sv.wd = save.width, sv.ht = save.height){
  if(save_or_not == FALSE){ # don't even bother 
    
  } else{
    parent.path.split <- strsplit(wd, '/')[[1]]
    parent.path.idx <- which(parent.path.split == 'CODE_DISTRIBUTE')
    parent.path <- paste(paste0(parent.path.split[c(1:parent.path.idx)], collapse='/'), results.path, sep='/')
    
    # create new save dir if it doesn't exist
    saveDir <- paste(parent.path, save.path, sep='/')
    if(!file.exists(saveDir)){dir.create(file.path(saveDir))}
    
    # correct size for crossplot if needed
    if(cp.size == TRUE){
      sv.wd <- save.width.cp
      sv.ht <- save.height.cp
    }
    
    # save the figure
    save.x <- paste(saveDir, save.name, sep='/')
    ggsave(save.x, this.plot, width=sv.wd, height=sv.ht, units='cm')
  }
}


## [3] INITIALIZE CLIMATE MODEL SIMULATIONS
prep.MEBM <- function(results.path, wd = getwd()){
  # find the results file
  p.path.split <- strsplit(wd, '/')[[1]]
  p.path.idx <- which(p.path.split == 'CODE_DISTRIBUTE')
  p.path <- paste(paste0(p.path.split[c(1:p.path.idx)], collapse='/'), results.path, sep='/')
  
  # find and source the forcing file
  forcing.iters.name <- list.files(p.path, pattern = 'forcing_iters', full.names = FALSE)[1]
  FORCE.DF <<- readRDS(paste(p.path, forcing.iters.name, sep='/'))
  
  # find and source the parameter lists
  par.list.names <- list.files(p.path, pattern = 'PAR-LIST.RDS', full.names = TRUE)
  parList.list <- list()
  for(j in 1:length(par.list.names)){
    parList.list[[j]] <- readRDS(par.list.names[j])
  }
  names(parList.list) <- dat$summary$RunName
  # send to global env. 
  parList.list <<- parList.list
  
  # find and source the config file
  # !! defaults to first config file if there are multiple !! 
  config.name <- list.files(p.path, pattern='copy_of_', full.names = FALSE)[1]
  source(paste(p.path, config.name, sep='/'))
  
  # read in MEBM code and c-cycle (for wx)
  source(code.init.fxn.loc)
  initialize.mebm()      # moist energy balance model functions
  initialize.c.cycle()   # carbon cycle and weathering functions
}


## [4] CALCULATE CLIMATE ANOMALIES
## ... anomaly relative to t == 0
clim.anoms <- function(zonal.clim, anoms = c('temp_C', 'Precip_m_yr', 'Discharge_m3_yr', 'Fwsil')){
  n.runs <- unique(zonal.clim$RUN.idx)
  for(i in n.runs){
    this.run <- zonal.clim[RUN.idx == n.runs[i]]
    for(j in 1:length(anoms)){
      this.var <- anoms[j]
      this.anom <- unlist(this.run[ , ..this.var]) - unlist(this.run[time == 0 , ..this.var])
      # bring together
      if(j == 1){
        out.anom <- as.data.table(this.anom)
      } else{
        out.anom <- as.data.table(cbind(out.anom, this.anom))
      }
    }
    # name the df
    colnames(out.anom) <- paste0(anoms, '_anom', each='')
    # combine with original dat
    if(i ==1){
      out.clim <- as.data.table(cbind(this.run, out.anom))
    } else{
      tmp.outclim <- as.data.table(cbind(this.run, out.anom))
      out.clim <- as.data.table(rbind(out.clim, tmp.outclim))
    }
  }
  # return result
  return(out.clim)
}


## [5] RUN THE MEBM FOR DEFINED TIMESTEPS
CLIM.RUN <- function(force.iters = FORCE.DF){
  # some code repeated from train.tracks
  # check whether geography is a list of many geogs, or a single geog
  if(is.null(nrow(geog.list))){  # if it doesn't have rows, then the vec is the list
    multi.geog <- TRUE
    geog.vec <- names(geog.list)
  } else{   # if it does have rows, then the vec is just this file
    multi.geog <- FALSE
    geog.vec <- LandFracFileName 
  }
  # number of iters
  n.idx <<- nrow(force.iters)
  
  # find initial clim
  ctrl.iter <- which(force.iters$CTRL.idx == 'Y')
  CTRL.CLIM <- dat$initial_clim[RUN.idx == ctrl.iter]
  runoff.INIT <<- round(CTRL.CLIM$budykoRunoff_m_yr, q_res)   # latitudinal runoff profile
  Temp.bin.INIT <<- round(CTRL.CLIM$temp_C, temp_res)
  
  # create a list to save results
  save.list <- list()
  
  # ... LOOP THROUGH
  for(i in 1:length(parList.list)){
    # ************************************ # 
    # TROUBLESHOOTING PRINT                #
    cat(paste("=======================================================",
              print(paste("SOLVING ITER:", i)),    #
              "=======================================================", 
              sep='\n'))
    # ************************************ # 
    # assign this parList 
    parList <<- parList.list[[i]]
    # get number of timesteps
    this.dat <- dat$timeseries[iter == i]
    this.initclim <- dat$initial_clim[RUN.idx == i]
    timesteps <- this.dat$time
    
    # assign weathering pars
    p <<- parList$ParamFunc_pars
    Ea.wx <- parList$ParamFunc_pars['Ea']
    Ts.wx <- parList$ParamFunc_pars['Ts']
    soil.CO2.flag.wx <- parList$ParamFunc_pars['soil.CO2.flag']
    pCO2.i <<- force.iters$co2.force[i]
    temp.a.i <<- round(weighted.mean(x = this.initclim$temp_C, w=this.initclim$Land_Area_m2), temp_res) 
    
    # change init if needed
    if(scale.to.init.co2 == FALSE){
      runoff.INIT <<- round(this.initclim$budykoRunoff_m_yr, q_res) #  this is used for MCfunc
      Temp.bin.INIT <<- round(this.initclim$temp_C, temp_res)
    }
    
    # ... LOOP timesteps
    for(j in 1:length(timesteps)){
      cat(paste("solving", j, "of", length(timesteps)))
      ## [1] Set geography
      if(multi.geog == TRUE){ # if we use more than one geography, test the one set by force.df
        thisPaleogeo <<- geog.list[[force.iters$geography.force[i]]]
        landfrac_index <<- geog.vec[i]
      } else{thisPaleogeo <<- geog.list; landfrac_index <<- geog.vec} # otherwise just use the single geog
      
      ## [2] set initial conditions
      pCO2 <<- this.dat$pCO2[j]
      # ... set correct sensitivity param
      parList[[par.here]][[SENSITIVITY.PAR]] <- force.iters$SENSITIVITY.force[i]
      parList <<- parList  # I reset this because there's an error when I try to change the global environment of parList (<<-) with the line above
      
      ## [3] Set climate state
      LBtemp <<- this.dat$LBtemp[j] ; RBtemp <<- this.dat$RBtemp[j]
      
      ## [4] initial mebm 
      df_MEBM <- try_with_time_limit(expr = MEBM_TRAIN.sensitivity(pCO2_in = round(pCO2, CO2_res), LandFrac_tx = thisPaleogeo, get_pars=FALSE ),
                                     cpu = cpu.timeMax, elapsed = elapsed.timeMax)
      # ... check if we have a result
      df_MEBM <- Snowball_NULL_navigator(thisdf = df_MEBM, thisCO2 = pCO2, SearchDir.BC = "R")
      # ... assign output 
      temp.lat <- round(df_MEBM$temp_C, temp_res)                 # [degC] MEBM temperature
      runoff.lat <-  round(df_MEBM$budykoRunoff_m_yr, q_res)   # [m yr-1] runoff calculated from Budyko framework (omega assigned in MEBM_TRAIN)
      # ... find iceline
      ice.latitudes <- which(df_MEBM$albedo >= iceAlf)
      if(length(ice.latitudes)== 0){iceline.lat <- 9999
      } else{iceline.lat <- min(abs(df_MEBM[ice.latitudes,]$Lat))}
      
      
      # ... set up vectors to fill in for loop
      Fwsil.vec <- numeric(length=length(xin))
      Fwcarb.vec <- numeric(length=length(xin))
      Conc.vec <- numeric(length=length(xin))
      Conc0.vec <- numeric(length=length(xin))
      Conc.carb <- numeric(length=length(xin))
      Conc0.carb <- numeric(length=length(xin))
      Dw.vec <- numeric(length=length(xin))
      area.frac.vec <- df_MEBM$Frac_of_tot_LandArea
      area.vec <- df_MEBM$Land_Area_m2
      # ------------------- DELETE?? -------------------- #
      # area.vec=rep(0.8,length(xin)) #this is the land area parameter -- need to fill in from input
      # total.area=sum(area.vec)#sum(area.vec)*length(xin)
      # ------------------------------------------------- #
      for(x.lat in 1:length(xin)){
        Temp.bin <- temp.lat[x.lat] # temperature in latitude bin
        
        # Runoff    
        q <- runoff.lat[x.lat]
        
        # Using a Maher and Chamberlain-style concentration feedback with a Winnick and Maher 
        # open-system calculation for Ceq
        Conc.results <- MCfunc.trans(pCO2 = pCO2, q0 = runoff.INIT[x.lat], q = q, Temperature = Temp.bin,
                                     Ea = Ea.wx, Ts=Ts.wx, soilCO2.use=soil.CO2.flag.wx,
                                     Temp.init = Temp.bin.INIT[x.lat])
        
        # Silicate Weathering
        Conc0.vec[x.lat] <- Conc.results$Granite$Conc0*parList$ParamFunc_pars['silicate.lith'] + Conc.results$Basalt$Conc0*(1-parList$ParamFunc_pars['silicate.lith'])
        Conc.vec[x.lat] <- Conc.results$Granite$Conc*parList$ParamFunc_pars['silicate.lith'] + Conc.results$Basalt$Conc*(1-parList$ParamFunc_pars['silicate.lith'])
        Dw.vec[x.lat] <- Conc.results$Granite$Dw*parList$ParamFunc_pars['silicate.lith'] + Conc.results$Basalt$Conc*(1-parList$ParamFunc_pars['silicate.lith'])
        
        Fwsil.vec[x.lat] <- q*Conc.vec[x.lat]*area.vec[x.lat]*this.dat$Fwsil_scalar[j]     # NOTE: q is not indexed by "x.lat" because it's assigned at beginning of loop
        
        # Carbonate Weathering
        Conc0.carb[x.lat] <- Conc.results$Carbonate$Conc0
        Conc.carb[x.lat] <- Conc.results$Carbonate$Conc
        Fwcarb.vec[x.lat] <- q*Conc.carb[x.lat]*area.vec[x.lat]*this.dat$Fwcarb_scalar[j]
      }
      # add weathering fluxes to MEBM output 
      df_MEBM$Fwsil <- Fwsil.vec
      df_MEBM$Fwcarb <- Fwcarb.vec
      # add discharge since we didn't earlier
      df_MEBM$Discharge_m3_yr <- df_MEBM$budykoRunoff_m_yr * df_MEBM$Land_Area_m2
      # add identifiers
      df_MEBM$Runname <- dat$summary$RunName[i]
      df_MEBM$RUN.idx <- i
      df_MEBM$time <- this.dat$time[i]
      
      # STORE RESULT
      if(j == 1){
        out.df <- df_MEBM
      } else{
        out.df <- as.data.table(rbind(out.df, df_MEBM))
      }
    }
    
    # save in list
    save.list[[i]] <- out.df
  }
  
  # name and return save.list
  names(save.list) <- dat$summary$RunName
  return(save.list)
}




