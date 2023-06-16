#-------------------------------------------------------------------------------#
#                Moist Energy Balance Model (MEBM)                              # 
#                                                                               #
# Multi-stability navigator --- this script includes functions that navigate    #
# the multi-stable solution space of MEBM by assigning different temperature    #
# boundary conditions                                                           #
#                                                                               #
# Following the solution of Gerard Roe, Univ. Washington (see Roe et al., 2015) #
# Hydroclimate parameterizations after Held, 2001; from Siler et al., (in Rev.) #
#      Coded and modified in R by Tyler Kukla, Stanford Univ. June, 2018        #
#-------------------------------------------------------------------------------#


# --- THE TIME LIMITER FOR MEBM 
try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf){
  y <- try({setTimeLimit(cpu, elapsed); expr}, silent = TRUE) 
  if(inherits(y, "try-error")) NULL else y 
}

# --- INSOLATION SELECTOR GIVEN TIME
insol_YearSelect <- function(t, insol_df = insol_df){
  insol_yrs <- unique(insol_df$yr_positive)
  yr.idx <- which.min(abs(t-insol_yrs))
  # pull out just the relevant yrs
  out.df <- insol_df[yr_positive == insol_yrs[yr.idx]]
  # return result
  return(out.df)
}


# --- GET THE MOST OFTEN OCCURRENCE IN VECTOR
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# --- GET COUNT OF INITIAL CLIM STATE (if you don't want to use mode but rather specificy count for selecting new state)
getcount <- function(iceDeal.count, iceDeal.last.count, test.vec.count){
  length_total <- length(test.vec.count)    # note: test.vec include original state plus all new states
  length_new <- length(which(test.vec.count == iceDeal.count))   # how many clim states equal the new one
  length_old <- length(which(test.vec.count == iceDeal.last.count))  # how many clim states equal the old one
  length_other <- length_total - length_new - length_old   # whatever is left (likely NAs for no solution)
  # output the result
  out.length <- c("total" = length_total, "new" = length_new, "old" = length_old, "other" = length_other) 
  return(out.length)
}

# --- PRESCRIBE THE NEW BOUNDARY CONDITION TEMPERATURES
TempBC_fun <- function(thisClim){
  if(thisClim == "south pole" | thisClim == "north pole"){
    LBtemp <<- LBtemp.options['one pole']
    RBtemp <<- RBtemp.options['one pole']
  } else if(thisClim == "both poles"){
    LBtemp <<- LBtemp.options['both poles']
    RBtemp <<- RBtemp.options['both poles']
  } else{
    LBtemp <<- LBtemp.options['no ice']
    RBtemp <<- RBtemp.options['no ice']
  }
}


# --- GET THE NEW ICE STATE
clim_state_ID <- function(thisdf){
  if(thisdf$temp_C[1] < ice_Threshold & thisdf$temp_C[length(xin)] < ice_Threshold){
    iceDeal <- "both poles"
  } else if(thisdf$temp_C[1] < ice_Threshold){
    iceDeal <- "south pole"
  } else if(thisdf$temp_C[length(xin)] < ice_Threshold){
    iceDeal <- "north pole"
  } else{iceDeal <- "no ice"}
  
  # spit out result
  return(iceDeal)
}

# --- TURN THE ICE STATE INTO A CODED INTEGER
clim_state_ID_INT <- function(this.state){
  out.int <- vector()  # initialize output vec
  for(i in 1:length(this.state)){
    if(this.state[i] == "no ice"){out.int[i] <- 0
    }else if(this.state[i] == "south pole"){out.int[i] <- 1
    }else if(this.state[i] == "north pole"){out.int[i] <- 2
    }else if(this.state[i] == "both poles"){out.int[i] <- 11
    }else{out.int[i] <- 9999}
  }
  return(out.int)
}

# --- TURN CODED STATE INTO STRING
clim_state_INT.2.Str <- function(this.int){
  out.str <- vector()   # initialize output vec
  for(i in 1:length(this.int)){
    if(this.int[i] == 0){out.str[i] <- "no_ice"
    }else if(this.int[i] == 1){out.str[i] <- "south_pole"
    }else if(this.int[i] == 2){out.str[i] <- "north_pole"
    }else if(this.int[i] == 11){out.str[i] <- "both poles"
    }else{out.str[i] <- "unknown_ERR"}
  }
  return(out.str)
}

## IF YOU GET A SNOWBALL OR NULL RESULT
# ... function steps through alternative solutions until
# ... an acceptable solution is found
# 
# ... df = output dataframe from MEBM
# ... SearchDir.BC = the search direction for temp boundary conditions after looking east and north (R is right/east; L is left/north)
# ... ExtraTimeFactor = the timeMax multiplier to allow trycatch to search a little longer 
Snowball_NULL_navigator <- function(thisdf, thisCO2, insol.t, SearchDir.BC="R", 
                                    step.size = 0.5, nsteps = 15, ExtraTimeFactor=3){
  # build function to test if the while loop is needed
  run.while_fun <- function(testdf){
    if(is.null(testdf)){
      out.while <- c("Y", "null")
    } else if(is.vector(testdf)){
      if(testdf == "break"){out.while <- "N"}
    } else if(max(testdf$temp_C, na.rm=T) <= ice_Threshold){
        out.while <- c("Y", "snowball")
    } else{out.while <- c("N")}
    
    return(out.while)
  }
  # what to do if you get a snowball or NULL result
  RBC.counter <- 1   # count the number of times we've iterated on the right boundary 
  LBC.counter <- 1   # count the number of times we've iterated on the left boundary 
  BOTH.counter <- 1  # count number of times we've iterated w/ both at once
  ALL.counter <- 1   # track all the attempts
  RBtemp.init <- RBtemp.save <- RBtemp  # save the initial conditions
  LBtemp.init <- LBtemp.save <- LBtemp
  # vector for directions left to analyze
  not.analyzed <- c("R", "L", "both")
  
  # test whether to begin the while loop
  run.while <- run.while_fun(thisdf)
  
  while(run.while[1] == "Y"){
    # -- FIRST ITERATION ONLY
    if(ALL.counter == 1){
      # [0] UPDATE USER 
      if(run.while[2] == "null"){print("NULL sol'n -- Searching for a real climate!")}
      if(run.while[2] == "snowball"){print("Fell into a snowball -- trying to get out!")}
      
      # [1] IF NULL, TRY JUST GIVING IT SOME MORE TIME
      if(run.while[2] == "null"){
        thisdf <- try_with_time_limit(expr = MEBM_TRAIN.sensitivity(pCO2_in = round(thisCO2, CO2_res), 
                                                                    LandFrac_tx = thisPaleogeo, get_pars=FALSE, 
                                                                    custom_insol = TRUE, insol_dat = insol.t),
                                      cpu = cpu.timeMax*ExtraTimeFactor, elapsed = elapsed.timeMax)
        LBtemp.save <- LBtemp # save in case this is the final result
        RBtemp.save <- RBtemp # save in case this is the final result
      }
      # update state
      run.while <- run.while_fun(thisdf)
      ALL.counter <- ALL.counter + 1 
    }
    
    # print("GOTCHA") # debug
    
    # if run.while still "Y", dive into the grid search
    if(run.while[1] == "Y"){
      # Set a counter so we only do each if statement once per loop
      # (this is critical because we update SearchDir.BC at end of "if"... 
      # we don't want it to cause us to slip into the next "if" and override progress)
      THIS.count <- ALL.counter
      # FIRST FOLLOW PRESCRIBED SEARCH DIR 
      if(SearchDir.BC == "R" & THIS.count == ALL.counter){
        RBtemp <<- RBtemp + step.size
        LBtemp <<- LBtemp.init
        thisdf <- try_with_time_limit(expr = MEBM_TRAIN.sensitivity(pCO2_in = round(thisCO2, CO2_res), 
                                                                    LandFrac_tx = thisPaleogeo, get_pars=FALSE, 
                                                                    custom_insol = TRUE, insol_dat = insol.t),
                                      cpu = cpu.timeMax, elapsed = elapsed.timeMax)
        RBC.counter <- RBC.counter + 1
        ALL.counter <- ALL.counter + 1
        LBtemp.save <- LBtemp # save in case this is the final result
        RBtemp.save <- RBtemp # save in case this is the final result
        # check whether to switch directions
        if(RBC.counter > nsteps){
          # update not.analyzed
          not.analyzed <- not.analyzed[not.analyzed != "R"]
          SearchDir.BC <- ifelse(length(not.analyzed) == 0, "DONE", not.analyzed[1]) # update search.dir to move along
          # reset temp BCs
          LBtemp <<- LBtemp.init
          RBtemp <<- RBtemp.init
        }
      } else if(SearchDir.BC == "L" & THIS.count == ALL.counter){
        RBtemp <<- RBtemp.init
        LBtemp <<- LBtemp + step.size
        thisdf <- try_with_time_limit(expr = MEBM_TRAIN.sensitivity(pCO2_in = round(thisCO2, CO2_res), 
                                                                    LandFrac_tx = thisPaleogeo, get_pars=FALSE, 
                                                                    custom_insol = TRUE, insol_dat = insol.t),
                                      cpu = cpu.timeMax, elapsed = elapsed.timeMax)
        LBC.counter <- LBC.counter + 1
        ALL.counter <- ALL.counter + 1
        LBtemp.save <- LBtemp # save in case this is the final result
        RBtemp.save <- RBtemp # save in case this is the final result
        # check whether to switch directions
        if(LBC.counter > nsteps){
          # update not.analyzed
          not.analyzed <- not.analyzed[not.analyzed != "L"]
          SearchDir.BC <- ifelse(length(not.analyzed) == 0, "DONE", not.analyzed[1]) # update search.dir to move along
          # reset temp BCs
          LBtemp <<- LBtemp.init
          RBtemp <<- RBtemp.init
        }
      } else if(SearchDir.BC == "both" & THIS.count == ALL.counter){
        RBtemp <<- RBtemp + step.size
        LBtemp <<- LBtemp + step.size
        thisdf <- try_with_time_limit(expr = MEBM_TRAIN.sensitivity(pCO2_in = round(thisCO2, CO2_res), 
                                                                    LandFrac_tx = thisPaleogeo, get_pars=FALSE, 
                                                                    custom_insol = TRUE, insol_dat = insol.t),
                                      cpu = cpu.timeMax, elapsed = elapsed.timeMax)
        ALL.counter <- ALL.counter + 1
        BOTH.counter <- BOTH.counter + 1
        LBtemp.save <- LBtemp # save in case this is the final result
        RBtemp.save <- RBtemp # save in case this is the final result
        # check whether to switch directions
        if(BOTH.counter > nsteps){
          # update not.analyzed
          not.analyzed <- not.analyzed[not.analyzed != "both"]
          SearchDir.BC <- ifelse(length(not.analyzed) == 0, "DONE", not.analyzed[1]) # update search.dir to move along
          # reset temp BCs
          LBtemp <<- LBtemp.init
          RBtemp <<- RBtemp.init
        }
      } else if(SearchDir.BC == "DONE" & THIS.count == ALL.counter){
        thisdf <- "break" # give up and select the last climate state
      }
    }
    
    # update run.while state
    run.while <- run.while_fun(thisdf)
    # print(paste0(c("STILL STUCK? ", "WHY? "),run.while, collapse = ' -- '))
    
  }
    
  # set LB and RB temp to the save solution
  LBtemp <<- LBtemp.save
  RBtemp <<- RBtemp.save
  # return the result
  # print("Out of the snowball!")
  return(thisdf)
}





## IF YOU GET A NEW CLIMATE STATE AND WANT TO MAKE SURE IT'S ROBUST
# ... function steps through adjacent solutions 
# ... other function takes result and decides what the climate state should be

## [1] STEP THROUGH NEW CLIM STATES (to be nested in other fxn)
# values defined in navigator fxn except: 
# cancel.early = determines whether we quit as soon as we violate nav.count.conditions
NewClimState_step <- function(in.df, in.CO2, in.iceDeal, in.iceDeal.last, in.insol.t, stepSize, 
                              nav.count.conditions, include.diagonals, cancel.early){
  # set the initial test (the result of recent clim state)
  LBtemp.init <- LBtemp; RBtemp.init <- RBtemp
  testing <- TRUE    # we are testing to make sure it's a robust result
  thisTest <- 1
  test.vec <- vector()
  test.vec[1] <- in.iceDeal
  outdf.list <- list()     # make a list for all the climate possibilities
  outdf.list[[1]] <- in.df
  while(testing == T){
    # set the direction of travel in the solution space 
    if(thisTest == 1){LBtemp <<- LBtemp.init ; RBtemp <<- RBtemp.init + stepSize   # increase RBtemp
    } else if(thisTest==2){LBtemp <<- LBtemp.init + stepSize; RBtemp <<- RBtemp.init  # increase LBtemp
    } else if(thisTest==3){LBtemp <<- LBtemp.init ; RBtemp <<- RBtemp.init - stepSize  # decrease RBtemp
    } else if(thisTest==4){LBtemp <<- LBtemp.init - stepSize; RBtemp <<- RBtemp.init # decrease LBtemp
    } else if(thisTest==5){LBtemp <<- LBtemp.init + stepSize ; RBtemp <<- RBtemp.init + stepSize  # increase both
    } else if(thisTest==6){LBtemp <<- LBtemp.init - stepSize ; RBtemp <<- RBtemp.init - stepSize  # decrease both
    } else if(thisTest==7){LBtemp <<- LBtemp.init + stepSize ; RBtemp <<- RBtemp.init - stepSize  # increase left decrease right
    } else if(thisTest==8){LBtemp <<- LBtemp.init - stepSize ; RBtemp <<- RBtemp.init + stepSize} # increase right decrease left
    # run the simulation
    thisdf.x <- try_with_time_limit(expr = MEBM_TRAIN.sensitivity(pCO2_in = round(in.CO2, CO2_res), 
                                                                  LandFrac_tx = thisPaleogeo, get_pars=FALSE, 
                                                                  custom_insol = TRUE, insol_dat = in.insol.t),
                                    cpu = cpu.timeMax, elapsed = elapsed.timeMax)
    # print("TEST X") # debug
    # thisdf.x <- try_with_time_limit(expr = MEBM_TRAIN(pCO2_in = round(in.CO2, CO2_res),  
    #                                                   LandFrac_tx = thisPaleogeo, custom_insol=FALSE, insol_dat=NULL),
    #                                 cpu = timeMax, elapsed = timeMax)
    outdf.list[[thisTest + 1]] <- thisdf.x
    # get the new ice deal
    if(!is.null(thisdf.x)){
      if(thisdf.x$temp_C[1] < ice_Threshold & thisdf.x$temp_C[length(xin)] < ice_Threshold){
        iceDeal.x <- "both poles"
      } else if(thisdf.x$temp_C[1] < ice_Threshold){
        iceDeal.x <- "south pole"
      } else if(thisdf.x$temp_C[length(xin)] < ice_Threshold){
        iceDeal.x <- "north pole"
      } else{iceDeal.x <- "no ice"}
    } else{iceDeal.x <- NA}
    # save the new ice deal
    test.vec[thisTest + 1] <- iceDeal.x
    thisTest <- thisTest + 1
    
    # test if we can abort the rest of the runs
    # do we have too many olds or too many others? 
    this.count <- getcount(iceDeal.count = in.iceDeal, iceDeal.last.count = in.iceDeal.last, test.vec.count = test.vec)
    if(cancel.early == TRUE){
      if((this.count["old"] > nav.count.conditions["old"]) | (this.count["other"] > nav.count.conditions["other"])){
        testing <- FALSE
      }
    }
    if(thisTest >= 5 & include.diagonals == FALSE){  # finish our testing once we've explored enough
      testing <- FALSE 
    } else if(thisTest >= 9 & include.diagonals == TRUE){ # finish our testing once we've explored enough
      testing <- FALSE
    }
  }
  
  # return results
  returnMe <- list(test.vec, outdf.list)
  return(returnMe)
}


# NAVIGATE TO A NEW CLIMATE STATE (or not)
# ... thisdf = output dataframe from MEBM
# ... thisCO2 = atmospheric CO2 for this step
# ... iceDeal.last = ice configuration at previous timestep
# ... iceDeal = proposed ice configuration for this timestep
# ... TieGoes2the = (only matters for nav.metric = "mode") "new" | "old" (whether a tie in the mode goes to the new or old state)
# ... clim.stepSize = [degrees C] step to search for other solutions in diff directions
# ... nav.metric = "count" | "mode" (whether we base deciding to change climate state on the count of new, old, other clim states, or on the most often occurring state)
# ... this.nav.count.conditions = [int vector] must have names "new" ; "old" ; "other" with values to compare test results. Clim state changes if more or equal to "new" & less or equal to "old" & less or equal to "other"
# ---   -- EXAMPLE: nav.count.conditions =  c("new"=3, "old"=0, "other"=1e3) ; here, "other" is really high so that it essentially doesn't matter in the decision (solutions will always be lower or equal to this number)
# ... this.include.diagonals = "TRUE" | "FALSE" (whether search should include the steps to the four diagonals (in addition to up down left right steps))
# ... CO2.check = "TRUE" | "FALSE" (whether to step through greater / lesser CO2 if the new clim state is accepted)
# ... CO2.step: [int, ppm CO2] how much to step up and down if CO2.check = TRUE
NewClimState_Navigator <- function(thisdf, thisCO2, iceDeal, iceDeal.last, insol.t, TieGoes2the = "new", 
                                   clim.stepSize = 1, nav.metric = "count", 
                                   this.nav.count.conditions = c("new"=4, "old"=0, "other"=1e3), 
                                   this.include.diagonals = TRUE, CO2.check = TRUE, CO2.step = 0.5){
  iceDeal.new <- iceDeal # define new state
  ## --- FIRST check that things will run properly --- ## 
  # if we were (accidentally) looking for more than 100% agreement, then set nav.count to total test count
  if(nav.metric == "count"){
    if(!all(c("new", "old", "other") %in% names(this.nav.count.conditions))){
      stop("this.nav.count.conditions must have names 'new', 'old', and 'other' (some or all missing)")
    }
    if(this.nav.count.conditions["old"] < 0){this.nav.count.conditions["old"] <- 0} # min "old" at 0
    if(this.nav.count.conditions["other"] < 0){this.nav.count.conditions["other"] <- 0} # min "old" at 0
    if(this.include.diagonals==TRUE & this.nav.count.conditions["new"] > 9){ # max new is 9
      this.nav.count.conditions["new"] <- 9
    } else if(this.include.diagonals==FALSE & this.nav.count.conditions["new"] > 5){
      this.nav.count.conditions["new"] <- 5
    }
  }
  
  ## --- CHECK IF THE CLIMATE STATE HAS CHANGED AND IF WE ACCEPT IT --- ##
  if(iceDeal != iceDeal.last){ # then we have a change in state and we wish to check that it's robust
    print(paste("New climate state! Testing if I should switch out of:", iceDeal.last, sep=' '))
    ## WALK THROUGH RBTEMP / LBTEMP SPACE
    result.step <- NewClimState_step(in.df = thisdf, in.CO2 = thisCO2, in.iceDeal = iceDeal, 
                                     in.iceDeal.last = iceDeal.last, in.insol.t = insol.t, stepSize = clim.stepSize, 
                                     nav.count.conditions = this.nav.count.conditions, 
                                     include.diagonals = this.include.diagonals, cancel.early = TRUE)
    test.vec <- result.step[[1]]
    outdf.list <- result.step[[2]]
    # find the most common solution
    iceDeal.mode <- getmode(na.omit(test.vec))
    iceDeal.count <- getcount(iceDeal.count = iceDeal, iceDeal.last.count = iceDeal.last, test.vec.count = test.vec)
    
    # if we test for new state by MODE ------ 
    if(nav.metric == "mode" & (iceDeal %in% iceDeal.mode)){
      if(length(iceDeal.mode > 1) & TieGoes2the=="old" & (iceDeal.last %in% iceDeal.mode)){ # if there's a tie and we decide to give the tie to the old state (thus the old state is part of the tie) then update to maintain the old state 
        iceDeal <- iceDeal.last
        thisdf <- outdf.list[[which(test.vec==iceDeal.last)[1]]]
        print(paste("Staying in old climate state:", iceDeal.last, sep=' '))
      } else{
        thisdf <- thisdf  # set df to the original value if initial iceDeal is the most common result (even if it's a tie)
        # print(paste("New climate state accepted:", iceDeal, sep=' '))
      }
      
    } else if(nav.metric == "mode" & (is.null(iceDeal.mode) | is.na(iceDeal.mode))){  # these shouldn't happen, but in case they do this is in here
      thisdf <- thisdf
    } else if(nav.metric == "mode" & (iceDeal %in% iceDeal.mode)==FALSE){   # then select the first one (whether there's one option or more)
      iceDeal <- iceDeal.mode[1]
      # ASSIGN THE NEW DATA
      thisdf <- outdf.list[[which(test.vec==iceDeal.mode[1])[1]]]
    }
    
    # if we test for new state by COUNT -----
    if(nav.metric == "count"){
      # test if we meet the conditions
      # same or more new than specified ; same or less old than specified ; same or less other than specified
      if((iceDeal.count["new"] >= this.nav.count.conditions["new"]) & (iceDeal.count["old"] <= this.nav.count.conditions["old"]) & (iceDeal.count["other"] <= this.nav.count.conditions["other"])){
        thisdf <- thisdf  # ...keep the new ice/clim state
        # print(paste("New climate state accepted:", iceDeal, sep=' '))
      } else if(length(which(test.vec==iceDeal.last)) >= 1){ # otherwise keep the old climate state if it's an option
        iceDeal <- iceDeal.last
        thisdf <- outdf.list[[which(test.vec==iceDeal.last)[1]]]  # revert to the old one
        print(paste("Staying in old climate state:", iceDeal.last, sep=' '))
      } else{  # if we don't meet conditions AND old climate state isn't possible, just take the first solution
        thisdf <- thisdf
      }
    }
  } else{thisdf <- thisdf}   # if iceDeal == iceDeal.last then we keep the solution
  # print(paste("All done! my clim state is now:", iceDeal, sep=' '))
  
  
  # # ---- LOOP OVER CO2? ---- # # 
  # Check if we should loop over CO2s 
  if(iceDeal != iceDeal.last & iceDeal == iceDeal.new & CO2.check == FALSE){
    print(paste("New climate state accepted:", iceDeal, sep=' '))
  } else if(iceDeal != iceDeal.last & iceDeal == iceDeal.new & CO2.check == TRUE){  # in this case, we have accepted the new clim state per the above test 
    # add CO2 step 
    result.step.1  <- NewClimState_step(in.df = thisdf, in.CO2 = (thisCO2 + CO2.step), in.iceDeal = iceDeal, 
                                       in.iceDeal.last = iceDeal.last, in.insol.t = insol.t, stepSize = clim.stepSize, 
                                       nav.count.conditions = this.nav.count.conditions, 
                                       include.diagonals = this.include.diagonals, cancel.early = FALSE)
    test.vec.1 <- result.step.1[[1]]
    outdf.list.1 <- result.step.1[[2]]
    # subtract CO2 step 
    result.step.2  <- NewClimState_step(in.df = thisdf, in.CO2 = (thisCO2 - CO2.step), in.iceDeal = iceDeal, 
                                        in.iceDeal.last = iceDeal.last, in.insol.t = insol.t, stepSize = clim.stepSize, 
                                        nav.count.conditions = this.nav.count.conditions, 
                                        include.diagonals = this.include.diagonals, cancel.early = FALSE)
    test.vec.2 <- result.step.1[[1]]
    outdf.list.2 <- result.step.1[[2]]
    
    # compare results
    test.vec.co2 <- c(test.vec.1, test.vec.2)
    outdf.list.co2 <- c(outdf.list.1, outdf.list.2)
    # find the most common solution
    iceDeal.mode.co2 <- getmode(na.omit(test.vec.co2))
    iceDeal.count.co2 <- getcount(iceDeal.count = iceDeal, iceDeal.last.count = iceDeal.last, test.vec.count = test.vec.co2)
    
    # just check to see if the olds show up 
    if(iceDeal.count.co2["old"] > (this.nav.count.conditions["old"] * 2)){ # reject new state
      iceDeal <- iceDeal.last
      thisdf <- outdf.list[[which(test.vec==iceDeal.last)[1]]]  # revert to the old one 
      if(is.null(thisdf)){ # if there was no old clim state in the original list, find one in the updated list (albeit w/ slightly diff co2)
        thisdf <- outdf.list.co2[[which(test.vec.co2 == iceDeal.last)[1]]]
      }
      print(paste("Staying in old climate state:", iceDeal.last, sep=' '))
    } else{ print(paste("New climate state accepted:", iceDeal, sep=' ')) } # keep the new state, do nothing
  } 
  
  # return the result
  return(thisdf)
}


