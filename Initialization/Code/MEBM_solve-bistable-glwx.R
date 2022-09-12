#-------------------------------------------------------------------------------#
#                Moist Energy Balance Model (MEBM) -- Solver                    # 
#                                                                               #
# Following the solution of Gerard Roe, Univ. Washington (see Roe et al., 2015) #
# Hydroclimate parameterizations after Held, 2001; from Siler et al., (2018)    #
#      Coded and modified in R by Tyler Kukla, Stanford Univ. June, 2018        #
#-------------------------------------------------------------------------------#


## Modified to read in glwx variable -- 
# ... glwx is a coefficient from 0 to 1 that determines the fraction of runoff (P minus E)
# ... that is allowed beneath an ice sheet (0 is no runoff, thus no weathering;
# ... and 1 is full runoff, thus weathering as if there's no ice sheet)

MEBM_solve <- function(xin_mod, LBC_in, RBC_in, budy_om_in, mypars){
#---------------------------------CALCULATE THE MODEL SOLUTION------------------------------------------#
  #... run the MEBM
        # yini ; yend -- have structure =c(y1=temp, y2=flux)
        # x -- grid
        # func -- bvp function
        # guess -- initial guess for temperature profile
        # atol -- calibrated based on Matlab output of similar model and computation time
        # method -- "lsoda" is fastest... if specifying a jacobian matrix this might not be the fastest 
  print(
    system.time(
      # sol <- bvpshoot(yini = c(NA, LBC), yend = c(NA, RBC), x = xin, func = MEBMfun,guess=c(5),
      #                 verbose=F, atol=1e-4, parms = pars_in, maxiter = 2000, method="lsode") # ,maxsteps=10^6)
      # sol <- bvptwp(yini = c(NA, LBC), yend = c(NA, RBC), x = xin, func = MEBMfun,atol=1e-14,
      #               nmax=500) #yguess=c(temp_guess),
      sol <- bvpcol(yini = c(NA, LBC_in), yend = c(NA, RBC_in), x = xin_mod, func = MEBMfun,
                    parms = mypars, atol=3, nmax=1000000, xguess=c(min(xin_mod), max(xin_mod)),
                    yguess= matrix(nrow = 2, ncol = 2, data=c(LBtemp,RBtemp,0, 0), byrow=T))
      
      # .. old
      #yguess=c(temp_guess),
      #                 verbose=F, atol=1e-14) #, maxiter = 1000)#, method="lsoda")#,maxstep=10^5)
      
    )
  )
  #diagnostics(sol)
  
  #... solve the hydrologic cycle 
  myTemp <- sol[,2] 
  hydro <- hydrofun(temp=myTemp, x=xin_mod)
  
  #... re-calculate MEBM variables
  alf_flag <- mypars['alf_flag']
  # Solve the Top of Atmosphere (TOA) fluxes
  # Source  (ASR)
  if(alf_flag == "D18"){  # calculate albedo using the hyperbolic formulation of Dortmans et al., 2018 (Clim of the Past)
    # omega value determines how smooth the ice - land transition 
    #... (larger = shallower; smaller = steeper)
    alf <- 0.5 * ( (afun(xin_mod) + iceAlf) + (afun(xin_mod) - iceAlf) * tanh((myTemp)/omega) )
    }
  
  if(alf_flag == "step"){ # if the temperature is below the ice-threshold, declare it ice
    # Classic ice-threshold value is -10 degC (goes back to calibrations of Budyko)
    ice_dex <- which(myTemp <= ice_Threshold)    # find the ice line
    alf <- afun(xin_mod) ; alf[ice_dex] <- iceAlf    # calculate albedo
  }
  
  # Calculate the source term
  A0_in <- as.numeric(mypars['A0'])
  B0_in <- as.numeric(mypars['B0'])
  Q0_in <- as.numeric(mypars['Q0'])
  
  if(insol_flag == "INSOL"){
    Src <- InsolFun(xin_mod) * (1-alf)
  } else if(insol_flag == "Traditional"){
    Src <- Q0_in*(1-0.241*(3*xin_mod^2-1))*(1-alf)   # [W m-2] follows the 2nd order Legendre polynomial approximation
  }
  
  # Sink
  Snk <- (A0_in+B0_in*myTemp)                     # [W m-2] OLR 
  
  # MSE divergence 
  divF  <- 1/(2*pi*Re^2)*pracma::gradient(hydro$mseFlux_W,xin)
  divF[1:2] <- divF[3]
  divF[(length(xin)-1):length(xin)] <- divF[length(xin)-2]
  
  # calculate continental runoff with Budyko 
  evap.temp <- hydro$Evap_m_yr ; precip.temp <- hydro$Precip_m_yr
  # .. if evap or precip are negative, set to basically zero (zero or less returns NAs)
  evap.temp[which(evap.temp < 0)] <- 1e-2 ; precip.temp[which(precip.temp < 0)] <- 1e-2
  # get krun and runoff
  krun <- ((evap.temp/precip.temp) - (1 + (evap.temp/precip.temp)^budy_om_in)^(1/budy_om_in))*-1    # budyko equation (Fu, 1981)
  runoff_out <- krun*precip.temp         # applies budyko to all land (alternatively, we could just apply it to E>P)
  runoff_out <- ifelse(alf==iceAlf, glwx*runoff_out, runoff_out)    # set runoff to glwx*runoff if there is an ice sheet
  
  # entropy production (integrate divergence)
  sdot = trapz(xin, divF)    
  
  #... generate an output dataframe 
  myLat <- rad2deg(asin(xin))
  mydf <- as.data.table(cbind(myLat, hydro, runoff_out, alf, Src, Snk, divF, sdot))
  colnames(mydf) <- c("Lat", "temp_C", "MSE_J_kg", "Precip_m_yr", "Evap_m_yr", "EminusP_m_yr", "mseFlux_W", 
                      "latentFlux_W", "sensibleFlux_W", "latentEddyFlux_W", "latentHadleyFlux_W", "LHflux_Div", 
                      "budykoRunoff_m_yr", "albedo", "ASR_W_m2", "OLR_W_m2", "MSEfluxDiv", "entropyProd")
  
  return(mydf)
}
