#-------------------------------------------------------------------------------#
#                Moist Energy Balance Model (MEBM) -- Initialize                # 
#                                                                               #
# Following the solution of Gerard Roe, Univ. Washington (see Roe et al., 2015) #
# Hydroclimate parameterizations after Held, 2001; from Siler et al., (2018)    #
#      Coded and modified in R by Tyler Kukla, Stanford Univ. June, 2018        #
#-------------------------------------------------------------------------------#
# req'd for hydro solve
library(tibble)  # use of tibble is since abandoned (09/2022)
library(pracma)
# req'd for main ODE 
library(bvpSolve)
# req'd for plotting
library(ggplot2)
library(ggpubr)


# ------------------------------------------------------------------------ #
#                     USED FOR SENSITIVITY ANALYSIS                        #
# ------------------------------------------------------------------------ #

# NOTE: get_pars = TRUE means the train does not run MEBM but returns the parameters
#                  in the function instead 
MEBM_TRAIN.sensitivity <- function(pCO2_in, custom_insol = FALSE, insol_dat = NULL, get_pars = FALSE, inputPars = parList,
                                   LandFrac_tx = parList$MEBMpars$LandFrac_tx, ice_Threshold = parList$TRAINpars$ice_Threshold, 
                                   landAlf = parList$TRAINpars$landAlf, oceanAlf = parList$TRAINpars$oceanAlf, 
                                   iceAlf = parList$TRAINpars$iceAlf, rel_hum = parList$TRAINpars$rel_hum, 
                                   budy_om = parList$MEBMpars$budy_om, modCO2 = parList$MEBMpars$modCO2, 
                                   myD = parList$MEBMpars$myD, Q0 = parList$MEBMpars$Q0, 
                                   A0.intercept = parList$MEBMpars$A0.intercept, A0.slope = parList$MEBMpars$A0.slope,
                                   B0 = parList$MEBMpars$B0, alf_flag = parList$MEBMpars$alf_flag){
  
  # -------------------------------------------------------------------------------------------------- #
  #                                     CONSTRUCT THE MEBM DOMAIN                                      #
  # -------------------------------------------------------------------------------------------------- #
  # ... Initialize MEBM
  # Domain grid and boundary conditions
  nNodes <- inputPars$MEBMpars$nNodes                     # number of nodes in the MEBM grid
  xin <<- inputPars$TRAINpars$xin     # x-values (0-1 is equator-pole); (-1 to 1 version [pole to pole] is not yet working)
  LBC <- RBC <- 0                   # zero-flux boundary condition
  # temp_guess <- 15                  # [degC] guess for initial temp value (no need to guess flux because we set BCs)
  # [temp_guess is currently fixed in MEBM_solve.R]
  
  #... Albedo and humidity parameterizations (sorry for sloppy coding here...)
  ice_Threshold <<- ice_Threshold               # [degC] threshold for ice formation
  landAlf <<- landAlf                    # [] average albedo of ice-free land surface
  oceanAlf <<- oceanAlf                # [] average albedo of ocean
  iceAlf <<- iceAlf                     # [] average albedo of ice 
  rel_hum <<- rel_hum                    # [] relative humidity
  # budyko omega
  budy_om <- budy_om
  
  #... Physical parameters input to MEBM
  # Diffusion coefficient
  D_HF10 <- 1.06e6                  # [m2 s-1] diffusivity in Hwang and Frierson 2010
  D_Sil <- 1.16e6                   # [m2 s-1] diffusivity in Siler et al., 2018
  Din <<- myD * p0/g                 # Here, F = -2pi*Re*cos(phi)*D*dh/dy. For DHF10 F = -2pi*Re*cos(phi)*D*(ps/g)*dh/dy    
  # Radiation parameterization
  Q0 <- Q0
  B0 <- B0
  #A0 <- 207.3                      # [W m-2] OLR constant - NO CO2 DEPENDENCY
  A0 <- A0.intercept-A0.slope*log(pCO2_in/modCO2) # 207.3-5.35*log(pCO2_in/420) # [W m-2] OLR constant depending on pCO2 [Myhre et al., 1998, North and Kim, 2017]
  # OPTIONS: "step" = discrete step from no ice to ice albedo
  #          "D18" = hyperbolic (smoothed) transition from no ice to ice albedo
  #                  from Dortmans et al., 2018 (Climate of the Past) -- this caused stability issues in the past
  #          "prescribed" = use calculated albedo from landfrac data
  # -------------------------------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------------------------------- #
  #                                   MODEL INPUT CALCULATIONS                                         #
  # -------------------------------------------------------------------------------------------------- #
  #... Continentality and calculate ice free albedo
  landFrac <- LandFrac_tx$LandFrac_inLatBelt_forAlbedo
  lat_orig <- LandFrac_tx$Lat
  lat_land_sinx <- sin(lat_orig * (pi/180))
  landFrac_int <- approx(x=lat_land_sinx, y=landFrac, xout=xin)[[2]]    # the interpolated fraction of land per lat
  
  alf_bylat <- (landFrac_int * landAlf) + ((1-landFrac_int) * oceanAlf)  # [] latitudinal albedo weighted by continent
  alf_noice <- (landFrac_int * landAlf) + ((1-landFrac_int) * oceanAlf)  # [] latitudinal albedo weighted by continent
  afun_df <- as.data.frame(list(x = xin, alf_noice = alf_noice))
  afun <<- approxfun(afun_df, rule=2)
  
  # Insolation flag
  if(custom_insol==FALSE){
    insol_flag <<- "Traditional"
  } else if(custom_insol==TRUE){
    insol_flag <<- "INSOL"}
  # OPTIONS: "INSOL" = get insolation from a given latitude / insol file (already loaded)
  #          "Traditional" = get insolation the old fashioned way (using Q0)
  
  
  #... Create a function returning insolation given some value of "x"
  if(insol_flag=="INSOL"){
    # isolate just this year
    yr_df <- insol_dat
    
    #... translate latitude to the sinx grid
    # first get these lats in the correct format
    yr_df$lat_x <- sin(yr_df$lat * (pi/180))
    # interpolate to new grid
    myInsol <- approx(x=yr_df$lat_x, y=yr_df$insol, xout=xin, rule=2)
    
    # create a function to call
    InsolFun <<- approxfun(myInsol, rule=2)
  }
  # -------------------------------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------------------------------- #
  
  
  
  #------------------------------------------GENERATE AND PLOT MODEL OUTPUT------------------------------------------#
  if(get_pars == FALSE){   # then run the normal MEBM; otherwise don't run mebm and just return parameters
    # output dataframe
    pars_in <- c('A0' = A0, 'B0' = B0, 'Q0' = Q0, 'Din' = Din, 'alf_flag' = alf_flag)
    df <- MEBM_solve(xin_mod = xin, LBC_in = LBC, RBC_in = RBC, budy_om_in = budy_om, mypars = pars_in)
    
    # add continentality to output dataframe 
    # -- Fraction of total land area
    tempFrac <- stats::approx(x=LandFrac_tx$Lat, y=LandFrac_tx$Frac_of_tot_LandArea, xout=df$Lat, rule=2)$y   # interpolated fraction yields sum > 1 if MEBM resolution > paleogeog resolution
    df$Frac_of_tot_LandArea <- tempFrac / sum(tempFrac)     # rescale so sums up to 100%
    # -- Land area [m2]
    tempArea <- stats::approx(x=LandFrac_tx$Lat, y=LandFrac_tx$Land_Area_m2, xout=df$Lat, rule=2)$y   # interpolated area yields too much area if MEBM resolution > paleogeog resolution
    df$Land_Area_m2 <- tempArea / sum(tempFrac)   # tempFrac scales MEBM resolution to paleogeography resolution
    
    return(df)
  } else if(get_pars == TRUE){  # return just a list of the variables
    
    return(sapply(ls(),function(x)get(x),simplify=F,USE.NAMES=T))
    
  }
  
}
