#-------------------------------------------------------------------------------#
#              Moist Energy Balance Model (MEBM) -- Hydro-solve                 # 
#                                                                               #
# Following the solution of Gerard Roe, Univ. Washington (see Roe et al., 2015) #
#                          ------------------------                             #
#              Formulation for the Hadley cell based on Held, 2001              #
#                 Weighting function from Siler et al., 2018                    #
#                          ------------------------                             #
#       Coded and modified in R by Tyler Kukla, Stanford Univ. June, 2018       #
#-------------------------------------------------------------------------------#

# function to convert from [W m-2] to [m yr-1]
Wm2flux <- function(Wm2){
  #... terms for conversion
  rho <- 1e3    # [kg m-3] density of water
  Lv <- 2.5e6   # [J kg-1] latent heat of vaporization
  spyr <- 60*60*24*365   # [s yr-1] seconds per year
  mm_m <- 1e3   # [mm m-1] millimeters per meter (multiply the result by this for mm/yr output)
  #... conversion
  ms <- Wm2*(1/rho)*(1/Lv)   # conversion to meters per second
  ms*spyr   # meters per year 
}

# function to solve hydrologic cycle from MEBM mse output
hydrofun <- function(temp, x, hadley_extent=15){
  # moist static energy
  qstar <- (eps/p0)*e0*exp((a*temp)/(b+temp))           # [g kg-1] saturation specific humidity 
  q <- rel_hum*qstar                                    # [g kg-1] specific humidity 
  h <- cp*(temp+273.15)+Lv*q                            # [J kg-1] Moist static energy (MSE)
  flux <- -2*pi*Din*(1-xin^2)*pracma::gradient(h, xin)       # [W] poleward MSE flux
  
  # hadley cell weighting function
  idx <- which(x>=0)                       # identify x values greater than or equal to zero... (seems like a formality)
  x15 <- sin(deg2rad(hadley_extent))         # find the x value equal to the latitude of prescribed hadley cell extent
  xw <- 0.15                                 # sets smoothness of transition from hadley to eddy transport
  wt <- 1-exp(-x^2/0.3^2) #0.5*(1+tanh((xin[idx]-x15)/xw))      # weighting function for eddy fluxes (from Nick Siler)
  # wt <- cos(pi*xin)*(1-abs(xin))
  
  # hadley cell partitioning
  heq <- max(h)                                # [J kg-1] moist static energy at the equator
  gms_factor <- 1.06                         # 1 + percentage increase in GMS per change in heq -- 1.06 in Siler et al. 2018; 1.08 in Bonan et al. 2023 
  ht <- gms_factor * heq                           # [J kg-1] mse in tropical upper troposphere (fit by Siler et al. 2018)
  gms <- ht - h                              # [J kg-1] gross moist stability
  F_hc <- (1-wt)*flux                        # Hadley cell flux
  V <- F_hc/gms                   # Mass transport of the hadley cell 
  #F_hc_approx = approxfun(F_hc,xin)
  #F_hc_approx(0)
  # which(V[500:1200]==min(abs(V[500:1200]-0)))
  # V[which(min(abs(V-0))]
  
  V[1] <- 0                                  # restore 0 flux BC
  V[length(V)] <-0
  
  F_LH_eddy <- wt*-2*pi*Din*(1-x^2)*pracma::gradient(Lv*q, x)   # [W] eddy latent heat fluxes including weighting function
  F_LH <- -Lv*V*q + F_LH_eddy                # [W] total latent heat flux, including equatorward component due to Hadley Cell
  F_LH_hadley <- -Lv*V*q                     # [W] latent heat flux from hadley cell transport
  F_DS <- flux - F_LH                        # [W] dry static heat flux
  divF_LH = 1/(2*pi*Re^2)*pracma::gradient(F_LH,x) # latent heat flux divergence
  divF_LH[1:2] <- divF_LH[3]
  divF_LH[(length(xin)-1):length(xin)] <- divF_LH[length(xin)-2]
  E_m_P = divF_LH/(Lv*rho)*pi*1e7;           # [m yr-1] E minus P 
  
  # Evaporation / Precipitation partitioning
  alpha <- Lv/461.5/(temp+273.15)^2                 # alpha parameter (from Nick Siler.. original ref?)
  beta <- cp/Lv/alpha/qstar                             # beta parameter
  RG <- 180*(1*(1-x^2)-0.4*exp(-(x/0.15)^2))    # [W m-2] made-up idealized pattern of R-G
  #%RG=170*(1*(1-x.^2));%-.4*exp(-(x./.15).^2))
  Ch <- 1.5e-3                                      # drag coefficient
  LWfb <- 0                                         # [W m-2 K-1] LW feedback at surface
  u <- 4+abs(sin(pi*x/1.5))*4                     # [m s-1] wind speed   [ #%-2.5*cos(3*asin(x)) ]
  rho_air <- 1.2                                    # air density   [ psfc./287./(T_ctrl+273.15) ]
  Evap <- (RG*alpha+rho_air*cp*(1-rel_hum)*Ch*u)/(alpha+cp/Lv/qstar)  # [W m-2] evaporation
  Prec <- Evap-divF_LH                                            # [W m-2] precipitation
  Evap <- Evap/(Lv*rho)*pi*1e7                      # [m yr-1] evaporation
  Prec <- Prec/(Lv*rho)*pi*1e7                      # [m yr-1] precipitation
  
  
  # generate an output dataframe 
  myout <- as.data.table(cbind(temp, h, Prec, Evap, E_m_P, flux, F_LH, F_DS, F_LH_eddy, F_LH_hadley, divF_LH))
  colnames(myout) <- c("temp_C", "MSE_Jkg", "Precip_m_yr", "Evap_m_yr", "EminusP_m_yr", "mseFlux_W", "latentFlux_W", "sensibleFlux_W",
                       "latentEddyFlux_W", "latentHadleyFlux_W", "LHflux_Div")

  return(myout)
}

