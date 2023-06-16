# Functions for CH2O-CHO
# For coupling with MEBM

# CH2OCHO coupling and train.tracks functions by T Kukla
# For more detailed explanation, please review Caves Rugenstein et al. 2019 (Nature). 

# Jeremy Rugenstein acknowledges Kimberly Lau and Adam Jost for help in 
# using package('deSolve').

# Created on: 04-03-2020
# Last Modified: 09-07-2022
# Last Modified by: TK

# Branch off of CLiBeSO-water-func.R

## --- function for trying to solve the conflicted file issue with saveRDS --- ## 
## ... https://stackoverflow.com/questions/44456298/issue-with-saverds 
## NOT USED
# saveRDS2 <- function(object,file){str(object);saveRDS(object,file)}

#### Model Parameter Function ####
paramfunc <- function(remove=NA) {
  p <- c( #parameter vector
    
    # Physical Constants
    Clim.Sens = Clim.Sens, # Earth System Sensitivity (°C/CO2 doubling) 
    temp.a.i = temp.a.i, # Temperature of the atmosphere
    q.sens = 0.04, # Sensitivity of runoff to warming
    
    # Initial Sulfur fluxes
    Fwsulf.i = 0.5e12, # from Calmels et al. 2007 (Geology)
    Fbsulf.i = 0.5e12, # To balance the initial sulfur cycle
    
    # Initial Carbon and Alkalinity fluxes (mol/yr)
    Fvolc.i = Fvolc.i, 
    Fbcarb.i = Fbcarb.i, 
    Fwsil.i = Fvolc.i, 
    Fwcarb.i = Fwcarb.i, 
    Fworg.i = Fworg.i, 
    Fborg.i = Fborg.i, 
    
    # Initial d13C Values
    Fwcarb.d13C = 0,
    Fworg.d13C = Fworg.d13C,
    Fvolc.d13C = -5, 
    DB = 27, # Cap-delta between inorganic carbon and organic carbon
    
    # Phosphorus Parameters (default parameters from Shield and Mills 2017)
    # 1 = no PO4 feedback
    # 2 = Van Cappellen and Ingall 1994 parameterization
    # 3 = Original COPSE formulation (Bergman et al., 2004)
    # 4 = Shields and Mills 2017
    # 5 = Payne and Kump 2007
    # 6 = arbitrary temp-dependent organic carbon burial feedback 
    #     [6 requires arbFB.n, positive values = negative feedback w/temp and vice versa]
    PO4.form = 1, # Determines what formulation of phosphorus burial/organic carbon burial to use
    Fwp.i = 2.4e10, # Assumes initial burial flux is 250x less than Fborg (initial value in CLiBeSO: 4.7e10)
    Fbp.i = 2.4e10, # Assumes initial burial flux is 250x less than Fborg (initial value in CLiBeSO: 4.7e10)
    sil.P = 0.58, # P derived from silicate weathering
    carb.P = 0.21, # P derived from carbonate weathering
    org.P = 0.21, # P derived from organic carbon weathering
    PO4.n = 3, # Exponential describing dependence of Fborg on [PO4] (range from 2-3) (Van Cappellen and Ingall 1994)
    arbFB.n = 0.2,  # (used w/ option 6) Exponent for Fborg feedback w/ temperature (positive is more Fborg w/ warming, negative is less)
    
    # Parameters used in MC 2014 Function
    keff.ref = 8.7e-6, # mol/m2/yr (reference reaction rate) (from Maher and Chamberlain SI Table S1)
    m = 270, # g/mol (molar mass) (from Maher and Chamberlain SI Table S1)
    A = 0.1, # m2/g (specific surface area) (from Maher and Chamberlain SI Table S1)
    Rn.max.ref = 1085, # umol SiO2/L/yr # Maher and Chamberlain 2014
    L.phi = 0.1, # Reactive length scale times effective porosity (Maher and Chamberlain 2014)
    Ts = 2e4, # Timescale for Solids (Maher and Chamberlain 2014)
    Ceq = 374, # Equilibrium Concentration for SiO2 (Maher and Chamberlain 2014)
    soil.CO2.flag = 1, # Flag to activate soil CO2 module
    # 1: Soil CO2 calculated using Volk 1987
    # 2: Soil CO2 equals atmospheric CO2
    # 3: Soil CO2 equals atmospheric CO2 multiplied by magnification factor (soil.CO2.mag)
    soil.CO2.mag = 10, # Soil CO2 magnification factor
    Ea = 38, # (kJ/mol) Activation Energy (Maher and Chamberlain 2014)
    
    # Other scaling exponents and factors
    borg.carb = 1, # Scaling between organic carbon burial and rain rate
    carb.clim.exp = 0.05, # Exponent carbonate weathering Arrhenius reaction (not implemented as a changeable parameter)
    #Rk.i = R.k.i, #1/(log2(pCO2.i/280)+1), # Silicate reactivity of land surface (k in Eq. 3) (Caves et al., 2016--EPSL)
    Trop.area = 1/3, # Effective "area" for tropical weathering
    Conc.sens = 0.01, # (% change per K)--only used in CLiBeSO.qC
    q0 = 0.3, # Initial runoff
    silicate.lith = 1, #1 = all granite, 0 = all basalt
    
    #Start.Ma = 70, # Starting time of the model (Ma)
    nbox = 200 # Number of latitudinal boxes
  )
  
  if (is.na(remove)==TRUE) {return.p <- p}
  else {for (i in 1:length(remove)){
    p <- p[-(which(names(p)==c(remove[i])))]}
    return.p <- p
  }
  return.p
}

#### Carbonate Solvers ####
# From Adam Jost

# CO2 Sys--pH and pCO2
# Used only to set initial conditions
co2sys.init <- function(pH, pCO2, temp=temp.o.i, sal=salinity, p=300, Ca2=Ca.i, Mg=Mg.i, SO42=0.0282){
  
  # Inputs should be in mol/kg
  
  Mg.Ca <- Mg/Ca2
  
  # From Wikipedia
  Ca2.m <- 0.01028
  Mg.m <- 0.05654
  SO42.m <- 0.0282
  
  # From Zeebe and Tyrrell (2019; GCA)--Table 2
  Ca.K1 <- 5e-3
  Mg.K1 <- 17e-3
  SO42.K1 <- 208e-3
  
  Ca.K2 <- 157e-3
  Mg.K2 <- 420e-3
  SO42.K2 <- 176e-3
  
  Ca.Ksp <- 185e-3
  Mg.Ksp <- 518e-3
  SO42.Ksp <- 106e-3
  
  t=273.15+temp #convert to K
  k0=0.02839188 #from seacarb, assumes sal=35, temp=25, p=0
  
  #calculate k1, k2, kw, and kb based on the method by DOE 1994 (as shown in Zeebe and Wolf-Gladrow)
  k1.m <- exp(2.83655 - 2307.1266/t - 1.5529413*log(t) - (0.207608410 + 4.0484/t)*sqrt(sal) + 0.0846834*sal - 0.00654208*sal^(3/2) + log(1-0.001005*sal))
  k1 <- k1.m*(1+((Ca.K1*(Ca2/Ca2.m-1))+(Mg.K1*(Mg/Mg.m-1))+(SO42.K1*(SO42/SO42.m-1))))
  
  k2.m = exp(-9.226508 - 3351.6106/t - 0.2005743*log(t) - (0.106901773 + 23.9722/t)*sqrt(sal) + 0.1130822*sal - 0.00846934*sal^(3/2) + log(1 - 0.001005*sal))
  k2 <- k2.m*(1+((Ca.K2*(Ca2/Ca2.m-1))+(Mg.K2*(Mg/Mg.m-1))+(SO42.K2*(SO42/SO42.m-1))))
  
  kw = exp(148.96502 - 13847.26/t - 23.6521*log(t) + (118.67/t - 5.977 + 1.0495*log(t))*sal^(1/2) - 0.01615*sal)
  kb = exp((-8966.90 - 2890.53*sal^(1/2) - 77.942*sal + 1.728*sal^(3/2) - 0.0996*sal^2)/t + 148.0248 +137.1942*sal^(1/2) + 1.62142*sal - (24.4344 + 25.085*sal^(1/2) + 0.2474*sal)*log(t) + 0.053105*sal^(1/2)*T)
  
  bt = 4.16e-4*sal/35 #eq. A.7.15 in Zeebe and Wolf-Gladrow 2001
  
  #calculate ksp for calcite and aragonite from Muci (1983) as described in Zeebe and Wolf-Gladrow, 2001
  ksp.c.m = 10^(-171.9065 - 0.077993*t + 2839.319/t + 71.595*log10(t) + (-0.77712 + 0.0028426*t + 178.34/t)*sal^(1/2) - 0.07711*sal + 0.0041249*sal^1.5)
  ksp.c <- ksp.c.m*(1+((Ca.Ksp*(Ca2/Ca2.m-1))+(Mg.Ksp*(Mg/Mg.m-1))+(SO42.Ksp*(SO42/SO42.m-1))))
  
  h <- 10^(-pH) # convert to [H+]
  CO2 <- pCO2*k0/1e6
  DIC <- CO2*(1+k1/h+k1*k2/h^2)
  HCO3 <- DIC/(1+h/k1+k2/h)
  CO32 <- DIC/(1+h/k2+h^2/k1*k2)
  OmegaCalcite = Ca2*CO32/ksp.c
  TA <- CO2*(k1/h+2*k1*k2/h^2)+bt*kb/(kb+h)+kw/h-h
  
  # list <- c(DIC=DIC,ALK=TA,CO2=CO2,HCO3=HCO3,CO32=CO32,OmegaCalcite=OmegaCalcite)
  results <- as.data.frame(array(data=c(DIC,TA,CO2,HCO3,CO32,OmegaCalcite),dim=c(1,6)))
  colnames(results) <- c("DIC","ALK","CO2","HCO3","CO32","OmegaCalcite")
  results
}

# CO2-Sys--accepts DIC and TA
# Used within the CH2O-CHO functions
co2sys <- function(DIC, TA, temp=temp.o.i, sal=salinity, p=300, Ca.in=Ca.i, Mg.in=Mg.i, SO42.in=0.0282){
  
  # Inputs should be in mol/kg
  
  Mg.Ca <- Mg.in/Ca.in
  
  # From Wikipedia
  Ca.m <- 0.01028
  Mg.m <- 0.05654
  SO42.m <- 0.0282
  
  # From Zeebe and Tyrrell (2019; GCA)--Table 2
  Ca.K1 <- 5e-3
  Mg.K1 <- 17e-3
  SO42.K1 <- 208e-3
  
  Ca.K2 <- 157e-3
  Mg.K2 <- 420e-3
  SO42.K2 <- 176e-3
  
  Ca.Ksp <- 185e-3
  Mg.Ksp <- 518e-3
  SO42.Ksp <- 106e-3
  
  t=273.15+temp #convert to K
  k0=0.02839188 #from seacarb, assumes sal=35, temp=25, p=0
  
  #calculate k1, k2, kw, and kb based on the method by DOE 1994 (as shown in Zeebe and Wolf-Gladrow)
  k1.m <- exp(2.83655 - 2307.1266/t - 1.5529413*log(t) - (0.207608410 + 4.0484/t)*sqrt(sal) + 0.0846834*sal - 0.00654208*sal^(3/2) + log(1-0.001005*sal))
  k1 <- k1.m*(1+((Ca.K1*(Ca.in/Ca.m-1))+(Mg.K1*(Mg.in/Mg.m-1))+(SO42.K1*(SO42.in/SO42.m-1))))
  
  k2.m = exp(-9.226508 - 3351.6106/t - 0.2005743*log(t) - (0.106901773 + 23.9722/t)*sqrt(sal) + 0.1130822*sal - 0.00846934*sal^(3/2) + log(1 - 0.001005*sal))
  k2 <- k2.m*(1+((Ca.K2*(Ca.in/Ca.m-1))+(Mg.K2*(Mg.in/Mg.m-1))+(SO42.K2*(SO42.in/SO42.m-1))))
  
  kw = exp(148.96502 - 13847.26/t - 23.6521*log(t) + (118.67/t - 5.977 + 1.0495*log(t))*sal^(1/2) - 0.01615*sal)
  kb = exp((-8966.90 - 2890.53*sal^(1/2) - 77.942*sal + 1.728*sal^(3/2) - 0.0996*sal^2)/t + 148.0248 +137.1942*sal^(1/2) + 1.62142*sal - (24.4344 + 25.085*sal^(1/2) + 0.2474*sal)*log(t) + 0.053105*sal^(1/2)*T)
  
  bt = 4.16e-4*sal/35 #eq. A.7.15 in Zeebe and Wolf-Gladrow 2001
  
  #calculate ksp for calcite and aragonite from Muci (1983) as described in Zeebe and Wolf-Gladrow, 2001
  ksp.c.m = 10^(-171.9065 - 0.077993*t + 2839.319/t + 71.595*log10(t) + (-0.77712 + 0.0028426*t + 178.34/t)*sal^(1/2) - 0.07711*sal + 0.0041249*sal^1.5)
  ksp.c <- ksp.c.m*(1+((Ca.Ksp*(Ca.in/Ca.m-1))+(Mg.Ksp*(Mg.in/Mg.m-1))+(SO42.Ksp*(SO42.in/SO42.m-1))))
  
  #creating a function to solve for [H+] (referred to as 'h' here) given DIC and TA. Refer to Zeebe and Wolf-Gladrow 2001 for instructions
  h.solve = function(DIC, TA){
    #polynomials
    a = -1
    b = -TA - kb - k1
    c = DIC*k1 - TA*kb - TA*k1 + kb*bt + kw - kb*k1 - k1*k2
    d = DIC*kb*k1 + 2*DIC*k1*k2 - TA*kb*k1 - TA*k1*k2 + kb*bt*k1 + kw*kb + kw*k1 - kb*k1*k2 
    e = 2*DIC*kb*k1*k2 - TA*kb*k1*k2 + kb*bt*k1*k2 + kw*kb*k1 + kw*k1*k2
    f = kw*kb*k1*k2
    p = c(f,e,d,c,b,a)
    r = polyroot(p)
    h=max(Re(r))
    return(h)
  }
  
  h = h.solve(DIC, TA)
  pH = -log10(h)
  CO2 = DIC/(1 + k1/h + k1*k2/(h^2))
  HCO3 = DIC/(1 + h/k1 + k2/h)
  CO32 = DIC/(1 + h/k2 + (h^2)/k1*k2)
  OmegaCalcite = Ca.in*CO32/ksp.c
  #omega.a = Ca2*CO32/ksp.a
  pCO2=CO2/k0*1e6
  
  results <- as.data.frame(array(data=c(DIC,TA,CO2,HCO3,CO32,OmegaCalcite,pH,pCO2),dim=c(1,8)))
  colnames(results) <- c("DIC","ALK","CO2","HCO3","CO32","OmegaCalcite","pH","pCO2")
  results
}

#### Weathering Functions ####
soilCO2func <- function(pCO2, pCO2.init, prod.max=2,pCO2.min=100) {
  
  # From Tyler Volk (1989)
  
  prod0 <- 1
  pCO2.soil.i <- pCO2.init*10
  pCO2.50 <- (prod.max/prod0-1)*(pCO2.init-pCO2.min)
  prod <- prod.max*(pCO2-pCO2.min)/(pCO2.50+(pCO2-pCO2.min))
  pCO2.soil <- (prod/prod0*(1-pCO2.init/pCO2.soil.i)+pCO2/pCO2.soil.i)*pCO2.soil.i # +(pCO2-pCO2.init)
  RCO2.soil <- pCO2.soil/pCO2.soil.i
  
  RCO2.soil
}

# Arrhenius Function 
arrhenius <- function(Temp.C,Ea.in){
  R <- 8.314 # J/K/mol
  Ea.used <- Ea.in*1000
  
  mkeffa <- exp((Ea.used/R)*((1/(p['temp.a.i']+273.15))-(1/(Temp.C+273.15))))
  mkeffa
}


# Maher and Chamberlain (2014) equation
MCfunc.trans <- function(pCO2,q0,q,Ea=p['Ea'],Ts=p['Ts'],keff.ref=p['keff.ref'],
                         Temperature=Temp.bin,soilCO2.use=p['soil.CO2.flag'],
                         mag.factor=p['soil.CO2.mag'], Temp.init=Temp.bin.i,
                         pCO2.modern = 280){
  
  # Ea is input in kJ/mol
  
  RCO2 <- pCO2/pCO2.i
  
  if (soilCO2.use==1){
    RCO2.soil <- soilCO2func(pCO2=pCO2, pCO2.init = pCO2.i)
    RCO2.init.mod.soil <- soilCO2func(pCO2=pCO2.i, pCO2.init = pCO2.modern)
  }
  if (soilCO2.use==2){
    RCO2.soil <- RCO2
    RCO2.init.mod.soil <- pCO2.i/pCO2.modern
  }
  if (soilCO2.use==3){
    RCO2.soil <- RCO2*mag.factor
    RCO2.init.mod.soil <- (pCO2.i/pCO2.modern) * mag.factor
  }
  
  m.init <- p['m'] # g/mol (molar mass) (from Maher and Chamberlain SI Table S1)
  A.init <- p['A'] # m2/g (specific surface area) (from Maher and Chamberlain SI Table S1)
  keff.init <- keff.ref*arrhenius(Temp.C=Temp.init,Ea.in=Ea)
  Rn.max.init <- p['Rn.max.ref']*(keff.init/keff.ref)
  fw.init <- 1/(1+m.init*A.init*keff.init*Ts)
  Ceq.init <- ((RCO2.init.mod.soil)^(0.316))*p['Ceq'] #1000
  
  # timeslice keff and rn.max
  keff <- keff.ref*arrhenius(Temp.C=Temperature,Ea.in=Ea)
  Rn.max <- p['Rn.max.ref']*(keff/keff.ref)
  fw <- 1/(1+p['m']*p['A']*keff*Ts)
  
  # Calculation Initial Concentration
  # changed Dw0 to be scaled to the initial temperature
  Dw0 <- p['L.phi']*Rn.max.init*fw.init/Ceq.init
  Conc0 <- p['Ceq']*(((exp(1)^2)*Dw0/q0)/(1+(exp(1)^2)*Dw0/q0))
  
  Ceq <- ((RCO2.soil)^(0.316))*p['Ceq'] #1000
  Dw <- p['L.phi']*Rn.max*fw/Ceq
  Conc <- Ceq*(((exp(1)^2)*Dw/q)/(1+(exp(1)^2)*Dw/q))
  
  # -- apply scaling for carbonate and basalt
  # Carbonate weathering
  carb.dw.scale <- 2.5    # following unpublished analysis of Ibarra and Winnick (along with power-law scaling of Bluth and Kump 1994)
  carb.Ceq.scale <- 2
  Dw0.carb <- Dw0 * carb.dw.scale
  Conc0.carb <- (p['Ceq']*carb.Ceq.scale) * (((exp(1)^2)*Dw0.carb/q0)/(1+(exp(1)^2)*Dw0.carb/q0))
  Ceq.carb <- Ceq*carb.Ceq.scale
  Dw.carb <- Dw * carb.dw.scale
  Conc.carb <- Ceq.carb*(((exp(1)^2)*Dw.carb/q)/(1+(exp(1)^2)*Dw.carb/q))
  
  # Basalt weathering
  ba.dw.scale <- 3.45    
  ba.Ceq.scale <- 1.3
  Dw0.ba <- Dw0 * ba.dw.scale
  Conc0.ba <- (p['Ceq']*ba.Ceq.scale) * (((exp(1)^2)*Dw0.ba/q0)/(1+(exp(1)^2)*Dw0.ba/q0))
  Ceq.ba <- Ceq*ba.Ceq.scale
  Dw.ba <- Dw * ba.dw.scale
  Conc.ba <- Ceq.ba*(((exp(1)^2)*Dw.ba/q)/(1+(exp(1)^2)*Dw.ba/q))
  
  
  results <- c(q0,q,Conc0,Conc,Ceq,Dw,RCO2.soil)
  results.carb <- c(q0,q,Conc0.carb,Conc.carb,Ceq.carb,Dw.carb,RCO2.soil)
  results.ba <- c(q0,q,Conc0.ba,Conc.ba,Ceq.ba,Dw.ba,RCO2.soil)
  names(results) <- c("q0","q","Conc0","Conc","Ceq","Dw","RCO2.soil")
  names(results.carb) <- c("q0","q","Conc0","Conc","Ceq","Dw","RCO2.soil")
  names(results.ba) <- c("q0","q","Conc0","Conc","Ceq","Dw","RCO2.soil")
  results <- as.data.frame(t(as.matrix(results,nrow=1,ncol=length(results))))
  results.carb <- as.data.frame(t(as.matrix(results.carb,nrow=1,ncol=length(results))))
  results.ba <- as.data.frame(t(as.matrix(results.ba,nrow=1,ncol=length(results))))
  total.results <- list(results,results.carb,results.ba)
  names(total.results) <- c("Granite","Carbonate","Basalt")
  
  # Return
  total.results
}

#### Base CH2O-CHO Functions ####
CH2OCHO = function(t,y,p){ 
  # run the loop
  with(as.list(c(p,y)),{
    #... identify whether this is the first loop or second loop of the numerical solver
    if(t==0){loop.idx <<- loop.idx + 1}
    if(loop.idx < 3){   # this counts the number of times t==t[1] occurs (twice in first loop, once in second)
      this.loop <<- 1
    } else if(loop.idx >= 3){
      this.loop <<- 2
    }
    
    # define forcing for this timestep
    Fextra = force.C(t)
    d13C.force = force.d13C(t)
    
    # Calculate carbon system parameters (two-step iterative procedure to get effect of temperature change correct)
    # CURRENTLY TURNED OFF AS TEMP CHANGES ARE SMALL
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=temp.o.i,
                         sal=salinity,p=300,Ca.in=y['Ca'],Mg.in=y['Mg'])
    pCO2 <- round(carb.parms$pCO2, CO2_res)
    omega <- carb.parms$OmegaCalcite
    RCO2 <- pCO2/pCO2.i
    
    # Volcanism
    Fvolc <- p['Fvolc.i'] + Fextra
    
    
    # ************************************************************************** #
    #                 MOIST ENERGY BALANCE MODEL CLIMATOLOGY                     #
    #   -- And track boundary condition temperatures from one loop to next --    #
    # ************************************************************************** #
    if(this.loop==1){
      # if at the start, restore initial conditions
      if(t==0){LBtemp <<- LBtemp.initial; RBtemp <<- RBtemp.initial; pCO2.mebm <<- pCO2.i; names(pCO2.mebm) <- NULL
        if(Bistable.Run == "Y"){TempBC_fun(thisClim = clim.state)}  # only update temperature BCs if bistable is on
      } else{pCO2.mebm <<- pCO2 ; names(pCO2.mebm) <- NULL}
    # run MEBM
    df_MEBM <- try_with_time_limit(expr = MEBM_TRAIN.sensitivity(pCO2_in = round(pCO2.mebm, CO2_res), LandFrac_tx = thisPaleogeo, get_pars=FALSE ),
                                   cpu = cpu.timeMax, elapsed = elapsed.timeMax)
    #... look for something else if solution is NULL or snowball
    if(SnowballFinder=="on"){
      df_MEBM <- Snowball_NULL_navigator(thisdf = df_MEBM, thisCO2 = pCO2, SearchDir.BC = "R")
    }
    
    #... get the new climate state
    clim.state <- clim_state_ID(thisdf = df_MEBM) 
    
    if(ClimStateTester=="on"){
      #... if we're in a new climate state, test if it's a robust change (and not just an outlier)
      df_MEBM <- NewClimState_Navigator(thisdf = df_MEBM, thisCO2 = pCO2, iceDeal = clim.state, 
                                        iceDeal.last = clim.state.last, this.nav.count.conditions = CLIM.NAV.CONDITIONS)
      #... update the clim state if needed
      clim.state <- clim_state_ID(thisdf = df_MEBM)
    }
    #... and set the old one 
    clim.state.last <<- clim.state
    
    #... save vectors --
    # [1] for the timestep 
    # [2] for the climate state if loop == 1
    # [3] for the LBtemp and RBtemp if loop == 1 
    if(t==0 & this.loop==1){clim.state.vec <<- clim.state
      timestep.vec <<- t
      pco2.vec <<- round(pCO2.mebm, CO2_res)
      LBtemp.vec <<- LBtemp
      RBtemp.vec <<- RBtemp
    } else if(this.loop==1){
      clim.state.vec[length(clim.state.vec)+1] <<- clim.state
      timestep.vec[length(timestep.vec)+1] <<- t
      pco2.vec[length(pco2.vec)+1] <<- round(pCO2.mebm, CO2_res)
      LBtemp.vec[length(LBtemp.vec)+1] <<- LBtemp
      RBtemp.vec[length(RBtemp.vec)+1] <<- RBtemp
    }
    # ... get the new temperature boundary conditions after we get a solution that works
    if(Bistable.Run == "Y"){TempBC_fun(thisClim = clim.state)}  # only update temperature BCs if bistable is on
    if(Bistable.Run == "N"){LBtemp <<- LBtemp.initial ; RBtemp <<- RBtemp.initial}
    
    
    } else{  # if second loop, assign values from first loop
      # assign the climate state based on the first loop
      clim.state.t <<- clim.state.vec[which(timestep.vec==t)][1]
      if(Bistable.Run == "Y"){TempBC_fun(thisClim = clim.state.t)}  # only update temperature BCs if bistable is on
      if(impose_loop1_BCs == "Y"){ # overwrite the boundary conditions to be consistent with the first loop (rather than just assigning based on the climate state)
        LBtemp <<- LBtemp.vec[which(timestep.vec==t)][1]
        RBtemp <<- RBtemp.vec[which(timestep.vec==t)][1]
      }
      # assign CO2 based on the first loop
      pCO2.mebm <<- pco2.vec[which(timestep.vec==t)][1] ; names(pCO2.mebm) <- NULL
      # override if t==1 (just to be sure we get the right initial value)
      if(t==0){if(Bistable.Run == "Y"){TempBC_fun(thisClim = clim.state)}  # only update temperature BCs if bistable is on
        if(impose_loop1_BCs == "Y"){ # overwrite the boundary conditions to be consistent with the first loop (rather than just assigning based on the climate state)
          LBtemp <<- LBtemp.vec[which(timestep.vec==t)][1]
          RBtemp <<- RBtemp.vec[which(timestep.vec==t)][1]
        }
        pCO2.mebm <<- pCO2.i}
      df_MEBM <- try_with_time_limit(expr = MEBM_TRAIN.sensitivity(pCO2_in = round(pCO2.mebm, CO2_res), LandFrac_tx = thisPaleogeo, get_pars=FALSE ),
                                     cpu = cpu.timeMax, elapsed = elapsed.timeMax)
      
      #... look for something else if solution is NULL or snowball
      if(SnowballFinder=="on"){
        df_MEBM <- Snowball_NULL_navigator(thisdf = df_MEBM, thisCO2 = pCO2, SearchDir.BC = "R")
      }
      #... get the new climate state
      clim.state <- clim_state_ID(thisdf = df_MEBM) 
      ## [NOTE: CLIMSTATETESTER is commented out because we assign based on the first loop here]
      # if(ClimStateTester=="on"){
      #   #... if we're in a new climate state, test if it's a robust change (and not just an outlier)
      #   df_MEBM <- NewClimState_Navigator(thisdf = df_MEBM, thisCO2 = pCO2, iceDeal = clim.state, iceDeal.last = clim.state.last)
      #   #... update the clim state if needed
      #   clim.state <- clim_state_ID(thisdf = df_MEBM)
      # }
    }
    
    # ... get integer code of the climate state
    clim.state.int <- clim_state_ID_INT(clim.state)
    # ... calculate global temp
    temp.i.wt <- round(weighted.mean(x = df_MEBM$temp_C, w=Lat_area_weight), temp_res)   # [degC] area-weighted mean global temperature
    runoff.i.wt <- round(weighted.mean(x=df_MEBM$budykoRunoff_m_yr, w=Lat_area_weight), q_res)
    runoff.i.wt.land <- round(weighted.mean(x=df_MEBM$budykoRunoff_m_yr, w=df_MEBM$Land_Area_m2), q_res)
    # ... assign new deep ocean temp
    temp.o.i <<- temp.i.wt - 10 ### assume a 10C gradient between air and deep ocean T
    # ... assign output 
    temp.lat <- round(df_MEBM$temp_C, temp_res)                 # [degC] MEBM temperature
    runoff.lat <-  round(df_MEBM$budykoRunoff_m_yr, q_res)   # [m yr-1] runoff calculated from Budyko framework (omega assigned in MEBM_TRAIN)
    # ... find iceline
    ice.latitudes <- which(df_MEBM$albedo >= iceAlf)
    if(length(ice.latitudes)== 0){iceline.lat <- 9999
    } else{iceline.lat <- min(abs(df_MEBM[ice.latitudes,]$Lat))}
    # ************************************************************************** #
    
    
    # ... set up vectors to fill in for loop
    Fwsil.vec <- numeric(length=p['nbox'])
    Fwcarb.vec <- numeric(length=p['nbox'])
    Conc.vec <- numeric(length=p['nbox'])
    Conc0.vec <- numeric(length=p['nbox'])
    Conc.carb <- numeric(length=p['nbox'])
    Conc0.carb <- numeric(length=p['nbox'])
    Dw.vec <- numeric(length=p['nbox'])
    area.frac.vec <- df_MEBM$Frac_of_tot_LandArea
    area.vec <- df_MEBM$Land_Area_m2
    # ------------------- DELETE?? -------------------- #
    # area.vec=rep(0.8,p['nbox']) #this is the land area parameter -- need to fill in from input
    # total.area=sum(area.vec)#sum(area.vec)*p['nbox']
    # ------------------------------------------------- #
    for (i in 1:p['nbox']){
      Temp.bin <- temp.lat[i] # temperature in latitude bin
      
      # Runoff    
      q0 <- runoff0[i]
      q <- runoff.lat[i]
      
      # Using a Maher and Chamberlain-style concentration feedback with a Winnick and Maher 
      # open-system calculation for Ceq
      Conc.results <- MCfunc.trans(pCO2 = pCO2, q0 = runoff.INIT[i], q = q, Temperature = Temp.bin,
                             Ea = p['Ea'], Ts=p['Ts'], soilCO2.use=p['soil.CO2.flag'],
                             Temp.init = Temp.bin.INIT[i])

      # Silicate Weathering
      Conc0.vec[i] <- Conc.results$Granite$Conc0*p['silicate.lith'] + Conc.results$Basalt$Conc0*(1-p['silicate.lith'])
      Conc.vec[i] <- Conc.results$Granite$Conc*p['silicate.lith'] + Conc.results$Basalt$Conc*(1-p['silicate.lith'])
      Dw.vec[i] <- Conc.results$Granite$Dw*p['silicate.lith'] + Conc.results$Basalt$Conc*(1-p['silicate.lith'])

      Fwsil.vec[i] <- q*Conc.vec[i]*area.vec[i]*Fwsil_scalar     # NOTE: q is not indexed by "i" because it's assigned at beginning of loop

      # Carbonate Weathering
      Conc0.carb[i] <- Conc.results$Carbonate$Conc0
      Conc.carb[i] <- Conc.results$Carbonate$Conc
      Fwcarb.vec[i] <- q*Conc.carb[i]*area.vec[i]*Fwcarb_scalar
    }
    # add weathering fluxes to MEBM output 
    df_MEBM$Fwsil <- Fwsil.vec
    df_MEBM$Fwcarb <- Fwcarb.vec
    # add discharge since we didn't earlier
    df_MEBM$Discharge_m3_yr <- df_MEBM$budykoRunoff_m_yr * df_MEBM$Land_Area_m2
    df_MEBM$time <- t
    df_MEBM$RUN.idx <- this.run
    # ... save MEBM results if end of second loop
    if(this.loop > 1){
      if(t == 0){  # if first timestep, create new save file
        save.MEBM <<- df_MEBM
      } else if(t %in% save.clim.here){   # if selected timestep, save this iter of clim
        save.MEBM <<- as.data.table(rbind(save.MEBM, df_MEBM))
      }
      # save result if last ts
      if(t == save.clim.here[length(save.clim.here)]){
        saveRDS(save.MEBM, file = paste(saveDir, paste(RunName_BASE, '-CLIM_TS-idx_', this.run, '.RDS', sep=''), sep='/'))
      }
    }
    
    # Sum weathering fluxes
    Fwsil <- sum(Fwsil.vec, na.rm=T)
    Fwcarb <- sum(Fwcarb.vec, na.rm=T)
    
    
    # Average global [C] (not weighted)
    Conc.avg <- mean(Conc.vec, na.rm=TRUE)
    
    Fbcarb <- p['Fbcarb.i']*(omega/omega.i) #carbonate burial, mol/yr
    
    # Alkalinity Fluxes    
    Falkin = 2*Fwcarb + 2*Fwsil
    Falkout = 2*Fbcarb
    
    # Organic weathering fluxes
    Fworg <- p['Fworg.i']
    
    # Phosphorus weathering fluxes
    Fwp <- p['Fwp.i']*(((p['sil.P'])*(Fwsil/p['Fwsil.i']))+((p['carb.P'])*(Fwcarb/p['Fwcarb.i']))+((p['org.P'])*(Fworg/p['Fworg.i'])))
    
    # Organic and Phosphorus Burial Fluxes
    if (p['PO4.form']==1) { # No coupling to phosphorus
      Fborg <- p['Fborg.i']*(Fbcarb/p['Fbcarb.i'])^p['borg.carb']
      Fbp <- Fwp
    }
    
    if (p['PO4.form']==2){ # Uses the Van Cappellen formulation
      Fborg <- p['Fborg.i']*(y['PO4']/PO4.i)^p['PO4.n']
      Fbp <- Fborg/250
    }
    
    if (p['PO4.form']==3){ # Uses the original COPSE formulation (Bergman et al., 2004)
      #NOTE: DYNAMIC OXYGEN IS NOT YET IMPLEMENTED IN THIS MODEL
      # This requires a prescribed anoxic fraction of 0.14
      Fborg <- p['Fborg.i']*(y['PO4']/PO4.i)^2
      Fbp <- Fborg/250 # Note, this would need to change if dynamic O2 is implemented
    }
    
    if (p['PO4.form']==4){ # Uses Shields and Mills 2017 formulation
      Fborg <- p['Fborg.i']*(y['PO4']/PO4.i)^2
      Fbp <- p['Fbp.i']*(Fborg/p['Fborg.i'])
    }
    
    if (p['PO4.form']==5){ # Uses Payne and Kump 2007 formulation
      Fbp <- p['Fbp.i']*(y['PO4']/PO4.i)
      Fborg <- p['Fborg.i']*(y['PO4']/PO4.i)
    }
    
    if (p['PO4.form']==6){  # Arbitrary temperature-dependent organic burial 
      Fborg <- Fworg * (temp.i.wt/temp.a.i)^arbFB.n # (note temp.a.i is weighted surface air temp from train.tracks fxn (first timestep))
      Fbp <- Fwp
    }
    
    # Sulfur Fluxes and associated Carbon fluxes
    Fwsulf <- p['Fwsulf.i']
    Fbsulf <- p['Fbsulf.i']*y['S']/S.i
    
    Fwsulf.C <- 2*Fwsulf # Assumes 1:2 molar ratio (S -> C)
    Fbsulf.C <- 2*Fbsulf # Assumes 1:2 molar ratio (S -> C)
    
    dDIC <- Fvolc + Fworg + Fwcarb - Fborg - Fbcarb + Fwsulf.C - Fbsulf.C
    dAlk <- Falkin - Falkout
    dd13C <- (Fworg*(p['Fworg.d13C']-y['d13C']) + Fwcarb*(p['Fwcarb.d13C']-y['d13C']) + Fextra*(d13C.force-y['d13C']) + Fvolc*(p['Fvolc.d13C']-y['d13C']) + p['DB']*Fborg)/y['DIC']
    dS <- Fwsulf - Fbsulf
    dPO4 <- Fwp - Fbp
    dMg <- 0
    dCa <- 0
    
    
    
    # ************************************ # 
    # TROUBLESHOOTING PRINT                #
    cat(paste("=======================================================",
              paste("time: ", t, sep=''),    #
              paste("CO2: ", round(pCO2,0), sep=''),  #
              paste("temp: ", round(temp.i.wt,1), sep=''),  #
              paste("Fwsil: ", round(Fwsil,1), sep=''),  #
              paste("Fvolc: ", round(Fvolc,1), sep=''),  #
              # paste("LBtemp: ", LBtemp, sep=''), #    
              if(t==0){"Time zero!!"}else{"********************"},  # debug...
              paste("This loop = ", this.loop, sep=''),
              paste("This run = ", this.run, " of ", n.idx, sep=''),
              "=======================================================", '\n',
              sep='\n'))
    # ************************************ # 
    
    
    
    # save and set results
    results <- c(dDIC, dAlk, dd13C, dS, dPO4, dMg, dCa)  
    
    res=list(results, pCO2, omega, Fwcarb, Fbcarb, Fborg, Fwsil, RCO2, Fvolc,
            Fworg, Fwsulf, Fbsulf, Conc.avg, Fwp, Fbp, temp.i.wt, runoff.i.wt, 
            runoff.i.wt.land, LBtemp, RBtemp, clim.state.int, iceline.lat)
  })
}

#### Uber CH2O-CHO Functions ####
# These functions are necessary for parallelizing the above CH2O-CHO functions
# [NOT ACTIVE IN THIS VERSION]
CH2OCHO.3.uber <- function(parameters){
  final.results <- as.data.frame(lsodes(y=y0,times=t,func=CH2OCHO,parms=parameters,rtol=1e-6))
  names(final.results) <- c('time','DIC', 'Alk','d13C','S','PO4','Mg','Ca','pCO2','omega', 'Fwcarb','Fbcarb', 
                            'Fborg','Fwsil','Fwsil.1','Fwsil.3','RCO2','Temp','Fworg',"Fwsulf","Fbsulf",
                            "q.1","q.2","q.3",
                            "q.avg","Conc.1","Conc.2","Conc.3","Conc.avg","Dw.1","Dw.2","Dw.3")
  list(parameters,final.results)
}




## FUNCTION: TRAIN TRACKS
## reads in the forcing and relevant files and runs a loop of simulations
train.tracks <- function(force.iters = FORCE.DF,
                         config.file = CONFIG.FN,
                         p.dir = parent.path){
  # ... START THE SIM
  # Modified from Donovan Bake: https://www.asciiart.eu/vehicles/trains 
  cat(paste("********************************************************************",
            "** ALL ABOARD!!!                                                  **",
            "                 _-====-__-======-__-========-_____-============-__",
            "               _(           _      __  _    _       _   _         _)",
            "            OO(            /   /_/ _/ / )_ /   /_/ / ) / )        )_",
            "           0  (_          (__ / / /_ (_/  (__ / / (_/ (_/          _)",
            "         o0     (_                                                _)",
            "        o         '=-___-===-_____-========-___________-===-dwb-='",
            "      .o                                _________                  ",
            "     . ______          ______________  |         |      _____      ",
            "   _()_||__|| ________ |            |  |_________|   __||___||__",
            "  (KLIR 2022| |      | |            | __Y______00_| |_         _|",
            " /-OO----OO''='OO--OO'='OO--------OO'='OO-------OO'='OO-------OO'=P",
            "===================================================================",
            sep='\n'))
  
  ## ... enjoy the graphic
  Sys.sleep(4)
  ## ... prompt on calc
  cat("Calculating climate at the first stop ... ", sep='\n')
  
  ## --------------------------------------------------------------- ## 
  ## FIRST - Copy the config.r file and save in the results dir      ## 
  old.path <- paste(p.dir, 'RUN', config.file, sep='/')
  new.path <- paste(saveDir, paste('copy_of_', config.file, sep=''), sep='/')
  # save result as a var to avoid TRUE/FALSE printing to console
  config.copy.result <- file.copy(from = old.path, to = new.path, overwrite = TRUE)
  ## --------------------------------------------------------------- ## 
  
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
  # loop through the forcing df
  for(RUN in 1:n.idx){
    # keep track of the run
    this.run <<- RUN
    ## [1] Set geography
    if(multi.geog == TRUE){ # if we use more than one geography, test the one set by force.df
      thisPaleogeo <<- geog.list[[force.iters$geography.force[RUN]]]
      landfrac_index <<- geog.vec[RUN]
      } else{thisPaleogeo <<- geog.list; landfrac_index <<- geog.vec} # otherwise just use the single geog
    
    ## [2] set other initial conditions
    pCO2 <- force.iters$co2.force[RUN]
    glwx <<- force.iters$glwx.force[RUN]  # goes to global environment because it's called in MEBM_solve-bistable-glwx.R script
    m <- force.iters$m.force[RUN]
    t.exc.start <- force.iters$t.start.force[RUN]
    t.exc.end <- force.iters$t.end.force[RUN]
    duration <- force.iters$duration.force[RUN]
    dt <- force.iters$dt.force[RUN]
    CONTROL.iter <- force.iters$CTRL.idx[RUN]
    # set to control iter for initialization
    SENSITIVITY.CTRL.iter <- force.iters[CTRL.idx == "Y"]$SENSITIVITY.force[1]
    # ... update the parList sensitivity vector with this iteration
    parList[[par.here]][[SENSITIVITY.PAR]] <- SENSITIVITY.CTRL.iter
    parList <<- parList  # I reset this because there's an error when I try to change the global environment of parList (<<-) with the line above
    # ******************************************** #
    # SET THE RUN NAME FOR THIS SPECIFIC OUTPUT    #
    # First remove odd chars from m
    msub1 <- gsub("[.]", "o", as.character(m)) # replace dot with small o
    msub2 <- gsub("[+]", "", msub1)  # replace plus with nothing
    RunName.idx <- paste(RunName_BASE, "co2-", pCO2, "_glwx-", glwx, "_m-", msub2, "_iter-", RUN, sep='')
    # ******************************************** #
    
    ## [3] Set climate state
    LBtemp <<- LBtemp.initial ; RBtemp <<- RBtemp.initial
    
    ## [4] initial mebm 
    dfMEBM_init <- MEBM_TRAIN.sensitivity(pCO2_in = pCO2, LandFrac_tx = thisPaleogeo, get_pars=FALSE )  
    Temp.bin.i <- round(dfMEBM_init$temp_C, temp_res)  # latitudinal temperature profile
    
    # set the new (last) clim state (i.e. that for the first time step)
    clim.state.last <<- clim.state <<- clim_state_ID(thisdf = dfMEBM_init)
    if(Bistable.Run=="Y"){TempBC_fun(thisClim = clim.state.last)
      LBtemp.initial <<- LBtemp ; RBtemp.initial <<- RBtemp}
    
    # ... calculate global mean temperature
    e_sq <- 0.00669437999014    # [m] eccentricity^2 (accounts for slight oblique shape of earth... negligible in calculation)
    deg2rad <- function(deg) {(deg * pi) / (180)}
    Lat_area_weight <<- (pi*Re*cos(deg2rad(dfMEBM_init$Lat)) / (180 * sqrt(1-e_sq * sin(deg2rad(dfMEBM_init$Lat))**2)))        # [m] meters per degree longitude (assigned to environment for CH2O-CHO)
    temp.i <- round(weighted.mean(x = dfMEBM_init$temp_C, w=Lat_area_weight), temp_res)   # [degC] area-weighted mean global temperature
    # ... set initial runoff value 
    runoff.i <- round(weighted.mean(x=dfMEBM_init$budykoRunoff_m_yr, w=dfMEBM_init$Land_Area_m2), q_res) # [m/yr] weighted by land area
    runoff0 <<- round(dfMEBM_init$budykoRunoff_m_yr, q_res)   # latitudinal runoff profile
    
    # --- UPDATE PARAMETER TO CURRENT ITER VALUE --- # 
    SENSITIVITY.iter <- force.iters$SENSITIVITY.force[RUN]
    parList[[par.here]][[SENSITIVITY.PAR]] <- SENSITIVITY.iter
    parList <<- parList  # I reset this because there's an error when I try to change the global environment of parList (<<-) with the line above
    
    ## [5] initialize geochem bcs
    # ... set boundary conditions
    bc = cbind(Age,Ca.i,Mg.i,k.run,clim.sens,degassing,pCO2,temp.i)  # bring in temp from MEBM initial output
    # get correct units
    Ca.i = bc[1,'Ca.i']/1000 # [mol/kg]
    Mg.i = bc[1,'Mg.i']/1000 # [mol/kg]
    # put in global env.
    bc <<- bc
    
    ## [6] Set initial CO2 and temp for weathering scaling
    if(scale.to.init.co2 == TRUE & RUN == 1){  # update only at first iter
      pCO2.i <<- bc[1,'pCO2']
      temp.a.i <<- as.vector(bc[1,'temp.i'])  # [degC] initial ocean surface temp
    } else if(scale.to.init.co2 == FALSE){  # update for every iter
      pCO2.i <<- bc[1,'pCO2']
      temp.a.i <<- as.vector(bc[1,'temp.i']) # [degC] initial ocean surface temp
    } 
    
    ## [7] Set c cycle perturbation and d13C
    tvolc <- m
    Finterp.volc = c(rep(0,(t.exc.start/dt+1)), rep(tvolc/(n*dt), n) ,rep(0,(duration-t.exc.end)/dt))
    d13Cinterp.pert = c(rep(0,(t.exc.start/dt+1)), rep(d13C.extra, n) ,rep(0,(duration-t.exc.end)/dt))
    # build function to interpolate perturbation
    event <<- data.frame(t,Finterp.volc)
    event.d13C <<- data.frame(t, d13Cinterp.pert)
    force.C <<- approxfun(x=event[,1],y=event[,2],method="linear",rule=2) # function for calculating Fextra
    force.d13C <<- approxfun(x=event.d13C[,1],y=event.d13C[,2],method="linear",rule=2)
    
    ## [8] Initialize carbonate system
    oceanV <<- 1.4e21 # L of ocean
    temp.o.i <<- temp.a.i - 10 ### assume a 10C gradient between air and deep ocean T
    salinity <<- 35 # salinity in salinity units
    sal <<- (salinity/1000)+1 # salinity in kg/L
    pH.i <<- 8.2
    carb.parms.i <<- co2sys.init(pH = pH.i,pCO2 = pCO2.i,temp=temp.o.i,sal=salinity,Ca2=Ca.i,Mg=Mg.i)
    omega.i <<- carb.parms.i$OmegaCalcite
    RCO2.i <<- pCO2.i/pCO2.i
    DIC.i <<- carb.parms.i$DIC*oceanV
    Alk.i <<- carb.parms.i$ALK*oceanV
    
    ## [9] initial fluxes and carbon isotope values
    # Initial carbon cycle fluxes and parameters (mol/yr)
    Fvolc.i <<- 8e12 # solid Earth degassing
    Fwcarb.i <<- 12e12 - (Fvolc.i-8e12)/2 # carbonate weathering
    Fbcarb.i <<- 20e12 + (Fvolc.i-8e12)/2 # carbonate burial
    Fworg.i <<- 8e12 # organic C weathering
    Fborg.i <<- 8e12 # burial of organic carbon
    # Initial d13C Values
    d13C.i <<- 0 # initial DIC permil
    Fwcarb.d13C <<- 0
    Fvolc.d13C <<- -5 
    DB <<- 27 # Cap-delta between inorganic carbon and organic carbon
    Fworg.d13C <<- d13C.i - (Fwcarb.i*(Fwcarb.d13C-d13C.i) + Fvolc.i*(Fvolc.d13C-d13C.i) + DB*Fborg.i)/(Fworg.i)
    S.i <<- 10e6*0.5e12 #p['Fwsulf.i'] # mol sulfur; assumes a 10 Ma residence time of S in the ocean; modern-day estimates are 5e19 (Broecker and Peng 1982), but this includes a substantial gypsum component
    PO4.i <<- 3.1e15 # mol P
    # (note, S is turned off in this version of the code)
    
    ## [10] Initialize paramfunc
    p <- paramfunc()
    y0 <- c(DIC=DIC.i, Alk=Alk.i, d13C=d13C.i, S=S.i, PO4=PO4.i, Mg=Mg.i, Ca=Ca.i) 
    names(y0) <- c(names(y0)[1:3],"S","PO4","Mg","Ca")
    # remove unnecessary param
    p <- paramfunc(remove=c("nbox"))  
    # update phosph.
    p['PO4.form'] <- 1  # 1=no coupling; 2=Van Cappellen formulation; 3=Bergman et al., 2004; 4=Shields and Mills 2017; 5=Payne and Kump 2007
    p <- c(p,nbox=length(xin))    # nbox set to length of nodes in MEBM
    # assign to global env
    p <<- p
    y0 <<- y0
    
    ## [11] initialize weathering
    Fwsil.vec.setup <<- numeric(length=length(xin))
    Fwcarb.vec.setup <<- numeric(length=length(xin))
    Conc.vec.setup <<- numeric(length=length(xin))
    Conc0.vec.setup <<- numeric(length=length(xin))
    Conc0.carb.setup <<- numeric(length=length(xin))
    Conc.carb.setup <<- numeric(length=length(xin))
    Dw.vec.setup <<- numeric(length=length(xin))
    area.frac.vec.setup <<- dfMEBM_init$Frac_of_tot_LandArea  # these go to global env. for the parList save at end
    area.vec.setup <<- dfMEBM_init$Land_Area_m2   # these go to global env. for the parList save at end
    # ... set initial runoff and temperature conditions
    if(scale.to.init.co2 == TRUE & CONTROL.iter == "Y"){
      q0.setup <<- runoff.i # this is just for the log file
      runoff.INIT <<- runoff0 #  this is used for MCfunc
      Temp.bin.INIT <<- Temp.bin.i
    } else if(scale.to.init.co2 == FALSE){
      q0.setup <<- runoff.i # this is just for the log file
      runoff.INIT <<- runoff0  # this is used for MCfunc
      Temp.bin.INIT <<- Temp.bin.i
    }
    # calculate the scalar to go from C*q to Fwsil-init
    for (i in 1:length(xin)){   # length(xin) is equal to nbox later on
      # Runoff  
      q.setup <- runoff0[i]  
      
      # Using a Maher and Chamberlain-style concentration feedback with a Winnick and Maher 
      # open-system calculation for Ceq
      Conc.results.setup <- MCfunc.trans(pCO2 = pCO2, q0 = runoff.INIT[i], q = q.setup, Temperature = Temp.bin.i[i],
                                         Ea = p['Ea'], Ts=p['Ts'], soilCO2.use=p['soil.CO2.flag'],
                                         Temp.init = Temp.bin.INIT[i])
      
      # Silicate Weathering
      Conc0.vec.setup[i] <- Conc.results.setup$Granite$Conc0*p['silicate.lith'] + Conc.results.setup$Basalt$Conc0*(1-p['silicate.lith'])
      Conc.vec.setup[i] <- Conc.results.setup$Granite$Conc*p['silicate.lith'] + Conc.results.setup$Basalt$Conc*(1-p['silicate.lith'])
      Dw.vec.setup[i] <- Conc.results.setup$Granite$Dw*p['silicate.lith'] + Conc.results.setup$Basalt$Conc*(1-p['silicate.lith'])
      # Carbonate Weathering
      Conc0.carb.setup[i] <- Conc.results.setup$Carbonate$Conc0
      Conc.carb.setup[i] <- Conc.results.setup$Carbonate$Conc
      
      # Fwsil vec
      Fwsil.vec.setup[i] <- q.setup*Conc.vec.setup[i]*area.vec.setup[i]
      Fwsil.vec.setup[which(is.na(Fwsil.vec.setup))] <- 0  # replace zero runoff NAs with zero wx
      # Fwcarb vec
      Fwcarb.vec.setup[i] <- q.setup*Conc.carb.setup[i]*area.vec.setup[i]
      Fwcarb.vec.setup[which(is.na(Fwcarb.vec.setup))] <- 0  # replace zero runoff NAs with zero wx
    }
    # Set the weathering scaler
    if((scale.to.init.co2 == TRUE & CONTROL.iter == "Y") | scale.to.init.co2 == FALSE){
      Fwsil.global.unscaled <- sum(Fwsil.vec.setup, na.rm=T)
      Fwcarb.global.unscaled <- sum(Fwcarb.vec.setup, na.rm=T)
      Fwsil_scalar <<-  Fvolc.i / Fwsil.global.unscaled
      Fwcarb_scalar <<- Fwcarb.i / Fwcarb.global.unscaled
    }
    
    Fwsil.global <- sum(Fwsil.vec.setup, na.rm=T) * Fwsil_scalar
    Fwcarb.global <- sum(Fwcarb.vec.setup, na.rm=T) * Fwcarb_scalar
    
    # -- RECORD SUMMARY RESULTS
    if(RUN == 1){
      outdf.sum1 <- force.iters[1,]
      outdf.sum2 <- as.data.table(cbind(pCO2, temp.i, runoff.i, glwx, RBtemp, LBtemp, RUN))
      colnames(outdf.sum2) <- c("pCO2", "temp_C", "runoff_land.wt", "glwx",
                                "RBtemp", "LBtemp", "iter")
      outdf.sum <- as.data.table(cbind(outdf.sum1, outdf.sum2))
      outdf.sum$SENSITIVITY_PAR <- SENSITIVITY.PAR
      outdf.sum$SENSITIVITY_VAL <- SENSITIVITY.iter
      outdf.sum$Fwsil_unequilibrated <- Fwsil.global
      outdf.sum$Fwsil_scalar <- Fwsil_scalar
      outdf.sum$climstate <- clim.state
      outdf.sum$geog <- landfrac_index
      outdf.sum$RunName <- RunName.idx
    } else{
      tmpdf.sum1 <- force.iters[RUN, ]
      tmpdf.sum2 <- as.data.table(cbind(pCO2, temp.i, runoff.i, glwx, RBtemp, LBtemp, RUN))
      colnames(tmpdf.sum2) <- c("pCO2", "temp_C", "runoff_land.wt", "glwx",
                                "RBtemp", "LBtemp", "iter")
      tmpdf.sum <- as.data.table(cbind(tmpdf.sum1, tmpdf.sum2))
      tmpdf.sum$SENSITIVITY_PAR <- SENSITIVITY.PAR
      tmpdf.sum$SENSITIVITY_VAL <- SENSITIVITY.iter
      tmpdf.sum$Fwsil_unequilibrated <- Fwsil.global
      tmpdf.sum$Fwsil_scalar <- Fwsil_scalar
      tmpdf.sum$climstate <- clim.state
      tmpdf.sum$geog <- landfrac_index
      tmpdf.sum$RunName <- RunName.idx
      # add to existing df
      outdf.sum <- as.data.table(rbind(outdf.sum, tmpdf.sum))
    }
    
    
    ## -- RUN THE TIME-TRANSIENT PART -- ## 
    # track where we are in numerical solver
    loop.idx <<- 0     # must start at zero [DO NOT CHANGE]
    # start run
    run.start <- proc.time()
    output <- as.data.frame(deSolve::rk4(y=y0, times=t, func=CH2OCHO, parms=p))
    run.end <- proc.time()
    elapsed.time <- (run.end[3] - run.start[3])
    simDuration <<- paste(round(elapsed.time/60,0),"min and", round((((elapsed.time-floor(elapsed.time))*100*60)/100),0),"sec")
    print(simDuration)
    colnames(output) <- c('time','DIC', 'Alk','d13C',"S","PO4",'Mg','Ca','pCO2','omega','Fwcarb','Fbcarb',
                          'Fborg','Fwsil','RCO2','Fvolc','Fworg', "Fwsulf","Fbsulf","Conc.avg","Fwp","Fbp", 
                          "temp_C", "runoff_budyko", "runoff_budyko_LANDwt", 
                          "LBtemp", "RBtemp", "climstate_int", "iceline_lowestLat")
    #... decode climate state
    output$climstate <- clim_state_INT.2.Str(output$climstate_int)
    #... add scalars
    output$Fwsil_scalar <- Fwsil_scalar
    output$Fwcarb_scalar <- Fwcarb_scalar
    #... save sensitivity sim info
    output$SENSITIVITY_PARAM <- SENSITIVITY.PAR
    output$SENSITIVITY_VALUE <- SENSITIVITY.iter
    output$pCO2.init <- pCO2
    output$iter <- RUN
    output$glwx <- glwx
    output$m <- m
    output$tstart <- t.exc.start
    output$tend <- t.exc.end
    output$duration <- duration
    output$dt <- dt
    output$CTRL <- CONTROL.iter
    
    # BUILD OUTPUT DF
    if(RUN==1){
      # add discharge and wx
      dfMEBM_init$Discharge_m3_yr <- dfMEBM_init$budykoRunoff_m_yr * dfMEBM_init$Land_Area_m2
      dfMEBM_init$Fwsil <- Fwsil.vec.setup * Fwsil_scalar
      dfMEBM_init$Fwcarb <- Fwcarb.vec.setup * Fwcarb_scalar
      # track initial climate
      dfMEBM_init$SENSITIVITY_PAR <- SENSITIVITY.PAR
      dfMEBM_init$SENSITIVITY_VEC <- SENSITIVITY.iter
      dfMEBM_init$RUN.idx <- RUN
      # store in output df
      output.final <- output
      output.initClim <- dfMEBM_init
    } else{
      # add discharge and wx
      dfMEBM_init$Discharge_m3_yr <- dfMEBM_init$budykoRunoff_m_yr * dfMEBM_init$Land_Area_m2
      dfMEBM_init$Fwsil <- Fwsil.vec.setup * Fwsil_scalar
      dfMEBM_init$Fwcarb <- Fwcarb.vec.setup * Fwcarb_scalar
      # track initial climate
      dfMEBM_init$SENSITIVITY_PAR <- SENSITIVITY.PAR
      dfMEBM_init$SENSITIVITY_VEC <- SENSITIVITY.iter
      dfMEBM_init$RUN.idx <- RUN
      # store in output df
      output.final <- as.data.table(rbind(output.final, output))
      output.initClim <- as.data.table(rbind(output.initClim, dfMEBM_init))
    }
    
    
    # ********************************************************************************** #
    # save all results at each step
    saveName.full <- paste(saveDir, '/', RunName_BASE, '-', 'FULL.RDS',  sep='')
    saveName.sum <- paste(saveDir, '/', RunName_BASE, '-', 'SUMMARY.RDS', sep='')
    saveName.initclim <- paste(saveDir, '/', RunName_BASE, '-', 'INITCLIM.RDS', sep='')
    # remove files if they exist (to prevent conflicted saves) (Sorry for the hacky solution)
    # ... in future, see also: https://stackoverflow.com/questions/44456298/issue-with-saverds
    if(file.exists(saveName.full)){file.remove(saveName.full)}
    if(file.exists(saveName.sum)){file.remove(saveName.sum)}
    if(file.exists(saveName.initclim)){file.remove(saveName.initclim)}
    # save new file
    saveRDS(output.final, file=saveName.full)
    saveRDS(outdf.sum, file=saveName.sum)
    saveRDS(output.initClim, file=saveName.initclim)
    # ********************************************************************************** #
    
    CH2O_CHO_mode <<- "transient"
    # ------------------------- #
    # WRITE LOG FILE            #
    MEBM_parameters <<- MEBM_TRAIN.sensitivity(pCO2_in = pCO2[1], LandFrac_tx = thisPaleogeo, get_pars=TRUE)  
    TRAIN_parameters <<- sapply(saveTRAIN, function(x)get(x),simplify=F,USE.NAMES=T)  # this is the same thing as PARlist
    paramfunc_parameters <<- p
    RunName <<- RunName.idx
    writeLogFile(userName = myName, userNotes = myNotes)
    # Save output file          #
    saveName <- paste(saveDir, '/', RunName.idx, '.RDS',  sep='')
    saveRDS(output, file=saveName)
    # save parameters 
    parameter_list <- list(MEBM_parameters, TRAIN_parameters, paramfunc_parameters) ; names(parameter_list) <- c("MEBMpars", "TRAINpars", "ParamFunc_pars")
    savePar <- paste(saveDir, '/', RunName.idx, '_PAR-LIST.RDS', sep='')
    saveRDS(parameter_list, savePar)
    # ------------------------- #
  }
  
}



