#-------------------------------------------------------------------------------#
#                CH2O-CHOO-TRAIN [MEBM + CLiBeSO] Log file writer               # 
#                                                                               #
# Script to write log files for coupled CH2O-CHOO-TRAIN simulations. Saves all  #
# CLiBeSO and MEBM parameters so the simulations should be fully reproducible   #
#-------------------------------------------------------------------------------#

# *************************************************************** #
# this script writes the log file for each simulation 
#... it works as a function where the only inputs are "user notes"
#    and user name (where notes is an opportunity for the user to
#    write whatever they wish and name is simply the name of the 
#    person running the simulation)
# *************************************************************** #


#... Terms that must be defined in the global environment:

# [GENERAL]
# RunName -- identifier for this run
# saveDir -- location where output is saved
# simDuration -- string of how long the run took (often "XX minutes and XX seconds")

# [MEBM]
# LandFracFileName -- name of the land fraction file
# interpCO2 -- name (and sometimes location) of CO2 file

# [CLiBeSO]




writeLogFile <- function(userName=NULL, userNotes=NULL){
  setwd(saveDir)
  # create file name
  logName <- paste('LogFile_', RunName,'_', Sys.Date(), ".txt", sep='')
  fileConn <- file(logName)
  writeLines(c("---------------------------------", 
               "- CH2O-CHO Simulation Log Sheet -",
               "---------------------------------",
               "MEBM Program: MEBM_v2.0",
               "Carbon cycle Program: ",
               paste("Run by:", userName),
               paste("Run name:", RunName),
               paste("Run Start Time:", Sys.time()),
               paste("Run duration:", simDuration),
               paste("Run file path:", saveDir),
               "",
               "---- CH20-CHO TRAIN mode:",
               paste("Simulation mode:", CH2O_CHO_mode),
               "",
               "---- Input files:",
               paste("Land fraction file:", LandFracFileName),
               paste("pCO2 file:", CO2.file),
               paste("Seawater values file:", sw.file),
               paste("Volcanism file:", Berner.volc.file),
               # paste("Clim sensitivity file:", sens.file),
               "",
               "",
               "---------------------------------", 
               "MEBM initialization -------------",
               "---------------------------------",
               "",
               "------------------ Grid setup ------------------",
               paste("nNodes = ", length(xin)),
               paste("xin min = ", min(xin)),
               paste("xin max = ", max(xin)),
               paste("Index of landfrac data from list:", landfrac_index),
               paste("Right boundary condition flag = ", MEBM_parameters$RBC),
               paste("Left boundary condition flag = ", MEBM_parameters$LBC),
               "",
               "--------------- Global constants ---------------",
               paste("Earth's radius [m] = ", Re),
               paste("Diffusivity coefficient [m2 s-1] = ", MEBM_parameters$myD),
               paste("Solar constant (Q0) [W m-2] = ", MEBM_parameters$Q0),
               paste("OLR response constant (B0) [W m-2 degC-1] = ", MEBM_parameters$B0),
               paste("OLR constant depending on pCO2 (A0) [W m-2 degC-1] = ", MEBM_parameters$A0),
               "",
               "-------------- Climate constants ---------------",
               paste("Modern pCO2 for A0 calc [ppm] = ", MEBM_parameters$modCO2),
               paste("Relative humidity [0-1] = ", rel_hum),
               paste("Budyko omega [non-dimensional; (global avg=2.6)] = ", MEBM_parameters$budy_om),
               "",
               "-------------------- Albedo --------------------",
               paste("Albedo formulation = ", MEBM_parameters$alf_flag),
               paste("Ice threshold [degC] = ", ice_Threshold),
               paste("Ice albedo [0-1] = ", iceAlf),
               paste("Land albedo [0-1] = ", landAlf),
               paste("Ocean albedo [0-1] = ", oceanAlf),
               "",
               "---------------- Fixed constants ---------------",
               "These are values that cannot be accessed by the log sheet, but the value at time of writing is listed for reference.",
               "Gross Moist Stability [J kg-1] = 1.5e4",
               "Approximate hadley cell extent = 15 degrees",
               "Smoothness of hadley transition coefficient (xw) = 0.15",
               "Constants for E/P partitioning --- ",
               "--- Drag coefficient (Ch) = 1.5e-3",
               "--- LW feedback at surface [W m-2 K-1] = 0",
               "--- Wind speed parameterization [m s-1] = 4+abs(sin(pi*x/1.5))*4",
               "",
               "",
               "---------------------------------", 
               "Carbon cycle initialization -----",
               "---------------------------------",
               "",
               "------------ Chemistry setup vectors -----------",
               paste("Age for geochem initialization = ", bc[,"Age"]),
               "",
               paste("Atmospheric pCO2 [ppm] = ", bc[,"pCO2"]),
               "",
               paste("Ca concentration marine [mol/kg] = ", bc[,"Ca.i"]/1000),
               "",
               paste("Mg concentration marine [mol/kg] = ", bc[,"Mg.i"]/1000),
               "",
               # paste("Clim sensitivity = ", as.vector(bc[,'clim.sens'])),
               # "",
               paste("Volcanic degassing scalar [ ] = ", degassing),
               "",
               paste("Global temperature [degC] = ", as.vector(bc[ ,'temp.i'])),
               "",
               "------------- CLiBeSO domain setup -------------",
               paste("Model run duration [years] = ", duration),
               paste("timestep [years] = ", dt),
               paste("Perturbation start time [year] = ", t.exc.start),
               paste("Perturbation end time [year] = ", t.exc.end),
               paste("Mass of carbon injected [moles C] = ", m),
               paste("Ocean volume [L] = ", oceanV),
               paste("Deep ocean temperature [degC] = ", temp.o.i),
               paste("Ocean salinity [kg/L] = ", sal),
               "",
               "--------------- Carbonate system ---------------",
               paste("Initial ocean pH = ", pH.i),
               paste("Initial calcite saturation state = ", omega.i),
               paste("Initial RCO2 = ", RCO2.i),
               paste("Initial DIC [moles] = ", DIC.i),
               paste("Initial Alkalinity [moles] = ", Alk.i),
               "",
               "----------------- Carbon fluxes ----------------",
               paste("Background volcanic flux [mol/yr] = ", Fvolc.i),
               paste("Carbonate weathering flux [mol/yr] = ", Fwcarb.i),
               paste("Carbonate burial flux [mol/yr] = ", Fbcarb.i),
               paste("Organic carbon weathering flux [mol/yr] = ", Fworg.i),
               paste("Organic carbon burial flux [mol/yr] = ", Fborg.i),
               "",
               "---------------- Carbon isotopes ---------------",
               paste("DIC d13C [per mille] = ", d13C.i),
               paste("d13C of weathered carbonate [per mille] = ", Fwcarb.d13C),
               paste("d13C of weathered organic carbon [per mille] = ", Fworg.d13C),
               paste("d13C of volcanism [per mille] = ", Fvolc.d13C),
               paste("Inorganic - organic cap delta [per mille] = ", DB),
               "",
               "-------------- P and S reservoirs --------------",
               paste("Initial sulfur reservoir size [mol] = ", S.i),
               paste("Initial phosphate reservoir size [mol] = ", PO4.i),
               "",
               "",
               "---------------------------------", 
               "User Notes ----------------------",
               "---------------------------------",
               userNotes,
               "",
               "",
               "************************************************",
               "************************************************",
               "---- end of log file"
  ), fileConn)
  close(fileConn)
  
  print(paste("---- Log file successfully written ----"))
}


