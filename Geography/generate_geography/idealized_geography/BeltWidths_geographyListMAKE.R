# --------------------------------------------------- # 
# Script to parameterize the southern hemisphere land # 
# planet so that total global land area is similar    #
# to the modern                                       #
# ---                                                 #
# T Kukla (ColoState Univ.)                           #
#                                                     #
# Created: 2/26/2022                                  #
# Last updated: 8/02/2022 (TK)                        #
# --------------------------------------------------- #
library(data.table)
library(pracma)
library(ggplot2)
library(ggthemes)

# generate belt-type geographies of different widths centered on some defined latitude

rm(list=ls())

# working dir is file dir
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# ... set width of the latitudinal belt
belt.widths <- c(20, 30, 40)  # degrees latitude

#... Bring in the land data
LandFracFileName <- "PaleoGeo_dataframes_2deg.RDS"
myLandFrac <- readRDS(LandFracFileName)

#... pull out modern land file
modland <- as.data.table(myLandFrac[["Ma_0"]])
#... create an idealized "cat-eye" planet
Re <- 6.37e6        # [m] radius of the earth
e_sq <- 0.00669437999014    # [m] eccentricity^2 (accounts for slight oblique shape of earth... negligible in calculation)
global_surface_area <- 4*pi*(Re**2)     # [m^2] surface area of Earth
# global_landArea <- sum(modland$Land_Area_m2)   # total land area

# ... distribute global land area so there is the same land frac in each lat bin above some minimum latitude
# ... start by getting latitudinal area and calculating the latitude that land reaches for a given fraction of land cover held constant
landFrac.vector <- seq(0.1, 1, length=100)
total.lat.area.wt <- (pi*Re*cos(deg2rad(modland$Lat)) / (180 * sqrt(1-e_sq * sin(deg2rad(modland$Lat))**2)))   
lat.area.frac.total <- total.lat.area.wt / sum(total.lat.area.wt)
lat.area <- global_surface_area * lat.area.frac.total
df.lat.area <- as.data.table(cbind(lat.area, modland$Lat)) ; colnames(df.lat.area) <- c("lat area", "lat deg")
thislat <- vector()

# ... loop through land-lines (lowest lat of land) then through latitudes
central.lat <- 20
# save result to list
out.list <- list()
for(land.line in 1:length(belt.widths)){
  # ... use the 100% area filled by land thing
  # calculate the inputs for the output df
  # find latitude of land extent
  land_alb_frac <- 1
  minlat <- central.lat - (belt.widths[land.line]/2) 
  newland <- modland # overwrite modland one at a time
  theserows <- newland[Lat >= minlat & Lat <= (minlat+belt.widths[land.line]), which = TRUE] # which statement returns row numbers: https://stackoverflow.com/questions/22408306/using-i-to-return-row-numbers-with-data-table-package
  # add land area
  newland[Lat >= minlat & Lat <= (minlat+belt.widths[land.line])]$Land_Area_m2 <- land_alb_frac * lat.area[theserows]
  newland[Lat < minlat | Lat > (minlat+belt.widths[land.line])]$Land_Area_m2 <- 0
  # add fraction of land in lat bin
  newland[Lat >= minlat & Lat <= (minlat+belt.widths[land.line])]$LandFrac_inLatBelt_forAlbedo <- land_alb_frac
  newland[Lat < minlat | Lat > (minlat+belt.widths[land.line])]$LandFrac_inLatBelt_forAlbedo <- 0
  # add frac of total land area
  newland$Frac_of_tot_LandArea <- newland$Land_Area_m2 / sum(newland$Land_Area_m2)
  
  # did we get all of the global land area?
  # sum(newland$Land_Area_m2) / global_landArea
  # does it look right?
  # plot(newland$Lat, newland$Land_Area_m2, type='l')
  
  # save to list
  out.list[[land.line]] <- newland
  
  # save to df
  if(land.line == 1){
    outdf <- newland
    outdf$minlat <- minlat
  } else{
    newland$minlat <- minlat
    outdf <- as.data.table(rbind(outdf, newland))
  }
}

names(out.list) <- paste("Belt_width_", belt.widths, sep='')


# check that it looks okay

ggplot(outdf) + 
  geom_line(aes(x=Lat, y=Land_Area_m2, color=minlat, group=minlat)) + 
  theme_few()


# ********************************************************************* # 
# SAVE RESULT --------------------------------------------------------- #
saveRDS(out.list, file='Subtropicbelt-Widths_geography_2deg-LIST.RDS')  #
# ********************************************************************* #


