# --------------------------------------------------- # 
# Script to parameterize the southern hemisphere land # 
# planet so that total global land area is similar    #
# to the modern                                       #
# ---                                                 #
# T Kukla (ColoState Univ.)                           #
#                                                     #
# Created: 2/26/2022                                  #
# Last updated: 2/26/2022 (TK)                        #
# --------------------------------------------------- #
library(data.table)

rm(list=ls())

# working dir is file dir
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#... Bring in the land data
LandFracFileName <- "PaleoGeo_dataframes_2deg.RDS"
myLandFrac <- readRDS(LandFracFileName)

deg2rad <- function(deg) {(deg * pi) / (180)}

#... pull out modern land file
modland <- as.data.table(myLandFrac[["Ma_0"]])
#... create an idealized "cat-eye" planet
Re <- 6.37e6        # [m] radius of the earth
e_sq <- 0.00669437999014    # [m] eccentricity^2 (accounts for slight oblique shape of earth... negligible in calculation)
global_surface_area <- 4*pi*(Re**2)     # [m^2] surface area of Earth
global_landArea <- sum(modland$Land_Area_m2)   # total land area

# ... distribute global land area into two mid-lat zones that are symmetrical. 
# ... start with a defined abs(latitude) then add land north/south (away from equator) 
# ... ... until we cover ~40% of Earth surface area (as modern world)
start.abs.lat <- 30
landFrac.vector <- seq(0.1, 1, length=100)
total.lat.area.wt <- (pi*Re*cos(deg2rad(modland$Lat)) / (180 * sqrt(1-e_sq * sin(deg2rad(modland$Lat))**2)))   
lat.area.frac.total <- total.lat.area.wt / sum(total.lat.area.wt)
lat.area <- global_surface_area * lat.area.frac.total
df.lat.area <- as.data.table(cbind(lat.area, modland$Lat)) ; colnames(df.lat.area) <- c("lat area", "lat deg")

newland <- modland # overwrite modland one at a time
newland$Land_Area_m2 <- 0 ; newland$LandFrac_inLatBelt_forAlbedo <- 0 ; newland$Frac_of_tot_LandArea <- 0

land.fill <- TRUE   # marks whether we keep adding land, or whether we're done
idx <- 1
while(land.fill == TRUE){
  # find the latitude columns
  if(idx == 1){
    this.idx <- which.min(abs(abs(newland$Lat) - start.abs.lat))
    reverse.idx <- nrow(newland) - this.idx + 1
    my.idx <- c(this.idx, reverse.idx)
    # newland[this.idx] ; newland[reverse.idx]
  }
  
  # add land area
  newland[my.idx]$Land_Area_m2 <- lat.area[my.idx]
  # add fraction of land in lat bin
  newland[my.idx]$LandFrac_inLatBelt_forAlbedo <- 1
  
  # update idx
  idx <- idx + 1
  my.idx <- c(my.idx[1] - 1, my.idx[2] + 1)
  # newland[my.idx]
  
  # update leftover land
  leftover_land <- global_landArea - sum(newland$Land_Area_m2)
  
  if(leftover_land < 0){
    land.fill <- FALSE
  }
}
newland$Frac_of_tot_LandArea <- newland$Land_Area_m2 / sum(newland$Land_Area_m2)

# see how it looks
plot(newland$Lat, newland$Land_Area_m2, type='l')
plot(newland$Lat, newland$LandFrac_inLatBelt_forAlbedo, type='l')
plot(newland$Lat, newland$Frac_of_tot_LandArea, type='l')

sum(newland$Land_Area_m2) / global_landArea
# ... distribute global land area so there is the same land frac in each lat bin above some minimum latitude


# ********************************************************* # 
# SAVE RESULT --------------------------------------------- #
saveRDS(newland, file='Midlat_geography_2deg.RDS')          #
# ********************************************************* #


