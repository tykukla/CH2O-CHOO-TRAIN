# --------------------------------------------------- # 
# Script to parameterize the cat-eye planet so that   #
# total global land area is similar to the modern     #
# ---                                                 #
# T Kukla (Stanford Univ.)                            #
#                                                     #
# Created: 2/26/2021                                  #
# Last updated: 2/26/2021 (TK)                        #
# --------------------------------------------------- #
library(data.table)

rm(list=ls())

# working dir is file dir
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#... Bring in the land data
LandFracFileName <- "PaleoGeo_dataframes_2deg.RDS"
myLandFrac <- readRDS(LandFracFileName)

#... pull out modern land file
modland <- as.data.table(myLandFrac[["Ma_0"]])
global_landArea <- sum(modland$Land_Area_m2)

#... create an idealized "cat-eye" planet
Re <- 6.37e6        # [m] radius of the earth
e_sq <- 0.00669437999014    # [m] eccentricity^2 (accounts for slight oblique shape of earth... negligible in calculation)
# set latitudinal land fraction equal to the fraction of global area that is land today
global_surface_area <- 4*pi*(Re**2)     # [m^2] surface area of Earth
latitudinal_landFrac <- global_landArea / global_surface_area
deg2rad <- function(deg) {(deg * pi) / (180)}
# set thisPaleogeo to modern (all modern data will be over-written)
thisPaleogeo <- modland
TEMP_lat_area_weight <- latitudinal_landFrac * (pi*Re*cos(deg2rad(modland$Lat)) / (180 * sqrt(1-e_sq * sin(deg2rad(modland$Lat))**2)))  # MULTIPLIED BY 0.4 BECAUSE IT'S THE LAND FRAC IN EACH LAT BELT
thisPaleogeo$Frac_of_tot_LandArea <- TEMP_lat_area_weight / sum(TEMP_lat_area_weight)  # rep(1/nrow(thisPaleogeo), nrow(thisPaleogeo))
thisPaleogeo$LandFrac_inLatBelt_forAlbedo <- latitudinal_landFrac  # constant if frac_of_tot is weighted by lat_area
thisPaleogeo$Land_Area_m2 <- global_landArea * thisPaleogeo$Frac_of_tot_LandArea

# quick visual check
plot(modland$Land_Area_m2, type='l')
lines(thisPaleogeo$Land_Area_m2, type='l', col='blue')


# ******************************************************* # 
# SAVE RESULT ------------------------------------------- #
saveRDS(thisPaleogeo, file='Cat-eye_geography_2deg.RDS')  #
# ******************************************************* #
