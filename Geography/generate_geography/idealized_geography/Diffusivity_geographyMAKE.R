# --------------------------------------------------- # 
# Script to parameterize the subtropic-limited        #
# cat-eye and a north-oriented geography for tests    # 
# with the diffusivity variations                     #
# ---                                                 #
# T Kukla (Stanford Univ.)                            #
#                                                     #
# Created: 8/04/2021                                  #
# Last updated: 8/04/2021 (TK)                        #
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
latitudinal_landFrac <- 0.5 # global_landArea / global_surface_area
deg2rad <- function(deg) {(deg * pi) / (180)}


## -- CAT-EYE from -25 to 25 Lat -- ## 
max.lat <- 28; min.lat <- -28
# set thisPaleogeo to modern (all modern data will be over-written)
thisPaleogeo <- modland
TEMP_lat_area_weight <- latitudinal_landFrac * (pi*Re*cos(deg2rad(modland$Lat)) / (180 * sqrt(1-e_sq * sin(deg2rad(modland$Lat))**2)))  # MULTIPLIED BY 0.4 BECAUSE IT'S THE LAND FRAC IN EACH LAT BELT
thisPaleogeo$Frac_of_tot_LandArea <- TEMP_lat_area_weight / sum(TEMP_lat_area_weight)  # rep(1/nrow(thisPaleogeo), nrow(thisPaleogeo))
thisPaleogeo$LandFrac_inLatBelt_forAlbedo <- latitudinal_landFrac  # constant if frac_of_tot is weighted by lat_area
thisPaleogeo$Land_Area_m2 <- global_landArea * thisPaleogeo$Frac_of_tot_LandArea

# restrict to lat.lims
thisPaleogeo$Frac_of_tot_LandArea <- ifelse(thisPaleogeo$Lat > max.lat | thisPaleogeo$Lat < min.lat, 
                                            0, thisPaleogeo$Frac_of_tot_LandArea)
thisPaleogeo$Land_Area_m2 <- ifelse(thisPaleogeo$Lat > max.lat | thisPaleogeo$Lat < min.lat, 
                                    0, thisPaleogeo$Land_Area_m2)
thisPaleogeo$LandFrac_inLatBelt_forAlbedo <- ifelse(thisPaleogeo$Lat > max.lat | thisPaleogeo$Lat < min.lat, 
                                                    0, thisPaleogeo$LandFrac_inLatBelt_forAlbedo)


# quick visual check
plot(modland$Land_Area_m2, type='l')
lines(thisPaleogeo$Land_Area_m2, type='l', col='blue')


## -- NORTHLAND north of 25 Lat -- ## 
# ... distribute global land area so there is the same land frac in each lat bin above some minimum latitude
# ... start by getting latitudinal area and calculating the latitude that land reaches for a given fraction of land cover held constant
landFrac.vector <- seq(0.1, 1, length=100)
total.lat.area.wt <- (pi*Re*cos(deg2rad(modland$Lat)) / (180 * sqrt(1-e_sq * sin(deg2rad(modland$Lat))**2)))   
lat.area.frac.total <- total.lat.area.wt / sum(total.lat.area.wt)
lat.area <- global_surface_area * lat.area.frac.total
df.lat.area <- as.data.table(cbind(lat.area, modland$Lat)) ; colnames(df.lat.area) <- c("lat area", "lat deg")
thislat <- vector()

# ... loop through land-lines (lowest lat of land) then through latitudes
lowest.lat <- c(40)
# ... use the 100% area filled by land thing
# calculate the inputs for the output df
# find latitude of land extent
land_alb_frac <- 1
minlat <- lowest.lat * -1
newland <- modland # overwrite modland one at a time
theserows <- newland[Lat <= minlat, which = TRUE] # which statement returns row numbers: https://stackoverflow.com/questions/22408306/using-i-to-return-row-numbers-with-data-table-package
# add land area
newland[Lat <= minlat]$Land_Area_m2 <- land_alb_frac * lat.area[theserows]
newland[Lat > minlat]$Land_Area_m2 <- 0
# add fraction of land in lat bin
newland[Lat <= minlat]$LandFrac_inLatBelt_forAlbedo <- land_alb_frac
newland[Lat > minlat]$LandFrac_inLatBelt_forAlbedo <- 0
# add frac of total land area
newland$Frac_of_tot_LandArea <- newland$Land_Area_m2 / sum(newland$Land_Area_m2)


# -- CHECK FOR CONSISTENCY -- # 
# did we get all of the global land area?
sum(thisPaleogeo$Land_Area_m2) / global_landArea
sum(newland$Land_Area_m2) / global_landArea
# does it look right?
plot(newland$Lat, newland$Land_Area_m2, type='l')


# ********************************************************* # 
# ## SAVE RESULT ------------------------------------------ #
outlist <- list(thisPaleogeo, newland)
names(outlist) <- c(paste("CatEye_Latfrac_0", latitudinal_landFrac*10, "_LatLim_", max.lat, sep=''), 
                    paste("Northland_LowLat_", lowest.lat, sep=''))
saveRDS(outlist, "CatEye_Northland_DiffusivGeog.RDS")
# ********************************************************* # 
