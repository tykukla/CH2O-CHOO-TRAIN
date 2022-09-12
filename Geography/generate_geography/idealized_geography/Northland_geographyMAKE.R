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
# ... Note -- I realized the formulation was S-hemisphere 
# ... after naming it northland... but the model runs 
# ... symmetrically and they're functionally identical
# (southland is changed to northland later by multiplying lat by -1)

library(data.table)

rm(list=ls())

# working dir is file dir
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#... Bring in the land data
LandFracFileName <- "PaleoGeo_dataframes_2deg.RDS"
myLandFrac <- readRDS(LandFracFileName)

#... pull out modern land file
modland <- as.data.table(myLandFrac[["Ma_0"]])
#... create an idealized "cat-eye" planet
Re <- 6.37e6        # [m] radius of the earth
e_sq <- 0.00669437999014    # [m] eccentricity^2 (accounts for slight oblique shape of earth... negligible in calculation)
global_surface_area <- 4*pi*(Re**2)     # [m^2] surface area of Earth
global_landArea <- sum(modland$Land_Area_m2)   # total land area

# ... distribute global land area so there is the same land frac in each lat bin above some minimum latitude
# ... start by getting latitudinal area and calculating the latitude that land reaches for a given fraction of land cover held constant
landFrac.vector <- seq(0.1, 1, length=100)
total.lat.area.wt <- (pi*Re*cos(deg2rad(modland$Lat)) / (180 * sqrt(1-e_sq * sin(deg2rad(modland$Lat))**2)))   
lat.area.frac.total <- total.lat.area.wt / sum(total.lat.area.wt)
lat.area <- global_surface_area * lat.area.frac.total
df.lat.area <- as.data.table(cbind(lat.area, modland$Lat)) ; colnames(df.lat.area) <- c("lat area", "lat deg")
thislat <- vector()
for(idx in 1:length(landFrac.vector)){
  thislandFrac <- landFrac.vector[idx]
  leftover_land_old <- global_landArea # Initialize the leftover land tracker
  leftover_land_new <- vector()
  for(latx in 1:length(total.lat.area.wt)){
    leftover_land_new[latx] <- leftover_land_old - (thislandFrac * lat.area[latx])
    leftover_land_old <- leftover_land_new[latx]
    print(leftover_land_old)
  }
  # put in df
  outdf <- as.data.table(cbind(leftover_land_new, modland$Lat)) ; colnames(outdf) <- c("land_used", "lat")
  thislat[idx] <- stats::approx(x=outdf$land_used, y=outdf$lat, xout=0)$y # find latitude where leftover land is zero
}

# lat.landfrac <- approx(thislat, landFrac.vector, xout=-10)$x

plot(landFrac.vector, thislat, type='l', col='darkred')

# ... okay let's just use the 100% area filled by land thing
# calculate the inputs for the output df
# find latitude of land extent
land_alb_frac <- 1
minlat.idx <- which.min(abs(outdf$land_used))
minlat <- outdf$lat[minlat.idx]
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
# change southland to northland
newland$Lat <- newland$Lat * -1

# did we get all of the global land area?
sum(newland$Land_Area_m2) / global_landArea
# does it look right?
plot(newland$Lat, newland$Land_Area_m2, type='l')



# ********************************************************* # 
# SAVE RESULT --------------------------------------------- #
saveRDS(newland, file='Northland_geography_2deg.RDS')  #
# ********************************************************* #


