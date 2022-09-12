# -------------------------------------- # 
# Function to build paleogeo data.tables # 
# from paleogeography shape files        #
# -------------------------------------- #
library(dplyr) 

#... reads in latitude, longitude, and land frac data
#... degSquared accounts for square grid cells
LandFrac_fun <- function(xyzLat, xyzLon, xyzFrac, degSquared = 2){
  Re <- 6378137.0    # [m] radius of earth in 
  e_sq <- 0.00669437999014    # [m] eccentricity^2 (accounts for slight oblique shape of earth... unimportant, but hey it's here)
  
  deg2rad <- function(deg) {(deg * pi) / (180)}
  #... get areas of each latitude bin
  # equations from: https://en.wikipedia.org/wiki/Latitude under the "Length of a degree of Latitude" section
  m_per_gridLat <- degSquared * (pi*Re*(1-e_sq) / (180 * ( (1-e_sq * sin(deg2rad(xyzLat))**2)**(3/2) )))   # [m] meters per degree latitude
  m_per_gridLon <- degSquared * (pi*Re*cos(deg2rad(xyzLat)) / (180 * sqrt(1-e_sq * sin(deg2rad(xyzLat))**2)))        # [m] meters per degree longitude 
  m_grid_area <- m_per_gridLat * m_per_gridLon   # [m^2] area of each grid cell in meters
  
  #... calculate the fraction of land area in each lat belt
  newxyz <- as.data.table(cbind(xyzLat, xyzLon, xyzFrac, m_grid_area))
  colnames(newxyz) <- c('Lat', 'Lon', 'Frac', 'AreaPerCell')
  # ... divide xyz frac by 100 or 1 depending on max xyzfrac
  frac.denom <- ifelse(max(newxyz$Frac) > 1.1, 100, 1)
  newxyz %>%
    group_by(Lat) %>%
    mutate(LandFrac = (sum(Frac)/frac.denom)/length(Lon)) %>%
    summarise_all(mean) -> newxyzOUT
  
  #... calculate the area at each latitudinal belt
  m_per_latBelt <- newxyzOUT$AreaPerCell * (360 / degSquared)  
  area_per_latBelt <- newxyzOUT$LandFrac * m_per_latBelt
  
  #... create the output file (fraction of total land per latitude ; 
  #                            and total land area)
  allmyland <- sum(area_per_latBelt)   # [m 2] land area on the planet
  land_FracOfTotal <- area_per_latBelt / allmyland
  myOut <- as.data.table(cbind(newxyzOUT$Lat, land_FracOfTotal, area_per_latBelt, newxyzOUT$LandFrac))
  colnames(myOut) <- c('Lat', 'Frac_of_tot_LandArea', 'Land_Area_m2', 'LandFrac_inLatBelt_forAlbedo')
  
  # return result
  return(myOut)
}