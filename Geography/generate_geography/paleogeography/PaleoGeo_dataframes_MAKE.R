# ---------------------------------------------------------------- # 
# ... script to calculate latitudinal continental fractions        # 
#     for different time steps to feed into the MEBM / C-cycle     # 
#     model                                                        # 
# ---------                                                        # 
# T. Kukla (Stanford Univ. | Colostate Univ.)                      # 
# Developed for use with C-cycle models                            # 
# [D. Ibarra, K. Lau, J. Rugenstein]                               #
# Modified Aug 2022 by T. Kukla                                    #
# ---------------------------------------------------------------- #
library(mapdata)
library(ggplot2)
library(velociraptr)
library(rgdal) 
library(RCurl)
library(raster)
library(data.table)
library(rstudioapi)  # for auto-setting working dir

rm(list=ls())   # clear global environment 

# working dir is file dir
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# ----------------------------------------------------------------------------------------------- #

### FORMAT OF SCRIPT ### 
## 1) bring in the function
## 2) bring in paleogeographic shape files
## 3) run the for loop which loops through each time interval and generates an output dataframe
## 4) these dataframes are then saved into a list

# ----------------------------------------------------------------------------------------------- #

#... select the times for which to generate a data.table
ages <- seq(0,545,5)      # [Ma] ages for which to download paleogeography

## 1) READ IN THE LAND FRAC FUNCTION
source('GeogConstructDF_Fxn.R')

## 2) DOWNLOAD PALEOGEOGRAPHIES FOR EASY OFF-LINE CALCULATION 
# PaleoGeoList <- list()
# for(i in 1:length(ages)){
#   myAge <- ages[i]                  # age of the paleogeographic map in millions of years ago
#   pmap_poly <- downloadPaleogeography(Age = myAge)
#   PaleoGeoList[[i]] <- pmap_poly
# 
#   # ------
#   print(paste("Age =", ages[i], 'Ma', sep=' '))
# }

# -- SAVE PALEOGEOGRAPHY FILES -- # 
# saveRDS(PaleoGeoList, file='SpatialPolygons_5MyrSteps.RDS')
# ------------------------------- # 

## 2) READ IN THE PALEOGEOGRAPHY DATA
PaleoGeoList <- readRDS('SpatialPolygons_5MyrSteps.RDS')

## 3) RUN THE FOR LOOP
myres <- 2     # [deg] resolution that is used for the paleogeography
outlist <- list()     # hold the output dataframes
for(i in 1:length(PaleoGeoList)){
  #... keep track
  print(paste("Age =", ages[i], 'Ma', sep=' '))
  
  #... bring in the map data 
  pmap_poly <- PaleoGeoList[[i]]
  
  #... change to raster
  r <- raster(ncols=(360/myres), nrows=(180/myres))
  x <- rasterize(pmap_poly,r,getCover=TRUE,progress="text")
  pts <- rasterToPoints(x=x)
  myxyz <- as.data.table(pts) ; colnames(myxyz) <- c('lon', 'lat', 'landsea')
  #... 100 = land ; 0 = ocean
  
  mydf <- LandFrac_fun(xyzLat = myxyz$lat, xyzLon = myxyz$lon, xyzFrac = myxyz$landsea)
  
  # --------------------- #
  # ADD IT TO THE LIST!   # 
  # --------------------- #
  
  outlist[[i]] <- mydf
}

#... name the list elements
names(outlist) <- paste('Ma', ages, sep='_')

# ---------------- SAVE THE FINAL LIST ----------------- # 
# saveRDS(outlist, file='PaleoGeo_dataframes_2deg.RDS')
# ------------------------------------------------------ # 






