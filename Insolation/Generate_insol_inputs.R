# ----------------------------------------------- #
# Pull out paleo-insolation data to make energy   #
# balance model input files for moist energy      #
# balance modeling                                #
# ---                                             #
# T Kukla (Stanford Univ. 2019)                   #
# ----------------------------------------------- #
library(ggplot2)
library(dplyr)
library(palinsol)
library(viridis)
library(data.table)
library(ggthemes)

rm(list=ls())

# find parent dir
insol.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)  # defaults to this file's path
parent.path.split <- strsplit(insol.dir, '/')[[1]]
parent.path.idx <- which(parent.path.split == 'CODE_DISTRIBUTE')
parent.path <- paste0(parent.path.split[c(1:parent.path.idx)], collapse='/')

# where to save result
save.2.run <- TRUE
if(save.2.run == TRUE){
  save.here <- 'Initialization/Insolation'  # insol files can be called from this folder to run
} else{
  save.here <- insol.dir                    # otherwise just save it in the insol folder
}

# ************************************* #
# Define the relevant days of year for  
# computation (normal [non-leap] year)
# (after: https://cals.arizona.edu/azmet/julian.html)
Jan.start <- 1; Jan.end <- 31
Feb.start <- 32; Feb.end <- 59
Mar.start <- 60; Mar.end <- 90
Apr.start <- 91; Apr.end <- 120
May.start <- 121; May.end <- 151
Jun.start <- 152; Jun.end <- 181
Jul.start <- 182; Jul.end <- 212
Aug.start <- 213; Aug.end <- 243
Sep.start <- 244; Sep.end <- 273
Oct.start <- 274; Oct.end <- 304
Nov.start <- 305; Nov.end <- 334
Dec.start <- 335; Dec.end <- 365

# ************************************* #
# Set the time-slices (years) over 
# which to calculate global insolation 
yrMin <- -10e6
yrMax <- -1
yrStep <- 2500
myYrs <- seq(yrMin, yrMax, yrStep)

# ************************************* #
# Set output grid 
myLats <- seq(-90, 89, 2)



# ****************************************************** #
# DEFINE THE MONTH-INTEGRATION FOR DATA COMPILATION      #
# EXAMPLES:                                              #
# [1] DJF: Month.start=Dec.start; Month.end=Feb.end      #
# [2] JJA: Month.start=Jun.start; Month.end=Aug.end      #
# [3] NDJFMAM: Month.start=Nov.start; Month.end=May.end  #
# [4] ANN: Month.start=Jul.start; Month.end=Jun.end      #
Month.start <- Jul.start                                 #
Month.end <- Jun.end                                     #
# ---- TIMESLICE NAME (for output file)                  #
tsl <- 'ANN'                                             #
# ****************************************************** #



# ====================== RUN THE FOR-LOOP ========================== #
for(i in 1:length(myYrs)){   # loop through each year
  #... extract the Laskar 2004 solution
  la_sol <- la04(t=myYrs[i])
  
  #... calculate the true solar longitude bounds
  TrueLon_min <- day2l(orbit=la_sol, day=Month.start)
  TrueLon_max <- day2l(orbit=la_sol, day=Month.end)
  
  #... vector to hold latitude outputs 
  Insol_by_lat <- vector()
  
  for(j in 1:length(myLats)){    # loop through each latitude
    #... get latitude in radians
    lat_rad <- myLats[j] * (pi/180)
    #... calculate insolation at this latitude and time
    thisInsol <- Insol_l1l2(orbit=la_sol, l1=TrueLon_min, l2=TrueLon_max, lat=lat_rad, avg=TRUE, ell=TRUE)
    
    #... save result
    Insol_by_lat[j] <- thisInsol
    
  }
  #... create the output dataframe if necessary
  if(i==1){
    outdf <- as.data.table(cbind(myLats, Insol_by_lat, myYrs[i])) ; colnames(outdf) <- c('lat', 'insol', 'yr')
  } else{
    ndf <- as.data.table(cbind(myLats, Insol_by_lat, myYrs[i])) ; colnames(ndf) <- c('lat', 'insol', 'yr')
    outdf <- as.data.table(rbind(outdf, ndf)) ; colnames(outdf) <- c('lat', 'insol', 'yr')
  }
  
  #... track progress
  print(paste('Calculating timeslice...', i, 'of', length(myYrs)), sep=' ')
}

# get positive years (easier for integrating with CH2O-CHO MEBM)
outdf$yr_positive <- -min(outdf$yr) + outdf$yr # shift everything up by most negative value 


# ******************************************************************************* # 
# -- UNCOMMENT TO MAKE AN IDENTICAL VERSION WITH CONSTANT INSOL AT YEAR 0 AS CTRL #
# outdf.initYr <- outdf[which(outdf$yr_positive==min(outdf$yr_positive))]  # line used to select initial year in INSOL script
# outdf.initYr.FILL <- as.data.table(cbind(outdf$yr_positive, outdf.initYr))
# outdf.initYr.FILL$yr_positive <- outdf.initYr.FILL$V1
# outdf.initYr.FILL[ , V1:=NULL]
# # check that it worked
# ggplot(outdf.initYr.FILL[yr_positive == 80000]) +
#   geom_line(aes(x=lat, y=insol), size=3) +
#   geom_line(data=outdf.initYr.FILL[yr_positive == 1147500],
#             aes(x=lat, y=insol), linetype='dashed', color='red') +
#   geom_line(data=outdf.initYr.FILL[yr_positive == 335000],
#             aes(x=lat, y=insol), linetype='dotdash', color='blue', size=2, alpha=0.6) +
#   lims(y=c(150, 450))
# ## ******** SAVE THE CONSTANT FILE ************ #
# saveName <- 'insol_10Myr-2500steps_ANN-CONSTANT_T0.RDS'
# saveRDS(outdf.initYr.FILL, paste(save.here, saveName, sep='/'))
# ## ******************************************* #


# ******** SAVE YOUR OUTPUT FILE ************ #
## -- UNCOMMENT THE LAST LINE TO SAVE
## -- (commented out to prevent accidental overwrite)
saveName <- 'insol_10Myr-2500steps_ANN.RDS'
# saveRDS(outdf, paste(save.here, saveName, sep='/'))     
# ******************************************* #



## -- VISUALIZE RESULTS ----------------------------------------------------------------------------
#... take a quick look at some of the output
ggplot(outdf) + 
  geom_line(aes(x=lat, y=insol, group=yr, color=yr)) + 
  scale_color_viridis_c() + 
  theme_linedraw()

#... raster of all lats
ggplot(outdf) + 
  geom_raster(aes(x=lat, y=yr_positive, fill=insol)) +
  scale_fill_viridis_c() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_few()

#... zoom in on a specific lat 
thislat <- 70 ; y.range <- 15
ggplot(outdf[lat == thislat]) + 
  geom_line(aes(x=yr, y=insol)) + 
  scale_y_continuous(limits=c(min(outdf[lat == thislat]$insol), min(outdf[lat == thislat]$insol) + y.range)) +
  theme_few()

## -- scratch calculations
# analyze magnitude of variability by latitude
these.lats <- unique(outdf$lat)
out.range <- out.factorchange <- out.lat <- vector()
for(j in 1:length(these.lats)){
  tempdf <- outdf[lat == these.lats[j]]
  # save results
  out.lat[j] <- these.lats[j]
  out.range[j] <- range(tempdf$insol)[2] - range(tempdf$insol)[1]
  out.factorchange[j] <- out.range[j] / tempdf[yr==max(tempdf$yr)]$insol[1]
}
# bring together
df.sum <- as.data.table(cbind(out.lat, out.range, out.factorchange))
colnames(df.sum) <- c('lat', 'range', 'factorchange')
df.sum$percentVar <- df.sum$factorchange * 100


ggplot(df.sum) + 
  geom_line(aes(x=lat, y=percentVar)) + 
  theme_few()

ggplot(df.sum) + 
  geom_line(aes(x=lat, y=range)) + 
  theme_few()



