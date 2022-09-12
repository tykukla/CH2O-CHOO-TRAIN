# ------------------------------------------------- # 
#                   ALL ABOARD!!                    # 
#                                                   #
#     Script to process and plot CH2O-CHOO TRAIN    #
#                     results                       #  
#                                                   #
#                                                   #
# *************************************             #
# T Kukla (Colostate Univ. 2022)                    #
# KV Lau (Penn State Univ. 2022)                    #
# DE Ibarra (Brown Univ. 2022)                      #
# JKC Rugenstein (Colostate Univ. 2022)             #
# *************************************             # 
#                                                   #
# ------------------------------------------------- # 
library(data.table)  # preferred data format
library(rstudioapi)  # for setting working dir to file path
library(stringr)     # for finding correct filenames
library(cmocean)     # colorbars
library(ggplot2)     # plotting
library(ggthemes)    # plotting
library(ggtext)      # markdown formatting
library(colorspace)  # for clim ts plotting

# clear environment
rm(list=ls())

# set working dir to file path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# read in the functions to help processing
source('processing_functions.R')
# and the plotting parameters
source('parameters_plotting.R')

# find results
RESULTS.FOLDER <- 'Example_Results'
dat <- collect.dat(results.path = RESULTS.FOLDER)

# where to save figures
save.folder <- 'figures'  # script will make dir if it doesn't already exist 
permit.save <- TRUE       # [TRUE ; FALSE] whether to save the plots that are generated

## --------------------------------------------------------------------------------
## RUN EVERYTHING BELOW THIS TO SAVE ALL FIGS IN RESULTS.FOLDER 
## ... (or go step-by-step to visualize in RStudio)


# ---------------------------------------------- #
# --- PLOT INITIAL CLIMATOLOGY (BY LATITUDE) --- #
# ... pull out clim dat
dat.clim <- dat$initial_clim

## [0] LAND AREA
p.la.name <- 'initClim_LandArea.png'
p.la <- ggplot(dat.clim) + 
  geom_line(aes(x=sin(Lat*(pi/180)), y=Land_Area_m2, group=RUN.idx, color=RUN.idx), size=ln.size) +
  scale_y_continuous(name = expression("Land Area (m"^"2"*")")) +
  scale_x_continuous(breaks=lat.breaks.sin, expand=c(0,0), labels = lat.breaks, name = "Latitude (°N)") + 
  scale_color_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  theme_few() +
  theme(axis.title = element_text(size=ax.lab.size), 
        axis.text = element_text(size=ax.txt.size))

p.la
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.la, save.name = p.la.name)

## [1] TEMPERATURE
p.t.clim.name <- 'initClim_temperature.png'
p.t.clim <- ggplot(dat.clim) + 
  geom_line(aes(x=sin(Lat*(pi/180)), y=temp_C, group=RUN.idx, color=RUN.idx), size=ln.size) +
  scale_y_continuous(name = "Temperature (°C)") +
  scale_x_continuous(breaks=lat.breaks.sin, expand=c(0,0), labels = lat.breaks, name = "Latitude (°N)") + 
  scale_color_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  theme_few() +
  theme(axis.title = element_text(size=ax.lab.size), 
        axis.text = element_text(size=ax.txt.size))

p.t.clim
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.t.clim, save.name = p.t.clim.name)


## [2] DISCHARGE (BUDYKO RUNOFF * LAND AREA)
p.Q.clim.name <- 'initClim_budykoDischarge.png'
p.Q.clim <- ggplot(dat.clim) + 
  geom_line(aes(x=sin(Lat*(pi/180)), y=Discharge_m3_yr, group=RUN.idx, color=RUN.idx), size=ln.size) +
  scale_y_continuous(name = expression("Discharge (m"^"3"*"yr"^"-1"*")")) +
  scale_x_continuous(breaks=lat.breaks.sin, expand=c(0,0), labels = lat.breaks, name = "Latitude (°N)") + 
  scale_color_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  theme_few() +
  theme(axis.title = element_text(size=ax.lab.size), 
        axis.text = element_text(size=ax.txt.size))

p.Q.clim
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.Q.clim, save.name = p.Q.clim.name)


## [3] DISCHARGE (BUDYKO RUNOFF * LAND AREA)
p.Fsil.clim.name <- 'initClim_Fwsil.png'
p.Fsil.clim <- ggplot(dat.clim) + 
  geom_line(aes(x=sin(Lat*(pi/180)), y=Fwsil, group=RUN.idx, color=RUN.idx), size=ln.size) +
  scale_y_continuous(name = expression("Silicate weathering flux (mol C yr"^" -1"*")")) +
  scale_x_continuous(breaks=lat.breaks.sin, expand=c(0,0), labels = lat.breaks, name = "Latitude (°N)") + 
  scale_color_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  theme_few() +
  theme(axis.title = element_text(size=ax.lab.size), 
        axis.text = element_text(size=ax.txt.size))

p.Fsil.clim
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.Fsil.clim, save.name = p.Fsil.clim.name)




# ------------------------------- #
# --- PLOT TIMESERIES RESULTS --- #
# ... pull out clim dat
dat.ts <- dat$timeseries

## [0] VOLCANIC FORCING (C emissions)
p.Fvolc.ts.name <- 'timeseries_pCO2.png'
p.Fvolc.ts <- ggplot(dat.ts) + 
  geom_line(aes(x=time/1e3, y=Fvolc, group=iter, color=iter), size=ln.size) +
  scale_y_continuous(name="Volcanic C emissions (mol C yr <sup>-1</sup>)") +
  scale_x_continuous(name="Time (kyr)", expand=c(0,0)) +
  scale_color_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.Fvolc.ts
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.Fvolc.ts, save.name = p.Fvolc.ts.name)


## [1] ATMOSPHERIC CO2
p.co2.ts.name <- 'timeseries_pCO2.png'
p.co2.ts <- ggplot(dat.ts) + 
  geom_line(aes(x=time/1e3, y=pCO2, group=iter, color=iter), size=ln.size) +
  scale_y_continuous(name="Atmospheric pCO<sub>2</sub> (ppmv)") +
  scale_x_continuous(name="Time (kyr)", expand=c(0,0)) +
  scale_color_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.co2.ts
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.co2.ts, save.name = p.co2.ts.name)


## [2] GLOBAL T
p.t.ts.name <- 'timeseries_temperature.png'
p.t.ts <- ggplot(dat.ts) + 
  geom_line(aes(x=time/1e3, y=temp_C, group=iter, color=iter), size=ln.size) +
  scale_y_continuous(name="Temperature (°C)") +
  scale_x_continuous(name="Time (kyr)", expand=c(0,0)) +
  scale_color_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.t.ts
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.t.ts, save.name = p.t.ts.name)


## [3] GLOBAL MEAN RUNOFF (LAND AREA-WTD)
p.q.ts.name <- 'timeseries_runoff_landWt.png'
p.q.ts <- ggplot(dat.ts) + 
  geom_line(aes(x=time/1e3, y=runoff_budyko_LANDwt, group=iter, color=iter), size=ln.size) +
  scale_y_continuous(name="Global mean runoff (m yr <sup>-1</sup>)") +
  scale_x_continuous(name="Time (kyr)", expand=c(0,0)) +
  scale_color_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.q.ts
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.q.ts, save.name = p.q.ts.name)


## [4] SILICATE WEATHERING FLUX
p.Fsil.ts.name <- 'timeseries_Fwsil.png'
p.Fsil.ts <- ggplot(dat.ts) + 
  geom_line(aes(x=time/1e3, y=Fwsil, group=iter, color=iter), size=ln.size) +
  scale_y_continuous(name="Fw,sil (mol C yr <sup>-1</sup>)") +
  scale_x_continuous(name="Time (kyr)", expand=c(0,0)) +
  scale_color_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.Fsil.ts
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.Fsil.ts, save.name = p.Fsil.ts.name)


## [5] CARBON ISOTOPES OF DIC
p.d13.ts.name <- 'timeseries_d13C.png'
p.d13.ts <- ggplot(dat.ts) + 
  geom_line(aes(x=time/1e3, y=d13C, group=iter, color=iter), size=ln.size) +
  scale_y_continuous(name="*&delta;*<sup>13</sup>C (&permil; VPDB)") +
  scale_x_continuous(name="Time (kyr)", expand=c(0,0)) +
  scale_color_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.d13.ts
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.d13.ts, save.name = p.d13.ts.name)



# ------------------------------------------ #
# --- PLOT TIMESERIES OF ZONAL MEAN CLIM --- #
# ... pull out clim dat
dat.ts.clim <- dat$timeseries_clim
# calculate some anomalies
anoms.2.calc <- c('temp_C', 'Precip_m_yr', 'Discharge_m3_yr', 'Fwsil')
dat.ts.clim <- clim.anoms(zonal.clim = dat.ts.clim, anoms = anoms.2.calc)
# pick a run to investigate
this.run.idx <- 3

## [1] TEMPERATURE ANOMS
p.t.tsclim.name <- paste('ClimTS_anom_tempC-idx_', this.run.idx, '.png', sep='')
p.t.tsclim <- ggplot(dat.ts.clim[RUN.idx == this.run.idx]) + 
  geom_tile(aes(x=time/1e3, y=sin(Lat*(pi/180)), fill=temp_C_anom)) +
  scale_y_continuous(breaks=lat.breaks.sin, expand=c(0,0), labels = lat.breaks, name = "Latitude (°N)") +
  scale_x_continuous(name = "Time (kyr)", expand=c(0,0)) +
  scale_fill_continuous_divergingx(palette=pal.name.clim, mid = 0, rev=TRUE) +
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.t.tsclim
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.t.tsclim, save.name = p.t.tsclim.name)


## [2] PRECIP ANOMS
p.Pr.tsclim.name <- paste('ClimTS_anom_Precip-idx_', this.run.idx, '.png', sep='')
p.Pr.tsclim <- ggplot(dat.ts.clim[RUN.idx == this.run.idx]) + 
  geom_tile(aes(x=time/1e3, y=sin(Lat*(pi/180)), fill=Precip_m_yr_anom)) +
  scale_y_continuous(breaks=lat.breaks.sin, expand=c(0,0), labels = lat.breaks, name = "Latitude (°N)") +
  scale_x_continuous(name = "Time (kyr)", expand=c(0,0)) +
  scale_fill_continuous_divergingx(palette=pal.name.clim, mid = 0, rev=TRUE) +
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.Pr.tsclim
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.Pr.tsclim, save.name = p.Pr.tsclim.name)

## [3] DISCHARGE ANOMS
p.Q.tsclim.name <- paste('ClimTS_anom_Discharge-idx_', this.run.idx, '.png', sep='')
p.Q.tsclim <- ggplot(dat.ts.clim[RUN.idx == this.run.idx]) + 
  geom_tile(aes(x=time/1e3, y=sin(Lat*(pi/180)), fill=Discharge_m3_yr_anom)) +
  scale_y_continuous(breaks=lat.breaks.sin, expand=c(0,0), labels = lat.breaks, name = "Latitude (°N)") +
  scale_x_continuous(name = "Time (kyr)", expand=c(0,0)) +
  scale_fill_continuous_divergingx(palette=pal.name.clim, mid = 0, rev=TRUE) +
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.Q.tsclim
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.Q.tsclim, save.name = p.Q.tsclim.name)


## [4] FWSIL ANOMS
p.Fsil.tsclim.name <- paste('ClimTS_anom_Fwsil-idx_', this.run.idx, '.png', sep='')
p.Fsil.tsclim <- ggplot(dat.ts.clim[RUN.idx == this.run.idx]) + 
  geom_tile(aes(x=time/1e3, y=sin(Lat*(pi/180)), fill=Fwsil_anom)) +
  scale_y_continuous(breaks=lat.breaks.sin, expand=c(0,0), labels = lat.breaks, name = "Latitude (°N)") +
  scale_x_continuous(name = "Time (kyr)", expand=c(0,0)) +
  scale_fill_continuous_divergingx(palette=pal.name.clim, mid = 0, rev=TRUE) +
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.Fsil.tsclim
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.Fsil.tsclim, save.name = p.Fsil.tsclim.name)




# ------------------------------- #
# --- CLIMATE / WX CROSSPLOTS --- #
## [1] TEMPERATURE - CO2 CROSSPLOT
p.tco2.cp.name <- 'crossplot_TempCO2.png'
p.tco2.cp <- ggplot(dat.ts) + 
  geom_point(aes(x=pCO2, y=temp_C, group=iter, fill=iter, shape=climstate), color='black', size=5) +
  scale_y_continuous(name="Temperature (°C)") +
  scale_x_continuous(name="Atmospheric pCO<sub>2</sub> (ppmv)") +
  scale_fill_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  scale_shape_manual(values=c("no_ice" = 21, "south_pole" = 22, "north_pole" = 23, "both_poles" = 24)) + 
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.tco2.cp
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.tco2.cp, save.name = p.tco2.cp.name, cp.size = TRUE)


## [2] TEMPERATURE - FVOLC CROSSPLOT
p.tFvolc.cp.name <- 'crossplot_TempFvolc.png'
p.tFvolc.cp <- ggplot(dat.ts) + 
  geom_point(aes(x=Fvolc, y=temp_C, group=iter, fill=iter, shape=climstate), color='black', size=5) +
  scale_y_continuous(name="Temperature (°C)") +
  scale_x_continuous(name="Volcanic C emissions (mol C yr <sup>-1</sup>)") +
  scale_fill_cmocean(name=pal.name, direction = pal.dir, start=pal.start, end=pal.end, guide='none') +
  scale_shape_manual(values=c("no_ice" = 21, "south_pole" = 22, "north_pole" = 23, "both_poles" = 24)) + 
  theme_few() +
  theme(axis.text=element_text(size=ax.txt.size), 
        axis.title.y=element_markdown(size=ax.lab.size),
        axis.title.x=element_markdown(size=ax.lab.size))

p.tFvolc.cp
# save result (won't save if permit.save == FALSE)
plot.save(this.plot = p.tFvolc.cp, save.name = p.tFvolc.cp.name, cp.size = TRUE)






