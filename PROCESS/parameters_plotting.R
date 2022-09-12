# --------------------------------------- # 
# -- PLOTTING PARAMETERS ---------------- #
# --------------------------------------- # 

# Set up x-ax transform
lat.breaks <- c(-90, -60, -30, 0, 30, 60, 90)    # which latitudes to label on x-ax
# lat.breaks <- c(-50, -30, -15, 0, 15, 30, 50)
lat.breaks.sin <- sin(lat.breaks * (pi/180))     # transform 
# set up colors
pal.name <- 'deep'       # name of palette (cmocean: https://matplotlib.org/cmocean/)
pal.dir <- -1            # direction of palette (forward = 1; backward = -1)
pal.start <- 0.1         # [0,1] start color value 
pal.end <- 0.9           # [0,1] end color value 
# other aes
ln.size <- 1.7           # line size
ax.lab.size <- 18        # axis label text size
ax.txt.size <- 14        # axis text text size
# xax lims [CURRENTLY NOT USED]
t.min <- 0/1e3 ; t.max <- 1e6/1e3   # [kyr] age limits for timeseries plots

# clim timeseries colors
pal.name.clim <- 'RdBu'  # name of palette (colorspace, see: https://stackoverflow.com/questions/58718527/setting-midpoint-for-continuous-diverging-color-scale-on-a-heatmap)

# save aes
# timeseries / lat sizes
save.width <- 18          # [cm] width of saved plot 
save.height <- 9          # [cm] height of saved plot 
# crossplot sizes
save.width.cp <- 15           # [cm] width of saved plot 
save.height.cp <- 12          # [cm] height of saved plot

