# The CH2O-CHOO TRAIN 
[![DOI](https://zenodo.org/badge/535833456.svg)](https://zenodo.org/badge/latestdoi/535833456)

The **C**arbon-**H<sub>2</sub>O** **C**oupled Hydr**O**l**O**gical model with **T**errestrial **R**unoff **A**nd **IN**solation (**CH2O-CHOO TRAIN**) is an Earth System model for climate, weathering, and the long-term carbon cycle. 

The model computes zonal-mean climate and mineral weathering fields and feeds these results to a global model for the long-term ocean-atmosphere inorganic carbon inventory.

## Running the model
The model is coded in R ([R download](https://cran.r-project.org/)) and designed to run in RStudio ([RStudio download](https://www.rstudio.com/products/rstudio/download/)). If you are new to R, consider searching for a guide to downloading R and RStudio on Google or YouTube.  

The easiest way to get started with the model is to click the green **CODE** button at the top of the page, then click **Download ZIP**. Unzip (or 'extract') the downloaded file. Then open 'CODE_DISTRIBUTE.Rproj' or any specific '.R' file in RStudio. Running the code through the .Rproj file ensures that the required packages (from the 'renv' directory) are loaded using versions that work with the code. At time of writing, we are not aware of any conflicts with more recent package updates.

### Primary user scripts
The two most important model scripts are found in the 'RUN' directory.
1. **config.R**
  - Walks through user-defined settings for the simulation (including which input files to call)
2. **CH2O-CHOO-TRAIN_RUN.R**
  - Reads the config file and begins the simulation. 
  - Automatically creates an output directory and saves results.

There is a README file in most subdirectories to explain each directory's contents. These README's follow the same format as the 'SUBDIR-ROADMAP.txt' file in the main directory. The subdirectories contain scripts to build new input files.

If the CH2O-CHOO-TRAIN_RUN.R file is executed with no other changes, the code will output the results in the 'Example_Results' directory. The 'figures' folder in this directory is created from 'PLOT_RESULTS.R' in the 'PROCESS' directory.

## CH2O-CHOO TRAIN Publications
1. Kukla, T., Ibarra, D.E., Lau, K.V., Rugenstein, J.K.C. (Submitted). All aboard! Earth system investigations with the CH2O-CHOO TRAIN.
2. Kukla, T., Lau, K.V., Ibarra, D.E., Rugenstein, J.K.C. (Upcoming). Deterministic icehouse and greenhouse climates on million-year timescales.

### Key references
This model is based on previous energy balance climate model, weathering, and carbon cycling work, including: 
- Flannery, B. P. (1984). Energy balance models incorporating transport of thermal and latent energy. Journal of Atmospheric Sciences, 41(3), 414-421.
- Maher, K., & Chamberlain, C. P. (2014). Hydrologic regulation of chemical weathering and the geologic carbon cycle. science, 343(6178), 1502-1504.
- Roe, G. H., Feldl, N., Armour, K. C., Hwang, Y. T., & Frierson, D. M. (2015). The remote impacts of climate feedbacks on regional climate predictability. Nature Geoscience, 8(2), 135-139.
- Shields, G. A., & Mills, B. J. (2017). Tectonic controls on the long-term carbon isotope mass balance. Proceedings of the National Academy of Sciences, 114(17), 4318-4323.
- Siler, N., Roe, G. H., & Armour, K. C. (2018). Insights into the zonal-mean response of the hydrologic cycle to global warming from a diffusive energy balance model. Journal of Climate, 31(18), 7481-7493.
- Winnick, M. J., & Maher, K. (2018). Relationships between CO2, thermodynamic limits on silicate weathering, and the strength of the silicate weathering feedback. Earth and Planetary Science Letters, 485, 111-120.
- Caves Rugenstein, J. K., Ibarra, D. E., & von Blanckenburg, F. (2019). Neogene cooling driven by land surface reactivity rather than increased weathering fluxes. Nature, 571(7763), 99-102.
