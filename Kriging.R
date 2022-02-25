# Article reference: https://swilke-geoscience.net/post/2020-09-10-kriging_with_r/kriging/

# We will need some packages for (spatial) data processing
library(tidyverse) # wrangling tabular data and plotting
library(sf) # processing spatial vector data - the easy way
library(sp) # processing spatial vector data - the way gstat needs it
library(raster) # processing spatial raster data. !!!overwrites dplyr::select!!!

# Packages for geostatistics
library(gstat)   # The most popular R-Package for Kriging (imho)
library(automap) # Automatize some (or all) parts of the gstat-workflow 

# Some packages to make pretty plots
library(patchwork)
library(viridis)

# Finally, Excel functions to store data and plots
library(XLConnect)

# Download the data for this tutorial from Github!
# The data for this tutorial is derived from a dataset published here
# https://www.opengeodata.nrw.de/produkte/umwelt_klima/wasser/flurabstandskarte_1988/
# licence information:
# Datenlizenz Deutschland – Flurabstandskarte NRW 1988 – Version 2.0
# see http://www.govdata.de/dl-de/by-2-0 for more details

grd_100_df <- readr::read_csv(
  "testdata1.csv",
) %>% 
  
# grd_100_df <- readr::read_csv(
#  "https://raw.githubusercontent.com/Ignimbrit/exchange/master/data/2020/grid_100.csv",
# ) %>% 
dplyr::select(-licence)

# The projection is EPSG 25832

head(grd_100_df)

grd_100_rstr <- raster::rasterFromXYZ(
  grd_100_df, 
  res = c(100, 100), # resolution in meter (see crs)
  crs = "+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
)

png("Rplot.png")
plot(grd_100_rstr)
dummy <- dev.off()

set.seed(42) # forreproducibility

# Simulate 100 random observation wells 
wellobs <- slice_sample(grd_100_df, n = 100) # 100 is original sample
png("Rplot1.png")

ggplot(
  data = wellobs,
  mapping = aes(x = X, y = Y, color = Z)
) +
  geom_point(size = 3) + 
  scale_color_viridis(option = "B") +
  theme_classic()
dummy <- dev.off()

# Convert to {sf} because that is the best way to store spatial points
wellobs_sf <- st_as_sf(wellobs, coords = c("X", "Y"), crs = 25832) %>% 
  cbind(st_coordinates(.))

# We will discuss later, what Z~1 does actually mean in this context
v_emp_OK <- gstat::variogram(
  Z~1,
  as(wellobs_sf, "Spatial") # switch from {sf} to {sp}
)

png("Rplot2.png")
plot(v_emp_OK)
dummy <- dev.off()

# automap's autofitVariogram actually produces more info than we need.
# I will only keep the var_model part.
v_mod_OK <- automap::autofitVariogram(Z~1, as(wellobs_sf, "Spatial"))$var_model

# To inspect the automatic fit that was chosen for us we can use
# automap's excellent build in methods for base::plot

png("Rplot3.png")
plot(automap::autofitVariogram(Z~1, as(wellobs_sf, "Spatial")))
dummy <- dev.off()

# technically we already have a grid from the initial dataset, but as we are 
# still working under the pretense that our only available data are the 
# simulated observation wells, we will construct our grid from that object.

# Step 1: define a grid based on the bounding box of our observations
grd_100_sf <- wellobs_sf %>% 
  st_bbox() %>% 
  st_as_sfc() %>% 
  st_make_grid(
    cellsize = c(100, 100), # 100m pixel size
    what = "centers"
  ) %>%
  st_as_sf() %>%
  cbind(., st_coordinates(.))

# Step 2: making our grid work for gstat
grd_100_sp <- as(grd_100_sf, "Spatial") # converting to {sp} format
gridded(grd_100_sp) <- TRUE             # informing the object that it is a grid
grd_100_sp <- as(grd_100_sp, "SpatialPixels") # specifying what kind of grid

# That second step there is based on a discussion I found on Stackoverflow
# https://stackoverflow.com/questions/43436466/create-grid-in-r-for-kriging-in-gstat

# Ordinary Kriging
OK <- krige(
  Z~1,                       # Z is our variable and "~1" means "depends on mean"
  as(wellobs_sf, "Spatial"), # input data in {sp} format
  grd_100_sp,                # locations to interpolate at
  model = v_mod_OK           # the variogram model fitted above
)
## [using ordinary kriging]

# Simple Kriging
SK <- krige(
  Z~1,                       # Z still depends on mean
  beta = mean(grd_100_df$Z), # but this time we know the mean's value
  as(wellobs_sf, "Spatial"), # input data in {sp} format
  grd_100_sp,                # locations to interpolate at
  model = v_mod_OK           # the variogram model fitted above
)
## [using simple kriging]

# Universal Kriging
# Implementing this method is somewhat different.
# we no longer assume that Z is essentially depending on a single mean but
# rather on the position of the interpolation location within our target grid
UK <- krige(
  Z~coords.x1+coords.x2, # Think "Z~X+Y" but {sp} conversion alters variable naming
  as(wellobs_sf, "Spatial"), # input data in {sp} format (`X` --> `coords.x1`)
  grd_100_sp,                # locations to interpolate at
  model = autofitVariogram(  # we need an appropriate variogram fit
    Z~X+Y,                   # here we can keep "X+Y" - it's just how it is
    as(wellobs_sf, "Spatial")
  )$var_model
)
## [using universal kriging]

# I'll also add an inverse distance weighted model to provide a baseline
# for model evaluation
# Note how the only difference to Ordinary Kriging is the absence of a
# fitted variogram model
idwres <- idw(
  Z~1,                       # idw also depends on mean
  as(wellobs_sf, "Spatial"), # input data in {sp} format
  grd_100_sp,                # locations to interpolate at
) 
## [inverse distance weighted interpolation]

# A function to plot rasters
plot_my_gstat_output <- function(raster_object, object_name){
  
  df <- rasterToPoints(raster_object) %>% as_tibble()
  colnames(df) <- c("X", "Y", "Z")
  
  ggplot(df, aes(x = X, y = Y, fill = Z)) +
    geom_raster() +
    ggtitle(label = object_name) +
    scale_fill_viridis(option = "B", limits = c(50, 100)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
}

p_orig <- plot_my_gstat_output(grd_100_rstr, "Original Raster")
p_idw <- plot_my_gstat_output(raster(idwres), "IDW")
p_SK <- plot_my_gstat_output(raster(SK), "Simple Kriging")
p_OK <- plot_my_gstat_output(raster(OK), "Ordinary Kriging")
p_UK <- plot_my_gstat_output(raster(UK), "Universal Kriging")

# I also want to display sampling locations
p_wellobs <- ggplot(
  data = wellobs,
  mapping = aes(x = X, y = Y, color = Z)
) +
  geom_point(size = 3) + 
  scale_color_viridis(option = "B",  limits = c(50, 100)) +
  ggtitle(label = "Observation Wells Sampled") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

## [Show Differences between Kriging Methods]
map(list(SK, OK, UK), raster) %>% 
  map(summary) %>%
  do.call("cbind", .) %>% 
  as.data.frame() %>% 
  setNames(c("SK", "OK", "UK"))

compMap <- map(list(SK, OK, UK), raster) %>% 
  map(summary) %>%
  do.call("cbind", .) %>% 
  as.data.frame() %>% 
  setNames(c("SK", "OK", "UK"))
# This works because of library(patchwork)
png("Rplot4.png")
(p_orig + p_wellobs + p_idw) / 
  (p_SK + p_OK + p_UK) + 
  plot_layout(guides = 'collect')
dummy <- dev.off()


# Save data and images in Excel
fn <- "RTest.xlsx"
if (file.exists(fn)) {
  file.remove(fn) 
}

wb <- loadWorkbook(fn, create=TRUE)
# Original Data
createSheet(wb,"Original Data")
writeWorksheet(wb,grd_100_df,sheet="Original Data",startRow=1,startCol=1)

# Plot 0
createSheet(wb, "Plots0")
createName(wb, name = "graph", formula = "Plots0!$A$1", overwrite=TRUE)
addImage(wb, filename = "Rplot.png", name="graph", originalSize=TRUE)

# Plot 1
createSheet(wb, "Plots1")
createName(wb, name = "graph", formula = "Plots1!$A$1", overwrite=TRUE)
addImage(wb, filename = "Rplot1.png", name="graph", originalSize=TRUE)

# Plot 2
createSheet(wb, "Plots2")
createName(wb, name = "graph", formula = "Plots2!$A$1", overwrite=TRUE)
addImage(wb, filename = "Rplot2.png", name="graph", originalSize=TRUE)

# Plot 3
createSheet(wb, "Plots3")
createName(wb, name = "graph", formula = "Plots3!$A$1", overwrite=TRUE)
addImage(wb, filename = "Rplot3.png", name="graph", originalSize=TRUE)

# Plot 4
createSheet(wb, "Plots4")
createName(wb, name = "graph", formula = "Plots4!$A$1", overwrite=TRUE)
addImage(wb, filename = "Rplot4.png", name="graph", originalSize=TRUE)

# Plot 4 Comparison
createSheet(wb, "Plots4")
writeWorksheet(wb,compMap,sheet="Plots4", startRow=1,startCol=12)

# Close workbook
dummy <- file.remove("Rplot.png")
dummy <- file.remove("Rplot1.png")
dummy <- file.remove("Rplot2.png")
dummy <- file.remove("Rplot3.png")
dummy <- file.remove("Rplot4.png")

saveWorkbook(wb)
