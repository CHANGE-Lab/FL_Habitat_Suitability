# working directories
setwd("Z:/Courtney/Stuart_MSc_Ch1/")
temp_wd = "Z:/Courtney/Stuart_MSc_Ch1/Temp/" # temporary files
source_wd = "Z:/Courtney/Stuart_MSc_Ch1/Source_Data/" # source data
dem_wd = "Z:/Courtney/Stuart_MSc_Ch1/Source_Data/Job606638_ncei_nintharcsec_dem/" # DEMs
csv_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Data/" # for writing tabular data
train_wd = "Z:/Courtney/Stuart_MSc_Ch1/Modeling_Data/Training/" # for training area data
test_wd = "Z:/Courtney/Stuart_MSc_Ch1/Modeling_Data/Testing/" # for testing area data

# libraries
library(easypackages)
libraries("rgdal", "gdalUtils", "raster", "sp", "sf", "tmap", "lwgeom", "rgeos",
          "cleangeo", "tidyverse", "stars", "fasterize", "PNWColors", "spex", "igraph", 
          "spatialEco", "tibble", "ncf", "spdep", "gstat", "geoR", "readxl",
          "tidyr", "dplyr")

# change where large temp rasters are saved
rasterOptions(tmpdir = "Z:/Courtney/Stuart_MSc_Ch1/Temp/")

# save PROJ.4 string for standard projection (ESPG:26958 NAD 83/Florida East) and source LAT/LON data (EPSG:4326 WGS 84/World Geodetic System 1984)
my_crs = CRS("+init=epsg:26958")
gcs = CRS("+init=epsg:4326")

# necessary data 
# training and testing areas, and full domain
train = st_read(paste0(train_wd, "GIS/Training_Area.shp")) # training area
test = st_read(paste0(test_wd, "GIS/Testing_Area.shp")) # testing area 

#### WATER QUALITY DATA ####
# full study domain for selecting water quality stations
domain = rbind(train, test)

# water quality data that is to be interpolated
wq = read_excel(paste0(source_wd, "Water_Conditions/WQFloridaKeys&Shelf (ppm) UPDATED 6-6-2020.xlsx"), 
                sheet = "Data in ppm") # from The SERC Water Quality Monitoring Network (http://serc.fiu.edu/wqmnetwork/)
head(wq, 5)

# now data prep and cleaning
wq = wq %>%
  mutate(DATE = as.Date(wq$DATE, origin = "1899-12-30")) %>% # this is the origin for excel sheets according to http://support.microsoft.com/kb/214330
  separate(DATE, into = c("YEAR", "MONTH", "DAY"), sep = "-", remove = F) %>% # separate out time data
  filter(YEAR %in% c("2014", "2015", "2016", "2017", "2018")) %>% # keeping only years of interest
  rename(TEMP_B = `TEMP-B`, SAL_B = `SAL-B`, DO_B = `DO-B`) %>% # rename variables because "-" causes problems
  st_as_sf(., coords = c("LONDEC", "LATDEC"), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) %>% # re-project to standard CRS
  mutate(LON_M = sf::st_coordinates(.)[,1], # save LON_M and LAT_M columns
         LAT_M = sf::st_coordinates(.)[,2]) %>%
  st_intersection(., st_make_valid(domain)) %>% # clip data to study domain 
  add_column(SEASON = NA) %>% # four seasons based on quarterly sampling 
  mutate(SEASON = cut(as.numeric(.$MONTH), breaks = c(0,3,6,9,Inf), labels = c("winter", "spring", "summer", "fall"))) %>%
  select(BASIN, STATION, SITE, YEAR, MONTH, SEASON, LON_M, LAT_M, TEMP_B, SAL_B, DO_B) 

#### WINTER DATA ####

# here are the stations that were sampled every winter for the full five year period
winter_sites = wq %>%
  select(STATION, SITE, YEAR, MONTH, SEASON, LON_M, LAT_M, TEMP_B, SAL_B, DO_B) %>%
  filter(SEASON == "winter", !is.na(TEMP_B), !is.na(SAL_B), !is.na(DO_B)) %>%
  group_by(SITE, LON_M, LAT_M) %>%
  count() %>%
  filter(n == 5) %>%
  ungroup()

# restricting the wq data to only those sites in the winter_sites data frame and calculating summary stats
winter_wq = wq %>%
  filter(SITE %in% winter_sites$SITE, SEASON == "winter") %>%
  group_by(SITE, LON_M, LAT_M) %>%
  mutate(MEAN_WIN_TEMP_B = mean(TEMP_B), SD_WIN_TEMP_B = sd(TEMP_B),
         MEAN_WIN_SAL_B = mean(SAL_B), SD_WIN_SAL_B = sd(SAL_B),
         MEAN_WIN_DO_B = mean(DO_B), SD_WIN_DO_B = sd(DO_B),
         N = as.integer(5)) %>%
  select(everything()) %>%
  distinct(BASIN, STATION, SITE, LON_M, LAT_M, N, MEAN_WIN_TEMP_B, SD_WIN_TEMP_B,
           MEAN_WIN_SAL_B, SD_WIN_SAL_B, MEAN_WIN_DO_B, SD_WIN_DO_B)

# for variogram modeling and spatial estimation we need SPDF, so convert sf object to SPDF
winter_wq_sp = winter_wq %>% st_drop_geometry()
coordinates(winter_wq_sp) = ~ LON_M + LAT_M
proj4string(winter_wq_sp) = proj4string(my_crs) # use my_crs CRS object to define proj4 string of the SPDF
summary(winter_wq_sp) # not much variation going on

# first taking a glimpse at correlograms and the Moran's I values for the data to ensure that there is spatial dependence
winter_coords = cbind(winter_wq_sp$LON_M, winter_wq_sp$LAT_M) # calculate a distance matrix
colnames(winter_coords) = c("LON_M", "LAT_M")
winter_distmat = as.matrix(dist(winter_coords))
winter_maxdist = 2/3 * max(winter_distmat) # maximum distance to consider in correlogram/variogram

# spline correlograms with 95% pointwise bootstrap CIs
w_temp_corr = spline.correlog(x = winter_wq_sp$LON_M, y = winter_wq_sp$LAT_M, z = winter_wq_sp$MEAN_WIN_TEMP_B, xmax = winter_maxdist, resamp = 100, type = "boot")
plot(w_temp_corr) # with 95% CIs from bootstrapping
w_sal_corr = spline.correlog(x = winter_wq_sp$LON_M, y = winter_wq_sp$LAT_M, z = winter_wq_sp$MEAN_WIN_SAL_B, xmax = winter_maxdist, resamp = 100, type = "boot")
plot(w_sal_corr) 
w_do_corr = spline.correlog(x = winter_wq_sp$LON_M, y = winter_wq_sp$LAT_M, z = winter_wq_sp$MEAN_WIN_DO_B, xmax = winter_maxdist, resamp = 100, type = "boot")
plot(w_do_corr)

# neighborhood list (neighbors within 20 km distance)
winter_neigh = dnearneigh(x = winter_coords, d1 = 0, d2 = 20000, longlat = F)
plot(winter_neigh, coordinates(winter_coords))
winter_wts = nb2listw(neighbours = winter_neigh, style = "W", zero.policy = T) # weights matrix for calculating Moran's I

# Moran's I with normal approximations
moran.test(winter_wq_sp$MEAN_WIN_TEMP_B, listw = winter_wts, randomisation = F, zero.policy = T)  # p < 0.01 sig. spatial dependence
moran.test(winter_wq_sp$MEAN_WIN_SAL_B, listw = winter_wts, randomisation = F, zero.policy = T)  # p < 0.01 sig. spatial dependence
moran.test(winter_wq_sp$MEAN_WIN_DO_B, listw = winter_wts, randomisation = F, zero.policy = T)  # p < 0.01 sig. spatial dependence

# Moran's I with Monte Carlo permutations, does everything match with the normal approximations?
moran.mc(winter_wq_sp$MEAN_WIN_TEMP_B, listw = winter_wts, nsim = 99, zero.policy = T) # p = 0.01 sig. spatial dependence
moran.mc(winter_wq_sp$MEAN_WIN_SAL_B, listw = winter_wts, nsim = 99, zero.policy = T) # p = 0.01 sig. spatial dependence
moran.mc(winter_wq_sp$MEAN_WIN_DO_B, listw = winter_wts, nsim = 99, zero.policy = T) # p = 0.01 sig. spatial dependence

# according to Moran tests with both normal and Monte Carlo approximations, there is significant spatial dependence in the winter temp, sal, and DO data.
# now  move on to variogram modeling to capture the spatial structure...
require(gstat)

w_temp_evgm = variogram(MEAN_WIN_TEMP_B ~ 1, winter_wq_sp, cutoff = winter_maxdist) # empirical variogram
plot(w_temp_evgm, xlab = "distance (m)", pch = 19)
w_temp_fvgm = fit.variogram(w_temp_evgm, vgm(psill = 0.25, model = "Sph", range = 50000, nugget = 0.1)) # fit variogram
w_temp_svgm_plot = plot(w_temp_evgm, model = w_temp_fvgm, xlab = "distance (m)", pch = 19)
w_temp_svgm_plot
print(w_temp_fvgm)

w_sal_evgm = variogram(MEAN_WIN_SAL_B ~ 1, winter_wq_sp, cutoff = winter_maxdist)
plot(w_sal_evgm, xlab = "distance (m)", pch = 19)
w_sal_fvgm = fit.variogram(w_sal_evgm, vgm(psill = 0.25, model = "Sph", range = 50000, nugget = 0.1))
w_sal_svgm_plot = plot(w_sal_evgm, model = w_sal_fvgm, xlab = "distance (m)", pch = 19)
w_sal_svgm_plot
print(w_sal_fvgm)

w_do_evgm = variogram(MEAN_WIN_DO_B ~ 1, winter_wq_sp, cutoff = winter_maxdist)
plot(w_do_evgm, xlab = "distance (m)", pch = 19)
w_do_fvgm = fit.variogram(w_do_evgm, vgm(psill = 0.005, model = "Sph", range = 50000, nugget = 0.005))
w_do_svgm_plot = plot(w_do_evgm, model = w_do_fvgm, xlab = "distance (m)", pch = 19)
w_do_svgm_plot
print(w_do_fvgm)

#### SUMMER DATA ####
# here are the stations that were sampled every summer for the full five year period
summer_sites = wq %>%
  select(STATION, SITE, YEAR, MONTH, SEASON, LON_M, LAT_M, TEMP_B, SAL_B, DO_B) %>%
  filter(SEASON == "summer", !is.na(TEMP_B), !is.na(SAL_B), !is.na(DO_B)) %>%
  group_by(SITE, LON_M, LAT_M) %>%
  count() %>%
  filter(n == 5) %>%
  ungroup()

# restricting the wq data frame to only those sites in the summer_sites data frame and calculating summary stats
summer_wq = wq %>%
  filter(SITE %in% summer_sites$SITE, SEASON == "summer") %>%
  group_by(SITE, LON_M, LAT_M) %>%
  mutate(MEAN_SUM_TEMP_B = mean(TEMP_B), SD_SUM_TEMP_B = sd(TEMP_B),
         MEAN_SUM_SAL_B = mean(SAL_B), SD_SUM_SAL_B = sd(SAL_B),
         MEAN_SUM_DO_B = mean(DO_B), SD_SUM_DO_B = sd(DO_B),
         N = as.integer(5)) %>%
  select(everything()) %>%
  distinct(BASIN, STATION, SITE, LON_M, LAT_M, N, MEAN_SUM_TEMP_B, SD_SUM_TEMP_B,
           MEAN_SUM_SAL_B, SD_SUM_SAL_B, MEAN_SUM_DO_B, SD_SUM_DO_B)

# for variogram modeling and spatial estimation we need SPDF, so convert sf object to SPDF
summer_wq_sp = summer_wq %>% st_drop_geometry()
coordinates(summer_wq_sp) = ~ LON_M + LAT_M
proj4string(summer_wq_sp) = proj4string(my_crs) # use my_crs CRS object to define proj4 string of the SPDF
summary(summer_wq_sp)

# first taking a glimpse at correlograms and the Moran's I values for the data to ensure that there is spatial dependence
summer_coords = cbind(summer_wq_sp$LON_M, summer_wq_sp$LAT_M) #calculate a distance matrix
colnames(summer_coords) = c("LON_M", "LAT_M")
summer_distmat = as.matrix(dist(summer_coords))
summer_maxdist = 2/3 * max(summer_distmat) #maximum distance to consider in correlogram/variogram

# spline correlograms with 95% pointwise bootstrap CIs
s_temp_corr = spline.correlog(x = summer_wq_sp$LON_M, y = summer_wq_sp$LAT_M, z = summer_wq_sp$MEAN_SUM_TEMP_B, xmax = summer_maxdist, resamp = 100, type = "boot")
plot(s_temp_corr) # with 95% CIs from bootstrapping
s_sal_corr = spline.correlog(x = summer_wq_sp$LON_M, y = summer_wq_sp$LAT_M, z = summer_wq_sp$MEAN_SUM_SAL_B, xmax = summer_maxdist, resamp = 100, type = "boot")
plot(s_sal_corr) 
s_do_corr = spline.correlog(x = summer_wq_sp$LON_M, y = summer_wq_sp$LAT_M, z = summer_wq_sp$MEAN_SUM_DO_B, xmax = summer_maxdist, resamp = 100, type = "boot")
plot(s_do_corr)

# neighborhood list (neighbors within 20 km distance)
summer_neigh = dnearneigh(x = summer_coords, d1 = 0, d2 = 20000, longlat = F)
plot(summer_neigh, coordinates(summer_coords))
summer_wts = nb2listw(neighbours = summer_neigh, style = "W", zero.policy = T) # weights matrix for calculating Moran's I


# Moran's I with normal approximations
moran.test(summer_wq_sp$MEAN_SUM_TEMP_B, listw = summer_wts, randomisation = F, zero.policy = T)  # est. Moran's I stat = 0.085335893, p = 0.1228 minimal spatial dependence
moran.test(summer_wq_sp$MEAN_SUM_SAL_B, listw = summer_wts, randomisation = F, zero.policy = T) # est. Moran's I stat = 0.564442066, p < 0.01 sig. spatial dependence
moran.test(summer_wq_sp$MEAN_SUM_DO_B, listw = summer_wts, randomisation = F, zero.policy = T) # est. Moran's I stat = 0.210497631, p < 0.01 sig. spatial dependence


# Moran's I with Monte Carlo permutations, does everything match with the normal approximations?
moran.mc(summer_wq_sp$MEAN_SUM_TEMP_B, listw = summer_wts, nsim = 99, zero.policy = T)  # est. Moran's I stat = 0.085336, p = 0.11 minimal spatial dependence
moran.mc(summer_wq_sp$MEAN_SUM_SAL_B, listw = summer_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.56444, p = 0.01 sig. spatial dependence
moran.mc(summer_wq_sp$MEAN_SUM_DO_B, listw = summer_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.2105, p = 0.01 sig. spatial dependence


# according to Moran tests with both normal and Monte Carlo estimations, 
# there is some (minimal) spatial dependence in the temp data,
# and significant spatial dependence in the sal and DO data.
# now we can move on to variogram modeling to capture the spatial structure...
require(gstat)

s_temp_evgm = variogram(MEAN_SUM_TEMP_B ~ 1, summer_wq_sp, cutoff = summer_maxdist) # empirical variogram
plot(s_temp_evgm, xlab = "distance (m)", pch = 19)
s_temp_fvgm = fit.variogram(s_temp_evgm, vgm(psill = 0.35, model = "Sph", range = 25000, nugget = 0.25))
s_temp_svgm_plot = plot(s_temp_evgm, model = s_temp_fvgm, xlab = "distance (m)", pch = 19)
s_temp_svgm_plot
print(s_temp_fvgm)

s_sal_evgm = variogram(MEAN_SUM_SAL_B ~ 1, summer_wq_sp, cutoff = summer_maxdist)
plot(s_sal_evgm, xlab = "distance (m)", pch = 19)
s_sal_fvgm = fit.variogram(s_sal_evgm, vgm(psill = 0.25, model = "Sph", range = 55000, nugget = 0.1))
s_sal_svgm_plot = plot(s_sal_evgm, model = s_sal_fvgm, xlab = "distance (m)", pch = 19)
s_sal_svgm_plot
print(s_sal_fvgm)

s_do_evgm = variogram(MEAN_SUM_DO_B ~ 1, summer_wq_sp, cutoff = summer_maxdist)
plot(s_do_evgm, xlab = "distance (m)", pch = 19)
s_do_fvgm = fit.variogram(s_do_evgm, vgm(psill = 0.002, model = "Sph", range = 50000, nugget = 0.003))
s_do_svgm_plot = plot(s_do_evgm, model = s_do_fvgm, xlab = "distance (m)", pch = 19)
s_do_svgm_plot
print(s_do_fvgm)

#### KRIGING: TRAINING AREA ####
# habitat raster as a guide for kriging
habitat_train = raster(paste0(train_wd, "Environmental/Habitat.asc")) 
crs(habitat_train) = my_crs

# prediction grid
train_grid = raster(ncol = ncol(habitat_train), nrow = nrow(habitat_train), xmn = xmin(habitat_train), 
                    xmx = xmax(habitat_train), ymn = ymin(habitat_train), ymx = ymax(habitat_train))
train_grid = as(train_grid, "SpatialPixels") # convert to spatial pixels object
proj4string(train_grid) = proj4string(my_crs) # assign projection 

# ordinary kriging (running in parallel) using empirical semivariograms calculated above
# starting with winter conditions in the training area

#### winter temperature ####
#Calculate the number of cores
library(doParallel)
no_cores = detectCores() - 2

# Initiate cluster 
cl = makeCluster(no_cores)

train_parts = split(x = 1:length(train_grid), f = 1:no_cores)

clusterExport(cl = cl, varlist = c("winter_wq_sp", "train_grid", "train_parts", "w_temp_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_temp_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_TEMP_B ~ 1, locations = winter_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = w_temp_fvgm))
stopCluster(cl)

# merge all the predictions
w_temp_train_merge = maptools::spRbind(w_temp_train_par[[1]], w_temp_train_par[[2]])
for (j in 3:length(w_temp_train_par)) {
  w_temp_train_merge = maptools::spRbind(w_temp_train_merge, w_temp_train_par[[j]])
}
w_temp_train_merge = SpatialPixelsDataFrame(points = w_temp_train_merge, data = w_temp_train_merge@data)

# save data
summary(w_temp_train_merge)
writeGDAL(w_temp_train_merge["var1.pred"], fname = paste0(temp_wd, "mean_win_temp_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_temp_train = raster(paste0(temp_wd, "mean_win_temp_train.tif"))
mean_win_temp_train = writeRaster(raster::mask(raster::crop(mean_win_temp_train, habitat_train), habitat_train),
                                  file = file.path(train_wd, "Environmental/Mean_Win_Temp.asc"), format = "ascii",
                                  overwrite = T)

rm(list = c("w_temp_train_merge", "w_temp_train_par", "mean_win_temp_train", "cl"))
showConnections()

#### winter salinity ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "train_grid", "train_parts", "w_sal_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_sal_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_SAL_B ~ 1, locations = winter_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = w_sal_fvgm))
stopCluster(cl)

# merge all the predictions
w_sal_train_merge = maptools::spRbind(w_sal_train_par[[1]], w_sal_train_par[[2]])
for (j in 3:length(w_sal_train_par)) {
  w_sal_train_merge = maptools::spRbind(w_sal_train_merge, w_sal_train_par[[j]])
}
w_sal_train_merge = SpatialPixelsDataFrame(points = w_sal_train_merge, data = w_sal_train_merge@data)

# save data
summary(w_sal_train_merge)
writeGDAL(w_sal_train_merge["var1.pred"], fname = paste0(temp_wd, "mean_win_sal_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_sal_train = raster(paste0(temp_wd, "mean_win_sal_train.tif"))
mean_win_sal_train = writeRaster(raster::mask(raster::crop(mean_win_sal_train, habitat_train), habitat_train),
                                  file = file.path(train_wd, "Environmental/Mean_Win_Sal.asc"), format = "ascii",
                                  overwrite = T)

rm(list = c("w_sal_train_merge", "w_sal_train_par", "mean_win_sal_train", "cl"))
showConnections()

#### winter DO ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "train_grid", "train_parts", "w_do_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_do_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_DO_B ~ 1, locations = winter_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = w_do_fvgm))
stopCluster(cl)

# merge all the predictions
w_do_train_merge = maptools::spRbind(w_do_train_par[[1]], w_do_train_par[[2]])
for (j in 3:length(w_do_train_par)) {
  w_do_train_merge = maptools::spRbind(w_do_train_merge, w_do_train_par[[j]])
}
w_do_train_merge = SpatialPixelsDataFrame(points = w_do_train_merge, data = w_do_train_merge@data)

# save data
summary(w_do_train_merge)
writeGDAL(w_do_train_merge["var1.pred"], fname = paste0(temp_wd, "mean_win_do_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_do_train = raster(paste0(temp_wd, "mean_win_do_train.tif"))
mean_win_do_train = writeRaster(raster::mask(raster::crop(mean_win_do_train, habitat_train), habitat_train),
                                 file = file.path(train_wd, "Environmental/Mean_Win_DO.asc"), format = "ascii",
                                 overwrite = T)

rm(list = c("w_do_train_merge", "w_do_train_par", "mean_win_do_train", "cl"))
showConnections()

#### summer temperature ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "train_grid", "train_parts", "s_temp_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_temp_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_TEMP_B ~ 1, locations = summer_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = s_temp_fvgm))
stopCluster(cl)

# merge all the predictions
s_temp_train_merge = maptools::spRbind(s_temp_train_par[[1]], s_temp_train_par[[2]])
for (j in 3:length(s_temp_train_par)) {
  s_temp_train_merge = maptools::spRbind(s_temp_train_merge, s_temp_train_par[[j]])
}
s_temp_train_merge = SpatialPixelsDataFrame(points = s_temp_train_merge, data = s_temp_train_merge@data)

# save data
summary(s_temp_train_merge)
writeGDAL(s_temp_train_merge["var1.pred"], fname = paste0(temp_wd, "mean_sum_temp_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_temp_train = raster(paste0(temp_wd, "mean_sum_temp_train.tif"))
mean_sum_temp_train = writeRaster(raster::mask(raster::crop(mean_sum_temp_train, habitat_train), habitat_train),
                                 file = file.path(train_wd, "Environmental/Mean_Sum_Temp.asc"), format = "ascii",
                                 overwrite = T)

rm(list = c("s_temp_train_merge", "s_temp_train_par", "mean_sum_temp_train", "cl"))
showConnections()

#### summer salinity ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "train_grid", "train_parts", "s_sal_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_sal_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_SAL_B ~ 1, locations = summer_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = s_sal_fvgm))
stopCluster(cl)

# merge all the predictions
s_sal_train_merge = maptools::spRbind(s_sal_train_par[[1]], s_sal_train_par[[2]])
for (j in 3:length(s_sal_train_par)) {
  s_sal_train_merge = maptools::spRbind(s_sal_train_merge, s_sal_train_par[[j]])
}
s_sal_train_merge = SpatialPixelsDataFrame(points = s_sal_train_merge, data = s_sal_train_merge@data)

# save data
summary(s_sal_train_merge)
writeGDAL(s_sal_train_merge["var1.pred"], fname = paste0(temp_wd, "mean_sum_sal_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_sal_train = raster(paste0(temp_wd, "mean_sum_sal_train.tif"))
mean_sum_sal_train = writeRaster(raster::mask(raster::crop(mean_sum_sal_train, habitat_train), habitat_train),
                                 file = file.path(train_wd, "Environmental/Mean_Sum_Sal.asc"), format = "ascii",
                                 overwrite = T)

rm(list = c("s_sal_train_merge", "s_sal_train_par", "mean_sum_sal_train", "cl"))
showConnections()

#### summer DO ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "train_grid", "train_parts", "s_do_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_do_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_DO_B ~ 1, locations = summer_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = s_do_fvgm))
stopCluster(cl)

# merge all the predictions
s_do_train_merge = maptools::spRbind(s_do_train_par[[1]], s_do_train_par[[2]])
for (j in 3:length(s_do_train_par)) {
  s_do_train_merge = maptools::spRbind(s_do_train_merge, s_do_train_par[[j]])
}
s_do_train_merge = SpatialPixelsDataFrame(points = s_do_train_merge, data = s_do_train_merge@data)

# save data
summary(s_do_train_merge)
writeGDAL(s_do_train_merge["var1.pred"], fname = paste0(temp_wd, "mean_sum_do_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_do_train = raster(paste0(temp_wd, "mean_sum_do_train.tif"))
mean_sum_do_train = writeRaster(raster::mask(raster::crop(mean_sum_do_train, habitat_train), habitat_train),
                                file = file.path(train_wd, "Environmental/Mean_Sum_DO.asc"), format = "ascii",
                                overwrite = T)

rm(list = c("s_do_train_merge", "s_do_train_par", "mean_sum_do_train", "cl"))
showConnections()

# remove the large training area grids before moving onto the testing area
rm(list = c("habitat_train", "train_grid", "train_parts"))

#### KRIGING: TESTING AREA ####
# habitat raster as a guide for kriging
habitat_test = raster(paste0(test_wd, "Environmental/Habitat.asc")) 
crs(habitat_test) = my_crs

# prediction grid
test_grid = raster(ncol = ncol(habitat_test), nrow = nrow(habitat_test), xmn = xmin(habitat_test), 
                    xmx = xmax(habitat_test), ymn = ymin(habitat_test), ymx = ymax(habitat_test))
train_grid = as(train_grid, "SpatialPixels") # convert to spatial pixels object
proj4string(train_grid) = proj4string(my_crs) # assign projection 

### winter temperature ####
# Initiate cluster 
cl = makeCluster(no_cores)

test_parts = split(x = 1:length(test_grid), f = 1:no_cores)

clusterExport(cl = cl, varlist = c("winter_wq_sp", "test_grid", "test_parts", "w_temp_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_temp_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_TEMP_B ~ 1, locations = winter_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = w_temp_fvgm))
stopCluster(cl)

# merge all the predictions
w_temp_test_merge = maptools::spRbind(w_temp_test_par[[1]], w_temp_test_par[[2]])
for (j in 3:length(w_temp_test_par)) {
  w_temp_test_merge = maptools::spRbind(w_temp_test_merge, w_temp_test_par[[j]])
}
w_temp_test_merge = SpatialPixelsDataFrame(points = w_temp_test_merge, data = w_temp_test_merge@data)

# save data
summary(w_temp_test_merge)
writeGDAL(w_temp_test_merge["var1.pred"], fname = paste0(temp_wd, "mean_win_temp_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_temp_test = raster(paste0(temp_wd, "mean_win_temp_test.tif"))
mean_win_temp_test = writeRaster(raster::mask(raster::crop(mean_win_temp_test, habitat_test), habitat_test),
                                  file = file.path(test_wd, "Environmental/Mean_Win_Temp.asc"), format = "ascii",
                                  overwrite = T)

rm(list = c("w_temp_test_merge", "w_temp_test_par", "mean_win_temp_test", "cl"))
showConnections()

#### winter salinity ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "test_grid", "test_parts", "w_sal_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_sal_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_SAL_B ~ 1, locations = winter_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = w_sal_fvgm))
stopCluster(cl)

# merge all the predictions
w_sal_test_merge = maptools::spRbind(w_sal_test_par[[1]], w_sal_test_par[[2]])
for (j in 3:length(w_sal_test_par)) {
  w_sal_test_merge = maptools::spRbind(w_sal_test_merge, w_sal_test_par[[j]])
}
w_sal_test_merge = SpatialPixelsDataFrame(points = w_sal_test_merge, data = w_sal_test_merge@data)

# save data
summary(w_sal_test_merge)
writeGDAL(w_sal_test_merge["var1.pred"], fname = paste0(temp_wd, "mean_win_sal_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_sal_test = raster(paste0(temp_wd, "mean_win_sal_test.tif"))
mean_win_sal_test = writeRaster(raster::mask(raster::crop(mean_win_sal_test, habitat_test), habitat_test),
                                 file = file.path(test_wd, "Environmental/Mean_Win_Sal.asc"), format = "ascii",
                                 overwrite = T)

rm(list = c("w_sal_test_merge", "w_sal_test_par", "mean_win_sal_test", "cl"))
showConnections()

#### winter DO ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "test_grid", "test_parts", "w_do_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_do_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_DO_B ~ 1, locations = winter_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = w_do_fvgm))
stopCluster(cl)

# merge all the predictions
w_do_test_merge = maptools::spRbind(w_do_test_par[[1]], w_do_test_par[[2]])
for (j in 3:length(w_do_test_par)) {
  w_do_test_merge = maptools::spRbind(w_do_test_merge, w_do_test_par[[j]])
}
w_do_test_merge = SpatialPixelsDataFrame(points = w_do_test_merge, data = w_do_test_merge@data)

# save data
summary(w_do_test_merge)
writeGDAL(w_do_test_merge["var1.pred"], fname = paste0(temp_wd, "mean_win_do_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_do_test = raster(paste0(temp_wd, "mean_win_do_test.tif"))
mean_win_do_test = writeRaster(raster::mask(raster::crop(mean_win_do_test, habitat_test), habitat_test),
                                file = file.path(test_wd, "Environmental/Mean_Win_DO.asc"), format = "ascii",
                                overwrite = T)

rm(list = c("w_do_test_merge", "w_do_test_par", "mean_win_do_test", "cl"))
showConnections()

#### summer temperature ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "test_grid", "test_parts", "s_temp_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_temp_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_TEMP_B ~ 1, locations = summer_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = s_temp_fvgm))
stopCluster(cl)

# merge all the predictions
s_temp_test_merge = maptools::spRbind(s_temp_test_par[[1]], s_temp_test_par[[2]])
for (j in 3:length(s_temp_test_par)) {
  s_temp_test_merge = maptools::spRbind(s_temp_test_merge, s_temp_test_par[[j]])
}
s_temp_test_merge = SpatialPixelsDataFrame(points = s_temp_test_merge, data = s_temp_test_merge@data)

# save data
summary(s_temp_test_merge)
writeGDAL(s_temp_test_merge["var1.pred"], fname = paste0(temp_wd, "mean_sum_temp_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_temp_test = raster(paste0(temp_wd, "mean_sum_temp_test.tif"))
mean_sum_temp_test = writeRaster(raster::mask(raster::crop(mean_sum_temp_test, habitat_test), habitat_test),
                                  file = file.path(test_wd, "Environmental/Mean_Sum_Temp.asc"), format = "ascii",
                                  overwrite = T)

rm(list = c("s_temp_test_merge", "s_temp_test_par", "mean_sum_temp_test", "cl"))
showConnections()

#### summer salinity ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "test_grid", "test_parts", "s_sal_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_sal_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_SAL_B ~ 1, locations = summer_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = s_sal_fvgm))
stopCluster(cl)

# merge all the predictions
s_sal_test_merge = maptools::spRbind(s_sal_test_par[[1]], s_sal_test_par[[2]])
for (j in 3:length(s_sal_test_par)) {
  s_sal_test_merge = maptools::spRbind(s_sal_test_merge, s_sal_test_par[[j]])
}
s_sal_test_merge = SpatialPixelsDataFrame(points = s_sal_test_merge, data = s_sal_test_merge@data)

# save data
summary(s_sal_test_merge)
writeGDAL(s_sal_test_merge["var1.pred"], fname = paste0(temp_wd, "mean_sum_sal_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_sal_test = raster(paste0(temp_wd, "mean_sum_sal_test.tif"))
mean_sum_sal_test = writeRaster(raster::mask(raster::crop(mean_sum_sal_test, habitat_test), habitat_test),
                                 file = file.path(test_wd, "Environmental/Mean_Sum_Sal.asc"), format = "ascii",
                                 overwrite = T)

rm(list = c("s_sal_test_merge", "s_sal_test_par", "mean_sum_sal_test", "cl"))
showConnections()

#### summer DO ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "test_grid", "test_parts", "s_do_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_do_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_DO_B ~ 1, locations = summer_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = s_do_fvgm))
stopCluster(cl)

# merge all the predictions
s_do_test_merge = maptools::spRbind(s_do_test_par[[1]], s_do_test_par[[2]])
for (j in 3:length(s_do_test_par)) {
  s_do_test_merge = maptools::spRbind(s_do_test_merge, s_do_test_par[[j]])
}
s_do_test_merge = SpatialPixelsDataFrame(points = s_do_test_merge, data = s_do_test_merge@data)

# save data
summary(s_do_test_merge)
writeGDAL(s_do_test_merge["var1.pred"], fname = paste0(temp_wd, "mean_sum_do_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_do_test = raster(paste0(temp_wd, "mean_sum_do_test.tif"))
mean_sum_do_test = writeRaster(raster::mask(raster::crop(mean_sum_do_test, habitat_test), habitat_test),
                                file = file.path(test_wd, "Environmental/Mean_Sum_DO.asc"), format = "ascii",
                                overwrite = T)

rm(list = c("s_do_test_merge", "s_do_test_par", "mean_sum_do_test", "cl"))
showConnections()

