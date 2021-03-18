# working directories
setwd("Z:/Courtney/Stuart_MSc_Ch1/") # main project folder
temp_wd = "Z:/Courtney/Stuart_MSc_Ch1/Temp/" # temporary files
source_wd = "Z:/Courtney/Stuart_MSc_Ch1/Source_Data/" # source data
dem_wd = "Z:/Courtney/Stuart_MSc_Ch1/Source_Data/Job606638_ncei_nintharcsec_dem/" # DEMs
csv_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Data/" # for writing tabular data
train_wd = "Z:/Courtney/Stuart_MSc_Ch1/Modeling_Data/Training/" # for training area data
test_wd = "Z:/Courtney/Stuart_MSc_Ch1/Modeling_Data/Testing/" # for testing area data

# libraries
library(easypackages)
libraries("tidyr", "rgdal", "gdalUtils", "raster", "sp", "sf", "tmap", "lwgeom", "rgeos",
          "cleangeo", "tidyverse", "stars", "fasterize", "PNWColors", "spex", "igraph", 
          "spatialEco", "tibble", "ncf", "spdep", "gstat", "geoR", "readxl", "dplyr",
          "parallel", "doParallel")

# change where large temp rasters are saved
rasterOptions(tmpdir = "Z:/Courtney/Stuart_MSc_Ch1/Temp/")

# save PROJ.4 string for standard projection (ESPG:26958 NAD 83/Florida East) and source LAT/LON data (EPSG:4326 WGS 84/World Geodetic System 1984)
my_crs = CRS("+init=epsg:26958")
gcs = CRS("+init=epsg:4326")
us_ft = CRS("+init=epsg:2236")

# full study domain for clipping water quality data
train = st_read(paste0(train_wd, "GIS/Training_Area.shp"))
test = st_read(paste0(test_wd, "GIS/Testing_Area.shp"))
domain = rbind(train, test)

### PREPARE WATER QUALITY DATA ####
# from The SERC Water Quality Monitoring Network (WQMN) (http://serc.fiu.edu/wqmnetwork/)
wqmn = read_excel(paste0(source_wd, "Water_Conditions/WQFloridaKeys&Shelf (ppm) UPDATED 6-6-2020.xlsx"), 
                sheet = "Data in ppm") 
head(wqmn, 5)

# pre-processing and cleaning 
wqmn = wqmn %>%
  mutate(DATE = as.Date(wqmn$DATE, origin = "1899-12-30")) %>% # this is the origin for excel sheets according to http://support.microsoft.com/kb/214330
  separate(DATE, into = c("YEAR", "MONTH", "DAY"), sep = "-", remove = F) %>% # separate out date data
  filter(YEAR %in% c("2014", "2015", "2016", "2017", "2018")) %>% # keeping only years of interest
  rename(TEMP = `TEMP-B`, SAL = `SAL-B`, DO = `DO-B`) %>% # rename variables because "-" causes problems
  st_as_sf(., coords = c("LONDEC", "LATDEC"), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) %>% # re-project to standard CRS
  mutate(LON_M = sf::st_coordinates(.)[,1], # save LON_M and LAT_M columns
         LAT_M = sf::st_coordinates(.)[,2]) %>%
  st_intersection(., st_make_valid(domain)) %>% # clip data to study domain 
  add_column(SEASON = NA) %>% # four seasons based on quarterly sampling 
  mutate(SEASON = cut(as.numeric(.$MONTH), breaks = c(0,3,6,9,Inf), labels = c("winter", "spring", "summer", "fall"))) %>%
  select(BASIN, SITE, YEAR, MONTH, DAY, SEASON, LON_M, LAT_M, TEMP, SAL, DO) 

# Biscayne Bay Water Quality (BBWQ) from Miami-Dade County Surface and Groundwater Quality Viewer (https://mdc.maps.arcgis.com/apps/webappviewer/index.html?id=3fd24515ee614f5db63924d7323a4ea7)
# BBWQ temperature data
bb_temp = read.csv(paste0(source_wd, "Water_Conditions/BBWQ_Temperature.csv"))
head(bb_temp, 5)

bb_temp = bb_temp %>%
  mutate(DATE = as.Date(bb_temp$DateCollected, format = '%d-%B-%y')) %>% 
  separate(DATE, into = c("YEAR", "MONTH", "DAY"), sep = "-", remove = F) %>% # separate out date data
  filter(YEAR %in% c("2014", "2015", "2016", "2017", "2018"), # keeping only years of interest
         DepthCode == "B") %>% # select bottom conditions
  rename(TEMP = Value, # rename columns to match WQMN data 
         SITE = StationUniqueID,
         BASIN = Water_Basin) %>%
  add_column(SEASON = NA) %>% # four seasons based on wqmn quarterly sampling 
  mutate(SEASON = cut(as.numeric(.$MONTH), breaks = c(0,3,6,9,Inf), labels = c("winter", "spring", "summer", "fall"))) %>%
  select(BASIN, SITE, YEAR, MONTH, DAY, SEASON, X_COORD, Y_COORD, TEMP)

# BBWQ saliniy data
bb_sal = read.csv(paste0(source_wd, "Water_Conditions/BBWQ_Salinity.csv"))
head(bb_sal, 5)

bb_sal = bb_sal %>%
  mutate(DATE = as.Date(bb_sal$DateCollected, format ='%B %d, %Y')) %>% 
  separate(DATE, into = c("YEAR", "MONTH", "DAY"), sep = "-", remove = F) %>% # separate out date data
  filter(YEAR %in% c("2014", "2015", "2016", "2017", "2018"), # keeping only years of interest
         DepthCode == "B") %>% # select bottom conditions
  rename(SAL = Value, # rename columns to match WQMN data 
         SITE = StationUniqueID,
         BASIN = Water_Basin) %>%
  add_column(SEASON = NA) %>% # four seasons based on WQMN quarterly sampling 
  mutate(SEASON = cut(as.numeric(.$MONTH), breaks = c(0,3,6,9,Inf), labels = c("winter", "spring", "summer", "fall"))) %>%
  select(BASIN, SITE, YEAR, MONTH, DAY, SEASON, X_COORD, Y_COORD, SAL)

# BBWQ dissolved oxygen data
bb_do = read.csv(paste0(source_wd, "Water_Conditions/BBWQ_Dissolved_Oxygen.csv"))
head(bb_do, 5)

bb_do = bb_do %>%
  mutate(DATE = as.Date(bb_do$DateCollected, format ='%B %d, %Y')) %>% 
  separate(DATE, into = c("YEAR", "MONTH", "DAY"), sep = "-", remove = F) %>% # separate out date data
  filter(YEAR %in% c("2014", "2015", "2016", "2017", "2018"), # keeping only years of interest
         DepthCode == "B") %>% # select bottom conditions
  rename(DO = Value, # rename columns to match WQMN data 
         SITE = StationUniqueID,
         BASIN = Water_Basin) %>%
  add_column(SEASON = NA) %>% # four seasons based on WQMN quarterly sampling 
  mutate(SEASON = cut(as.numeric(.$MONTH), breaks = c(0,3,6,9,Inf), labels = c("winter", "spring", "summer", "fall"))) %>%
  select(BASIN, SITE, YEAR, MONTH, DAY, SEASON, X_COORD, Y_COORD, DO)

#### WINTER DATA ####
# here are the wqmn stations that were sampled every winter for the full five year period
win_wqmn_sites = wqmn %>%
  select(SITE, YEAR, MONTH, DAY, SEASON, LON_M, LAT_M, TEMP, SAL, DO) %>%
  filter(SEASON == "winter", !is.na(TEMP), !is.na(SAL), !is.na(DO)) %>%
  group_by(SITE, LON_M, LAT_M) %>%
  count() %>%
  filter(n == 5) %>% # five years (2014, 2015, 2016, 2017, 2018)
  ungroup()

# restricting the wqmn data to only those sites in the win_wqmn_sites data frame and calculating summary stats
win_wqmn = wqmn %>%
  filter(SITE %in% win_wqmn_sites$SITE, SEASON == "winter") %>%
  group_by(SITE, LON_M, LAT_M) %>%
  mutate(MEAN_WIN_TEMP = mean(TEMP), SD_WIN_TEMP = sd(TEMP), # mean temp over the five years
         MEAN_WIN_SAL = mean(SAL), SD_WIN_SAL = sd(SAL), # mean sal over the five years
         MEAN_WIN_DO = mean(DO), SD_WIN_DO = sd(DO)) %>% # mean DO over the five years
  select(everything()) %>%
  distinct(BASIN, SITE, LON_M, LAT_M, MEAN_WIN_TEMP, SD_WIN_TEMP,
           MEAN_WIN_SAL, SD_WIN_SAL, MEAN_WIN_DO, SD_WIN_DO)


# annual mean temps, sal, do for biscayne bay
win_bb_temp = bb_temp %>%
  filter(SEASON == "winter") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_WIN_TEMP = mean(TEMP), SD_ANN_WIN_TEMP = sd(TEMP)) %>% # mean DO for each of the five years
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, ANN_MEAN_WIN_TEMP, SD_ANN_WIN_TEMP)

win_bb_sal = bb_sal %>%
  filter(SEASON == "winter") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_WIN_SAL = mean(SAL), SD_ANN_WIN_SAL = sd(SAL)) %>% # mean DO for each of the five years
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, ANN_MEAN_WIN_SAL, SD_ANN_WIN_SAL)

win_bb_do = bb_do %>%
  filter(SEASON == "winter") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_WIN_DO = mean(DO), SD_ANN_WIN_DO = sd(DO)) %>% # mean DO for each of the five years
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, ANN_MEAN_WIN_DO, SD_ANN_WIN_DO)

win_bb_temp_sites = win_bb_temp %>%
  group_by(SITE, X_COORD, Y_COORD) %>%
  count() %>%
  filter(n == 5) %>% # five years (2014, 2015, 2016, 2017, 2018)
  ungroup()

win_bb_sal_sites = win_bb_sal %>%
  group_by(SITE, X_COORD, Y_COORD) %>%
  count() %>%
  filter(n == 5) %>% # five years (2014, 2015, 2016, 2017, 2018)
  ungroup()

win_bb_do_sites = win_bb_do %>%
  group_by(SITE, X_COORD, Y_COORD) %>%
  count() %>%
  filter(n == 5) %>% # five years (2014, 2015, 2016, 2017, 2018)
  ungroup()

# restricting the bbwq data to only those sites that were sampled in each of the five years
# calculating summary stats from the annual means 
win_bb_temp = win_bb_temp %>%
  filter(SITE %in% win_bb_temp_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_WIN_TEMP = mean(ANN_MEAN_WIN_TEMP), SD_WIN_TEMP = sd(ANN_MEAN_WIN_TEMP)) %>% # overall mean temp
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_WIN_TEMP, SD_WIN_TEMP)

win_bb_sal = win_bb_sal %>%
  filter(SITE %in% win_bb_sal_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_WIN_SAL = mean(ANN_MEAN_WIN_SAL), SD_WIN_SAL = sd(ANN_MEAN_WIN_SAL)) %>% # overall mean sal
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_WIN_SAL, SD_WIN_SAL)

win_bb_do = win_bb_do %>%
  filter(SITE %in% win_bb_do_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_WIN_DO = mean(ANN_MEAN_WIN_DO), SD_WIN_DO = sd(ANN_MEAN_WIN_DO)) %>% # overall mean DO
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_WIN_DO, SD_WIN_DO)

# combine the temperature, salinity, and DO data from Biscayne Bay and clip to study domain
win_bbwq = left_join(win_bb_temp, win_bb_sal)
win_bbwq = left_join(win_bbwq, win_bb_do)

win_bbwq = win_bbwq %>%
  st_as_sf(., coords = c("X_COORD", "Y_COORD"), crs = us_ft) %>% # convert to spatial
  st_transform(., my_crs) %>% # re-project to standard CRS
  add_column(LON_M = sf::st_coordinates(.)[,1], # save LON_M and LAT_M columns
         LAT_M = sf::st_coordinates(.)[,2]) %>%
  st_intersection(., st_make_valid(domain)) %>% # clip data to study domain 
  distinct(BASIN, SITE, LON_M, LAT_M, MEAN_WIN_TEMP, SD_WIN_TEMP,
           MEAN_WIN_SAL, SD_WIN_SAL, MEAN_WIN_DO, SD_WIN_DO)

# remove intermediate biscayne bay data tables
rm(list = c("win_bb_temp", "win_bb_temp_sites" ,"win_bb_sal", "win_bb_sal_sites",
            "win_bb_do", "win_bb_do_sites"))

# combining the winter water quality data from biscayne bay and the fknms
winter_wq = rbind(win_wqmn, win_bbwq)

# for variogram modeling and spatial estimation we need SPDF, so convert sf object to SPDF
winter_wq_sp = winter_wq %>% st_drop_geometry()
coordinates(winter_wq_sp) = ~ LON_M + LAT_M
proj4string(winter_wq_sp) = proj4string(my_crs) # use my_crs CRS object to define proj4 string of the SPDF
summary(winter_wq_sp)
write.csv(winter_wq_sp, paste0(csv_wd, "Winter_Water_Conditions.csv")) 

# first taking a glimpse at correlograms and the Moran's I values for the data to ensure that there is spatial dependence
winter_coords = cbind(winter_wq_sp$LON_M, winter_wq_sp$LAT_M) # calculate a distance matrix
colnames(winter_coords) = c("LON_M", "LAT_M")
winter_distmat = as.matrix(dist(winter_coords))
winter_maxdist = 2/3 * max(winter_distmat) # maximum distance to consider in correlogram/variogram

# spline correlograms with 95% pointwise bootstrap CIs
w_temp_corr = spline.correlog(x = winter_wq_sp$LON_M, y = winter_wq_sp$LAT_M, z = winter_wq_sp$MEAN_WIN_TEMP, xmax = winter_maxdist, resamp = 100, type = "boot")
plot(w_temp_corr) # with 95% CIs from bootstrapping
w_sal_corr = spline.correlog(x = winter_wq_sp$LON_M, y = winter_wq_sp$LAT_M, z = winter_wq_sp$MEAN_WIN_SAL, xmax = winter_maxdist, resamp = 100, type = "boot")
plot(w_sal_corr) 
w_do_corr = spline.correlog(x = winter_wq_sp$LON_M, y = winter_wq_sp$LAT_M, z = winter_wq_sp$MEAN_WIN_DO, xmax = winter_maxdist, resamp = 100, type = "boot")
plot(w_do_corr)

# neighborhood list (neighbors within 20 km distance so every site has at least one neighbor)
winter_neigh = dnearneigh(x = winter_coords, d1 = 0, d2 = 20000, longlat = F)
plot(winter_neigh, coordinates(winter_coords))
winter_wts = nb2listw(neighbours = winter_neigh, style = "W", zero.policy = T) # weights matrix for calculating Moran's I

# Moran's I with normal approximations
moran.test(winter_wq_sp$MEAN_WIN_TEMP, listw = winter_wts, randomisation = F, zero.policy = T) # est. Moran's I stat = 0.585978759, p < 0.01 sig. spatial dependence
moran.test(winter_wq_sp$MEAN_WIN_SAL, listw = winter_wts, randomisation = F, zero.policy = T) # est. Moran's I stat = 0.288572497, p < 0.01 sig. spatial dependence
moran.test(winter_wq_sp$MEAN_WIN_DO, listw = winter_wts, randomisation = F, zero.policy = T) # est. Moran's I stat = 0.244604117, p < 0.01 sig. spatial dependence

# Moran's I with Monte Carlo permutations, does everything match with the normal approximations?
moran.mc(winter_wq_sp$MEAN_WIN_TEMP, listw = winter_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.58598, p = 0.01 sig. spatial dependence
moran.mc(winter_wq_sp$MEAN_WIN_SAL, listw = winter_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.28857, p = 0.01 sig. spatial dependence
moran.mc(winter_wq_sp$MEAN_WIN_DO, listw = winter_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.2446, p = 0.01 sig. spatial dependence

# according to Moran tests with both normal and Monte Carlo approximations, there is significant spatial dependence 
# in the winter temp, sal, and DO data.
# now  move on to variogram modeling to capture the spatial structure...
require(gstat)

w_temp_evgm = variogram(MEAN_WIN_TEMP ~ 1, winter_wq_sp, cutoff = winter_maxdist) # empirical variogram
plot(w_temp_evgm, xlab = "distance (m)", pch = 19)
w_temp_fvgm = fit.variogram(w_temp_evgm, vgm(psill = 0.8, model = "Sph", range = 70000, nugget = 0)) # fit variogram
w_temp_svgm_plot = plot(w_temp_evgm, model = w_temp_fvgm, xlab = "distance (m)", pch = 19)
w_temp_svgm_plot
print(w_temp_fvgm)

w_sal_evgm = variogram(MEAN_WIN_SAL ~ 1, winter_wq_sp, cutoff = winter_maxdist)
plot(w_sal_evgm, xlab = "distance (m)", pch = 19)
w_sal_fvgm = fit.variogram(w_sal_evgm, vgm(psill = 5, model = "Sph", range = 20000, nugget = 3))
w_sal_svgm_plot = plot(w_sal_evgm, model = w_sal_fvgm, xlab = "distance (m)", pch = 19)
w_sal_svgm_plot
print(w_sal_fvgm)

w_do_evgm = variogram(MEAN_WIN_DO ~ 1, winter_wq_sp, cutoff = winter_maxdist)
plot(w_do_evgm, xlab = "distance (m)", pch = 19)
w_do_fvgm = fit.variogram(w_do_evgm, vgm(model = "Sph"), fit.sills = TRUE, fit.ranges = TRUE)
w_do_fvgm_plot = plot(w_do_evgm, model = w_do_fvgm, xlab = "distance (m)", pch = 19)
w_do_fvgm_plot
print(w_do_fvgm)


#### SUMMER DATA ####
# here are the stations that were sampled every summer for the full five year period
sum_wqmn_sites = wqmn %>%
  select(SITE, YEAR, MONTH, SEASON, LON_M, LAT_M, TEMP, SAL, DO) %>%
  filter(SEASON == "summer", !is.na(TEMP), !is.na(SAL), !is.na(DO)) %>%
  group_by(SITE, LON_M, LAT_M) %>%
  count() %>%
  filter(n == 5) %>%
  ungroup()

# restricting the wqmn data frame to only those sites in the sum_wqmn_sites data frame and calculating summary stats
sum_wqmn = wqmn %>%
  filter(SITE %in% sum_wqmn_sites$SITE, SEASON == "summer") %>%
  group_by(SITE, LON_M, LAT_M) %>%
  mutate(MEAN_SUM_TEMP = mean(TEMP), SD_SUM_TEMP = sd(TEMP),
         MEAN_SUM_SAL = mean(SAL), SD_SUM_SAL = sd(SAL),
         MEAN_SUM_DO = mean(DO), SD_SUM_DO = sd(DO)) %>%
  select(everything()) %>%
  distinct(BASIN, SITE, LON_M, LAT_M, MEAN_SUM_TEMP, SD_SUM_TEMP,
           MEAN_SUM_SAL, SD_SUM_SAL, MEAN_SUM_DO, SD_SUM_DO)

# annual mean temps, sal, do for biscayne bay
sum_bb_temp = bb_temp %>%
  filter(SEASON == "summer") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_SUM_TEMP = mean(TEMP), SD_ANN_SUM_TEMP = sd(TEMP)) %>% # mean temp for each of the five years
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, ANN_MEAN_SUM_TEMP, SD_ANN_SUM_TEMP)

sum_bb_sal = bb_sal %>%
  filter(SEASON == "summer") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_SUM_SAL = mean(SAL), SD_ANN_SUM_SAL = sd(SAL)) %>% # mean sal for each of the five years
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, ANN_MEAN_SUM_SAL, SD_ANN_SUM_SAL)

sum_bb_do = bb_do %>%
  filter(SEASON == "summer") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_SUM_DO = mean(DO), SD_ANN_SUM_DO = sd(DO)) %>% # mean DO for each of the five years
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, ANN_MEAN_SUM_DO, SD_ANN_SUM_DO)

sum_bb_temp_sites = sum_bb_temp %>%
  group_by(SITE, X_COORD, Y_COORD) %>%
  count() %>%
  filter(n == 5) %>% # five years (2014, 2015, 2016, 2017, 2018)
  ungroup()

sum_bb_sal_sites = sum_bb_sal %>%
  group_by(SITE, X_COORD, Y_COORD) %>%
  count() %>%
  filter(n == 5) %>% # five years (2014, 2015, 2016, 2017, 2018)
  ungroup()

sum_bb_do_sites = sum_bb_do %>%
  group_by(SITE, X_COORD, Y_COORD) %>%
  count() %>%
  filter(n == 5) %>% # five years (2014, 2015, 2016, 2017, 2018)
  ungroup()

# restricting the bbwq data to only those sites that were sampled in each of the five years
# calculating summary stats from the annual means 
sum_bb_temp = sum_bb_temp %>%
  filter(SITE %in% sum_bb_temp_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_SUM_TEMP = mean(ANN_MEAN_SUM_TEMP), SD_SUM_TEMP = sd(ANN_MEAN_SUM_TEMP)) %>% # overall mean temp
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_SUM_TEMP, SD_SUM_TEMP)

sum_bb_sal = sum_bb_sal %>%
  filter(SITE %in% sum_bb_sal_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_SUM_SAL = mean(ANN_MEAN_SUM_SAL), SD_SUM_SAL = sd(ANN_MEAN_SUM_SAL)) %>% # overall mean sal
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_SUM_SAL, SD_SUM_SAL)

sum_bb_do = sum_bb_do %>%
  filter(SITE %in% sum_bb_do_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_SUM_DO = mean(ANN_MEAN_SUM_DO), SD_SUM_DO = sd(ANN_MEAN_SUM_DO)) %>% # overall mean DO
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_SUM_DO, SD_SUM_DO)

# combine the temperature, salinity, and DO data from Biscayne Bay and clip to study domain
sum_bbwq = left_join(sum_bb_temp, sum_bb_sal)
sum_bbwq = left_join(sum_bbwq, sum_bb_do)

sum_bbwq = sum_bbwq %>%
  st_as_sf(., coords = c("X_COORD", "Y_COORD"), crs = us_ft) %>% # convert to spatial
  st_transform(., my_crs) %>% # re-project to standard CRS
  add_column(LON_M = sf::st_coordinates(.)[,1], # save LON_M and LAT_M columns
             LAT_M = sf::st_coordinates(.)[,2]) %>%
  st_intersection(., st_make_valid(domain)) %>% # clip data to study domain 
  distinct(BASIN, SITE, LON_M, LAT_M, MEAN_SUM_TEMP, SD_SUM_TEMP,
           MEAN_SUM_SAL, SD_SUM_SAL, MEAN_SUM_DO, SD_SUM_DO)

# remove intermediate biscayne bay data tables
rm(list = c("sum_bb_temp", "sum_bb_temp_sites", "sum_bb_sal", "sum_bb_sal_sites",
            "sum_bb_do", "sum_bb_do_sites"))

# combining the summer water quality data from biscayne bay and the fknms
summer_wq = rbind(sum_wqmn, sum_bbwq)

# for variogram modeling and spatial estimation we need SPDF, so convert sf object to SPDF
summer_wq_sp = summer_wq %>% st_drop_geometry()
coordinates(summer_wq_sp) = ~ LON_M + LAT_M
proj4string(summer_wq_sp) = proj4string(my_crs) # use my_crs CRS object to define proj4 string of the SPDF
summary(summer_wq_sp)
write.csv(summer_wq_sp, paste0(csv_wd, "Summer_Water_Conditions.csv"))


# first taking a glimpse at correlograms and the Moran's I values for the data to ensure that there is spatial dependence
summer_coords = cbind(summer_wq_sp$LON_M, summer_wq_sp$LAT_M) #calculate a distance matrix
colnames(summer_coords) = c("LON_M", "LAT_M")
summer_distmat = as.matrix(dist(summer_coords))
summer_maxdist = 2/3 * max(summer_distmat) #maximum distance to consider in correlogram/variogram

# spline correlograms with 95% pointwise bootstrap CIs
s_temp_corr = spline.correlog(x = summer_wq_sp$LON_M, y = summer_wq_sp$LAT_M, z = summer_wq_sp$MEAN_SUM_TEMP, xmax = summer_maxdist, resamp = 100, type = "boot")
plot(s_temp_corr) # with 95% CIs from bootstrapping
s_sal_corr = spline.correlog(x = summer_wq_sp$LON_M, y = summer_wq_sp$LAT_M, z = summer_wq_sp$MEAN_SUM_SAL, xmax = summer_maxdist, resamp = 100, type = "boot")
plot(s_sal_corr) 
s_do_corr = spline.correlog(x = summer_wq_sp$LON_M, y = summer_wq_sp$LAT_M, z = summer_wq_sp$MEAN_SUM_DO, xmax = summer_maxdist, resamp = 100, type = "boot")
plot(s_do_corr)

# neighborhood list (neighbors within 20 km distance so each site has at least one neighbor)
summer_neigh = dnearneigh(x = summer_coords, d1 = 0, d2 = 20000, longlat = F)
plot(summer_neigh, coordinates(summer_coords))
summer_wts = nb2listw(neighbours = summer_neigh, style = "W", zero.policy = T) # weights matrix for calculating Moran's I


# Moran's I with normal approximations
moran.test(summer_wq_sp$MEAN_SUM_TEMP, listw = summer_wts, randomisation = F, zero.policy = T)  # est. Moran's I stat = 0.155536517, p < 0.01 sig. spatial dependence
moran.test(summer_wq_sp$MEAN_SUM_SAL, listw = summer_wts, randomisation = F, zero.policy = T) # est. Moran's I stat =  0.186665202, p < 0.01 sig. spatial dependence
moran.test(summer_wq_sp$MEAN_SUM_DO, listw = summer_wts, randomisation = F, zero.policy = T) # est. Moran's I stat = 0.254710245, p < 0.01 sig. spatial dependence


# Moran's I with Monte Carlo permutations, does everything match with the normal approximations?
moran.mc(summer_wq_sp$MEAN_SUM_TEMP, listw = summer_wts, nsim = 99, zero.policy = T)  # est. Moran's I stat = 0.15554, p = 0.01 sig. spatial dependence
moran.mc(summer_wq_sp$MEAN_SUM_SAL, listw = summer_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.18667, p = 0.01 sig. spatial dependence
moran.mc(summer_wq_sp$MEAN_SUM_DO, listw = summer_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.25471, p = 0.01 sig. spatial dependence


# according to Moran tests with both normal and Monte Carlo estimations, 
# there is significant spatial dependence in the summer temp, sal, and DO data.
# now move on to variogram modeling to capture the spatial structure...
require(gstat)

s_temp_evgm = variogram(MEAN_SUM_TEMP ~ 1, summer_wq_sp, cutoff = summer_maxdist) # empirical variogram
plot(s_temp_evgm, xlab = "distance (m)", pch = 19)
s_temp_fvgm = fit.variogram(s_temp_evgm, vgm(psill = 0.1, model = "Sph", range = 30000, nugget = 0.2))
s_temp_svgm_plot = plot(s_temp_evgm, model = s_temp_fvgm, xlab = "distance (m)", pch = 19)
s_temp_svgm_plot
print(s_temp_fvgm)

s_sal_evgm = variogram(MEAN_SUM_SAL ~ 1, summer_wq_sp, cutoff = summer_maxdist)
plot(s_sal_evgm, xlab = "distance (m)", pch = 19)
s_sal_fvgm = fit.variogram(s_sal_evgm, vgm(model = "Sph"), fit.sills = TRUE, fit.ranges = TRUE)
s_sal_fvgm_plot = plot(s_sal_evgm, model = s_sal_fvgm, xlab = "distance (m)", pch = 19)
s_sal_fvgm_plot
print(s_sal_fvgm)

s_do_evgm = variogram(MEAN_SUM_DO ~ 1, summer_wq_sp, cutoff = summer_maxdist)
plot(s_do_evgm, xlab = "distance (m)", pch = 19)
s_do_fvgm = fit.variogram(s_do_evgm, vgm(psill = 0.05, model = "Sph", range = 20000, nugget = 0.01))
s_do_svgm_plot = plot(s_do_evgm, model = s_do_fvgm, xlab = "distance (m)", pch = 19)
s_do_svgm_plot
print(s_do_fvgm)

#### KRIGING: TRAINING AREA ####
# prediction grid
habitat_train = raster(paste0(train_wd, "Environmental/Habitat.asc")) 
crs(habitat_train) = my_crs

train_grid = raster(ncol = ncol(habitat_train), 
                    nrow = nrow(habitat_train),
                    xmn = xmin(habitat_train),
                    xmx = xmax(habitat_train),
                    ymn = ymin(habitat_train), 
                    ymx = ymax(habitat_train)) %>%
  as(., "SpatialPixels")
proj4string(train_grid) = proj4string(my_crs) # assign projection 

# ordinary kriging (running in parallel) using empirical semivariograms calculated above
# starting with winter conditions in the training area
library(doParallel)
library(parallel)

#### winter temperature ####
# calculate the number of cores
no_cores = detectCores() - 2

# initiate cluster 
cl = makeCluster(no_cores)

# split training area into pieces for each core
train_parts = split(x = 1:length(train_grid), f = 1:no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "train_grid", "train_parts", "w_temp_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_temp_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_TEMP ~ 1, locations = winter_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = w_temp_fvgm))
stopCluster(cl)
showConnections()

# merge all the predictions
w_temp_train_merge = maptools::spRbind(w_temp_train_par[[1]], w_temp_train_par[[2]])
for (j in 3:length(w_temp_train_par)) {
  w_temp_train_merge = maptools::spRbind(w_temp_train_merge, w_temp_train_par[[j]])
}
w_temp_train_merge = SpatialPixelsDataFrame(points = w_temp_train_merge, data = w_temp_train_merge@data)

# save new surface
summary(w_temp_train_merge)
writeGDAL(w_temp_train_merge["var1.pred"], fname = paste0(temp_wd, "win_temp_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_temp_train = raster(paste0(temp_wd, "win_temp_train.tif"))
mean_win_temp_train = writeRaster(raster::mask(raster::crop(mean_win_temp_train, habitat_train), habitat_train),
                                  file = file.path(train_wd, "Environmental/Win_Temp.asc"), format = "ascii",
                                  overwrite = T)

rm(list = c("w_temp_train_merge", "w_temp_train_par", "mean_win_temp_train", "cl"))
showConnections()

#### winter salinity ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "train_grid", "train_parts", "w_sal_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_sal_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_SAL ~ 1, locations = winter_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = w_sal_fvgm))
stopCluster(cl)

# merge all the predictions
w_sal_train_merge = maptools::spRbind(w_sal_train_par[[1]], w_sal_train_par[[2]])
for (j in 3:length(w_sal_train_par)) {
  w_sal_train_merge = maptools::spRbind(w_sal_train_merge, w_sal_train_par[[j]])
}
w_sal_train_merge = SpatialPixelsDataFrame(points = w_sal_train_merge, data = w_sal_train_merge@data)

# save new surface
summary(w_sal_train_merge)
writeGDAL(w_sal_train_merge["var1.pred"], fname = paste0(temp_wd, "win_sal_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_sal_train = raster(paste0(temp_wd, "win_sal_train.tif"))
mean_win_sal_train = writeRaster(raster::mask(raster::crop(mean_win_sal_train, habitat_train), habitat_train),
                                 file = file.path(train_wd, "Environmental/Win_Sal.asc"), format = "ascii",
                                 overwrite = T)

rm(list = c("w_sal_train_merge", "w_sal_train_par", "mean_win_sal_train", "cl"))
showConnections()

#### winter DO ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "train_grid", "train_parts", "w_do_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_do_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_DO ~ 1, locations = winter_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = w_do_fvgm))
stopCluster(cl)

# merge all the predictions
w_do_train_merge = maptools::spRbind(w_do_train_par[[1]], w_do_train_par[[2]])
for (j in 3:length(w_do_train_par)) {
  w_do_train_merge = maptools::spRbind(w_do_train_merge, w_do_train_par[[j]])
}
w_do_train_merge = SpatialPixelsDataFrame(points = w_do_train_merge, data = w_do_train_merge@data)

# save new surface
summary(w_do_train_merge)
writeGDAL(w_do_train_merge["var1.pred"], fname = paste0(temp_wd, "win_do_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_do_train = raster(paste0(temp_wd, "win_do_train.tif"))
mean_win_do_train = writeRaster(raster::mask(raster::crop(mean_win_do_train, habitat_train), habitat_train),
                                file = file.path(train_wd, "Environmental/Win_DO.asc"), format = "ascii",
                                overwrite = T)

rm(list = c("w_do_train_merge", "w_do_train_par", "mean_win_do_train", "cl"))
showConnections()

#### summer temperature ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "train_grid", "train_parts", "s_temp_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_temp_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_TEMP ~ 1, locations = summer_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = s_temp_fvgm))
stopCluster(cl)

# merge all the predictions
s_temp_train_merge = maptools::spRbind(s_temp_train_par[[1]], s_temp_train_par[[2]])
for (j in 3:length(s_temp_train_par)) {
  s_temp_train_merge = maptools::spRbind(s_temp_train_merge, s_temp_train_par[[j]])
}
s_temp_train_merge = SpatialPixelsDataFrame(points = s_temp_train_merge, data = s_temp_train_merge@data)

# save new surface
summary(s_temp_train_merge)
writeGDAL(s_temp_train_merge["var1.pred"], fname = paste0(temp_wd, "sum_temp_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_temp_train = raster(paste0(temp_wd, "sum_temp_train.tif"))
mean_sum_temp_train = writeRaster(raster::mask(raster::crop(mean_sum_temp_train, habitat_train), habitat_train),
                                  file = file.path(train_wd, "Environmental/Sum_Temp.asc"), format = "ascii",
                                  overwrite = T)

rm(list = c("s_temp_train_merge", "s_temp_train_par", "mean_sum_temp_train", "cl"))
showConnections()


#### summer salinity ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "train_grid", "train_parts", "s_sal_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_sal_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_SAL ~ 1, locations = summer_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = s_sal_fvgm))
stopCluster(cl)

# merge all the predictions
s_sal_train_merge = maptools::spRbind(s_sal_train_par[[1]], s_sal_train_par[[2]])
for (j in 3:length(s_sal_train_par)) {
  s_sal_train_merge = maptools::spRbind(s_sal_train_merge, s_sal_train_par[[j]])
}
s_sal_train_merge = SpatialPixelsDataFrame(points = s_sal_train_merge, data = s_sal_train_merge@data)

# save new surface
summary(s_sal_train_merge)
writeGDAL(s_sal_train_merge["var1.pred"], fname = paste0(temp_wd, "sum_sal_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_sal_train = raster(paste0(temp_wd, "sum_sal_train.tif"))
mean_sum_sal_train = writeRaster(raster::mask(raster::crop(mean_sum_sal_train, habitat_train), habitat_train),
                                 file = file.path(train_wd, "Environmental/Sum_Sal.asc"), format = "ascii",
                                 overwrite = T)

rm(list = c("s_sal_train_merge", "s_sal_train_par", "mean_sum_sal_train", "cl"))
showConnections()

#### summer DO ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "train_grid", "train_parts", "s_do_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_do_train_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_DO ~ 1, locations = summer_wq_sp,
        newdata = train_grid[train_parts[[x]],], model = s_do_fvgm))
stopCluster(cl)

# merge all the predictions
s_do_train_merge = maptools::spRbind(s_do_train_par[[1]], s_do_train_par[[2]])
for (j in 3:length(s_do_train_par)) {
  s_do_train_merge = maptools::spRbind(s_do_train_merge, s_do_train_par[[j]])
}
s_do_train_merge = SpatialPixelsDataFrame(points = s_do_train_merge, data = s_do_train_merge@data)

# save new surface
summary(s_do_train_merge)
writeGDAL(s_do_train_merge["var1.pred"], fname = paste0(temp_wd, "sum_do_train.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_do_train = raster(paste0(temp_wd, "sum_do_train.tif"))
mean_sum_do_train = writeRaster(raster::mask(raster::crop(mean_sum_do_train, habitat_train), habitat_train),
                                file = file.path(train_wd, "Environmental/Sum_DO.asc"), format = "ascii",
                                overwrite = T)

rm(list = c("s_do_train_merge", "s_do_train_par", "mean_sum_do_train", "cl"))
showConnections()

# remove the large training area grids before moving onto the testing area
rm(list = c("habitat_train", "train_grid", "train_parts"))


#### KRIGING: TESTING AREA ####
# habitat raster as a guide for kriging
# here it is if you need to reload it
habitat_test = raster(paste0(test_wd, "Environmental/Habitat.asc")) 
crs(habitat_test) = my_crs

# create prediction grid for testing area
test_grid = raster(ncol = ncol(habitat_test),
                   nrow = nrow(habitat_test),
                   xmn = xmin(habitat_test),
                   xmx = xmax(habitat_test),
                   ymn = ymin(habitat_test),
                   ymx = ymax(habitat_test)) %>%
  as(., "SpatialPixels")
proj4string(test_grid) = proj4string(my_crs) # assign projection 

### winter temperature ####
# initiate cluster 
cl = makeCluster(no_cores)

# split testing area into pieces for each of the cores
test_parts = split(x = 1:length(test_grid), f = 1:no_cores)

clusterExport(cl = cl, varlist = c("winter_wq_sp", "test_grid", "test_parts", "w_temp_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_temp_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_TEMP ~ 1, locations = winter_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = w_temp_fvgm))
stopCluster(cl)

# merge all the predictions
w_temp_test_merge = maptools::spRbind(w_temp_test_par[[1]], w_temp_test_par[[2]])
for (j in 3:length(w_temp_test_par)) {
  w_temp_test_merge = maptools::spRbind(w_temp_test_merge, w_temp_test_par[[j]])
}
w_temp_test_merge = SpatialPixelsDataFrame(points = w_temp_test_merge, data = w_temp_test_merge@data)

# save new surface
summary(w_temp_test_merge)
writeGDAL(w_temp_test_merge["var1.pred"], fname = paste0(temp_wd, "win_temp_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_temp_test = raster(paste0(temp_wd, "win_temp_test.tif"))
mean_win_temp_test = writeRaster(raster::mask(raster::crop(mean_win_temp_test, habitat_test), habitat_test),
                                 file = file.path(test_wd, "Environmental/Win_Temp.asc"), format = "ascii",
                                 overwrite = T)
compareRaster(mean_win_temp_test, habitat_test, res = T, extent = T, rowcol = T)
rm(list = c("w_temp_test_merge", "w_temp_test_par", "mean_win_temp_test", "cl"))
showConnections()

#### winter salinity ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "test_grid", "test_parts", "w_sal_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_sal_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_SAL ~ 1, locations = winter_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = w_sal_fvgm))
stopCluster(cl)

# merge all the predictions
w_sal_test_merge = maptools::spRbind(w_sal_test_par[[1]], w_sal_test_par[[2]])
for (j in 3:length(w_sal_test_par)) {
  w_sal_test_merge = maptools::spRbind(w_sal_test_merge, w_sal_test_par[[j]])
}
w_sal_test_merge = SpatialPixelsDataFrame(points = w_sal_test_merge, data = w_sal_test_merge@data)

# save new surface
summary(w_sal_test_merge)
writeGDAL(w_sal_test_merge["var1.pred"], fname = paste0(temp_wd, "win_sal_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_sal_test = raster(paste0(temp_wd, "win_sal_test.tif"))
mean_win_sal_test = writeRaster(raster::mask(raster::crop(mean_win_sal_test, habitat_test), habitat_test),
                                file = file.path(test_wd, "Environmental/Win_Sal.asc"), format = "ascii",
                                overwrite = T)

rm(list = c("w_sal_test_merge", "w_sal_test_par", "mean_win_sal_test", "cl"))
showConnections()

#### winter DO ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "test_grid", "test_parts", "w_do_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_do_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_DO ~ 1, locations = winter_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = w_do_fvgm))
stopCluster(cl)

# merge all the predictions
w_do_test_merge = maptools::spRbind(w_do_test_par[[1]], w_do_test_par[[2]])
for (j in 3:length(w_do_test_par)) {
  w_do_test_merge = maptools::spRbind(w_do_test_merge, w_do_test_par[[j]])
}
w_do_test_merge = SpatialPixelsDataFrame(points = w_do_test_merge, data = w_do_test_merge@data)

# save new surface
summary(w_do_test_merge)
writeGDAL(w_do_test_merge["var1.pred"], fname = paste0(temp_wd, "win_do_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_do_test = raster(paste0(temp_wd, "win_do_test.tif"))
mean_win_do_test = writeRaster(raster::mask(raster::crop(mean_win_do_test, habitat_test), habitat_test),
                               file = file.path(test_wd, "Environmental/Win_DO.asc"), format = "ascii",
                               overwrite = T)

rm(list = c("w_do_test_merge", "w_do_test_par", "mean_win_do_test", "cl"))
showConnections()


#### summer temperature ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "test_grid", "test_parts", "s_temp_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_temp_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_TEMP ~ 1, locations = summer_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = s_temp_fvgm))
stopCluster(cl)

# merge all the predictions
s_temp_test_merge = maptools::spRbind(s_temp_test_par[[1]], s_temp_test_par[[2]])
for (j in 3:length(s_temp_test_par)) {
  s_temp_test_merge = maptools::spRbind(s_temp_test_merge, s_temp_test_par[[j]])
}
s_temp_test_merge = SpatialPixelsDataFrame(points = s_temp_test_merge, data = s_temp_test_merge@data)

# save new surface
summary(s_temp_test_merge)
writeGDAL(s_temp_test_merge["var1.pred"], fname = paste0(temp_wd, "sum_temp_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_temp_test = raster(paste0(temp_wd, "sum_temp_test.tif"))
mean_sum_temp_test = writeRaster(raster::mask(raster::crop(mean_sum_temp_test, habitat_test), habitat_test),
                                 file = file.path(test_wd, "Environmental/Sum_Temp.asc"), format = "ascii",
                                 overwrite = T)

rm(list = c("s_temp_test_merge", "s_temp_test_par", "mean_sum_temp_test", "cl"))
showConnections()

#### summer salinity ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "test_grid", "test_parts", "s_sal_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_sal_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_SAL ~ 1, locations = summer_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = s_sal_fvgm))
stopCluster(cl)

# merge all the predictions
s_sal_test_merge = maptools::spRbind(s_sal_test_par[[1]], s_sal_test_par[[2]])
for (j in 3:length(s_sal_test_par)) {
  s_sal_test_merge = maptools::spRbind(s_sal_test_merge, s_sal_test_par[[j]])
}
s_sal_test_merge = SpatialPixelsDataFrame(points = s_sal_test_merge, data = s_sal_test_merge@data)

# save new surface
summary(s_sal_test_merge)
writeGDAL(s_sal_test_merge["var1.pred"], fname = paste0(temp_wd, "sum_sal_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_sal_test = raster(paste0(temp_wd, "sum_sal_test.tif"))
mean_sum_sal_test = writeRaster(raster::mask(raster::crop(mean_sum_sal_test, habitat_test), habitat_test),
                                file = file.path(test_wd, "Environmental/Sum_Sal.asc"), format = "ascii",
                                overwrite = T)

rm(list = c("s_sal_test_merge", "s_sal_test_par", "mean_sum_sal_test", "cl"))
showConnections()

#### summer DO ####
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "test_grid", "test_parts", "s_do_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_do_test_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_DO ~ 1, locations = summer_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = s_do_fvgm))
stopCluster(cl)

# merge all the predictions
s_do_test_merge = maptools::spRbind(s_do_test_par[[1]], s_do_test_par[[2]])
for (j in 3:length(s_do_test_par)) {
  s_do_test_merge = maptools::spRbind(s_do_test_merge, s_do_test_par[[j]])
}
s_do_test_merge = SpatialPixelsDataFrame(points = s_do_test_merge, data = s_do_test_merge@data)

# save new surface
summary(s_do_test_merge)
writeGDAL(s_do_test_merge["var1.pred"], fname = paste0(temp_wd, "sum_do_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_do_test = raster(paste0(temp_wd, "sum_do_test.tif"))
mean_sum_do_test = writeRaster(raster::mask(raster::crop(mean_sum_do_test, habitat_test), habitat_test),
                               file = file.path(test_wd, "Environmental/Sum_DO.asc"), format = "ascii",
                               overwrite = T)

rm(list = c("s_do_test_merge", "s_do_test_par", "mean_sum_do_test", "cl"))
showConnections()
