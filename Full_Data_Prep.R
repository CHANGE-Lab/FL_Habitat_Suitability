#### SET-UP ####
# This script was used to prepare data for Courtney Stuart's first MSc chapter 
# in the lab of Dr. Stephanie Green at the University of Alberta (2019-2021). 
# Data are specific to southern Florida and include: benthic habitat 
# classifications, bathymetric and topographic surfaces, georeferenced reef fish
# occurrence records, and bottom water conditions. These data were used for habitat
# suitability modeling, using penalized logistic regression and maximum entropy
# techniques.

# AUTHOR: Courtney Stuart
# DATE: May 13, 2021

#### TO USE THIS FILE ####
# This is script 1 of 4 in Courtney's data analysis pipeline

# working directories
setwd("Z:/Courtney/Stuart_MSc_Ch1/") # main project folder

# data directories 
temp_wd = "Z:/Courtney/Stuart_MSc_Ch1/Temporary/" # temporary files
source_wd = "Z:/Courtney/Stuart_MSc_Ch1/Source_Data/" # source data
dem_wd = "Z:/Courtney/Stuart_MSc_Ch1/Source_Data/Job606638_ncei_nintharcsec_dem/" # DEMs
csv_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Data/" # for writing tabular data
fish_wd = "Z:/Courtney/Stuart_MSc_Ch1/Species_Occurrence/" # for fish data
spatial_wd = "Z:/Courtney/Stuart_MSc_Ch1/Spatial_Predictors/" # for spatial predictor rasters
gis_wd = "Z:/Courtney/Stuart_MSc_Ch1/GIS_Files/" # for any GIS shapefiles

# libraries
library(easypackages)
libraries("conflicted", "rgdal", "gdalUtils", "raster", "sp", "sf", "tmap", "dplyr", 
          "lwgeom", "rgeos", "cleangeo", "tidyverse", "stars", "fasterize", 
          "PNWColors", "spex", "igraph", "spatialEco")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("xlim", "spex")

# change where large temp rasters are saved
rasterOptions(tmpdir = "Z:/Courtney/Stuart_MSc_Ch1/Temp/")

# save PROJ.4 string for standard projection (ESPG:26958 NAD 83/Florida East) 
# and source LAT/LON data (EPSG:4326 WGS 84/World Geodetic System 1984)
my_crs = CRS("+init=epsg:26958")
gcs = CRS("+init=epsg:4326")
us_ft = CRS("+init=epsg:2236")

#### HABITAT & DEPTH DATA ####
# read in park/marine sanctuary polygons
# Florida Keys National Marine Sanctuary (FKNMS) shapefile 
# (https://sanctuaries.noaa.gov/library/imast_gis.html)

fknms = st_read(dsn = paste0(source_wd, "Parks/fknms_py.shp")) %>%
  st_transform(., my_crs)
compareCRS(fknms, my_crs) # check projection
st_is_valid(fknms, reason = T) # check geometry

# National Park Service shapefile (for Biscayne National Park (BNP))
# (https://public-nps.opendata.arcgis.com/datasets/nps-boundary-1/data)
nps = st_read(dsn = paste0(source_wd, "Parks/NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service.shp")) %>%
  st_transform(., my_crs)
bnp = nps[nps$UNIT_NAME == "Biscayne National Park",] # extract BNP
compareCRS(bnp, my_crs)
st_is_valid(bnp, reason = T)

# how do the shapefiles line up
tm_shape(bnp) +
  tm_fill("navy") +
  tm_shape(fknms) +
  tm_fill("gray")

# there are gaps at the interior seams where FKNMS and BNP should meet, here is 
# an arbitrary polygon that I constructed in a GIS to fill the interior gaps 
# while still respecting the outer boundaries of the two parks
fill_gaps = data.frame(
  lon = c(291433.147, 292034.632, 292232.734, 292246.008, 285348.704, 279047.242,
          275554.735, 275369.526, 275054.142, 274774.742, 274393.741, 272105.355,
          270054.830, 269509.126, 268907.462, 268944.004, 269174.427, 270765.899,
          272195.446, 274481.451, 275137.619, 275518.619, 276132.454, 279709.628,
          284056.211, 290049.036, 290684.037),
  lat = c(145564.616, 146169.449, 146170.187, 132196.017, 103890.571, 110144.275, 
          110223.650, 111678.862, 112384.771, 112841.972, 113121.373, 114417.569, 
          114864.055, 115029.419, 114478.159, 114673.867, 115390.709, 115045.427,
          114723.957, 113390.455, 113030.621, 112247.452, 111570.118, 111675.951, 
          107164.796, 137644.857, 145185.497)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = my_crs) %>%
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# take a look now
tm_shape(bnp) +
  tm_fill("navy") +
  tm_shape(fknms) +
  tm_fill("gray") +
  tm_shape(fill_gaps) +
  tm_fill("orange")

# union the parks and gap fill data to produce one complete polygon
parks_union = st_union(fknms, bnp)
parks_union = st_union(parks_union, fill_gaps)

tm_shape(parks_union) +
  tm_fill("navy") 

# have an outline now of the parks' outer borders, so remove temporary data
rm(list = c("fknms", "bnp", "nps", "fill_gaps"))

# load habitat data from Unified Reef Tract Map (URM)
# (https://myfwc.com/research/gis/regional-projects/unified-reef-map/)
require(rgdal)
fgdb = paste0(source_wd, "Unified_Reef_Map/FWC_UnifiedFloridaReefMap_v2.0.gdb")

# List all feature classes in file geodatabase
subset(ogrDrivers(), grepl("GDB", name)) 
fc_list = ogrListLayers(fgdb) 
print(fc_list) 

# Read in URM feature class
reef_map = st_read(dsn = fgdb,layer = "UnifiedReefMap") %>%
  filter(!st_is_empty(.)) %>%
  st_transform(., my_crs)
compareCRS(reef_map, my_crs) # check projection to be sure

# check benthic habitat classes
unique(reef_map$ClassLv1) 

# create IDs because MaxEnt requires categorical variables to be defined 
# numerically, not with words
ClassLv1 = unique(reef_map$ClassLv1) # list categories
# create a data.frame of IDs
ClassLv1_df = data.frame(ID = 1:length(ClassLv1), ClassLv1 = ClassLv1) 
# match IDs
reef_map$ClassLv1_ID = ClassLv1_df$ID[match(reef_map$ClassLv1, ClassLv1_df$ClassLv1)] 
unique(reef_map$ClassLv1_ID) # did it work?  
# save for later
write.csv(ClassLv1_df, paste0(csv_wd, "URM_ClassLv1_IDs.csv"), row.names = F)

# clip unified reef map data to parks
reef_clip = st_intersection(st_make_valid(reef_map), parks_union)

# create a palette for plotting benthic habitat classes
pal = pnw_palette("Bay", 14, type = "continuous") 

tm_shape(reef_clip, projection = my_crs) +
  tm_fill("ClassLv1_ID", palette = pal, style = "cat") 

# supplementary shoreline mangrove habitat data
# (https://geodata.myfwc.com/datasets/mangrove-habitat-in-florida-1/explore)
mg_shore = 
  st_read(dsn = paste0(source_wd, 
                       "Mangrove_Habitat/Mangrove_Habitat_in_Florida.shp")) %>%
  filter(!st_is_empty(.)) %>%
  st_transform(., my_crs) %>%
  st_cast(., "MULTIPOLYGON")
compareCRS(mg_shore, my_crs)

# there is a small gap between the mangrove data and the reef tract map along 
# the mainland coast, and these missing areas are important coastline mangroves 
# (visible in satellite imagery) add a 10 m buffer around the mangrove data and 
# then use gDifference to keep only the non-overlapping areas (AKA, respect the 
# boundaries of the reef map).
mg_buff = st_buffer(mg_shore, dist = 10)

# now keep only the non-overlapping regions of the mangrove data 
# (FYI: time-consuming step)
mg_clip = rgeos::gDifference(as_Spatial(st_intersection(mg_buff, parks_union)),
                             as_Spatial(reef_clip)) %>% st_as_sf()

# double check the ID assigned to mangroves from the reef map and add it to the
# mangrove data
ClassLv1_df
mg_clip$ClassLv1 = rep(as.character("Mangrove"), nrow(mg_clip))
mg_clip$ClassLv1_ID = rep(as.integer(11), nrow(mg_clip))

# save mangrove and reef data for mapping
st_write(mg_clip, dsn = paste0(gis_wd, "Mangrove_Habitat.shp"), 
         driver = "ESRI Shapefile", append = F)
st_write(reef_clip, dsn = paste0(gis_wd, "Reef_Map_Habitat.shp"), 
         driver = "ESRI Shapefile", append = F)

#  check geometry for both habitat layers
head(mg_clip, 1) # sf column: geometry 
head(reef_clip, 1) # sf column: Shape
reef_clip = reef_clip %>% rename(geometry = Shape)
head(reef_clip, 1) # fixed

# combining the reef and mangrove habitat layers
hab_clip = rbind(mg_clip, select(reef_clip, geometry, ClassLv1, ClassLv1_ID)) %>%
  filter(!ClassLv1 %in% c("Land", "Not Classified")) %>% # only want real benthic habitats
  st_cast("MULTIPOLYGON")

tm_shape(hab_clip, projection = my_crs) +
  tm_fill("ClassLv1_ID", palette = pal, style = "cat") 

# template raster for habitat w/ 5 x 5 m res and standard CRS
ras = raster(ext = extent(hab_clip), res = c(5,5), crs = my_crs)

# convert benthic habitat polygons to a temporary raster
hab_ras = writeRaster((fasterize(hab_clip, ras, field = "ClassLv1_ID", fun = "max")), 
                      file = file.path(temp_wd, "hab_ras.tif"), format = "GTiff", 
                      overwrite = T)

# pause on habitat for now so we can bring in the depth data

# load Continuously Updated DEM - 1/9 Arc-Second Resolution Bathymetric-Topographic 
# Tiles  from NOAA-NCEI
# (https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ngdc.mgg.dem:999919)

setwd(dem_wd) # temporarily set WD to source DEMs folder (downloaded tiles)
dems = c('Job606638_ncei_nintharcsec_dem_002_001.tif', 
         'Job606638_ncei_nintharcsec_dem_002_000.tif',
         'Job606638_ncei_nintharcsec_dem_001_001.tif',
         'Job606638_ncei_nintharcsec_dem_001_000.tif',
         'Job606638_ncei_nintharcsec_dem_000_000.tif')
#'Job606638_ncei_nintharcsec_dem_000_001.tif' # this tile is for Gulf side, we don't need it

# overall extent of bathy-topo tiles
e = extent(618778.56, 988742.88, 61706.64, 640446.24)
template = raster(e) # use extent to create empty template raster

# set CRS of template using info in metadata of tiles
# Horizontal CS: EPSG:3512 NAD83(NSRS2007) / Florida East (US foot)
template_crs = CRS("+init=epsg:3512")
crs(template) = template_crs

# GDAL version 2.1.3
# write the empty template raster 
writeRaster(template, file = file.path(temp_wd, "dem_mosaic.tif"), 
            format = "GTiff", overwrite = T)

# mosaic tiles and save them to the empty template raster
gdalUtils::mosaic_rasters(gdalfile = dems, 
                          dst_dataset = file.path(temp_wd, "dem_mosaic.tif"),
                          of = "GTiff", overwrite = T)

# re-project DEM mosaic to standard CRS and use resolution of 5 x 5 m
gdalUtils::gdalwarp(srcfile = paste0(temp_wd, "dem_mosaic.tif"), # mosaic DEM
                    dstfile = paste0(temp_wd, "dem_mosaic_proj.tif"), # destination projected DEM
                    s_srs = template_crs, # source crs
                    t_srs = my_crs, # destination crs
                    tr = c(5, 5), # destination cell resolution in meters
                    overwrite = T,
                    verbose = T)

# re-load the projected DEM
dem_proj = raster(paste0(temp_wd, "dem_mosaic_proj.tif"))

# the projected DEM still stores depth data in units US feet, convert to meters
dem_proj = dem_proj * 0.304801  # (1 ft = 0.304801 m)
dem_proj = raster::setMinMax(dem_proj) # calculate and save min-max values
dem_proj # check attributes

# remove any cells that have a value > 1 m (aka sites above ground)
# keep 0 to 1 m because these areas may be partially inundated by tides
dem_proj1 = dem_proj
dem_proj1[dem_proj1 > 1] = NA

# plot the depth data using new palette
pal2 = pnw_palette("Winter", n = 100, "continuous")
plot(dem_proj1, col = pal2)

# the DEM data will be used in ArcGIS with the Benthic Terrain Modeler
# extension to calculate seafloor surface morphometrics -- the rugosity and slope
# tools will require the FULL DEM (projected mosaic) rather than the clipped 
# DEMs because these tools do NOT respect ArcGIS' environment settings
writeRaster(dem_proj1, paste0(dem_wd, "Full_Mosaic_DEM.tif"), format = "GTiff",
            overwrite = T)

# remove temporary products
rm(list = c("dems", "e", "template", "template_crs")) 

# back to normal working directory
setwd("Z:/Courtney/Stuart_MSc_Ch1/")

# clip habitat data to DEM extent first, then clip DEM back to habitat; 
# mask and match extents
hab_crop_dem = writeRaster((raster::crop(hab_ras, dem_proj1)), 
                           file = file.path(temp_wd, "hab_crop_dem.tif"), 
                           format = "GTiff", overwrite = T)
dem_crop_hab = writeRaster((raster::crop(dem_proj1, hab_ras)), 
                           file = file.path(temp_wd, "dem_crop_hab.tif"), 
                           format = "GTiff", overwrite = T)
extent(hab_crop_dem) = extent(dem_crop_hab) # match extents
habitat = raster::mask(hab_crop_dem, dem_crop_hab) # habitat data within DEM extent
depth = raster::mask(dem_crop_hab, hab_crop_dem) # depth data matching habitat extent 
extent(habitat) = extent(depth)
compareRaster(depth, habitat, extent = T, res = T, crs = T, rowcol = T)

# they should line up perfectly 
tm_shape(depth) +
  tm_raster(palette = pal2, style = "cont") +
  tm_shape(habitat) +
  tm_raster(palette = pal, n = 14) +
  tm_layout(frame = F)

# there is a lot of blank space (NAs) in the lower left portion of the rasters,
# clip both raster data layers to remove some of this blank space (save storage)
clip = data.frame(
  lon = c(148729.419, 148729.419, 292477.832, 292477.832), # bbox values in meters
  lat = c(20120.419, 148774.322, 148774.322, 20120.419)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = my_crs) %>%
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

depth = raster::mask(raster::crop(depth, extent(clip)), clip)
crs(depth) = my_crs
habitat = raster::mask(raster::crop(habitat, extent(clip)), clip)
crs(habitat) = my_crs

# write out data - need to use ArcGIS' Benthic Terrain Modeler extension to 
# calculate seascape topography metrics from the depth raster (*remember, the
# slope and rugosity tools will require the FULL DEM in the source DEM folder*)
# MaxEnt requires ascii format
writeRaster(depth, file = file.path(spatial_wd, "Depth.asc"), 
            format = "ascii", overwrite = T)
writeRaster(habitat, file = file.path(spatial_wd, "Habitat.asc"), 
            format = "ascii", overwrite = T)

# create a simple, constant-value raster from the habitat raster for clipping 
writeRaster(((habitat * 0) + 1), file = file.path(temp_wd, "Constant_Hab_Ras.tif"), 
            format = "GTiff", overwrite = T)

# define function for converting raster to polygons 
# credit: John Baumgartner (johnbaums / polygonize.R on github)
polygonize <- function(srcfile, dstfile, mask, ogr_format, fieldname, band, connect8 = FALSE) {
  options <- if(isTRUE(connect8)) 'CONNECT8=8' else character(0)
  contour_options <- character(0)
  if(missing(mask)) mask <- character(0)
  .Call("_sf_CPL_polygonize", PACKAGE = "sf", srcfile, mask, 
        'GTiff', ogr_format, layer = dstfile, options, 0, #iPixValField, 
        contour_options, use_contours = FALSE, use_integer = TRUE)
}

# now convert the constant-value habitat raster to simple polygon
polygonize(srcfile = paste0(temp_wd, "Constant_Hab_Ras.tif"), 
           dstfile = file.path(temp_wd, "Constant_Hab_Shp.shp"), 
           ogr_format ='ESRI Shapefile', connect8 = F)

# This simple, constant-value habitat raster converted to a shapefile 
# represents our study domain
domain = st_read(dsn = paste0(temp_wd, "Constant_Hab_Shp.shp")) %>%
  st_transform(., my_crs)
domain = domain[domain$Value == 1,] # selecting only areas that had habitat data
plot(domain)
# saving domain area as a shapefile
st_write(domain, dsn = paste0(gis_wd, "Study_Domain.shp"), append = F)

# lastly, create "proximity to mangrove" surface by calculating Euclidean distance 
# from each cell to the nearest mangrove cell
ClassLv1_df # again finding mangrove ID --> 11 is what we want to keep
m = c(1, NA, 2, NA, 3, NA, 4, NA, 5, NA, 6, NA, 7, NA, 8, NA,
      9, NA, 10, NA, 11, 11, 12, NA, 13, NA, 14, NA) # cell IDs 1:10 & 12-14 become NA
mat = matrix(m, ncol = 2, byrow = T)
mg_dist = writeRaster(raster::mask(raster::crop(gridDistance(habitat, origin = 11),
                                                habitat), habitat),
                      file = file.path(spatial_wd, "Mangrove_Dist.asc"),
                      format = "ascii", overwrite = T)

##### FISH DATA #####
# reef visual census data in the study domain
library(rvc)
library(snow)
require(rvc)
rvc = getSampleData(years = c(2014, 2016, 2018), regions = "FLA KEYS")
head(rvc, 1)

# convert the fork length column (LEN) to total length (TOT_LEN) using FishBase 
# length-length conversion for gray snapper specifically (TL = 0 + 1.049 x FL)
# will repeat later for bluestriped grunt using their specific FL-TL conversion
rvc = rvc %>% mutate(LG_TOT_LEN = (LEN*1.049)) 

# how many unique sites were sampled in 2014, 2016, and 2018?
rvc_sites = rvc %>% distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, 
                             STATION_NR, LAT_DEGREES, LON_DEGREES, 
                             .keep_all = T)

#### gray snapper ####
# start with presence and absence records specific to the adult gray snapper 
# stage (> 24.71 cm TL (size at maturity))
lg_adult = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, 
           LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% # group by primary and second-stage sample units
  filter(SPECIES_CD == "LUT GRIS" & LG_TOT_LEN > 24.71) %>% 
  summarize(N = sum(NUM)) %>% # was an adult gray snapper ever seen at this second-stage sample unit (SSU)?
  ungroup() %>%
  distinct()
lg_adult$LIFE_STAGE = rep("ADULT", nrow(lg_adult)) # specify that these are all adult records
lg_adult$PRES = ifelse(lg_adult$N > 0, 1, 0) # if at least one adult was seen at the SSU, it's a presence [1], else an absence [0]
lg_adult$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_adult)) # source column for when RVC and MG data are compiled

# now the inferred absences sites (sites where either no gray snapper were seen 
# or only those of another age class were seen)
lg_adult_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, 
           LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "LUT GRIS" & !LG_TOT_LEN > 24.71) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 b/c these are adult absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, 
           LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
lg_adult_abs$LIFE_STAGE = rep("ADULT", nrow(lg_adult_abs))
lg_adult_abs$PRES = ifelse(lg_adult_abs$N > 0, 1, 0)
lg_adult_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_adult_abs))

# combining adult absence and presence data and using distinct because some of 
# the absence sites might be shared across both the adult stage-specific data 
# and the inferred absence data
lg_adult_rvc = rbind(lg_adult, lg_adult_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, 
           LAT_DEGREES, LON_DEGREES, .keep_all = T)

# repeat for the subadult gray snapper stage (9.51 cm <= TOT LEN <= 24.71 cm 
# (between size at 1 YR and size at maturation))
lg_subadult = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, 
           LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% 
  filter(SPECIES_CD == "LUT GRIS" & (LG_TOT_LEN >= 9.51 & LG_TOT_LEN <= 24.71)) %>% 
  summarize(N = sum(NUM)) %>% # was a subadult gray snapper ever seen at this SSU?
  ungroup() %>%
  distinct()
lg_subadult$LIFE_STAGE = rep("SUBADULT", nrow(lg_subadult)) 
lg_subadult$PRES = ifelse(lg_subadult$N > 0, 1, 0)
lg_subadult$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_subadult))

# now the inferred absences sites (sites where either no gray snapper were seen
# or only those of another age class were seen)
lg_subadult_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, 
           LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "LUT GRIS" & !(LG_TOT_LEN >= 9.51 & LG_TOT_LEN <= 24.71)) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 b/c these are subadult absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, 
           LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
lg_subadult_abs$LIFE_STAGE = rep("SUBADULT", nrow(lg_subadult_abs))
lg_subadult_abs$PRES = ifelse(lg_subadult_abs$N > 0, 1, 0)
lg_subadult_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_subadult_abs))

# combining subadult absence and presence data and using distinct because some 
# of the absence sites might be shared across both the subadult stage-specific 
# data and the inferred absence data
lg_subadult_rvc = rbind(lg_subadult, lg_subadult_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR,
           LAT_DEGREES, LON_DEGREES, .keep_all = T)

# finally, repeat for the juvenile gray snapper stage (TOT LEN < 9.51 cm 
# (smaller than size at 1 YR), but not equal to 0!)
lg_juvenile = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR,
           LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% 
  filter(SPECIES_CD == "LUT GRIS" & (LG_TOT_LEN != 0 & LG_TOT_LEN < 9.51)) %>% 
  summarize(N = sum(NUM)) %>% # was a juvenile gray snapper ever seen at this SSU?
  ungroup() %>%
  distinct()
lg_juvenile$LIFE_STAGE = rep("JUVENILE", nrow(lg_juvenile)) 
lg_juvenile$PRES = ifelse(lg_juvenile$N > 0, 1, 0)
lg_juvenile$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_juvenile))

# now the inferred absences sites (sites where either no gray snapper were seen 
# or only those of another age class were seen)
lg_juvenile_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, 
           LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "LUT GRIS" & !(LG_TOT_LEN != 0 & LG_TOT_LEN < 9.51)) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 b/c these are juvenile absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, 
           LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
lg_juvenile_abs$LIFE_STAGE = rep("JUVENILE", nrow(lg_juvenile_abs))
lg_juvenile_abs$PRES = ifelse(lg_juvenile_abs$N > 0, 1, 0)
lg_juvenile_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_juvenile_abs))

# combining juvenile absence and presence data and using distinct because some 
# of the absence sites might be shared across both the juvenile stage-specific 
# data and the inferred absence data
lg_juvenile_rvc = rbind(lg_juvenile, lg_juvenile_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR,
           LAT_DEGREES, LON_DEGREES, .keep_all = T)

# compiling all RVC datasets just in case
lg_rvc = do.call("rbind", list(lg_adult_rvc, lg_subadult_rvc, lg_juvenile_rvc))

# and getting rid of all intermediate tables
rm(list = c("lg_adult", "lg_adult_abs", "lg_subadult", "lg_subadult_abs", 
            "lg_juvenile", "lg_juvenile_abs"))

# RVC data done, now time for Mangrove Visual Survey (MVS) data
mg = read.csv(paste0(source_wd, "Fish_Data/Shoreline_Mangrove_Surveys/DATA_Serafy_BB_MANGROVE_FISH_DATA_1998W_2019W.csv"), stringsAsFactors = F)
head(mg, 1)

# filter out years of interest (2014, 2016, 2018)
# convert min, avg, and max total length values from inches to centimeters 
# (1 in = 2.54 cm), rename lat/long and species name columns to match RVC dataset 
# (to join tables later)
mg = mg %>%
  filter(YR %in% c(2014, 2016, 2018)) %>%
  mutate(MIN_TOT_LEN = (as.numeric(MIN_IN))*2.54) %>%
  mutate(AVE_TOT_LEN = (as.numeric(AVE_IN))*2.54) %>%
  mutate(MAX_TOT_LEN = (as.numeric(MAX_IN))*2.54) %>%
  rename(YEAR = YR, MONTH = MO, DAY = DY, SPECIES_CD = SP, LAT_DEGREES = LAT, 
         LON_DEGREES = LON)

# first find absences shared across life stages by isolating sites where no gray
# snapper were observed
mg_sites = mg %>% distinct(Site, LON_DEGREES, LAT_DEGREES)

# now isolate gray snapper presence sites
lg_pres_sites = mg %>%
  filter(SPECIES_CD == "LUT_GRI") %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# these absence sites apply to all gray snapper
lg_abs = mg %>%
  filter(!Site %in% lg_pres_sites$Site) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0, "SPECIES_CD" = "LUT GRIS")
lg_abs$LIFE_STAGE = rep("ADULT", nrow(lg_abs)) # life stage column (temporarily ADULT)
lg_abs$PRES = ifelse(lg_abs$N > 0, 1, 0) # was a gray snapper ever seen here?
lg_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(lg_abs)) # source column for when data are compiled

# now considering only sites where gray snapper were observed
lg = mg %>% filter(SPECIES_CD == "LUT_GRI")
view(lg)

# need to parse out life history stages using minimum and maximum body lengths --
# looks like there are "." entries for min and max lengths when only one fish
# was caught at a site, change that to the actual length values using
# data.frame[row_number, column_number] = new_value 
# (might be an easier way but I'm not sure how)
lg[46, 11] = as.numeric(3.00)
lg[46, 13] = as.numeric(3.00)
lg[50, 11] = as.numeric(4.00)
lg[50, 13] = as.numeric(4.00)
lg[67, 11] = as.numeric(4.00)
lg[67, 13] = as.numeric(4.00)
lg[69, 11] = as.numeric(11.00)
lg[69, 13] = as.numeric(11.00)
lg[70, 11] = as.numeric(13.00)
lg[70, 13] = as.numeric(13.00)
lg[98, 11] = as.numeric(1.00)
lg[98, 13] = as.numeric(1.00)
lg[101, 11] = as.numeric(6.00)
lg[101, 13] = as.numeric(6.00)
lg[250, 13] = as.numeric(8.00)


view(lg) # looks better, now fix the min, ave, and max total length columns again
lg = lg %>% 
  mutate(MIN_TOT_LEN = (as.numeric(MIN_IN)*2.54)) %>%
  mutate(AVE_TOT_LEN = (as.numeric(AVE_IN)*2.54)) %>%
  mutate(MAX_TOT_LEN = (as.numeric(MAX_IN)*2.54))
summary(lg$MIN_TOT_LEN) # check for NAs
summary(lg$MAX_TOT_LEN)

# adult stage-specific records where EITHER min or max total lengths exceed the 
# size at maturity (> 24.71 cm TL)
lg_adult = lg %>%   
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(MIN_TOT_LEN > 24.71 | MAX_TOT_LEN > 24.71) %>% 
  summarize(N = sum(as.numeric(NO))) %>%
  mutate(SPECIES_CD = "LUT GRIS") %>% # to match RVC
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE) %>%
  ungroup()
lg_adult$LIFE_STAGE = rep("ADULT",nrow(lg_adult))
lg_adult$PRES = ifelse(lg_adult$N>0, 1, 0)
lg_adult$SOURCE = rep("MANGROVE VISUAL SURVEY",nrow(lg_adult))

# now the inferred absences sites (sites where either no gray snapper were seen 
# or only those of another age class were seen)
lg_adult_abs = lg %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN > 24.71 | MAX_TOT_LEN > 24.71)) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT adults
  ungroup()
lg_adult_abs$LIFE_STAGE = rep("ADULT", nrow(lg_adult_abs))
lg_adult_abs$PRES = ifelse(lg_adult_abs$N > 0, 1, 0) 
lg_adult_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(lg_adult_abs))

# combining the adult stage-specific data and inferred absence data, 
# again using distinct to remove repeating values
lg_adult_mg = do.call("rbind", list(lg_adult, lg_adult_abs, lg_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)


# repeat for the subadult gray snapper stage (9.51 cm <= TOT LEN <= 24.71 cm 
# (between size at 1 YR and size at maturation)). remember first to change the 
# lifestage column of the shared absence dataframe
lg_abs$LIFE_STAGE = rep("SUBADULT", nrow(lg_abs))

lg_subadult = lg %>%   
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter((MIN_TOT_LEN >= 9.51 & MIN_TOT_LEN <= 24.71) | (MAX_TOT_LEN >= 9.51 & MAX_TOT_LEN <= 24.71)) %>% # if either of these are true, then a subadult was present
  summarize(N = sum(as.numeric(NO))) %>%
  mutate(SPECIES_CD = "LUT GRIS") %>% # to match RVC
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE) %>%
  ungroup()
lg_subadult$LIFE_STAGE = rep("SUBADULT", nrow(lg_subadult))
lg_subadult$PRES = ifelse(lg_subadult$N > 0, 1, 0)
lg_subadult$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(lg_subadult))

# # now the inferred absences sites (sites where either no gray snapper were 
# seen or only those of another age class were seen)
lg_subadult_abs = lg %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN >= 9.51 & MIN_TOT_LEN <= 24.71) | (MAX_TOT_LEN >= 9.51 & MAX_TOT_LEN <= 24.71)) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT sub-adults
  ungroup()
lg_subadult_abs$LIFE_STAGE = rep("SUBADULT", nrow(lg_subadult_abs))
lg_subadult_abs$PRES = ifelse(lg_subadult_abs$N > 0, 1, 0) 
lg_subadult_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(lg_subadult_abs))

# combining the subadult stage-specific data and inferred absence data, 
# again using distinct to remove repeating values
lg_subadult_mg = do.call("rbind", list(lg_subadult, lg_subadult_abs, lg_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# finally, repeat for juvenile gray snapper stage (TOT LEN < 9.51 cm (smaller 
# than size at 1 YR), but not equal to 0!)
lg_abs$LIFE_STAGE = rep("JUVENILE", nrow(lg_abs))

lg_juvenile = lg %>%   
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(MIN_TOT_LEN != 0 & (MIN_TOT_LEN < 9.51 | MAX_TOT_LEN < 9.51)) %>% # if either of these are true, then a juvenile was present
  summarize(N = sum(as.numeric(NO))) %>%
  mutate(SPECIES_CD = "LUT GRIS") %>% # to match RVC
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE) %>%
  ungroup()
lg_juvenile$LIFE_STAGE = rep("JUVENILE", nrow(lg_juvenile))
lg_juvenile$PRES = ifelse(lg_juvenile$N > 0, 1, 0)
lg_juvenile$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(lg_juvenile))

# # now the inferred absences sites (sites where either no gray snapper were 
# seen or only those of another age class were seen)
lg_juvenile_abs = lg %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN != 0 & (MIN_TOT_LEN < 9.51 | MAX_TOT_LEN < 9.51))) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT juveniles
  ungroup()
lg_juvenile_abs$LIFE_STAGE = rep("JUVENILE", nrow(lg_juvenile_abs))
lg_juvenile_abs$PRES = ifelse(lg_juvenile_abs$N > 0, 1, 0) 
lg_juvenile_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(lg_juvenile_abs))

# combining the juvenile stage-specific data and inferred absence data, again 
# using distinct to remove repeating values
lg_juvenile_mg = do.call("rbind", list(lg_juvenile, lg_juvenile_abs, lg_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# compiling all MG datasets just in case
lg_mg = do.call("rbind", list(lg_adult_mg, lg_subadult_mg, lg_juvenile_mg))

# and getting rid of all intermediate tables
rm(list = c("lg_adult", "lg_adult_abs", "lg_subadult", "lg_subadult_abs",
            "lg_juvenile", "lg_juvenile_abs", "lg", "lg_abs", "lg_pres_sites"))

# now combining the two data sources - RVC and MG - and clip to study domain
# adults
lg_adult_clip = full_join(lg_adult_rvc, lg_adult_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) # re-project 
lg_adult_clip = lg_adult_clip[domain, ]
compareCRS(lg_adult_clip, my_crs) # just to be sure

rm(list = c("lg_adult_mg", "lg_adult_rvc"))

# subadults
lg_subadult_clip = full_join(lg_subadult_rvc, lg_subadult_mg)  %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) # re-project 
lg_subadult_clip = lg_subadult_clip[domain, ] # clip
compareCRS(lg_subadult_clip, my_crs) # just to be sure

rm(list = c("lg_subadult_mg", "lg_subadult_rvc"))

# juveniles
lg_juvenile_clip = full_join(lg_juvenile_rvc, lg_juvenile_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) # re-project 
lg_juvenile_clip = lg_juvenile_clip[domain, ] # clip
compareCRS(lg_juvenile_clip, my_crs) # just to be sure

rm(list = c("lg_juvenile_mg", "lg_juvenile_rvc"))

# save presence-absence data for logistic regression modeling
st_write(lg_adult_clip %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(fish_wd, "Presence_Absence/Adult_Gray_Snapper_PA_Full.csv"), 
         append = FALSE)

st_write(lg_subadult_clip %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(fish_wd, "Presence_Absence/Subadult_Gray_Snapper_PA_Full.csv"), 
         append = FALSE)

st_write(lg_juvenile_clip %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(fish_wd, "Presence_Absence/Juvenile_Gray_Snapper_PA_Full.csv"), 
         append = FALSE)

# now filtering out only the presence records for MaxEnt modeling
# MaxEnt only needs species name, lon, lat (in that order)
st_write(lg_adult_clip %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), 
         dsn = paste0(fish_wd, "Presence_Only/Adult_Gray_Snapper_PO_Full.csv"), 
         append = FALSE)

st_write(lg_subadult_clip %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), 
         dsn = paste0(fish_wd, "Presence_Only/Subadult_Gray_Snapper_PO_Full.csv"), 
         append = FALSE)

st_write(lg_juvenile_clip %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), 
         dsn = paste0(fish_wd, "Presence_Only/Juvenile_Gray_Snapper_PO_Full.csv"), 
         append = FALSE)

#### bluestriped grunts ####

# first convert fork length to total length using the bluestriped grunt
# length-length conversion from FishBase  TL = 0 + 1.034 x FL
rvc = rvc %>% mutate(HS_TOT_LEN = (LEN*1.034)) 

# start with presence and absence records specific to the adult bluestriped grunt 
# stage (> 25.33 cm TL (size at maturity))
hs_adult = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% # group by both primary and second-stage sample units
  filter(SPECIES_CD == "HAE SCIU" & HS_TOT_LEN > 25.33) %>% 
  summarize(N = sum(NUM)) %>% # was an adult bluestriped grunt ever seen at this SSU? 
  ungroup() %>%
  distinct()
hs_adult$LIFE_STAGE = rep("ADULT", nrow(hs_adult)) # specify that these are all adult records
hs_adult$PRES = ifelse(hs_adult$N > 0, 1, 0) # if at least one adult was seen at the SSU, it's a presence [1], else an absence [0]
hs_adult$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_adult)) # source column for when RVC and MG data are compiled

# now the inferred absences sites (sites where either no gray bluestriped grunts 
# were seen or only those of another age class were seen)
hs_adult_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "HAE SCIU" & !HS_TOT_LEN > 25.33) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 b/c these are adult absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
hs_adult_abs$LIFE_STAGE = rep("ADULT", nrow(hs_adult_abs))
hs_adult_abs$PRES = ifelse(hs_adult_abs$N > 0, 1, 0)
hs_adult_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_adult_abs))

# combining adult absence and presence data and using distinct because some of 
# the absence sites might be shared across both the adult stage-specific data 
# and the inferred absence data
hs_adult_rvc = rbind(hs_adult, hs_adult_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, .keep_all = T)

# repeat for the subadult bluestriped grunt stage (11.90 cm <= TOT LEN <= 25.33 cm 
# (between size at 1 YR and size at maturation))
hs_subadult = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% 
  filter(SPECIES_CD == "HAE SCIU" & (HS_TOT_LEN >= 11.90 & HS_TOT_LEN <= 25.33)) %>% 
  summarize(N = sum(NUM)) %>% # was a subadult bluestriped grunt ever seen at this SSU?
  ungroup() %>%
  distinct()
hs_subadult$LIFE_STAGE = rep("SUBADULT", nrow(hs_subadult)) 
hs_subadult$PRES = ifelse(hs_subadult$N > 0, 1, 0)
hs_subadult$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_subadult))

# now the inferred absences sites (sites where either no bluestriped grunt were 
# seen or only those of another age class were seen)
hs_subadult_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "HAE SCIU" & !(HS_TOT_LEN >= 11.90 & HS_TOT_LEN <= 25.33)) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 b/c these are subadult absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
hs_subadult_abs$LIFE_STAGE = rep("SUBADULT", nrow(hs_subadult_abs))
hs_subadult_abs$PRES = ifelse(hs_subadult_abs$N > 0, 1, 0)
hs_subadult_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_subadult_abs))

# combining subadult absence and presence data and using distinct because some 
# of the absence sites might be shared across both the subadult stage-specific 
# data and the inferred absence data
hs_subadult_rvc = rbind(hs_subadult, hs_subadult_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, .keep_all = T)

# finally, repeat for the juvenile bluestriped grunt stage (TOT LEN < 11.90 cm 
# (smaller than size at 1 YR), but not equal to 0!)
hs_juvenile = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% 
  filter(SPECIES_CD == "HAE SCIU" & (HS_TOT_LEN != 0 & HS_TOT_LEN < 11.90)) %>% 
  summarize(N = sum(NUM)) %>% # was a juvenile bluestriped grunt ever seen at this SSU?
  ungroup() %>%
  distinct()
hs_juvenile$LIFE_STAGE = rep("JUVENILE", nrow(hs_juvenile)) 
hs_juvenile$PRES = ifelse(hs_juvenile$N > 0, 1, 0)
hs_juvenile$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_juvenile))

# now the inferred absences sites (sites where either no bluestriped grunts were
# seen or only those of another age class were seen)
hs_juvenile_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "HAE SCIU" & !(HS_TOT_LEN != 0 & HS_TOT_LEN < 11.90)) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 b/c these are juvenile absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
hs_juvenile_abs$LIFE_STAGE = rep("JUVENILE", nrow(hs_juvenile_abs))
hs_juvenile_abs$PRES = ifelse(hs_juvenile_abs$N > 0, 1, 0)
hs_juvenile_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_juvenile_abs))

# combining juvenile absence and presence data and using distinct because some 
# of the absence sites might be shared across both the juvenile stage-specific 
# and the inferred absence data
hs_juvenile_rvc = rbind(hs_juvenile, hs_juvenile_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, .keep_all = T)


# compiling all RVC datasets just in case
hs_rvc = do.call("rbind", list(hs_adult_rvc, hs_subadult_rvc, hs_juvenile_rvc))

# and getting rid of all intermediate tables
rm(list = c("hs_adult", "hs_adult_abs", "hs_subadult", "hs_subadult_abs", 
            "hs_juvenile", "hs_juvenile_abs"))

# mangrove data
# bluestriped grunt presence sites
hs_pres_sites = mg %>%
  filter(SPECIES_CD == "HAE_SCI") %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# these absence sites apply to all bluestriped grunts
hs_abs = mg %>%
  filter(!Site %in% hs_pres_sites$Site) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0, "SPECIES_CD" = "HAE SCIU")
hs_abs$LIFE_STAGE = rep("ADULT", nrow(hs_abs)) # life stage column (temporarily ADULT)
hs_abs$PRES = ifelse(hs_abs$N > 0, 1, 0) # was a bluestriped grunt ever seen here?
hs_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(hs_abs)) # source column for when data are compiled

# now considering only sites where bluestriped grunt were observed
hs = mg %>% filter(SPECIES_CD == "HAE_SCI")
view(hs)

# for min and max lengths when only one fish was caught at a site, change that to the actual length values using
# data.frame[row_number, column_number] = new_value (might be an easier way but I'm not sure how)
hs[29, 11] = as.numeric(6.00)
hs[29, 13] = as.numeric(6.00)
hs[30, 11] = as.numeric(4.00)
hs[30, 13] = as.numeric(4.00)
hs[32, 11] = as.numeric(3.00)
hs[32, 13] = as.numeric(3.00)

# now fix the min, ave, and max total length columns again
hs = hs %>% 
  mutate(MIN_TOT_LEN = (as.numeric(MIN_IN)*2.54)) %>%
  mutate(AVE_TOT_LEN = (as.numeric(AVE_IN)*2.54)) %>%
  mutate(MAX_TOT_LEN = (as.numeric(MAX_IN)*2.54))
summary(hs$MIN_TOT_LEN) # check for NAs
summary(hs$MAX_TOT_LEN)


# adult stage-specific records where EITHER min or max total lengths exceed the
# size at maturity (> 24.71 cm TL)
hs_adult = hs %>%   
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(MIN_TOT_LEN > 25.33 | MAX_TOT_LEN > 25.33) %>% 
  summarize(N = sum(as.numeric(NO))) %>%
  mutate(SPECIES_CD = "HAE SCIU") %>% # to match RVC
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE) %>%
  ungroup()
hs_adult$LIFE_STAGE = rep("ADULT",nrow(hs_adult))
hs_adult$PRES = ifelse(hs_adult$N>0, 1, 0)
hs_adult$SOURCE = rep("MANGROVE VISUAL SURVEY",nrow(hs_adult))

# now the inferred absences sites (sites where either no bluestriped grunts were
# seen or only those of another age class were seen)
hs_adult_abs = hs %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN > 25.33 | MAX_TOT_LEN > 25.33)) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT adults
  ungroup()
hs_adult_abs$LIFE_STAGE = rep("ADULT", nrow(hs_adult_abs))
hs_adult_abs$PRES = ifelse(hs_adult_abs$N > 0, 1, 0) 
hs_adult_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(hs_adult_abs))

# combining the adult stage-specific data and inferred absence data, again using
# distinct to remove repeating values
hs_adult_mg = do.call("rbind", list(hs_adult, hs_adult_abs, hs_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# repeat for the subadult bluestriped grunt stage (11.90 cm <= TOT LEN <= 25.33 cm 
# (between size at 1 YR and size at maturation)) remember first to change the 
# lifestage column of the shared absence dataframe
hs_abs$LIFE_STAGE = rep("SUBADULT", nrow(hs_abs))

hs_subadult = hs %>%   
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter((MIN_TOT_LEN >= 11.90 & MIN_TOT_LEN <= 25.33) | (MAX_TOT_LEN >= 11.90 & MAX_TOT_LEN <= 25.33)) %>% # if either of these are true, then a subadult was present
  summarize(N = sum(as.numeric(NO))) %>%
  mutate(SPECIES_CD = "HAE SCIU") %>% # to match RVC
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE) %>%
  ungroup()
hs_subadult$LIFE_STAGE = rep("SUBADULT", nrow(hs_subadult))
hs_subadult$PRES = ifelse(hs_subadult$N > 0, 1, 0)
hs_subadult$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(hs_subadult))

# now the inferred absences sites (sites where either no bluestriped grunts were 
# seen or only those of another age class were seen)
hs_subadult_abs = hs %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN >= 11.90 & MIN_TOT_LEN <= 25.33) | (MAX_TOT_LEN >= 11.90 & MAX_TOT_LEN <= 25.33)) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 b/c these are NOT subadults
  ungroup()
hs_subadult_abs$LIFE_STAGE = rep("SUBADULT", nrow(hs_subadult_abs))
hs_subadult_abs$PRES = ifelse(hs_subadult_abs$N > 0, 1, 0) 
hs_subadult_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(hs_subadult_abs))

# combining the subadult stage-specific data and inferred absence data, a
# gain using distinct to remove repeating values
hs_subadult_mg = do.call("rbind", list(hs_subadult, hs_subadult_abs, hs_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# finally, repeat for juvenile bluestriped grunt stage (TOT LEN < 11.90 cm 
# (smaller than size at 1 YR), but not equal to 0!)
hs_abs$LIFE_STAGE = rep("JUVENILE", nrow(hs_abs))

hs_juvenile = hs %>%   
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(MIN_TOT_LEN != 0 & (MIN_TOT_LEN < 11.90 | MAX_TOT_LEN < 11.90)) %>% # if either of these are true, then a juvenile was present
  summarize(N = sum(as.numeric(NO))) %>%
  mutate(SPECIES_CD = "HAE SCIU") %>% # to match RVC
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE) %>%
  ungroup()
hs_juvenile$LIFE_STAGE = rep("JUVENILE", nrow(hs_juvenile))
hs_juvenile$PRES = ifelse(hs_juvenile$N > 0, 1, 0)
hs_juvenile$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(hs_juvenile))

# now the inferred absences sites (sites where either no bluestriped grunt were 
# seen or only those of another age class were seen)
hs_juvenile_abs = hs %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN != 0 & (MIN_TOT_LEN < 11.90 | MAX_TOT_LEN < 11.90))) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT juveniles
  ungroup()
hs_juvenile_abs$LIFE_STAGE = rep("JUVENILE", nrow(hs_juvenile_abs))
hs_juvenile_abs$PRES = ifelse(hs_juvenile_abs$N > 0, 1, 0) 
hs_juvenile_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(hs_juvenile_abs))

# combining the juvenile stage-specific data and inferred absence data, again 
# using distinct to remove repeating values
hs_juvenile_mg = do.call("rbind", list(hs_juvenile, hs_juvenile_abs, hs_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# compiling all MG datasets just in case
hs_mg = do.call("rbind", list(hs_adult_mg, hs_subadult_mg, hs_juvenile_mg))

# and getting rid of all intermediate tables
rm(list = c("hs_adult", "hs_adult_abs", "hs_subadult", "hs_subadult_abs", 
            "hs_juvenile", "hs_juvenile_abs", "hs", "hs_abs", "hs_pres_sites"))

# now combining the two data sources - RVC and MG - and clipping to domain
# adults
hs_adult_clip = full_join(hs_adult_rvc, hs_adult_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) # re-project 
hs_adult_clip = hs_adult_clip[domain, ]
compareCRS(hs_adult_clip, my_crs) # just to be sure

# sub-adults
hs_subadult_clip = full_join(hs_subadult_rvc, hs_subadult_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) # re-project 
hs_subadult_clip = hs_subadult_clip[domain, ]
compareCRS(hs_subadult_clip, my_crs) # just to be sure

# juvenile
hs_juvenile_clip = full_join(hs_juvenile_rvc, hs_juvenile_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) # re-project 
hs_juvenile_clip = hs_juvenile_clip[domain, ]
compareCRS(hs_juvenile_clip, my_crs) # just to be sure

# save presence-absence data for logistic regression modeling
st_write(hs_adult_clip %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(fish_wd, "Presence_Absence/Adult_Bluestriped_Grunt_PA_Full.csv"), 
         append = FALSE)

st_write(hs_subadult_clip %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(fish_wd, "Presence_Absence/Subadult_Bluestriped_Grunt_PA_Full.csv"), 
         append = FALSE)

st_write(hs_juvenile_clip %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(fish_wd, "Presence_Absence/Juvenile_Bluestriped_Grunt_PA_Full.csv"), 
         append = FALSE)

# now filtering out only the presence records for MaxEnt modeling
# MaxEnt only needs species name, lon, lat (in that order)
st_write(hs_adult_clip %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), 
         dsn = paste0(fish_wd, "Presence_Only/Adult_Bluestriped_Grunt_PO_Full.csv"), 
         append = FALSE)

st_write(hs_subadult_clip %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M),
         dsn = paste0(fish_wd, "Presence_Only/Subadult_Bluestriped_Grunt_PO_Full.csv"), 
         append = FALSE)

st_write(hs_juvenile_clip %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), 
         dsn = paste0(fish_wd, "Presence_Only/Juvenile_Bluestriped_Grunt_PO_Full.csv"), 
         append = FALSE)

# create sampling effort raster to parse out sampling bias: rule of thumb for 
# selecting bandwidth according to Scott (1992) and Bowman and Azzalini (1997)
choose_bw = function(spdf) {
  X = coordinates(spdf)
  sigma = c(sd(X[,1]), sd(X[,2])) * (2 / (3 * nrow(X))) ^ (1/6)
}

sampling_effort = full_join(rvc_sites, mg_sites) %>%
  select(LON_DEGREES, LAT_DEGREES) %>%
  st_as_sf(., coords = c(1, 2), crs = gcs) %>% 
  st_transform(., my_crs) %>% 
  add_column("count" = 1) 
sampling_effort = sampling_effort[domain, ] %>%
  as(., "Spatial")

# create a template grid using the habitat raster as a guide
domain_grid = raster(ncol = ncol(habitat), nrow = nrow(habitat), 
                     xmn = xmin(habitat), xmx = xmax(habitat), 
                     ymn = ymin(habitat), ymx = ymax(habitat))
domain_grid[] <- rep(1,ncell(domain_grid))
domain_bw = choose_bw(sampling_effort)

# calculate kernel density surface
domain_kde = sp.kde(x = sampling_effort,
                    bw = domain_bw,
                    newdata = domain_grid,
                    standardize = T)
tm_shape(domain_kde) + tm_raster() 
domain_kde = raster::mask(raster::crop(domain_kde, habitat), habitat)
compareRaster(domain_kde, habitat, extent = T, crs = T, rowcol = T)
# add small constant value because bias grid can't have 0 in MaxEnt
domain_kde = domain_kde + (1/10000) 
tm_shape(domain_kde) + tm_raster(n = 5, palette = "-RdBu") 
domain_kde = raster::setMinMax(domain_kde) 
writeRaster(domain_kde, filename = paste0(fish_wd, "Presence_Only/Bias.asc"),
            format = "ascii", overwrite = T)

# free up some space in the global environment
rm(list = c("lg_adult_clip", "lg_subadult_clip", "lg_juvenile_clip",
            "hs_adult_clip", "hs_subadult_clip", "hs_juvenile_clip",
            "rvc", "mg", "rvc_sites", "mg_sites", "lg_mg", "lg_rvc",
            "hs_mg", "hs_rvc", "mg_dist", "domain_kde", "fknms", "gcs"))


#### WATER QUALITY DATA ####
libraries("tidyr", "rgdal", "gdalUtils", "raster", "sp", "sf", "tmap", "lwgeom", "rgeos",
          "cleangeo", "tidyverse", "stars", "fasterize", "PNWColors", "spex", "igraph", 
          "spatialEco", "tibble", "ncf", "spdep", "gstat", "geoR", "readxl", "dplyr",
          "parallel", "doParallel", "conflicted") 
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# from The SERC Water Quality Monitoring Network (WQMN) 
# (http://serc.fiu.edu/wqmnetwork/)
wqmn = read_excel(paste0(source_wd, "Water_Conditions/WQFloridaKeys&Shelf (ppm) UPDATED 6-6-2020.xlsx"), 
                  sheet = "Data in ppm") 
head(wqmn, 5)

# pre-processing and cleaning 
wqmn = wqmn %>%
  mutate(DATE = as.Date(wqmn$DATE, origin = "1899-12-30")) %>% # excel origin http://support.microsoft.com/kb/214330
  separate(DATE, into = c("YEAR", "MONTH", "DAY"), sep = "-", remove = F) %>% # separate out date data
  filter(YEAR %in% c("2014", "2015", "2016", "2017", "2018")) %>% # keeping only years of interest
  rename(TEMP = `TEMP-B`, SAL = `SAL-B`, DO = `DO-B`) %>% # rename variables because "-" causes problems
  st_as_sf(., coords = c("LONDEC", "LATDEC"), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) %>% # re-project to standard CRS
  mutate(LON_M = sf::st_coordinates(.)[,1], # save LON_M and LAT_M columns
         LAT_M = sf::st_coordinates(.)[,2]) %>%
  st_intersection(., st_make_valid(domain)) %>% # clip data to study domain 
  add_column(SEASON = NA) %>% # four seasons based on quarterly sampling 
  mutate(SEASON = cut(as.numeric(.$MONTH), breaks = c(0,3,6,9,Inf), 
                      labels = c("winter", "spring", "summer", "fall"))) %>%
  select(BASIN, SITE, YEAR, MONTH, DAY, SEASON, LON_M, LAT_M, TEMP, SAL, DO) 

# Biscayne Bay Water Quality (BBWQ) from Miami-Dade County Surface and Groundwater 
# Quality Viewer (https://mdc.maps.arcgis.com/apps/webappviewer/index.html?id=3fd24515ee614f5db63924d7323a4ea7)
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
  add_column(SEASON = NA) %>% # four seasons based on WQMN quarterly sampling 
  mutate(SEASON = cut(as.numeric(.$MONTH), breaks = c(0,3,6,9,Inf), 
                      labels = c("winter", "spring", "summer", "fall"))) %>%
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
  mutate(SEASON = cut(as.numeric(.$MONTH), breaks = c(0,3,6,9,Inf), 
                      labels = c("winter", "spring", "summer", "fall"))) %>%
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
  mutate(SEASON = cut(as.numeric(.$MONTH), breaks = c(0,3,6,9,Inf), 
                      labels = c("winter", "spring", "summer", "fall"))) %>%
  select(BASIN, SITE, YEAR, MONTH, DAY, SEASON, X_COORD, Y_COORD, DO)

#### WINTER DATA ####
# here are the WQMN stations that were sampled every winter for the full 5yr period
win_wqmn_sites = wqmn %>%
  select(SITE, YEAR, MONTH, DAY, SEASON, LON_M, LAT_M, TEMP, SAL, DO) %>%
  filter(SEASON == "winter", !is.na(TEMP), !is.na(SAL), !is.na(DO)) %>%
  group_by(SITE, LON_M, LAT_M) %>%
  count() %>%
  filter(n == 5) %>% # five years (2014, 2015, 2016, 2017, 2018)
  ungroup()

# restricting the WQMN data to only those sites in the win_wqmn_sites data frame
# and calculating summary stats
win_wqmn = wqmn %>%
  filter(SITE %in% win_wqmn_sites$SITE, SEASON == "winter") %>%
  group_by(SITE, LON_M, LAT_M) %>%
  mutate(MEAN_WIN_TEMP = mean(TEMP), SD_WIN_TEMP = sd(TEMP), # mean temp 
         MEAN_WIN_SAL = mean(SAL), SD_WIN_SAL = sd(SAL), # mean sal 
         MEAN_WIN_DO = mean(DO), SD_WIN_DO = sd(DO)) %>% # mean DO 
  select(everything()) %>%
  distinct(BASIN, SITE, LON_M, LAT_M, MEAN_WIN_TEMP, SD_WIN_TEMP,
           MEAN_WIN_SAL, SD_WIN_SAL, MEAN_WIN_DO, SD_WIN_DO)


# annual mean temps, sal, do for Biscayne Bay
win_bb_temp = bb_temp %>%
  filter(SEASON == "winter") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_WIN_TEMP = mean(TEMP), SD_ANN_WIN_TEMP = sd(TEMP)) %>%
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, ANN_MEAN_WIN_TEMP, SD_ANN_WIN_TEMP)

win_bb_sal = bb_sal %>%
  filter(SEASON == "winter") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_WIN_SAL = mean(SAL), SD_ANN_WIN_SAL = sd(SAL)) %>% 
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, ANN_MEAN_WIN_SAL, SD_ANN_WIN_SAL)

win_bb_do = bb_do %>%
  filter(SEASON == "winter") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_WIN_DO = mean(DO), SD_ANN_WIN_DO = sd(DO)) %>% 
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

# restricting the BBWQ data to only those sites that were sampled in each of the five years
# calculating summary stats from the annual means 
win_bb_temp = win_bb_temp %>%
  filter(SITE %in% win_bb_temp_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_WIN_TEMP = mean(ANN_MEAN_WIN_TEMP), SD_WIN_TEMP = sd(ANN_MEAN_WIN_TEMP)) %>% 
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_WIN_TEMP, SD_WIN_TEMP)

win_bb_sal = win_bb_sal %>%
  filter(SITE %in% win_bb_sal_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_WIN_SAL = mean(ANN_MEAN_WIN_SAL), SD_WIN_SAL = sd(ANN_MEAN_WIN_SAL)) %>% 
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_WIN_SAL, SD_WIN_SAL)

win_bb_do = win_bb_do %>%
  filter(SITE %in% win_bb_do_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_WIN_DO = mean(ANN_MEAN_WIN_DO), SD_WIN_DO = sd(ANN_MEAN_WIN_DO)) %>% 
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_WIN_DO, SD_WIN_DO)

# combine the temperature, salinity, and DO data from Biscayne Bay
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

# remove intermediate Biscayne Bay data tables
rm(list = c("win_bb_temp", "win_bb_temp_sites" ,"win_bb_sal", "win_bb_sal_sites",
            "win_bb_do", "win_bb_do_sites"))

# combining the winter water quality data from Biscayne Bay and the FKNMS
winter_wq = rbind(win_wqmn, win_bbwq)

# for variogram modeling and spatial estimation we need SPDF, so convert from sf 
winter_wq_sp = winter_wq %>% st_drop_geometry()
coordinates(winter_wq_sp) = ~ LON_M + LAT_M
proj4string(winter_wq_sp) = proj4string(my_crs)
summary(winter_wq_sp)
write.csv(winter_wq_sp, paste0(csv_wd, "Winter_Water_Conditions.csv"), row.names = F) 


# first taking a glimpse at correlograms and the Moran's I values for the data 
# to ensure that there is spatial dependence
winter_coords = cbind(winter_wq_sp$LON_M, winter_wq_sp$LAT_M) # distance matrix
colnames(winter_coords) = c("LON_M", "LAT_M")
winter_distmat = as.matrix(dist(winter_coords))
winter_maxdist = 2/3 * max(winter_distmat) # maximum distance to consider 

# spline correlograms with 95% pointwise bootstrap CIs
w_temp_corr = spline.correlog(x = winter_wq_sp$LON_M, y = winter_wq_sp$LAT_M, z = winter_wq_sp$MEAN_WIN_TEMP, xmax = winter_maxdist, resamp = 100, type = "boot")
plot(w_temp_corr) 
w_sal_corr = spline.correlog(x = winter_wq_sp$LON_M, y = winter_wq_sp$LAT_M, z = winter_wq_sp$MEAN_WIN_SAL, xmax = winter_maxdist, resamp = 100, type = "boot")
plot(w_sal_corr) 
w_do_corr = spline.correlog(x = winter_wq_sp$LON_M, y = winter_wq_sp$LAT_M, z = winter_wq_sp$MEAN_WIN_DO, xmax = winter_maxdist, resamp = 100, type = "boot")
plot(w_do_corr)

# neighborhood list (neighbors within 15 km distance so every site has at least one neighbor)
winter_neigh = dnearneigh(x = winter_coords, d1 = 0, d2 = 15000, longlat = F)
plot(winter_neigh, coordinates(winter_coords))
winter_wts = nb2listw(neighbours = winter_neigh, style = "W", zero.policy = T) # weights matrix for calculating Moran's I

# Moran's I with normal approximations
moran.test(winter_wq_sp$MEAN_WIN_TEMP, listw = winter_wts, randomisation = F, zero.policy = T) # est. Moran's I stat = 0.65, p < 0.01 sig. spatial dependence
moran.test(winter_wq_sp$MEAN_WIN_SAL, listw = winter_wts, randomisation = F, zero.policy = T) # est. Moran's I stat = 0.45, p < 0.01 sig. spatial dependence
moran.test(winter_wq_sp$MEAN_WIN_DO, listw = winter_wts, randomisation = F, zero.policy = T) # est. Moran's I stat = 0.38, p < 0.01 sig. spatial dependence

# Moran's I with Monte Carlo permutations
moran.mc(winter_wq_sp$MEAN_WIN_TEMP, listw = winter_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.65, p = 0.01 sig. spatial dependence
moran.mc(winter_wq_sp$MEAN_WIN_SAL, listw = winter_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.45, p = 0.01 sig. spatial dependence
moran.mc(winter_wq_sp$MEAN_WIN_DO, listw = winter_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.38, p = 0.01 sig. spatial dependence

# according to Moran tests with both normal and Monte Carlo approximations, there 
# is significant spatial dependence in the winter temp, sal, and DO data.
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
w_sal_fvgm = fit.variogram(w_sal_evgm, vgm(psill = 10, model = "Sph", range = 20000, nugget = 3))
w_sal_svgm_plot = plot(w_sal_evgm, model = w_sal_fvgm, xlab = "distance (m)", pch = 19)
w_sal_svgm_plot
print(w_sal_fvgm)

w_do_evgm = variogram(MEAN_WIN_DO ~ 1, winter_wq_sp, cutoff = winter_maxdist)
plot(w_do_evgm, xlab = "distance (m)", pch = 19)
w_do_fvgm = fit.variogram(w_do_evgm, vgm(psill = 0.08, model = "Sph", range = 10000, nugget = 0.02))
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

# restricting the WQMN data frame to only those sites in the sum_wqmn_sites data 
# frame and calculating summary stats
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
  mutate(ANN_MEAN_SUM_TEMP = mean(TEMP), SD_ANN_SUM_TEMP = sd(TEMP)) %>% 
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, ANN_MEAN_SUM_TEMP, SD_ANN_SUM_TEMP)

sum_bb_sal = bb_sal %>%
  filter(SEASON == "summer") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_SUM_SAL = mean(SAL), SD_ANN_SUM_SAL = sd(SAL)) %>% 
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, ANN_MEAN_SUM_SAL, SD_ANN_SUM_SAL)

sum_bb_do = bb_do %>%
  filter(SEASON == "summer") %>%
  group_by(SITE, X_COORD, Y_COORD, YEAR) %>% # annual means (a mean for each year)
  mutate(ANN_MEAN_SUM_DO = mean(DO), SD_ANN_SUM_DO = sd(DO)) %>% 
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

# restricting the BBWQ data to only those sites that were sampled in each of the
# five years and calculating summary stats from the annual means 
sum_bb_temp = sum_bb_temp %>%
  filter(SITE %in% sum_bb_temp_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_SUM_TEMP = mean(ANN_MEAN_SUM_TEMP), SD_SUM_TEMP = sd(ANN_MEAN_SUM_TEMP)) %>% 
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_SUM_TEMP, SD_SUM_TEMP)

sum_bb_sal = sum_bb_sal %>%
  filter(SITE %in% sum_bb_sal_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_SUM_SAL = mean(ANN_MEAN_SUM_SAL), SD_SUM_SAL = sd(ANN_MEAN_SUM_SAL)) %>% 
  select(everything()) %>%
  distinct(BASIN, SITE, X_COORD, Y_COORD, MEAN_SUM_SAL, SD_SUM_SAL)

sum_bb_do = sum_bb_do %>%
  filter(SITE %in% sum_bb_do_sites$SITE) %>%
  group_by(SITE, X_COORD, Y_COORD) %>% # overall mean (mean of the five annual means)
  mutate(MEAN_SUM_DO = mean(ANN_MEAN_SUM_DO), SD_SUM_DO = sd(ANN_MEAN_SUM_DO)) %>% 
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

# remove intermediate Biscayne Bay data tables
rm(list = c("sum_bb_temp", "sum_bb_temp_sites", "sum_bb_sal", "sum_bb_sal_sites",
            "sum_bb_do", "sum_bb_do_sites"))

# combining the summer water quality data from Biscayne Bay and the FKNMS
summer_wq = rbind(sum_wqmn, sum_bbwq)

# for variogram modeling and spatial estimation we need SPDF, so convert from sf 
summer_wq_sp = summer_wq %>% st_drop_geometry()
coordinates(summer_wq_sp) = ~ LON_M + LAT_M
proj4string(summer_wq_sp) = proj4string(my_crs) 
summary(summer_wq_sp)
write.csv(summer_wq_sp, paste0(csv_wd, "Summer_Water_Conditions.csv"), row.names = FALSE)


# first taking a glimpse at correlograms and the Moran's I values for the data 
# to ensure that there is spatial dependence
summer_coords = cbind(summer_wq_sp$LON_M, summer_wq_sp$LAT_M) # distance matrix
colnames(summer_coords) = c("LON_M", "LAT_M")
summer_distmat = as.matrix(dist(summer_coords))
summer_maxdist = 2/3 * max(summer_distmat) # maximum distance to consider

# spline correlograms with 95% pointwise bootstrap CIs
s_temp_corr = spline.correlog(x = summer_wq_sp$LON_M, y = summer_wq_sp$LAT_M, z = summer_wq_sp$MEAN_SUM_TEMP, xmax = summer_maxdist, resamp = 100, type = "boot")
plot(s_temp_corr) 
s_sal_corr = spline.correlog(x = summer_wq_sp$LON_M, y = summer_wq_sp$LAT_M, z = summer_wq_sp$MEAN_SUM_SAL, xmax = summer_maxdist, resamp = 100, type = "boot")
plot(s_sal_corr) 
s_do_corr = spline.correlog(x = summer_wq_sp$LON_M, y = summer_wq_sp$LAT_M, z = summer_wq_sp$MEAN_SUM_DO, xmax = summer_maxdist, resamp = 100, type = "boot")
plot(s_do_corr)

# neighborhood list (neighbors within 15 km distance so each site has at least one neighbor)
summer_neigh = dnearneigh(x = summer_coords, d1 = 0, d2 = 15000, longlat = F)
plot(summer_neigh, coordinates(summer_coords))
summer_wts = nb2listw(neighbours = summer_neigh, style = "W", zero.policy = T) # weights matrix for calculating Moran's I


# Moran's I with normal approximations
moran.test(summer_wq_sp$MEAN_SUM_TEMP, listw = summer_wts, randomisation = F, zero.policy = T)  # est. Moran's I stat = 0.18, p < 0.01 sig. spatial dependence
moran.test(summer_wq_sp$MEAN_SUM_SAL, listw = summer_wts, randomisation = F, zero.policy = T) # est. Moran's I stat =  0.27, p < 0.01 sig. spatial dependence
moran.test(summer_wq_sp$MEAN_SUM_DO, listw = summer_wts, randomisation = F, zero.policy = T) # est. Moran's I stat = 0.37, p < 0.01 sig. spatial dependence


# Moran's I with Monte Carlo permutations
moran.mc(summer_wq_sp$MEAN_SUM_TEMP, listw = summer_wts, nsim = 99, zero.policy = T)  # est. Moran's I stat = 0.18, p = 0.01 sig. spatial dependence
moran.mc(summer_wq_sp$MEAN_SUM_SAL, listw = summer_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.27, p = 0.01 sig. spatial dependence
moran.mc(summer_wq_sp$MEAN_SUM_DO, listw = summer_wts, nsim = 99, zero.policy = T) # est. Moran's I stat = 0.37, p = 0.01 sig. spatial dependence


# according to Moran tests with both normal and Monte Carlo estimations, 
# there is significant spatial dependence in the summer temp, sal, and DO data.
# now move on to variogram modeling to capture the spatial structure...
require(gstat)

s_temp_evgm = variogram(MEAN_SUM_TEMP ~ 1, summer_wq_sp, cutoff = summer_maxdist) # empirical variogram
plot(s_temp_evgm, xlab = "distance (m)", pch = 19)
s_temp_fvgm = fit.variogram(s_temp_evgm, vgm(psill = 0.15, model = "Sph", range = 30000, nugget = 0.2))
s_temp_svgm_plot = plot(s_temp_evgm, model = s_temp_fvgm, xlab = "distance (m)", pch = 19)
s_temp_svgm_plot
print(s_temp_fvgm)

s_sal_evgm = variogram(MEAN_SUM_SAL ~ 1, summer_wq_sp, cutoff = summer_maxdist)
plot(s_sal_evgm, xlab = "distance (m)", pch = 19)
# this is an odd one, but there is clear spatial dependence according to the Moran's 
# I tests and correlograms, let fit.variogram find the appropriate fitted parameters
s_sal_fvgm = fit.variogram(s_sal_evgm, vgm(model = "Sph"), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 1)
s_sal_fvgm_plot = plot(s_sal_evgm, model = s_sal_fvgm, xlab = "distance (m)", pch = 19)
s_sal_fvgm_plot
print(s_sal_fvgm)

s_do_evgm = variogram(MEAN_SUM_DO ~ 1, summer_wq_sp, cutoff = summer_maxdist)
plot(s_do_evgm, xlab = "distance (m)", pch = 19)
s_do_fvgm = fit.variogram(s_do_evgm, vgm(psill = 0.05, model = "Sph", range = 20000, nugget = 0.01))
s_do_svgm_plot = plot(s_do_evgm, model = s_do_fvgm, xlab = "distance (m)", pch = 19)
s_do_svgm_plot
print(s_do_fvgm)

#### KRIGING ####
# create a prediction grid
grid = (raster(paste0(spatial_wd, "Habitat.asc"))*0)
crs(grid) = my_crs
plot(grid)
grid = grid %>% as(., "SpatialPixels")
proj4string(grid) = proj4string(my_crs)

# initiate cluster and divvy up prediction grid for cores
library(parallel)
no_cores = 10 # this depends on your computer system 
cl = makeCluster(no_cores)
parts = split(x = 1:length(grid), f = 1:no_cores)
stopCluster(cl)

# winter temperature
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "grid", "parts", "w_temp_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_temp_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_TEMP ~ 1, locations = winter_wq_sp,
        newdata = grid[parts[[x]],], model = w_temp_fvgm))
stopCluster(cl)
showConnections()

w_temp_merge = maptools::spRbind(w_temp_par[[1]], w_temp_par[[2]])
for (j in 3:length(w_temp_par)) {
  w_temp_merge = maptools::spRbind(w_temp_merge, w_temp_par[[j]])
}
w_temp_merge = SpatialPixelsDataFrame(points = w_temp_merge, data = w_temp_merge@data)
# save new surface
summary(w_temp_merge)
writeGDAL(w_temp_merge["var1.pred"], fname = paste0(temp_wd, "Mean_Win_Temp.tif"),
          drivername = "GTiff", type = "Float32")
writeGDAL(w_temp_merge["var1.var"], fname = paste0(temp_wd, "Mean_Win_Temp_Var.tif"),
          drivername = "GTiff", type = "Float32")

w_temp = raster(paste0(temp_wd, "Mean_Win_Temp.tif"))
crs(w_temp) = my_crs
extent = extent(raster(paste0(spatial_wd, "Habitat.asc")))
w_temp2 = extend(w_temp, extent)
w_temp2
writeRaster(w_temp2, filename = paste0(spatial_wd, "Mean_win_Temp.asc"), 
            format = "ascii", overwrite = T)

# winter salinity
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "grid", "parts", "w_sal_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_sal_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_SAL ~ 1, locations = winter_wq_sp,
        newdata = grid[parts[[x]],], model = w_sal_fvgm))
stopCluster(cl)
showConnections()

w_sal_merge = maptools::spRbind(w_sal_par[[1]], w_sal_par[[2]])
for (j in 3:length(w_sal_par)) {
  w_sal_merge = maptools::spRbind(w_sal_merge, w_sal_par[[j]])
}
w_sal_merge = SpatialPixelsDataFrame(points = w_sal_merge, data = w_sal_merge@data)
# save new surface
summary(w_sal_merge)
writeGDAL(w_sal_merge["var1.pred"], fname = paste0(temp_wd, "Mean_Win_Sal.tif"),
          drivername = "GTiff", type = "Float32")
writeGDAL(w_sal_merge["var1.var"], fname = paste0(temp_wd, "Mean_Win_Sal_Var.tif"),
          drivername = "GTiff", type = "Float32")

w_sal = raster(paste0(temp_wd, "Mean_Win_Sal.tif"))
crs(w_sal) = my_crs
extent = extent(raster(paste0(spatial_wd, "Habitat.asc")))
w_sal2 = extend(w_sal, extent)
w_sal2
writeRaster(w_sal2, filename = paste0(spatial_wd, "Mean_Win_Sal.asc"), 
            format = "ascii", overwrite = T)

# winter dissolved oxygen
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("winter_wq_sp", "grid", "parts", "w_do_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
w_do_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_WIN_DO ~ 1, locations = winter_wq_sp,
        newdata = grid[parts[[x]],], model = w_do_fvgm))
stopCluster(cl)
showConnections()

w_do_merge = maptools::spRbind(w_do_par[[1]], w_do_par[[2]])
for (j in 3:length(w_do_par)) {
  w_do_merge = maptools::spRbind(w_do_merge, w_do_par[[j]])
}
w_do_merge = SpatialPixelsDataFrame(points = w_do_merge, data = w_do_merge@data)
# save new surface
summary(w_do_merge)
writeGDAL(w_do_merge["var1.pred"], fname = paste0(temp_wd, "Mean_Win_DO.tif"),
          drivername = "GTiff", type = "Float32")
writeGDAL(w_do_merge["var1.var"], fname = paste0(temp_wd, "Mean_Win_DO_Var.tif"),
          drivername = "GTiff", type = "Float32")

w_do = raster(paste0(temp_wd, "Mean_Win_DO.tif"))
crs(w_do) = my_crs
extent = extent(raster(paste0(spatial_wd, "Habitat.asc")))
w_do2 = extend(w_do, extent)
w_do2
writeRaster(w_do2, filename = paste0(spatial_wd, "Mean_Win_DO.asc"), 
            format = "ascii", overwrite = T)

# summer temperature
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "grid", "parts", "s_temp_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_temp_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_TEMP ~ 1, locations = summer_wq_sp,
        newdata = grid[parts[[x]],], model = s_temp_fvgm))
stopCluster(cl)
showConnections()

s_temp_merge = maptools::spRbind(s_temp_par[[1]], s_temp_par[[2]])
for (j in 3:length(s_temp_par)) {
  s_temp_merge = maptools::spRbind(s_temp_merge, s_temp_par[[j]])
}
s_temp_merge = SpatialPixelsDataFrame(points = s_temp_merge, data = s_temp_merge@data)
# save new surface
summary(s_temp_merge)
writeGDAL(s_temp_merge["var1.pred"], fname = paste0(temp_wd, "Mean_Sum_Temp.tif"),
          drivername = "GTiff", type = "Float32")
writeGDAL(s_temp_merge["var1.var"], fname = paste0(temp_wd, "Mean_Sum_Temp_Var.tif"),
          drivername = "GTiff", type = "Float32")

s_temp = raster(paste0(temp_wd, "Mean_Sum_Temp.tif"))
crs(s_temp) = my_crs
extent = extent(raster(paste0(spatial_wd, "Habitat.asc")))
s_temp2 = extend(s_temp, extent)
s_temp2
writeRaster(s_temp2, filename = paste0(spatial_wd, "Mean_Sum_Temp.asc"), 
            format = "ascii", overwrite = T)

# summer salinity
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "grid", "parts", "s_sal_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_sal_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_SAL ~ 1, locations = summer_wq_sp,
        newdata = grid[parts[[x]],], model = s_sal_fvgm))
stopCluster(cl)
showConnections()

s_sal_merge = maptools::spRbind(s_sal_par[[1]], s_sal_par[[2]])
for (j in 3:length(s_sal_par)) {
  s_sal_merge = maptools::spRbind(s_sal_merge, s_sal_par[[j]])
}
s_sal_merge = SpatialPixelsDataFrame(points = s_sal_merge, data = s_sal_merge@data)
# save new surface
summary(s_sal_merge)
writeGDAL(s_sal_merge["var1.pred"], fname = paste0(temp_wd, "Mean_Sum_Sal.tif"),
          drivername = "GTiff", type = "Float32")
writeGDAL(s_sal_merge["var1.var"], fname = paste0(temp_wd, "Mean_Sum_Sal_Var.tif"),
          drivername = "GTiff", type = "Float32")

s_sal = raster(paste0(temp_wd, "Mean_Sum_Sal.tif"))
crs(s_sal) = my_crs
extent = extent(raster(paste0(spatial_wd, "Habitat.asc")))
s_sal2 = extend(s_sal, extent)
s_sal2
writeRaster(s_sal2, filename = paste0(spatial_wd, "Mean_Sum_Sal.asc"), 
            format = "ascii", overwrite = T)

# summer dissolved oxygen
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = c("summer_wq_sp", "grid", "parts", "s_do_fvgm"),
              envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
s_do_par = parLapply(cl = cl, X = 1:no_cores, fun = function(x)
  krige(formula = MEAN_SUM_DO ~ 1, locations = summer_wq_sp,
        newdata = grid[parts[[x]],], model = s_do_fvgm))
stopCluster(cl)
showConnections()

s_do_merge = maptools::spRbind(s_do_par[[1]], s_do_par[[2]])
for (j in 3:length(s_do_par)) {
  s_do_merge = maptools::spRbind(s_do_merge, s_do_par[[j]])
}
s_do_merge = SpatialPixelsDataFrame(points = s_do_merge, data = s_do_merge@data)
# save new surface
summary(s_do_merge)
writeGDAL(s_do_merge["var1.pred"], fname = paste0(temp_wd, "Mean_Sum_DO.tif"),
          drivername = "GTiff", type = "Float32")
writeGDAL(s_do_merge["var1.var"], fname = paste0(temp_wd, "Mean_Sum_DO_Var.tif"),
          drivername = "GTiff", type = "Float32")

s_do = raster(paste0(temp_wd, "Mean_Sum_DO.tif"))
crs(s_do) = my_crs
extent = extent(raster(paste0(spatial_wd, "Habitat.asc")))
s_do2 = extend(s_do, extent)
s_do2
writeRaster(s_do2, filename = paste0(spatial_wd, "Mean_Sum_DO.asc"), 
            format = "ascii", overwrite = T)

#### SEAFLOOR MORPHOLOGY ####
# Though the seafloor surface morphology rasters exported from ArcGIS were
# calculated using the depth raster, ArcGIS and RStudio export ASCII rasters with
# a slightly different level of precision, causing their extents to vary a bit. 
# This is a problem, as MaxEnt requires all ASCII files to have the SAME EXACT
# origin and extent. To resolve this, read in the tiff files exported from ArcGIS
# and export them to the spatial predictors folder as ASCIIs.  

# BPI Broad
bpi_broad = raster(paste0(temp_wd, "BPI_Broad.tif"))
crs(bpi_broad) = my_crs
writeRaster(bpi_broad, paste0(spatial_wd, "BPI_Broad.asc"),
            format = "ascii", overwrite = TRUE)

# BPI Fine
bpi_fine = raster(paste0(temp_wd, "BPI_Fine.tif"))
crs(bpi_fine) = my_crs
writeRaster(bpi_fine, paste0(spatial_wd, "BPI_Fine.asc"),
            format = "ascii", overwrite = TRUE)

# Curvature
curvature = raster(paste0(temp_wd, "Curvature.tif"))
crs(curvature) = my_crs
writeRaster(curvature, paste0(spatial_wd, "Curvature.asc"),
            format = "ascii", overwrite = TRUE)

# Plan Curvature
plan_curve = raster(paste0(temp_wd, "Plan_Curve.tif"))
crs(plan_curve) = my_crs
writeRaster(plan_curve, paste0(spatial_wd, "Plan_Curve.asc"),
            format = "ascii", overwrite = TRUE)

# Slope
slope = raster(paste0(temp_wd, "Slope.tif"))
crs(slope) = my_crs
writeRaster(slope, paste0(spatial_wd, "Slope.asc"),
            format = "ascii", overwrite = TRUE)

# Rugosity
rugosity = raster(paste0(temp_wd, "Rugosity.tif"))
crs(rugosity) = my_crs
writeRaster(rugosity, paste0(spatial_wd, "Rugosity.asc"),
            format = "ascii", overwrite = TRUE)

# Standard deviation of depth
StDev_Depth = raster(paste0(temp_wd, "Depth_sdev_003.tif"))
crs(rugosity) = my_crs
writeRaster(StDev_Depth, paste0(spatial_wd, "StDev_Depth.asc"),
            format = "ascii", overwrite = TRUE)
