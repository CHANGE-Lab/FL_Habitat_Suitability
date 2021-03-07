#### SET-UP ####
# This script was used to prepare data for Courtney Stuart's first MSc chapter in the lab of 
# Dr. Stephanie Green at the University of Alberta. Data are specific to southern Florida and 
# include: benthic habitat classifications, bathymetric and topographic surfaces, georeferenced 
# reef fish occurrence records, and bottom water conditions. 

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
libraries("rgdal", "gdalUtils", "raster", "sp", "sf", "tmap", "dplyr", "lwgeom", "rgeos",
          "cleangeo", "tidyverse", "stars", "fasterize", "PNWColors", "spex", "igraph", 
          "spatialEco")


# change where large temp rasters are saved
rasterOptions(tmpdir = "Z:/Courtney/Stuart_MSc_Ch1/Temp/")

# save PROJ.4 string for standard projection (ESPG:26958 NAD 83/Florida East) and source LAT/LON data (EPSG:4326 WGS 84/World Geodetic System 1984)
my_crs = CRS("+init=epsg:26958")
gcs = CRS("+init=epsg:4326")

#### HABITAT & DEPTH DATA ####

# read in park/marine sanctuary polygons
# Florida Keys National Marine Sanctuary shapefile (https://sanctuaries.noaa.gov/library/imast_gis.html)

fknms = st_read(dsn = paste0(source_wd, "Parks/fknms_py.shp")) %>%
  st_transform(., my_crs)
compareCRS(fknms, my_crs) # check projection
st_is_valid(fknms, reason = T) # check geometry

# National Park Service shapefile (https://public-nps.opendata.arcgis.com/datasets/nps-boundary-1/data)
nps = st_read(dsn = paste0(source_wd, "Parks/NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service.shp")) %>%
  st_transform(., my_crs)
bnp = nps[nps$UNIT_NAME == "Biscayne National Park",] # extract Biscayne Bay
compareCRS(bnp, my_crs)
st_is_valid(bnp)

# how do the shapefiles line up
tm_shape(bnp) +
  tm_fill("navy") +
  tm_shape(fknms) +
  tm_fill("gray")

# there are gaps at the seams where FKNMS and BNP should meet, here is an arbitrary
# polygon that I constructed in a GIS to fill the interior gaps while still respecting
# the outer boundaries of the two parks
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

# have an outline now of the parks' outer borders, so remove temporary data
rm(list = c("fknms", "bnp", "nps", "fill_gaps"))

# make simple polygon feature for clipping out the training area 
# (constrained by Tavernier Creek to the southwest & Key Biscayne to the northeast)
train_poly = data.frame(
  lon = c(236904.809, 275906.500, 308636.563, 272163.225), # bbox values in meters
  lat = c(78590.323, 159992.442, 142403.545, 62074.390)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = my_crs) %>%
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# make a simple polygon feature for clipping out the testing area 
# (Tavernier Creek to the northeast & Cudjoe Key to the southwest)
test_poly = data.frame(
  lon = c(149311.942, 244481.888, 250880.718, 250698.734, 149245.026), # bbox values in meters
  lat = c(75028.405, 75041.028, 72043.654, 20743.032, 20841.839)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = my_crs) %>%
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# clip marine parks geometry to training and testing areas 
parks_train = st_intersection(train_poly, parks_union)
parks_test = st_intersection(test_poly, parks_union)

tm_shape(parks_train) +
  tm_fill("navy") +
  tm_shape(parks_test) +
  tm_fill("gray")

st_write(parks_train, paste0(temp_wd, "parks_train.shp"))
st_write(parks_test, paste0(temp_wd, "parks_test.shp"))

# load habitat data from Unified Reef Tract Map 
# (https://myfwc.com/research/gis/regional-projects/unified-reef-map/)
require(rgdal)
fgdb = paste0(source_wd, "Unified_Reef_Map/FWC_UnifiedFloridaReefMap_v2.0.gdb")

# List all feature classes in file gdb
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

# create IDs because MaxEnt requires categorical variables to be defined numerically, not with words
ClassLv1 = unique(reef_map$ClassLv1) # list categories
ClassLv1_df = data.frame(ID = 1:length(ClassLv1), ClassLv1 = ClassLv1) # create a data.frame of IDs
reef_map$ClassLv1_ID = ClassLv1_df$ID[match(reef_map$ClassLv1, ClassLv1_df$ClassLv1)] # match IDs
unique(reef_map$ClassLv1_ID) # did it work? yes. 
write.csv(ClassLv1_df, paste0(csv_wd, "URM_ClassLv1_IDs.csv"), row.names = F) # save for metadata

# clip unified reef map data to training and testing areas
reef_train = st_intersection(st_make_valid(reef_map), parks_train)
reef_test = st_intersection(st_make_valid(reef_map), parks_test)

# save for mapping
# ignore warnings about abbreviating field names and having too many decimals in the field "Shp_Ar", that's OK
st_write(reef_train, dsn = paste0(train_wd, "GIS/Reef_Train.shp"), driver = "ESRI Shapefile", append = F)
st_write(reef_test, dsn = paste0(test_wd, "GIS/Reef_Test.shp"), driver = "ESRI Shapefile", append = F)

# create a palette for plotting benthic habitat classes
pal = pnw_palette("Bay", 14, type = "continuous") 

# supplementary shoreline mangrove habitat data
mg_shore = 
  st_read(dsn = paste0(source_wd, "Mangrove_Habitat/Mangrove_Habitat_in_Florida.shp")) %>%
  filter(!st_is_empty(.)) %>%
  st_transform(., my_crs) %>%
  st_cast(., "MULTIPOLYGON")
compareCRS(mg_shore, my_crs)

# there is a small gap between the mangrove data and the reef tract map along the mainland
# coast, and these missing areas are important coastline mangroves (visible in satellite 
# imagery) add a 10 m buffer around the mangrove data and then use gDifference to keep only
# the non-overlapping areas (AKA, respect the boundaries of the reef map).
mg_buff = st_buffer(mg_shore, dist = 10)

# now keep only the non-overlapping regions of the mangrove data (FYI: time-consuming step)
mg_train = rgeos::gDifference(as_Spatial(st_intersection(mg_buff, parks_train)),
                              as_Spatial(reef_train)) %>% st_as_sf()

mg_test = rgeos::gDifference(as_Spatial(st_intersection(mg_buff, parks_test)), 
                             as_Spatial(reef_test)) %>% st_as_sf()

tm_shape(parks_train) + tm_borders("gray") + 
  tm_shape(mg_train) + tm_fill("orange")

# double check the ID assigned to mangroves from the reef map and add it to the mangrove data
ClassLv1_df
mg_train$ClassLv1 = rep(as.character("Mangrove"), nrow(mg_train))
mg_train$ClassLv1_ID = rep(as.integer(11), nrow(mg_train))
mg_test$ClassLv1 = rep(as.character("Mangrove"), nrow(mg_test))
mg_test$ClassLv1_ID = rep(as.integer(11), nrow(mg_test))

# save for mapping
st_write(mg_train, dsn = paste0(train_wd, "GIS/Mangrove_Train.shp"), 
         driver = "ESRI Shapefile", append = F)
st_write(mg_test, dsn = paste0(test_wd, "GIS/Mangrove_Test.shp"), 
         driver = "ESRI Shapefile", append = F)

#  check geometry for both habitat layers
head(mg_train, 1) # sf column: geometry 
head(reef_train, 1) # sf column: Shape

reef_train = reef_train %>% rename(geometry = Shape)
head(reef_train, 1) # fixed

reef_test = reef_test %>% rename(geometry = Shape)
head(reef_test, 1) # fixed

# combining all habitats in the training area
hab_train = rbind(mg_train, select(reef_train, geometry, ClassLv1, ClassLv1_ID)) %>%
  filter(!ClassLv1 %in% c("Land", "Not Classified")) %>% # only want real benthic habitats
  st_cast("MULTIPOLYGON")

tm_shape(hab_train, projection = my_crs) +
  tm_fill("ClassLv1_ID", palette = pal, style = "cat") 

# combining all habitats in the testing area
hab_test = rbind(mg_test, select(reef_test, geometry, ClassLv1, ClassLv1_ID)) %>%
  filter(!ClassLv1 %in% c("Land", "Not Classified")) %>% # only want real benthic habitats
  st_cast("MULTIPOLYGON")

tm_shape(hab_test, projection = my_crs) +
  tm_fill("ClassLv1_ID", palette = pal, style = "cat") 

# template raster for habitat training area w/ 5 x 5 m res and standard CRS
train_ras = raster(ext = extent(hab_train), res = c(5,5), crs = my_crs)

# convert benthic habitat polygons to a temporary raster
hab_train_ras = writeRaster((fasterize(hab_train, train_ras, field = "ClassLv1_ID", fun = "max")), 
                            file = file.path(temp_wd, "hab_train_ras.tif"), format = "GTiff", 
                            overwrite = T)

plot(hab_train_ras, col = pal)

# repeat for testing area
test_ras = raster(ext = extent(hab_test), res = c(5,5))
crs(test_ras) = my_crs
hab_test_ras = writeRaster((fasterize(hab_test, test_ras, field = "ClassLv1_ID", fun = "max")), 
                           file = file.path(temp_wd, "hab_test_ras.tif"), format = "GTiff", 
                           overwrite = T)

plot(hab_test_ras, col = pal)

# remove temporary files 
rm(list = c("reef_map", "reef_train", "reef_test", "fgdb", "mg_buff", 
            "mg_buff_train", "mg_buff_test", "mg_shore", "mg_train", "mg_test")) 

# load Continuously Updated DEM - 1/9 Arc-Second Resolution Bathymetric-Topographic Tiles  from NOAA-NCEI
# (https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ngdc.mgg.dem:999919)

setwd(dem_wd) # temporarily set WD to source DEMs folder (downloaded tiles)
dems = c('Job606638_ncei_nintharcsec_dem_002_001.tif', 
         'Job606638_ncei_nintharcsec_dem_002_000.tif',
         'Job606638_ncei_nintharcsec_dem_001_001.tif',
         'Job606638_ncei_nintharcsec_dem_001_000.tif',
         'Job606638_ncei_nintharcsec_dem_000_000.tif')
         #'Job606638_ncei_nintharcsec_dem_000_001.tif' # this tile is for Gulf of Mexico, don't need it

e = extent(618778.56, 988742.88, 61706.64, 640446.24)

template = raster(e) # use extent to create empty template raster

# set CRS of template using info in metadata of tiles
# Horizontal CS: EPSG:3512 NAD83(NSRS2007) / Florida East (US foot)
# Vertical CS: 6360 NAVD88_height_(ftUS)
template_crs = CRS("+init=epsg:3512")
crs(template) = template_crs

# GDAL version 2.1.3
# write the empty template raster 
writeRaster(template, file = file.path(temp_wd, "dem_mosaic.tif"), format = "GTiff",
            overwrite = T)

# mosaic tiles and save them to the empty template raster
gdalUtils::mosaic_rasters(gdalfile = dems, dst_dataset = file.path(temp_wd, "dem_mosaic.tif"),
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

# the projected DEM still stores depth data in units US feet, convert to meters instead
dem_proj = dem_proj * 0.304801  # (1 ft = 0.304801 m)
dem_proj = raster::setMinMax(dem_proj) # calculate and save the min and max values of the raster
dem_proj # check attributes

# the DEM data will be used in ArcGIS with the Benthic Terrain Modeler
# extension to calculate seafloor morphometrics -- the rugosity and slope
# tools will require the FULL DEM (projected mosaic) rather than the clipped (train
# & test) DEMs because these tools do NOT respect ArcGIS' environment settings
writeRaster(dem_proj, paste0(dem_wd, "Full_Mosaic_DEM.tif"), format = "GTiff",
            overwrite = T)

pal2 = pnw_palette("Winter", n = 100, "continuous") # palette for plotting depth
plot(dem_proj, col = pal2)

# remove temporary products
rm(list = c("dems", "e", "template", "template_crs")) 

# back to normal working directory
setwd("Z:/Courtney/Stuart_MSc_Ch1/")

# clip habitat data to DEM extent first, then clip DEM back to habitat; mask and match extents. 
# Starting with the training area:
hab_train_crop_dem = writeRaster((raster::crop(hab_train_ras, dem_proj)), 
                                 file = file.path(temp_wd, "hab_train_crop_dem.tif"), 
                                 format = "GTiff", overwrite = T)
dem_crop_hab_train = writeRaster((raster::crop(dem_proj, hab_train_ras)), 
                                 file = file.path(temp_wd, "dem_crop_hab_train.tif"), 
                                 format = "GTiff", overwrite = T)
extent(hab_train_crop_dem) = extent(dem_crop_hab_train) # match extents
habitat_train = raster::mask(hab_train_crop_dem, dem_crop_hab_train) # habitat data within DEM extent
depth_train = raster::mask(dem_crop_hab_train, hab_train_crop_dem) # depth data matching habitat extent 
extent(habitat_train) = extent(depth_train)
compareRaster(depth_train, habitat_train, extent = T, res = T, crs = T, rowcol = T)

# they should line up perfectly 
tm_shape(depth_train) +
  tm_raster(palette = pal2, style = "cont") +
  tm_shape(habitat_train) +
  tm_raster(palette = pal, n = 14) +
  tm_layout(frame = F)

# write out data - need to use ArcGIS' Benthic Terrain Modeler extension to 
# calculate seascape topography metrics from the depth raster (*remember, the
# slope and rugosity tools will require the FULL DEM in the source DEM folder)
# MaxEnt requires ascii format
writeRaster(depth_train, file = file.path(train_wd, "Environmental/Depth.asc"), 
            format = "ascii", overwrite = T)
writeRaster(habitat_train, file = file.path(train_wd, "Environmental/Habitat.asc"), 
            format = "ascii", overwrite = T)

# create a simple, constant-value raster from the habitat training raster for clipping purposes
writeRaster(((habitat_train * 0) + 1), file = file.path(temp_wd, "simp_train.tif"), 
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

# now convert simple training raster to simple polygon
polygonize(srcfile = paste0(temp_wd, "simp_train.tif"), 
           dstfile = file.path(temp_wd, "simp_train.shp"), 
           ogr_format ='ESRI Shapefile', connect8 = F)

train = st_read(dsn = paste0(temp_wd, "simp_train.shp")) %>%
  st_transform(., my_crs)
train = train[train$Value == 1,] # selecting only areas that had training habitat data
plot(train)
# saving training area as a shapefile
st_write(train, dsn = paste0(train_wd, "GIS/Training_Area.shp"), append = F)

# repeat steps for testing area
# clip habitat data to DEM extent first, then clip DEM back to habitat; mask and match extents
hab_test_crop_dem = writeRaster((raster::crop(hab_test_ras, dem_proj)), 
                                file = file.path(temp_wd, "hab_test_crop_dem.tif"), 
                                format = "GTiff", overwrite = T)
dem_crop_hab_test = writeRaster((raster::crop(dem_proj, hab_test_ras)), 
                                file = file.path(temp_wd, "dem_crop_hab_test.tif"), 
                                format = "GTiff", overwrite = T)
extent(hab_test_crop_dem) = extent(dem_crop_hab_test) # match extents
crs(hab_test_crop_dem) = my_crs
crs(dem_crop_hab_test) = my_crs

# write out data - need to use ArcGIS' Benthic Terrain Modeler extension to calculate seascape topography 
# metrics from the depth raster; MaxEnt requires ascii format
habitat_test = writeRaster((raster::mask(hab_test_crop_dem, dem_crop_hab_test)),
                           file = file.path(test_wd, "Environmental/Habitat.asc"),
                           format = "ascii", overwrite = T)# habitat data within DEM extent
depth_test = writeRaster((raster::mask(dem_crop_hab_test, hab_test_crop_dem)),
                         file = file.path(test_wd, "Environmental/Depth.asc"),
                         format = "ascii", overwrite = T)# depth data matching habitat extent
extent(habitat_test) = extent(depth_test)
compareRaster(depth_test, habitat_test, extent = T, res = T, crs = T, rowcol = T)
writeRaster(((habitat_test * 0) + 1), file = file.path(temp_wd, "simp_test.tif"), 
            format = "GTiff", overwrite = T)


# should overlap perfectly 
tm_shape(depth_test) +
  tm_raster(palette = pal2, style = "cont") +
  tm_shape(habitat_test) +
  tm_raster(palette = pal) +
  tm_layout(frame = F)

# now convert simple testing raster to simple polygon
polygonize(srcfile = paste0(temp_wd, "simp_test.tif"), 
           dstfile = file.path(temp_wd, "simp_test.shp"), 
           ogr_format='ESRI Shapefile', connect8 = F)

test = st_read(dsn = paste0(temp_wd, "simp_test.shp")) %>%
  st_transform(., my_crs)
test = test[test$Value == 1,] # selecting only areas that had testing habitat data
plot(test)
# saving the testing area as a shapefile
st_write(test, dsn = paste0(test_wd, "GIS/Testing_Area.shp"), append = F)

# clearing out some space in the environment
rm(list = c("hab_test", "hab_test_crop_dem", "hab_test_ras", "hab_train_crop_dem", 
            "hab_train_ras", "hab_train", "dem_crop_hab_test", "dem_crop_hab_train", "dem_proj"))

# lastly, create "proximity to mangrove" surfaces by calculating Euclidean distance from each 
# cell to the nearest mangrove cell
ClassLv1_df # again finding mangrove ID --> 11 is what we want to keep
m = c(1, NA, 2, NA, 3, NA, 4, NA, 5, NA, 6, NA, 7, NA, 8, NA,
      9, NA, 10, NA, 11, 11, 12, NA, 13, NA, 14, NA) # cells with IDs 1:10 & 12-14 become NA
mat = matrix(m, ncol = 2, byrow = T)
mg_euc_train = writeRaster(raster::mask(raster::crop(gridDistance(habitat_train, origin = 11),
                                                     habitat_train), habitat_train),
                           file = file.path(train_wd, "Environmental/Mangrove_Dist.asc"),
                           format = "ascii", overwrite = T)
mg_euc_test = writeRaster(raster::mask(raster::crop(gridDistance(habitat_test, origin = 11),
                                                    habitat_test), habitat_test),
                          file = file.path(test_wd, "Environmental/Mangrove_Dist.asc"),
                          format = "ascii", overwrite = T) 


##### FISH DATA #####
library(rvc)
require(rvc)
rvc = getSampleData(years = c(2014, 2016, 2018), regions = "FLA KEYS")
head(rvc, 1)

# convert the fork length column (LEN) to total length (TOT_LEN) using FishBase 
# length-length conversion for gray snapper specifically TL = 0 + 1.049 x FL
# this is temporary, will repeat later for bluestriped grunt using their specific FL-TL conversion)
rvc = rvc %>% mutate(TOT_LEN = (LEN*1.049)) 

# how many unique sites were sampled in 2014, 2016, and 2018?
rvc_sites = rvc %>% distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, .keep_all = T)

#### gray snapper ####
# start with presence and absence records specific to the adult gray snapper stage (> 24.71 cm TL (size at maturity))
lg_adult = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% # group by both primary and second-stage sample units
  filter(SPECIES_CD == "LUT GRIS" & TOT_LEN > 24.71) %>% 
  summarize(N = sum(NUM)) %>% # was an adult gray snapper ever seen at this second-stage sample unit (SSU)?
  ungroup() %>%
  distinct()
lg_adult$LIFE_STAGE = rep("ADULT", nrow(lg_adult)) # specify that these are all adult records
lg_adult$PRES = ifelse(lg_adult$N > 0, 1, 0) # if at least one adult was seen at the SSU, it's a presence [1], else an absence [0]
lg_adult$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_adult)) # source column for when RVC and MG data are compiled

# now the inferred absences sites (sites where either no gray snapper were seen or only those of another age class were seen)
lg_adult_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "LUT GRIS" & !TOT_LEN > 24.71) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 because these are adult absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
lg_adult_abs$LIFE_STAGE = rep("ADULT", nrow(lg_adult_abs))
lg_adult_abs$PRES = ifelse(lg_adult_abs$N > 0, 1, 0)
lg_adult_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_adult_abs))

# combining adult absence and presence data and using distinct because some of the absence sites might be shared across
# both the adult stage-specific data and the inferred absence data
lg_adult_rvc = rbind(lg_adult, lg_adult_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, .keep_all = T)

# repeat for the subadult gray snapper stage (9.51 cm <= TOT LEN <= 24.71 cm (between size at 1 YR and size at maturation))
lg_subadult = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% 
  filter(SPECIES_CD == "LUT GRIS" & (TOT_LEN >= 9.51 & TOT_LEN <= 24.71)) %>% 
  summarize(N = sum(NUM)) %>% # was a subadult gray snapper ever seen at this second-stage sample unit (SSU)?
  ungroup() %>%
  distinct()
lg_subadult$LIFE_STAGE = rep("SUBADULT", nrow(lg_subadult)) 
lg_subadult$PRES = ifelse(lg_subadult$N > 0, 1, 0)
lg_subadult$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_subadult))

# now the inferred absences sites (sites where either no gray snapper were seen or only those of another age class were seen)
lg_subadult_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "LUT GRIS" & !(TOT_LEN >= 9.51 & TOT_LEN <= 24.71)) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 because these are subadult absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
lg_subadult_abs$LIFE_STAGE = rep("SUBADULT", nrow(lg_subadult_abs))
lg_subadult_abs$PRES = ifelse(lg_subadult_abs$N > 0, 1, 0)
lg_subadult_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_subadult_abs))

# combining subadult absence and presence data and using distinct because some of the absence sites might be shared across
# both the subadult stage-specific data and the inferred absence data
lg_subadult_rvc = rbind(lg_subadult, lg_subadult_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, .keep_all = T)

# finally, repeat for the juvenile gray snapper stage (TOT LEN < 9.51 cm (smaller than size at 1 YR), but not equal to 0!)
lg_juvenile = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% 
  filter(SPECIES_CD == "LUT GRIS" & (TOT_LEN != 0 & TOT_LEN < 9.51)) %>% 
  summarize(N = sum(NUM)) %>% # was a juvenile gray snapper ever seen at this second-stage sample unit (SSU)?
  ungroup() %>%
  distinct()
lg_juvenile$LIFE_STAGE = rep("JUVENILE", nrow(lg_juvenile)) 
lg_juvenile$PRES = ifelse(lg_juvenile$N > 0, 1, 0)
lg_juvenile$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_juvenile))

# now the inferred absences sites (sites where either no gray snapper were seen or only those of another age class were seen)
lg_juvenile_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "LUT GRIS" & !(TOT_LEN != 0 & TOT_LEN < 9.51)) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 because these are juvenile absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
lg_juvenile_abs$LIFE_STAGE = rep("JUVENILE", nrow(lg_juvenile_abs))
lg_juvenile_abs$PRES = ifelse(lg_juvenile_abs$N > 0, 1, 0)
lg_juvenile_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_juvenile_abs))

# combining juvenile absence and presence data and using distinct because some of the absence sites might be shared across
# both the juvenile stage-specific data and the inferred absence data
lg_juvenile_rvc = rbind(lg_juvenile, lg_juvenile_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, .keep_all = T)

# compiling all RVC datasets just in case
lg_rvc = do.call("rbind", list(lg_adult_rvc,
                               lg_subadult_rvc,
                               lg_juvenile_rvc))

# and getting rid of all intermediate tables
rm(list = c("lg_adult", "lg_adult_abs", "lg_subadult", "lg_subadult_abs", "lg_juvenile", "lg_juvenile_abs"))

# Reef Visual Census (RVC) data done, now time for Southeast Fisheries Science Center Mangrove Study Data (MG)
mg = read.csv(paste0(source_wd, "Fish_Data/Shoreline_Mangrove_Surveys/DATA_Serafy_BB_MANGROVE_FISH_DATA_1998W_2019W.csv"), stringsAsFactors = F)
head(mg, 1)

# filter out years of interest (2014, 2016, 2018)
# convert min, avg, and max total length values from inches to centimeters (1 in = 2.54 cm)
# rename lat/long and species name columns to match RVC dataset (to join tables later)
mg = mg %>%
  filter(YR %in% c(2014, 2016, 2018)) %>%
  mutate(MIN_TOT_LEN = (as.numeric(MIN_IN))*2.54) %>%
  mutate(AVE_TOT_LEN = (as.numeric(AVE_IN))*2.54) %>%
  mutate(MAX_TOT_LEN = (as.numeric(MAX_IN))*2.54) %>%
  rename(YEAR = YR, MONTH = MO, DAY = DY, SPECIES_CD = SP, LAT_DEGREES = LAT, LON_DEGREES = LON)

# first find absences shared across life stages by isolating sites where no gray snapper were observed
mg_sites = mg %>% distinct(Site, LON_DEGREES, LAT_DEGREES)

# gray snapper presence sites
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

# need to parse out life history stages using minimum and maximum body lengths -- looks like there are "." entries
# for min and max lengths when only one fish was caught at a site, change that to the actual length values using
# data.frame[row_number, column_number] = new_value (might be an easier way but I'm not sure how)
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

# adult stage-specific records where EITHER min or max total lengths exceed the size at maturity (> 24.71 cm TL)
lg_adult = lg %>%   
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(MIN_TOT_LEN > 24.71 | MAX_TOT_LEN > 24.71) %>% # if either of these are true, then an adult was present
  summarize(N = sum(as.numeric(NO))) %>%
  mutate(SPECIES_CD = "LUT GRIS") %>% # to match RVC
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE) %>%
  ungroup()
lg_adult$LIFE_STAGE = rep("ADULT",nrow(lg_adult))
lg_adult$PRES = ifelse(lg_adult$N>0, 1, 0)
lg_adult$SOURCE = rep("MANGROVE VISUAL SURVEY",nrow(lg_adult))

# # now the inferred absences sites (sites where either no gray snapper were seen or only those of another age class were seen)
lg_adult_abs = lg %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN > 24.71 | MAX_TOT_LEN > 24.71)) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT adults
  ungroup()
lg_adult_abs$LIFE_STAGE = rep("ADULT", nrow(lg_adult_abs))
lg_adult_abs$PRES = ifelse(lg_adult_abs$N > 0, 1, 0) 
lg_adult_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(lg_adult_abs))

# combining the adult stage-specific data and inferred absence data, again using distinct to remove repeating values
lg_adult_mg = do.call("rbind", list(lg_adult, lg_adult_abs, lg_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)


# repeat for the subadult gray snapper stage (9.51 cm <= TOT LEN <= 24.71 cm (between size at 1 YR and size at maturation))
# remember first to change the lifestage column of the shared absence dataframe
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

# # now the inferred absences sites (sites where either no gray snapper were seen or only those of another age class were seen)
lg_subadult_abs = lg %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN >= 9.51 & MIN_TOT_LEN <= 24.71) | (MAX_TOT_LEN >= 9.51 & MAX_TOT_LEN <= 24.71)) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT adults
  ungroup()
lg_subadult_abs$LIFE_STAGE = rep("SUBADULT", nrow(lg_subadult_abs))
lg_subadult_abs$PRES = ifelse(lg_subadult_abs$N > 0, 1, 0) 
lg_subadult_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(lg_subadult_abs))

# combining the subadult stage-specific data and inferred absence data, again using distinct to remove repeating values
lg_subadult_mg = do.call("rbind", list(lg_subadult, lg_subadult_abs, lg_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# finally, repeat for juvenile gray snapper stage (TOT LEN < 9.51 cm (smaller than size at 1 YR), but not equal to 0!)
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

# # now the inferred absences sites (sites where either no gray snapper were seen or only those of another age class were seen)
lg_juvenile_abs = lg %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN != 0 & (MIN_TOT_LEN < 9.51 | MAX_TOT_LEN < 9.51))) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT adults
  ungroup()
lg_juvenile_abs$LIFE_STAGE = rep("JUVENILE", nrow(lg_juvenile_abs))
lg_juvenile_abs$PRES = ifelse(lg_juvenile_abs$N > 0, 1, 0) 
lg_juvenile_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(lg_juvenile_abs))

# combining the juvenile stage-specific data and inferred absence data, again using distinct to remove repeating values
lg_juvenile_mg = do.call("rbind", list(lg_juvenile, lg_juvenile_abs, lg_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# compiling all MG datasets just in case
lg_mg = do.call("rbind", list(lg_adult_mg,
                              lg_subadult_mg,
                              lg_juvenile_mg))

# and getting rid of all intermediate tables
rm(list = c("lg_adult", "lg_adult_abs", "lg_subadult", "lg_subadult_abs", "lg_juvenile", 
            "lg_juvenile_abs", "lg", "lg_abs", "lg_pres_sites"))

# now combining the two data sources - RVC and MG - and clipping to training and testing areas
# adults
train_lg_adult = full_join(lg_adult_rvc, lg_adult_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) %>% # re-project 
  st_intersection(., st_make_valid(train)) # clip
compareCRS(train_lg_adult, my_crs) # just to be sure

test_lg_adult = full_join(lg_adult_rvc, lg_adult_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% 
  st_transform(., my_crs) %>% 
  st_intersection(., st_make_valid(test)) 
compareCRS(test_lg_adult, my_crs) 

rm(list = c("lg_adult_mg", "lg_adult_rvc"))

# subadults
train_lg_subadult = full_join(lg_subadult_rvc, lg_subadult_mg)  %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) %>% # re-project 
  st_intersection(., st_make_valid(train)) # clip
compareCRS(train_lg_subadult, my_crs) # just to be sure

test_lg_subadult = full_join(lg_subadult_rvc, lg_subadult_mg)  %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% 
  st_transform(., my_crs) %>% 
  st_intersection(., st_make_valid(test)) 
compareCRS(test_lg_subadult, my_crs) 

rm(list = c("lg_subadult_mg", "lg_subadult_rvc"))

# juveniles
train_lg_juvenile = full_join(lg_juvenile_rvc, lg_juvenile_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) %>% # re-project 
  st_intersection(., st_make_valid(train)) # clip
compareCRS(train_lg_juvenile, my_crs) # just to be sure

test_lg_juvenile = full_join(lg_juvenile_rvc, lg_juvenile_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% 
  st_transform(., my_crs) %>% 
  st_intersection(., st_make_valid(test)) 
compareCRS(test_lg_juvenile, my_crs)

rm(list = c("lg_juvenile_mg", "lg_juvenile_rvc"))

# write out presence-absence data from training and testing sites, respectively 
st_write(train_lg_adult %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(train_wd, "Occurrences/Presence_Absence/Adult_Gray_Snapper_PA_Train.csv"), 
         append = FALSE)

st_write(train_lg_subadult %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(train_wd, "Occurrences/Presence_Absence/Subadult_Gray_Snapper_PA_Train.csv"),
         append = FALSE)

st_write(train_lg_juvenile %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(train_wd, "Occurrences/Presence_Absence/Juvenile_Gray_Snapper_PA_Train.csv"), 
         append = FALSE)

st_write(test_lg_adult %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(test_wd, "Occurrences/Presence_Absence/Adult_Gray_Snapper_PA_Test.csv"), 
         append = FALSE)

st_write(test_lg_subadult %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(test_wd, "Occurrences/Presence_Absence/Subadult_Gray_Snapper_PA_Test.csv"),
         append = FALSE)

st_write(test_lg_juvenile %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(test_wd, "Occurrences/Presence_Absence/Juvenile_Gray_Snapper_PA_Test.csv"), 
         append = FALSE)

# now filtering out only the presence records for MaxEnt modeling
st_write(train_lg_adult %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(train_wd, "Occurrences/Presence_Only/Adult_Gray_Snapper_PO_Train.csv"), 
         append = FALSE)

st_write(train_lg_subadult %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(train_wd, "Occurrences/Presence_Only/Subadult_Gray_Snapper_PO_Train.csv"), 
         append = FALSE)

st_write(train_lg_juvenile %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(train_wd, "Occurrences/Presence_Only/Juvenile_Gray_Snapper_PO_Train.csv"), 
         append = FALSE)

st_write(test_lg_adult %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(test_wd, "Occurrences/Presence_Only/Adult_Gray_Snapper_PO_Test.csv"), 
         append = FALSE)

st_write(test_lg_subadult %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(test_wd, "Occurrences/Presence_Only/Subadult_Gray_Snapper_PO_Test.csv"), 
         append = FALSE)

st_write(test_lg_juvenile %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(test_wd, "Occurrences/Presence_Only/Juvenile_Gray_Snapper_PO_Test.csv"), 
         append = FALSE)

# cleaning out the global environmental before moving on to bluestriped grunts
rm(list = c("adult_gray_snapper", "juvenile_gray_snapper", "subadult_gray_snapper"))

#### bluestriped grunts ####

# first convert fork length to total length using the bluestriped grunt
# length-length conversion from FishBase  TL = 0 + 1.034 x FL
rvc = rvc %>% mutate(TOT_LEN = (LEN*1.034)) 

# start with presence and absence records specific to the adult bluestriped grunt 
# stage (> 25.33 cm TL (size at maturity))
hs_adult = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% # group by both primary and second-stage sample units
  filter(SPECIES_CD == "HAE SCIU" & TOT_LEN > 25.33) %>% 
  summarize(N = sum(NUM)) %>% # was an adult bluestriped grunt ever seen at this second-stage sample unit (SSU)?
  ungroup() %>%
  distinct()
hs_adult$LIFE_STAGE = rep("ADULT", nrow(hs_adult)) # specify that these are all adult records
hs_adult$PRES = ifelse(hs_adult$N > 0, 1, 0) # if at least one adult was seen at the SSU, it's a presence [1], else an absence [0]
hs_adult$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_adult)) # source column for when RVC and MG data are compiled

# now the inferred absences sites (sites where either no gray bluestriped grunts 
# were seen or only those of another age class were seen)
hs_adult_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "HAE SCIU" & !TOT_LEN > 25.33) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 because these are adult absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
hs_adult_abs$LIFE_STAGE = rep("ADULT", nrow(hs_adult_abs))
hs_adult_abs$PRES = ifelse(hs_adult_abs$N > 0, 1, 0)
hs_adult_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_adult_abs))

# combining adult absence and presence data and using distinct because some of the absence sites might be shared across
# both the adult stage-specific data and the inferred absence data
hs_adult_rvc = rbind(hs_adult, hs_adult_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, .keep_all = T)

# repeat for the subadult bluestriped grunt stage (11.90 cm <= TOT LEN <= 25.33 cm (between size at 1 YR and size at maturation))
hs_subadult = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% 
  filter(SPECIES_CD == "HAE SCIU" & (TOT_LEN >= 11.90 & TOT_LEN <= 25.33)) %>% 
  summarize(N = sum(NUM)) %>% # was a subadult bluestriped grunt ever seen at this second-stage sample unit (SSU)?
  ungroup() %>%
  distinct()
hs_subadult$LIFE_STAGE = rep("SUBADULT", nrow(hs_subadult)) 
hs_subadult$PRES = ifelse(hs_subadult$N > 0, 1, 0)
hs_subadult$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_subadult))

# now the inferred absences sites (sites where either no bluestriped grunt were seen or only those of another age class were seen)
hs_subadult_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "HAE SCIU" & !(TOT_LEN >= 11.90 & TOT_LEN <= 25.33)) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 because these are subadult absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
hs_subadult_abs$LIFE_STAGE = rep("SUBADULT", nrow(hs_subadult_abs))
hs_subadult_abs$PRES = ifelse(hs_subadult_abs$N > 0, 1, 0)
hs_subadult_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_subadult_abs))

# combining subadult absence and presence data and using distinct because some of the absence sites might be shared across
# both the subadult stage-specific data and the inferred absence data
hs_subadult_rvc = rbind(hs_subadult, hs_subadult_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, .keep_all = T)

# finally, repeat for the juvenile bluestriped grunt stage (TOT LEN < 11.90 cm (smaller than size at 1 YR), but not equal to 0!)
hs_juvenile = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>% 
  filter(SPECIES_CD == "HAE SCIU" & (TOT_LEN != 0 & TOT_LEN < 11.90)) %>% 
  summarize(N = sum(NUM)) %>% # was a juvenile bluestriped grunt ever seen at this second-stage sample unit (SSU)?
  ungroup() %>%
  distinct()
hs_juvenile$LIFE_STAGE = rep("JUVENILE", nrow(hs_juvenile)) 
hs_juvenile$PRES = ifelse(hs_juvenile$N > 0, 1, 0)
hs_juvenile$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_juvenile))

# now the inferred absences sites (sites where either no bluestriped grunts were seen or only those of another age class were seen)
hs_juvenile_abs = rvc %>%
  group_by(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "HAE SCIU" & !(TOT_LEN != 0 & TOT_LEN < 11.90)) %>% 
  add_column(N = 0) %>% # manually assign a value of 0 because these are juvenile absence sites
  ungroup() %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, SPECIES_CD, N)
hs_juvenile_abs$LIFE_STAGE = rep("JUVENILE", nrow(hs_juvenile_abs))
hs_juvenile_abs$PRES = ifelse(hs_juvenile_abs$N > 0, 1, 0)
hs_juvenile_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_juvenile_abs))

# combining juvenile absence and presence data and using distinct because some of the absence sites might be shared across
# both the juvenile stage-specific data and the inferred absence data
hs_juvenile_rvc = rbind(hs_juvenile, hs_juvenile_abs) %>%
  distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, .keep_all = T)


# compiling all RVC datasets just in case
hs_rvc = do.call("rbind", list(hs_adult_rvc,
                               hs_subadult_rvc,
                               hs_juvenile_rvc))

# and getting rid of all intermediate tables
rm(list = c("hs_adult", "hs_adult_abs", "hs_subadult", "hs_subadult_abs", "hs_juvenile", 
            "hs_juvenile_abs"))

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


# adult stage-specific records where EITHER min or max total lengths exceed the size at maturity (> 24.71 cm TL)
hs_adult = hs %>%   
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(MIN_TOT_LEN > 25.33 | MAX_TOT_LEN > 25.33) %>% # if either of these are true, then an adult was present
  summarize(N = sum(as.numeric(NO))) %>%
  mutate(SPECIES_CD = "HAE SCIU") %>% # to match RVC
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE) %>%
  ungroup()
hs_adult$LIFE_STAGE = rep("ADULT",nrow(hs_adult))
hs_adult$PRES = ifelse(hs_adult$N>0, 1, 0)
hs_adult$SOURCE = rep("MANGROVE VISUAL SURVEY",nrow(hs_adult))

# # now the inferred absences sites (sites where either no bluestriped grunts were seen or only those of another age class were seen)
hs_adult_abs = hs %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN > 25.33 | MAX_TOT_LEN > 25.33)) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT adults
  ungroup()
hs_adult_abs$LIFE_STAGE = rep("ADULT", nrow(hs_adult_abs))
hs_adult_abs$PRES = ifelse(hs_adult_abs$N > 0, 1, 0) 
hs_adult_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(hs_adult_abs))

# combining the adult stage-specific data and inferred absence data, again using distinct to remove repeating values
hs_adult_mg = do.call("rbind", list(hs_adult, hs_adult_abs, hs_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# repeat for the subadult bluestriped grunt stage (11.90 cm <= TOT LEN <= 25.33 cm (between size at 1 YR and size at maturation))
# remember first to change the lifestage column of the shared absence dataframe
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

# # now the inferred absences sites (sites where either no bluestriped grunts were seen or only those of another age class were seen)
hs_subadult_abs = hs %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN >= 11.90 & MIN_TOT_LEN <= 25.33) | (MAX_TOT_LEN >= 11.90 & MAX_TOT_LEN <= 25.33)) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT adults
  ungroup()
hs_subadult_abs$LIFE_STAGE = rep("SUBADULT", nrow(hs_subadult_abs))
hs_subadult_abs$PRES = ifelse(hs_subadult_abs$N > 0, 1, 0) 
hs_subadult_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(hs_subadult_abs))

# combining the subadult stage-specific data and inferred absence data, again using distinct to remove repeating values
hs_subadult_mg = do.call("rbind", list(hs_subadult, hs_subadult_abs, hs_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# finally, repeat for juvenile bluestriped grunt stage (TOT LEN < 11.90 cm (smaller than size at 1 YR), but not equal to 0!)
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

# # now the inferred absences sites (sites where either no bluestriped grunt were seen or only those of another age class were seen)
hs_juvenile_abs = hs %>%  
  group_by(Site, LAT_DEGREES, LON_DEGREES, SPECIES_CD) %>%
  filter(!(MIN_TOT_LEN != 0 & (MIN_TOT_LEN < 11.90 | MAX_TOT_LEN < 11.90))) %>%
  distinct(Site, LAT_DEGREES, LON_DEGREES) %>%
  add_column("N" = 0) %>% # assign 0 because these are NOT adults
  ungroup()
hs_juvenile_abs$LIFE_STAGE = rep("JUVENILE", nrow(hs_juvenile_abs))
hs_juvenile_abs$PRES = ifelse(hs_juvenile_abs$N > 0, 1, 0) 
hs_juvenile_abs$SOURCE = rep("MANGROVE VISUAL SURVEY", nrow(hs_juvenile_abs))

# combining the juvenile stage-specific data and inferred absence data, again using distinct to remove repeating values
hs_juvenile_mg = do.call("rbind", list(hs_juvenile, hs_juvenile_abs, hs_abs)) %>%
  distinct(Site, LON_DEGREES, LAT_DEGREES, .keep_all = TRUE)

# compiling all MG datasets just in case
hs_mg = do.call("rbind", list(hs_adult_mg,
                              hs_subadult_mg,
                              hs_juvenile_mg))

# and getting rid of all intermediate tables
rm(list = c("hs_adult", "hs_adult_abs", "hs_subadult", "hs_subadult_abs", "hs_juvenile", 
            "hs_juvenile_abs", "hs", "hs_abs", "hs_pres_sites"))

# now combining the two data sources - RVC and MG - and clipping to training and testing areas
# adults
train_hs_adult = full_join(hs_adult_rvc, hs_adult_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) %>% # re-project 
  st_intersection(., st_make_valid(train)) # clip
compareCRS(train_hs_adult, my_crs) # just to be sure

test_hs_adult = full_join(hs_adult_rvc, hs_adult_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% 
  st_transform(., my_crs) %>% 
  st_intersection(., st_make_valid(test)) 
compareCRS(test_hs_adult, my_crs) 

rm(list = c("hs_adult_mg", "hs_adult_rvc"))

# subadults
train_hs_subadult = full_join(hs_subadult_rvc, hs_subadult_mg)  %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) %>% # re-project 
  st_intersection(., st_make_valid(train)) # clip
compareCRS(train_hs_adult, my_crs) # just to be sure

test_hs_subadult = full_join(hs_subadult_rvc, hs_subadult_mg)  %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% 
  st_transform(., my_crs) %>% 
  st_intersection(., st_make_valid(test)) 
compareCRS(test_hs_adult, my_crs) 

rm(list = c("hs_subadult_mg", "hs_subadult_rvc"))

# juveniles
train_hs_juvenile = full_join(hs_juvenile_rvc, hs_juvenile_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% # convert to spatial (lon, lat)
  st_transform(., my_crs) %>% # re-project 
  st_intersection(., st_make_valid(train)) # clip
compareCRS(train_hs_juvenile, my_crs) # just to be sure

test_hs_juvenile = full_join(hs_juvenile_rvc, hs_juvenile_mg) %>%
  select(SPECIES_CD, LON_DEGREES, LAT_DEGREES, LIFE_STAGE, PRES, N, SOURCE) %>%
  st_as_sf(., coords = c(2, 3), crs = gcs) %>% 
  st_transform(., my_crs) %>% 
  st_intersection(., st_make_valid(test)) 
compareCRS(test_hs_juvenile, my_crs)

rm(list = c("hs_juvenile_mg", "hs_juvenile_rvc"))

# write out presence-absence data from training and testing sites, respectively 
st_write(train_hs_adult %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(train_wd, "Occurrences/Presence_Absence/Adult_Bluestriped_Grunt_PA_Train.csv"), 
         append = FALSE)

st_write(train_hs_subadult %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(train_wd, "Occurrences/Presence_Absence/Subadult_Bluestriped_Grunt_PA_Train.csv"),
         append = FALSE)

st_write(train_hs_juvenile %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(train_wd, "Occurrences/Presence_Absence/Juvenile_Bluestriped_Grunt_PA_Train.csv"), 
         append = FALSE)

st_write(test_hs_adult %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(test_wd, "Occurrences/Presence_Absence/Adult_Bluestriped_Grunt_PA_Test.csv"), 
         append = FALSE)

st_write(test_hs_subadult %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(test_wd, "Occurrences/Presence_Absence/Subadult_Bluestriped_Grunt_PA_Test.csv"),
         append = FALSE)

st_write(test_hs_juvenile %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]),
         dsn = paste0(test_wd, "Occurrences/Presence_Absence/Juvenile_Bluestriped_Grunt_PA_Test.csv"),
         append = FALSE)

# now filtering out only the presence records for MaxEnt modeling
st_write(train_hs_adult %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(train_wd, "Occurrences/Presence_Only/Adult_Bluestriped_Grunt_PO_Train.csv"), 
         append = FALSE)

st_write(train_hs_subadult %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(train_wd, "Occurrences/Presence_Only/Subadult_Bluestriped_Grunt_PO_Train.csv"),
         append = FALSE)

st_write(train_hs_juvenile %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(train_wd, "Occurrences/Presence_Only/Juvenile_Bluestriped_Grunt_PO_Train.csv"), 
         append = FALSE)

st_write(test_hs_adult %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(test_wd, "Occurrences/Presence_Only/Adult_Bluestriped_Grunt_PO_Test.csv"), 
         append = FALSE)

st_write(test_hs_subadult %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(test_wd, "Occurrences/Presence_Only/Subadult_Bluestriped_Grunt_PO_Test.csv"), 
         append = FALSE)

st_write(test_hs_juvenile %>%
           filter(PRES == 1) %>%
           mutate(LON_M = sf::st_coordinates(.)[,1],
                  LAT_M = sf::st_coordinates(.)[,2]) %>%
           select(SPECIES_CD, LON_M, LAT_M), # maxent only needs species name, lon, lat (in that order)
         dsn = paste0(test_wd, "Occurrences/Presence_Only/Juvenile_Bluestriped_Grunt_PO_Test.csv"), 
         append = FALSE)

# cleaning out the global environmental before moving on to bluestriped grunts
rm(list = c("adult_bluestriped_grunt", "juvenile_bluestriped_grunt", "subadult_bluestriped_grunt"))


# create sampling effort raster for training area to parse out sampling bias
library(spatialEco)

# rule of thumb for selecting bandwidth according to Scott (1992) and Bowman and Azzalini (1997)
choose_bw = function(spdf) {
  X = coordinates(spdf)
  sigma = c(sd(X[,1]), sd(X[,2])) * (2 / (3 * nrow(X))) ^ (1/6)
}

sampling_effort_train = full_join(rvc_sites, mg_sites) %>%
  select(LON_DEGREES, LAT_DEGREES) %>%
  st_as_sf(., coords = c(1, 2), crs = gcs) %>% 
  st_transform(., my_crs) %>% 
  st_intersection(., st_make_valid(train)) %>%
  add_column("count" = 1) %>%
  as(., "Spatial")

train_grid = raster(ncol = ncol(habitat_train), nrow = nrow(habitat_train), 
                    xmn = xmin(habitat_train), xmx = xmax(habitat_train), 
                    ymn = ymin(habitat_train), ymx = ymax(habitat_train))
train_grid[] <- rep(1,ncell(train_grid))
train_bw = choose_bw(sampling_effort_train)

train_kde = sp.kde(x = sampling_effort_train,
                   bw = train_bw,
                   newdata = train_grid,
                   standardize = T)
tm_shape(train_kde) + tm_raster() 
train_kde = raster::mask(raster::crop(train_kde, habitat_train), habitat_train)
compareRaster(train_kde, habitat_train, extent = T, crs = T, rowcol = T)
tm_shape(train_kde) + tm_raster(n = 5, palette = "-RdBu") #+ tm_shape(sampling_effort_train) + tm_dots(size = 0.01)
writeRaster(train_kde, filename = paste0(train_wd, "Occurrences/bias.asc"),
            format = "ascii", overwrite = T)

# free up some space
rm(list = c("train_lg_adult", "train_lg_subadult", "train_lg_juvenile",
            "test_lg_adult", "test_lg_subadult", "test_lg_juvenile",
            "train_hs_adult", "train_hs_subadult", "train_hs_juvenile",
            "test_hs_adult", "test_hs_subadult", "test_hs_juvenile",
            "rvc", "mg", "rvc_sites", "mg_sites", "lg_mg", "lg_rvc",
            "hs_mg", "hs_rvc", "mg_euc_train", "mg_euc_test", "train_kde"))

#### WATER QUALITY DATA ####
libraries("sp", "tibble", "ncf", "spdep", "gstat", "geoR", "readxl", "tidyr")
rasterOptions(tmpdir = "Z:/Courtney/GIS/Temp/") # temporary working directory for raster package

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
summary(winter_wq_sp)
write.csv(winter_wq_sp, paste0(csv_wd, "Winter_Water_Conditions.csv"))

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

# according to Moran tests with both normal and Monte Carlo approximations, there is significant spatial dependence 
# in the winter temp, sal, and DO data.
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
write.csv(summer_wq_sp, paste0(csv_wd, "Summer_Water_Conditions.csv"))

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
# now move on to variogram modeling to capture the spatial structure...
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
# prediction grid (based on habitat raster, here it is if you need to reload)
# habitat_train = raster(paste0(train_wd, "Environmental/Habitat.asc")) 
# crs(habitat_train) = my_crs

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

#### winter temperature ####
# Calculate the number of cores
no_cores = detectCores() - 2

# Initiate cluster 
cl = makeCluster(no_cores)

# split training area into pieces for each core
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

# save the new surface
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

# save the new surface
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

# save the new surface
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

# save the new surface
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

# save the new surface
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

# save the new surface
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
  krige(formula = MEAN_WIN_TEMP_B ~ 1, locations = winter_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = w_temp_fvgm))
stopCluster(cl)

# merge all the predictions
w_temp_test_merge = maptools::spRbind(w_temp_test_par[[1]], w_temp_test_par[[2]])
for (j in 3:length(w_temp_test_par)) {
  w_temp_test_merge = maptools::spRbind(w_temp_test_merge, w_temp_test_par[[j]])
}
w_temp_test_merge = SpatialPixelsDataFrame(points = w_temp_test_merge, data = w_temp_test_merge@data)

# save the new surface
summary(w_temp_test_merge)
writeGDAL(w_temp_test_merge["var1.pred"], fname = paste0(temp_wd, "mean_win_temp_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_win_temp_test = raster(paste0(temp_wd, "mean_win_temp_test.tif"))
mean_win_temp_test = writeRaster(raster::mask(raster::crop(mean_win_temp_test, habitat_test), habitat_test),
                                 file = file.path(test_wd, "Environmental/Mean_Win_Temp.asc"), format = "ascii",
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
  krige(formula = MEAN_WIN_SAL_B ~ 1, locations = winter_wq_sp,
        newdata = test_grid[test_parts[[x]],], model = w_sal_fvgm))
stopCluster(cl)

# merge all the predictions
w_sal_test_merge = maptools::spRbind(w_sal_test_par[[1]], w_sal_test_par[[2]])
for (j in 3:length(w_sal_test_par)) {
  w_sal_test_merge = maptools::spRbind(w_sal_test_merge, w_sal_test_par[[j]])
}
w_sal_test_merge = SpatialPixelsDataFrame(points = w_sal_test_merge, data = w_sal_test_merge@data)

# save the new surface
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

# save the new surface
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

# save the new surface
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

# save the new surface
summary(s_do_test_merge)
writeGDAL(s_do_test_merge["var1.pred"], fname = paste0(temp_wd, "mean_sum_do_test.tif"),
          drivername = "GTiff", type = "Float32")
mean_sum_do_test = raster(paste0(temp_wd, "mean_sum_do_test.tif"))
mean_sum_do_test = writeRaster(raster::mask(raster::crop(mean_sum_do_test, habitat_test), habitat_test),
                               file = file.path(test_wd, "Environmental/Mean_Sum_DO.asc"), format = "ascii",
                               overwrite = T)

rm(list = c("s_do_test_merge", "s_do_test_par", "mean_sum_do_test", "cl"))
showConnections()
