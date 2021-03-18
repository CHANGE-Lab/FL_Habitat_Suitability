# Logistic Regression modeling of habitat suitability for sub-adult gray snapper and
# bluestriped grunts in the Florida Keys. For this research, the Florida Keys region
# was divided into two coastal zones, a training area and a testing area. The training
# area included Biscayne National Park (BNP) and a portion of the Florida Keys National
# Marine Sanctuary (FKNMS) extending from its northernmost boundary to Tavernier Creek. 
# The testing area was entirely encompassed by the FKNMS, bounded by Tavernier Creek
# to the northeast and Cudjoe Key to the southwest. 

# working directories
setwd("Z:/Courtney/Stuart_MSc_Ch1/") # main project folder
temp_wd = "Z:/Courtney/Stuart_MSc_Ch1/Temp/" # temporary files
source_wd = "Z:/Courtney/Stuart_MSc_Ch1/Source_Data/" # source data
dem_wd = "Z:/Courtney/Stuart_MSc_Ch1/Source_Data/Job606638_ncei_nintharcsec_dem/" # DEMs
csv_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Data/" # for writing tabular data
train_wd = "Z:/Courtney/Stuart_MSc_Ch1/Modeling_Data/Training/" # for training area data
test_wd = "Z:/Courtney/Stuart_MSc_Ch1/Modeling_Data/Testing/" # for testing area data
plots_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Figures/"


# libraries
library(easypackages)
libraries("rgdal", "gdalUtils", "raster", "sp", "sf", "tmap", "dplyr", "lwgeom", "rgeos",
          "cleangeo", "tidyverse", "stars", "fasterize", "PNWColors", "spex", "igraph", 
          "spatialEco", "MLeval")

# change where large temp rasters are saved
rasterOptions(tmpdir = "Z:/Courtney/Stuart_MSc_Ch1/Temp/")

# save PROJ.4 string for standard projection (ESPG:26958 NAD 83/Florida East) and 
# source LAT/LON data (EPSG:4326 WGS 84/World Geodetic System 1984)
my_crs = CRS("+init=epsg:26958")
gcs = CRS("+init=epsg:4326")

#### SET-UP: TRAINING AREA ####

# presence-absence sub-adult gray snapper (Lutjanus griseus) and bluestriped grunt 
# (Haemulon sciurus) records in the training area
lg_occ_train = read.csv(paste0(train_wd, "Occurrences/Presence_Absence/Subadult_Gray_Snapper_PA_Train.csv"))[,-6]
lg_coords_train = lg_occ_train %>% select(LON_M, LAT_M) %>%
  st_as_sf(., coords = c("LON_M", "LAT_M"), crs = my_crs)

hs_occ_train = read.csv(paste0(train_wd, "Occurrences/Presence_Absence/Subadult_Bluestriped_Grunt_PA_Train.csv"))[,-6]
hs_coords_train = hs_occ_train %>% select(LON_M, LAT_M) %>%
  st_as_sf(., coords = c("LON_M", "LAT_M"), crs = my_crs)

# spatial predictors in training area
habitat_train = raster(paste0(train_wd, "Environmental/Habitat.asc"))
mg_dist_train = raster(paste0(train_wd, "Environmental/Mangrove_Dist.asc"))
depth_train = raster(paste0(train_wd, "Environmental/Depth.asc"))
sd_depth_train = raster(paste0(train_wd, "Environmental/SD_Depth.asc"))
slope_train = raster(paste0(train_wd, "Environmental/Slope.asc"))
curvature_train = raster(paste0(train_wd, "Environmental/Curvature.asc"))
plan_curve_train = raster(paste0(train_wd, "Environmental/Plan_Curve.asc"))
bpi_fine_train = raster(paste0(train_wd, "Environmental/BPI_Fine.asc"))
bpi_broad_train = raster(paste0(train_wd, "Environmental/BPI_Broad.asc"))
rugosity_train = raster(paste0(train_wd, "Environmental/Rugosity.asc"))
sum_temp_train = raster(paste0(train_wd, "Environmental/Sum_Temp.asc"))
sum_do_train = raster(paste0(train_wd, "Environmental/Sum_DO.asc"))
sum_sal_train = raster(paste0(train_wd, "Environmental/Sum_Sal.asc"))
win_temp_train = raster(paste0(train_wd, "Environmental/Win_Temp.asc"))
win_do_train = raster(paste0(train_wd, "Environmental/Win_DO.asc"))
win_sal_train = raster(paste0(train_wd, "Environmental/Win_Sal.asc"))

# define crs
crs(habitat_train) = my_crs
crs(mg_dist_train) = my_crs
crs(depth_train) = my_crs
crs(sd_depth_train) = my_crs
crs(slope_train) = my_crs
crs(curvature_train) = my_crs
crs(plan_curve_train) = my_crs
crs(bpi_fine_train) = my_crs
crs(bpi_broad_train) = my_crs
crs(rugosity_train) = my_crs
crs(sum_temp_train) = my_crs
crs(sum_do_train) = my_crs
crs(sum_sal_train) = my_crs
crs(win_temp_train) = my_crs
crs(win_do_train) = my_crs
crs(win_sal_train) = my_crs

# create raster stack for training area
env_train = stack(x = c(habitat_train, mg_dist_train, depth_train, sd_depth_train,
                        slope_train, curvature_train, plan_curve_train, bpi_fine_train,
                        bpi_broad_train, rugosity_train, sum_temp_train, sum_do_train, 
                        sum_sal_train, win_temp_train, win_do_train, win_sal_train))


# combine fish and environmental data in the training area and save the final data frames
lg_train = cbind(lg_occ_train, raster::extract(env_train, lg_coords_train)) 
lg_train$Habitat = factor(lg_train$Habitat) # convert habitat to factor
# re-label response variable (1 = presence, 0 = absence) & convert to
# factor, otherwise caret will return errors later
lg_train$PRES[lg_train$PRES == 1] = "PRESENCE"
lg_train$PRES[lg_train$PRES == 0] = "ABSENCE"
lg_train$PRES = as.factor(lg_train$PRES)
write.csv(lg_train, paste0(csv_wd, "Gray_Snapper_Training_Data.csv"), row.names = F)


hs_train = cbind(hs_occ_train, raster::extract(env_train, hs_coords_train)) 
hs_train$Habitat = factor(hs_train$Habitat) # convert habitat to factor
# re-label response variable (1 = presence, 0 = absence) & convert to
# factor, otherwise caret will return errors later
hs_train$PRES[hs_train$PRES == 1] = "PRESENCE"
hs_train$PRES[hs_train$PRES == 0] = "ABSENCE"
hs_train$PRES = as.factor(hs_train$PRES)
write.csv(hs_train, paste0(csv_wd, "Bluestriped_Grunt_Training_Data.csv"), row.names = F)


#### SET-UP: TESTING AREA ####

# presence-absence sub-adult gray snapper (Lutjanus griseus) and bluestriped grunt 
# (Haemulon sciurus) records in the testing area
lg_occ_test = read.csv(paste0(test_wd, "Occurrences/Presence_Absence/Subadult_Gray_Snapper_PA_Test.csv"))[,-6]
lg_coords_test = lg_occ_test %>% select(LON_M, LAT_M) %>%
  st_as_sf(., coords = c("LON_M", "LAT_M"), crs = my_crs)

hs_occ_test = read.csv(paste0(test_wd, "Occurrences/Presence_Absence/Subadult_Bluestriped_Grunt_PA_Test.csv"))[,-6]
hs_coords_test = hs_occ_test %>% select(LON_M, LAT_M) %>%
  st_as_sf(., coords = c("LON_M", "LAT_M"), crs = my_crs)

# spatial predictors in testing area
habitat_test = raster(paste0(test_wd, "Environmental/Habitat.asc"))
mg_dist_test = raster(paste0(test_wd, "Environmental/Mangrove_Dist.asc"))
depth_test = raster(paste0(test_wd, "Environmental/Depth.asc"))
sd_depth_test = raster(paste0(test_wd, "Environmental/SD_Depth.asc"))
slope_test = raster(paste0(test_wd, "Environmental/Slope.asc"))
curvature_test = raster(paste0(test_wd, "Environmental/Curvature.asc"))
plan_curve_test = raster(paste0(test_wd, "Environmental/Plan_Curve.asc"))
bpi_fine_test = raster(paste0(test_wd, "Environmental/BPI_Fine.asc"))
bpi_broad_test = raster(paste0(test_wd, "Environmental/BPI_Broad.asc"))
rugosity_test = raster(paste0(test_wd, "Environmental/Rugosity.asc"))
sum_temp_test = raster(paste0(test_wd, "Environmental/Sum_Temp.asc"))
sum_do_test = raster(paste0(test_wd, "Environmental/Sum_DO.asc"))
sum_sal_test = raster(paste0(test_wd, "Environmental/Sum_Sal.asc"))
win_temp_test = raster(paste0(test_wd, "Environmental/Win_Temp.asc"))
win_do_test = raster(paste0(test_wd, "Environmental/Win_DO.asc"))
win_sal_test = raster(paste0(test_wd, "Environmental/Win_Sal.asc"))

# define crs
crs(habitat_test) = my_crs
crs(mg_dist_test) = my_crs
crs(depth_test) = my_crs
crs(sd_depth_test) = my_crs
crs(slope_test) = my_crs
crs(curvature_test) = my_crs
crs(plan_curve_test) = my_crs
crs(bpi_fine_test) = my_crs
crs(bpi_broad_test) = my_crs
crs(rugosity_test) = my_crs
crs(sum_temp_test) = my_crs
crs(sum_do_test) = my_crs
crs(sum_sal_test) = my_crs
crs(win_temp_test) = my_crs
crs(win_do_test) = my_crs
crs(win_sal_test) = my_crs

# create raster stack for testing area
env_test = stack(x = c(habitat_test, mg_dist_test, depth_test, sd_depth_test,
                        slope_test, curvature_test, plan_curve_test, bpi_fine_test,
                        bpi_broad_test, rugosity_test, sum_temp_test, sum_do_test, 
                        sum_sal_test, win_temp_test, win_do_test, win_sal_test))


# combine fish and environmental data in the testing area
lg_test = cbind(lg_occ_test, raster::extract(env_test, lg_coords_test)) 
lg_test$Habitat = factor(lg_test$Habitat) # convert habitat to factor
# re-label response variable (1 = presence, 0 = absence) & convert to
# factor, otherwise caret will return errors later
lg_test$PRES[lg_test$PRES == 1] = "PRESENCE"
lg_test$PRES[lg_test$PRES == 0] = "ABSENCE"
lg_test$PRES = as.factor(lg_test$PRES)
write.csv(lg_test, paste0(csv_wd, "Gray_Snapper_Testing_Data.csv"), row.names = F)


hs_test = cbind(hs_occ_test, raster::extract(env_test, hs_coords_test)) 
hs_test$Habitat = factor(hs_test$Habitat) # convert habitat to factor
# re-label response variable (1 = presence, 0 = absence) & convert to
# factor, otherwise caret will return errors later
hs_test$PRES[hs_test$PRES == 1] = "PRESENCE"
hs_test$PRES[hs_test$PRES == 0] = "ABSENCE"
hs_test$PRES = as.factor(hs_test$PRES)
write.csv(hs_test, paste0(csv_wd, "Bluestriped_Grunt_Testing_Data.csv"), row.names = F)


#### CARET: TRAIN GRAY SNAPPER MODEL ####
library(caret)
set.seed(42)

# set training method for caret (10 fold cross-validation)
trControl = trainControl(method = "cv", number = 10, savePredictions = TRUE,
                         summaryFunction = twoClassSummary, classProbs = TRUE)

# full model with all 11 predictors
lg_fit1 = caret::train(PRES ~ Habitat + Mangrove_Dist + Depth + Slope + 
                         Curvature + BPI_Fine + BPI_Broad + Rugosity + Sum_Temp +
                         Sum_Sal + Win_Sal, data = lg_train,
                       method = "glm", family = "binomial", metric = "ROC",
                       trControl = trControl, na.action = na.pass)
lg_fit1
summary(lg_fit1)

# AUC-ROC Curve for training CV
lg_eval = evalm(lg_fit1, percent = 95, gnames = "Lujanus griseus", 
                rlinethick = 0.75, fsize = 10)
lg_eval$roc
png(filename = paste0(plots_wd, "Gray_Snapper_LogReg_AUC.png"), width = 5, height = 3.15, units = "in", 
    res = 300)
lg_eval$roc
dev.off()

# variable importance
varImp(lg_fit1, useModel = T, scale = F)
plot((varImp(lg_fit1, useModel = T, scale = F)))

# look at variable importance in terms of AUC
lg_df = lg_train %>% select(PRES, Habitat, Mangrove_Dist, Depth, Slope,
                            Rugosity, Curvature, BPI_Fine, BPI_Broad, 
                            Sum_Temp, Sum_Sal, Win_Sal)

lg_AUC_imp = filterVarImp(lg_df[,2-12], lg_df[,1])
lg_AUC_imp

# create predicted HSM raster across entire training area
beginCluster(n = 10)
lg_train_HSM = raster::predict(env_train, lg_fit1)
endCluster()

# binary presence/absence predictions for testing area
lg_test_pred = predict(lg_fit1, lg_test)
confusionMatrix(lg_test_pred, lg_test$PRES, positive = "PRESENCE")

# predictions if using cut-off of 0.5 for suitable 
lg_test_pred2 = predict(lg_fit1, lg_test, type = "prob")
new_lg_test_pred = ifelse(lg_test_pred2$PRESENCE >= 0.4, "PRESENCE", "ABSENCE")
new_lg_test_pred = as.factor(new_lg_test_pred)
confusionMatrix(new_lg_test_pred, lg_test$PRES, positive = "PRESENCE")



# full model with all 12 predictors
hs_fit1 = caret::train(PRES ~ Habitat + Mangrove_Dist + Depth + Slope + 
                         Curvature + BPI_Fine + BPI_Broad + Rugosity + Sum_Temp +
                         Sum_Sal + Win_Sal, data = hs_train,
                       method = "glm", family = "binomial", metric = "ROC",
                       trControl = trControl)
hs_fit1
summary(hs_fit1)



# AUC-ROC Curve for training CV
hs_eval = evalm(hs_fit1, percent = 95, gnames = "Haemulon sciurus", 
                rlinethick = 0.75, fsize = 10)
hs_eval$roc
png(filename = paste0(plots_wd, "Bluestriped_Grunt_LogReg_AUC.png"), width = 5, height = 3.15, units = "in", 
    res = 300)
hs_eval$roc
dev.off()

# variable importance
varImp(hs_fit1, useModel = TRUE, scale = F)

hs_df = hs_train %>% select(PRES, Habitat, Mangrove_Dist, Depth, Slope,
                            Rugosity, Curvature, BPI_Fine, BPI_Broad, 
                            Sum_Temp, Sum_Sal, Win_Sal)

hs_AUC_imp = filterVarImp(hs_df[,2-12], hs_df[,1])
hs_AUC_imp

# create predicted HSM raster across entire training area
beginCluster(n = 10)
hs_train_HSM = raster::predict(env_train, hs_fit1)
endCluster()

# predictions for testing area (as is)
hs_test_pred = predict(hs_fit1, hs_test)
confusionMatrix(hs_test_pred, hs_test$PRES, positive = "PRESENCE")
view(hs_test_pred)
# predictions if using cut-off of 0.4 for suitable 
hs_test_pred2 = predict(hs_fit1, hs_test, type = "prob")
new_hs_test_pred = ifelse(hs_test_pred2$PRESENCE >= 0.4, "PRESENCE", "ABSENCE")
new_hs_test_pred = as.factor(new_hs_test_pred)
confusionMatrix(new_hs_test_pred, hs_test$PRES, positive = "PRESENCE")
