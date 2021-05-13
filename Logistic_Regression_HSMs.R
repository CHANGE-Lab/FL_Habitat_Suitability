# This script is part of Courtney Stuart's first MSc Chapter in the Lab of 
# Dr. Stephanie Green at the University of Alberta. Data are specific to southern
# Florida, including the georeferenced occurrence records of two reef fish species
# (gray snapper Lutjanus griseus (lg) and bluestriped grunt Haemulon sciurus (hs))
# and rasters of various environmental covariates. Three logistic regression models
# were fitted for each species to create habitat suitability predictions: standard 
# logistic regression, lasso-regularized regression, and ridge-regularized regression. 


#### SET-UP ####

# working directories
setwd("Z:/Courtney/Stuart_MSc_Ch1/") # main project folder

# data directories 
temp_wd = "Z:/Courtney/Stuart_MSc_Ch1/Temporary/" # temporary files
csv_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Data/" # for writing tabular data
fish_wd = "Z:/Courtney/Stuart_MSc_Ch1/Species_Occurrence/" # for fish data
spatial_wd = "Z:/Courtney/Stuart_MSc_Ch1/Spatial_Predictors/" # for spatial predictor rasters
gis_wd = "Z:/Courtney/Stuart_MSc_Ch1/GIS_Files/" # for any GIS shapefiles
plots_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Figures/"
temp_plots = "Z:/Courtney/Stuart_MSc_Ch1/Plots/" # temp plots for Courtney

# libraries
library(easypackages)
libraries("glmnet", "tidyr", "caret", "MLeval", "dplyr", "raster", 
          "PNWColors", "cowplot", "sf", "conflicted")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# save PROJ.4 string for standard projection (ESPG:26958 NAD 83/Florida East) and source LAT/LON data (EPSG:4326 WGS 84/World Geodetic System 1984)
my_crs = CRS("+init=epsg:26958")

# extract data values from the 12 spatial predictors at the point locations
# of each presence-absence record (following multicollinearity assessment)
# subadult gray snapper (Lutjanus griseus)
lg_occ = read.csv(paste0(fish_wd, "Presence_Absence/Subadult_Gray_Snapper_PA.csv"))
lg_coords = lg_occ %>% dplyr::select(LON_M, LAT_M) %>%
  st_as_sf(., coords = c("LON_M", "LAT_M"), crs = my_crs)

# subadult bluestriped grunt
hs_occ = read.csv(paste0(fish_wd, "Presence_Absence/Subadult_Bluestriped_Grunt_PA.csv"))
hs_coords = hs_occ %>% dplyr::select(LON_M, LAT_M) %>%
  st_as_sf(., coords = c("LON_M", "LAT_M"), crs = my_crs)

# load environmental rasters 
habitat = raster(paste0(spatial_wd, "Habitat.asc"))
mg_dist = raster(paste0(spatial_wd, "Mangrove_Dist.asc"))
depth = raster(paste0(spatial_wd, "Depth.asc"))
slope = raster(paste0(spatial_wd, "Slope.asc"))
curvature = raster(paste0(spatial_wd, "Curvature.asc"))
bpi_fine = raster(paste0(spatial_wd, "BPI_Fine.asc"))
bpi_broad = raster(paste0(spatial_wd, "BPI_Broad.asc"))
rugosity = raster(paste0(spatial_wd, "Rugosity.asc"))
sum_temp = raster(paste0(spatial_wd, "Mean_Sum_Temp.asc"))
sum_sal = raster(paste0(spatial_wd, "Mean_Sum_Sal.asc"))
win_temp = raster(paste0(spatial_wd, "Mean_Win_Temp.asc"))
win_sal = raster(paste0(spatial_wd, "Mean_Win_Sal.asc"))

# define crs
crs(habitat) = my_crs
crs(mg_dist) = my_crs
crs(depth) = my_crs
crs(slope) = my_crs
crs(curvature) = my_crs
crs(bpi_fine) = my_crs
crs(bpi_broad) = my_crs
crs(rugosity) = my_crs
crs(sum_temp) = my_crs
crs(sum_sal) = my_crs
crs(win_temp) = my_crs
crs(win_sal) = my_crs

# create raster stack 
env = stack(x = c(habitat, mg_dist, depth, slope, curvature, bpi_fine, 
            bpi_broad, rugosity, sum_temp, sum_sal, win_temp, win_sal))

# combine fish records and environmental data
lg_full = cbind(lg_occ, raster::extract(env, lg_coords)) 
lg_full = as.data.frame(lg_full) %>%
  mutate(PRES = as.factor(PRES),
         Habitat = as.factor(Habitat))
names(lg_full)[names(lg_full) == "slope"] = "Slope"
names(lg_full)[names(lg_full) == "bpi_fine"] = "BPI_Fine"
names(lg_full)[names(lg_full) == "rugosity"] = "Rugosity"
# re-label response variable (1 = presence, 0 = absence) & convert to
# factor, otherwise caret will return errors later
lg_full$PRES2 = as.factor(ifelse(lg_full$PRES == 1, "PRESENCE", "ABSENCE"))
lg_full = lg_full %>% relocate(PRES2, .after = PRES)
# save data frame to csv
write.csv(lg_full, paste0(csv_wd, "Subadult_Gray_Snapper_Full_Dataset.csv"), row.names = F)

# repeat for subadult bluestriped grunts (Haemulon sciurus)
hs_full = cbind(hs_occ, raster::extract(env, hs_coords))
hs_full = as.data.frame(hs_full) %>%
  mutate(PRES = as.factor(PRES),
         Habitat = as.factor(Habitat))
names(hs_full)[names(hs_full) == "slope"] = "Slope"
names(hs_full)[names(hs_full) == "bpi_fine"] = "BPI_Fine"
names(hs_full)[names(hs_full) == "rugosity"] = "Rugosity"
# re-label response variable (1 = presence, 0 = absence) & convert to
# factor, otherwise caret will return errors later
hs_full$PRES2 = as.factor(ifelse(hs_full$PRES == 1, "PRESENCE", "ABSENCE"))
hs_full = hs_full %>% relocate(PRES2, .after = PRES)
# save data frame to csv
write.csv(hs_full, paste0(csv_wd, "Subadult_Bluestriped_Grunt_Full_Dataset.csv"), row.names = F)

# random 70-30% data split for model training and evaluation
library(ISLR)
set.seed(123)   # set seed to ensure you always have same random numbers generated
lg_train_id = sample(seq_len(nrow(lg_full)),size = floor(0.70*nrow(lg_full)))  # Randomly identifies therows equal to sample size ( defined in previous instruction) from  all the rows of Smarket dataset and stores the row number in train_ind
lg_train = lg_full[lg_train_id,] #creates the training dataset (70%)
lg_test = lg_full[-lg_train_id,]  # creates the test dataset (30%)

hs_train_id = sample(seq_len(nrow(hs_full)),size = floor(0.70*nrow(hs_full)))  # Randomly identifies therows equal to sample size ( defined in previous instruction) from  all the rows of Smarket dataset and stores the row number in train_ind
hs_train = hs_full[hs_train_id,] #creates the training dataset (70%)
hs_test = hs_full[-hs_train_id,]  # creates the test dataset (30%)

# save training records (presence-only) for MaxEnt
write.csv((lg_train %>% filter(PRES == 1) %>% select(SPECIES_CD, LON_M, LAT_M)), 
          paste0(fish_wd, "Presence_Only/Subadult_Gray_Snapper_Train.csv"), row.names = F)

write.csv((hs_train %>% filter(PRES == 1) %>% select(SPECIES_CD, LON_M, LAT_M)), 
          paste0(fish_wd, "Presence_Only/Subadult_Bluestriped_Grunt_Train.csv"), row.names = F)

# save testing records (presence-only) for MaxEnt
write.csv((lg_test %>% filter(PRES == 1) %>% select(SPECIES_CD, LON_M, LAT_M)), 
          paste0(fish_wd, "Presence_Only/Subadult_Gray_Snapper_Test.csv"), row.names = F)

write.csv((hs_test %>% filter(PRES == 1) %>% select(SPECIES_CD, LON_M, LAT_M)), 
          paste0(fish_wd, "Presence_Only/Subadult_Bluestriped_Grunt_Test.csv"), row.names = F)


# keep what is needed for logistic regression in caret
lg_train = lg_train %>% select(PRES2, Habitat, Mangrove_Dist, Depth, Slope, 
                         Curvature, BPI_Fine, BPI_Broad, Rugosity,
                         Mean_Sum_Temp, Mean_Sum_Sal, Mean_Win_Temp, 
                         Mean_Win_Sal)
lg_train$PRES2 = relevel(lg_train$PRES2, ref = "ABSENCE") #  ensure that absence is the reference category

lg_test = lg_test %>% select(PRES2, Habitat, Mangrove_Dist, Depth, Slope, 
                               Curvature, BPI_Fine, BPI_Broad, Rugosity,
                               Mean_Sum_Temp, Mean_Sum_Sal, Mean_Win_Temp, 
                               Mean_Win_Sal)
lg_test$PRES2 = relevel(lg_test$PRES2, ref = "ABSENCE") #  ensure that absence is the reference category


hs_train = hs_train %>% select(PRES2, Habitat, Mangrove_Dist, Depth, Slope, 
                         Curvature, BPI_Fine, BPI_Broad, Rugosity,
                         Mean_Sum_Temp, Mean_Sum_Sal, Mean_Win_Temp, 
                         Mean_Win_Sal)
hs_train$PRES2 = relevel(hs_train$PRES2, ref = "ABSENCE") #  ensure that absence is the reference category

hs_test = hs_test %>% select(PRES2, Habitat, Mangrove_Dist, Depth, Slope, 
                               Curvature, BPI_Fine, BPI_Broad, Rugosity,
                               Mean_Sum_Temp, Mean_Sum_Sal, Mean_Win_Temp, 
                               Mean_Win_Sal)
hs_test$PRES2 = relevel(hs_test$PRES2, ref = "ABSENCE") #  ensure that absence is the reference category

#### GRAY SNAPPER ####
#### Standard Logistic Regression ####
trainControl = trainControl(method = "cv", number = 10,
                            summaryFunction = twoClassSummary, 
                            classProbs = TRUE, savePredictions = TRUE)
lg_glm = caret::train(PRES2 ~ ., lg_train, method = "glm", family = "binomial",
                      trControl = trainControl, metric = "ROC", maximize = TRUE,
                      preProcess = c("center", "scale"), na.action = na.pass)

# ROC, sensitivity, and specificity (RECALL caret uses absence as a reference class,
# so sensitivity is for absence & specificity is for presence)
lg_glm 

# coefficients, deviance, and AIC, where coefficients characterize the relationship 
# between the predictors and species presence on a log-odds scale
summary(lg_glm) 

exp(coef(lg_glm$finalModel)) # or convert coefficients to odds ratio instead

# plot AUC-ROC Curve for training CV
lg_glm_eval = evalm(lg_glm, percent = 95, gnames = "Lujanus griseus", 
                    rlinethick = 0.75, fsize = 8, showplots = FALSE)
lg_glm_eval$roc

# predict probability of presence / absence across testing area validation data
lg_glm_probs = predict(lg_glm, newdata = lg_test, type = "prob") # get rid of last empty row

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
lg_glm_pred = predict(lg_glm, lg_test)
lg_glm_CM = confusionMatrix(lg_glm_pred, lg_test$PRES2,
                            positive = "PRESENCE") # ignore the empty final row that stores only habitat data

# IMPORTANT: remember, sensitivity refers to the reference class (absence) and specificity refers to the
# predicted class (presence). Swap the values if you want presence to be the positive class
lg_glm_CM 

#### Lasso Regularization ####
# LASSO = Least Absolute Shrinkage and Selection Operator, parsimonious model that
# performs L1 regularization (alpha = 1) by minimizing the absolute magnitude of 
# the regression coefficients . L1 is also used by MaxEnt. Coefficients that are 
# responsible for large variance are shrunk to ZERO.

# convert training data to matrix format for glmnet
lg_train_mat = model.matrix(PRES2 ~ ., lg_train)
# convert class to numerical variable
lg_train_y = ifelse(lg_train$PRES2 == "PRESENCE", 1, 0)

# convert testing data to matrix format
lg_test_mat = model.matrix(PRES2 ~ ., lg_test)
# convert class to numerical variable
lg_test_y = ifelse(lg_test$PRES2 == "PRESENCE", 1, 0)

# use 10-fold cross validation to find appropriate lambda
lg_lasso_cv = cv.glmnet(lg_train_mat, lg_train_y, alpha = 1, nfolds = 10,
                        family = "binomial", standardize = T, 
                        type.measure = "deviance")
plot(lg_lasso_cv)

# fitted coefficients, using minimum lambda
coef(lg_lasso_cv, s = "lambda.min")
# penalty term using minimum lambda
sum(coef(lg_lasso_cv, s = "lambda.min")[-1] ^ 2)

# run lasso regression in caret using minimum lambda value from cross validation
lg_lasso_grid = expand.grid(alpha = 1, lambda = lg_lasso_cv$lambda.min)
lg_lasso = train(PRES2 ~ ., lg_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = lg_lasso_grid, 
                 metric = "ROC", maximize = TRUE,  preProcess = c("center", "scale"))
lg_lasso # ROC, sensitivity, specificity
max(lg_lasso[["results"]]$ROC) # max AUC-ROC
coef(lg_lasso$finalModel, lg_lasso$bestTune$lambda) # coefficients log-odds
exp(coef(lg_lasso$finalModel, lg_lasso$bestTune$lambda)) # odds ratio


# plot AUC-ROC Curve from training CV
lg_lasso_eval = evalm(lg_lasso, percent = 95, gnames = "Lasso Regression", 
                      rlinethick = 0.75, fsize = 8, showplots = FALSE)
plot(lg_lasso_eval$roc)
#auc1 = plot(lg_lasso_eval$roc) + ggtitle("Lutjanus griseus") + 
 # theme(plot.title = element_text(size = 14, face = "italic")) + 
  #theme(legend.position = c(0.66, 0.14))
#png(filename = paste0(temp_plots, "Gray_Snapper_Lasso_AUC.png"), 
 #   width = 3.15, height = 3.15, units = "in", res = 300)
#auc1
#dev.off()

# predict probability of presence / absence across testing area validation data
lg_lasso_probs = predict(lg_lasso, newdata = lg_test, type = "prob")

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
lg_lasso_pred = predict(lg_lasso, newdata = lg_test)
lg_lasso_CM = confusionMatrix(data = lg_lasso_pred, reference = lg_test$PRES,
                              positive = "PRESENCE")
lg_lasso_CM 


#### Ridge Regularization ####
# L2 regularization (alpha = 0 ) that performs shrinkage by minimizing the sum of
# the squared coefficients. Unlike lasso, the coefficients in ridge regression can
# only asymptotically approach a value of zero.

# use 10-fold cross validation to find appropriate lambda
lg_ridge_cv = cv.glmnet(lg_train_mat, lg_train_y, alpha = 0, nfolds = 10,
                        family = "binomial", standardize = T, 
                        type.measure = "deviance")
plot(lg_ridge_cv)

# fitted coefficients, using minimum lambda
coef(lg_ridge_cv, s = "lambda.min")
# penalty term using minimum lambda
sum(coef(lg_ridge_cv, s = "lambda.min")[-1] ^ 2)

# run lasso regression in caret using minimum lambda value from cross validation
lg_ridge_grid = expand.grid(alpha = 0, lambda = lg_ridge_cv$lambda.min)
lg_ridge = train(PRES2 ~ ., lg_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = lg_ridge_grid, 
                 metric = "ROC", maximize = TRUE, preProcess = c("center", "scale"))
lg_ridge # ROC, sensitivity, specificity
max(lg_ridge[["results"]]$ROC) # max AUC-ROC
coef(lg_ridge$finalModel, lg_ridge$bestTune$lambda) # coefficients log-odds
exp(coef(lg_ridge$finalModel, lg_ridge$bestTune$lambda)) # odds ratio

# AUC-ROC Curve from training CV
lg_ridge_eval = evalm(lg_ridge, percent = 95, gnames = "Ridge Regression", 
                      rlinethick = 0.75, fsize = 8, showplots = FALSE)
plot(lg_ridge_eval$roc)
#auc2 = plot(lg_ridge_eval$roc) + ggtitle("") + 
 # theme(plot.title = element_text(size = 14, face = "italic")) + 
  #theme(legend.position = c(0.66, 0.14))
#png(filename = paste0(temp_plots, "Gray_Snapper_Ridge_AUC.png"), 
 #   width = 3.15, height = 3.15, units = "in", res = 300)
#auc2
#dev.off()

# predict probability of presence / absence across testing area validation data
lg_ridge_probs = predict(lg_ridge, newdata = lg_test, type = "prob")

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
lg_ridge_pred = predict(lg_ridge, newdata = lg_test)
lg_ridge_CM = confusionMatrix(data = lg_ridge_pred, reference = lg_test$PRES,
                              positive = "PRESENCE")
lg_ridge_CM 


#### BLUESTRIPED GRUNT ####
#### Standard Logistic Regression ####
hs_glm = caret::train(PRES2 ~ ., hs_train, method = "glm", family = "binomial",
                      trControl = trainControl, metric = "ROC", maximize = TRUE,
                      preProcess = c("center", "scale"))
hs_glm # ROC, sensitivity, specificity
summary(hs_glm) # coefficients, deviance, AIC

# plot AUC-ROC Curve from training CV
hs_glm_eval = evalm(hs_glm, percent = 95, gnames = "Haemulon sciurus", 
                    rlinethick = 0.75, fsize = 8, showplots = FALSE)
hs_glm_eval$roc

# predict probability of presence / absence across testing area validation data
hs_glm_probs = predict(hs_glm, newdata = hs_test, type = "prob")

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
hs_glm_pred = predict(hs_glm, hs_test)
hs_glm_CM = confusionMatrix(hs_glm_pred, hs_test$PRES,
                            positive = "PRESENCE")
hs_glm_CM 


#### Lasso Regularization ####
# convert training data to matrix format for glmnet
hs_train_mat = model.matrix(PRES2 ~ ., hs_train)
#convert class to numerical variable
hs_train_y = ifelse(hs_train$PRES2 == "PRESENCE", 1, 0)

# convert testing data to matrix format
hs_test_mat = model.matrix(PRES2 ~ ., hs_test)
#convert class to numerical variable
hs_test_y = ifelse(hs_test$PRES2 == "PRESENCE", 1, 0)

# use 10-fold cross validation to find appropriate lambda
hs_lasso_cv = cv.glmnet(hs_train_mat, hs_train_y, alpha = 1, nfolds = 10,
                        family = "binomial", standardize = T, 
                        type.measure = "deviance")
plot(hs_lasso_cv)

# fitted coefficients, using minimum lambda
coef(hs_lasso_cv, s = "lambda.min")
# penalty term using minimum lambda
sum(coef(hs_lasso_cv, s = "lambda.min")[-1] ^ 2)

# run logistic regression in caret using minimum lambda value from cross validation
hs_lasso_grid = expand.grid(alpha = 1, lambda = hs_lasso_cv$lambda.min)
hs_lasso = train(PRES2 ~ ., hs_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = hs_lasso_grid, 
                 metric = "ROC", maximize = TRUE, preProcess = c("center", "scale"))
hs_lasso # ROC, sensitivity, specificity
max(hs_lasso[["results"]]$ROC) # max AUC-ROC
coef(hs_lasso$finalModel, hs_lasso$bestTune$lambda) # coefficients log-odds
exp(coef(hs_lasso$finalModel, hs_lasso$bestTune$lambda)) # odds-ratio

# plot AUC-ROC Curve from training CV
hs_lasso_eval = evalm(hs_lasso, percent = 95, gnames = "Lasso Regression", 
                      rlinethick = 0.75, fsize = 8, showplots = FALSE)
plot(hs_lasso_eval$roc)
#auc3 = plot(hs_lasso_eval$roc) + ggtitle("Haemulon sciurus") + 
  #theme(plot.title = element_text(size = 14, face = "italic")) + 
  #theme(legend.position = c(0.66, 0.14))
#png(filename = paste0(temp_plots, "Bluestriped_Grunt_Lasso_AUC.png"), 
#    width = 3.15, height = 3.15, units = "in", res = 300)
#auc3
#dev.off()

# predict probability of presence / absence across testing area validation data
hs_lasso_probs = predict(hs_lasso, newdata = hs_test, type = "prob")

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
hs_lasso_pred = predict(hs_lasso, newdata = hs_test)
hs_lasso_CM = confusionMatrix(data = hs_lasso_pred, reference = hs_test$PRES,
                              positive = "PRESENCE")
hs_lasso_CM


#### Ridge Regularization ####
# use 10-fold cross validation to find appropriate lambda
hs_ridge_cv = cv.glmnet(hs_train_mat, hs_train_y, alpha = 0, nfolds = 10,
                        family = "binomial", standardize = T, 
                        type.measure = "deviance")
plot(hs_ridge_cv)

# fitted coefficients, using minimum lambda
coef(hs_ridge_cv, s = "lambda.min")
# penalty term using minimum lambda
sum(coef(hs_ridge_cv, s = "lambda.min")[-1] ^ 2)

# run ridge regression in caret using minimum lambda value from cross validation
hs_ridge_grid = expand.grid(alpha = 0, lambda = hs_ridge_cv$lambda.min)
hs_ridge = train(PRES2 ~ ., hs_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = hs_ridge_grid, 
                 metric = "ROC", maximize = TRUE, preProcess = c("center", "scale"))
hs_ridge # ROC, sensitivity, specificity
max(hs_ridge[["results"]]$ROC) # max AUC-ROC
coef(hs_ridge$finalModel, hs_ridge$bestTune$lambda) # coefficients log-odds
exp(coef(hs_ridge$finalModel, hs_ridge$bestTune$lambda)) # odds ratio

# plot AUC-ROC Curve from training CV
hs_ridge_eval = evalm(hs_lasso, percent = 95, gnames = "Ridge Regression", 
                      rlinethick = 0.75, fsize = 8, showplots = FALSE)
plot(hs_ridge_eval$roc)
#auc4 = plot(hs_ridge_eval$roc) + ggtitle(" ") + 
  #theme(plot.title = element_text(size = 14, face = "italic")) + 
  #theme(legend.position = c(0.66, 0.14))
#png(filename = paste0(temp_plots, "Bluestriped_Grunt_Ridge_AUC.png"), 
  #  width = 3.15, height = 3.15, units = "in", res = 300)
#auc4
#dev.off()

# predict probability of presence / absence across testing area validation data
hs_ridge_probs = predict(hs_ridge, newdata = hs_test, type = "prob")

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
hs_ridge_pred = predict(hs_ridge, newdata = hs_test)
hs_ridge_CM = confusionMatrix(data = hs_ridge_pred, reference = hs_test$PRES,
                              positive = "PRESENCE")
hs_ridge_CM

#### COMPARE ACCURACIES ####
# gray snapper
lg_glm_CM # 77.71 %
lg_lasso_CM # 76.39 %
lg_ridge_CM # 77.09 %

# bluestriped grunt
hs_glm_CM # 74.61 %
hs_lasso_CM # 74.61 %
hs_ridge_CM # 73.37 %

#### COMPARE AUC-ROC ####
# gray snapper
max(lg_glm[["results"]]$ROC) # ~0.7326
max(lg_lasso[["results"]]$ROC) #  ~ 0.7347
max(lg_ridge[["results"]]$ROC) #  ~ 0.7413

# bluestriped grunt
max(hs_glm[["results"]]$ROC) # ~0.7586
max(hs_lasso[["results"]]$ROC) #  ~ 0.7605
max(hs_ridge[["results"]]$ROC) #  ~ 0.7557

#### REGRESSION COEFFICIENTS ####
# look at one of the dgC matrices to find out the order of predictors
coef(lg_lasso$finalModel, lg_lasso$bestTune$lambda) 
coef_names = as.data.frame(c("Intercept", "Scattered Coral/Rock", "Continuous Seagrass",
                             "Discontinuous Seagrass", "Unconsolidated Sediment",
                             "Aggregate Reef", "Pavement", "Reef Rubble", "Mangrove",
                             "Mangrove Distance", "Depth", "Slope", "Curvature",
                             "BPI Fine", "BPI Broad", "Rugosity", "Summer Temperature",
                             "Summer Salinity", "Winter Temperature", "Winter Salinity"))

lg_lasso_coef = as.data.frame(as.vector(coef(lg_lasso$finalModel, lg_lasso$bestTune$lambda)))
lg_ridge_coef = as.data.frame(as.vector(coef(lg_ridge$finalModel, lg_ridge$bestTune$lambda)))
hs_lasso_coef = as.data.frame(as.vector(coef(hs_lasso$finalModel, hs_lasso$bestTune$lambda)))
hs_ridge_coef = as.data.frame(as.vector(coef(hs_ridge$finalModel, hs_ridge$bestTune$lambda)))

# merge and add species info
lg_lasso_coef = cbind(coef_names, lg_lasso_coef)
colnames(lg_lasso_coef) = c("Variable", "Coefficient")
lg_lasso_coef$Species = rep("Lutjanus griseus", nrow(lg_lasso_coef))
lg_lasso_coef$Life_Stage = rep("Subadult", nrow(lg_lasso_coef))

lg_ridge_coef = cbind(coef_names, lg_ridge_coef)
colnames(lg_ridge_coef) = c("Variable", "Coefficient")
lg_ridge_coef$Species = rep("Lutjanus griseus", nrow(lg_ridge_coef))
lg_ridge_coef$Life_Stage = rep("Subadult", nrow(lg_ridge_coef))


hs_lasso_coef = cbind(coef_names, hs_lasso_coef)
colnames(hs_lasso_coef) = c("Variable", "Coefficient")
hs_lasso_coef$Species = rep("Haemulon sciurus", nrow(hs_lasso_coef))
hs_lasso_coef$Life_Stage = rep("Subadult", nrow(hs_lasso_coef))


hs_ridge_coef = cbind(coef_names, hs_ridge_coef)
colnames(hs_ridge_coef) = c("Variable", "Coefficient")
hs_ridge_coef$Species = rep("Haemulon sciurus", nrow(hs_ridge_coef))
hs_ridge_coef$Life_Stage = rep("Subadult", nrow(hs_ridge_coef))

# combine lasso coefficients for both species, repeat for ridge and save
lasso_coef = rbind(lg_lasso_coef, hs_lasso_coef)
ridge_coef = rbind(lg_ridge_coef, hs_ridge_coef)

# the original depth data are negative, meaning larger negative numbers (deeper water)
# are technically "smaller" than less negative (shallower) values. This is a bit confusing 
# when looking at coefficients, because technically a positive depth coefficient means that
# habitat suitability increases as depth values get larger (shallower). Instead, flip
# the sign of the coefficient so that they represent the change in suitability at increasing depths
lasso_coef[11, 2] = as.numeric(-1*lasso_coef[11, 2])
lasso_coef[31, 2] = as.numeric(-1*lasso_coef[31, 2])

ridge_coef[11, 2] = as.numeric(-1*ridge_coef[11, 2])
ridge_coef[31, 2] = as.numeric(-1*ridge_coef[31, 2])

# finally add the odds ratios to each table and save
lasso_coef$Odds_Ratio = exp(lasso_coef$Coefficient)
ridge_coef$Odds_Ratio = exp(ridge_coef$Coefficient)
lasso_coef = lasso_coef %>% relocate(Odds_Ratio, .after = Coefficient)
ridge_coef = ridge_coef %>% relocate(Odds_Ratio, .after = Coefficient)

write.csv(lasso_coef, paste0(csv_wd, "Lasso_Coefficients.csv"), row.names = FALSE)
write.csv(ridge_coef, paste0(csv_wd, "Ridge_Coefficients.csv"), row.names = FALSE)

# plotting coefficients
my_pal = pnw_palette("Bay",8)
my_pal

lasso_coef_plot = ggplot(data = lasso_coef, aes(x = Coefficient, y = Variable, fill = Species)) +
  geom_vline(xintercept = 0, color = "gray") + geom_point(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  labs(y = "", x = "Lasso regression coefficient") + xlim(-2,2)  +
  scale_y_discrete(limits = c("Intercept", "Winter Salinity", "Winter Temperature", "Summer Salinity", 
                              "Summer Temperature", "Rugosity", "BPI Broad", "BPI Fine",
                              "Curvature", "Slope", "Depth", "Mangrove Distance", "Mangrove",
                              "Reef Rubble", "Pavement", "Aggregate Reef", "Unconsolidated Sediment",
                              "Discontinuous Seagrass", "Continuous Seagrass", "Scattered Coral/Rock")) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_rect(color = "black"), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), legend.text = element_text(face = "italic", size = 8),
                     legend.margin = margin(0,0,0,0), plot.margin = margin(c(0,0,0,0)),
                     axis.text = element_text(size = 8), axis.title.x = element_text(size = 10),
                     legend.box.margin=margin(t = 0, r = 0, b = -10, l = 0))
lasso_coef_plot

ridge_coef_plot = ggplot(data = ridge_coef, aes(x = Coefficient, y = Variable, fill = Species)) +
  geom_vline(xintercept = 0, color = "gray") + geom_point(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  labs(y = "", x = "Ridge regression coefficient") + xlim(-2,2) +
  scale_y_discrete(limits = c("Intercept", "Winter Salinity", "Winter Temperature", "Summer Salinity", 
                              "Summer Temperature", "Rugosity", "BPI Broad", "BPI Fine",
                              "Curvature", "Slope", "Depth", "Mangrove Distance", "Mangrove",
                              "Reef Rubble", "Pavement", "Aggregate Reef", "Unconsolidated Sediment",
                              "Discontinuous Seagrass", "Continuous Seagrass", "Scattered Coral/Rock")) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_rect(color = "black"), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), legend.text = element_text(face = "italic", size = 8),
                     legend.margin = margin(0,0,0,0), plot.margin = margin(c(0,0,0,0)),
                     axis.text = element_text(size = 8), axis.title.x = element_text(size = 10),
                     legend.box.margin = margin(t = 0, r = 0, b = -10, l = 0))
ridge_coef_plot

ridge_coef_plot2 = ridge_coef_plot + theme(legend.position = "none",  plot.margin = margin(c(t = 8,0,b = 5,0)))

coef_grid = plot_grid(lasso_coef_plot, ridge_coef_plot2, nrow = 2)

ggsave(plot = lasso_coef_plot, filename = paste0(temp_plots, "Lasso_Coefficients.png"), height = 3.15, 
       width = 5, units = "in", dpi = 450)

ggsave(plot = ridge_coef_plot, filename = paste0(temp_plots, "Ridge_Coefficients.png"), height = 3.15, 
       width = 5, units = "in", dpi = 450)

ggsave(plot = coef_grid, filename = paste0(temp_plots, "Regression_Coefficients_Grid.png"), height = 5, 
       width = 5, units = "in", dpi = 450)

#### VARIABLE IMPORTANCE ####

lg_lasso_vimp = caret::varImp(lg_lasso, lambda = lg_lasso$bestTune$lambda)$importance
lg_lasso_vimp # check order of variables 
names = data.frame("Variable" = c("Scattered Coral/Rock", "Continuous Seagrass",
                                  "Discontinuous Seagrass", "Unconsolidated Sediment",
                                  "Aggregate Reef", "Pavement", "Reef Rubble", "Mangrove",
                                  "Mangrove Distance", "Depth", "Slope", "Curvature",
                                  "BPI Fine", "BPI Broad", "Rugosity", "Summer Temperature",
                                  "Summer Salinity", "Winter Temperature", "Winter Salinity"))
lg_lasso_vimp = cbind(lg_lasso_vimp, names)

lg_ridge_vimp = cbind(caret::varImp(lg_ridge, lambda = lg_ridge$bestTune$lambda)$importance, names)
hs_lasso_vimp = cbind(caret::varImp(hs_lasso, lambda = hs_lasso$bestTune$lambda)$importance, names)
hs_ridge_vimp = cbind(caret::varImp(hs_ridge, lambda = hs_ridge$bestTune$lambda)$importance, names)

lg_lasso_vimp$Species = rep("Lutjanus griseus", nrow(lg_lasso_vimp))
lg_ridge_vimp$Species = rep("Lutjanus griseus", nrow(lg_ridge_vimp))
hs_lasso_vimp$Species = rep("Haemulon sciurus", nrow(hs_lasso_vimp))
hs_ridge_vimp$Species = rep("Haemulon sciurus", nrow(hs_ridge_vimp))

lg_lasso_vimp$Life_Stage = rep("Subadult", nrow(lg_lasso_vimp))
lg_ridge_vimp$Life_Stage = rep("Subadult", nrow(lg_ridge_vimp))
hs_lasso_vimp$Life_Stage = rep("Subadult", nrow(hs_lasso_vimp))
hs_ridge_vimp$Life_Stage = rep("Subadult", nrow(hs_ridge_vimp))

lasso_vimp = rbind(lg_lasso_vimp, hs_lasso_vimp)
ridge_vimp = rbind(lg_ridge_vimp, hs_ridge_vimp)

# save variable importance info
write.csv(lasso_vimp, paste0(csv_wd, "Lasso_Variable_Importance.csv"))
write.csv(ridge_vimp, paste0(csv_wd, "Ridge_Variable_Importance.csv"))


lasso_vimp_plot = 
  ggplot(lasso_vimp, aes(x = Variable, y = Overall, fill = Species)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  xlab(" ") + ylab("Lasso Regression Variable Importance") +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  scale_x_discrete(limits = c("Winter Salinity", "Winter Temperature", "Summer Salinity", 
                              "Summer Temperature", "Rugosity", "BPI Broad", "BPI Fine",
                              "Curvature", "Slope", "Depth", "Mangrove Distance", "Mangrove",
                              "Reef Rubble", "Pavement", "Aggregate Reef", "Unconsolidated Sediment",
                              "Discontinuous Seagrass", "Continuous Seagrass", "Scattered Coral/Rock")) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(color = "black"),
                     axis.text.y = element_text(color = "black"),
                     legend.text = element_text(face = "italic"),
                     legend.box.margin = margin(t = 0, r = 0, b = -10, l = 0)) +
  geom_text(aes(label = round(Overall, digits = 2)), size = 3, 
            position = position_dodge(width = 1), hjust = -0.25, vjust = 0.45) + 
  coord_flip() 
lasso_vimp_plot

ridge_vimp_plot = 
  ggplot(ridge_vimp, aes(x = Variable, y = Overall, fill = Species)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  xlab(" ") + ylab("Ridge Regression Variable Importance") +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) + 
  scale_x_discrete(limits = c("Winter Salinity", "Winter Temperature", "Summer Salinity", 
                              "Summer Temperature", "Rugosity", "BPI Broad", "BPI Fine",
                              "Curvature", "Slope", "Depth", "Mangrove Distance", "Mangrove",
                              "Reef Rubble", "Pavement", "Aggregate Reef", "Unconsolidated Sediment",
                              "Discontinuous Seagrass", "Continuous Seagrass", "Scattered Coral/Rock")) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(color = "black"),
                     axis.text.y = element_text(color = "black"),
                     legend.text = element_text(face = "italic"),
                     legend.box.margin = margin(t = 0, r = 0, b = -10, l = 0)) +
  geom_text(aes(label = round(Overall, digits = 2)), size = 3, 
            position = position_dodge(width = 1), hjust = -0.25, vjust = 0.45) + 
  coord_flip()
ridge_vimp_plot

ggsave(plot = lasso_vimp_plot, filename = paste0(temp_plots, "Lasso_Variable_Importance.png"), height = 5, 
       width = 8.5, units = "in", dpi = 450)

ggsave(plot = ridge_vimp_plot, filename = paste0(temp_plots, "Ridge_Variable_Importance.png"), height = 5, 
       width = 8.5, units = "in", dpi = 450)

#### TOP 5 VARIABLES ####
# five most important predictors across species-model combinations
lasso_top5 = lasso_vimp %>% filter((Species == "Lutjanus griseus" & Variable %in% c("Depth", "Mangrove Distance", "Reef Rubble",
                                                                       "Pavement", "Summer Temperature")) |
                                     (Species == "Haemulon sciurus" & Variable %in% c("Winter Salinity", "Depth", "Pavement",
                                                                         "Summer Salinity", "Unconsolidated Sediment")))

ridge_top5 = ridge_vimp %>% filter((Species == "Lutjanus griseus" & Variable %in% c("Depth", "Mangrove Distance", "Pavement",
                                                                                    "Reef Rubble", "Winter Temperature")) |
                                     (Species == "Haemulon sciurus" & Variable %in% c("Winter Salinity", "Depth", "Pavement",
                                                                                      "Summer Salinity", "Unconsolidated Sediment")))

# plotting
lasso_top5_plot = ggplot((lasso_top5 %>% group_by(Species)), 
                         aes(x = interaction(Variable, Overall, Species), 
                             y = Overall, fill = Species)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge(), group = lasso_top5$Species) +
  xlab(" ") + ylab("Lasso regression variable importance") +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) + 
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(color = "black"),
                     axis.text.y = element_text(color = "black"),
                     legend.text = element_text(face = "italic"),
                     legend.box.margin = margin(t = 0, r = 0, b = -10, l = 0)) +
  geom_text(aes(label = round(Overall, digits = 2)), size = 3, 
            position = position_dodge(width = 1), hjust = -0.25, vjust = 0.45) + 
  coord_flip()
lasso_top5_plot
lasso_top5_plot = lasso_top5_plot + scale_x_discrete(labels = c("Depth.100.Lutjanus griseus" = "Depth",
                                                                "Mangrove Distance.57.1531077410818.Lutjanus griseus" = "Mangrove Distance",
                                                                "Reef Rubble.49.6410017869303.Lutjanus griseus" = "Reef Rubble",
                                                                "Pavement.47.2176003335797.Lutjanus griseus" = "Pavement",
                                                                "Summer Temperature.38.4700757087235.Lutjanus griseus" = "Summer Temperature",
                                                                "Winter Salinity.100.Haemulon sciurus" = "Winter Salinity",
                                                                "Depth.50.9633789015947.Haemulon sciurus" = "Depth",
                                                                "Pavement.30.1289151728483.Haemulon sciurus" = "Pavement",
                                                                "Summer Salinity.26.6816671874025.Haemulon sciurus" = "Summer Salinity",
                                                                "Unconsolidated Sediment.20.420931627249.Haemulon sciurus" = "Unconsolidated Sediment"))
lasso_top5_plot

ridge_top5_plot = ggplot((ridge_top5 %>% group_by(Species)), 
                         aes(x = interaction(Variable, Overall, Species), 
                             y = Overall, fill = Species)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge(), group = ridge_top5$Species) +
  xlab(" ") + ylab("Ridge regression variable importance") +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) + 
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(color = "black"),
                     axis.text.y = element_text(color = "black"),
                     legend.text = element_text(face = "italic"),
                     legend.box.margin = margin(t = 0, r = 0, b = -10, l = 0)) +
  geom_text(aes(label = round(Overall, digits = 2)), size = 3, 
            position = position_dodge(width = 1), hjust = -0.25, vjust = 0.45) + 
  coord_flip()
ridge_top5_plot
ridge_top5_plot = ridge_top5_plot + scale_x_discrete(labels = c("Depth.100.Lutjanus griseus" = "Depth",
                                                                "Mangrove Distance.90.4765641574874.Lutjanus griseus" = "Mangrove Distance",
                                                                "Pavement.77.1229623162004.Lutjanus griseus" = "Pavement",
                                                                "Reef Rubble.72.767215856622.Lutjanus griseus" = "Reef Rubble",
                                                                "Winter Temperature.66.1729254744799.Lutjanus griseus" = "Winter Temperature",
                                                                "Winter Salinity.100.Haemulon sciurus" = "Winter Salinity",
                                                                "Depth.48.0209817702403.Haemulon sciurus" = "Depth",
                                                                "Pavement.42.6318832181052.Haemulon sciurus" = "Pavement",
                                                                "Summer Salinity.38.9856234349051.Haemulon sciurus" = "Summer Salinity",
                                                                "Unconsolidated Sediment.28.7635412734945.Haemulon sciurus" = "Unconsolidated Sediment"))
ridge_top5_plot


ggsave(plot = lasso_top5_plot, filename = paste0(temp_plots, "Lasso_Top5_Variables.png"), height = 5, 
       width = 8.5, units = "in", dpi = 450)

ggsave(plot = ridge_top5_plot, filename = paste0(temp_plots, "Ridge_Top5_Variables.png"), height = 5, 
       width = 8.5, units = "in", dpi = 450)


#### CONTINUOUS (RASTER) SUITABILITY ####
# save folder path
HSMs = "Z:/Courtney/Stuart_MSc_Ch1/HSMs/Logistic_Regression/"

# make predictions across training area 
# FYI: very time consuming step, each model takes up to 10 hours to build
library(snow)

# make sure raster stack layer names match lasso model variables
env
lg_lasso$coefnames
names(env) <- c("Habitat", "Mangrove_Dist", "Depth", "Slope", "Curvature", "BPI_Fine", "BPI_Broad", "Rugosity",
                "Mean_Sum_Temp", "Mean_Sum_Sal", "Mean_Win_Temp", "Mean_Win_Sal")

# IMPORTANT: specify (1 - predict) because the raster package will use the first class as its
# output, and that is the absence reference class in this case. We want probability of presence!
# you'll notice this small caveat if you open the lg_lasso_probs df - absence is first.
beginCluster(n = 10)
lg_lasso_HSM = writeRaster((1-(raster::predict(env, lg_lasso, type = "prob"))),
                                 filename = paste0(HSMs, 
                                 "Subadult_Gray_Snapper_Lasso.asc"),
                                 format = "ascii", overwrite = TRUE)
endCluster()

# repeat for lg ridge training
beginCluster(n = 10)
lg_train_ridge_HSM = writeRaster((1-(raster::predict(env_train, lg_ridge, type = "prob"))),
                                 filename = paste0(HSMs, 
                                 "Logistic_Regression/Subadult_Gray_Snapper_Ridge_Train.asc"),
                                 format = "ascii", overwrite = TRUE)
endCluster()

# repeat for hs lasso training
beginCluster(n = 10)
hs_train_lasso_HSM = writeRaster((1-(raster::predict(env_train, hs_lasso, type = "prob"))),
                                 filename = paste0(HSMs, 
                                 "Logistic_Regression/Subadult_Bluestriped_Grunt_Lasso_Train.asc"),
                                 format = "ascii", overwrite = TRUE)
endCluster()

# repeat for hs ridge training
beginCluster(n = 10)
hs_train_ridge_HSM = writeRaster((1-(raster::predict(env_train, hs_ridge, type = "prob"))),
                                 filename = paste0(HSMs, 
                                 "Logistic_Regression/Subadult_Bluestriped_Grunt_Ridge_Train.asc"),
                                 format = "ascii", overwrite = TRUE)
endCluster()

#### Testing Area ####
# spatial predictors in testing area
habitat_test = raster(paste0(test_wd, "Environmental/Habitat.asc"))
mg_dist_test = raster(paste0(test_wd, "Environmental/Mangrove_Dist.asc"))
depth_test = raster(paste0(test_wd, "Environmental/Depth.asc"))
slope_test = raster(paste0(test_wd, "Environmental/Slope.asc"))
curvature_test = raster(paste0(test_wd, "Environmental/Curvature.asc"))
bpi_fine_test = raster(paste0(test_wd, "Environmental/BPI_Fine.asc"))
bpi_broad_test = raster(paste0(test_wd, "Environmental/BPI_Broad.asc"))
rugosity_test = raster(paste0(test_wd, "Environmental/Rugosity.asc"))
sum_temp_test = raster(paste0(test_wd, "Environmental/Mean_Sum_Temp.asc"))
sum_sal_test = raster(paste0(test_wd, "Environmental/Mean_Sum_Sal.asc"))
win_sal_test = raster(paste0(test_wd, "Environmental/Mean_Win_Sal.asc"))

# define crs
crs(habitat_test) = my_crs
crs(mg_dist_test) = my_crs
crs(depth_test) = my_crs
crs(slope_test) = my_crs
crs(curvature_test) = my_crs
crs(bpi_fine_test) = my_crs
crs(bpi_broad_test) = my_crs
crs(rugosity_test) = my_crs
crs(sum_temp_test) = my_crs
crs(sum_sal_test) = my_crs
crs(win_sal_test) = my_crs

# create raster stack for testing area
env_test = stack(x = c(habitat_test, mg_dist_test, depth_test,
                       slope_test, curvature_test, bpi_fine_test,
                       bpi_broad_test, rugosity_test, sum_temp_test, 
                       sum_sal_test, win_sal_test))

beginCluster(n = 10)
lg_test_lasso_HSM = writeRaster((1-(raster::predict(env_test, lg_lasso, type = "prob"))),
                                filename = paste0(HSMs, 
                                "Logistic_Regression/Subadult_Gray_Snapper_Lasso_Test.asc"),
                                format = "ascii", overwrite = TRUE)
endCluster()

# repeat for lg ridge testing
beginCluster(n = 10)
lg_test_ridge_HSM = writeRaster((1-(raster::predict(env_test, lg_ridge, type = "prob"))),
                                filename = paste0(HSMs, 
                                "Logistic_Regression/Subadult_Gray_Snapper_Ridge_Test.asc"),
                                format = "ascii", overwrite = TRUE)
endCluster()

# repeat for hs lasso testing
beginCluster(n = 10)
hs_test_lasso_HSM = writeRaster((1-(raster::predict(env_test, hs_lasso, type = "prob"))),
                                filename = paste0(HSMs, 
                                "Logistic_Regression/Subadult_Bluestriped_Grunt_Lasso_Test.asc"),
                                format = "ascii", overwrite = TRUE)
endCluster()

# repeat for hs ridge testing
beginCluster(n = 10)
hs_test_ridge_HSM = writeRaster((1-(raster::predict(env_test, hs_ridge, type = "prob"))),
                                filename = paste0(HSMs, 
                                "Logistic_Regression/Subadult_Bluestriped_Grunt_Ridge_Test.asc"),
                                format = "ascii", overwrite = TRUE)
endCluster()

#### MAX SSS THRESHOLD ####
# maximize the sum of the training sensitivity + specificity (max SSS)
lg_lasso_th = thresholder(lg_lasso, statistics = c("Sensitivity", "Specificity"),
                          threshold = seq(0, 1, by = 0.01))
lg_lasso_th$SSS = lg_lasso_th$Sensitivity + lg_lasso_th$Specificity
# max SSS occurs at a threshold of 70 (remember this is referring to absence),
# so anything at or above 30% is considered suitable

lg_ridge_th = thresholder(lg_ridge, statistics = c("Sensitivity", "Specificity"),
                          threshold = seq(0, 1, by = 0.01))
lg_ridge_th$SSS = lg_ridge_th$Sensitivity + lg_ridge_th$Specificity
# max SSS occurs at a threshold of 70 (remember this is referring to absence),
# so anything at or above 30% is considered suitable

hs_lasso_th = thresholder(hs_lasso, statistics = c("Sensitivity", "Specificity"),
                          threshold = seq(0, 1, by = 0.01))
hs_lasso_th$SSS = hs_lasso_th$Sensitivity + hs_lasso_th$Specificity
# max SSS occurs at a threshold of 70 (remember this is referring to absence),
# so anything at or above 30% is considered suitable

hs_ridge_th = thresholder(hs_ridge, statistics = c("Sensitivity", "Specificity"),
                          threshold = seq(0, 1, by = 0.01))
hs_ridge_th$SSS = hs_ridge_th$Sensitivity + hs_ridge_th$Specificity
# max SSS occurs at a threshold of 73 (remember this is referring to absence),
# so anything at or above 27% is considered suitable

# predict binary response (presence/absence) across testing area and calculate map accuracy
# max SSS threshold for suitability is 0.30
lg_lasso_pred30 = as.factor(ifelse(lg_lasso_probs$PRESENCE >= 0.3, "PRESENCE", "ABSENCE"))
lg_lasso_CM30 = confusionMatrix(data = lg_lasso_pred30, reference = lg_test[-663,]$PRES,
                              positive = "PRESENCE")
lg_lasso_CM30 


lg_ridge_pred30 = as.factor(ifelse(lg_ridge_probs$PRESENCE >= 0.3, "PRESENCE", "ABSENCE"))
lg_ridge_CM30 = confusionMatrix(data = lg_ridge_pred30, reference = lg_test[-663,]$PRES,
                                positive = "PRESENCE")
lg_ridge_CM30 


hs_lasso_pred30 = as.factor(ifelse(hs_lasso_probs$PRESENCE >= 0.3, "PRESENCE", "ABSENCE"))
hs_lasso_CM30 = confusionMatrix(data = hs_lasso_pred30, reference = hs_test[-663,]$PRES,
                                positive = "PRESENCE")
hs_lasso_CM30 


hs_ridge_pred30 = as.factor(ifelse(hs_ridge_probs$PRESENCE >= 0.3, "PRESENCE", "ABSENCE"))
hs_ridge_CM30 = confusionMatrix(data = hs_ridge_pred30, reference = hs_test[-663,]$PRES,
                                positive = "PRESENCE")
hs_ridge_CM30 
