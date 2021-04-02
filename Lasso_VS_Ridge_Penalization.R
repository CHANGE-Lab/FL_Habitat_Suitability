# This script is part of Courtney Stuart's first MSc Chapter in the Lab of 
# Dr. Stephanie Green at the University of Alberta. Data are specific to southern
# Florida, including the georeferenced occurrence records of two reef fish species
# (gray snapper Lutjanus griseus (lg) and bluestriped grunt Haemulon sciurus (hs))
# and rasters of various environmental covariates. Three logistic regression models
# were fitted for each species to create habitat suitability predictions: standard 
# logistic regression, lasso-regularized regression, and ridge-regularized regression. 


#### SET-UP ####
# data directories
setwd("Z:/Courtney/Stuart_MSc_Ch1/") # main project folder
temp_wd = "Z:/Courtney/Stuart_MSc_Ch1/Temp/" # temporary files
csv_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Data/" # for writing tabular data
train_wd = "Z:/Courtney/Stuart_MSc_Ch1/Modeling_Data/Training/" # for training area data
test_wd = "Z:/Courtney/Stuart_MSc_Ch1/Modeling_Data/Testing/" # for testing area data
plots_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Figures/" # for figures
temp_plots = "Z:/Courtney/Stuart_MSc_Ch1/Plots/" # temporary plots for Courtney

# libraries
library(easypackages)
libraries("glmnet", "tidyr", "caret", "MLeval", "dplyr", "raster")


# presence-absence datasets from model training area and model testing area
# keep only the reduced covariate set (following multicollinearity assessment)
lg_train = read.csv(paste0(csv_wd, "Gray_Snapper_Training_Data.csv")) %>%
  dplyr::select(PRES, Habitat, Mangrove_Dist, Depth, Slope, Curvature, Rugosity,
         BPI_Broad, BPI_Fine, Sum_Temp, Sum_Sal, Win_Sal) %>%
  dplyr::mutate(PRES = as.factor(PRES),
         Habitat = as.factor(Habitat))
levels(lg_train$PRES) # reverse levels so "PRESENCE" is the "positive" class 
lg_train$PRES = factor(lg_train$PRES, levels=rev(levels(lg_train$PRES)))

lg_test = read.csv(paste0(csv_wd, "Gray_Snapper_Testing_Data.csv")) %>%
  dplyr::select(PRES, Habitat, Mangrove_Dist, Depth, Slope, Curvature, Rugosity,
         BPI_Broad, BPI_Fine, Sum_Temp, Sum_Sal, Win_Sal) %>%
  dplyr::mutate(PRES = as.factor(PRES),
         Habitat = as.factor(Habitat))
levels(lg_test$PRES) # reverse levels so "PRESENCE" is the "positive" class 
lg_test$PRES = factor(lg_test$PRES, levels=rev(levels(lg_test$PRES)))


hs_train = read.csv(paste0(csv_wd, "Bluestriped_Grunt_Training_Data.csv"))%>%
  dplyr::select(PRES, Habitat, Mangrove_Dist, Depth, Slope, Curvature, Rugosity,
         BPI_Broad, BPI_Fine, Sum_Temp, Sum_Sal, Win_Sal) %>%
  dplyr::mutate(PRES = as.factor(PRES),
         Habitat = as.factor(Habitat))
levels(hs_train$PRES) # reverse levels so "PRESENCE" is the "positive" class 
hs_train$PRES = factor(hs_train$PRES, levels=rev(levels(hs_train$PRES)))

hs_test = read.csv(paste0(csv_wd, "Bluestriped_Grunt_Testing_Data.csv"))%>%
  dplyr::select(PRES, Habitat, Mangrove_Dist, Depth, Slope, Curvature, Rugosity,
         BPI_Broad, BPI_Fine, Sum_Temp, Sum_Sal, Win_Sal) %>%
  dplyr::mutate(PRES = as.factor(PRES),
         Habitat = as.factor(Habitat))
levels(hs_test$PRES) # reverse levels so "PRESENCE" is the "positive" class 
hs_test$PRES = factor(hs_test$PRES, levels=rev(levels(hs_test$PRES)))


# add an empty row with mangrove habitat (ID# 11) to the test data so an equal 
# number of habitat levels are present in the training and testing data
add_mg = data.frame(NA, as.factor(11), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
names(add_mg) = c("PRES", "Habitat", "Mangrove_Dist", "Depth", "Slope", "Curvature",
                 "Rugosity", "BPI_Broad", "BPI_Fine", "Sum_Temp", "Sum_Sal", "Win_Sal")
lg_test = rbind(lg_test, add_mg) # adds an extra row that is empty other than 11 for habitat
hs_test = rbind(hs_test, add_mg)

# before getting started
set.seed(42)

#### GRAY SNAPPER ####
#### Standard Logistic Regression ####
# build standard logistic regression model (not penalized)
trainControl = trainControl(method = "cv", number = 10,
                          summaryFunction = twoClassSummary, 
                          classProbs = TRUE, savePredictions = TRUE)
lg_glm = caret::train(PRES ~ ., lg_train, method = "glm", family = "binomial",
                      trControl = trainControl, metric = "ROC", maximize = TRUE)
lg_glm # ROC, sensitivity, specificity
summary(lg_glm) # coefficients, deviance, AIC

# plot AUC-ROC Curve for training CV
lg_glm_eval = evalm(lg_glm, percent = 95, gnames = "Lujanus griseus", 
      rlinethick = 0.75, fsize = 8, showplots = FALSE)
lg_glm_eval$roc

# binary presence/absence predictions for testing area
lg_glm_pred = predict(lg_glm, lg_test)
lg_glm_CM = confusionMatrix(lg_glm_pred, lg_test[-663,]$PRES, 
                            positive = "PRESENCE") # ignore the empty final row that stores only habitat data
lg_glm_CM # accuracy ~ 84% (83.84)

lg_glm_imp = varImp(lg_glm, useModel = TRUE, nonpara = TRUE, scale = TRUE)
View(lg_glm_imp$importance)

# plot AUC-ROC Curve for testing area
lg_glm_pred2 = predict(lg_glm, newdata = lg_test[-663,], type = "prob")
lg_glm_test_AUC = evalm(data.frame(lg_glm_pred2, lg_test[-663,]$PRES), 
                        showplots = FALSE)
lg_glm_test_AUC$roc

#### Lasso Regularization ####
# LASSO = Least Absolute Shrinkage and Selection Operator, parsimonious model that
# performs L1 regularization (alpha = 1) by minimizing the absolute magnitude of 
# the regression coefficients . L1 is also used by MaxEnt. Coefficients that are 
# responsible for large variance are shrunk to ZERO.

# convert training data to matrix format for glmnet
lg_train_mat = model.matrix(PRES ~ ., lg_train)
#convert class to numerical variable
lg_train_y = ifelse(lg_train$PRES == "PRESENCE", 1, 0)

# convert testing data to matrix format
lg_test_mat = model.matrix(PRES ~ ., lg_test)
#convert class to numerical variable
lg_test_y = ifelse(lg_test$PRES == "PRESENCE", 1, 0)

# fit lasso using 10-fold cross validation
lg_lasso_cv = cv.glmnet(lg_train_mat, lg_train_y, alpha = 1, nfolds = 10,
                        family = "binomial", standardize = T, 
                        type.measure = "deviance")
plot(lg_lasso_cv)

# fitted coefficients, using minimum lambda
coef(lg_lasso_cv, s = "lambda.min")
# penalty term using minimum lambda
sum(coef(lg_lasso_cv, s = "lambda.min")[-1] ^ 2)

# run penalized model in caret using minimum lambda value from cross validation
lg_lasso_grid = expand.grid(alpha = 1, lambda = lg_lasso_cv$lambda.min)
lg_lasso = train(PRES ~ ., lg_train, method = "glmnet", family = "binomial",
              trControl = trainControl, tuneGrid = lg_lasso_grid, 
              metric = "ROC", maximize = TRUE)
lg_lasso # ROC, sensitivity, specificity
max(lg_lasso[["results"]]$ROC) # max AUC-ROC
coef(lg_lasso$finalModel, lg_lasso$bestTune$lambda) # coefficients 
# Slope, curvature, BPI Broad, Sum Sal, and a few habitat classes driven to 0 

# predictions and confusion matrix
lg_lasso_pred = predict(lg_lasso, newdata = lg_test[-663,])
lg_lasso_CM = confusionMatrix(data = lg_lasso_pred, reference = lg_test[-663,]$PRES, 
                              positive = "PRESENCE")
lg_lasso_CM # accuracy ~ 84% (84.14)

lg_lasso_imp = varImp(lg_lasso, useModel = TRUE, nonpara = TRUE, scale = TRUE)
View(lg_lasso_imp$importance)

# plot AUC-ROC Curve for training CV
lg_lasso_eval = evalm(lg_lasso, percent = 95, gnames = "Lasso Regression", 
                    rlinethick = 0.75, fsize = 8, showplots = FALSE)
auc1 = plot(lg_lasso_eval$roc) + ggtitle("Lutjanus griseus") + 
  theme(plot.title = element_text(size = 14, face = "italic")) + 
  theme(legend.position = c(0.66, 0.14))
png(filename = paste0(temp_plots, "Gray_Snapper_Lasso_AUC.png"), 
    width = 3.15, height = 3.15, units = "in", res = 300)
auc1
dev.off()

# plot AUC-ROC Curve for testing area
lg_lasso_pred2 = predict(lg_lasso, newdata = lg_test[-663,], type = "prob")
lg_lasso_test_AUC = evalm(data.frame(lg_lasso_pred2, lg_test[-663,]$PRES),
                          showplots = FALSE)
lg_lasso_test_AUC$roc

#### Ridge Regularization ####
# L2 regularization (aplha = 0 ) that performs shrinkage by minimizing the sum of
# the squared coefficients. Unlike lasso, the coefficients in ridge regression can
# only asymptotically approach a value of zero.

# fit ridge using 10-fold cross validation
lg_ridge_cv = cv.glmnet(lg_train_mat, lg_train_y, alpha = 0, nfolds = 10,
                        family = "binomial", standardize = T, 
                        type.measure = "deviance")
plot(lg_ridge_cv)

# fitted coefficients, using minimum lambda
coef(lg_ridge_cv, s = "lambda.min")
# penalty term using minimum lambda
sum(coef(lg_ridge_cv, s = "lambda.min")[-1] ^ 2)

# run penalized model in caret using minimum lambda value from cross validation
lg_ridge_grid = expand.grid(alpha = 0, lambda = lg_ridge_cv$lambda.min)
lg_ridge = train(PRES ~ ., lg_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = lg_ridge_grid, 
                 metric = "ROC", maximize = TRUE)
lg_ridge # ROC, sensitivity, specificity
max(lg_ridge[["results"]]$ROC) # max AUC-ROC
coef(lg_ridge$finalModel, lg_ridge$bestTune$lambda) # coefficients 
# Mangrove distance, BPI Broad, and BPI Fine pretty close to 0

# predictions and confusion matrix
lg_ridge_pred = predict(lg_ridge, newdata = lg_test[-663,])
lg_ridge_CM = confusionMatrix(data = lg_ridge_pred, reference = lg_test[-663,]$PRES, 
                              positive = "PRESENCE")
lg_ridge_CM # accuracy ~ 85% (84.59)

lg_ridge_imp = varImp(lg_ridge, useModel = TRUE, nonpara = TRUE, scale = TRUE)
View(lg_ridge_imp$importance)

# AUC-ROC Curve for training CV
lg_ridge_eval = evalm(lg_ridge, percent = 95, gnames = "Ridge Regression", 
                rlinethick = 0.75, fsize = 8, showplots = FALSE)
auc2 = plot(lg_ridge_eval$roc) + ggtitle("") + 
  theme(plot.title = element_text(size = 14, face = "italic")) + 
  theme(legend.position = c(0.66, 0.14))
png(filename = paste0(temp_plots, "Gray_Snapper_Ridge_AUC.png"), 
    width = 3.15, height = 3.15, units = "in", res = 300)
auc2
dev.off()

# plot AUC-ROC Curve for testing area
lg_ridge_pred2 = predict(lg_ridge, newdata = lg_test[-663,], type = "prob")
lg_ridge_test_AUC = evalm(data.frame(lg_ridge_pred2, lg_test[-663,]$PRES),
                          showplots = FALSE)
lg_ridge_test_AUC$roc


#### BLUESTRIPED GRUNT ####
#### Standard Logistic Regression ####
# build standard logistic regression model (not penalized)
hs_glm = caret::train(PRES ~ ., hs_train, method = "glm", family = "binomial",
                      trControl = trainControl, metric = "ROC", maximize = TRUE)
hs_glm # ROC, sensitivity, specificity
summary(hs_glm) # coefficients, deviance, AIC

# plot AUC-ROC Curve for training CV
hs_glm_eval = evalm(hs_glm, percent = 95, gnames = "Haemulon sciurus", 
                    rlinethick = 0.75, fsize = 8, showplots = FALSE)
hs_glm_eval$roc

# binary presence/absence predictions for testing area
hs_glm_pred = predict(hs_glm, hs_test)
hs_glm_CM = confusionMatrix(hs_glm_pred, hs_test[-663,]$PRES, 
                            positive = "PRESENCE") # ignore the empty final row that stores only habitat data
hs_glm_CM # accuracy ~ 68% (68.43)

# variable importance
hs_glm_imp = varImp(hs_glm, useModel = TRUE, nonpara = TRUE, scale = TRUE)
View(hs_glm_imp$importance)

# plot AUC-ROC Curve for testing area
hs_glm_pred2 = predict(hs_glm, newdata = hs_test[-663,], type = "prob")
hs_glm_test_AUC = evalm(data.frame(hs_glm_pred2, hs_test[-663,]$PRES),
                        showplots = FALSE)
hs_glm_test_AUC$roc

#### Lasso Regularization ####
# convert training data to matrix format for glmnet
hs_train_mat = model.matrix(PRES ~ ., hs_train)
#convert class to numerical variable
hs_train_y = ifelse(hs_train$PRES == "PRESENCE", 1, 0)

# convert testing data to matrix format
hs_test_mat = model.matrix(PRES ~ ., hs_test)
#convert class to numerical variable
hs_test_y = ifelse(hs_test$PRES == "PRESENCE", 1, 0)

# fit lasso using 10-fold cross validation
hs_lasso_cv = cv.glmnet(hs_train_mat, hs_train_y, alpha = 1, nfolds = 10,
                        family = "binomial", standardize = T, 
                        type.measure = "deviance")
plot(hs_lasso_cv)

# fitted coefficients, using minimum lambda
coef(hs_lasso_cv, s = "lambda.min")
# penalty term using minimum lambda
sum(coef(hs_lasso_cv, s = "lambda.min")[-1] ^ 2)

# run penalized model in caret using minimum lambda value from cross validation
hs_lasso_grid = expand.grid(alpha = 1, lambda = hs_lasso_cv$lambda.min)
hs_lasso = train(PRES ~ ., hs_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = hs_lasso_grid, 
                 metric = "ROC", maximize = TRUE)
hs_lasso # ROC, sensitivity, specificity
max(hs_lasso[["results"]]$ROC) # max AUC-ROC
coef(hs_lasso$finalModel, hs_lasso$bestTune$lambda) # coefficients 

# predictions and confusion matrix
hs_lasso_pred = predict(hs_lasso, newdata = hs_test[-663,])
hs_lasso_CM = confusionMatrix(data = hs_lasso_pred, reference = hs_test[-663,]$PRES, 
                              positive = "PRESENCE")
hs_lasso_CM # accuracy ~ 69% (68.58)

hs_lasso_imp = varImp(hs_lasso, useModel = TRUE, nonpara = TRUE, scale = TRUE)
View(hs_lasso_imp$importance)

# plot AUC-ROC Curve for training CV
hs_lasso_eval = evalm(hs_lasso, percent = 95, gnames = "Lasso Regression", 
                      rlinethick = 0.75, fsize = 8, showplots = FALSE)
auc3 = plot(hs_lasso_eval$roc) + ggtitle("Haemulon sciurus") + 
  theme(plot.title = element_text(size = 14, face = "italic")) + 
  theme(legend.position = c(0.66, 0.14))
png(filename = paste0(temp_plots, "Bluestriped_Grunt_Lasso_AUC.png"), 
    width = 3.15, height = 3.15, units = "in", res = 300)
auc3
dev.off()

# plot AUC-ROC Curve for testing area
hs_lasso_pred2 = predict(hs_lasso, newdata = hs_test[-663,], type = "prob")
hs_lasso_test_AUC = evalm(data.frame(hs_lasso_pred2, hs_test[-663,]$PRES),
                          showplots = FALSE)
hs_lasso_test_AUC$roc

#### Ridge Regularization ####
# fit ridge using 10-fold cross validation
hs_ridge_cv = cv.glmnet(hs_train_mat, hs_train_y, alpha = 0, nfolds = 10,
                        family = "binomial", standardize = T, 
                        type.measure = "deviance")
plot(hs_ridge_cv)

# fitted coefficients, using minimum lambda
coef(hs_ridge_cv, s = "lambda.min")
# penalty term using minimum lambda
sum(coef(hs_ridge_cv, s = "lambda.min")[-1] ^ 2)

# run penalized model in caret using minimum lambda value from cross validation
hs_ridge_grid = expand.grid(alpha = 0, lambda = hs_ridge_cv$lambda.min)
hs_ridge = train(PRES ~ ., hs_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = hs_ridge_grid, 
                 metric = "ROC", maximize = TRUE)
hs_ridge # ROC, sensitivity, specificity
max(hs_ridge[["results"]]$ROC) # max AUC-ROC
coef(hs_ridge$finalModel, hs_ridge$bestTune$lambda) # coefficients 
# Mangrove distance, BPI Broad, and BPI Fine pretty close to 0

# predictions and confusion matrix
hs_ridge_pred = predict(hs_ridge, newdata = hs_test[-663,])
hs_ridge_CM = confusionMatrix(data = hs_ridge_pred, reference = hs_test[-663,]$PRES, 
                              positive = "PRESENCE")
hs_ridge_CM # accuracy ~ 70% (70.24)

hs_ridge_imp = varImp(hs_ridge, useModel = TRUE, nonpara = TRUE, scale = TRUE)
View(hs_ridge_imp$importance)

# plot AUC-ROC Curve for training CV
hs_ridge_eval = evalm(hs_lasso, percent = 95, gnames = "Ridge Regression", 
                      rlinethick = 0.75, fsize = 8, showplots = FALSE)
auc4 = plot(hs_ridge_eval$roc) + ggtitle(" ") + 
  theme(plot.title = element_text(size = 14, face = "italic")) + 
  theme(legend.position = c(0.66, 0.14))
png(filename = paste0(temp_plots, "Bluestriped_Grunt_Ridge_AUC.png"), 
    width = 3.15, height = 3.15, units = "in", res = 300)
auc4
dev.off()

# plot AUC-ROC Curve for testing area
hs_ridge_pred2 = predict(hs_ridge, newdata = hs_test[-663,], type = "prob")
hs_ridge_test_AUC = evalm(data.frame(hs_ridge_pred2, hs_test[-663,]$PRES),
                          showplots = FALSE)
hs_ridge_test_AUC$roc

#### COMPARE ACCURACIES ####
# gray snapper
lg_glm_CM # 83.84 %
lg_lasso_CM # 84.14 %
lg_ridge_CM # 84.59 %

# bluestriped grunt
hs_glm_CM # 68.43 %
hs_lasso_CM # 68.58 %
hs_ridge_CM # 70.24 %

#### COMPARE AUC-ROC AND PLOT ####
# gray snapper
max(lg_glm[["results"]]$ROC) # ~0.7482
max(lg_lasso[["results"]]$ROC) #  ~ 0.7517
max(lg_ridge[["results"]]$ROC) #  ~ 0.7498

# bluestriped grunt
max(hs_glm[["results"]]$ROC) # ~0.7987
max(hs_lasso[["results"]]$ROC) #  ~ 0.7995
max(hs_ridge[["results"]]$ROC) #  ~ 0.7939


#### CONTINUOUS (RASTER) SUITABILITY ####
### Training Area ####

# save PROJ.4 string for standard projection (ESPG:26958 NAD 83/Florida East)
my_crs = CRS("+init=epsg:26958")

# spatial predictors in training area
habitat_train = raster(paste0(train_wd, "Environmental/Habitat.asc"))
mg_dist_train = raster(paste0(train_wd, "Environmental/Mangrove_Dist.asc"))
depth_train = raster(paste0(train_wd, "Environmental/Depth.asc"))
slope_train = raster(paste0(train_wd, "Environmental/Slope.asc"))
curvature_train = raster(paste0(train_wd, "Environmental/Curvature.asc"))
bpi_fine_train = raster(paste0(train_wd, "Environmental/BPI_Fine.asc"))
bpi_broad_train = raster(paste0(train_wd, "Environmental/BPI_Broad.asc"))
rugosity_train = raster(paste0(train_wd, "Environmental/Rugosity.asc"))
sum_temp_train = raster(paste0(train_wd, "Environmental/Sum_Temp.asc"))
sum_sal_train = raster(paste0(train_wd, "Environmental/Sum_Sal.asc"))
win_sal_train = raster(paste0(train_wd, "Environmental/Win_Sal.asc"))

# define crs
crs(habitat_train) = my_crs
crs(mg_dist_train) = my_crs
crs(depth_train) = my_crs
crs(slope_train) = my_crs
crs(curvature_train) = my_crs
crs(bpi_fine_train) = my_crs
crs(bpi_broad_train) = my_crs
crs(rugosity_train) = my_crs
crs(sum_temp_train) = my_crs
crs(sum_sal_train) = my_crs
crs(win_sal_train) = my_crs

# create raster stack for training area
env_train = stack(x = c(habitat_train, mg_dist_train, depth_train,
                        slope_train, curvature_train, bpi_fine_train,
                        bpi_broad_train, rugosity_train, sum_temp_train, 
                        sum_sal_train, win_sal_train))

rm(list = c("habitat_train", "mg_dist_train", "depth_train",
            "slope_train", "curvature_train", "bpi_fine_train",
            "bpi_broad_train", "rugosity_train", "sum_temp_train", 
            "sum_sal_train", "win_sal_train"))


# make predictions across training area
# FYI: very time consuming step, each model takes up to 10 hours to build
library(snow)

beginCluster(n = 10)
lg_train_lasso_HSM = raster::predict(env_train, lg_lasso, fun = predict,
                                     type = "prob", filename = paste0(train_wd,
                                     "HSMs/Logistic_Regresssion/lg_lasso_HSM.asc"),
                                     format = "ascii", progress = "window", 
                                     overwrite = TRUE)
endCluster()


beginCluster(n = 10)
lg_train_ridge_HSM = raster::predict(env_train, lg_ridge, fun = predict,
                                     type = "prob", filename = paste0(train_wd,
                                    "HSMs/Logistic_Regresssion/lg_ridge_HSM.asc"),
                                    format = "ascii", progress = "window")
endCluster()

beginCluster(n = 10)
hs_train_lasso_HSM = raster::predict(env_train, hs_lasso, fun = predict,
                                     type = "prob", filename = paste0(train_wd,
                                    "HSMs/Logistic_Regresssion/hs_lasso_HSM.asc"),
                                     format = "ascii", progress = "window", 
                                     overwrite = TRUE)
endCluster()


beginCluster(n = 10)
hs_train_ridge_HSM = raster::predict(env_train, hs_ridge, fun = predict,
                                     type = "prob", filename = paste0(train_wd,
                                    "HSMs/Logistic_Regresssion/hs_ridge_HSM.asc"),
                                     format = "ascii", progress = "window")
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
sum_temp_test = raster(paste0(test_wd, "Environmental/Sum_Temp.asc"))
sum_sal_test = raster(paste0(test_wd, "Environmental/Sum_Sal.asc"))
win_sal_test = raster(paste0(test_wd, "Environmental/Win_Sal.asc"))

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

rm(list = c("habitat_test", "mg_dist_test", "depth_test",
            "slope_test", "curvature_test", "bpi_fine_test",
            "bpi_broad_test", "rugosity_test", "sum_temp_test", 
            "sum_sal_test", "win_sal_test"))

# make predictions across testing area
beginCluster(n = 10)
lg_test_lasso_HSM = raster::predict(env_test, lg_lasso, fun = predict,
                                     type = "prob", filename = paste0(test_wd,
                                     "HSMs/Logistic_Regresssion/lg_lasso_HSM.asc"),
                                     format = "ascii", progress = "window", 
                                     overwrite = TRUE)
endCluster()


beginCluster(n = 10)
lg_test_ridge_HSM = raster::predict(env_test, lg_ridge, fun = predict,
                                     type = "prob", filename = paste0(test_wd,
                                     "HSMs/Logistic_Regresssion/lg_ridge_HSM.asc"),
                                     format = "ascii", progress = "window",
                                     overwrite = TRUE)
endCluster()

# make predictions across testing area
beginCluster(n = 10)
hs_test_lasso_HSM = raster::predict(env_test, hs_lasso, fun = predict,
                                    type = "prob", filename = paste0(test_wd,
                                    "HSMs/Logistic_Regresssion/hs_lasso_HSM.asc"),
                                    format = "ascii", progress = "window", 
                                    overwrite = TRUE)
endCluster()

beginCluster(n = 10)
hs_test_ridge_HSM = raster::predict(env_test, hs_ridge, fun = predict,
                                    type = "prob", filename = paste0(test_wd,
                                    "HSMs/Logistic_Regresssion/hs_ridge_HSM.asc"),
                                    format = "ascii", progress = "window",
                                    overwrite = TRUE)
endCluster()
