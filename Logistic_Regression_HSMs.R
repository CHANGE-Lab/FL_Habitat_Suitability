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
libraries("glmnet", "tidyr", "caret", "MLeval", "dplyr", "raster", "PNWColors", "cowplot")


# presence-absence datasets from model training area and model testing area
# keep only the reduced covariate set (following multicollinearity assessment)

# gray snapper (Lutjanus griseus)
lg_train = read.csv(paste0(csv_wd, "Subadult_Gray_Snapper_Training_Data.csv")) %>%
  dplyr::select(PRES, Habitat, Mangrove_Dist, Depth, Slope, Curvature, Rugosity,
                BPI_Broad, BPI_Fine, Mean_Sum_Temp, Mean_Sum_Sal, Mean_Win_Sal) %>%
  dplyr::mutate(PRES = as.factor(PRES),
                Habitat = as.factor(Habitat))
lg_train$PRES = relevel(lg_train$PRES, ref = "ABSENCE") #  ensure that absence is the reference category

lg_test = read.csv(paste0(csv_wd, "Subadult_Gray_Snapper_Testing_Data.csv")) %>%
  dplyr::select(PRES, Habitat, Mangrove_Dist, Depth, Slope, Curvature, Rugosity,
                BPI_Broad, BPI_Fine, Mean_Sum_Temp, Mean_Sum_Sal, Mean_Win_Sal) %>%
  dplyr::mutate(PRES = as.factor(PRES),
                Habitat = as.factor(Habitat))
lg_test$PRES = relevel(lg_test$PRES, ref = "ABSENCE") #  ensure that absence is the reference category


# bluestriped grunt (Haemulon sciurus)
hs_train = read.csv(paste0(csv_wd, "Subadult_Bluestriped_Grunt_Training_Data.csv"))%>%
  dplyr::select(PRES, Habitat, Mangrove_Dist, Depth, Slope, Curvature, Rugosity,
                BPI_Broad, BPI_Fine, Mean_Sum_Temp, Mean_Sum_Sal, Mean_Win_Sal) %>%
  dplyr::mutate(PRES = as.factor(PRES),
                Habitat = as.factor(Habitat))
hs_train$PRES = relevel(hs_train$PRES, ref = "ABSENCE") # ensure that absence is the reference category

hs_test = read.csv(paste0(csv_wd, "Subadult_Bluestriped_Grunt_Testing_Data.csv"))%>%
  dplyr::select(PRES, Habitat, Mangrove_Dist, Depth, Slope, Curvature, Rugosity,
                BPI_Broad, BPI_Fine, Mean_Sum_Temp, Mean_Sum_Sal, Mean_Win_Sal) %>%
  dplyr::mutate(PRES = as.factor(PRES),
                Habitat = as.factor(Habitat))
hs_test$PRES = relevel(hs_test$PRES, ref = "ABSENCE") # ensure that absence is the reference category


# add an empty row with mangrove habitat (ID# 11) to the test data so an equal 
# number of habitat levels are present in the training and testing data
add_mg = data.frame(NA, as.factor(11), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
names(add_mg) = c("PRES", "Habitat", "Mangrove_Dist", "Depth", "Slope", "Curvature",
                  "Rugosity", "BPI_Broad", "BPI_Fine", "Mean_Sum_Temp", "Mean_Sum_Sal", "Mean_Win_Sal")
lg_test = rbind(lg_test, add_mg) # adds an extra row that is empty other than 11 for habitat
hs_test = rbind(hs_test, add_mg)

# before getting started
set.seed(42)

#### GRAY SNAPPER ####
#### Standard Logistic Regression ####
trainControl = trainControl(method = "cv", number = 10,
                            summaryFunction = twoClassSummary, 
                            classProbs = TRUE, savePredictions = TRUE)
lg_glm = caret::train(PRES ~ ., lg_train, method = "glm", family = "binomial",
                      trControl = trainControl, metric = "ROC", maximize = TRUE,
                      preProcess = c("center", "scale"))

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
lg_glm_probs = predict(lg_glm, newdata = lg_test[-663,], type = "prob") # get rid of last empty row

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
lg_glm_pred = predict(lg_glm, lg_test[-663,])
lg_glm_CM = confusionMatrix(lg_glm_pred, lg_test[-663,]$PRES,
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
lg_train_mat = model.matrix(PRES ~ ., lg_train)
# convert class to numerical variable
lg_train_y = ifelse(lg_train$PRES == "PRESENCE", 1, 0)

# convert testing data to matrix format
lg_test_mat = model.matrix(PRES ~ ., lg_test)
# convert class to numerical variable
lg_test_y = ifelse(lg_test$PRES == "PRESENCE", 1, 0)

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
lg_lasso = train(PRES ~ ., lg_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = lg_lasso_grid, 
                 metric = "ROC", maximize = TRUE,  preProcess = c("center", "scale"))
lg_lasso # ROC, sensitivity, specificity
max(lg_lasso[["results"]]$ROC) # max AUC-ROC
coef(lg_lasso$finalModel, lg_lasso$bestTune$lambda) # coefficients log-odds
exp(coef(lg_lasso$finalModel, lg_lasso$bestTune$lambda)) # odds ratio
# Slope, curvature, BPI Broad, Sum Sal, and a few habitat classes driven to 0


# plot AUC-ROC Curve from training CV
lg_lasso_eval = evalm(lg_lasso, percent = 95, gnames = "Lasso Regression", 
                      rlinethick = 0.75, fsize = 8, showplots = FALSE)
auc1 = plot(lg_lasso_eval$roc) + ggtitle("Lutjanus griseus") + 
  theme(plot.title = element_text(size = 14, face = "italic")) + 
  theme(legend.position = c(0.66, 0.14))
png(filename = paste0(temp_plots, "Gray_Snapper_Lasso_AUC.png"), 
    width = 3.15, height = 3.15, units = "in", res = 300)
auc1
dev.off()

# predict probability of presence / absence across testing area validation data
lg_lasso_probs = predict(lg_lasso, newdata = lg_test[-663,], type = "prob")

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
lg_lasso_pred = predict(lg_lasso, newdata = lg_test[-663,])
lg_lasso_CM = confusionMatrix(data = lg_lasso_pred, reference = lg_test[-663,]$PRES,
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
lg_ridge = train(PRES ~ ., lg_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = lg_ridge_grid, 
                 metric = "ROC", maximize = TRUE, preProcess = c("center", "scale"))
lg_ridge # ROC, sensitivity, specificity
max(lg_ridge[["results"]]$ROC) # max AUC-ROC
coef(lg_ridge$finalModel, lg_ridge$bestTune$lambda) # coefficients log-odds
exp(coef(lg_ridge$finalModel, lg_ridge$bestTune$lambda)) # odds ratio
# Mangrove distance, BPI Broad, and BPI Fine pretty close to 0

# AUC-ROC Curve from training CV
lg_ridge_eval = evalm(lg_ridge, percent = 95, gnames = "Ridge Regression", 
                      rlinethick = 0.75, fsize = 8, showplots = FALSE)
auc2 = plot(lg_ridge_eval$roc) + ggtitle("") + 
  theme(plot.title = element_text(size = 14, face = "italic")) + 
  theme(legend.position = c(0.66, 0.14))
png(filename = paste0(temp_plots, "Gray_Snapper_Ridge_AUC.png"), 
    width = 3.15, height = 3.15, units = "in", res = 300)
auc2
dev.off()

# predict probability of presence / absence across testing area validation data
lg_ridge_probs = predict(lg_ridge, newdata = lg_test[-663,], type = "prob")

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
lg_ridge_pred = predict(lg_ridge, newdata = lg_test[-663,])
lg_ridge_CM = confusionMatrix(data = lg_ridge_pred, reference = lg_test[-663,]$PRES,
                              positive = "PRESENCE")
lg_ridge_CM 


#### BLUESTRIPED GRUNT ####
#### Standard Logistic Regression ####
hs_glm = caret::train(PRES ~ ., hs_train, method = "glm", family = "binomial",
                      trControl = trainControl, metric = "ROC", maximize = TRUE,
                      preProcess = c("center", "scale"))
hs_glm # ROC, sensitivity, specificity
summary(hs_glm) # coefficients, deviance, AIC

# plot AUC-ROC Curve from training CV
hs_glm_eval = evalm(hs_glm, percent = 95, gnames = "Haemulon sciurus", 
                    rlinethick = 0.75, fsize = 8, showplots = FALSE)
hs_glm_eval$roc

# predict probability of presence / absence across testing area validation data
hs_glm_probs = predict(hs_glm, newdata = hs_test[-663,], type = "prob")

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
hs_glm_pred = predict(hs_glm, hs_test)
hs_glm_CM = confusionMatrix(hs_glm_pred, hs_test[-663,]$PRES,
                            positive = "PRESENCE")
hs_glm_CM 


#### Lasso Regularization ####
# convert training data to matrix format for glmnet
hs_train_mat = model.matrix(PRES ~ ., hs_train)
#convert class to numerical variable
hs_train_y = ifelse(hs_train$PRES == "PRESENCE", 1, 0)

# convert testing data to matrix format
hs_test_mat = model.matrix(PRES ~ ., hs_test)
#convert class to numerical variable
hs_test_y = ifelse(hs_test$PRES == "PRESENCE", 1, 0)

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
hs_lasso = train(PRES ~ ., hs_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = hs_lasso_grid, 
                 metric = "ROC", maximize = TRUE, preProcess = c("center", "scale"))
hs_lasso # ROC, sensitivity, specificity
max(hs_lasso[["results"]]$ROC) # max AUC-ROC
coef(hs_lasso$finalModel, hs_lasso$bestTune$lambda) # coefficients log-odds
exp(coef(hs_lasso$finalModel, hs_lasso$bestTune$lambda)) # odds-ratio

# plot AUC-ROC Curve from training CV
hs_lasso_eval = evalm(hs_lasso, percent = 95, gnames = "Lasso Regression", 
                      rlinethick = 0.75, fsize = 8, showplots = FALSE)
auc3 = plot(hs_lasso_eval$roc) + ggtitle("Haemulon sciurus") + 
  theme(plot.title = element_text(size = 14, face = "italic")) + 
  theme(legend.position = c(0.66, 0.14))
png(filename = paste0(temp_plots, "Bluestriped_Grunt_Lasso_AUC.png"), 
    width = 3.15, height = 3.15, units = "in", res = 300)
auc3
dev.off()

# predict probability of presence / absence across testing area validation data
hs_lasso_probs = predict(hs_lasso, newdata = hs_test[-663,], type = "prob")

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
hs_lasso_pred = predict(hs_lasso, newdata = hs_test[-663,])
hs_lasso_CM = confusionMatrix(data = hs_lasso_pred, reference = hs_test[-663,]$PRES,
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
hs_ridge = train(PRES ~ ., hs_train, method = "glmnet", family = "binomial",
                 trControl = trainControl, tuneGrid = hs_ridge_grid, 
                 metric = "ROC", maximize = TRUE, preProcess = c("center", "scale"))
hs_ridge # ROC, sensitivity, specificity
max(hs_ridge[["results"]]$ROC) # max AUC-ROC
coef(hs_ridge$finalModel, hs_ridge$bestTune$lambda) # coefficients log-odds
exp(coef(hs_ridge$finalModel, hs_ridge$bestTune$lambda)) # odds ratio
# Mangrove distance, BPI Broad, and BPI Fine pretty close to 0

# plot AUC-ROC Curve from training CV
hs_ridge_eval = evalm(hs_lasso, percent = 95, gnames = "Ridge Regression", 
                      rlinethick = 0.75, fsize = 8, showplots = FALSE)
auc4 = plot(hs_ridge_eval$roc) + ggtitle(" ") + 
  theme(plot.title = element_text(size = 14, face = "italic")) + 
  theme(legend.position = c(0.66, 0.14))
png(filename = paste0(temp_plots, "Bluestriped_Grunt_Ridge_AUC.png"), 
    width = 3.15, height = 3.15, units = "in", res = 300)
auc4
dev.off()

# predict probability of presence / absence across testing area validation data
hs_ridge_probs = predict(hs_ridge, newdata = hs_test[-663,], type = "prob")

# predict binary response (presence/absence) across testing area and calculate map accuracy
# default threshold for suitability is 0.5
hs_ridge_pred = predict(hs_ridge, newdata = hs_test[-663,])
hs_ridge_CM = confusionMatrix(data = hs_ridge_pred, reference = hs_test[-663,]$PRES,
                              positive = "PRESENCE")
hs_ridge_CM

#### COMPARE ACCURACIES ####
# gray snapper
lg_glm_CM # 83.84 %
lg_lasso_CM # 84.14 %
lg_ridge_CM # 84.59 %

# bluestriped grunt
hs_glm_CM # 68.43 %
hs_lasso_CM # 68.58 %
hs_ridge_CM # 70.24 %

#### COMPARE AUC-ROC ####
# gray snapper
max(lg_glm[["results"]]$ROC) # ~0.7482
max(lg_lasso[["results"]]$ROC) #  ~ 0.7517
max(lg_ridge[["results"]]$ROC) #  ~ 0.7503

# bluestriped grunt
max(hs_glm[["results"]]$ROC) # ~0.7991
max(hs_lasso[["results"]]$ROC) #  ~ 0.7996
max(hs_ridge[["results"]]$ROC) #  ~ 0.7937

#### REGRESSION COEFFICIENTS ####
# look at one of the dgC matrices to find out the order of predictors
coef(lg_lasso$finalModel, lg_lasso$bestTune$lambda) 
coef_names = as.data.frame(c("Intercept", "Scattered Coral/Rock", "Continuous Seagrass",
                             "Discontinuous Seagrass", "Unconsolidated Sediment",
                             "Aggregate Reef", "Pavement", "Reef Rubble", "Mangrove",
                             "Mangrove Distance", "Depth", "Slope", "Curvature",
                             "Rugosity", "BPI Broad", "BPI Fine", "Summer Temperature",
                             "Summer Salinity", "Winter Salinity"))

lg_lasso_coef = as.data.frame(as.vector(coef(lg_lasso$finalModel, lg_lasso$bestTune$lambda)))
lg_ridge_coef = as.data.frame(as.vector(coef(lg_ridge$finalModel, lg_ridge$bestTune$lambda)))
hs_lasso_coef = as.data.frame(as.vector(coef(hs_lasso$finalModel, hs_lasso$bestTune$lambda)))
hs_ridge_coef = as.data.frame(as.vector(coef(hs_ridge$finalModel, hs_ridge$bestTune$lambda)))

# merge and add species info
lg_lasso_coef = cbind(coef_names, lg_lasso_coef)
colnames(lg_lasso_coef) = c("Variable", "Coefficient")
lg_lasso_coef$Species = rep("Lutjanus griseus", nrow(lg_lasso_coef))

lg_ridge_coef = cbind(coef_names, lg_ridge_coef)
colnames(lg_ridge_coef) = c("Variable", "Coefficient")
lg_ridge_coef$Species = rep("Lutjanus griseus", nrow(lg_ridge_coef))

hs_lasso_coef = cbind(coef_names, hs_lasso_coef)
colnames(hs_lasso_coef) = c("Variable", "Coefficient")
hs_lasso_coef$Species = rep("Haemulon sciurus", nrow(hs_lasso_coef))

hs_ridge_coef = cbind(coef_names, hs_ridge_coef)
colnames(hs_ridge_coef) = c("Variable", "Coefficient")
hs_ridge_coef$Species = rep("Haemulon sciurus", nrow(hs_ridge_coef))

# combine lasso coefficients for both species, repeat for ridge and save
lasso_coef = rbind(lg_lasso_coef, hs_lasso_coef)
ridge_coef = rbind(lg_ridge_coef, hs_ridge_coef)

# the original depth data are negative, meaning larger negative numbers (deeper water)
# are technically "smaller" than less negative (shallower) values. This is a bit confusing 
# when looking at coefficients, because technically a positive depth coefficient means that
# habitat suitability increases as depth values get larger (shallower). Instead, flip
# the sign of the coefficient so that they represent the change in suitability at increasing depths
lasso_coef[11, 2] = as.numeric(-1*lasso_coef[11, 2])
lasso_coef[30, 2] = as.numeric(-1*lasso_coef[30, 2])

ridge_coef[11, 2] = as.numeric(-1*ridge_coef[11, 2])
ridge_coef[30, 2] = as.numeric(-1*ridge_coef[30, 2])

# finally add the odds ratios to each table and save
lasso_coef$Odds_Ratio = exp(lasso_coef$Coefficient)
ridge_coef$Odds_Ratio = exp(ridge_coef$Coefficient)
write.csv(lasso_coef, paste0(csv_wd, "Lasso_Coefficients.csv"))
write.csv(ridge_coef, paste0(csv_wd, "Ridge_Coefficients.csv"))

# plotting coefficients
my_pal = pnw_palette("Bay",8)
my_pal

lasso_coef_plot = ggplot(data = lasso_coef, aes(x = Coefficient, y = Variable, fill = Species)) +
  geom_vline(xintercept = 0, color = "gray") + geom_point(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  labs(y = "", x = "Lasso regression coefficient") + xlim(-2,2)  +
  scale_y_discrete(limits = c("Intercept", "Winter Salinity", "Summer Salinity", "Summer Temperature", "BPI Broad", "BPI Fine",
                              "Rugosity", "Curvature", "Slope", "Depth", "Mangrove Distance", "Discontinuous Seagrass", 
                              "Continuous Seagrass", "Mangrove", "Pavement", "Aggregate Reef","Reef Rubble", 
                              "Scattered Coral/Rock", "Unconsolidated Sediment")) +
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
  scale_y_discrete(limits = c("Intercept", "Winter Salinity", "Summer Salinity", "Summer Temperature", "BPI Broad", "BPI Fine",
                              "Rugosity", "Curvature", "Slope", "Depth", "Mangrove Distance", "Discontinuous Seagrass", 
                              "Continuous Seagrass", "Mangrove", "Pavement", "Aggregate Reef","Reef Rubble", 
                              "Scattered Coral/Rock", "Unconsolidated Sediment")) +
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
names = data.frame("Variable" = c("Scattered Coral/Rock", "Continuous Seagrass", "Discontinuous Seagrass",
                                  "Unconsolidated Sediment", "Aggregate Reef", "Pavement",
                                  "Reef Rubble", "Mangrove", "Mangrove Distance", "Depth", 
                                  "Slope", "Curvature", "Rugosity", "BPI Broad", "BPI Fine", 
                                  "Summer Temperature", "Summer Salinity", "Winter Salinity"))
lg_lasso_vimp = cbind(lg_lasso_vimp, names)

lg_ridge_vimp = cbind(caret::varImp(lg_ridge, lambda = lg_ridge$bestTune$lambda)$importance, names)
hs_lasso_vimp = cbind(caret::varImp(hs_lasso, lambda = hs_lasso$bestTune$lambda)$importance, names)
hs_ridge_vimp = cbind(caret::varImp(hs_ridge, lambda = hs_ridge$bestTune$lambda)$importance, names)

lg_lasso_vimp$Species = rep("Lutjanus griseus", nrow(lg_lasso_vimp))
lg_ridge_vimp$Species = rep("Lutjanus griseus", nrow(lg_ridge_vimp))
hs_lasso_vimp$Species = rep("Haemulon sciurus", nrow(hs_lasso_vimp))
hs_ridge_vimp$Species = rep("Haemulon sciurus", nrow(hs_ridge_vimp))

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
  scale_x_discrete(limits = c("Winter Salinity", "Summer Salinity", "Summer Temperature", "BPI Broad", "BPI Fine",
                              "Rugosity", "Curvature", "Slope", "Depth", "Mangrove Distance", "Discontinuous Seagrass", 
                              "Continuous Seagrass", "Mangrove", "Pavement", "Aggregate Reef","Reef Rubble", 
                              "Scattered Coral/Rock", "Unconsolidated Sediment")) +
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
  scale_x_discrete(limits = c("Winter Salinity", "Summer Salinity", "Summer Temperature", "BPI Broad", "BPI Fine",
                              "Rugosity", "Curvature", "Slope", "Depth", "Mangrove Distance", "Discontinuous Seagrass", 
                              "Continuous Seagrass", "Mangrove", "Pavement", "Aggregate Reef","Reef Rubble", 
                              "Scattered Coral/Rock", "Unconsolidated Sediment")) +
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
sum_temp_train = raster(paste0(train_wd, "Environmental/Mean_Sum_Temp.asc"))
sum_sal_train = raster(paste0(train_wd, "Environmental/Mean_Sum_Sal.asc"))
win_sal_train = raster(paste0(train_wd, "Environmental/Mean_Win_Sal.asc"))

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

# save folder path
HSMs = "Z:/Courtney/Stuart_MSc_Ch1/HSMs/"

# make predictions across training area 
# FYI: very time consuming step, each model takes up to 10 hours to build
library(snow)

# IMPORTANT: specify (1 - predict) because the raster package will use the first class as its
# output, and that is the absence reference class in this case. We want probability of presence!
# you'll notice this small caveat if you open the lg_lasso_probs df - absence is first.
beginCluster(n = 10)
lg_train_lasso_HSM = writeRaster((1-(raster::predict(env_train, lg_lasso, type = "prob"))),
                                 filename = paste0(HSMs, 
                                 "Logistic_Regression/Subadult_Gray_Snapper_Lasso_Train.asc"),
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
