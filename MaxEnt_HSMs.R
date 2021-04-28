# MaxEnt predictions and response curves


#### SET-UP ####
# data directories
setwd("Z:/Courtney/Stuart_MSc_Ch1/") # main project folder
temp_wd = "Z:/Courtney/Stuart_MSc_Ch1/Temp/" # temporary files
csv_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Data/" # for writing tabular data
train_wd = "Z:/Courtney/Stuart_MSc_Ch1/Modeling_Data/Training/" # for training area data
test_wd = "Z:/Courtney/Stuart_MSc_Ch1/Modeling_Data/Testing/" # for testing area data
plots_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Figures/" # for figures
temp_plots = "Z:/Courtney/Stuart_MSc_Ch1/Plots/" # temporary plots for Courtney
HSMs = "Z:/Courtney/Stuart_MSc_Ch1/HSMs/"


# libraries
library(easypackages)
libraries("glmnet", "tidyverse", "tidyr", "caret", "MLeval", "dplyr", "raster", 
          "sf", "ggplot2", "fs", "purrr", "stringr", "PNWColors", "cowplot",
          "dismo", "conflicted")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# save PROJ.4 string for standard projection (ESPG:26958 NAD 83/Florida East)
my_crs = CRS("+init=epsg:26958")

#### GRAY SNAPPER ####

# presence-absence datasets from model training area and model testing area
lg_train_pa = read.csv(paste0(train_wd, "Occurrences/Presence_Absence/Subadult_Gray_Snapper_PA_Train.csv"))
lg_train_pa$PRES = as.factor(ifelse(lg_train_pa$PRES == 1, "PRESENCE", "ABSENCE"))
lg_train_coords = lg_train_pa %>% dplyr::select(LON_M, LAT_M)

lg_test_pa = read.csv(paste0(test_wd, "Occurrences/Presence_Absence/Subadult_Gray_Snapper_PA_Test.csv"))
lg_test_pa$PRES = as.factor(ifelse(lg_test_pa$PRES == 1, "PRESENCE", "ABSENCE"))
lg_test_coords = lg_test_pa %>% dplyr::select(LON_M, LAT_M)


# lg maxents
lg_maxent_train = raster(paste0(HSMs, "MaxEnt/Subadult_Gray_Snapper/LUT_GRIS_avg.asc"))
lg_maxent_test = raster(paste0(HSMs, "MaxEnt/Subadult_Gray_Snapper/LUT_GRIS_Environmental_avg.asc"))

# maxent predictions (50% suitability threshold)
lg_maxent_train_prob = cbind(lg_train_coords, raster::extract(lg_maxent_train, lg_train_coords)) 
lg_maxent_test_prob = cbind(lg_test_coords, raster::extract(lg_maxent_test, lg_test_coords)) 

# default class predictions at suitability threshold 0.5
lg_maxent_train_pred = as.factor(ifelse(lg_maxent_train_prob$`raster::extract(lg_maxent_train, lg_train_coords)`>= 0.5,
                              "PRESENCE", "ABSENCE"))
lg_maxent_test_pred = as.factor(ifelse(lg_maxent_test_prob$`raster::extract(lg_maxent_test, lg_test_coords)`>= 0.5, 
                             "PRESENCE", "ABSENCE"))

# confusion matrix
lg_maxent_train_CM = confusionMatrix(data = lg_maxent_train_pred, reference = lg_train_pa$PRES, 
                                     positive = "PRESENCE")
lg_maxent_train_CM 
  
lg_maxent_test_CM = confusionMatrix(data = lg_maxent_test_pred, reference = lg_test_pa$PRES, 
                                    positive = "PRESENCE")
lg_maxent_test_CM 


# AUC for testing area
# probability of presence
lg_test_prob_pres = lg_maxent_test_prob$`raster::extract(lg_maxent_test, lg_test_coords)`
# probability of absence
lg_test_prob_abs = rep(1 - lg_test_prob_pres)
lg_test_probs = as.data.frame(cbind(lg_test_prob_pres, lg_test_prob_abs)) %>%
  dplyr::rename(PRESENCE = lg_test_prob_pres,
                ABSENCE = lg_test_prob_abs)

lg_maxent_test_AUC = evalm(data.frame(lg_test_probs, lg_test_pa$PRES),
                          showplots = FALSE)

# maxent predictions (Maximize sum of training sensitivity + specificity)
lg_maxent_train_predSSS = as.factor(ifelse(lg_maxent_train_prob$`raster::extract(lg_maxent_train, lg_train_coords)`>= 0.6459,
                                        "PRESENCE", "ABSENCE"))
lg_maxent_test_predSSS = as.factor(ifelse(lg_maxent_test_prob$`raster::extract(lg_maxent_test, lg_test_coords)`>= 0.6459, 
                                       "PRESENCE", "ABSENCE"))

# confusion matrix
lg_maxent_train_CM_SSS = confusionMatrix(data = lg_maxent_train_predSSS, reference = lg_train_pa$PRES, 
                                         positive = "PRESENCE")
lg_maxent_train_CM_SSS 

lg_maxent_test_CM_SSS = confusionMatrix(data = lg_maxent_test_predSSS, reference = lg_test_pa$PRES, 
                                        positive = "PRESENCE")
lg_maxent_test_CM_SSS 


#### BLUESTRIPED GRUNT ####

# presence-absence datasets from model training area and model testing area
hs_train_pa = read.csv(paste0(train_wd, "Occurrences/Presence_Absence/Subadult_Bluestriped_Grunt_PA_Train.csv"))
hs_train_pa$PRES = as.factor(ifelse(hs_train_pa$PRES == 1, "PRESENCE", "ABSENCE"))
hs_train_coords = hs_train_pa %>% dplyr::select(LON_M, LAT_M)

hs_test_pa = read.csv(paste0(test_wd, "Occurrences/Presence_Absence/Subadult_Bluestriped_Grunt_PA_Test.csv"))
hs_test_pa$PRES = as.factor(ifelse(hs_test_pa$PRES == 1, "PRESENCE", "ABSENCE"))
hs_test_coords = hs_test_pa %>% dplyr::select(LON_M, LAT_M)

# hs maxents
hs_maxent_train = raster(paste0(HSMs, "MaxEnt/Subadult_Bluestriped_Grunt/HAE_SCIU_avg.asc"))
hs_maxent_test = raster(paste0(HSMs, "MaxEnt/Subadult_Bluestriped_Grunt/HAE_SCIU_Environmental_avg.asc"))

# maxent predictions (50% suitability threshold)
hs_maxent_train_prob = cbind(hs_train_coords, raster::extract(hs_maxent_train, hs_train_coords)) 
hs_maxent_test_prob = cbind(hs_test_coords, raster::extract(hs_maxent_test, hs_test_coords)) 

# default class predictions at suitability threshold 0.5
hs_maxent_train_pred = as.factor(ifelse(hs_maxent_train_prob$`raster::extract(hs_maxent_train, hs_train_coords)`>=0.5,
                                        "PRESENCE", "ABSENCE"))
hs_maxent_test_pred = as.factor(ifelse(hs_maxent_test_prob$`raster::extract(hs_maxent_test, hs_test_coords)`>=0.5, 
                                       "PRESENCE", "ABSENCE"))

# confusion matrix
hs_maxent_train_CM = confusionMatrix(data = hs_maxent_train_pred, reference = hs_train_pa$PRES, 
                                     positive = "PRESENCE")
hs_maxent_train_CM 

hs_maxent_test_CM = confusionMatrix(data = hs_maxent_test_pred, reference = hs_test_pa$PRES, 
                                    positive = "PRESENCE")
hs_maxent_test_CM 

# AUC for testing area
# probability of presence
hs_test_prob_pres = hs_maxent_test_prob$`raster::extract(hs_maxent_test, hs_test_coords)`
# probability of absence
hs_test_prob_abs = rep(1 - hs_test_prob_pres)
hs_test_probs = as.data.frame(cbind(hs_test_prob_pres, hs_test_prob_abs)) %>%
  dplyr::rename(PRESENCE = hs_test_prob_pres,
                ABSENCE = hs_test_prob_abs)

hs_maxent_test_AUC = evalm(data.frame(hs_test_probs, hs_test_pa$PRES),
                           showplots = FALSE)

# maxent predictions (Maximize sum of training sensitivity + specificity)
hs_maxent_train_predSSS = as.factor(ifelse(hs_maxent_train_prob$`raster::extract(hs_maxent_train, hs_train_coords)`>= 0.5892,
                                           "PRESENCE", "ABSENCE"))
hs_maxent_test_predSSS = as.factor(ifelse(hs_maxent_test_prob$`raster::extract(hs_maxent_test, hs_test_coords)`>= 0.5892, 
                                          "PRESENCE", "ABSENCE"))

# confusion matrix
hs_maxent_train_CM_SSS = confusionMatrix(data = hs_maxent_train_predSSS, reference = hs_train_pa$PRES, 
                                         positive = "PRESENCE")
hs_maxent_train_CM_SSS 

hs_maxent_test_CM_SSS = confusionMatrix(data = hs_maxent_test_predSSS, reference = hs_test_pa$PRES, 
                                        positive = "PRESENCE")
hs_maxent_test_CM_SSS 


#### RESPONSE CURVES ####
### HABITAT ####
# gray snapper
URM = read.csv(paste0(csv_wd, "URM_ClassLv1_IDs.csv"))

lg_hab = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                    pattern = "*Habitat_only.csv")

for(i in 1:length(lg_hab)) {                              
  assign(paste0("lg_hab", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                          lg_hab[i])))

}

hab_lg_list = list(lg_hab1, lg_hab2, lg_hab3, lg_hab4, lg_hab5, lg_hab6,
                    lg_hab7, lg_hab8, lg_hab9, lg_hab10)
# make a list that corresponds with the ID number you want for each dataframe
hab_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# Iterate across each dataframe in the list
hab_map = Map(cbind, hab_lg_list, ID = hab_lg_IDs)

lg_hab_df = do.call("rbind", list(hab_map[[1]], hab_map[[2]], hab_map[[3]],
                                  hab_map[[4]], hab_map[[5]], hab_map[[6]],
                                  hab_map[[7]], hab_map[[8]], hab_map[[9]],
                                  hab_map[[10]]))

rm(list = c("lg_hab1", "lg_hab2", "lg_hab3", "lg_hab4", "lg_hab5", "lg_hab6", "lg_hab7",
            "lg_hab8", "lg_hab9", "lg_hab10", "hab_lg_list", "hab_lg_IDs", "hab_map"))

lg_hab_df$ClassLv1_ID = URM$ClassLv1[match(lg_hab_df$x, URM$ID)] # match IDs
lg_hab_df$ID = as.factor(lg_hab_df$ID)
lg_hab_df$Species = rep("Lutjanus griseus", nrow(lg_hab_df))

lg_hab_df = lg_hab_df %>%
  rename(Variable = variable,
         Hab_ID = x,
         HSI = y,
         Fold = ID,
         Hab_Name = ClassLv1_ID)

lg_hab_summary = lg_hab_df %>%
  group_by(Hab_Name) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>%
  ungroup() %>%
  dplyr::select(Species, Variable, Hab_ID, Hab_Name, Mean_HSI, SD_HSI) %>%
  distinct()
                                                                                                           
#  bluestriped grunts
hs_hab = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                    pattern = "*Habitat_only.csv")

for(i in 1:length(hs_hab)) {                              
  assign(paste0("hs_hab", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_hab[i])))
  
}

hab_hs_list = list(hs_hab1, hs_hab2, hs_hab3, hs_hab4, hs_hab5, hs_hab6,
                   hs_hab7, hs_hab8, hs_hab9, hs_hab10)
# make a list that corresponds with the ID number you want for each dataframe
hab_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# Iterate across each dataframe in the list
hab_map = Map(cbind, hab_hs_list, ID = hab_hs_IDs)

hs_hab_df = do.call("rbind", list(hab_map[[1]], hab_map[[2]], hab_map[[3]],
                                  hab_map[[4]], hab_map[[5]], hab_map[[6]],
                                  hab_map[[7]], hab_map[[8]], hab_map[[9]],
                                  hab_map[[10]]))

rm(list = c("hs_hab1", "hs_hab2", "hs_hab3", "hs_hab4", "hs_hab5", "hs_hab6", "hs_hab7",
            "hs_hab8", "hs_hab9", "hs_hab10", "hab_hs_list", "hab_hs_IDs", "hab_map"))

hs_hab_df$ClassLv1_ID = URM$ClassLv1[match(hs_hab_df$x, URM$ID)] # match IDs
hs_hab_df$ID = as.factor(hs_hab_df$ID)
hs_hab_df$Species = rep("Haemulon sciurus", nrow(hs_hab_df))

hs_hab_df = hs_hab_df %>%
  rename(Variable = variable,
         Hab_ID = x,
         HSI = y,
         Fold = ID,
         Hab_Name = ClassLv1_ID)

hs_hab_summary = hs_hab_df %>%
  group_by(Hab_Name) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>%
  ungroup() %>%
  dplyr::select(Species, Variable, Hab_ID, Hab_Name, Mean_HSI, SD_HSI) %>%
  distinct()

# combined habitat results
hab_summary = rbind(lg_hab_summary, hs_hab_summary)

# palette
my_pal = pnw_palette("Bay",8)
my_pal

hab_plot = ggplot(hab_summary, aes(x = Hab_Name, y = Mean_HSI, fill = Species, width = 0.5)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  xlab(" ") + ylab("Predicted Habitat Suitability") +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_errorbar(aes(ymin = Mean_HSI-SD_HSI, ymax = Mean_HSI+SD_HSI), width = 0.3,
                position = position_dodge(0.5)) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(), 
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.text.y = element_text(color = "black"),
                     axis.line = element_line(color = "black"),
                     axis.text.x = element_text(color = "black", angle = 30, vjust = 1,
                                                hjust = 1, size = 7)) + 
  theme(legend.text = element_text(face = "italic")) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "),
                                                 width = 25))
hab_plot

ggsave(filename = paste0(temp_plots, "Habitat_Curve.png"), hab_plot, width = 5,
       height = 3.15, units = "in", dpi = 450)

#### MANGROVE DISTANCE ####
# gray snapper
lg_mg = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                    pattern = "*Mangrove_Dist_only.csv")

for(i in 1:length(lg_mg)) {                              
  assign(paste0("lg_mg", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                         lg_mg[i])))
  
}

mg_lg_list = list(lg_mg1, lg_mg2, lg_mg3, lg_mg4, lg_mg5, lg_mg6,
                   lg_mg7, lg_mg8, lg_mg9, lg_mg10)
mg_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
mg_map = Map(cbind, mg_lg_list, ID = mg_lg_IDs)

lg_mg_df = do.call("rbind", list(mg_map[[1]], mg_map[[2]], mg_map[[3]],
                                  mg_map[[4]], mg_map[[5]], mg_map[[6]],
                                  mg_map[[7]], mg_map[[8]], mg_map[[9]],
                                  mg_map[[10]]))
lg_mg_df$ID = as.factor(lg_mg_df$ID)


rm(list = c("lg_mg1", "lg_mg2", "lg_mg3", "lg_mg4", "lg_mg5", "lg_mg6", "lg_mg7", 
            "lg_mg8","lg_mg9", "lg_mg10", "mg_lg_list", "mg_lg_IDs", "mg_map"))

lg_mg_df$Species = rep("Lutjanus griseus", nrow(lg_mg_df))

lg_mg_df = lg_mg_df %>%
  rename(Variable = variable,
         MG_Dist = x,
         HSI = y,
         Fold = ID) %>%
  filter(MG_Dist >= 0)

lg_mg_summary = lg_mg_df %>%
  group_by(MG_Dist) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>%
  ungroup() %>%
  dplyr::select(Species, Variable, MG_Dist, Mean_HSI, SD_HSI) %>%
  distinct()

#  bluestriped grunts
hs_mg = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                   pattern = "*Mangrove_Dist_only.csv")

for(i in 1:length(hs_mg)) {                              
  assign(paste0("hs_mg", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_mg[i])))
  
}

mg_hs_list = list(hs_mg1, hs_mg2, hs_mg3, hs_mg4, hs_mg5, hs_mg6,
                  hs_mg7, hs_mg8, hs_mg9, hs_mg10)
mg_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
mg_map = Map(cbind, mg_hs_list, ID = mg_hs_IDs)

hs_mg_df = do.call("rbind", list(mg_map[[1]], mg_map[[2]], mg_map[[3]],
                                 mg_map[[4]], mg_map[[5]], mg_map[[6]],
                                 mg_map[[7]], mg_map[[8]], mg_map[[9]],
                                 mg_map[[10]]))
hs_mg_df$ID = as.factor(hs_mg_df$ID)


rm(list = c("hs_mg1", "hs_mg2", "hs_mg3", "hs_mg4", "hs_mg5", "hs_mg6", "hs_mg7", 
            "hs_mg8","hs_mg9", "hs_mg10", "mg_hs_list", "mg_hs_IDs", "mg_map"))

hs_mg_df$Species = rep("Haemulon sciurus", nrow(hs_mg_df))

hs_mg_df = hs_mg_df %>%
  rename(Variable = variable,
         MG_Dist = x,
         HSI = y,
         Fold = ID) %>%
  filter(MG_Dist >= 0)

hs_mg_summary = hs_mg_df %>%
  group_by(MG_Dist) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>%
  ungroup() %>%
  dplyr::select(Species, Variable, MG_Dist, Mean_HSI, SD_HSI) %>%
  distinct()

# combined mangrove distance results
mg_summary = rbind(lg_mg_summary, hs_mg_summary)

mg_plot = 
  ggplot(mg_summary, aes(y = Mean_HSI, x = MG_Dist, group = Species, color = Species, fill = Species)) +
  geom_point(aes(color = Species)) +
  geom_line(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_ribbon(aes(ymin = Mean_HSI - SD_HSI, ymax = Mean_HSI + SD_HSI),
              alpha = 0.3, show.legend = FALSE) +
  xlab("Distance to Nearest Mangrove (m)") +
  ylab("Predicted Habitat Suitability") + ylim(0, 1) + 
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"))
mg_plot

ggsave(filename = paste0(temp_plots, "Mangrove_Dist_Curve.png"), mg_plot, width = 3.15,
       height = 3.15, units = "in", dpi = 450)

#### DEPTH ####
# gray snapper
lg_depth = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                   pattern = "*Depth_only.csv")

for(i in 1:length(lg_depth)) {                              
  assign(paste0("lg_depth", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                         lg_depth[i])))
  
}

depth_lg_list = list(lg_depth1, lg_depth2, lg_depth3, lg_depth4, lg_depth5, lg_depth6,
                  lg_depth7, lg_depth8, lg_depth9, lg_depth10)
depth_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
depth_map = Map(cbind, depth_lg_list, ID = depth_lg_IDs)

lg_depth_df = do.call("rbind", list(depth_map[[1]], depth_map[[2]], depth_map[[3]],
                                 depth_map[[4]], depth_map[[5]], depth_map[[6]],
                                 depth_map[[7]], depth_map[[8]], depth_map[[9]],
                                 depth_map[[10]]))
lg_depth_df$ID = as.factor(lg_depth_df$ID)


rm(list = c("lg_depth1", "lg_depth2", "lg_depth3", "lg_depth4", "lg_depth5", "lg_depth6", "lg_depth7", 
            "lg_depth8","lg_depth9", "lg_depth10", "depth_lg_list", "depth_lg_IDs", "depth_map"))

lg_depth_df$Species = rep("Lutjanus griseus", nrow(lg_depth_df))

lg_depth_df = lg_depth_df %>%
  rename(Variable = variable,
         Depth = x,
         HSI = y,
         Fold = ID) %>%
  filter(Depth <= 0) %>% # remove extrapolations to 0+ depth (above ground)
  mutate(Depth = Depth*-1) # remove negative sign, assume negative considering this is depth


lg_depth_summary = lg_depth_df %>%
  group_by(Depth) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>%
  ungroup() %>%
  dplyr::select(Species, Variable, Depth, Mean_HSI, SD_HSI) %>%
  distinct()

# bluestriped grunt
hs_depth = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                      pattern = "*Depth_only.csv")

for(i in 1:length(hs_depth)) {                              
  assign(paste0("hs_depth", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_depth[i])))
  
}

depth_hs_list = list(hs_depth1, hs_depth2, hs_depth3, hs_depth4, hs_depth5, hs_depth6,
                     hs_depth7, hs_depth8, hs_depth9, hs_depth10)
depth_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
depth_map = Map(cbind, depth_hs_list, ID = depth_hs_IDs)

hs_depth_df = do.call("rbind", list(depth_map[[1]], depth_map[[2]], depth_map[[3]],
                                    depth_map[[4]], depth_map[[5]], depth_map[[6]],
                                    depth_map[[7]], depth_map[[8]], depth_map[[9]],
                                    depth_map[[10]]))
hs_depth_df$ID = as.factor(hs_depth_df$ID)


rm(list = c("hs_depth1", "hs_depth2", "hs_depth3", "hs_depth4", "hs_depth5", "hs_depth6", "hs_depth7", 
            "hs_depth8","hs_depth9", "hs_depth10", "depth_hs_list", "depth_hs_IDs", "depth_map"))

hs_depth_df$Species = rep("Haemulon sciurus", nrow(hs_depth_df))

hs_depth_df = hs_depth_df %>%
  rename(Variable = variable,
         Depth = x,
         HSI = y,
         Fold = ID) %>%
  filter(Depth <= 0) %>% # remove extrapolations to 0+ depth (above ground)
  mutate(Depth = Depth*-1) # remove negative sign, assume negative considering this is depth

hs_depth_summary = hs_depth_df %>%
  group_by(Depth) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>%
  ungroup() %>%
  dplyr::select(Species, Variable, Depth, Mean_HSI, SD_HSI) %>%
  distinct()

# combined depth results
depth_summary = rbind(lg_depth_summary, hs_depth_summary)

depth_plot = ggplot(depth_summary, aes(y = Mean_HSI, x = Depth, group = Species, color = Species, fill = Species)) +
  geom_point(aes(color = Species)) +
  geom_line(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_ribbon(aes(ymin = Mean_HSI - SD_HSI, ymax = Mean_HSI + SD_HSI),
              alpha = 0.3, show.legend = FALSE) +
  xlab("Depth (m)") +
  ylab("Predicted Habitat Suitability") + ylim(0, 1) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"))
depth_plot

ggsave(filename = paste0(temp_plots, "Depth_Curve.png"), depth_plot, width = 3.15,
       height = 3.15, units = "in", dpi = 450)

#### SLOPE ####
# gray snapper
lg_slope = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                      pattern = "*Slope_only.csv")

for(i in 1:length(lg_slope)) {                              
  assign(paste0("lg_slope", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                         lg_slope[i])))
  
}

slope_lg_list = list(lg_slope1, lg_slope2, lg_slope3, lg_slope4, lg_slope5, lg_slope6,
                     lg_slope7, lg_slope8, lg_slope9, lg_slope10)
slope_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
slope_map = Map(cbind, slope_lg_list, ID = slope_lg_IDs)

lg_slope_df = do.call("rbind", list(slope_map[[1]], slope_map[[2]], slope_map[[3]],
                                    slope_map[[4]], slope_map[[5]], slope_map[[6]],
                                    slope_map[[7]], slope_map[[8]], slope_map[[9]],
                                    slope_map[[10]]))
lg_slope_df$ID = as.factor(lg_slope_df$ID)


rm(list = c("lg_slope1", "lg_slope2", "lg_slope3", "lg_slope4", "lg_slope5", "lg_slope6", "lg_slope7", 
            "lg_slope8","lg_slope9", "lg_slope10", "slope_lg_list", "slope_lg_IDs", "slope_map"))

lg_slope_df$Species = rep("Lutjanus griseus", nrow(lg_slope_df))

lg_slope_df = lg_slope_df %>%
  rename(Variable = variable,
         Slope = x,
         HSI = y,
         Fold = ID) %>%
  filter(Slope >= 0)

lg_slope_summary = lg_slope_df %>%
  group_by(Slope) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>%
  ungroup() %>%
  dplyr::select(Species, Variable, Slope, Mean_HSI, SD_HSI) %>%
  distinct()

# bluestriped grunt
hs_slope = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                      pattern = "*Slope_only.csv")

for(i in 1:length(hs_slope)) {                              
  assign(paste0("hs_slope", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_slope[i])))
  
}

slope_hs_list = list(hs_slope1, hs_slope2, hs_slope3, hs_slope4, hs_slope5, hs_slope6,
                     hs_slope7, hs_slope8, hs_slope9, hs_slope10)
slope_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
slope_map = Map(cbind, slope_hs_list, ID = slope_hs_IDs)

hs_slope_df = do.call("rbind", list(slope_map[[1]], slope_map[[2]], slope_map[[3]],
                                    slope_map[[4]], slope_map[[5]], slope_map[[6]],
                                    slope_map[[7]], slope_map[[8]], slope_map[[9]],
                                    slope_map[[10]]))
hs_slope_df$ID = as.factor(hs_slope_df$ID)


rm(list = c("hs_slope1", "hs_slope2", "hs_slope3", "hs_slope4", "hs_slope5", "hs_slope6", "hs_slope7", 
            "hs_slope8","hs_slope9", "hs_slope10", "slope_hs_list", "slope_hs_IDs", "slope_map"))

hs_slope_df$Species = rep("Haemulon sciurus", nrow(hs_slope_df))

hs_slope_df = hs_slope_df %>%
  rename(Variable = variable,
         Slope = x,
         HSI = y,
         Fold = ID)%>%
  filter(Slope >= 0) 
 
hs_slope_summary = hs_slope_df %>%
  group_by(Slope) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>%
  ungroup() %>%
  dplyr::select(Species, Variable, Slope, Mean_HSI, SD_HSI) %>%
  distinct()

# combined slope results
slope_summary = rbind(lg_slope_summary, hs_slope_summary)

slope_plot = ggplot(slope_summary, aes(y = Mean_HSI, x = Slope, group = Species, color = Species, fill = Species)) +
  geom_point(aes(color = Species)) +
  geom_line(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_ribbon(aes(ymin = Mean_HSI - SD_HSI, ymax = Mean_HSI + SD_HSI),
              alpha = 0.3, show.legend = FALSE) +
  xlab("Slope (Degrees)") +
  ylab("Predicted Habitat Suitability") + ylim(0, 1) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"))
slope_plot


ggsave(filename = paste0(temp_plots, "Slope_Curve.png"), slope_plot, width = 3.15,
       height = 3.15, units = "in", dpi = 450)

#### RUGOSITY ####
# gray snapper
lg_rug = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                      pattern = "*Rugosity_only.csv")

for(i in 1:length(lg_rug)) {                              
  assign(paste0("lg_rug", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                         lg_rug[i])))
  
}

rug_lg_list = list(lg_rug1, lg_rug2, lg_rug3, lg_rug4, lg_rug5, lg_rug6,
                     lg_rug7, lg_rug8, lg_rug9, lg_rug10)
rug_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
rug_map = Map(cbind, rug_lg_list, ID = rug_lg_IDs)

lg_rug_df = do.call("rbind", list(rug_map[[1]], rug_map[[2]], rug_map[[3]],
                                    rug_map[[4]], rug_map[[5]], rug_map[[6]],
                                    rug_map[[7]], rug_map[[8]], rug_map[[9]],
                                    rug_map[[10]]))
lg_rug_df$ID = as.factor(lg_rug_df$ID)


rm(list = c("lg_rug1", "lg_rug2", "lg_rug3", "lg_rug4", "lg_rug5", "lg_rug6", "lg_rug7", 
            "lg_rug8","lg_rug9", "lg_rug10", "rug_lg_list", "rug_lg_IDs", "rug_map"))

lg_rug_df$Species = rep("Lutjanus griseus", nrow(lg_rug_df))

lg_rug_df = lg_rug_df %>%
  rename(Variable = variable,
         Rug = x,
         HSI = y,
         Fold = ID)

lg_rug_summary = lg_rug_df %>%
  group_by(Rug) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% ### NOTE, these are all zero because rugosity was set to 0 by L1 reg
  ungroup() %>%
  dplyr::select(Species, Variable, Rug, Mean_HSI, SD_HSI) %>%
  distinct()

# bluestriped grunt
hs_rug = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                    pattern = "*Rugosity_only.csv")

for(i in 1:length(hs_rug)) {                              
  assign(paste0("hs_rug", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_rug[i])))
  
}

rug_hs_list = list(hs_rug1, hs_rug2, hs_rug3, hs_rug4, hs_rug5, hs_rug6,
                   hs_rug7, hs_rug8, hs_rug9, hs_rug10)
rug_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
rug_map = Map(cbind, rug_hs_list, ID = rug_hs_IDs)

hs_rug_df = do.call("rbind", list(rug_map[[1]], rug_map[[2]], rug_map[[3]],
                                  rug_map[[4]], rug_map[[5]], rug_map[[6]],
                                  rug_map[[7]], rug_map[[8]], rug_map[[9]],
                                  rug_map[[10]]))
hs_rug_df$ID = as.factor(hs_rug_df$ID)


rm(list = c("hs_rug1", "hs_rug2", "hs_rug3", "hs_rug4", "hs_rug5", "hs_rug6", "hs_rug7", 
            "hs_rug8","hs_rug9", "hs_rug10", "rug_hs_list", "rug_hs_IDs", "rug_map"))

hs_rug_df$Species = rep("Haemulon sciurus", nrow(hs_rug_df))

hs_rug_df = hs_rug_df %>%
  rename(Variable = variable,
         Rug = x,
         HSI = y,
         Fold = ID)

hs_rug_summary = hs_rug_df %>%
  group_by(Rug) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% ### NOTE, these are all zero because rugosity was set to 0 by L1 reg
  ungroup() %>%
  dplyr::select(Species, Variable, Rug, Mean_HSI, SD_HSI) %>%
  distinct()

# combined rugosity results
rug_summary = rbind(lg_rug_summary, hs_rug_summary)

rug_plot = ggplot(rug_summary, aes(y = Mean_HSI, x = Rug, group = Species, color = Species, fill = Species)) +
  geom_point(aes(color = Species)) +
  geom_line(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_ribbon(aes(ymin = Mean_HSI - SD_HSI, ymax = Mean_HSI + SD_HSI),
              alpha = 0.3, show.legend = FALSE) +
  xlab("Rugosity (Slope-Corrected SA:PA)") +
  ylab("Predicted Habitat Suitability") + ylim(0, 1) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"))
rug_plot

ggsave(filename = paste0(temp_plots, "Rugosity_Curve.png"), rug_plot, width = 3.15,
       height = 3.15, units = "in", dpi = 450)

#### CURVATURE ####
# gray snapper
lg_curv = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                    pattern = "*Curvature_only.csv")

for(i in 1:length(lg_curv)) {                              
  assign(paste0("lg_curv", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                         lg_curv[i])))
  
}

curv_lg_list = list(lg_curv1, lg_curv2, lg_curv3, lg_curv4, lg_curv5, lg_curv6,
                   lg_curv7, lg_curv8, lg_curv9, lg_curv10)
curv_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
curv_map = Map(cbind, curv_lg_list, ID = curv_lg_IDs)

lg_curv_df = do.call("rbind", list(curv_map[[1]], curv_map[[2]], curv_map[[3]],
                                  curv_map[[4]], curv_map[[5]], curv_map[[6]],
                                  curv_map[[7]], curv_map[[8]], curv_map[[9]],
                                  curv_map[[10]]))
lg_curv_df$ID = as.factor(lg_curv_df$ID)


rm(list = c("lg_curv1", "lg_curv2", "lg_curv3", "lg_curv4", "lg_curv5", "lg_curv6", "lg_curv7", 
            "lg_curv8","lg_curv9", "lg_curv10", "curv_lg_list", "curv_lg_IDs", "curv_map"))

lg_curv_df$Species = rep("Lutjanus griseus", nrow(lg_curv_df))

lg_curv_df = lg_curv_df %>%
  rename(Variable = variable,
         Curvature = x,
         HSI = y,
         Fold = ID)

lg_curv_summary = lg_curv_df %>%
  group_by(Curvature) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, Curvature, Mean_HSI, SD_HSI) %>%
  distinct()

# bluestriped grunt
hs_curv = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                     pattern = "*Curvature_only.csv")

for(i in 1:length(hs_curv)) {                              
  assign(paste0("hs_curv", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_curv[i])))
  
}

curv_hs_list = list(hs_curv1, hs_curv2, hs_curv3, hs_curv4, hs_curv5, hs_curv6,
                    hs_curv7, hs_curv8, hs_curv9, hs_curv10)
curv_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
curv_map = Map(cbind, curv_hs_list, ID = curv_hs_IDs)

hs_curv_df = do.call("rbind", list(curv_map[[1]], curv_map[[2]], curv_map[[3]],
                                   curv_map[[4]], curv_map[[5]], curv_map[[6]],
                                   curv_map[[7]], curv_map[[8]], curv_map[[9]],
                                   curv_map[[10]]))
hs_curv_df$ID = as.factor(hs_curv_df$ID)


rm(list = c("hs_curv1", "hs_curv2", "hs_curv3", "hs_curv4", "hs_curv5", "hs_curv6", "hs_curv7", 
            "hs_curv8","hs_curv9", "hs_curv10", "curv_hs_list", "curv_hs_IDs", "curv_map"))

hs_curv_df$Species = rep("Haemulon sciurus", nrow(hs_curv_df))

hs_curv_df = hs_curv_df %>%
  rename(Variable = variable,
         Curvature = x,
         HSI = y,
         Fold = ID)

hs_curv_summary = hs_curv_df %>%
  group_by(Curvature) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, Curvature, Mean_HSI, SD_HSI) %>%
  distinct()

# combined curvature results
curv_summary = rbind(lg_curv_summary, hs_curv_summary)

curv_plot = ggplot(curv_summary, aes(y = Mean_HSI, x = Curvature, group = Species, color = Species, fill = Species)) +
  geom_point(aes(color = Species)) +
  geom_line(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_ribbon(aes(ymin = Mean_HSI - SD_HSI, ymax = Mean_HSI + SD_HSI),
              alpha = 0.3, show.legend = FALSE) +
  xlab("Curvature") +
  ylab("Predicted Habitat Suitability") + ylim(0, 1) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"))
curv_plot

ggsave(filename = paste0(temp_plots, "Curvature_Curve.png"), curv_plot, width = 3.15,
       height = 3.15, units = "in", dpi = 450)

#### BPI FINE ####
# gray snapper
lg_bpiF = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                     pattern = "*BPI_Fine_only.csv")

for(i in 1:length(lg_bpiF)) {                              
  assign(paste0("lg_bpiF", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                         lg_bpiF[i])))
  
}

bpiF_lg_list = list(lg_bpiF1, lg_bpiF2, lg_bpiF3, lg_bpiF4, lg_bpiF5, lg_bpiF6,
                    lg_bpiF7, lg_bpiF8, lg_bpiF9, lg_bpiF10)
bpiF_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
bpiF_map = Map(cbind, bpiF_lg_list, ID = bpiF_lg_IDs)

lg_bpiF_df = do.call("rbind", list(bpiF_map[[1]], bpiF_map[[2]], bpiF_map[[3]],
                                   bpiF_map[[4]], bpiF_map[[5]], bpiF_map[[6]],
                                   bpiF_map[[7]], bpiF_map[[8]], bpiF_map[[9]],
                                   bpiF_map[[10]]))
lg_bpiF_df$ID = as.factor(lg_bpiF_df$ID)


rm(list = c("lg_bpiF1", "lg_bpiF2", "lg_bpiF3", "lg_bpiF4", "lg_bpiF5", "lg_bpiF6", "lg_bpiF7", 
            "lg_bpiF8","lg_bpiF9", "lg_bpiF10", "bpiF_lg_list", "bpiF_lg_IDs", "bpiF_map"))

lg_bpiF_df$Species = rep("Lutjanus griseus", nrow(lg_bpiF_df))

lg_bpiF_df = lg_bpiF_df %>%
  rename(Variable = variable,
         BPI_Fine = x,
         HSI = y,
         Fold = ID)

lg_bpiF_summary = lg_bpiF_df %>%
  group_by(BPI_Fine) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, BPI_Fine, Mean_HSI, SD_HSI) %>%
  distinct()

# bluestriped grunt
hs_bpiF = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                     pattern = "*BPI_Fine_only.csv")

for(i in 1:length(hs_bpiF)) {                              
  assign(paste0("hs_bpiF", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_bpiF[i])))
  
}

bpiF_hs_list = list(hs_bpiF1, hs_bpiF2, hs_bpiF3, hs_bpiF4, hs_bpiF5, hs_bpiF6,
                    hs_bpiF7, hs_bpiF8, hs_bpiF9, hs_bpiF10)
bpiF_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
bpiF_map = Map(cbind, bpiF_hs_list, ID = bpiF_hs_IDs)

hs_bpiF_df = do.call("rbind", list(bpiF_map[[1]], bpiF_map[[2]], bpiF_map[[3]],
                                   bpiF_map[[4]], bpiF_map[[5]], bpiF_map[[6]],
                                   bpiF_map[[7]], bpiF_map[[8]], bpiF_map[[9]],
                                   bpiF_map[[10]]))
hs_bpiF_df$ID = as.factor(hs_bpiF_df$ID)


rm(list = c("hs_bpiF1", "hs_bpiF2", "hs_bpiF3", "hs_bpiF4", "hs_bpiF5", "hs_bpiF6", "hs_bpiF7", 
            "hs_bpiF8","hs_bpiF9", "hs_bpiF10", "bpiF_hs_list", "bpiF_hs_IDs", "bpiF_map"))

hs_bpiF_df$Species = rep("Haemulon sciurus", nrow(hs_bpiF_df))

hs_bpiF_df = hs_bpiF_df %>%
  rename(Variable = variable,
         BPI_Fine = x,
         HSI = y,
         Fold = ID)

hs_bpiF_summary = hs_bpiF_df %>%
  group_by(BPI_Fine) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, BPI_Fine, Mean_HSI, SD_HSI) %>%
  distinct()

bpiF_summary = rbind(lg_bpiF_summary, hs_bpiF_summary)

bpiF_plot = ggplot(bpiF_summary, aes(y = Mean_HSI, x = BPI_Fine, group = Species, color = Species, fill = Species)) +
  geom_point(aes(color = Species)) +
  geom_line(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_ribbon(aes(ymin = Mean_HSI - SD_HSI, ymax = Mean_HSI + SD_HSI),
              alpha = 0.3, show.legend = FALSE) +
  xlab("Fine-Scale BPI") +
  ylab("Predicted Habitat Suitability") + ylim(0, 1) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"))
bpiF_plot

ggsave(filename = paste0(temp_plots, "BPI_Fine_Curve.png"), bpiF_plot, width = 3.15,
       height = 3.15, units = "in", dpi = 450)

#### BPI BROAD ####
# gray snapper
lg_bpiB = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                     pattern = "*BPI_Broad_only.csv")

for(i in 1:length(lg_bpiB)) {                              
  assign(paste0("lg_bpiB", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                         lg_bpiB[i])))
  
}

bpiB_lg_list = list(lg_bpiB1, lg_bpiB2, lg_bpiB3, lg_bpiB4, lg_bpiB5, lg_bpiB6,
                    lg_bpiB7, lg_bpiB8, lg_bpiB9, lg_bpiB10)
bpiB_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
bpiB_map = Map(cbind, bpiB_lg_list, ID = bpiB_lg_IDs)

lg_bpiB_df = do.call("rbind", list(bpiB_map[[1]], bpiB_map[[2]], bpiB_map[[3]],
                                   bpiB_map[[4]], bpiB_map[[5]], bpiB_map[[6]],
                                   bpiB_map[[7]], bpiB_map[[8]], bpiB_map[[9]],
                                   bpiB_map[[10]]))
lg_bpiB_df$ID = as.factor(lg_bpiB_df$ID)


rm(list = c("lg_bpiB1", "lg_bpiB2", "lg_bpiB3", "lg_bpiB4", "lg_bpiB5", "lg_bpiB6", "lg_bpiB7", 
            "lg_bpiB8","lg_bpiB9", "lg_bpiB10", "bpiB_lg_list", "bpiB_lg_IDs", "bpiB_map"))

lg_bpiB_df$Species = rep("Lutjanus griseus", nrow(lg_bpiB_df))

lg_bpiB_df = lg_bpiB_df %>%
  rename(Variable = variable,
         BPI_Broad = x,
         HSI = y,
         Fold = ID)

lg_bpiB_summary = lg_bpiB_df %>%
  group_by(BPI_Broad) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, BPI_Broad, Mean_HSI, SD_HSI) %>%
  distinct()

# bluestriped grunt
hs_bpiB = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                     pattern = "*BPI_Broad_only.csv")

for(i in 1:length(hs_bpiB)) {                              
  assign(paste0("hs_bpiB", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_bpiB[i])))
  
}

bpiB_hs_list = list(hs_bpiB1, hs_bpiB2, hs_bpiB3, hs_bpiB4, hs_bpiB5, hs_bpiB6,
                    hs_bpiB7, hs_bpiB8, hs_bpiB9, hs_bpiB10)
bpiB_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
bpiB_map = Map(cbind, bpiB_hs_list, ID = bpiB_hs_IDs)

hs_bpiB_df = do.call("rbind", list(bpiB_map[[1]], bpiB_map[[2]], bpiB_map[[3]],
                                   bpiB_map[[4]], bpiB_map[[5]], bpiB_map[[6]],
                                   bpiB_map[[7]], bpiB_map[[8]], bpiB_map[[9]],
                                   bpiB_map[[10]]))
hs_bpiB_df$ID = as.factor(hs_bpiB_df$ID)


rm(list = c("hs_bpiB1", "hs_bpiB2", "hs_bpiB3", "hs_bpiB4", "hs_bpiB5", "hs_bpiB6", "hs_bpiB7", 
            "hs_bpiB8","hs_bpiB9", "hs_bpiB10", "bpiB_hs_list", "bpiB_hs_IDs", "bpiB_map"))

hs_bpiB_df$Species = rep("Haemulon sciurus", nrow(hs_bpiB_df))

hs_bpiB_df = hs_bpiB_df %>%
  rename(Variable = variable,
         BPI_Broad = x,
         HSI = y,
         Fold = ID)

hs_bpiB_summary = hs_bpiB_df %>%
  group_by(BPI_Broad) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, BPI_Broad, Mean_HSI, SD_HSI) %>%
  distinct()

bpiB_summary = rbind(lg_bpiB_summary, hs_bpiB_summary)

bpiB_plot = ggplot(bpiB_summary, aes(y = Mean_HSI, x = BPI_Broad, group = Species, color = Species, fill = Species)) +
  geom_point(aes(color = Species)) +
  geom_line(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_ribbon(aes(ymin = Mean_HSI - SD_HSI, ymax = Mean_HSI + SD_HSI),
              alpha = 0.3, show.legend = FALSE) +
  xlab("Broad-Scale BPI") +
  ylab("Predicted Habitat Suitability") + ylim(0, 1) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"))
bpiB_plot

ggsave(filename = paste0(temp_plots, "BPI_Broad_Curve.png"), bpiB_plot, width = 3.15,
       height = 3.15, units = "in", dpi = 450)

#### WINTER SALINITY ####
# gray snapper
lg_wsal = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                     pattern = "*Win_Sal_only.csv")

for(i in 1:length(lg_wsal)) {                              
  assign(paste0("lg_wsal", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                         lg_wsal[i])))
  
}

wsal_lg_list = list(lg_wsal1, lg_wsal2, lg_wsal3, lg_wsal4, lg_wsal5, lg_wsal6,
                    lg_wsal7, lg_wsal8, lg_wsal9, lg_wsal10)
wsal_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
wsal_map = Map(cbind, wsal_lg_list, ID = wsal_lg_IDs)

lg_wsal_df = do.call("rbind", list(wsal_map[[1]], wsal_map[[2]], wsal_map[[3]],
                                   wsal_map[[4]], wsal_map[[5]], wsal_map[[6]],
                                   wsal_map[[7]], wsal_map[[8]], wsal_map[[9]],
                                   wsal_map[[10]]))
lg_wsal_df$ID = as.factor(lg_wsal_df$ID)


rm(list = c("lg_wsal1", "lg_wsal2", "lg_wsal3", "lg_wsal4", "lg_wsal5", "lg_wsal6", "lg_wsal7", 
            "lg_wsal8","lg_wsal9", "lg_wsal10", "wsal_lg_list", "wsal_lg_IDs", "wsal_map"))

lg_wsal_df$Species = rep("Lutjanus griseus", nrow(lg_wsal_df))

lg_wsal_df = lg_wsal_df %>%
  rename(Variable = variable,
         Win_Sal = x,
         HSI = y,
         Fold = ID)

lg_wsal_summary = lg_wsal_df %>%
  group_by(Win_Sal) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, Win_Sal, Mean_HSI, SD_HSI) %>%
  distinct()

# bluestriped grunt
hs_wsal = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                     pattern = "*Win_Sal_only.csv")

for(i in 1:length(hs_wsal)) {                              
  assign(paste0("hs_wsal", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_wsal[i])))
  
}

wsal_hs_list = list(hs_wsal1, hs_wsal2, hs_wsal3, hs_wsal4, hs_wsal5, hs_wsal6,
                    hs_wsal7, hs_wsal8, hs_wsal9, hs_wsal10)
wsal_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
wsal_map = Map(cbind, wsal_hs_list, ID = wsal_hs_IDs)

hs_wsal_df = do.call("rbind", list(wsal_map[[1]], wsal_map[[2]], wsal_map[[3]],
                                   wsal_map[[4]], wsal_map[[5]], wsal_map[[6]],
                                   wsal_map[[7]], wsal_map[[8]], wsal_map[[9]],
                                   wsal_map[[10]]))
hs_wsal_df$ID = as.factor(hs_wsal_df$ID)


rm(list = c("hs_wsal1", "hs_wsal2", "hs_wsal3", "hs_wsal4", "hs_wsal5", "hs_wsal6", "hs_wsal7", 
            "hs_wsal8","hs_wsal9", "hs_wsal10", "wsal_hs_list", "wsal_hs_IDs", "wsal_map"))

hs_wsal_df$Species = rep("Haemulon sciurus", nrow(hs_wsal_df))

hs_wsal_df = hs_wsal_df %>%
  rename(Variable = variable,
         Win_Sal = x,
         HSI = y,
         Fold = ID)

hs_wsal_summary = hs_wsal_df %>%
  group_by(Win_Sal) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, Win_Sal, Mean_HSI, SD_HSI) %>%
  distinct()

wsal_summary = rbind(lg_wsal_summary, hs_wsal_summary)

wsal_plot = ggplot(wsal_summary, aes(y = Mean_HSI, x = Win_Sal, group = Species, color = Species, fill = Species)) +
  geom_point(aes(color = Species)) +
  geom_line(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_ribbon(aes(ymin = Mean_HSI - SD_HSI, ymax = Mean_HSI + SD_HSI),
              alpha = 0.3, show.legend = FALSE) +
  xlab("Winter Salinity (PSU)") +
  ylab("Predicted Habitat Suitability") + ylim(0, 1) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"))
wsal_plot

ggsave(filename = paste0(temp_plots, "Win_Sal_Curve.png"), wsal_plot, width = 3.15,
       height = 3.15, units = "in", dpi = 450)


#### SUMMER SALINITY ####
# gray snapper
lg_ssal = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                     pattern = "*Sum_Sal_only.csv")

for(i in 1:length(lg_ssal)) {                              
  assign(paste0("lg_ssal", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                         lg_ssal[i])))
  
}

ssal_lg_list = list(lg_ssal1, lg_ssal2, lg_ssal3, lg_ssal4, lg_ssal5, lg_ssal6,
                    lg_ssal7, lg_ssal8, lg_ssal9, lg_ssal10)
ssal_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
ssal_map = Map(cbind, ssal_lg_list, ID = ssal_lg_IDs)

lg_ssal_df = do.call("rbind", list(ssal_map[[1]], ssal_map[[2]], ssal_map[[3]],
                                   ssal_map[[4]], ssal_map[[5]], ssal_map[[6]],
                                   ssal_map[[7]], ssal_map[[8]], ssal_map[[9]],
                                   ssal_map[[10]]))
lg_ssal_df$ID = as.factor(lg_ssal_df$ID)


rm(list = c("lg_ssal1", "lg_ssal2", "lg_ssal3", "lg_ssal4", "lg_ssal5", "lg_ssal6", "lg_ssal7", 
            "lg_ssal8","lg_ssal9", "lg_ssal10", "ssal_lg_list", "ssal_lg_IDs", "ssal_map"))

lg_ssal_df$Species = rep("Lutjanus griseus", nrow(lg_ssal_df))

lg_ssal_df = lg_ssal_df %>%
  rename(Variable = variable,
         Sum_Sal = x,
         HSI = y,
         Fold = ID)

lg_ssal_summary = lg_ssal_df %>%
  group_by(Sum_Sal) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, Sum_Sal, Mean_HSI, SD_HSI) %>%
  distinct()

# bluestriped grunt
hs_ssal = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                     pattern = "*Sum_Sal_only.csv")

for(i in 1:length(hs_ssal)) {                              
  assign(paste0("hs_ssal", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_ssal[i])))
  
}

ssal_hs_list = list(hs_ssal1, hs_ssal2, hs_ssal3, hs_ssal4, hs_ssal5, hs_ssal6,
                    hs_ssal7, hs_ssal8, hs_ssal9, hs_ssal10)
ssal_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
ssal_map = Map(cbind, ssal_hs_list, ID = ssal_hs_IDs)

hs_ssal_df = do.call("rbind", list(ssal_map[[1]], ssal_map[[2]], ssal_map[[3]],
                                   ssal_map[[4]], ssal_map[[5]], ssal_map[[6]],
                                   ssal_map[[7]], ssal_map[[8]], ssal_map[[9]],
                                   ssal_map[[10]]))
hs_ssal_df$ID = as.factor(hs_ssal_df$ID)


rm(list = c("hs_ssal1", "hs_ssal2", "hs_ssal3", "hs_ssal4", "hs_ssal5", "hs_ssal6", "hs_ssal7", 
            "hs_ssal8","hs_ssal9", "hs_ssal10", "ssal_hs_list", "ssal_hs_IDs", "ssal_map"))

hs_ssal_df$Species = rep("Haemulon sciurus", nrow(hs_ssal_df))

hs_ssal_df = hs_ssal_df %>%
  rename(Variable = variable,
         Sum_Sal = x,
         HSI = y,
         Fold = ID)

hs_ssal_summary = hs_ssal_df %>%
  group_by(Sum_Sal) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, Sum_Sal, Mean_HSI, SD_HSI) %>%
  distinct()

ssal_summary = rbind(lg_ssal_summary, hs_ssal_summary)

ssal_plot = ggplot(ssal_summary, aes(y = Mean_HSI, x = Sum_Sal, group = Species, color = Species, fill = Species)) +
  geom_point(aes(color = Species)) +
  geom_line(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_ribbon(aes(ymin = Mean_HSI - SD_HSI, ymax = Mean_HSI + SD_HSI),
              alpha = 0.3, show.legend = FALSE) +
  xlab("Summer Salinity (PSU)") +
  ylab("Predicted Habitat Suitability") + ylim(0, 1) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"))
ssal_plot

ggsave(filename = paste0(temp_plots, "Sum_Sal_Curve.png"), ssal_plot, width = 3.15,
       height = 3.15, units = "in", dpi = 450)

#### SUMMER TEMPERATURE ####
# gray snapper
lg_stemp = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                     pattern = "*Sum_Temp_only.csv")

for(i in 1:length(lg_stemp)) {                              
  assign(paste0("lg_stemp", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/plots/",
                         lg_stemp[i])))
  
}

stemp_lg_list = list(lg_stemp1, lg_stemp2, lg_stemp3, lg_stemp4, lg_stemp5, lg_stemp6,
                    lg_stemp7, lg_stemp8, lg_stemp9, lg_stemp10)
stemp_lg_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
stemp_map = Map(cbind, stemp_lg_list, ID = stemp_lg_IDs)

lg_stemp_df = do.call("rbind", list(stemp_map[[1]], stemp_map[[2]], stemp_map[[3]],
                                   stemp_map[[4]], stemp_map[[5]], stemp_map[[6]],
                                   stemp_map[[7]], stemp_map[[8]], stemp_map[[9]],
                                   stemp_map[[10]]))
lg_stemp_df$ID = as.factor(lg_stemp_df$ID)


rm(list = c("lg_stemp1", "lg_stemp2", "lg_stemp3", "lg_stemp4", "lg_stemp5", "lg_stemp6", "lg_stemp7", 
            "lg_stemp8","lg_stemp9", "lg_stemp10", "stemp_lg_list", "stemp_lg_IDs", "stemp_map"))

lg_stemp_df$Species = rep("Lutjanus griseus", nrow(lg_stemp_df))

lg_stemp_df = lg_stemp_df %>%
  rename(Variable = variable,
         Sum_Temp = x,
         HSI = y,
         Fold = ID)

lg_stemp_summary = lg_stemp_df %>%
  group_by(Sum_Temp) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, Sum_Temp, Mean_HSI, SD_HSI) %>%
  distinct()

# bluestriped grunt
hs_stemp = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                      pattern = "*Sum_Temp_only.csv")

for(i in 1:length(hs_stemp)) {                              
  assign(paste0("hs_stemp", i),                                   
         read.csv(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/plots/",
                         hs_stemp[i])))
  
}

stemp_hs_list = list(hs_stemp1, hs_stemp2, hs_stemp3, hs_stemp4, hs_stemp5, hs_stemp6,
                     hs_stemp7, hs_stemp8, hs_stemp9, hs_stemp10)
stemp_hs_IDs = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
stemp_map = Map(cbind, stemp_hs_list, ID = stemp_hs_IDs)

hs_stemp_df = do.call("rbind", list(stemp_map[[1]], stemp_map[[2]], stemp_map[[3]],
                                    stemp_map[[4]], stemp_map[[5]], stemp_map[[6]],
                                    stemp_map[[7]], stemp_map[[8]], stemp_map[[9]],
                                    stemp_map[[10]]))
hs_stemp_df$ID = as.factor(hs_stemp_df$ID)


rm(list = c("hs_stemp1", "hs_stemp2", "hs_stemp3", "hs_stemp4", "hs_stemp5", "hs_stemp6", "hs_stemp7", 
            "hs_stemp8","hs_stemp9", "hs_stemp10", "stemp_hs_list", "stemp_hs_IDs", "stemp_map"))

hs_stemp_df$Species = rep("Haemulon sciurus", nrow(hs_stemp_df))

hs_stemp_df = hs_stemp_df %>%
  rename(Variable = variable,
         Sum_Temp = x,
         HSI = y,
         Fold = ID)

hs_stemp_summary = hs_stemp_df %>%
  group_by(Sum_Temp) %>%
  mutate(Mean_HSI = mean(HSI),
         SD_HSI = sd(HSI)) %>% 
  ungroup() %>%
  dplyr::select(Species, Variable, Sum_Temp, Mean_HSI, SD_HSI) %>%
  distinct()

stemp_summary = rbind(lg_stemp_summary, hs_stemp_summary)

stemp_plot = ggplot(stemp_summary, aes(y = Mean_HSI, x = Sum_Temp, group = Species, color = Species, fill = Species)) +
  geom_point(aes(color = Species)) +
  geom_line(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  geom_ribbon(aes(ymin = Mean_HSI - SD_HSI, ymax = Mean_HSI + SD_HSI),
              alpha = 0.3, show.legend = FALSE) +
  xlab("Summer Temperature (C)") +
  ylab("Predicted Habitat Suitability") + ylim(0, 1) + xlim(29.65, 31) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.text = element_text(face = "italic"))
stemp_plot

ggsave(filename = paste0(temp_plots, "Sum_Temp_Curve.png"), stemp_plot, width = 3.15,
       height = 3.15, units = "in", dpi = 450)

#### VARIABLE CONTRIBUTIONS ####
lg_results = read.csv(paste0(HSMs, "MaxEnt/Subadult_Gray_Snapper/maxentResults.csv"))
lg_results = lg_results[11,] %>% # keep only average from 10 CV runs
  mutate(Species = "Lutjanus griseus")

hs_results = read.csv(paste0(HSMs, "MaxEnt/Subadult_Bluestriped_Grunt/maxentResults.csv"))
hs_results = hs_results[11,] %>%
  mutate(Species = "Haemulon sciurus")

results = rbind(lg_results, hs_results)

perc = results %>%
  dplyr::select(1, 12:22) %>%
  rename(Habitat = Habitat.contribution, Mangrove_Dist = Mangrove_Dist.contribution,
         Depth = Depth.contribution, Slope = Slope.contribution,
         Curvature = Curvature.contribution, Rugosity = Rugosity.contribution,
         BPI_Fine = BPI_Fine.contribution, BPI_Broad = BPI_Broad.contribution,
         Summer_Salinity = Mean_Sum_Sal.contribution, Summer_Temperature = Mean_Sum_Temp.contribution,
         Winter_Salinity = Mean_Win_Sal.contribution)

perc = gather(perc, Variable, Contribution, 2:12, factor_key = TRUE)
perc = cbind(perc, data.frame("Variable2" = c("BPI Broad", "BPI Broad", "BPI Fine",
                                          "BPI Fine", "Curvature", "Curvature",
                                          "Depth", "Depth", "Habitat", "Habitat",
                                          "Mangrove Distance", "Mangrove Distance",
                                          "Summer Salinity", "Summer Salinity",
                                          "Summer Temperature", "Summer Temperature",
                                          "Winter Salinity", "Winter Salinity",
                                          "Rugosity", "Rugosity", "Slope", "Slope"))) %>%
  select(Contribution, Variable2, Species) %>%
  rename(Overall = Contribution, Variable = Variable2)

# save variable importance info
write.csv(perc, paste0(csv_wd, "MaxEnt_Variable_Importance.csv"))

perc_cont = 
  ggplot(perc, aes(x = Variable, y = Overall, fill = Species)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  xlab(" ") + ylab("MaxEnt Variable Importance (%)") +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  scale_x_discrete(limits = c("Winter Salinity", "Summer Salinity", "Summer Temperature", 
                              "BPI Broad", "BPI Fine", "Rugosity", "Curvature", "Slope",
                              "Depth", "Mangrove Distance", "Habitat")) +
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
  coord_flip() + theme(plot.margin = margin(0, 1, 0, 0, "cm"))

perc_cont
ggsave(filename = paste0(temp_plots, "MaxEnt_Variable_Importance.png"), perc_cont, width = 8.5,
       height = 3.5, units = "in", dpi = 400)

# permutation importance
perm = results %>%
  dplyr::select(1, 23:33) %>%
  rename(Habitat = Habitat.permutation.importance, Mangrove_Dist = Mangrove_Dist.permutation.importance,
         Depth = Depth.permutation.importance, Slope = Slope.permutation.importance,
         Curvature = Curvature.permutation.importance, Rugosity = Rugosity.permutation.importance,
         BPI_Fine = BPI_Fine.permutation.importance, BPI_Broad = BPI_Broad.permutation.importance,
         Summer_Salinity = Mean_Sum_Sal.permutation.importance, Summer_Temperature = Mean_Sum_Temp.permutation.importance,
         Winter_Salinity = Mean_Win_Sal.permutation.importance)

perm = gather(perm, Variable, Contribution, 2:12, factor_key = TRUE)

perm_imp = 
  ggplot(perm, aes(x = Variable, y = Contribution, fill = Species)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  xlab(" ") + ylab("Permutation Importance (%)") +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(color = "black"),
                     axis.text.y = element_text(color = "black"),
                     axis.line = element_line(color = "black"),
                     legend.text = element_text(face = "italic")) +
  geom_text(aes(label = round(Contribution, digits = 2)), size = 2.5, 
            position = position_dodge(width = 1), hjust = -0.25, vjust = 0.35) + coord_flip()

perm_imp

ggsave(filename = paste0(temp_plots, "Permutation_Importance.png"), perm_imp, width = 8.5,
       height = 5, units = "in", dpi = 400)

#### SPECIES OVERLAP ####
# overlap on a scale of 0 (no overlap) to 1 (identical)
lg_mosaic = mosaic(lg_maxent_train, lg_maxent_test, fun = mean)
hs_mosaic = mosaic(hs_maxent_train, hs_maxent_test, fun = mean)
overlap = nicheOverlap(lg_mosaic, hs_mosaic, stat = 'I', 
                       mask = FALSE, checkNegatives = FALSE)

# PARSE LAMBDAS ###
library(devtools)
# install_github('johnbaums/rmaxent')
library(rmaxent)

lg_lambdas = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/",
                      pattern = "*.lambdas")

for(i in 1:length(lg_lambdas)) {                              
  assign(paste0("lg_lambdas", i),                                   
         parse_lambdas(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Gray_Snapper/",
                         lg_lambdas[i])))
  
}

lg_lambdas_linear = as.data.frame(do.call("rbind", list(as.data.frame(lg_lambdas1$lambdas), 
                                                        as.data.frame(lg_lambdas2$lambdas), 
                                                        as.data.frame(lg_lambdas3$lambdas), 
                                                        as.data.frame(lg_lambdas4$lambdas), 
                                                        as.data.frame(lg_lambdas5$lambdas), 
                                                        as.data.frame(lg_lambdas6$lambdas), 
                                                        as.data.frame(lg_lambdas7$lambdas),
                                                        as.data.frame(lg_lambdas8$lambdas), 
                                                        as.data.frame(lg_lambdas9$lambdas), 
                                                        as.data.frame(lg_lambdas10$lambdas)))) 
lg_lambdas_linear = lg_lambdas_linear %>% 
  filter(type == "categorical" | type == "linear") %>%
  group_by(feature) %>%
  mutate(Avg_Coef = mean(lambda)) %>%
  distinct(feature, var, Avg_Coef) %>%
  ungroup()
lg_lambdas_linear$Species = rep("Lutjanus griseus", nrow(lg_lambdas_linear))
lg_vars = as.data.frame(c("Individual/Aggregate Patch Reef", "Discontinuous Seagrass", 
                          "Unconsolidated Sediment", "Aggregate Reef", "Mangrove", 
                          "BPI Broad", "BPI Fine", "Curvature", "Depth", "Mangrove Distance",
                          "Rugosity", "Slope", "Summer Salinity", "Summer Temperature",
                          "Winter Salinity", "Pavement"))
colnames(lg_vars) = "Variable"
lg_lambdas_linear = cbind(lg_lambdas_linear, lg_vars)

hs_lambdas = list.files("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/",
                        pattern = "*.lambdas")

for(i in 1:length(hs_lambdas)) {                              
  assign(paste0("hs_lambdas", i),                                   
         parse_lambdas(paste0("Z:/Courtney/Stuart_MSc_Ch1/HSMs/MaxEnt/Subadult_Bluestriped_Grunt/",
                              hs_lambdas[i])))
  
}

hs_lambdas_linear = as.data.frame(do.call("rbind", list(as.data.frame(hs_lambdas1$lambdas), 
                                                        as.data.frame(hs_lambdas2$lambdas), 
                                                        as.data.frame(hs_lambdas3$lambdas), 
                                                        as.data.frame(hs_lambdas4$lambdas), 
                                                        as.data.frame(hs_lambdas5$lambdas), 
                                                        as.data.frame(hs_lambdas6$lambdas), 
                                                        as.data.frame(hs_lambdas7$lambdas),
                                                        as.data.frame(hs_lambdas8$lambdas), 
                                                        as.data.frame(hs_lambdas9$lambdas), 
                                                        as.data.frame(hs_lambdas10$lambdas)))) 

hs_lambdas_linear = hs_lambdas_linear %>% 
  filter(type == "categorical" | type == "linear") %>%
  group_by(feature) %>%
  mutate(Avg_Coef = mean(lambda))%>%
  distinct(feature, var, Avg_Coef) %>%
  ungroup()
hs_lambdas_linear$Species = rep("Haemulon sciurus", nrow(hs_lambdas_linear))
hs_vars = as.data.frame(c("Individual/Aggregate Patch Reef", "Continuous Seagrass", 
                          "Discontinuous Seagrass", "Unconsolidated Sediment", 
                          "Reef Rubble", "BPI Broad", "BPI Fine", "Curvature", "Depth", 
                          "Mangrove Distance", "Rugosity", "Slope", "Summer Salinity",
                          "Summer Temperature", "Winter Salinity", "Aggregate Reef"))
colnames(hs_vars) = "Variable"
hs_lambdas_linear = cbind(hs_lambdas_linear, hs_vars)

lambdas_linear = rbind(lg_lambdas_linear, hs_lambdas_linear)

# the original depth data are negative, meaning larger negative numbers (deeper water)
# are technically "smaller" than less negative (shallower) values. This is a bit confusing 
# when looking at coefficients, because technically a positive depth coefficient means that
# habitat suitability increases as depth values get larger (shallower). Instead, flip
# the sign of the coefficient so that they represent the change in suitability at increasing depths
lambdas_linear[9, 3] = as.numeric(-1*lambdas_linear[9, 3])
lambdas_linear[25, 3] = as.numeric(-1*lambdas_linear[25, 3])
levels(lambdas_linear$Variable)[levels(lambdas_linear$Variable) =="Individual/Aggregate Patch Reef"] = "Patch Reef"

lambdas_plot = ggplot(data = lambdas_linear, aes(x = Avg_Coef, y = Variable, fill = Species)) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(color = Species)) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  labs(y = "",
       x = "MaxEnt coefficients") + xlim(-7.5, 7.5) +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_rect(color = "black"), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), legend.text = element_text(face = "italic", size = 8),
                     legend.margin = margin(0,0,0,0), plot.margin = margin(c(0,0,0,0)),
                     axis.text = element_text(size = 8), axis.title.x = element_text(size = 10),
                     legend.box.margin=margin(t = 0, r = 0, b = -10, l = 0)) 
lambdas_plot 

lambdas_plot = lambdas_plot + scale_y_discrete(
  labels = str_wrap(c("Winter Salinity", "Summer Salinity", "Summer Temperature", "BPI Broad", "BPI Fine",
                      "Curvature", "Rugosity", "Slope", "Depth", "Mangrove Distance", "Discontinuous Seagrass", 
                      "Continuous Seagrass", "Mangrove", "Pavement", "Aggregate Reef","Reef Rubble", 
                      "Patch Reef", "Unconsolidated Sediment"), width = 25))
lambdas_plot

ggsave(plot = lambdas_plot, filename = paste0(temp_plots, "MaxEnt_Coefficients.png"), height = 2.6, 
       width = 5, units = "in", dpi = 450)


