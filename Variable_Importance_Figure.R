# Script for plotting the top five most important variables for each model type

#### SET-UP ####
# data directories
setwd("Z:/Courtney/Stuart_MSc_Ch1/") # main project folder
temp_wd = "Z:/Courtney/Stuart_MSc_Ch1/Temp/" # temporary files
csv_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Data/" # for writing tabular data
plots_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Figures/" # for figures
temp_plots = "Z:/Courtney/Stuart_MSc_Ch1/Plots/" # temporary plots for Courtney

# libraries
library(easypackages)
libraries("tidyverse", "tidyr", "dplyr", "ggplot2",  "PNWColors", "cowplot", "conflicted")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# palette
my_pal = pnw_palette("Bay",8)
my_pal

# variable importance csv files
lasso_vimp = read.csv(paste0(csv_wd, "Lasso_Variable_Importance.csv"))
ridge_vimp = read.csv(paste0(csv_wd, "Ridge_Variable_Importance.csv"))
maxent_vimp = read.csv(paste0(csv_wd, "MaxEnt_Variable_Importance.csv"))

# top 5 variables from each model type
# lasso
lg_lasso_top5 = lasso_vimp %>%
  select(Overall, Variable, Species) %>%
  filter(Species == "Lutjanus griseus" & 
           Variable %in% c("Depth", "Mangrove Distance", "Summer Temperature",
                           "Reef Rubble", "Pavement"))
hs_lasso_top5 = lasso_vimp %>%
  select(Overall, Variable, Species) %>%
  filter(Species == "Haemulon sciurus" & 
           Variable %in% c("Winter Salinity", "Depth", "Pavement",
                           "Discontinuous Seagrass", "BPI Broad"))

lasso_top5 = rbind(lg_lasso_top5, hs_lasso_top5)

# ridge
lg_ridge_top5 = ridge_vimp %>%
  select(Overall, Variable, Species) %>%
  filter(Species == "Lutjanus griseus" & 
           Variable %in% c("Depth", "Mangrove Distance", "Summer Temperature",
                           "Reef Rubble", "Pavement"))
hs_ridge_top5 = ridge_vimp %>%
  select(Overall, Variable, Species) %>%
  filter(Species == "Haemulon sciurus" & 
           Variable %in% c("Winter Salinity", "Depth", "Summer Salinity",
                           "Pavement", "Discontinuous Seagrass"))

ridge_top5 = rbind(lg_ridge_top5, hs_ridge_top5)

# maxent
lg_maxent_top5 = maxent_vimp %>%
  select(Overall, Variable, Species) %>%
  filter(Species == "Lutjanus griseus" & 
           Variable %in% c("Habitat", "Mangrove Distance", "Depth",
                           "Slope", "Summer Salinity"))
hs_maxent_top5 = maxent_vimp %>%
  select(Overall, Variable, Species) %>%
  filter(Species == "Haemulon sciurus" & 
           Variable %in% c("Habitat", "Mangrove Distance", "Slope",
                           "Winter Salinity", "Depth"))

maxent_top5 = rbind(lg_maxent_top5, hs_maxent_top5)

# lasso plot
lasso_plot = ggplot(lasso_top5, aes(x = interaction(Variable, Species), 
                                    y = Overall, fill = Species, group = Species)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  scale_x_discrete(breaks = interaction(lasso_top5$Variable, lasso_top5$Species),
                   labels = lasso_top5$Variable) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 105)) +
  xlab(" ") + ylab("Lasso Regression Variable Importance") +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.margin = unit(c(0,0,0,0),"cm"),
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(color = "black"),
                     axis.text.y = element_text(color = "black"),
                     legend.text = element_text(face = "italic"),
                     legend.box.margin = margin(t = 0, r = 0, b = -10, l = 0)) +
  geom_text(aes(label = round(Overall, digits = 2)), size = 3, 
            position = position_dodge(width = 1), hjust = -0.15, vjust = 0.45) + 
  coord_flip() 
lasso_plot

ggsave(plot = lasso_plot, filename = paste0(temp_plots, "Lasso_Top5_Variables.png"), 
       height = 3.5, width = 7, units = "in", dpi = 450)

# ridge plot
ridge_plot = ggplot(ridge_top5, aes(x = interaction(Variable, Species), 
                                    y = Overall, fill = Species, group = Species)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  scale_x_discrete(breaks = interaction(ridge_top5$Variable, ridge_top5$Species),
                   labels = ridge_top5$Variable) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 105)) +
  xlab(" ") + ylab("Ridge Regression Variable Importance") +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     plot.margin = unit(c(0,0,0,0),"cm"),
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(color = "black"),
                     axis.text.y = element_text(color = "black"),
                     legend.text = element_text(face = "italic"),
                     legend.box.margin = margin(t = 0, r = 0, b = -10, l = 0)) +
  geom_text(aes(label = round(Overall, digits = 2)), size = 3, 
            position = position_dodge(width = 1), hjust = -0.15, vjust = 0.45) + 
  coord_flip() 
ridge_plot

ggsave(plot = ridge_plot, filename = paste0(temp_plots, "Ridge_Top5_Variables.png"), 
       height = 3.5, width = 7, units = "in", dpi = 450)

# maxent plot
maxent_plot = ggplot(maxent_top5, aes(x = interaction(Variable, Species), 
                                    y = Overall, fill = Species, group = Species)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  scale_x_discrete(breaks = interaction(maxent_top5$Variable, maxent_top5$Species),
                   labels = maxent_top5$Variable) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 65)) +
  xlab(" ") + ylab("MaxEnt Variable Importance") +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank(),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     plot.margin = unit(c(0,0,0,0),"cm"),
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(color = "black"),
                     axis.text.y = element_text(color = "black"),
                     legend.text = element_text(face = "italic"),
                     legend.box.margin = margin(t = 0, r = 0, b = -10, l = 0)) +
  geom_text(aes(label = round(Overall, digits = 2)), size = 3, 
            position = position_dodge(width = 1), hjust = -0.15, vjust = 0.45) + 
  coord_flip() 
maxent_plot

ggsave(plot = maxent_plot, filename = paste0(temp_plots, "MaxEnt_Top5_Variables.png"), 
       height = 3.5, width = 7, units = "in", dpi = 450)
