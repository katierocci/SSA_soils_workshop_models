# PCA space for sensitivity analysis
# Sophie von Fromm
# 2025-09-23

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggpubr)

## Load and prepare AfSIS data
afsis <- read.csv("forcing_data/afsis_ref_updated8.csv", as.is = T)

afsis_data <- afsis %>% 
  filter(Depth == "Topsoil") %>% 
  filter(CORG <= 20) %>% 
  #Calculate SOC stocks for observational data (kg/m2)
  mutate(C_stock = ((CORG/100)*bd_extracted*20)*10) %>% 
  drop_na(npp_modis, pH, Clay_2um, LIG_N, stemp, sm, Clay_63um, bd_extracted,
          C_stock) %>% 
  mutate(npp_modis.gC.m2.yr = npp_modis*1000) %>% 
  mutate(DataType = "observed") %>% 
  rename(Soil_Organic_Carbon_kg_m2 = C_stock,
         ANPP = npp_modis.gC.m2.yr,
         CLAY = Clay_2um,
         TSOI = stemp,
         THETA_LIQ = sm) %>% 
  dplyr::select(ANPP, CLAY, LIG_N, TSOI, THETA_LIQ,
                DataType) 

head(afsis_data)

## MIMICS ##
mimics_sens <- read.csv("./MIMICS_SensitivityAnalysisOutput_091725.csv", as.is = T) %>% 
  mutate(DataType = "MIMICS") %>% 
  dplyr::select(ANPP, CLAY, LIG_N, TSOI, THETA_LIQ,
                DataType)
head(mimics_sens)

## Merge both dataset
afsis_mimics <- rbind(afsis_data, mimics_sens) %>% 
  tibble()
head(afsis_mimics)

## PCA analysis
pca_all <- afsis_mimics %>% 
  dplyr::select(-DataType) %>% 
  #generate PCA space based on AfSIS data
  PCA(graph = FALSE, 
      ind.sup = c((nrow(afsis_data) + 1):(nrow(afsis_data) + nrow(mimics_sens))))

factoextra::fviz_pca_biplot(pca_all, col.var = "black", col.ind.sup = NA,
                            col.ind  = "grey", geom.ind = "point") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none")

## Extract PCA data for easier plotting
ind <- data.frame(pca_all$ind$coord)
ind.sup <- data.frame(pca_all$ind.sup$coord)
pca_all_df <- rbind(ind, ind.sup)

var <- facto_summarize(pca_all, element = "var", 
                       result = c("coord", "contrib", "cos2"))

var %>% 
  dplyr::select(name, Dim.1, Dim.2)

#factor for scaling variables to space of individuals (for Dim.1 and Dim.2)
r <- min((max(ind[, "Dim.1"]) - min(ind[, "Dim.1"])/(max(var[, "Dim.1"]) - 
                                                       min(var[, "Dim.1"]))), 
         (max(ind[, "Dim.2"]) - min(ind[, "Dim.2"])/(max(var[, "Dim.2"]) - 
                                                       min(var[, "Dim.2"]))))
#r <- 4.606812

pca_data <- cbind(afsis_mimics, pca_all_df) %>% 
  tibble()

ggplot() +
  # PCA coordinate lies
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # Sensitivity analysis PCA space
  geom_point(data = pca_data %>% 
               filter(DataType == "MIMICS"),
             aes(x = Dim.1, y = Dim.2, color = DataType)) +
  # AfSIS PCA space
  geom_point(data = pca_data %>% 
               filter(DataType == "observed"),
             aes(x = Dim.1, y = Dim.2, color = DataType)) +
  # PCA arrows 
  geom_segment(data = var, aes(x = 0, xend = Dim.1*r*0.7, y = 0, yend = Dim.2*r*0.7),
               linewidth = 0.75, arrow = arrow(length = unit(0.02, "npc"))) +
  # PCA arrow labels (climate data)
  ggrepel::geom_text_repel(data = var, aes(x = Dim.1*r*0.75, y = Dim.2*r*0.75, 
                                           label = name),
                           size = 4, color = "black", fontface = "bold", 
                           nudge_x = 0.4) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank())

