## AfSIS statistical model runs ##
## Sophie von Fromm ##
## 2025-08-27 ##

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(mlr3)
library(mlr3learners)
library(mlr3viz)
library(mlr3spatial)
library(mlr3spatiotempcv) 
library(iml)
library(sf)
library(scales)

## Random forest

# Load and prepare AfSIS data
afsis <- read.csv("forcing_data/afsis_ref_updated8.csv", as.is = T)

afsis_data <- afsis %>% 
  filter(Depth == "Topsoil") %>% 
  filter(CORG <= 20) %>% 
  #Calculate SOC stocks for observational data (kg/m2)
  mutate(C_stock = ((CORG/100)*bd_extracted*20)*10) %>% 
  drop_na(npp_modis, pH, Clay_2um, LIG_N, stemp, sm, Clay_63um, bd_extracted,
          C_stock)

afsis_data$Plot <- as.character(afsis_data$Plot)

afsis_data$Cluster <- as.character(afsis_data$Cluster)

# Select variables (based on model forcing patameters)
afsis_mimics_rf <- afsis_data %>% 
  unite("id", Country:Cluster) %>%
  mutate(npp_modis.gC.m2.yr = npp_modis*1000) %>% 
  dplyr::select(id, npp_modis.gC.m2.yr, Clay_2um, LIG_N, stemp, sm, C_stock)

afsis_millennial_rf <- afsis_data %>% 
  unite("id", Country:Cluster) %>%
  mutate(npp_modis.gC.m2.d = npp_modis*1000/365) %>% 
  dplyr::select(id, npp_modis.gC.m2.d, Clay_63um, pH, stemp, sm, C_stock)

## MIMICS
mimics_task_rf <- as_task_regr(x = afsis_mimics_rf,
                               target = "C_stock")

mimics_lrn_rf <- lrn("regr.ranger", importance = "permutation",
                     num.trees = 1000)

# Add id as group for CV (same id kept together)
mimics_task_rf$set_col_roles("id", roles = "group")
print(mimics_task_rf)

# cross-validation
set.seed(42)
resampling <- rsmp("cv", folds = 10)
resampling$instantiate(mimics_task_rf)

## Train model & check performance
mimics_rf <- mlr3::resample(task = mimics_task_rf, learner = mimics_lrn_rf, 
                            resampling = resampling, store_models = TRUE)

# R2 = 0.43, mae = 1.98, rmse = 2.89
mimics_rf$aggregate(measures = msrs(c("regr.rsq", "regr.mae", "regr.rmse")))

mimics_rf_pred <- mimics_rf$prediction(predict_sets = "test")
mimics_rf_pred_df <- data.frame(truth = mimics_rf_pred$truth,
                                response = mimics_rf_pred$response)

mimics_rf_pred_df %>% 
  ggplot(aes(x = response, y = truth)) +
  geom_point(size = 1) +
  geom_rug(length = unit(0.25, "cm")) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Observed SOC stocks [kg/m2]", 
                     limits = c(0,25), expand = c(0,0)) +
  scale_x_continuous("Predicted SOC stocks [kg/m2]", 
                     limits = c(0,25), expand = c(0,0))
ggsave(paste0("./model_output/RF_MIMICS_obs_pred_cv_10f_",
              Sys.Date(), ".jpeg"), width = 6, height = 6)

mimics_vi <- lapply(mimics_rf$learners, function(x) x$model$variable.importance)

mimics_vi %>%
  plyr::ldply() %>%
  pivot_longer(everything(), names_to = "variable", values_to = "x") %>%
  summarise(median = median(x, na.rm = TRUE),
            mad    = mad(x, na.rm = TRUE),
            .by = variable) %>%
  arrange(median)

## Partial dependence plots
mimics_task_rf_pdp <- as_task_regr(x = afsis_mimics_rf %>% 
                                     dplyr::select(-id),
                                   target = "C_stock")

mimics_lrn_rf_pdp <- lrn("regr.ranger", importance = "permutation",
                         num.trees = 1000) 

set.seed(42)
mimics_lrn_rf_pdp$train(mimics_task_rf_pdp)

mimics_model_rf <- Predictor$new(mimics_lrn_rf_pdp, data = afsis_mimics_rf %>% 
                                   dplyr::select(-id))

#PDP
mimics_effect_rf_pdp <- FeatureEffects$new(mimics_model_rf, method = "pdp",
                                           features = c("npp_modis.gC.m2.yr", 
                                                        "Clay_2um",
                                                        "LIG_N", "stemp", "sm"))

plot(mimics_effect_rf_pdp)


