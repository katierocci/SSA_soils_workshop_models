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
afsis <- read.csv("forcing_data/afsis_ref_updated9.csv", as.is = T)

afsis_data <- afsis %>% 
  tibble::rowid_to_column("Set") %>% 
  dplyr::filter(Depth == "Topsoil") %>% 
  dplyr::filter(CORG <= 20) %>% 
  #Calculate SOC stocks for observational data (kg/m2)
  dplyr::mutate(C_stock = ((CORG/100)*bd_extracted*20)*10) %>% 
  drop_na(npp_modis, pH, Clay_2um, LIG_N, stemp, sm, Clay_63um, bd_extracted,
          C_stock)

afsis_data$Plot <- as.character(afsis_data$Plot)

afsis_data$Cluster <- as.character(afsis_data$Cluster)

# Select variables (based on model forcing parameters)
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

mimics_vi_df <- mimics_vi %>%
  plyr::ldply() %>%
  pivot_longer(everything(), names_to = "variable", values_to = "x") %>%
  summarise(median = median(x, na.rm = TRUE),
            mad    = mad(x, na.rm = TRUE),
            .by = variable) %>%
  arrange(median) %>% 
  mutate(median_pct = (median / sum(median, na.rm = TRUE)) * 100,
         mad_pct = (mad / sum(median, na.rm = TRUE)) * 100) 

mimics_vi_df %>% 
  ggplot(aes(x = reorder(variable, -median_pct), y = median_pct)) +
  geom_col() +
  geom_errorbar(aes(ymin = median_pct - mad_pct,
                    ymax = median_pct + mad_pct),
                width = 0.15) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_discrete("") +
  scale_y_continuous("Relative explained variation (%)", expand = c(0,0),
                     limits = c(0,35))
ggsave(paste0("./model_output/RF_MIMICS_vi_cv_10f_",
              Sys.Date(), ".jpeg"), width = 8, height = 6)

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
ggsave(paste0("./model_output/RF_MIMICS_pdp_",
              Sys.Date(), ".jpeg"), width = 10, height = 6)

#PDP - 2D
mimics_effect_rf_pdp_2d <- FeatureEffect$new(mimics_model_rf, method = "pdp",
                                             feature = c("stemp", "sm"))

plot(mimics_effect_rf_pdp_2d)
ggsave(paste0("./model_output/RF_MIMICS_pdp_2d_sm_stemp",
              Sys.Date(), ".jpeg"), width = 10, height = 6)

## Millennial
millennial_task_rf <- as_task_regr(x = afsis_millennial_rf,
                                   target = "C_stock")

millennial_lrn_rf <- lrn("regr.ranger", importance = "permutation",
                         num.trees = 1000)

# Add id as group for CV (same id kept together)
millennial_task_rf$set_col_roles("id", roles = "group")
print(millennial_task_rf)

# cross-validation
set.seed(42)
resampling_millennial <- rsmp("cv", folds = 10)
resampling_millennial$instantiate(millennial_task_rf)

## Train model & check performance
millennial_rf <- mlr3::resample(task = millennial_task_rf, learner = millennial_lrn_rf, 
                                resampling = resampling_millennial, store_models = TRUE)

# R2 = 0.47, mae = 1.87, rmse = 2.78
millennial_rf$aggregate(measures = msrs(c("regr.rsq", "regr.mae", "regr.rmse")))

millennial_rf_pred <- millennial_rf$prediction(predict_sets = "test")
millennial_rf_pred_df <- data.frame(truth = millennial_rf_pred$truth,
                                    response = millennial_rf_pred$response)

millennial_rf_pred_df %>% 
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
ggsave(paste0("./model_output/RF_Millennial_obs_pred_cv_10f_",
              Sys.Date(), ".jpeg"), width = 6, height = 6)

millennial_vi <- lapply(millennial_rf$learners, function(x) x$model$variable.importance)

millennial_vi_df <- millennial_vi %>%
  plyr::ldply() %>%
  pivot_longer(everything(), names_to = "variable", values_to = "x") %>%
  summarise(median = median(x, na.rm = TRUE),
            mad    = mad(x, na.rm = TRUE),
            .by = variable) %>%
  arrange(median) %>% 
  mutate(median_pct = (median / sum(median, na.rm = TRUE)) * 100,
         mad_pct = (mad / sum(median, na.rm = TRUE)) * 100) 

millennial_vi_df %>% 
  ggplot(aes(x = reorder(variable, -median_pct), y = median_pct)) +
  geom_col() +
  geom_errorbar(aes(ymin = median_pct - mad_pct,
                    ymax = median_pct + mad_pct),
                width = 0.15) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_discrete("") +
  scale_y_continuous("Relative explained variation (%)", expand = c(0,0),
                     limits = c(0,35))
ggsave(paste0("./model_output/RF_Millennial_vi_cv_10f_",
              Sys.Date(), ".jpeg"), width = 8, height = 6)

## Partial dependence plots
millennial_task_rf_pdp <- as_task_regr(x = afsis_millennial_rf %>% 
                                         dplyr::select(-id),
                                       target = "C_stock")

millennial_lrn_rf_pdp <- lrn("regr.ranger", importance = "permutation",
                             num.trees = 1000) 

set.seed(42)
millennial_lrn_rf_pdp$train(millennial_task_rf_pdp)

millennial_model_rf <- Predictor$new(millennial_lrn_rf_pdp, data = afsis_millennial_rf %>% 
                                       dplyr::select(-id))

#PDP
millennial_effect_rf_pdp <- FeatureEffects$new(millennial_model_rf, method = "pdp",
                                               features = c("npp_modis.gC.m2.d", 
                                                            "Clay_63um",
                                                            "pH", "stemp", "sm"))

plot(millennial_effect_rf_pdp)
ggsave(paste0("./model_output/RF_Millennial_pdp_",
              Sys.Date(), ".jpeg"), width = 10, height = 6)

#PDP - 2D
millennial_effect_rf_pdp_2d <- FeatureEffect$new(millennial_model_rf, method = "pdp",
                                                 feature = c("stemp", "sm"))

plot(millennial_effect_rf_pdp_2d)
ggsave(paste0("./model_output/RF_Millennial_pdp_2d_sm_stemp",
              Sys.Date(), ".jpeg"), width = 10, height = 6)

####################### Model runs #####################

## Load fitted model runs
fitted_run <- read.csv("model_output/both_model_results_fitted_run_2025-10-10.csv") 

afsis_red <- afsis_data %>% 
  dplyr::select(Set, Longitude, Latitude, Region, Country, Site, Cluster, Plot,
                Clay_2um, Clay_63um, pH, bd_extracted, npp_modis, sm, stemp)

## merge both model and afsis datasets
models_afsis <-  fitted_run %>% 
  dplyr::filter(Type %in% c("MIMICS", "Millennial")) %>%
  rename(Soil_Organic_Carbon_kg_m2_mod = Soil_Organic_Carbon_kg_m2) %>%
  left_join(afsis_red, by = "Set") 

head(models_afsis)

# Select variables (based on model forcing parameters)
model_afsis_mimics_rf <- models_afsis %>% 
  filter(Type == "MIMICS") %>% 
  unite("id", Country:Cluster) %>%
  mutate(npp_modis.gC.m2.yr = npp_modis*1000) %>% 
  dplyr::select(id, npp_modis.gC.m2.yr, Clay_2um, LIG_N, stemp, sm, 
                Soil_Organic_Carbon_kg_m2_mod)

model_afsis_millennial_rf <- models_afsis %>%
  filter(Type == "Millennial") %>% 
  unite("id", Country:Cluster) %>%
  mutate(npp_modis.gC.m2.d = npp_modis*1000/365) %>%
  dplyr::select(id, npp_modis.gC.m2.d, Clay_63um, pH, stemp, sm, 
                Soil_Organic_Carbon_kg_m2_mod)

## MIMICS
model_mimics_task_rf <- as_task_regr(x = model_afsis_mimics_rf,
                                     target = "Soil_Organic_Carbon_kg_m2_mod")

model_mimics_lrn_rf <- lrn("regr.ranger", importance = "permutation",
                           num.trees = 1000)

# Add id as group for CV (same id kept together)
model_mimics_task_rf$set_col_roles("id", roles = "group")
print(model_mimics_task_rf)

# cross-validation
set.seed(42)
resampling <- rsmp("cv", folds = 10)
resampling$instantiate(model_mimics_task_rf)

## Train model & check performance
model_mimics_rf <- mlr3::resample(task = model_mimics_task_rf, 
                                  learner = model_mimics_lrn_rf, 
                                  resampling = resampling, store_models = TRUE)

# R2 = 0.94, mae = 0.35, rmse = 0.79
model_mimics_rf$aggregate(measures = msrs(c("regr.rsq", "regr.mae", "regr.rmse")))

model_mimics_rf_pred <- model_mimics_rf$prediction(predict_sets = "test")
model_mimics_rf_pred_df <- data.frame(truth = model_mimics_rf_pred$truth,
                                      response = model_mimics_rf_pred$response)

model_mimics_rf_pred_df %>% 
  ggplot(aes(x = response, y = truth)) +
  geom_point(size = 1) +
  geom_rug(length = unit(0.25, "cm")) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Observed SOC stocks [kg/m2]", 
                     limits = c(-2,43), expand = c(0,0)) +
  scale_x_continuous("Predicted SOC stocks [kg/m2]", 
                     limits = c(-2,43), expand = c(0,0))
ggsave(paste0("./model_output/RF_model_MIMICS_obs_pred_cv_10f_",
              Sys.Date(), ".jpeg"), width = 6, height = 6)

model_mimics_vi <- lapply(model_mimics_rf$learners, function(x) x$model$variable.importance)

model_mimics_vi_df <- model_mimics_vi %>%
  plyr::ldply() %>%
  pivot_longer(everything(), names_to = "variable", values_to = "x") %>%
  dplyr::summarise(median = median(x, na.rm = TRUE),
                   mad    = mad(x, na.rm = TRUE),
            .by = variable) %>%
  arrange(median) %>% 
  mutate(median_pct = (median / sum(median, na.rm = TRUE)) * 100,
         mad_pct = (mad / sum(median, na.rm = TRUE)) * 100) 

model_mimics_vi_df %>% 
  ggplot(aes(x = reorder(variable, -median_pct), y = median_pct)) +
  geom_col() +
  geom_errorbar(aes(ymin = median_pct - mad_pct,
                    ymax = median_pct + mad_pct),
                width = 0.15) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_discrete("") +
  scale_y_continuous("Relative explained variation (%)", expand = c(0,0),
                     limits = c(0,45))
ggsave(paste0("./model_output/RF_model_MIMICS_vi_cv_10f_",
              Sys.Date(), ".jpeg"), width = 8, height = 6)

## Partial dependence plots
model_mimics_task_rf_pdp <- as_task_regr(x = model_afsis_mimics_rf %>% 
                                           dplyr::select(-id),
                                         target = "Soil_Organic_Carbon_kg_m2_mod")

model_mimics_lrn_rf_pdp <- lrn("regr.ranger", importance = "permutation",
                               num.trees = 1000) 

set.seed(42)
model_mimics_lrn_rf_pdp$train(model_mimics_task_rf_pdp)

model_mimics_model_rf <- Predictor$new(model_mimics_lrn_rf_pdp, 
                                       data = model_afsis_mimics_rf %>% 
                                         dplyr::select(-id))

#PDP
model_mimics_effect_rf_pdp <- FeatureEffects$new(model_mimics_model_rf, 
                                                 method = "pdp",
                                                 features = c("npp_modis.gC.m2.yr", 
                                                              "Clay_2um",
                                                              "LIG_N", "stemp", "sm"))

plot(model_mimics_effect_rf_pdp)
ggsave(paste0("./model_output/RF_model_MIMICS_pdp_",
              Sys.Date(), ".jpeg"), width = 10, height = 6)

## Millennial
model_millennial_task_rf <- as_task_regr(x = model_afsis_millennial_rf,
                                     target = "Soil_Organic_Carbon_kg_m2_mod")

model_millennial_lrn_rf <- lrn("regr.ranger", importance = "permutation",
                           num.trees = 1000)

# Add id as group for CV (same id kept together)
model_millennial_task_rf$set_col_roles("id", roles = "group")
print(model_millennial_task_rf)

# cross-validation
set.seed(42)
resampling <- rsmp("cv", folds = 10)
resampling$instantiate(model_millennial_task_rf)

## Train model & check performance
model_millennial_rf <- mlr3::resample(task = model_millennial_task_rf, 
                                  learner = model_millennial_lrn_rf, 
                                  resampling = resampling, store_models = TRUE)

# R2 = 0.98, mae = 0.23, rmse = 0.35
model_millennial_rf$aggregate(measures = msrs(c("regr.rsq", "regr.mae", "regr.rmse")))

model_millennial_rf_pred <- model_millennial_rf$prediction(predict_sets = "test")
model_millennial_rf_pred_df <- data.frame(truth = model_millennial_rf_pred$truth,
                                      response = model_millennial_rf_pred$response)

model_millennial_rf_pred_df %>% 
  ggplot(aes(x = response, y = truth)) +
  geom_point(size = 1) +
  geom_rug(length = unit(0.25, "cm")) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Observed SOC stocks [kg/m2]", 
                     limits = c(-1,15), expand = c(0,0)) +
  scale_x_continuous("Predicted SOC stocks [kg/m2]", 
                     limits = c(-1,15), expand = c(0,0))
ggsave(paste0("./model_output/RF_model_millennial_obs_pred_cv_10f_",
              Sys.Date(), ".jpeg"), width = 6, height = 6)

model_millennial_vi <- lapply(model_millennial_rf$learners, function(x) x$model$variable.importance)

model_millennial_vi_df <- model_millennial_vi %>%
  plyr::ldply() %>%
  pivot_longer(everything(), names_to = "variable", values_to = "x") %>%
  dplyr::summarise(median = median(x, na.rm = TRUE),
                   mad    = mad(x, na.rm = TRUE),
                   .by = variable) %>%
  arrange(median) %>% 
  mutate(median_pct = (median / sum(median, na.rm = TRUE)) * 100,
         mad_pct = (mad / sum(median, na.rm = TRUE)) * 100) 

model_millennial_vi_df %>% 
  ggplot(aes(x = reorder(variable, -median_pct), y = median_pct)) +
  geom_col() +
  geom_errorbar(aes(ymin = median_pct - mad_pct,
                    ymax = median_pct + mad_pct),
                width = 0.15) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_discrete("") +
  scale_y_continuous("Relative explained variation (%)", expand = c(0,0),
                     limits = c(0,50))
ggsave(paste0("./model_output/RF_model_millennial_vi_cv_10f_",
              Sys.Date(), ".jpeg"), width = 8, height = 6)

## Partial dependence plots
model_millennial_task_rf_pdp <- as_task_regr(x = model_afsis_millennial_rf %>% 
                                           dplyr::select(-id),
                                         target = "Soil_Organic_Carbon_kg_m2_mod")

model_millennial_lrn_rf_pdp <- lrn("regr.ranger", importance = "permutation",
                               num.trees = 1000) 

set.seed(42)
model_millennial_lrn_rf_pdp$train(model_millennial_task_rf_pdp)

model_millennial_model_rf <- Predictor$new(model_millennial_lrn_rf_pdp, 
                                       data = model_afsis_millennial_rf %>% 
                                         dplyr::select(-id))

#PDP
model_millennial_effect_rf_pdp <- FeatureEffects$new(model_millennial_model_rf, 
                                                 method = "pdp",
                                                 features = c("npp_modis.gC.m2.d", 
                                                              "Clay_63um", "pH",
                                                              "stemp", "sm"))

plot(model_millennial_effect_rf_pdp)
ggsave(paste0("./model_output/RF_model_millennial_pdp_",
              Sys.Date(), ".jpeg"), width = 10, height = 6)





####################### Sensitivity Model runs #####################

## MIMICS ##
mimics_sens <- read.csv("./model_output/MIMICS_SensitivityAnalysisOutput_110525.csv", as.is = T)
head(mimics_sens)

# Check data range
p1 <- mimics_sens %>% 
  ggplot(aes(x = ANPP, y = Soil_Organic_Carbon_kg_m2)) +
  geom_point() +
  geom_rug()
p2 <- mimics_sens %>% 
  ggplot(aes(x = CLAY, y = Soil_Organic_Carbon_kg_m2)) +
  geom_point() +
  geom_rug()
p3 <- mimics_sens %>% 
  ggplot(aes(x = LIG_N, y = Soil_Organic_Carbon_kg_m2)) +
  geom_point() +
  geom_rug()
p4 <- mimics_sens %>% 
  ggplot(aes(x = TSOI, y = Soil_Organic_Carbon_kg_m2)) +
  geom_point() +
  geom_rug()
p5 <- mimics_sens %>% 
  ggplot(aes(x = THETA_LIQ, y = Soil_Organic_Carbon_kg_m2)) +
  geom_point() +
  geom_rug()
ggarrange(p1, p2, p3, p4, p5)
ggsave(paste0("./model_output/MIMICS_sens_data_range_",
              Sys.Date(), ".jpeg"), width = 10, height = 6)

mimics_sens_df_rf <- mimics_sens %>% 
  #select input variables: SOC stocks, NPP, clay 2um, lignin N, stemp, sm
  dplyr::select(Soil_Organic_Carbon_kg_m2, ANPP, CLAY, LIG_N, TSOI, THETA_LIQ)

summary(mimics_sens_df_rf)

mimics_sens_task_rf <- as_task_regr(x = mimics_sens_df_rf,
                               target = "Soil_Organic_Carbon_kg_m2")

mimics_sens_lrn_rf <- lrn("regr.ranger", importance = "permutation",
                          num.trees = 1000)

# Add id as group for CV (same id kept together)
# mimics_task_rf$set_col_roles("id", roles = "group")
# print(mimics_task_rf)

# cross-validation
set.seed(42)
resampling <- rsmp("cv", folds = 10)
resampling$instantiate(mimics_sens_task_rf)

## Train model & check performance
mimics_sens_rf <- mlr3::resample(task = mimics_sens_task_rf, learner = mimics_sens_lrn_rf, 
                                 resampling = resampling, store_models = TRUE)

# R2 = 0.99, mae = 0.04, rmse = 0.18
mimics_sens_rf$aggregate(measures = msrs(c("regr.rsq", "regr.mae", "regr.rmse")))

mimics_sens_rf_pred <- mimics_sens_rf$prediction(predict_sets = "test")
mimics_sens_rf_pred_df <- data.frame(truth = mimics_sens_rf_pred$truth,
                                     response = mimics_sens_rf_pred$response)

mimics_sens_rf_pred_df %>% 
  ggplot(aes(x = response, y = truth)) +
  geom_point(size = 1) +
  geom_rug(length = unit(0.25, "cm")) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Observed SOC stocks [kg/m2]", 
                     limits = c(-2,45), expand = c(0,0)) +
  scale_x_continuous("Predicted SOC stocks [kg/m2]", 
                     limits = c(-2,45), expand = c(0,0))
ggsave(paste0("./model_output/RF_MIMICS_sens_obs_pred_cv_10f_",
              Sys.Date(), ".jpeg"), width = 6, height = 6)

mimics_sens_vi <- lapply(mimics_sens_rf$learners, function(x) x$model$variable.importance)

mimics_sens_vi_df <- mimics_sens_vi %>%
  plyr::ldply() %>%
  pivot_longer(everything(), names_to = "variable", values_to = "x") %>%
  dplyr::summarise(median = median(x, na.rm = TRUE),
                   mad    = mad(x, na.rm = TRUE),
                   .by = variable) %>%
  arrange(median) %>% 
  mutate(median_pct = (abs(median) / sum(abs(median), na.rm = TRUE)) * 100,
         mad_pct = (mad / sum(abs(median), na.rm = TRUE)) * 100) 

mimics_sens_vi_df %>% 
  ggplot(aes(x = reorder(variable, -median_pct), y = median_pct)) +
  geom_col() +
  geom_errorbar(aes(ymin = median_pct - mad_pct,
                    ymax = median_pct + mad_pct),
                width = 0.15) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_discrete("") +
  scale_y_continuous("Relative explained variation (%)", expand = c(0,0),
                     limits = c(0,70))
ggsave(paste0("./model_output/RF_MIMICS_sens_vi_cv_10f_",
              Sys.Date(), ".jpeg"), width = 8, height = 6)

## Partial dependence plots
mimics_sens_task_rf_pdp <- as_task_regr(x = mimics_sens_df_rf,
                                        target = "Soil_Organic_Carbon_kg_m2")

mimics_sens_lrn_rf_pdp <- lrn("regr.ranger", importance = "permutation",
                              num.trees = 1000) 

set.seed(42)
mimics_sens_lrn_rf_pdp$train(mimics_sens_task_rf_pdp)
mimics_sens_model_rf <- Predictor$new(mimics_sens_lrn_rf_pdp, data = mimics_sens_df_rf)

# Increase the maximum allowed size for parallel processing
options(future.globals.maxSize = 4000 * 1024^2)  # ~4 GB

# Create a smaller dataset for PDP calculation
set.seed(42)
n_sample <- min(1000, nrow(mimics_sens_df_rf))  # Use max 1000 rows
sample_idx <- sample(nrow(mimics_sens_df_rf), n_sample)
mimics_sens_df_rf_small <- mimics_sens_df_rf[sample_idx, ]

# Create new predictor with smaller dataset
mimics_sens_model_rf_small <- Predictor$new(mimics_sens_lrn_rf_pdp, data = mimics_sens_df_rf_small)

# Calculate PDP
mimics_sens_effect_rf_pdp <- FeatureEffects$new(mimics_sens_model_rf_small, method = "pdp",
                                                features = c("ANPP", "CLAY", "LIG_N",
                                                             "TSOI", "THETA_LIQ"))
plot(mimics_sens_effect_rf_pdp)
ggsave(paste0("./model_output/RF_MIMICS_sens_pdp_",
              Sys.Date(), ".jpeg"), width = 10, height = 6)

## Millennial ##
millennial_sens <- read.csv("./model_output/Millennial_SensitivityAnalysisOutput_110525.csv", as.is = T) %>% 
  #remove model runs that did not reach steady state and have 0 C in MAOM fraction
  filter(MAOM > 0)
head(millennial_sens)
summary(millennial_sens$MAOM)

# # Create a smaller dataset for faster computing
# set.seed(42)
# n_sample <- min(100000, nrow(millennial_sens))  # Use max 100000 rows
# sample_idx <- sample(nrow(millennial_sens), n_sample)
# millennial_sens_small <- millennial_sens[sample_idx, ]

# Check data range
p1 <- millennial_sens %>% 
  ggplot(aes(x = forc_npp, y = Soil_Organic_Carbon_kg_m2)) +
  geom_point() +
  geom_rug()
p2 <- millennial_sens %>% 
  ggplot(aes(x = param_claysilt, y = Soil_Organic_Carbon_kg_m2)) +
  geom_point() +
  geom_rug()
p3 <- millennial_sens %>% 
  ggplot(aes(x = param_pH, y = Soil_Organic_Carbon_kg_m2)) +
  geom_point() +
  geom_rug()
p4 <- millennial_sens %>% 
  ggplot(aes(x = forc_st, y = Soil_Organic_Carbon_kg_m2)) +
  geom_point() +
  geom_rug()
p5 <- millennial_sens %>% 
  ggplot(aes(x = forc_sw, y = Soil_Organic_Carbon_kg_m2)) +
  geom_point() +
  geom_rug()
ggarrange(p1, p2, p3, p4, p5)
ggsave(paste0("./model_output/Millennial_sens_data_range_",
              Sys.Date(), ".jpeg"), width = 10, height = 6)

millennial_sens_df_rf <- millennial_sens %>% 
  #select input variables: SOC stocks, NPP, clay 2um, lignin N, stemp, sm
  dplyr::select(Soil_Organic_Carbon_kg_m2, forc_npp, param_claysilt, param_pH, 
                forc_st, forc_sw)

summary(millennial_sens_df_rf)

millennial_sens_task_rf <- as_task_regr(x = millennial_sens_df_rf,
                                        target = "Soil_Organic_Carbon_kg_m2")

millennial_sens_lrn_rf <- lrn("regr.ranger", importance = "permutation",
                              num.trees = 1000)

# Add id as group for CV (same id kept together)
# millennial_task_rf$set_col_roles("id", roles = "group")
# print(millennial_task_rf)

# cross-validation
set.seed(42)
resampling <- rsmp("cv", folds = 10)
resampling$instantiate(millennial_sens_task_rf)

## Train model & check performance
millennial_sens_rf <- mlr3::resample(task = millennial_sens_task_rf, 
                                     learner = millennial_sens_lrn_rf, 
                                     resampling = resampling, store_models = TRUE)

# R2 = 0.99, mae = 0.11, rmse = 0.21
millennial_sens_rf$aggregate(measures = msrs(c("regr.rsq", "regr.mae", "regr.rmse")))

millennial_sens_rf_pred <- millennial_sens_rf$prediction(predict_sets = "test")
millennial_sens_rf_pred_df <- data.frame(truth = millennial_sens_rf_pred$truth,
                                         response = millennial_sens_rf_pred$response)

millennial_sens_rf_pred_df %>% 
  ggplot(aes(x = response, y = truth)) +
  geom_point(size = 1) +
  geom_rug(length = unit(0.25, "cm")) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Observed SOC stocks [kg/m2]", 
                     limits = c(-2,40), expand = c(0,0)) +
  scale_x_continuous("Predicted SOC stocks [kg/m2]", 
                     limits = c(-2,40), expand = c(0,0))
ggsave(paste0("./model_output/RF_Millennial_sens_obs_pred_cv_10f_",
              Sys.Date(), ".jpeg"), width = 6, height = 6)

millennial_sens_vi <- lapply(millennial_sens_rf$learners, function(x) x$model$variable.importance)

millennial_sens_vi_df <- millennial_sens_vi %>%
  plyr::ldply() %>%
  pivot_longer(everything(), names_to = "variable", values_to = "x") %>%
  group_by(variable) %>% 
  dplyr::summarise(median = median(x, na.rm = TRUE),
                   mad    = mad(x, na.rm = TRUE)) %>%
  arrange(median) %>% 
  mutate(median_pct = (abs(median) / sum(abs(median), na.rm = TRUE)) * 100,
         mad_pct = (mad / sum(abs(median), na.rm = TRUE)) * 100) 

millennial_sens_vi_df %>% 
  ggplot(aes(x = reorder(variable, -median_pct), y = median_pct)) +
  geom_col() +
  geom_errorbar(aes(ymin = median_pct - mad_pct,
                    ymax = median_pct + mad_pct),
                width = 0.15) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_x_discrete("") +
  scale_y_continuous("Relative explained variation (%)", expand = c(0,0),
                     limits = c(0,70))
ggsave(paste0("./model_output/RF_Millennial_sens_vi_cv_10f_",
              Sys.Date(), ".jpeg"), width = 8, height = 6)

## Partial dependence plots
millennial_sens_task_rf_pdp <- as_task_regr(x = millennial_sens_df_rf,
                                            target = "Soil_Organic_Carbon_kg_m2")

millennial_sens_lrn_rf_pdp <- lrn("regr.ranger", importance = "permutation",
                                  num.trees = 1000) 

set.seed(42)
millennial_sens_lrn_rf_pdp$train(millennial_sens_task_rf_pdp)

millennial_sens_model_rf <- Predictor$new(millennial_sens_lrn_rf_pdp, 
                                          data = millennial_sens_df_rf)

# Increase the maximum allowed size for parallel processing
options(future.globals.maxSize = 4000 * 1024^2)  # ~4 GB

# Create a smaller dataset for PDP calculation
set.seed(42)
n_sample <- min(1000, nrow(millennial_sens_df_rf))  # Use max 1000 rows
sample_idx <- sample(nrow(millennial_sens_df_rf), n_sample)
millennial_sens_df_rf_small <- millennial_sens_df_rf[sample_idx, ]

# Create new predictor with smaller dataset
millennial_sens_model_rf_small <- Predictor$new(millennial_sens_lrn_rf_pdp, 
                                                data = millennial_sens_df_rf_small)

# Calculate PDP
millennial_sens_effect_rf_pdp <- FeatureEffects$new(millennial_sens_model_rf_small, 
                                                    method = "pdp",
                                                    features = c("forc_npp", "param_claysilt", 
                                                                 "param_pH",
                                                                 "forc_st", "forc_sw"))
plot(millennial_sens_effect_rf_pdp)
ggsave(paste0("./model_output/RF_Millennial_sens_pdp_",
              Sys.Date(), ".jpeg"), width = 10, height = 6)
