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

# Define color scheme
colors_models <- c("MIMICS" = "#1b9e77", 
                   "Millennial" = "#d95f02", 
                   "Century" = "#7570b3")

# Load and prepare AfSIS data
afsis <- read.csv("forcing_data/afsis_ref_updated9.csv", as.is = T)

afsis_data <- afsis %>% 
  tibble::rowid_to_column("Set") %>% 
  dplyr::filter(Depth == "Topsoil") %>% 
  dplyr::filter(CORG <= 20) %>% 
  #Calculate SOC stocks for observational data (kg/m2)
  dplyr::mutate(C_stock = ((CORG/100)*bd_extracted*20)*10) %>% 
  drop_na(npp_modis, pH, Clay_2um, LIG_N, stemp, sm, Clay_63um, bd_extracted,
          C_stock, Litterfall.gC.m2.yr)

afsis_data$Plot <- as.character(afsis_data$Plot)

afsis_data$Cluster <- as.character(afsis_data$Cluster)

########## RANDOM FOREST WITH OBSERVED AFSIS DATA ##########

# Define model configurations
model_configs <- list(
  MIMICS = list(
    name = "MIMICS",
    data_prep = function(df) {
      df %>%
        unite("id", Country:Cluster) %>%
        mutate(npp_modis.gC.m2.yr = npp_modis*1000) %>%
        dplyr::select(id, npp_modis.gC.m2.yr, Clay_2um, LIG_N, stemp, sm, C_stock)
    },
    features = c("npp_modis.gC.m2.yr", "Clay_2um", "LIG_N", "stemp", "sm")
  ),
  Millennial = list(
    name = "Millennial",
    data_prep = function(df) {
      df %>%
        unite("id", Country:Cluster) %>%
        mutate(npp_modis.gC.m2.d = npp_modis*1000/365) %>%
        dplyr::select(id, npp_modis.gC.m2.d, Clay_63um, pH, stemp, sm, C_stock)
    },
    features = c("npp_modis.gC.m2.d", "Clay_63um", "pH", "stemp", "sm")
  ),
  Century = list(
    name = "Century",
    data_prep = function(df) {
      df %>%
        unite("id", Country:Cluster) %>%
        mutate(npp_modis.gC.m2.d = npp_modis*1000/365,
               Lig_frc = LIG/100) %>%
        dplyr::select(id, npp_modis.gC.m2.d, Clay_63um, LIG_N, Lig_frc, stemp, sm, C_stock)
    },
    features = c("npp_modis.gC.m2.d", "Clay_63um", "stemp", "sm", "LIG_N", "Lig_frc")
  )
)

# Function to train and evaluate RF model
train_rf_model <- function(data, model_name) {
  task <- as_task_regr(x = data, target = "C_stock")
  task$set_col_roles("id", roles = "group")
  
  learner <- lrn("regr.ranger", importance = "permutation", num.trees = 1000)
  
  set.seed(42)
  resampling <- rsmp("cv", folds = 10)
  resampling$instantiate(task)
  
  rf_result <- mlr3::resample(task = task, learner = learner,
                              resampling = resampling, store_models = TRUE)
  
  # Get performance metrics
  metrics_raw <- rf_result$aggregate(measures = msrs(c("regr.rsq", "regr.mae", 
                                                       "regr.rmse")))
  
  # Convert to data.frame with model name
  metrics_df <- data.frame(
    model = model_name,
    rsq = as.numeric(metrics_raw["rsq"]),
    mae = as.numeric(metrics_raw["regr.mae"]),
    rmse = as.numeric(metrics_raw["regr.rmse"])
  )
  
  # Get predictions
  rf_pred <- rf_result$prediction(predict_sets = "test")
  pred_df <- data.frame(
    model = model_name,
    truth = rf_pred$truth,
    response = rf_pred$response
  )
  
  # Get variable importance
  vi <- lapply(rf_result$learners, function(x) x$model$variable.importance)
  vi_df <- vi %>%
    plyr::ldply() %>%
    pivot_longer(everything(), names_to = "variable", values_to = "x") %>%
    dplyr::summarise(median = median(x, na.rm = TRUE),
                     mad = mad(x, na.rm = TRUE),
                     .by = variable) %>%
    dplyr::arrange(median) %>%
    dplyr::mutate(median_pct = (median / sum(median, na.rm = TRUE)) * 100,
                  mad_pct = (mad / sum(median, na.rm = TRUE)) * 100,
                  model = model_name)
  
  list(
    metrics = metrics_df,
    predictions = pred_df,
    vi = vi_df,
    model = rf_result,
    task_pdp = task
  )
}

# Train all models
all_models <- purrr::map(model_configs, ~train_rf_model(.$data_prep(afsis_data), 
                                                        .$name))

# Extract metrics
all_metrics <- map_df(all_models, ~.$metrics)

# Extract and combine predictions
all_predictions <- map_df(all_models, ~.$predictions)

# Extract and combine variable importance
all_vi <- map_df(all_models, ~.$vi)

# Create metrics labels for the plot
metrics_for_plot <- all_metrics %>%
  mutate(
    label = paste0("R² = ", round(rsq, 2), "\n",
                   "MAE = ", round(mae, 2), "\n",
                   "RMSE = ", round(rmse, 2))
  )

# Create labels for facets
facet_labels <- c(
  MIMICS = "c) MIMICS",
  Millennial = "b) Millennial",
  Century = "a) Century"
)

# COMBINED OBSERVED VS PREDICTED PLOT
combined_pred_plot <- all_predictions %>%
  ggplot(aes(x = response, y = truth, color = model)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~model, labeller = labeller(model = facet_labels)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = colors_models, guide = "none") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        strip.text = element_text(size = 12, face = "bold", hjust = 0),
        strip.background = element_blank()) +
  scale_y_continuous("Observed SOC stocks [kg/m²]",
                     limits = c(0, 25), expand = c(0, 0)) +
  scale_x_continuous("Predicted SOC stocks [kg/m²]",
                     limits = c(0, 25), expand = c(0, 0)) +
  # Add performance metrics as text annotations
  geom_text(data = metrics_for_plot, 
            aes(x = 19, y = 6, label = label),
            inherit.aes = FALSE,
            vjust = "top", hjust = "left",
            size = 3)

print(combined_pred_plot)

# COMBINED VARIABLE IMPORTANCE PLOT
# Normalize variable names across models for comparison
all_vi_normalized <- all_vi %>%
  dplyr::mutate(variable = case_when(
    variable == "npp_modis.gC.m2.yr" ~ "NPP",
    variable == "npp_modis.gC.m2.d" ~ "NPP",
    variable == "Clay_2um" ~ "Clay",
    variable == "Clay_63um" ~ "Clay",
    variable == "stemp" ~ "Soil temperature",
    variable == "sm" ~ "Soil moisture",
    variable == "LIG_N" ~ "Lignin:N ratio",
    variable == "Lig_frc" ~ "Lignin",
    TRUE ~ variable
  )) %>%
  group_by(model, variable) %>%
  dplyr::summarise(across(c(median_pct, mad_pct), mean),
                   .groups = "drop")

combined_vi_plot <- all_vi_normalized %>%
  ggplot(aes(x = reorder(variable, -median_pct), y = median_pct, 
             fill = model, color = model)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(ymin = median_pct - mad_pct,
                    ymax = median_pct + mad_pct),
                position = position_dodge(width = 0.9),
                width = 0.15, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = colors_models, name = "Model") +
  scale_color_manual(values = colors_models, guide = "none") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top") +
  scale_x_discrete("") +
  scale_y_continuous("Relative explained variation (%)", 
                     expand = c(0, 0), limits = c(0, 35))

print(combined_vi_plot)

# PARTIAL DEPENDENCE PLOTS 
# Create a function for PDP
create_pdp <- function(model_list, model_name, config, data) {
  task_pdp <- as_task_regr(x = config$data_prep(data) %>%
                             dplyr::select(-id),
                           target = "C_stock")
  
  learner_pdp <- lrn("regr.ranger", importance = "permutation", num.trees = 1000)
  set.seed(42)
  learner_pdp$train(task_pdp)
  
  model_rf <- Predictor$new(learner_pdp,
                            data = config$data_prep(data) %>%
                              dplyr::select(-id))
  
  effect_pdp <- FeatureEffects$new(model_rf, method = "pdp",
                                   features = config$features)
  return(effect_pdp)
}

pdp_results <- purrr::map2(
  seq_along(model_configs),
  names(model_configs),
  ~create_pdp(all_models[[.y]], .y, model_configs[[.y]], afsis_data)
)

# Function to extract and process PDP results
extract_pdp_data <- function(pdp_results, model_names) {
  map_df(
    seq_along(pdp_results),
    ~rbindlist(pdp_results[[.x]]$results) %>%
      tibble() %>%
      dplyr::select(-.type) %>%
      dplyr::rename(
        feature_value = .borders,
        predicted_value = .value,
        feature_name = .feature
      ) %>%
      dplyr::mutate(model = model_names[.x])
  )
}

# Extract PDP data
pdp_obs_all <- extract_pdp_data(pdp_results, names(model_configs)) %>% 
  dplyr::mutate(data_source = "Observed")

head(pdp_obs_all)

#Function to plot PDP for each model
plot_pdp_fun <- function(pdp_data, model_name) {
  pdp_data %>%
    dplyr::filter(model == model_name) %>%
    ggplot(aes(y = predicted_value, x = feature_value)) +
    geom_path() +
    facet_wrap(~feature_name, scales = "free_x") +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(color = "black"),
          strip.text = element_text(size = 12, face = "bold")) +
    scale_y_continuous("Predicted SOC stock [kg m-2]", expand = c(0,0),
                       limits = c(0,10)) +
    scale_x_continuous("Predictor range", expand = c(0,0))
}

# Plot PDP's
plot_pdp_fun(pdp_obs_all, "MIMICS")

plot_pdp_fun(pdp_obs_all, "Millennial")

plot_pdp_fun(pdp_obs_all, "Century")

####################### RANDOM FOREST WITH FITTED MODEL RUNS #####################

# Load fitted model runs
fitted_run <- read.csv("model_output/all_model_results_fitted_run_2025-12-03.csv")

fitted_run %>% 
  group_by(Type) %>% 
  reframe(mean = mean(Soil_Organic_Carbon_kg_m2))

afsis_red <- afsis_data %>%
  dplyr::select(SSN, Longitude, Latitude, Region, Country, Site, Cluster, Plot,
                Clay_2um, Clay_63um, pH, bd_extracted, npp_modis, sm, stemp, 
                LIG, LIG_N)

# Merge datasets
models_afsis <- fitted_run %>%
  dplyr::rename(C_stock = Soil_Organic_Carbon_kg_m2) %>%
  left_join(afsis_red, by = "SSN")

models_afsis %>% 
  group_by(Type) %>% 
  reframe(mean = mean(C_stock))

# Define model configurations
model_run_configs <- list(
  MIMICS = list(
    name = "MIMICS",
    type = "MIMICS",
    data_prep = function(df) {
      df %>%
        filter(Type == "MIMICS") %>%
        unite("id", Country:Cluster) %>%
        mutate(npp_modis.gC.m2.yr = npp_modis*1000) %>%
        dplyr::select(id, npp_modis.gC.m2.yr, Clay_2um, LIG_N, stemp, sm,
                      C_stock)
    },
    features = c("npp_modis.gC.m2.yr", "Clay_2um", "LIG_N", "stemp", "sm")
  ),
  Millennial = list(
    name = "Millennial",
    type = "Millennial",
    data_prep = function(df) {
      df %>%
        filter(Type == "Millennial") %>%
        unite("id", Country:Cluster) %>%
        mutate(npp_modis.gC.m2.d = npp_modis*1000/365) %>%
        dplyr::select(id, npp_modis.gC.m2.d, Clay_63um, pH, stemp, sm,
                      C_stock)
    },
    features = c("npp_modis.gC.m2.d", "Clay_63um", "pH", "stemp", "sm")
  ),
  Century = list(
    name = "Century",
    type = "Century",
    data_prep = function(df) {
      df %>%
        filter(Type == "Century") %>%
        unite("id", Country:Cluster) %>%
        mutate(npp_modis.gC.m2.d = npp_modis*1000/365,
               Lig_frc = LIG/100) %>%
        dplyr::select(id, npp_modis.gC.m2.d, Clay_63um, LIG_N, Lig_frc, stemp, sm,
                      C_stock)
    },
    features = c("npp_modis.gC.m2.d", "Clay_63um", "stemp", "sm", "LIG_N", "Lig_frc")
  )
)

# Train all models
all_models_fitted <- purrr::map(model_run_configs, 
                         ~train_rf_model(.$data_prep(models_afsis), .$name))

# Extract results
all_metrics_fitted <- map_df(all_models_fitted, ~.$metrics)
all_predictions_fitted <- map_df(all_models_fitted, ~.$predictions)
all_vi_fitted <- map_df(all_models_fitted, ~.$vi)

# Create metrics labels
metrics_for_plot_fitted <- all_metrics_fitted %>%
  mutate(
    label = paste0("R² = ", round(rsq, 2), "\n",
                   "MAE = ", round(mae, 2), "\n",
                   "RMSE = ", round(rmse, 2))
  )

# Facet labels
facet_labels_fitted <- c(
  MIMICS = "c) MIMICS",
  Millennial = "b) Millennial",
  Century = "a) Century"
)

# COMBINED OBSERVED VS PREDICTED PLOT (Model runs)
combined_pred_plot_fitted <- all_predictions_fitted %>%
  ggplot(aes(x = response, y = truth, color = model)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~model, labeller = labeller(model = facet_labels_fitted)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = colors_models, guide = "none") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        strip.text = element_text(size = 12, face = "bold", hjust = 0),
        strip.background = element_blank()) +
  scale_y_continuous("Modeled SOC stocks [kg/m²]",
                     limits = c(-1, 43), expand = c(0, 0)) +
  scale_x_continuous("Predicted SOC stocks [kg/m²]",
                     limits = c(-1, 25), expand = c(0, 0)) +
  geom_text(data = metrics_for_plot_fitted,
            aes(x = 15, y = 8, label = label),
            inherit.aes = FALSE,
            vjust = "top", hjust = "left",
            size = 3)

print(combined_pred_plot_fitted)

# COMBINED VARIABLE IMPORTANCE PLOT (Model runs)
all_vi_normalized_fitted <- all_vi_fitted %>%
  mutate(variable = case_when(
    variable == "npp_modis.gC.m2.yr" ~ "NPP",
    variable == "npp_modis.gC.m2.d" ~ "NPP",
    variable == "Clay_2um" ~ "Clay",
    variable == "Clay_63um" ~ "Clay",
    variable == "stemp" ~ "Soil temperature",
    variable == "sm" ~ "Soil moisture",
    variable == "LIG_N" ~ "Lignin:N ratio",
    variable == "Lig_frc" ~ "Lignin",
    TRUE ~ variable
  )) %>%
  group_by(model, variable) %>%
  dplyr::summarise(across(c(median_pct, mad_pct), mean),
                   .groups = "drop")

combined_vi_plot_fitted <- all_vi_normalized_fitted %>%
  ggplot(aes(x = reorder(variable, -median_pct), y = median_pct,
             fill = model, color = model)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(ymin = median_pct - mad_pct,
                    ymax = median_pct + mad_pct),
                position = position_dodge(width = 0.9),
                width = 0.15, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = colors_models, name = "Model") +
  scale_color_manual(values = colors_models, guide = "none") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top") +
  scale_x_discrete("") +
  scale_y_continuous("Relative explained variation (%)",
                     expand = c(0, 0), limits = c(0, 55))

print(combined_vi_plot_fitted)

# PARTIAL DEPENDENCE PLOTS
# Run PDP function
pdp_results_fitted <- purrr::map2(
  seq_along(model_run_configs),
  names(model_run_configs),
  ~create_pdp(all_models_fitted[[.y]], .y, model_run_configs[[.y]], models_afsis)
)

# Extract PDP data
pdp_fit_all <- extract_pdp_data(pdp_results_fitted, names(model_run_configs)) %>% 
  mutate(data_source = "Fitted")

head(pdp_fit_all)

# Plot PDP's
plot_pdp_fun(pdp_fit_all, "MIMICS") +
  scale_y_continuous("Predicted SOC stock [kg m-2]", expand = c(0,0),
                     limits = c(0,20))

plot_pdp_fun(pdp_fit_all, "Millennial") +
  scale_y_continuous("Predicted SOC stock [kg m-2]", expand = c(0,0),
                     limits = c(0,20))

plot_pdp_fun(pdp_fit_all, "Century") +
  scale_y_continuous("Predicted SOC stock [kg m-2]", expand = c(0,0),
                     limits = c(0,20))

####################### RANDOM FOREST WITH SENSITIVITY MODEL RUNS #####################

# Define sensitivity model configurations WITH data_prep functions
sens_model_configs <- list(
  MIMICS = list(
    name = "MIMICS",
    data_prep = function(df) {
      df %>%
        dplyr::select(Soil_Organic_Carbon_kg_m2, ANPP, CLAY, LIG_N, TSOI, THETA_LIQ) %>%
        setNames(c("C_stock", "npp_modis.gC.m2.yr", "Clay_2um", "LIG_N", "stemp", "sm")) %>% 
        #REDUCE NUMBER FOR FASTER PROCESSING
        slice_sample(n = 5000)
    },
    features = c("npp_modis.gC.m2.yr", "Clay_2um", "LIG_N", "stemp", "sm")
  ),
  Millennial = list(
    name = "Millennial",
    data_prep = function(df) {
      df %>%
        filter(MAOM > 0) %>%
        dplyr::select(Soil_Organic_Carbon_kg_m2, forc_npp, param_claysilt, param_pH, forc_st, forc_sw) %>%
        setNames(c("C_stock", "npp_modis.gC.m2.d", "Clay_63um", "pH", "stemp", "sm")) %>% 
        #REDUCE NUMBER FOR FASTER PROCESSING
        slice_sample(n = 5000) 
    },
    features = c("npp_modis.gC.m2.d", "Clay_63um", "pH", "stemp", "sm")
  ),
  Century = list(
    name = "Century",
    data_prep = function(df) {
      df %>%
        dplyr::select(Soil_Organic_Carbon_kg_m2, forc_npp, param_claysilt, LigFrac, LN, forc_st, forc_sw) %>%
        setNames(c("C_stock", "npp_modis.gC.m2.d", "Clay_63um", "Lig_frc", "LIG_N", "stemp", "sm")) %>% 
        #REDUCE NUMBER FOR FASTER PROCESSING
        slice_sample(n = 5000)
    },
    features = c("npp_modis.gC.m2.d", "Clay_63um", "stemp", "sm", "LIG_N", "Lig_frc")
  )
)

# Read the sensitivity files
sens_mimics_raw <- read.csv("./model_output/MIMICS_SensitivityAnalysisOutput_2025-11-18.csv")
mean(sens_mimics_raw$Soil_Organic_Carbon_kg_m2)
sens_millennial_raw <- read.csv("./model_output/Millennial_SensitivityAnalysisOutput_2025-11-18.csv") %>% 
  dplyr::filter(MAOM > 0)
mean(sens_millennial_raw$Soil_Organic_Carbon_kg_m2)
sens_century_raw <- read.csv("./model_output/Century_SensitivityAnalysisOutput_2025-12-11.csv") %>% 
  mutate(param_claysilt = param_claysilt*100)
mean(sens_century_raw$Soil_Organic_Carbon_kg_m2)

# Function to train RF model on sensitivity data
train_sens_rf_model <- function(data, model_name, features) {
  task <- as_task_regr(x = data, target = "C_stock")
  learner <- lrn("regr.ranger", importance = "permutation", num.trees = 1000)
  set.seed(42)
  resampling <- rsmp("cv", folds = 10)
  resampling$instantiate(task)
  rf_result <- mlr3::resample(task = task, learner = learner,
                              resampling = resampling, store_models = TRUE)
  # Get performance metrics
  metrics_raw <- rf_result$aggregate(measures = msrs(c("regr.rsq", "regr.mae",
                                                       "regr.rmse")))
  metrics_df <- data.frame(
    model = model_name,
    rsq = as.numeric(metrics_raw["rsq"]),
    mae = as.numeric(metrics_raw["regr.mae"]),
    rmse = as.numeric(metrics_raw["regr.rmse"])
  )
  # Get predictions
  rf_pred <- rf_result$prediction(predict_sets = "test")
  pred_df <- data.frame(
    model = model_name,
    truth = rf_pred$truth,
    response = rf_pred$response
  )
  # Get variable importance
  vi <- lapply(rf_result$learners, function(x) x$model$variable.importance)
  vi_df <- vi %>%
    plyr::ldply() %>%
    pivot_longer(everything(), names_to = "variable", values_to = "x") %>%
    dplyr::summarise(median = median(x, na.rm = TRUE),
                     mad = mad(x, na.rm = TRUE),
                     .by = variable) %>%
    dplyr::arrange(median) %>%
    dplyr::mutate(median_pct = (abs(median) / sum(abs(median), na.rm = TRUE)) * 100,
                  mad_pct = (mad / sum(abs(median), na.rm = TRUE)) * 100,
                  model = model_name)
  list(
    metrics = metrics_df,
    predictions = pred_df,
    vi = vi_df,
    model = rf_result,
    task = task
  )
}

# Load and train sensitivity models
all_models_sens <- purrr::map(sens_model_configs,
                              ~train_sens_rf_model(.$data_prep(get(paste0("sens_", 
                                                                          tolower(.$name), "_raw"))), .$name))

# Extract results
all_metrics_sens <- map_df(all_models_sens, ~.$metrics)
all_predictions_sens <- map_df(all_models_sens, ~.$predictions)
all_vi_sens <- map_df(all_models_sens, ~.$vi)

# Create metrics labels
metrics_for_plot_sens <- all_metrics_sens %>%
  mutate(
    label = paste0("R² = ", round(rsq, 2), "\n",
                   "MAE = ", round(mae, 2), "\n",
                   "RMSE = ", round(rmse, 2))
  )

# Facet labels
facet_labels_sens <- c(
  MIMICS = "c) MIMICS",
  Millennial = "b) Millennial",
  Century = "a) Century"
)

# COMBINED OBSERVED VS PREDICTED PLOT
combined_pred_plot_sens <- all_predictions_sens %>%
  ggplot(aes(x = response, y = truth, color = model)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~model, labeller = labeller(model = facet_labels_sens)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = colors_models, guide = "none") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        strip.text = element_text(size = 12, face = "bold", hjust = 0),
        strip.background = element_blank()) +
  scale_y_continuous("Modeled SOC stocks [kg/m²]",
                     limits = c(-2, 41), expand = c(0, 0)) +
  scale_x_continuous("Predicted SOC stocks [kg/m²]",
                     limits = c(-2, 41), expand = c(0, 0)) +
  geom_text(data = metrics_for_plot_sens,
            aes(x = 20, y = 8, label = label),
            inherit.aes = FALSE,
            vjust = "top", hjust = "left",
            size = 3)

print(combined_pred_plot_sens)

# COMBINED VARIABLE IMPORTANCE PLOT 
all_vi_normalized_sens <- all_vi_sens %>%
  mutate(variable = case_when(
    variable == "npp_modis.gC.m2.yr" ~ "NPP",
    variable == "npp_modis.gC.m2.d" ~ "NPP",
    variable == "Clay_2um" ~ "Clay",
    variable == "Clay_63um" ~ "Clay",
    variable == "stemp" ~ "Soil temperature",
    variable == "sm" ~ "Soil moisture",
    variable == "LIG_N" ~ "Lignin:N ratio",
    variable == "Lig_frc" ~ "Lignin",
    TRUE ~ variable
  )) %>% 
  group_by(model, variable) %>%
  dplyr::summarise(across(c(median_pct, mad_pct), mean),
                   .groups = "drop")

combined_vi_plot_sens <- all_vi_normalized_sens %>%
  ggplot(aes(x = reorder(variable, -median_pct), y = median_pct,
             fill = model, color = model)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(ymin = median_pct - mad_pct,
                    ymax = median_pct + mad_pct),
                position = position_dodge(width = 0.9),
                width = 0.15, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = colors_models, name = "Model") +
  scale_color_manual(values = colors_models, guide = "none") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "top") +
  scale_x_discrete("") +
  scale_y_continuous("Relative explained variation (%)",
                     expand = c(0, 0), limits = c(0, 80))

print(combined_vi_plot_sens)

# PARTIAL DEPENDENCE PLOTS 
# Create a function for PDP
options(future.globals.maxSize = 6 * 1024^3) 
create_pdp_sens <- function(model_list, model_name, config, data) {
  task_pdp <- as_task_regr(x = config$data_prep(data),
                           target = "C_stock")
  
  learner_pdp <- lrn("regr.ranger", importance = "permutation", num.trees = 1000)
  set.seed(42)
  learner_pdp$train(task_pdp)
  
  model_rf <- Predictor$new(learner_pdp,
                            data = config$data_prep(data))
  
  effect_pdp <- FeatureEffects$new(model_rf, method = "pdp",
                                   features = config$features)
  return(effect_pdp)
}

# Run PDP function
pdp_results_sens <- purrr::map2(
  seq_along(sens_model_configs),
  names(sens_model_configs),
  ~create_pdp_sens(all_models_sens[[.y]], .y, sens_model_configs[[.y]], 
                   get(paste0("sens_", tolower(.y), "_raw")))
)

# Extract PDP data
pdp_sens_all <- extract_pdp_data(pdp_results_sens, names(sens_model_configs)) %>% 
  mutate(data_source = "Sensitivity")

head(pdp_sens_all)

# Plot PDP's
plot_pdp_fun(pdp_sens_all, "MIMICS") +
  scale_y_continuous("Predicted SOC stock [kg m-2]", expand = c(0,0),
                     limits = c(0,15))

plot_pdp_fun(pdp_sens_all, "Millennial") +
  scale_y_continuous("Predicted SOC stock [kg m-2]", expand = c(0,0),
                     limits = c(0,15))

plot_pdp_fun(pdp_sens_all, "Century") +
  scale_y_continuous("Predicted SOC stock [kg m-2]", expand = c(0,0),
                     limits = c(0,15))

####################### COMPARE PROCESS-BASED MODELS ACROSS RF APPROACHES #######################

# Combine all predictions with their data source
all_predictions_combined <- bind_rows(
  all_predictions %>% 
    mutate(data_source = "Observed"),
  all_predictions_fitted %>% 
    mutate(data_source = "Fitted"),
  all_predictions_sens %>% 
    mutate(data_source = "Sensitivity")
)

# Create mapping for model names
model_order <- c("Century", "Millennial", "MIMICS")
source_order <- c("Observed", "Fitted", "Sensitivity")

# Define facet labels columns (data sources)
facet_labels_sources <- c(
  Observed = "Observed Data",
  Fitted = "Fitted Model Runs",
  Sensitivity = "Sensitivity Analysis"
)

# Combine all metrics
all_metrics_comparison <- bind_rows(
  all_metrics %>% mutate(data_source = "Observed"),
  all_metrics_fitted %>% mutate(data_source = "Fitted"),
  all_metrics_sens %>% mutate(data_source = "Sensitivity")
)

# Prepare metrics for annotation
metrics_for_annotation <- all_metrics_comparison %>%
  mutate(model = factor(model, levels = model_order),
         data_source = factor(data_source, levels = source_order)) %>%
  mutate(label = paste0("R² = ", round(rsq, 2), "\n",
                        "MAE = ", round(mae, 2), "\n",
                        "RMSE = ", round(rmse, 2)))

# OBSERVED VS PREDICTED BY PROCESS-BASED MODEL (across RF approaches)
comparison_pred_plot <- all_predictions_combined %>%
  mutate(model = factor(model, levels = model_order),
         data_source = factor(data_source, levels = source_order)) %>%
  ggplot(aes(x = response, y = truth, color = model)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_grid(data_source ~ model, 
             labeller = labeller(data_source = facet_labels_sources)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = colors_models, 
                     guide = "none") +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_blank()) +
  scale_y_continuous("Observed/Modeled SOC stocks [kg/m²]", expand = c(0,0),
                     limits = c(-1,45)) +
  scale_x_continuous("Predicted SOC stocks [kg/m²]", expand = c(0,0),
                     limits = c(-1,30)) +
  # Add metrics as text annotations
  geom_text(data = metrics_for_annotation,
            aes(x = 1, y = 42, label = label),
            inherit.aes = FALSE,
            vjust = "top", hjust = "left",
            size = 3, color = "black")

print(comparison_pred_plot)
ggsave("./figures/all_RF_models_comparison_obs_pred.jpeg",
       comparison_pred_plot, width = 8, height = 6)

### VARIABLE IMPORTANCE COMPARISON
all_vi_comparison <- bind_rows(
  all_vi %>%
    mutate(data_source = "Observed"),
  all_vi_fitted %>%
    mutate(data_source = "Fitted"),
  all_vi_sens %>%
    mutate(data_source = "Sensitivity")) 

# Create variable label mapping
variable_labels <- c(
  "npp_modis.gC.m2.yr" = "NPP [gC m-2 yr-1]",
  "npp_modis.gC.m2.d" = "NPP [gC m-2 d-1]",
  "stemp" = "Soil temperature [°C]",
  "sm" = "Soil moisture [m3 m-3]",
  "Clay_2um" = "Clay < 2 μm [%]",
  "Clay_63um" = "Clay + silt < 63 μm [%]",
  "LIG_N" = "Lignin:N ratio",
  "Lig_frc" = "Lignin [%]",
  "pH" = "pH"
)
  
# Create the plot
vi_comparison_plot <- all_vi_comparison %>%
  mutate(model = factor(model, levels = model_order),
         data_source = factor(data_source, levels = source_order),
         variable_label = dplyr::recode(variable, !!!variable_labels)) %>%
  ggplot(aes(x = reorder(variable_label, -median_pct), y = median_pct,
             fill = data_source)) +
  geom_col(position = "dodge", alpha = 0.85, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, median_pct - mad_pct),
                    ymax = median_pct + mad_pct),
                position = position_dodge(width = 0.9),
                width = 0.2, color = "black", linewidth = 0.4) +
  facet_wrap(~model, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = c("Observed" = "#d8b365",
                               "Fitted" = "#5ab4ac",
                               "Sensitivity" = "#af8dc3"),
                    name = "Random forest model") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", 
                                   size = 10),
        axis.text.y = element_text(color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        legend.position = "top",
        panel.spacing.x = unit(0.8, "lines")) +
  scale_y_continuous("Relative explained variation (%)",
                     expand = c(0, 0), limits = c(0, 65)) +
  scale_x_discrete("")

print(vi_comparison_plot)
ggsave("./figures/all_RF_models_vi_comparison.jpeg",
       vi_comparison_plot, width = 11, height = 5)

# Create PDP plots comparing the three RF approaches for each model

pdp_all <- rbind(pdp_obs_all, pdp_fit_all, pdp_sens_all)

pdp_all$data_source <- factor(pdp_all$data_source, levels = c("Observed", 
                                                              "Fitted", "Sensitivity"),
                              ordered = TRUE)

#Function to plot PDP for each model
plot_pdp_all_fun <- function(pdp_data, model_name) {
  pdp_data %>%
    dplyr::filter(model == model_name) %>%
    ggplot(aes(y = predicted_value, x = feature_value, color = data_source)) +
    geom_path(linewidth = 0.9) +
    facet_wrap(~feature_name, scales = "free_x") +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          legend.position = "top") +
    scale_y_continuous("Predicted SOC stock [kg m-2]", expand = c(0,0),
                       limits = c(0,20)) +
    scale_x_continuous("Predictor range", expand = c(0,0)) +
    scale_color_manual(values = c("Observed" = "#d8b365",
                                  "Fitted" = "#5ab4ac",
                                  "Sensitivity" = "#af8dc3"),
                       name = "Random forest model")
}

# Plot PDP's
plot_pdp_all_fun(pdp_all, "MIMICS")
ggsave("./figures/RF_all_MIMICS_pdp.jpeg",
       width = 10, height = 6)

plot_pdp_all_fun(pdp_all, "Millennial")
ggsave("./figures/RF_all_Millennial_pdp.jpeg",
       width = 10, height = 6)

plot_pdp_all_fun(pdp_all, "Century")
ggsave("./figures/RF_all_Century_pdp.jpeg",
       width = 10, height = 6)

#Function to plot PDP for each model - for Cornell talk
p_clay <- pdp_all %>%
  filter(data_source != "Sensitivity") %>%
  dplyr::filter(feature_name == "Clay_63um"|
                  feature_name == "Clay_2um") %>%
  mutate(feature_name = "Clay content [%]") %>%
  ggplot(aes(y = predicted_value, x = feature_value, color = data_source)) +
  geom_path(linewidth = 1) +
  facet_grid(~ model, labeller = labeller(model = facet_labels)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        strip.text = element_text(size = 12, face = "bold", hjust = 0),
        legend.position = "top",
        strip.background = element_rect(fill = "white", color = NA),
        panel.spacing.x = unit(0.5, "cm"),
        plot.margin = margin(r = 10, b = 10, l = 10),
        axis.title.y = element_blank()) +
  scale_y_continuous("Predicted SOC stock [kg m-2]", expand = c(0,0),
                     limits = c(0,20)) +
  scale_x_continuous("Clay content [%]", expand = c(0,0), limits = c(0,100)) +
  scale_color_manual(values = c("Observed" = "#d8b365",
                                "Fitted" = "#5ab4ac"),
                     name = "Random forest model")

p_npp <- pdp_all %>%
  filter(data_source != "Sensitivity") %>%
  dplyr::filter(feature_name == "npp_modis.gC.m2.yr"|
                  feature_name == "npp_modis.gC.m2.d") %>%
  mutate(feature_value = ifelse(model == "MIMICS" &
                                  (feature_name == "npp_modis.gC.m2.yr"),
                                feature_value / 365,
                                feature_value),
         feature_name = ifelse(feature_name == "npp_modis.gC.m2.yr",
                               "npp_modis.gC.m2.d",
                               feature_name),
         feature_name = "NPP [g C m-2 d-1]") %>%
  ggplot(aes(y = predicted_value, x = feature_value, color = data_source)) +
  geom_path(linewidth = 1) +
  facet_grid(~ model) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        strip.text = element_blank(),
        legend.position = "top",
        strip.background = element_blank(),
        panel.spacing.x = unit(0.5, "cm"),
        plot.margin = margin(r = 10, t = 0, l = 10, b = 10),
        axis.title.y = element_blank()) +
  scale_y_continuous("Predicted SOC stock [kg m-2]", expand = c(0,0),
                     limits = c(0,20)) +
  scale_x_continuous("NPP [gC m-2 d-1]", expand = c(0,0), limits = c(0,6)) +
  scale_color_manual(values = c("Observed" = "#d8b365",
                                "Fitted" = "#5ab4ac"),
                     name = "Random forest model")

annotate_figure(
  ggarrange(p_clay, p_npp, common.legend = TRUE, nrow = 2),
  left = text_grob("Predicted SOC stock [kg m-2]",
                   rot = 90, vjust = 1, hjust = 0.5))
ggsave("./figures/RF_Century_MIMICS_Millennial_pdp_red.jpeg",
       height = 5, width = 9)





