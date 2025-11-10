## AfSIS process-based model analysis ##
## Sophie von Fromm ##
## 2025-10-08 ##

library(tidyverse)
library(ggpubr)

## Load default and fitted model runs
default_run <- read.csv("model_output/both_model_results_default_run_2025-11-10.csv") %>% 
  mutate(ModelRun = "default") 
fitted_run <- read.csv("model_output/both_model_results_fitted_run_2025-11-10.csv") %>% 
  mutate(ModelRun = "fitted")

## Load AfSIS data
afsis_data <- read.csv("forcing_data/afsis_ref_updated9.csv", as.is = T) %>% 
  tibble::rowid_to_column("Set")

## merge both model datasets
both_runs <- default_run %>% 
  full_join(fitted_run)

both_runs_obs <-  both_runs %>% 
  filter(Type %in% c("MIMICS", "Millennial", "Century")) %>%
  rename(Soil_Organic_Carbon_kg_m2_mod = Soil_Organic_Carbon_kg_m2) %>%
  left_join(
    both_runs %>% 
      filter(Type == "Observed") %>% 
      select(SSN, ModelRun, Soil_Organic_Carbon_kg_m2) %>% 
      rename(Soil_Organic_Carbon_kg_m2_obs = Soil_Organic_Carbon_kg_m2),
    by = c("SSN", "ModelRun")
  ) %>%
  select(SSN, Type, ModelRun, Soil_Organic_Carbon_kg_m2_mod, Soil_Organic_Carbon_kg_m2_obs)

median(both_runs_obs$Soil_Organic_Carbon_kg_m2_obs)

## Plot data and calculate model performance
# Function that returns a tibble with desired statistics for specific subset of data
extract_lm_stats <- function(df, type_filter, modelrun_filter) {
  data_sub <- df %>% filter(Type == type_filter, ModelRun == modelrun_filter)
  
  lm_fit <- lm(Soil_Organic_Carbon_kg_m2_obs ~ Soil_Organic_Carbon_kg_m2_mod, 
               data = data_sub)
  
  summary_lm <- summary(lm_fit)
  glance_lm <- broom::glance(lm_fit)
  
  # RMSE calculation
  rmse <- sqrt(mean((data_sub$Soil_Organic_Carbon_kg_m2_obs - data_sub$Soil_Organic_Carbon_kg_m2_mod)^2))
  
  # p-value of slope (coefficient on Soil_Organic_Carbon_kg_m2_obs)
  p_value <- coef(summary_lm)[2, "Pr(>|t|)"]
  
  tibble(
    Type = type_filter,
    ModelRun = modelrun_filter,
    adj_r_squared = summary_lm$adj.r.squared,
    p_value = p_value,
    rmse = rmse
  )
}

existing_combos <- both_runs_obs %>%
  distinct(Type, ModelRun)

stats_table <- existing_combos %>%
  rowwise() %>%
  do(extract_lm_stats(both_runs_obs, .$Type, .$ModelRun)) %>%
  ungroup() %>% 
  mutate(
    label = sprintf(
      "adj. R² = %.2f\np = %.3g\nRMSE = %.2f",
      adj_r_squared, p_value, rmse
    )
  )

both_runs_obs %>% 
  ggplot(aes(y = Soil_Organic_Carbon_kg_m2_obs, x = Soil_Organic_Carbon_kg_m2_mod,
             color = ModelRun)) +
  geom_abline(slope = 1) +
  geom_point(shape = 21) +
  facet_wrap(~Type) +
  geom_smooth(method = "lm", aes(fill = ModelRun)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_y_continuous("Observed SOC stocks [kg m-2]", limits = c(0,24),
                     expand = c(0,0)) +
  scale_x_continuous("Modeled SOC stocks [kg m-2]", limits = c(0,41),
                     expand = c(0,0)) +
  scale_color_manual("Model run", values = c("#f03b20", "#feb24c")) +
  scale_fill_manual("Model run", values = c("#f03b20", "#feb24c")) +
  geom_text(
    data = stats_table %>% filter(ModelRun == "default"),
    aes(x = 25, y = 22, label = label), color = "#f03b20",
    inherit.aes = FALSE,
    size = 4,
    hjust = 0
  ) +
  geom_text(
    data = stats_table %>% filter(ModelRun == "fitted"),
    aes(x = 25, y = 18.5, label = label), color = "#feb24c",
    inherit.aes = FALSE,
    size = 4,
    hjust = 0
  )
ggsave(paste0("./model_output/ModelFits_Century_MIMICS_Millennial_DefaultFitted_",
              Sys.Date(), ".jpeg"), height = 6, width = 12)  

## Bias plots
afsis_red <- afsis_data %>% 
  filter(Depth == "Topsoil") %>% 
  filter(CORG <= 20) %>% 
  dplyr::select(SSN, Longitude, Latitude, Region, Country, Site, Cluster, Plot,
                Clay_2um, Clay_63um, pH, Am_Ox_Al, Am_Ox_Fe, Caex, Clay_1_1,
                Clay_2_1, LIG_N, bd_extracted, npp_modis, sm, stemp, MAT, PET,
                MAP, LIG)

model_afsis <- both_runs_obs %>% 
  left_join(afsis_red, by = "SSN") %>% 
  mutate(SOC_bias_kg_m2 = Soil_Organic_Carbon_kg_m2_mod - Soil_Organic_Carbon_kg_m2_obs)

# function for bias plots with forcing parameters
bias_plot_fun_forc_prop <- function(xvar, model){
  model_afsis %>% 
    filter(Type == {{model}}) %>% 
    ggplot(aes(y = SOC_bias_kg_m2, x = {{xvar}}, color = ModelRun)) +
    geom_hline(yintercept = 0) +
    geom_point(shape = 21) +
    scale_y_continuous("Bias in SOC stocks [kg m-2]") +
    geom_smooth(aes(fill = ModelRun), method = "lm") +
    theme_classic(base_size = 14) +
    theme(axis.text = element_text(color = "black")) +
    scale_color_manual("Model run", values = c("#f03b20", "#feb24c")) +
    scale_fill_manual("Model run", values = c("#f03b20", "#feb24c"))
}

# Century 
f1_cen <- bias_plot_fun_forc_prop(model = "Century", xvar = npp_modis*1000/52) +
  scale_x_continuous("NPP [gC m-2 w-1]") +
  scale_color_manual("Model run", values = c("#feb24c")) +
  scale_fill_manual("Model run", values = c("#feb24c"))

f2_cen <- bias_plot_fun_forc_prop(model = "Century", xvar = MAP) +
  scale_x_continuous("MAP [mm]") +
  scale_color_manual("Model run", values = c("#feb24c")) +
  scale_fill_manual("Model run", values = c("#feb24c"))

f3_cen <- bias_plot_fun_forc_prop(model = "Century", xvar = PET) +
  scale_x_continuous("PET [mm]") +
  scale_color_manual("Model run", values = c("#feb24c")) +
  scale_fill_manual("Model run", values = c("#feb24c"))

f4_cen <- bias_plot_fun_forc_prop(model = "Century", xvar = MAT) +
  scale_x_continuous("MAT [C]") +
  scale_color_manual("Model run", values = c("#feb24c")) +
  scale_fill_manual("Model run", values = c("#feb24c"))

f5_cen <- bias_plot_fun_forc_prop(model = "Century", xvar = Clay_2um/100) +
  scale_x_continuous("Clay content < 2 um") +
  scale_color_manual("Model run", values = c("#feb24c")) +
  scale_fill_manual("Model run", values = c("#feb24c"))

f6_cen <- bias_plot_fun_forc_prop(model = "Century", xvar = (Clay_63um - Clay_2um)/100) +
  scale_x_continuous("Silt content > 2 um & < 63 um") +
  scale_color_manual("Model run", values = c("#feb24c")) +
  scale_fill_manual("Model run", values = c("#feb24c")) 

f7_cen <- bias_plot_fun_forc_prop(model = "Century", xvar = LIG_N/100) +
  scale_x_continuous("Lignin:N ratio")  +
  scale_color_manual("Model run", values = c("#feb24c")) +
  scale_fill_manual("Model run", values = c("#feb24c"))

f8_cen <- bias_plot_fun_forc_prop(model = "Century", xvar = LIG/100) +
  scale_x_continuous("Lignin")  +
  scale_color_manual("Model run", values = c("#feb24c")) +
  scale_fill_manual("Model run", values = c("#feb24c"))

ggarrange(f1_cen, f2_cen, f3_cen, f4_cen, f5_cen, f6_cen, f7_cen, f8_cen,
          common.legend = TRUE, nrow = 2, ncol = 4)

ggsave(paste0("./model_output/BiasPlots_Century_Fitted_forc_",
              Sys.Date(), ".jpeg"), height = 6, width = 12) 

# Millennial 
f1_mil <- bias_plot_fun_forc_prop(model = "Millennial", xvar = npp_modis*1000/365) +
  scale_x_continuous("NPP [gC m-2 d-1]")

f2_mil <- bias_plot_fun_forc_prop(model = "Millennial", xvar = pH) +
  scale_x_continuous("pH")

f3_mil <- bias_plot_fun_forc_prop(model = "Millennial", xvar = Clay_63um) +
  scale_x_continuous("Clay content < 63 um [%]")

f4_mil <- bias_plot_fun_forc_prop(model = "Millennial", xvar = sm) +
  scale_x_continuous("Soil moisture [m-3 m-3]")

f5_mil <- bias_plot_fun_forc_prop(model = "Millennial", xvar = stemp) +
  scale_x_continuous("Soil temperature [C]")

ggarrange(f1_mil, f2_mil, f3_mil, f4_mil, f5_mil, common.legend = TRUE)

ggsave(paste0("./model_output/BiasPlots_Millennial_DefaultFitted_forc_",
              Sys.Date(), ".jpeg"), height = 6, width = 12) 

# MIMICS
f1_mic <- bias_plot_fun_forc_prop(model = "MIMICS", xvar = npp_modis*1000) +
  scale_x_continuous("NPP [gC m-2 yr-1]", expand = c(0,0))

f2_mic <- bias_plot_fun_forc_prop(model = "MIMICS", xvar = LIG_N) +
  scale_x_continuous("Lignin:N ratio")

f3_mic <- bias_plot_fun_forc_prop(model = "MIMICS", xvar = Clay_2um) +
  scale_x_continuous("Clay content < 2 um [%]", expand = c(0,0))

f4_mic <- bias_plot_fun_forc_prop(model = "MIMICS", xvar = sm) +
  scale_x_continuous("Soil moisture [m-3 m-3]", expand = c(0,0))

f5_mic <- bias_plot_fun_forc_prop(model = "MIMICS", xvar = stemp) +
  scale_x_continuous("Soil temperature [C]", expand = c(0,0))

ggarrange(f1_mic, f2_mic, f3_mic, f4_mic, f5_mic, common.legend = TRUE)

ggsave(paste0("./model_output/BiasPlots_MIMICS_DefaultFitted_forc_",
              Sys.Date(), ".jpeg"), height = 6, width = 12) 

## function for bias plots with soil parameters
# Function to extract linear model statistics
extract_bias_lm_stats <- function(df, xvar_name, type_filter) {
  data_sub <- df %>% 
    filter(Type == type_filter, ModelRun == "fitted")
  
  # Create formula dynamically
  formula <- as.formula(paste("SOC_bias_kg_m2 ~", xvar_name))
  lm_fit <- lm(formula, data = data_sub)
  
  summary_lm <- summary(lm_fit)
  
  # RMSE calculation
  rmse <- sqrt(mean(lm_fit$residuals^2))
  
  # p-value of slope (coefficient on x variable)
  p_value <- coef(summary_lm)[2, "Pr(>|t|)"]
  
  tibble(
    Type = type_filter,
    adj_r_squared = summary_lm$adj.r.squared,
    p_value = p_value,
    rmse = rmse
  )
}

# Plotting function
bias_plot_fun_forc_prop <- function(xvar, xvar_name){  
  
  # Get types
  types <- model_afsis %>% 
    filter(ModelRun == "fitted") %>% 
    pull(Type) %>% 
    unique()
  
  # Calculate stats for both types
  stats_table <- tibble(Type = types) %>%
    rowwise() %>%
    do(extract_bias_lm_stats(model_afsis, xvar_name, .$Type)) %>%
    ungroup() %>%
    mutate(
      label = sprintf(
        "adj. R² = %.2f\np = %.3g\nRMSE = %.2f",
        adj_r_squared, p_value, rmse
      )
    )
  
  # Create the plot
  model_afsis %>%  
    filter(ModelRun == "fitted") %>%  
    ggplot(aes(y = SOC_bias_kg_m2, x = {{xvar}}, color = Type)) +  
    geom_hline(yintercept = 0) +  
    geom_point(shape = 21) +  
    geom_smooth(aes(fill = Type), method = "lm") +  
    theme_classic(base_size = 12) +  
    theme(axis.text = element_text(color = "black"),  
          axis.title.y = element_blank()) +  
    scale_color_manual("Model", values = c("#33a02c",  "#a6cee3", "#1f78b4")) +  
    scale_fill_manual("Model", values = c("#33a02c",  "#a6cee3", "#1f78b4")) +
    geom_text(
      data = stats_table %>% filter(Type == "Century"),
      aes(x = Inf, y = Inf, label = label), 
      color = "#33a02c",
      inherit.aes = FALSE,
      size = 2.5,
      hjust = 3.5,
      vjust = 1.5
    ) +
    geom_text(
      data = stats_table %>% filter(Type == "Millennial"),
      aes(x = Inf, y = Inf, label = label), 
      color = "#a6cee3",
      inherit.aes = FALSE,
      size = 2.5,
      hjust = 2.4,
      vjust = 1.5
    ) +
    geom_text(
      data = stats_table %>% filter(Type == "MIMICS"),
      aes(x = Inf, y = Inf, label = label), 
      color = "#1f78b4",
      inherit.aes = FALSE,
      size = 2.5,
      hjust = 1.1,
      vjust = 1.5
    )
}

# Create plots with the updated function
s1 <- bias_plot_fun_forc_prop(Am_Ox_Al/1000, "Am_Ox_Al") +  
  scale_x_continuous("Oxalate-extractable Al [g kg-1]")  
s2 <- bias_plot_fun_forc_prop(Am_Ox_Fe/1000, "Am_Ox_Fe") +  
  scale_x_continuous("Oxalate-extractable Fe [g kg-1]")  
s3 <- bias_plot_fun_forc_prop(pH, "pH") +  
  scale_x_continuous("pH")  
s4 <- bias_plot_fun_forc_prop(Caex, "Caex") +  
  scale_x_continuous("exchangeable Ca [cmol kg-1 soil-1]")  
s5 <- bias_plot_fun_forc_prop(Clay_2um, "Clay_2um") +  
  scale_x_continuous("Clay content < 2 um [%]")  
s6 <- bias_plot_fun_forc_prop(Clay_63um, "Clay_63um") +  
  scale_x_continuous("Clay + fine silt content < 63 um [%]")  
s7 <- bias_plot_fun_forc_prop(Clay_1_1, "Clay_1_1") +  
  scale_x_continuous("1:1 clay minerals [%]")  
s8 <- bias_plot_fun_forc_prop(Clay_2_1, "Clay_2_1") +  
  scale_x_continuous("2:1 clay minerals [%]")  

# Combine plots
annotate_figure(  
  ggarrange(s1, s2, s3, s4, s5, s6, s7, s8, common.legend = TRUE,  
            ncol = 4, nrow = 2),  
  left = text_grob("Bias in SOC stocks [kg m-2]", rot = 90, size = 14))

ggsave(paste0("./model_output/BiasPlots_Century_MIMICS_Millennial_Fitted_soil_",
              Sys.Date(), ".jpeg"), height = 6, width = 12) 

# Spatial plotting
s9 <- bias_plot_fun_forc_prop(Longitude, "Longitude")
s10 <- bias_plot_fun_forc_prop(Latitude, "Latitude")

annotate_figure(  
  ggarrange(s9, s10, common.legend = TRUE),  
  left = text_grob("Bias in SOC stocks [kg m-2]", rot = 90, size = 14))

ggsave(paste0("./model_output/BiasPlots_Century_MIMICS_Millennial_Fitted_long_lat_",
              Sys.Date(), ".jpeg"), height = 6, width = 12) 
