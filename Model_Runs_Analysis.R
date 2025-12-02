## AfSIS process-based model analysis ##
## Sophie von Fromm ##
## 2025-10-08 ##

library(tidyverse)
library(ggpubr)

## Load default and fitted model runs
default_run <- read.csv("model_output/all_model_results_default_run_2025-12-01.csv") %>% 
  mutate(ModelRun = "default") 
fitted_run <- read.csv("model_output/all_model_results_fitted_run_2025-12-01.csv") %>% 
  mutate(ModelRun = "fitted")

## Load AfSIS data
afsis_data <- read.csv("forcing_data/afsis_ref_updated9.csv", as.is = T) 

# Check data range for AfSIS data (used in this study)
facet_labels <- c(
  Clay_2um = "a) Clay content < 2 um [%]",
  Clay_63um = "b) Clay content < 63 um [%]",
  LIG = "c) Lignin content [%]",
  LIG_N = "d) Lignin to N ratio",
  npp_modis.gC.m2.d = "e) NPP [gC m-2 d-1]",
  pH = "f) pH",
  sm = "g) Soil moisture [m-3 m-3]",
  Soil_Organic_Carbon_kg_m2 = "h) SOC stocks [kg m-2]",
  stemp = "i) Soil temperature [C]"
)
  
afsis_data %>% 
  filter(Depth == "Topsoil") %>% 
  filter(CORG <= 20) %>%
  mutate(Soil_Organic_Carbon_kg_m2 = CORG*10*bd_extracted*20/100,
         npp_modis.gC.m2.d = npp_modis*1000/365) %>% 
  drop_na(npp_modis.gC.m2.d, pH, Clay_2um, LIG_N, stemp, sm, Clay_63um, bd_extracted,
          Litterfall.gC.m2.yr) %>% 
  dplyr::select(SSN, Clay_2um, Clay_63um, pH, npp_modis.gC.m2.d, sm, stemp, LIG, 
                LIG_N, Soil_Organic_Carbon_kg_m2) %>% 
  pivot_longer(!SSN) %>% 
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free", labeller = labeller(name = facet_labels)) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        strip.text = element_text(size = 12, face = "bold", hjust = 0),
        strip.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA)) +
  scale_x_continuous("Variable range") +
  scale_y_continuous("Density")
ggsave("./figures/AfSIS_DataRange.jpeg", height = 6, width = 12)  

## merge both model datasets
both_runs <- default_run %>% 
  full_join(fitted_run)

both_runs_obs <-  both_runs %>% 
  filter(Type %in% c("MIMICS", "Millennial", "Century")) %>%
  dplyr::rename(Soil_Organic_Carbon_kg_m2_mod = Soil_Organic_Carbon_kg_m2) %>%
  left_join(
    both_runs %>% 
      filter(Type == "Observed") %>% 
      select(SSN, ModelRun, Soil_Organic_Carbon_kg_m2) %>% 
      dplyr::rename(Soil_Organic_Carbon_kg_m2_obs = Soil_Organic_Carbon_kg_m2),
    by = c("SSN", "ModelRun")
  ) %>%
  dplyr::select(SSN, Type, ModelRun, Soil_Organic_Carbon_kg_m2_mod, Soil_Organic_Carbon_kg_m2_obs)

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

facet_labels_mod <- c(
  Century = "a) Century",
  Millennial = "b) Millennial",
  MIMICS = "c) MIMICS"
)

both_runs_obs %>% 
  ggplot(aes(y = Soil_Organic_Carbon_kg_m2_obs, x = Soil_Organic_Carbon_kg_m2_mod,
             color = ModelRun)) +
  geom_abline(slope = 1) +
  geom_point(shape = 21) +
  facet_wrap(~Type, labeller = labeller(Type = facet_labels_mod)) +
  geom_smooth(method = "lm", aes(fill = ModelRun)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        strip.text = element_text(face = "bold", hjust = 0),
        strip.background = element_rect(fill = "transparent", color = "white"),
        plot.background = element_rect(fill = "transparent", color = "white"),
        panel.background = element_rect(fill = "transparent", color = "white"),
       legend.background = element_rect(fill = "transparent", color = "white"),
       legend.position = "top") +
  scale_y_continuous("Observed SOC stocks [kg m-2]", limits = c(NA,24),
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(0,24)) +
  scale_x_continuous("Modeled SOC stocks [kg m-2]", limits = c(0,41),
                     expand = c(0,0)) +
  scale_color_manual("Model run", values = c("#f03b20", "#feb24c")) +
  scale_fill_manual("Model run", values = c("#f03b20", "#feb24c")) +
  geom_text(
    data = stats_table %>% filter(ModelRun == "default"),
    aes(x = 24.5, y = 22, label = label), color = "#f03b20",
    inherit.aes = FALSE,
    size = 4,
    hjust = 0
  ) +
  geom_text(
    data = stats_table %>% filter(ModelRun == "fitted"),
    aes(x = 24.5, y = 18, label = label), color = "#feb24c",
    inherit.aes = FALSE,
    size = 4,
    hjust = 0
  )
ggsave("./figures/ModelFits_Century_MIMICS_Millennial_DefaultFitted.jpeg", 
       height = 6, width = 12)  

### Bias plots
afsis_red <- afsis_data %>% 
  filter(Depth == "Topsoil") %>% 
  filter(CORG <= 20) %>% 
  dplyr::select(SSN, Longitude, Latitude, Region, Country, Site, Cluster, Plot,
                Clay_2um, Clay_63um, pH, Am_Ox_Al, Am_Ox_Fe, Caex, Clay_1_1,
                Clay_2_1, LIG_N, bd_extracted, npp_modis, sm, stemp, LIG)

model_afsis <- both_runs_obs %>% 
  left_join(afsis_red, by = "SSN") %>% 
  mutate(SOC_bias_kg_m2 = Soil_Organic_Carbon_kg_m2_mod - Soil_Organic_Carbon_kg_m2_obs)

## Forcing parameters
# Function to extract linear model statistics
extract_bias_forc_lm_stats <- function(df, xvar_name, model_filter, modelrun_filter) {
  data_sub <- df %>%
    filter(Type == model_filter, ModelRun == modelrun_filter)
  
  # Create formula dynamically
  formula <- as.formula(paste("SOC_bias_kg_m2 ~", xvar_name))
  lm_fit <- lm(formula, data = data_sub)
  summary_lm <- summary(lm_fit)
  
  # RMSE calculation
  rmse <- sqrt(mean(lm_fit$residuals^2))
  
  # p-value of slope
  p_value <- coef(summary_lm)[2, "Pr(>|t|)"]
  
  tibble(
    ModelRun = modelrun_filter,
    adj_r_squared = summary_lm$adj.r.squared,
    p_value = p_value,
    rmse = rmse
  )
}

# function for bias plots with forcing parameters
bias_plot_fun_forc <- function(xvar, xvar_name, model){
  # Get the model runs for this model
  modelruns <- model_afsis %>%
    filter(Type == model) %>%
    pull(ModelRun) %>%
    unique()
  
  # Calculate stats for both model runs
  stats_table <- tibble(ModelRun = modelruns) %>%
    rowwise() %>%
    do(extract_bias_forc_lm_stats(model_afsis, xvar_name, model, .$ModelRun)) %>%
    ungroup() %>%
    mutate(
      label = sprintf(
        "adj. R² = %.2f\np = %.3e\nRMSE = %.2f",
        adj_r_squared, p_value, rmse
      )
    )
  
  # Create the plot
  model_afsis %>%
    filter(Type == model) %>%
    ggplot(aes(y = SOC_bias_kg_m2, x = {{xvar}}, color = ModelRun)) +
    geom_hline(yintercept = 0) +
    geom_point(shape = 21, alpha = 0.3) +
    geom_smooth(aes(fill = ModelRun), method = "lm") +
    theme_classic(base_size = 14) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_blank(),
          plot.background = element_rect(fill = "transparent", color = "white"),
          legend.background = element_rect(fill = "transparent", color = "white"),
          panel.background = element_rect(fill = "transparent", color = "white"),
          title = element_text(size = 10)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(-20, 38)) +
    scale_color_manual("Model run", values = c("#f03b20", "#feb24c")) +
    scale_fill_manual("Model run", values = c("#f03b20", "#feb24c")) +
    # Add text annotations for stats
    geom_text(
      data = stats_table %>% filter(ModelRun == "default"),
      aes(x = Inf, y = 30, label = label),
      color = "#f03b20",
      inherit.aes = FALSE,
      size = 2.5,
      hjust = 3.0
    ) +
    geom_text(
      data = stats_table %>% filter(ModelRun == "fitted"),
      aes(x = Inf, y = 30, label = label),
      color = "#feb24c",
      inherit.aes = FALSE,
      size = 2.5,
      hjust = 1.8
    )
}

# Century 
f1_cen <- bias_plot_fun_forc(model = "Century", xvar = npp_modis*1000/365,
                             xvar_name = "npp_modis") +
  ggtitle("a) NPP [gC m-2 d-1]") 

f2_cen <- bias_plot_fun_forc(model = "Century", xvar = Clay_63um, 
                             xvar_name = "Clay_63um") +
  ggtitle("b) Clay content < 63 um [%]")

f3_cen <- bias_plot_fun_forc(model = "Century", xvar = sm, 
                             xvar_name = "sm") +
  ggtitle("c) Soil moisture [m-3 m-3]")

f4_cen <- bias_plot_fun_forc(model = "Century", xvar = stemp,
                             xvar_name = "stemp") +
  ggtitle("d) Soil temperature [C]")

f5_cen <- bias_plot_fun_forc(model = "Century", xvar = LIG, 
                             xvar_name = "LIG") +
  ggtitle("e) Lignin [%]")

f6_cen <- bias_plot_fun_forc(model = "Century", xvar = LIG_N, 
                             xvar_name = "LIG_N") +
  ggtitle("f) Lignin:N ratio")

# Combine plots
annotate_figure(  
  ggarrange(f1_cen, f2_cen, f3_cen, f4_cen, f5_cen, f6_cen,
            common.legend = TRUE, nrow = 2, ncol = 3),  
  left = text_grob("Bias in SOC stocks [kg m-2]", rot = 90, size = 13))

ggsave("./figures/BiasPlots_Century_DefaultFitted_forc.jpeg", 
       height = 5, width = 9) 

# Millennial 
f1_mil <- bias_plot_fun_forc(model = "Millennial", xvar = npp_modis*1000/365,
                             xvar_name = "npp_modis") +
  ggtitle("a) NPP [gC m-2 d-1]")

f2_mil <- bias_plot_fun_forc(model = "Millennial", xvar = Clay_63um,
                             xvar_name = "Clay_63um") +
  ggtitle("b) Clay content < 63 um [%]")

f3_mil <- bias_plot_fun_forc(model = "Millennial", xvar = sm,
                             xvar_name = "sm") +
  ggtitle("c) Soil moisture [m-3 m-3]")

f4_mil <- bias_plot_fun_forc(model = "Millennial", xvar = stemp,
                             xvar_name = "stemp") +
  ggtitle("d) Soil temperature [C]")

f5_mil <- bias_plot_fun_forc(model = "Millennial", xvar = pH,
                             xvar_name = "pH") +
  ggtitle("e) pH")

# Combine plots
annotate_figure(  
  ggarrange(f1_mil, f2_mil, f3_mil, f4_mil, f5_mil, common.legend = TRUE),  
  left = text_grob("Bias in SOC stocks [kg m-2]", rot = 90, size = 13))

ggsave("./figures/BiasPlots_Millennial_DefaultFitted_forc.jpeg", 
       height = 5, width = 9) 

# MIMICS
f1_mic <- bias_plot_fun_forc(model = "MIMICS", xvar = npp_modis*1000,
                             xvar_name = "npp_modis") +
  ggtitle("a) NPP [gC m-2 yr-1]")

f2_mic <- bias_plot_fun_forc(model = "MIMICS", xvar = Clay_2um,
                             xvar_name = "Clay_2um") +
  ggtitle("b) Clay content < 2 um [%]")

f3_mic <- bias_plot_fun_forc(model = "MIMICS", xvar = sm,
                             xvar_name = "sm") +
  ggtitle("c) Soil moisture [m-3 m-3]")

f4_mic <- bias_plot_fun_forc(model = "MIMICS", xvar = stemp,
                             xvar_name = "stemp") +
  ggtitle("d) Soil temperature [C]")

f5_mic <- bias_plot_fun_forc(model = "MIMICS", xvar = LIG_N,
                             xvar_name = "LIG_N") +
  ggtitle("e) Lignin:N ratio")

# Combine plots
annotate_figure(  
  ggarrange(f1_mic, f2_mic, f3_mic, f4_mic, f5_mic, common.legend = TRUE),  
  left = text_grob("Bias in SOC stocks [kg m-2]", rot = 90, size = 13))

ggsave("./figures/BiasPlots_MIMICS_DefaultFitted_forc.jpeg", 
       height = 5, width = 9) 

### soil parameters
# Function to extract linear model statistics
extract_bias_prop_lm_stats <- function(df, xvar_name, type_filter) {
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
bias_plot_fun_prop <- function(xvar, xvar_name){  
  
  # Get types
  types <- model_afsis %>% 
    filter(ModelRun == "fitted") %>% 
    pull(Type) %>% 
    unique()
  
  # Calculate stats for both types
  stats_table <- tibble(Type = types) %>%
    rowwise() %>%
    do(extract_bias_prop_lm_stats(model_afsis, xvar_name, .$Type)) %>%
    ungroup() %>%
    mutate(
      label = sprintf(
        "adj. R² = %.2f\np = %.3e\nRMSE = %.2f",
        adj_r_squared, p_value, rmse
      )
    )
  
  # Create the plot
  model_afsis %>%  
    filter(ModelRun == "fitted") %>%  
    ggplot(aes(y = SOC_bias_kg_m2, x = {{xvar}}, color = Type)) +  
    geom_hline(yintercept = 0) +  
    geom_point(shape = 21, alpha = 0.2) +  
    geom_smooth(aes(fill = Type), method = "lm") +  
    theme_classic(base_size = 13) +  
    theme(axis.text = element_text(color = "black"),  
          axis.title = element_blank(),
          plot.background = element_rect(fill = "transparent", color = "white"),
          legend.background = element_rect(fill = "transparent", color = "white"),
          panel.background = element_rect(fill = "transparent", color = "white"),
          title = element_text(size = 10)) +  
    scale_color_manual("Model", values = c("#1b9e77",  "#d95f02", "#7570b3")) +  
    scale_x_continuous(expand = c(0,0)) +
    scale_fill_manual("Model", values = c("#1b9e77",  "#d95f02", "#7570b3")) +
    geom_text(
      data = stats_table %>% filter(Type == "Century"),
      aes(x = Inf, y = 25, label = label), 
      color = "#1b9e77",
      inherit.aes = FALSE,
      size = 2.5,
      hjust = 3.5
    ) +
    geom_text(
      data = stats_table %>% filter(Type == "Millennial"),
      aes(x = Inf, y = 25, label = label), 
      color = "#d95f02",
      inherit.aes = FALSE,
      size = 2.5,
      hjust = 2.4
    ) +
    geom_text(
      data = stats_table %>% filter(Type == "MIMICS"),
      aes(x = Inf, y = 25, label = label), 
      color = "#7570b3",
      inherit.aes = FALSE,
      size = 2.5,
      hjust = 1.1
    )
}

# Create plots with the updated function
s1 <- bias_plot_fun_prop(Am_Ox_Al/1000, "Am_Ox_Al") +  
  ggtitle("a) Oxalate-extractable Al [g kg-1]") 
s2 <- bias_plot_fun_prop(Am_Ox_Fe/1000, "Am_Ox_Fe") +  
  ggtitle("b) Oxalate-extractable Fe [g kg-1]")  
s3 <- bias_plot_fun_prop(pH, "pH") +  
  ggtitle("c) pH")  
s4 <- bias_plot_fun_prop(Caex, "Caex") +  
  ggtitle("d) exchangeable Ca [cmol kg-1]")  
s5 <- bias_plot_fun_prop(Clay_2um, "Clay_2um") +  
  ggtitle("e) Clay content <2 um [%]")  
s6 <- bias_plot_fun_prop(Clay_63um, "Clay_63um") +  
  ggtitle("f) Clay + fine silt content <63 um [%]")  
s7 <- bias_plot_fun_prop(Clay_1_1, "Clay_1_1") +  
  ggtitle("g) 1:1 clay minerals [%]")  
s8 <- bias_plot_fun_prop(Clay_2_1, "Clay_2_1") +  
  ggtitle("h) 2:1 clay minerals [%]")  

# Combine plots
annotate_figure(  
  ggarrange(s1, s2, s3, s4, s5, s6, s7, s8, common.legend = TRUE,  
            ncol = 4, nrow = 2),  
  left = text_grob("Bias in SOC stocks [kg m-2]", rot = 90, size = 13))

ggsave("./figures/BiasPlots_Century_MIMICS_Millennial_Fitted_soil.jpeg", 
       height = 6, width = 11.5)

# Figure for Cornell talk:
# annotate_figure(  
#   ggarrange(s4, s5, s7, s8, common.legend = TRUE,  
#             ncol = 4, nrow = 1),  
#   left = text_grob("Bias in SOC stocks [kg m-2]", rot = 90, size = 12))
# 
# ggsave("./figures/BiasPlots_Century_MIMICS_Millennial_Fitted_soil_red.jpeg", 
#        height = 4, width = 12)


# Spatial plotting
s9 <- bias_plot_fun_prop(Longitude, "Longitude") +
  ggtitle("a) Longitude")
s10 <- bias_plot_fun_prop(Latitude, "Latitude") +
  ggtitle("b) Latitude")

annotate_figure(  
  ggarrange(s9, s10, common.legend = TRUE),  
  left = text_grob("Bias in SOC stocks [kg m-2]", rot = 90, size = 12))

ggsave("./figures/BiasPlots_Century_MIMICS_Millennial_Fitted_long_lat.jpeg", 
       height = 6, width = 12) 
