## AfSIS process-based model analysis ##
## Sophie von Fromm ##
## 2025-10-08 ##

library(tidyverse)

## Load default and fitted model runs
default_run <- read.csv("model_output/both_model_results_default_run_2025-10-08.csv") %>% 
  mutate(ModelRun = "default") 
fitted_run <- read.csv("model_output/both_model_results_fitted_run_2025-10-08.csv") %>% 
  mutate(ModelRun = "fitted")

## Load AfSIS data
afsis_data <- read.csv("forcing_data/afsis_ref_updated8.csv", as.is = T) %>% 
  filter(Depth == "Topsoil")

## merge both model datasets
both_runs <- default_run %>% 
  full_join(fitted_run)

both_runs_obs <-  both_runs %>% 
  filter(Type %in% c("MIMICS", "Millennial")) %>%
  rename(Soil_Organic_Carbon_kg_m2_mod = Soil_Organic_Carbon_kg_m2) %>%
  left_join(
    both_runs %>% 
      filter(Type == "Observed") %>% 
      select(Set, ModelRun, Soil_Organic_Carbon_kg_m2) %>% 
      rename(Soil_Organic_Carbon_kg_m2_obs = Soil_Organic_Carbon_kg_m2),
    by = c("Set", "ModelRun")
  ) %>%
  select(Set, Type, ModelRun, Soil_Organic_Carbon_kg_m2_mod, Soil_Organic_Carbon_kg_m2_obs)

## Plot data and calculate model performance
# Function that returns a tibble with desired statistics for specific subset of data
extract_lm_stats <- function(df, type_filter, modelrun_filter) {
  data_sub <- df %>% filter(Type == type_filter, ModelRun == modelrun_filter)
  
  lm_fit <- lm(Soil_Organic_Carbon_kg_m2_obs ~ Soil_Organic_Carbon_kg_m2_mod, 
               data = data_sub)
  
  summary_lm <- summary(lm_fit)
  glance_lm <- broom::glance(lm_fit)
  
  # RMSE calculation
  rmse <- sqrt(mean(lm_fit$residuals^2))
  
  # p-value of slope (coefficient on Soil_Organic_Carbon_kg_m2_obs)
  # It's the 2nd coefficient
  p_value <- coef(summary_lm)[2, "Pr(>|t|)"]
  
  tibble(
    Type = type_filter,
    ModelRun = modelrun_filter,
    adj_r_squared = summary_lm$adj.r.squared,
    p_value = p_value,
    rmse = rmse
  )
}


lm_fit <- lm(Soil_Organic_Carbon_kg_m2_obs ~ Soil_Organic_Carbon_kg_m2_mod, 
             data = both_runs_obs %>% 
               filter(Type == "MIMICS",
                      ModelRun == "fitted"))

summary(lm_fit)
sqrt(mean(lm_fit$residuals^2))

types <- unique(both_runs_obs$Type)
modelruns <- unique(both_runs_obs$ModelRun)

stats_table <- tidyr::crossing(Type = types, ModelRun = modelruns) %>%
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
ggsave(paste0("./model_output/ModelFits_MIMICS_Millennial_DefaultFitted_",
              Sys.Date(), ".jpeg"), height = 6, width = 9)  



