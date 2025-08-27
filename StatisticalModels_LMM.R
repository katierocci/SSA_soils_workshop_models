## AfSIS statistical model runs ##
## Sophie von Fromm ##
## 2025-08-27 ##

library(tidyverse)
library(ggpubr)
library(nlme)
library(MuMIn)

### Linear mixed-effects models

# Load and prepare AfSIS data
afsis <- read.csv("forcing_data/afsis_ref_updated8.csv", as.is = T)

afsis_data <- afsis %>% 
  filter(Depth == "Topsoil") %>% 
  filter(CORG <= 20) %>% 
  #Calculate SOC stocks for observational data (kg/m2)
  mutate(C_stock = ((CORG/100)*bd_extracted*20)*10) %>% 
  drop_na(npp_modis, pH, Clay_2um, LIG_N, stemp, sm, Clay_63um, bd_extracted,
          C_stock)

afsis_data$Plot <- as.factor(afsis_data$Plot)

afsis_data$Cluster <- as.factor(afsis_data$Cluster)

# Select variables (based on model forcing patameters)
afsis_mimics <- afsis_data %>% 
  mutate(npp_modis.gC.m2.yr = npp_modis*1000) %>% 
  dplyr::select(SSN, Longitude, Latitude, Site, Cluster, Plot, Depth,
                npp_modis.gC.m2.yr, Clay_2um, LIG_N, stemp, sm, C_stock)

afsis_millennial <- afsis_data %>% 
  mutate(npp_modis.gC.m2.d = npp_modis*1000/365) %>% 
  dplyr::select(SSN, Longitude, Latitude, Site, Cluster, Plot, Depth,
                npp_modis.gC.m2.d, Clay_63um, pH, stemp, sm, C_stock)

#### Linear mixed-effects model ####
## MIMICS ##
afsis_mimics_lmm <- afsis_mimics %>% 
  # normalize and scale data
  dplyr::select(SSN, Site, Cluster, Plot,
                npp_modis.gC.m2.yr, Clay_2um, LIG_N, stemp, sm, C_stock) %>% 
  # add small value to avoid 0's for transformation
  dplyr::mutate(npp_modis.gC.m2.yr = npp_modis.gC.m2.yr + 0.0001) %>% 
  dplyr::mutate_if(is.numeric, ~predict(bestNormalize::boxcox(.x)))

lmim0 <- lme(C_stock ~ 1,
            random = ~1|Site/Cluster/Plot, method = "ML", data = afsis_mimics_lmm)

lmim1 <- update(lmim0, ~. + npp_modis.gC.m2.yr)
lmim2 <- update(lmim1, ~. + Clay_2um)
lmim3 <- update(lmim2, ~. + LIG_N)
lmim4 <- update(lmim3, ~. + sm)
lmim5 <- update(lmim4, ~. + stemp)

# no-autocorrelation
as.data.frame(car::vif(lmim5))
as.data.frame(car::vif(lmim5)) %>% 
  filter(car::vif(lmim5) > 3)

# Summary output and diagnostic plot of full model (all parameters)
summary(lmim5)
# Fitted vs residuals
plot(lmim5, main = "Residuals vs Fitted")
# scale-location plot
plot(lmim5, sqrt(abs(resid(.))) ~ fitted(.))
# Q-Q-Plot
qqnorm(lmim5, abline = c(0,1),
       main = "qqnorm Plot")

# For each predictor
E2 <- resid(lmim5, type = "normalized")
F2 <- fitted(lmim5)
plot(x = afsis_mimics_lmm$Clay_2um,
     y = E2, ylab = "Residuals",
     xlab = "Clay_2um")
plot(x = afsis_mimics_lmm$npp_modis.gC.m2.yr,
     y = E2, ylab = "Residuals",
     xlab = "npp_modis.gC.m2.yr")
plot(x = afsis_mimics_lmm$stemp,
     y = E2, ylab = "Residuals",
     xlab = "stemp")
plot(x = afsis_mimics_lmm$sm,
     y = E2, ylab = "Residuals",
     xlab = "sm")
plot(x = afsis_mimics_lmm$LIG_N,
     y = E2, ylab = "Residuals",
     xlab = "LIG_N")

# Anova output for all models (step-wise) and full model
ava_mim <- anova(lmim0, lmim1, lmim2, lmim3, lmim4, lmim5)
ava_mim

lmim5_final <- update(lmim5, method = "REML")
summary(lmim5_final)

#NPP
r.squaredGLMM(lmim1)
NPP <- r.squaredGLMM(lmim1)[1,1]
#Clay
r.squaredGLMM(lmim2)
Clay_2um <- r.squaredGLMM(lmim2)[1,1]-r.squaredGLMM(lmim1)[1,1]
#lig_n
r.squaredGLMM(lmim3)
lig_n <- r.squaredGLMM(lmim3)[1,1]-r.squaredGLMM(lmim2)[1,1]
#sm
r.squaredGLMM(lmim4)
sm <- r.squaredGLMM(lmim4)[1,1]-r.squaredGLMM(lmim3)[1,1]
#stemp
r.squaredGLMM(lmim5)
stemp <- r.squaredGLMM(lmim5)[1,1]-r.squaredGLMM(lmim4)[1,1]

# Plot model results
R2m.mim <- tibble(NPP, Clay_2um, lig_n, sm, stemp) %>% 
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 2)) %>% 
  mutate(Models = "MIMICS (R² = 0.45)") %>% 
  dplyr::select(Models, everything()) %>% 
  mutate(across(NPP:stemp, ~ ifelse(.x < 0, 0, .x))) %>% 
  pivot_longer(!Models, names_to = "predictor", values_to = "values")

R2m.mim$predictor <- factor(R2m.mim$predictor,
                            levels = c("NPP", "Clay_2um", "lig_n", "sm", "stemp"),
                            ordered = TRUE)

p_R2m.mim <- R2m.mim %>% 
  ggplot(aes(y = Models, x = values*100, fill = predictor)) +
  geom_bar(stat = "identity", color = "black", 
           position = position_stack(reverse = TRUE)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  scale_x_continuous("Explained variation (%)",
                     expand = c(0,0), limits = c(0,55)) +
  scale_y_discrete(expand = c(0,0))

## Millennial ##
afsis_millennia_lmm <- afsis_millennial %>% 
  # normalize and scale data
  dplyr::select(SSN, Site, Cluster, Plot,
                npp_modis.gC.m2.d, Clay_63um, pH, stemp, sm, C_stock) %>% 
  # add small value to avoid 0's for transformation
  dplyr::mutate(npp_modis.gC.m2.d = npp_modis.gC.m2.d + 0.0001) %>% 
  dplyr::mutate_if(is.numeric, ~predict(bestNormalize::boxcox(.x)))

lmil0 <- lme(C_stock ~ 1,
             random = ~1|Site/Cluster/Plot, method = "ML", data = afsis_millennia_lmm)

lmil1 <- update(lmil0, ~. + npp_modis.gC.m2.d)
lmil2 <- update(lmil1, ~. + Clay_63um)
lmil3 <- update(lmil2, ~. + pH)
lmil4 <- update(lmil3, ~. + sm)
lmil5 <- update(lmil4, ~. + stemp)

# no-autocorrelation
as.data.frame(car::vif(lmil5))
as.data.frame(car::vif(lmil5)) %>% 
  filter(car::vif(lmil5) > 3)

# Summary output and diagnostic plot of full model (all parameters)
summary(lmil5)
# Fitted vs residuals
plot(lmil5, main = "Residuals vs Fitted")
# scale-location plot
plot(lmil5, sqrt(abs(resid(.))) ~ fitted(.))
# Q-Q-Plot
qqnorm(lmil5, abline = c(0,1),
       main = "qqnorm Plot")

# For each predictor
E2 <- resid(lmil5, type = "normalized")
F2 <- fitted(lmil5)
plot(x = afsis_millennia_lmm$Clay_63um,
     y = E2, ylab = "Residuals",
     xlab = "Clay_63um")
plot(x = afsis_millennia_lmm$npp_modis.gC.m2.d,
     y = E2, ylab = "Residuals",
     xlab = "npp_modis.gC.m2.d")
plot(x = afsis_millennia_lmm$stemp,
     y = E2, ylab = "Residuals",
     xlab = "stemp")
plot(x = afsis_millennia_lmm$sm,
     y = E2, ylab = "Residuals",
     xlab = "sm")
plot(x = afsis_millennia_lmm$pH,
     y = E2, ylab = "Residuals",
     xlab = "pH")

# Anova output for all models (step-wise) and full model
ava_mil <- anova(lmil0, lmil1, lmil2, lmil3, lmil4, lmil5)
ava_mil

lmil5_final <- update(lmil5, method = "REML")
summary(lmil5_final)

#NPP
r.squaredGLMM(lmil1)
NPP <- r.squaredGLMM(lmil1)[1,1]
#Clay
r.squaredGLMM(lmil2)
Clay_63um <- r.squaredGLMM(lmil2)[1,1]-r.squaredGLMM(lmil1)[1,1]
#pH
r.squaredGLMM(lmil3)
pH <- r.squaredGLMM(lmil3)[1,1]-r.squaredGLMM(lmil2)[1,1]
#sm
r.squaredGLMM(lmil4)
sm <- r.squaredGLMM(lmil4)[1,1]-r.squaredGLMM(lmil3)[1,1]
#stemp
r.squaredGLMM(lmil5)
stemp <- r.squaredGLMM(lmil5)[1,1]-r.squaredGLMM(lmil4)[1,1]

# Plot model results
R2m.mil <- tibble(NPP, Clay_63um, pH, sm, stemp) %>% 
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 2)) %>% 
  mutate(Models = "Millennial (R² = 0.54)") %>% 
  dplyr::select(Models, everything()) %>% 
  mutate(across(NPP:stemp, ~ ifelse(.x < 0, 0, .x))) %>% 
  pivot_longer(!Models, names_to = "predictor", values_to = "values")

R2m.mil$predictor <- factor(R2m.mil$predictor,
                            levels = c("NPP", "Clay_63um", "pH", "sm", "stemp"),
                            ordered = TRUE)

p_R2m.mil <- R2m.mil %>% 
  ggplot(aes(y = Models, x = values*100, fill = predictor)) +
  geom_bar(stat = "identity", color = "black", 
           position = position_stack(reverse = TRUE)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_blank()) +
  scale_x_continuous("Explained variation (%)",
                     expand = c(0,0), limits = c(0,55)) +
  scale_y_discrete(expand = c(0,0))

ggarrange(p_R2m.mim, p_R2m.mil, nrow = 2)
ggsave(paste0("./model_output/LMM_MIMICS_Millennial_R2m_",
              Sys.Date(), ".jpeg"), width = 12, height = 6)
