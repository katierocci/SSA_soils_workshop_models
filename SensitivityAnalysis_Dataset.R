## AfSIS process-based model - sensitivity analysis preparation ##
## Sophie von Fromm ##
## 2025-10-09 ##

library(tidyverse)
library(terra)
library(sf)
library(tmap)

## Load AfSIS data
afsis_data <- read.csv("forcing_data/afsis_ref_updated9.csv", as.is = T) %>% 
  filter(Depth == "Topsoil")

## Load africa map
# Boundaries for Africa
data("World")
africa_map <- World %>%
  dplyr::filter(continent == "Africa") %>% 
  st_transform(4326)

# Sub-Saharan Africa (based on definition of the United Nations geoscheme for Africa)
ssa_map_sf <- africa_map %>% 
  filter(name != "Algeria",
         name != "Egypt",
         name != "Libya",
         name != "W. Sahara",
         name != "Tunisia",
         name != "Morocco",
         name != "Djibouti",
         name != "Sudan") 

## Extract random sampling locations
# Generate random points within the SSA boundaries
set.seed(42)
random_points <- st_sample(
  ssa_map_sf, 
  size = 150000,  # Number of points to generate
  type = "random"  # Random sampling
)

# Convert to dataframe with longitude and latitude
random_points_df <- st_coordinates(random_points) %>%
  as.data.frame() %>%
  dplyr::rename(Longitude = X, Latitude = Y)

random_points_df %>% 
  ggplot(aes(x = Longitude, y = Latitude)) +
  geom_point(size = 0.1)

## Load spatial products
# soil moisture
sm_dir <- "GlobalData/sm.tif"
sm_raster <- terra::rast(sm_dir)
sm_raster
# plot(sm_raster)

# soil temperature
stemp_dir <- "GlobalData/stemp.tif"
stemp_raster <- terra::rast(stemp_dir)
stemp_raster
# plot(stemp_raster)

# NPP
npp_dir <- "GlobalData/NPP_Africa_5km.tif"
npp_raster <- terra::rast(npp_dir)
npp_raster
# plot(npp_raster)

# bulk density
bd_dir <- "GlobalData/sol_db_od_m_30m_0..20cm_2001..2017_v0.13_wgs84.tif"
bd_raster <- terra::rast(bd_dir)
bd_raster
# plot(bd_raster)

## Extract values from spatial products
random_points_sf <- st_as_sf(random_points_df, coords = c("Longitude", "Latitude"), 
                             crs = 4326)

# Extract all spatial products
random_spatial_points_df <- random_points_df %>%
  mutate(
    sm = terra::extract(sm_raster, vect(random_points_sf))[, 2],
    stemp = terra::extract(stemp_raster, vect(random_points_sf))[, 2],
    npp_modis = terra::extract(npp_raster, vect(random_points_sf))[, 2],
    bd_extracted = terra::extract(bd_raster, vect(random_points_sf))[, 2]
  )

# Check for NAs (points that fall outside raster extents)
summary(random_spatial_points_df)

# Remove rows with any NAs and reduce to 100,000 samples
set.seed(123)
random_spatial_points_NA_df <- random_spatial_points_df %>%
  filter(complete.cases(.)) %>% 
  # Convert bulk density from 10×kg/m³ to g/cm³
  mutate(bd_extracted = bd_extracted / 100) %>% 
  slice_sample(n = 100000)

## Model the remaining variables based on observational relationships
# Prepare observational data (only complete cases)
afsis_model_data <- afsis_data %>%
  filter(CORG <= 20) %>% 
  select(Clay_2um, Clay_63um, pH, LIG_N, LIG, sm, stemp, npp_modis, bd_extracted) %>%
  filter(complete.cases(.))

# Fit models for variables without spatial products
model_clay_2um <- lm(Clay_2um ~ sm + stemp + npp_modis + bd_extracted, 
                     data = afsis_model_data)
model_clay_63um <- lm(Clay_63um ~ sm + stemp + npp_modis + bd_extracted, 
                      data = afsis_model_data)
model_pH <- lm(pH ~ sm + stemp + npp_modis + bd_extracted, 
               data = afsis_model_data)
model_LIG_N <- lm(LIG_N ~ sm + stemp + npp_modis + bd_extracted, 
                  data = afsis_model_data)
model_LIG <- lm(LIG ~ sm + stemp + npp_modis + bd_extracted, 
                data = afsis_model_data)

# Check model summaries
summary(model_clay_2um)
summary(model_clay_63um)
summary(model_pH)
summary(model_LIG_N)
summary(model_LIG)

# Predict values for random points with realistic residual variation
set.seed(42)
pred_random_spatial_points_df <- random_spatial_points_NA_df %>%
  mutate(
    Clay_2um = predict(model_clay_2um, .) + rnorm(n(), 0, sigma(model_clay_2um)),
    Clay_63um = predict(model_clay_63um, .) + rnorm(n(), 0, sigma(model_clay_63um)),
    pH = predict(model_pH, .) + rnorm(n(), 0, sigma(model_pH)),
    LIG_N = predict(model_LIG_N, .) + rnorm(n(), 0, sigma(model_LIG_N)),
    LIG = predict(model_LIG, .) + rnorm(n(), 0, sigma(model_LIG))
  )

# Apply constraints based on observed ranges
pred_point_constr <- pred_random_spatial_points_df %>%
  mutate(
    Clay_2um = pmax(min(afsis_model_data$Clay_2um), 
                    pmin(Clay_2um, max(afsis_model_data$Clay_2um))),
    Clay_63um = pmax(min(afsis_model_data$Clay_63um), 
                     pmin(Clay_63um, max(afsis_model_data$Clay_63um))),
    # Ensure Clay_2um <= Clay_63um
    Clay_2um = pmin(Clay_2um, Clay_63um),
    pH = pmax(min(afsis_model_data$pH), 
              pmin(pH, max(afsis_model_data$pH))),
    LIG_N = pmax(min(afsis_model_data$LIG_N), 
                 pmin(LIG_N, max(afsis_model_data$LIG_N))),
    LIG = pmax(min(afsis_model_data$LIG), 
               pmin(LIG, max(afsis_model_data$LIG)))
  )

summary(pred_point_constr)
summary(afsis_model_data)

## Verify the results look reasonable
# Compare distributions to observational data
par(mfrow = c(2, 5))
for(var in c("Clay_2um", "Clay_63um", "pH", "LIG_N", "LIG", "sm", "stemp", "npp_modis", "bd_extracted")) {
  hist(afsis_model_data[[var]], main = paste("Observed", var), col = "lightblue", breaks = 30)
  hist(pred_point_constr[[var]], main = paste("Simulated", var), col = "lightcoral", breaks = 30)
}


# Check correlations are maintained
cor(afsis_model_data[, c("Clay_2um", "Clay_63um", "pH", "LIG_N", "LIG", "sm", 
                         "stemp", "npp_modis", "bd_extracted")])
cor(pred_point_constr[, c("Clay_2um", "Clay_63um", "pH", "LIG_N", "LIG", "sm", 
                          "stemp", "npp_modis", "bd_extracted")])

write.csv(pred_point_constr, row.names = FALSE,
          paste0("./forcing_data/sensitivity_analysis_simulated_input_data_",
                 Sys.Date(), ".csv"))


