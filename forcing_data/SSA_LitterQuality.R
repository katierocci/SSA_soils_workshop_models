#code for extracting litter quality data for AfSIS sites

#load libraries
library("ggplot2")
library("dplyr")
library("tidyr")

#load pre-processed TRY data so manageable size for R processing
TRY_filtered <- read.csv("TRY_LigN_filtered_data.csv")

#summarize lignin and N by species
TRY_LN <- TRY_filtered %>% mutate(Trait = case_when(TraitID==70 ~ 'LitN',
                                                    TraitID==73 ~ 'LitLig')) %>%
  select(AccSpeciesName, Trait, OrigValueStr) %>% group_by(AccSpeciesName, Trait) %>%
  summarise(Value = mean(OrigValueStr)) %>%
  pivot_wider(names_from = Trait, values_from = Value) %>% mutate(LitLigN = LitLig/LitN)

#read in plant growth form data and join by species name
TRY_PGF <- read.csv("TRY_Categorical_Traits.csv")
TRY_PGF_filter <- TRY_PGF %>% select(AccSpeciesName, PlantGrowthForm)
TRY_all <- TRY_LN %>% inner_join(TRY_PGF_filter, by="AccSpeciesName")
TRY_all_sum <- TRY_all %>% group_by(PlantGrowthForm) %>%
  summarise(LitLigN_mean = mean(LitLigN, na.rm=TRUE), LitLig_mean = mean(LitLig, na.rm=TRUE), count = n())

#converting summary data into individual values for each plant growth form
#lignin:N
forb_ligN <- as.numeric(TRY_all_sum[4,2]) #no forb so using herb since forbs are herbaceous (non-woody)
gram_ligN <- as.numeric(TRY_all_sum[3,2])
shrub_ligN <- as.numeric(TRY_all_sum[7,2])
tree_ligN <- as.numeric(TRY_all_sum[9,2])
#lignin
forb_lig <- as.numeric(TRY_all_sum[4,3]) #no forb so using herb since forbs are herbaceous (non-woody)
gram_lig <- as.numeric(TRY_all_sum[3,3])
shrub_lig <- as.numeric(TRY_all_sum[7,3])
tree_lig <- as.numeric(TRY_all_sum[9,3])

#load most recent afsis data
afsis <- read.csv("afsis_ref_updated9.csv")
#calculate lignin and lignin:N for each afsis observation based on vegetation structure
afsis_plant_trait <- afsis %>% mutate(VegTot = Shrubs+Graminoids+Forbs+Trees,
                                      PerShrub = Shrubs/VegTot,
                                      PerGram = Graminoids/VegTot,
                                      PerForb = Forbs/VegTot,
                                      PerTree = Trees/VegTot,
                                      LIG_N=case_when(VegStructure=="Forest" ~ tree_ligN,
                                                      VegStructure=="Woodland" ~ tree_ligN*0.4+gram_ligN*0.6,
                                                      VegStructure=="Bushland" ~ shrub_ligN*0.2+tree_ligN*0.2+gram_ligN*0.6,
                                                      VegStructure=="Thicket" ~ shrub_ligN,
                                                      VegStructure=="Shrubland" ~ shrub_ligN,
                                                      VegStructure=="Grassland" ~ gram_ligN*0.5+forb_ligN*0.5,
                                                      VegStructure=="Wooded Grassland" ~ gram_ligN*0.375+forb_ligN*0.375+tree_ligN*0.125+shrub_ligN*0.125,
                                                      VegStructure=="Cropland" ~ forb_ligN,
                                                      VegStructure=="Mangrove" ~ shrub_ligN,
                                                      TRUE ~ PerShrub*shrub_ligN + PerGram*gram_ligN+PerForb*forb_ligN+PerTree*tree_ligN),
                                      LIG=case_when(VegStructure=="Forest" ~ tree_lig,
                                                    VegStructure=="Woodland" ~ tree_lig*0.4+gram_lig*0.6,
                                                    VegStructure=="Bushland" ~ shrub_lig*0.2+tree_lig*0.2+gram_lig*0.6,
                                                    VegStructure=="Thicket" ~ shrub_lig,
                                                    VegStructure=="Shrubland" ~ shrub_lig,
                                                    VegStructure=="Grassland" ~ gram_lig*0.5+forb_lig*0.5,
                                                    VegStructure=="Wooded Grassland" ~ gram_lig*0.375+forb_lig*0.375+tree_lig*0.125+shrub_lig*0.125,
                                                    VegStructure=="Cropland" ~ forb_lig,
                                                    VegStructure=="Mangrove" ~ shrub_lig,
                                                    TRUE ~ PerShrub*shrub_lig + PerGram*gram_lig+PerForb*forb_lig+PerTree*tree_lig)) %>%
  select(SSN, LIG_N, LIG)

#join with afsis data
afsis2 <- afsis %>% inner_join(afsis_plant_trait, by="SSN")
