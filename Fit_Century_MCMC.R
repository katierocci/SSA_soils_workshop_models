##MCMC fitting for Century , adapted from MIMICS_HiRes repo by Derek Pierson

## Set working drive
# FYI: Packages are loaded in MIMICS_base_ftn.R

# For local run
#setwd("C:/Users/krocci/OneDrive - UCB-O365/Documents/MIMICS_HiRes")

#################################################
#load libraries needed outside of model code
################################################
library(tidyr)
library(dplyr)
library(caret)
library(here)
library(rootSolve)
library(parallel)
library(future)
library(gridExtra)

########################################
# Load forcing data
########################################
afsis_data <- read.csv("forcing_data/afsis_ref_updated9.csv", as.is=T)

## Prepare forcing data
forcing_data <- afsis_data %>% 
  tibble::rowid_to_column("Set") %>% 
  filter(Depth == 'Topsoil') %>% 
  filter(CORG <= 20)

data <- forcing_data %>%
  drop_na(Litterfall.gC.m2.yr, NPP.gC.m2.d, Clay_2um, Clay_63um, stemp, 
          sm, pH, bd_extracted, LIG_N) %>%
  mutate(npp_modis.gC.m2.d = npp_modis*1000/365, #Century units now same as Millennial
         param_claysilt = Clay_63um/100, #fraction in Century
         forc_st = stemp, 
         forc_sw = sm, 
         forc_npp = npp_modis.gC.m2.d,
         LN = LIG_N, 
         LigFrac = LIG/100,
         SOC = CORG*10*bd_extracted*20/100,
         input_to_metlitter = 0.85 - 0.013 * LN,
         X=SSN,
         site_id = row_number()) %>%
  dplyr::select(X, forc_st, forc_sw, forc_npp, LN, param_claysilt, LigFrac, 
                SOC, input_to_metlitter, site_id)

#partition training data
proptrain = 0.8

#Partitions data into training and test datasets
set.seed(1)
trainIndex <- createDataPartition(data$SOC, p = proptrain, list = FALSE, times = 1)
obsTrain <- data[trainIndex,] %>% select(X,SOC)
obsTest <- data[-trainIndex,] %>% select(X,SOC)
obsFull <- data %>% select(X,SOC)

dataTrain <- data[trainIndex,]
dataTest <- data[-trainIndex,] %>% select(-SOC)


########################################
# Load Century parameters and function
########################################

#load functions
source("functions/run_functions_Century.R")
source("ODEs/derivs_Century.R")

# Read in parameters - Described in Table A1 of Abramoff et al. (2021)
parameters.file <- read.table(here("parameters/soilpara_in_century.txt")) #default
parameters <- as.list(parameters.file$V2)
names(parameters) <- parameters.file$V1

# Function to run model for one row of site parameters
run_model_for_site <- function(site_row, base_params) {
  # Copy full parameter list
  site_params <- as.list(base_params)
  # Override default site-specific value with values from AfSIS
  site_params$input_to_metlitter <- as.numeric(site_row$input_to_metlitter) #doesn;t seem to exist????
  #print(site_params$input_to_metlitter)
  site_params$param_claysilt <- as.numeric(site_row$param_claysilt)
  site_params$LigFrac <- as.numeric(site_row$LigFrac)
  
  # Format input data
  input_for_site <- data.frame(
    forc_st  = rep(as.numeric(site_row$forc_st), 365),
    forc_sw  = rep(as.numeric(site_row$forc_sw), 365),
    forc_npp = rep(as.numeric(site_row$forc_npp), 365),
    stringsAsFactors = FALSE
  )
  
  # Run the model with an error-catching wrapper
  result <- tryCatch({
    Solve_Model_Century(input_for_site, derivs_Century, site_params)
  }, error = function(e) {
    warning(paste("Model failed for site:", site_row$site_id, ":", e$message))
    return(NULL)
  })
  
  if (is.null(result)) return(NULL)
  
  pool_df <- as.data.frame(t(result$y))
  #print(pool_df)
  pool_df$site_id <- site_row$site_id
  return(pool_df)
}

############################################
#determine inital RMSE for training dataset
#############################################
# Apply model across all rows of the site parameter dataframe
results_list_train <- apply(dataTrain, 1, function(row) {
  row <- as.list(row)  # Convert row to list for easy referencing
  run_model_for_site(row, parameters)
})
# Combine all into one dataframe
Century_ss_train <- dplyr::bind_rows(results_list_train)
Century_ss_train$site_id <- as.integer(Century_ss_train$site_id)
Cent_ids <- dataTrain %>% select(X, site_id)
Century_ss_train <- Century_ss_train %>% inner_join(Cent_ids, by="site_id")

#combine with obs and assess trained fit
Train_ModObs <- Century_ss_train %>% mutate(X, SOC.mod = (ACTIVE + SLOW + PASSIVE) / 1000) %>% select(X,SOC.mod) %>% inner_join(obsTrain, by="X") #StrLitter + MetLitter + 
RMSE_train <- sqrt(mean((Train_ModObs$SOC - Train_ModObs$SOC.mod)^2)) #5.13
lm_train <- lm(Train_ModObs$SOC ~ Train_ModObs$SOC.mod, data = Train_ModObs)
lm_train_sum <- summary(lm_train)
r2_train <- lm_train_sum$adj.r.squared #0.13
max_RMSE=100

########################################
# Set allowable min/max range for Century parameters
## Values are multipliers of default parameters for Century parameters
########################################

p_rng <- data.frame(Parameter = c("w1", "w2", "t3", "t4", "k_active", "k_slow", "k_passive", 
                                  "metlitter_to_active", "strlitter_to_active", "strlitter_to_slow"),
                    P_min = c(0.8, 0.8, 0.8, 0.8, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25),
                    P_max = c(1.2, 1.2, 1.2, 1.2, 4, 4, 4, 1.4, 1.4, 1.4))

#defaults
w1_default = parameters.file[1,2]
w2_default = parameters.file[3,2]
t3_default = parameters.file[5,2]
t4_default = parameters.file[6,2]
k_active_default = parameters.file[9,2]
k_slow_default = parameters.file[10,2]
k_passive_default = parameters.file[11,2]
metlitter_to_active_default = parameters.file[18,2]
strlitter_to_active_default = parameters.file[19,2]
strlitter_to_slow_default = parameters.file[20,2]

########################################
# Create dataframe to store MCMC steps
########################################
MCMC_out <- data.frame(i=0,
                       iter=1,
                       w1_x=1,
                       w2_x=1,
                       t3_x=1,
                       t4_x=1,
                       k_active_x=1,
                       k_slow_x=1,
                       k_passive_x=1,
                       metlitter_to_active_x=1,
                       strlitter_to_active_x=1,
                       strlitter_to_slow_x=1,
                       r2=r2_train,
                       RMSE=RMSE_train,
                       improve=1)

########################################
# ### Allow multi-core use (not sure this is helpful for a loop)
# Set number of cores to use
########################################

# Set number of CPU cores to use
nbr_cores <- detectCores(all.tests = FALSE, logical = TRUE)-1

plan(multicore, gc = FALSE, workers = nbr_cores)

# Begin track time required
start_time <- Sys.time()

########################################
### Run MCMC
################################

#set initial loop values (aka initial priors)
curr_p <- data.frame(w1_x=1,
                     w2_x=1,
                     t3_x=1,
                     t4_x=1,
                     k_active_x=1,
                     k_slow_x=1,
                     k_passive_x=1,
                     metlitter_to_active_x=1,
                     strlitter_to_active_x=1,
                     strlitter_to_slow_x=1,
                     run_num=NA)

# Set initial cost value (RMSE value to improve from)
#curr_cost <- (RMSE_train/max_RMSE) - r2_train
curr_RMSE <- RMSE_train
curr_r2 <- r2_train
#testing
#curr_RMSE <- 100
#curr_r2 <- 0.02

#Set trackers
iters_wo_improve = 0

#Set number of iterations for each parameter proposal 
Cen_runs <- 100 #originally 200

# Send progress statement to console
print(paste0("Running ", as.character(Cen_runs), " MCMC iterations"))

parameters_MCMC <- parameters.file
#Run MCMC loop
for(i in 1:Cen_runs) {
  
  print(paste0("Running proposal set #", as.character(i)))
  
  for(j in 1:10) {
    
    #Set new parameter value
    test_p <- curr_p
    #print(test_p)
    
    #Get random parameters to test, in groups
    if(j == 1) {test_p[1,1] <- runif(1, p_rng[1,2], p_rng[1,3])} #w1
    if(j == 2) {test_p[1,2] <- runif(1, p_rng[2,2], p_rng[2,3])} #w2
    if(j == 3) {test_p[1,3] <- runif(1, p_rng[3,2], p_rng[3,3])} #t3
    if(j == 4) {test_p[1,4] <- runif(1, p_rng[4,2], p_rng[4,3])} #t4
    if(j == 5) {test_p[1,5] <- runif(1, p_rng[5,2], p_rng[5,3])} #k_active
    if(j == 6) {test_p[1,6] <- runif(1, p_rng[6,2], p_rng[6,3])} #k_slow
    if(j == 7) {test_p[1,7] <- runif(1, p_rng[7,2], p_rng[7,3])} #k_passive
    if(j == 8) {test_p[1,8] <- runif(1, p_rng[8,2], p_rng[8,3])} #metlitter_to_active
    if(j == 9) {test_p[1,9] <- runif(1, p_rng[9,2], p_rng[9,3])} #strlitter_to_active
    if(j == 10) {test_p[1,10] <- runif(1, p_rng[10,2], p_rng[10,3])} #strlitter_to_slow
    
    #apply random parameters
    parameters_MCMC[1,2] = w1_default * test_p$w1_x[1]
    parameters_MCMC[2,2] = w2_default * test_p$w2_x[1]
    parameters_MCMC[3,2] = t3_default * test_p$t3_x[1]
    parameters_MCMC[4,2] = t4_default * test_p$t4_x[1]
    parameters_MCMC[5,2] = k_active_default * test_p$k_active_x[1]
    parameters_MCMC[6,2] = k_slow_default * test_p$k_slow_x[1]
    parameters_MCMC[7,2] = k_passive_default * test_p$k_passive_x[1]
    parameters_MCMC[8,2] = metlitter_to_active_default * test_p$metlitter_to_active_x[1]
    parameters_MCMC[9,2] = strlitter_to_active_default * test_p$strlitter_to_active_x[1]
    parameters_MCMC[10,2] = strlitter_to_slow_default * test_p$strlitter_to_slow_x[1]
    parameters_list <- as.list(parameters_MCMC$V2)
    names(parameters_list) <- parameters_MCMC$V1
    #print(parameters_MCMC)
    
    #run model
    Cenout_results <- apply(dataTrain, 1, function(row) {
      row <- as.list(row)  # Convert row to list for easy referencing
      run_model_for_site(row, parameters_list)
    })
    # Combine all into one dataframe
    Cenout <- dplyr::bind_rows(Cenout_results)
    Cenout$site_id <- as.integer(Cenout$site_id)
    Cent_ids <- dataTrain %>% select(X, site_id)
    Cenout <- Cenout %>% inner_join(Cent_ids, by="site_id")
    #r2 and RMSE
    Train_ModObs_MCMC <- Cenout %>% mutate(X, SOC.mod = (ACTIVE + SLOW + PASSIVE) / 1000) %>% select(X,SOC.mod) %>% inner_join(obsTrain, by="X") #StrLitter + MetLitter + 
    RMSE_MCMC <- sqrt(mean((Train_ModObs_MCMC$SOC - Train_ModObs_MCMC$SOC.mod)^2)) #5.13
    lm_MCMC <- lm(SOC ~ SOC.mod, data = Train_ModObs_MCMC)
    lm_sum_MCMC <- summary(lm_MCMC)
    r2_MCMC <- lm_sum_MCMC$adj.r.squared #0.13
    #print(RMSE_MCMC)
    #print(r2_MCMC)
    
    #log parameter updates in dataframe
    iter_out <- data.frame(i=i,
                           iter = ((i-1)*3) + j, 
                           w1_x=test_p[1],
                           w2_x=test_p[2],
                           t3_x=test_p[3],
                           t4_x=test_p[4],
                           k_active_x=test_p[5],
                           k_slow_x=test_p[6],
                           k_passive_x=test_p[7],
                           metlitter_to_active_x=test_p[8],
                           strlitter_to_active_x=test_p[9],
                           strlitter_to_slow_x=test_p[10],
                           r2=r2_MCMC,
                           RMSE=RMSE_MCMC,
                           improve=0)
    
    #Make decision based on cost outcome
    if(RMSE_MCMC < curr_RMSE &&
       r2_MCMC > curr_r2) 
    {
      
      #Update targets
      curr_p <- test_p
      #curr_cost <- (MIMout$RMSE/max_RMSE)-MIMout$r2
      curr_RMSE = RMSE_MCMC
      curr_r2 = r2_MCMC 
      iter_out$improve <- 1
      iters_wo_improve <- 0
      
      # Print to console
      print(paste0("IMPROVED RMSE and R2 TO ", round(curr_RMSE,4), " AND ", round(curr_r2,3)))
      
      ## Walk proposal distributions 
      # ONLY USEFUL IF COMPUTATIONAL POWER IS LIMITED, comment out if not
      #######################################################################
      # Set walk rate
      walk_rt = 2 # Set the parameter range min to be the current value divided by
      # this number, and the max to the current value multiplied
      # by this number
      
      # New proposal distributions
      ####################################
      p_rng[1,2] <- iter_out$w1_x / walk_rt # V_slope min
      p_rng[1,3] <- iter_out$w1_x +(iter_out$w1_x-(iter_out$w1_x/walk_rt)) # V_slope max
      
      p_rng[2,2] <- iter_out$w2_x / walk_rt # V_int min
      p_rng[2,3] <- iter_out$w2_x +(iter_out$w2_x-(iter_out$w2_x/walk_rt)) # V_int max
      
      p_rng[3,2] <- iter_out$t3_x / walk_rt # K_slope min
      p_rng[3,3] <- iter_out$t3_x +(iter_out$t3_x-(iter_out$t3_x/walk_rt)) # K_slope max
      
      p_rng[3,2] <- iter_out$t4_x / walk_rt # K_int min
      p_rng[3,3] <- iter_out$t4_x +(iter_out$t4_x-(iter_out$t4_x/walk_rt)) # K_int max
      
      p_rng[5,2] <- iter_out$k_active_x / walk_rt # Tau min
      p_rng[5,3] <- iter_out$k_active_x +(iter_out$k_active_x-(iter_out$k_active/walk_rt)) # Tau max
      
      p_rng[6,2] <- iter_out$k_slow_x / walk_rt # CUE min
      p_rng[6,3] <- iter_out$k_slow_x +(iter_out$k_slow_x-(iter_out$k_slow_x/walk_rt)) # CUE max
      
      p_rng[7,2] <- iter_out$k_passive_x / walk_rt # desorb min
      p_rng[7,3] <- iter_out$k_passive_x +(iter_out$k_passive_x-(iter_out$k_passive_x/walk_rt)) # desorb max
      
      p_rng[8,2] <- iter_out$metlitter_to_active_x / walk_rt # fPHYS min
      p_rng[8,3] <- iter_out$metlitter_to_active_x +(iter_out$metlitter_to_active_x-(iter_out$metlitter_to_active_x/walk_rt)) # fPHYS max
      
      p_rng[9,2] <- iter_out$strlitter_to_active_x / walk_rt # fPHYS min
      p_rng[9,3] <- iter_out$strlitter_to_active_x +(iter_out$strlitter_to_active_x-(iter_out$strlitter_to_active_x/walk_rt)) # fPHYS max
      
      p_rng[10,2] <- iter_out$strlitter_to_slow_x / walk_rt # fPHYS min
      p_rng[10,3] <- iter_out$strlitter_to_slow_x +(iter_out$strlitter_to_slow_x-(iter_out$strlitter_to_slow_x/walk_rt)) # fPHYS max
      
    } else {
      #update tracker for number of iterations without improvement
      iters_wo_improve <- iters_wo_improve + 1
      
      # Slowly tighten or expand distributions when, over many iterations,
      # no RMSE improvement is found
      #######################################################
      # Would such a process be an improvement?
      
    }
    
    # Export MCMC data
    MCMC_out <- rbind(MCMC_out, iter_out)
    
  }
}

#Print time required
Sys.time() - start_time

# Release CPU cores
plan(sequential)

nbrOfWorkers()


#######################
# Plot MCMC walk
#######################

write.csv(MCMC_out, "Century_MCMC_out_112525.csv")

pRMSE <- ggplot(MCMC_out, aes(x=iter, y=RMSE)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pr2 <- ggplot(MCMC_out, aes(x=iter, y=r2)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pTau_x <- ggplot(MCMC_out, aes(x=iter, y=w1_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pCUE_x <-ggplot(MCMC_out, aes(x=iter, y=w2_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pDesorb_x <- ggplot(MCMC_out, aes(x=iter, y=t3_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pFPHYS_x <- ggplot(MCMC_out, aes(x=iter, y=t4_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pVslope_x <- ggplot(MCMC_out, aes(x=iter, y=k_active_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pVint_x <- ggplot(MCMC_out, aes(x=iter, y=k_passive_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pKslope_x <- ggplot(MCMC_out, aes(x=iter, y=k_slow_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pKint_x <- ggplot(MCMC_out, aes(x=iter, y=metlitter_to_active_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pSLA_x <- ggplot(MCMC_out, aes(x=iter, y=strlitter_to_active_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pSLS_x <- ggplot(MCMC_out, aes(x=iter, y=strlitter_to_slow_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")

walk_plot <- grid.arrange(pRMSE, pr2, pTau_x, pCUE_x, pDesorb_x, pFPHYS_x, pVslope_x, pVint_x, pKslope_x, ncol = 2)

MCMC_out_improved <- MCMC_out%>%filter(improve == 1) %>%mutate(cost = (RMSE/100)-r2)
hist(MCMC_out_improved$RMSE)
hist(MCMC_out_improved$r2)
MCMC_SingleBest <- MCMC_out %>% filter(iter==5 & i==1)
write.csv(MCMC_SingleBest, "parameters/Century_MCMC_out_SingleBest_112625.csv")


#filtering later
MCMC_out <- read.csv("Century_MCMC_out_112525.csv")
MCMC_out_filt3 <- MCMC_out %>% filter(w1_x>0.5, w1_x<2,w2_x>0.5, w2_x<2,t3_x>0.5, t3_x<2,t4_x>0.5, t4_x<2,k_active_x>0.5, k_active_x<2,k_slow_x>0.5, k_slow_x<2,k_passive_x>0.5, k_passive_x<2,
                                      metlitter_to_active_x>0.5, metlitter_to_active_x<2,strlitter_to_active_x>0.5, strlitter_to_active_x<2, strlitter_to_slow_x>0.5, strlitter_to_slow_x<2) 
MCMC_SingleBest2 <- MCMC_out %>% filter(iter==5 & i==1)
write.csv(MCMC_SingleBest2, "parameters/Century_MCMC_out_SingleBest2_112625.csv")

#save plot
ggsave(file=paste0("MCMC/Output/", format(Sys.time(), "%Y%m%d_%H%M%S_"), "MIM_MCMC_pCombos-", as.character(MIM_runs),"_walk_plot", ".jpeg"), 
       plot=walk_plot,
       width=10,
       height=8)


##################################
#try parameters in test dataset and with full dataset - still MIMICS, fix!
#################################
#default pset - before running, re-load source data to get orginal paramters back!
MIMICS_afsis_test <- dataTest %>% split(1:nrow(dataTest)) %>% map(~ MIMICS_SP(df=.)) %>% bind_rows()
Test_ModObs <- MIMICS_afsis_test %>% mutate(X=Site, SOC.mod = MICr+MICK+SOMp+SOMc+SOMa) %>% select(X,SOC.mod) %>% inner_join(obsTest, by="X")
RMSE_test <- sqrt(mean((Test_ModObs$SOC - Test_ModObs$SOC.mod)^2)) #4.79
lm_test <- lm(Test_ModObs$SOC ~ Test_ModObs$SOC.mod, data = Test_ModObs)
lm_test_sum <- summary(lm_test)
r2_test <- lm_test_sum$adj.r.squared

#parameterized
Tau_MULT <<- Tau_MULT*MCMC_SingleBest$Tau_x
desorb_MULT <<- desorb_MULT*MCMC_SingleBest$desorb_x
fPHYS_MULT <<- fPHYS_MULT*MCMC_SingleBest$fPHYS_x
Vslope_MULT <<- Vslope_MULT*MCMC_SingleBest$Vslope_x
Vint_MULT <<- Vint_MULT*MCMC_SingleBest$Vint_x
Kslope_MULT <<- Kslope_MULT*MCMC_SingleBest$Kslope_x
Kint_MULT <<- Kint_MULT*MCMC_SingleBest$Kint_x
CUE_MULT <<- CUE_MULT*MCMC_SingleBest$CUE_x
MIMICS_afsis_test_paramed <- dataTest %>% split(1:nrow(dataTest)) %>% map(~ MIMICS_SP(df=.)) %>% bind_rows()
Test_ModObs_paramed <- MIMICS_afsis_test_paramed %>% mutate(X=Site, SOC.mod = MICr+MICK+SOMp+SOMc+SOMa) %>% select(X,SOC.mod) %>% inner_join(obsTest, by="X")
RMSE_test_paramed <- sqrt(mean((Test_ModObs_paramed$SOC - Test_ModObs_paramed$SOC.mod)^2)) #3.69
lm_test_paramed <- lm(Test_ModObs_paramed$SOC ~ Test_ModObs_paramed$SOC.mod, data = Test_ModObs_paramed)
lm_test_sum_paramed <- summary(lm_test_paramed)
r2_test_paramed <- lm_test_sum_paramed$adj.r.squared

#full dataset
MIMICS_afsis_full <- data %>% split(1:nrow(data)) %>% map(~ MIMICS_SP(df=.)) %>% bind_rows()
ModObs_full <- MIMICS_afsis_full %>% mutate(X=Site, SOC.mod = MICr+MICK+SOMp+SOMc+SOMa) %>% select(X,SOC.mod) %>% inner_join(obsFull, by="X")
RMSE_full_paramed <- sqrt(mean((ModObs_full$SOC - ModObs_full$SOC.mod)^2)) #4.62
lm_mimics_paramed <- lm(SOC ~ SOC.mod, data = ModObs_full)
summary(lm_mimics_paramed) #r2=0.066
ggplot(data=ModObs_full, aes(x = SOC.mod, y = SOC)) +
  geom_point() +
  scale_x_continuous("Predicted SOC (kg/m2)", expand = c(0,0)) +
  scale_y_continuous("Observed SOC (kg/m2)", expand = c(0,0)) +
  geom_abline(slope = 1) +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 16)

