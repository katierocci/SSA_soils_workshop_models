##MCMC fitting for MIMICS, adapted from MIMICS_HiRes repo by Derek Pierson

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

########################################
# Load forcing data
########################################
forcing_data <- read.csv("forcing_data/afsis_ref_updated8.csv", as.is=T)
data <- forcing_data %>% 
  drop_na(npp_modis, pH, Clay_2um, LIG_N, stemp, sm, NPP.gC.m2.d, Clay_63um, bd_extracted)  %>% 
  filter(CORG <=20) %>%
  mutate(SOC = CORG*10*bd_extracted*20/100,npp_modis.gC.m2.yr=npp_modis*1000) %>%
  select(X, ANPP=npp_modis.gC.m2.yr,
         CLAY = Clay_2um,
         lig_N = LIG_N,
         TSOI = stemp,
         theta_liq = sm,
         Depth = Depth, SOC) %>%
  filter(Depth=='Topsoil')

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
# Load MIMICS data and ftns from Brute Forcing script
########################################
source("functions/MIMICS_MCMC_funcs/MIMICS_repeat_base.R")

############################################
#determine inital RMSE for training dataset
#############################################
MIMICS_afsis_train <- dataTrain %>% split(1:nrow(dataTrain)) %>% map(~ MIMICS_SP(df=.)) %>% bind_rows()
Train_ModObs <- MIMICS_afsis_train %>% mutate(X=Site, SOC.mod = MICr+MICK+SOMp+SOMc+SOMa) %>% select(X,SOC.mod) %>% inner_join(obsTrain, by="X")
RMSE_train <- sqrt(mean((Train_ModObs$SOC - Train_ModObs$SOC.mod)^2))
lm_train <- lm(Train_ModObs$SOC ~ Train_ModObs$SOC.mod, data = Train_ModObs)
lm_train_sum <- summary(lm_train)
r2_train <- lm_train_sum$adj.r.squared
max_RMSE=100

########################################
# Set allowable min/max range for each MIMICS parameter
## Values are multipliers of default parameters in MIMICS sandbox
########################################

p_rng <- data.frame(Parameter = c("Vslope", "Vint", "Kslope", "Kint", "Tau", "CUE", "desorb", "fPHYS"),
                    P_min = c(0.4, 0.3, 0.4, 0.3, 0.3, 0.2, 0.001, 0.01),
                    P_max = c(4, 3, 4, 3, 3, 1.33, 1, 4))

########################################
# Create dataframe to store MCMC steps
########################################
MCMC_out <- data.frame(i=0,
                       iter=1,
                       Vslope_x=1,
                       Vint_x=1,
                       Kslope_x=1,
                       Kint_x=1,
                       Tau_x=1,
                       CUE_x=1,
                       desorb_x=0.17,
                       fPHYS_x=0.22,
                       slope=0,
                       r2=r2_train,
                       RMSE=RMSE_train,
                       #testing
                       #r2=0.02,
                       #RMSE=100,
                       MICpropSOC=0,
                       LITpropSOC=0,
                       MIM_CO_Avg=0,
                       SOMpTOvAvg=0,
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

#set initial loop values (aka initial piors)
curr_p <- data.frame(Vslope_x = 1,   
                     Vint_x = 1,
                     Kslope_x = 1,
                     Kint_x = 1,
                     Tau_x = 1,
                     CUE_x = 1,
                     desorb_x = 0.17,
                     fPHYS_x = 0.22,
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
MIM_runs <- 10 #originally 200

# Send progress statement to console
print(paste0("Running ", as.character(MIM_runs), " MCMC iterations"))

#Run MCMC loop
for(i in 1:MIM_runs) {
  
  print(paste0("Running proposal set #", as.character(i)))
  
  for(j in 1:8) {
    
    #Set new parameter value
    test_p <- curr_p
    
    #Get random parameters to test, in groups
    if(j == 1) {test_p[1,1] <- runif(1, p_rng[1,2], p_rng[1,3])} #Vslope
    if(j == 2) {test_p[1,2] <- runif(1, p_rng[2,2], p_rng[2,3])} #Vint
    if(j == 3) {test_p[1,3] <- runif(1, p_rng[3,2], p_rng[3,3])} #Kslope
    if(j == 4) {test_p[1,4] <- runif(1, p_rng[4,2], p_rng[4,3])} #Kint
    if(j == 5) {test_p[1,5] <- runif(1, p_rng[5,2], p_rng[5,3])} #Tau
    if(j == 6) {test_p[1,6] <- runif(1, p_rng[6,2], p_rng[6,3])} #CUE 
    if(j == 7) {test_p[1,7] <- runif(1, p_rng[7,2], p_rng[7,3])} #desorb
    if(j == 8) {test_p[1,8] <- runif(1, p_rng[8,2], p_rng[8,3])} #fPHYS
    
    print(test_p)
    
    #Run MIMICS ftn with test parameters
    MIMout <- MIMrepeat(forcing_df = dataTrain, rparams = test_p)
    
    #log parameter updates in dataframe
    iter_out <- data.frame(i=i,
                           iter = ((i-1)*3) + j, 
                           Vslope_x=test_p[1],
                           Vint_x=test_p[2],
                           Kslope_x=test_p[3],
                           Kint_x=test_p[4],
                           Tau_x=test_p[5],
                           CUE_x=test_p[6],
                           desorb_x=test_p[7],
                           fPHYS_x=test_p[8],
                           slope=MIMout$slope,
                           r2=MIMout$r2,
                           RMSE=MIMout$RMSE,
                           MICpropSOC=MIMout$MICpropSOC,
                           LITpropSOC=MIMout$LITpropSOC,
                           MIM_CO_Avg=MIMout$MIM_CO_mn,
                           SOMpTOvAvg=MIMout$SOMpTO,
                           improve=0)
    
    #Make decision based on cost outcome
    if(MIMout$RMSE < curr_RMSE &&
       MIMout$r2 > curr_r2 &&
       #MIMout$RMSE/max_RMSE)-MIMout$r2) < curr_cost &&
       MIMout$MICpropSOC > 0.008 &&
       MIMout$MICpropSOC < 0.08 && #informed by Xu et al., 2013, https://doi.org/10.1111/geb.12029
       MIMout$LITpropSOC > 0.05 &&
       MIMout$LITpropSOC < 0.50 &&
       MIMout$MIM_CO_mn > 0.01 &&
       MIMout$MIM_CO_mn < 100 &&
       MIMout$SOMpTO > 50 &&
       MIMout$SOMpTO < 1000) 
      {
      
      #Update targets
      curr_p <- test_p
      #curr_cost <- (MIMout$RMSE/max_RMSE)-MIMout$r2
      curr_RMSE = MIMout$RMSE
      curr_r2 = MIMout$r2
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
      p_rng[1,2] <- iter_out$Vslope_x / walk_rt # V_slope min
      p_rng[1,3] <- iter_out$Vslope_x +(iter_out$Vslope_x-(iter_out$Vslope_x/walk_rt)) # V_slope max
      
      p_rng[2,2] <- iter_out$Vint_x / walk_rt # V_int min
      p_rng[2,3] <- iter_out$Vint_x +(iter_out$Vint_x-(iter_out$Vint_x/walk_rt)) # V_int max
      
      p_rng[3,2] <- iter_out$Kslope_x / walk_rt # K_slope min
      p_rng[3,3] <- iter_out$Kslope_x +(iter_out$Kslope_x-(iter_out$Kslope_x/walk_rt)) # K_slope max
      
      p_rng[3,2] <- iter_out$Kint_x / walk_rt # K_int min
      p_rng[3,3] <- iter_out$Kint_x +(iter_out$Kint_x-(iter_out$Kint_x/walk_rt)) # K_int max
      
      p_rng[5,2] <- iter_out$Tau_x / walk_rt # Tau min
      p_rng[5,3] <- iter_out$Tau_x +(iter_out$Tau_x-(iter_out$Tau_x/walk_rt)) # Tau max
      
      p_rng[6,2] <- iter_out$CUE_x / walk_rt # CUE min
      p_rng[6,3] <- iter_out$CUE_x +(iter_out$CUE_x-(iter_out$CUE_x/walk_rt)) # CUE max
      
      p_rng[7,2] <- iter_out$desorb_x / walk_rt # desorb min
      p_rng[7,3] <- iter_out$desorb_x +(iter_out$desorb_x-(iter_out$desorb_x/walk_rt)) # desorb max
      
      p_rng[8,2] <- iter_out$fPHYS_x / walk_rt # fPHYS min
      p_rng[8,3] <- iter_out$fPHYS_x +(iter_out$fPHYS_x-(iter_out$fPHYS_x/walk_rt)) # fPHYS max
      
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

pRMSE <- ggplot(MCMC_out, aes(x=iter, y=RMSE)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pr2 <- ggplot(MCMC_out, aes(x=iter, y=r2)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pTau_x <- ggplot(MCMC_out, aes(x=iter, y=Tau_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pCUE_x <-ggplot(MCMC_out, aes(x=iter, y=CUE_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pDesorb_x <- ggplot(MCMC_out, aes(x=iter, y=desorb_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pFPHYS_x <- ggplot(MCMC_out, aes(x=iter, y=fPHYS_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pVslope_x <- ggplot(MCMC_out, aes(x=iter, y=Vslope_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pVint_x <- ggplot(MCMC_out, aes(x=iter, y=Vint_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pKslope_x <- ggplot(MCMC_out, aes(x=iter, y=Kslope_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")
pKint_x <- ggplot(MCMC_out, aes(x=iter, y=Kint_x)) + geom_line(color="grey50", alpha=0.5) + geom_point(size=3, color="grey50", alpha=0.5)  + geom_line(data=MCMC_out %>% filter(improve > 0), color="red", linewidth=1) + geom_point(data=MCMC_out %>% filter(improve > 0), color="red", size=4) + theme_minimal() +theme(legend.position = "none")

walk_plot <- grid.arrange(pRMSE, pr2, pTau_x, pCUE_x, pDesorb_x, pFPHYS_x, pVslope_x, pVint_x, pKslope_x, ncol = 2)

MCMC_out_check <- MCMC_out %>% filter(MIM_CO_Avg<100 & LITpropSOC<1 & RMSE<11) #only one paramter changing for best pset - seems messed up... and other issues with data - might try just r2 and see if that helps?

MCMC_out_improved <- MCMC_out%>%filter(improve == 1) %>%mutate(cost = (RMSE/100)-r2)
hist(MCMC_out_improved$RMSE)
hist(MCMC_out_improved$r2)
MCMC_SingleBest <- MCMC_out %>% filter(iter==24 & i==7)

#save plot
ggsave(file=paste0("MCMC/Output/", format(Sys.time(), "%Y%m%d_%H%M%S_"), "MIM_MCMC_pCombos-", as.character(MIM_runs),"_walk_plot", ".jpeg"), 
       plot=walk_plot,
       width=10,
       height=8)

#######################
# Export MCMC run data
#######################
write.csv(MCMC_out_check, "AfSIS_MCMC_Training_BestRMSEr2combo.csv")
write.csv(MCMC_SingleBest, "AfSIS_MCMC_Training_Best.csv")

######################
#import other good psets

#updated inputs - current best
MCMC_OldBest <- read.csv("AfSIS_MCMC_Training_UpdatedInputs.csv")
MCMC_improved <- MCMC_OldBest %>% filter(improve == 1)
MCMC_SingleBest <- MCMC_OldBest %>% filter(iter==351 & i==117)

#prioritizing r2
MCMC_r2 <- read.csv("AfSIS_MCMC_Training_SlightlyWorseRMSE.csv")
MCMC_out_check_r2 <- MCMC_r2 %>% filter(MIM_CO_Avg<100 & LITpropSOC<1 & RMSE<11)
MCMC_SingleBest <- MCMC_r2 %>% filter(iter==28 & i==9)
write.csv(MCMC_SingleBest, "AfSIS_MIMICS_MCMC_Parameterized.csv")

#original data
#MCMC_SingleBest <- read.csv("AfSIS_MCMC_Training_Best.csv")


##################################
#try parameters in test dataset and with full dataset
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

