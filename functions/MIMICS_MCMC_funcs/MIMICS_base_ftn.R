## Set working drive
#setwd("C:/github/MIMICS_HiRes")

#Libraries
library(rootSolve)
library(boot)
library(ggplot2)
library(tidyverse)
library(Metrics) 
library(parallel)
library(furrr)
library(purrr)
library(grid)
library(gridExtra)

#bring in RXEQ function
source("functions/MIMICS_MCMC_funcs/RXEQ_ftn.R")

### "Jitter step"
#Ftn to re-run stode ftn if solve is unsuccessful (mostly necessary to counter MICk crash to 0)
stode_jitter <- function(stode_y = Ty, stode_time = 1e6, stode_fun = RXEQ, stode_parms = Tpars, stode_pos = TRUE, run_i = 0) {
  success <- FALSE
  while (!success) {
    run_i <- run_i + 1
    #do something
    test  <- quiet(stode(y = stode_y, time = stode_time, fun = stode_fun, parms = stode_parms, positive = stode_pos)) #Suppress: "diagonal element is zero"
    tbl <- as.numeric(test[[1]])
    
    # Repeat stode ftn if the r or K microbial pools crash below 1e-10
    success <- tbl[3] > 1e-10 & tbl[4] > 1e-10
    
    if(!success) {
      #Add 1% on to Ty$mic2 if no success
      stode_y['MIC_2'] = stode_y['MIC_2'] * 1.01
      #print(stode_y['MIC_2'])
    }
    
    if(run_i > 5) {
      success <- TRUE
      #print(paste0("Jitter cap: ", as.character(run_i)))
    }
  }
  
  return(c(test,run_i))
}


########################################
# Set MIMICS default parameters
########################################
Vslope  <- rep(0.063, 6)
Vint    <- rep(5.47, 6)
aV      <- rep(0.000008, 6)  
Kslope  <- rep(c(0.025, 0.035, 0.025),2)
Kint    <- rep(3.19, 6)
aK      <- rep(10, 6)
vMOD    <- c(10, 2, 10, 3, 3, 2)
kMOD    <- c(8, 2, 4, 2, 4, 6)
KO      <- c(6, 6)
CUE     <- c(0.55, 0.25, 0.75, 0.35)
tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
Tau_MOD <- c(100, 0.8, 1.2, 2)
Tau_MULT <- 1
fPHYS_r <- c(0.3, 1.3)
fPHYS_K <- c(0.2, 0.8)
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(0.000015, -1.5)
PHYS_scalar <- c(2, -2, NA, NA, NA, NA)
FI      <- c(0.05, 0.05)
fmet_p <- c(1, 0.85, 0.013)
depth <- 30 ###
h2y        <- 24*365
MICROtoECO <- depth * 1e4 * 1e-3  # mgC/cm3 to g/m2

#Set default multipliers
Tau_MULT = 1
desorb_MULT = 1
fPHYS_MULT = 1


########################################
# Apply parameter multipliers
########################################
# Vslope = Vslope * 1.693578
# Vint = Vint * 0.633318
# Kslope = Kslope * 1.782366
# Kint = Kint * 0.3609913
# CUE = CUE * 1
# Tau_MULT = 1
# desorb_MULT = 2.3635554
# fPHYS_MULT = 2.0716163

###########################################
# MIMICS single point function
###########################################
MIMICS_SS <- function(df){
  
  
  # Convert all column names to upper case
  colnames(df) <- toupper(colnames(df))
  
  #print(df)
  
  ### Setup a var to collect run notes
  note <- ""
  
  ### Bring in forcing ANPP value, convert gDW to gC
  ANPP <- df$ANPP/2
  
  ### Bring in CLAY value, convert from percent to decimal
  fCLAY <- df$CLAY/100
  
  ### Bring in metals value
  AlFe <- df$ALFE
  
  ### Bring in TSOI value
  TSOI <- df$TSOI
  
  ### Bring in lig:N forcing data
  LIG_N <- df$LIG_N
  LIG <- df$LIG
  CN <- df$CN
  
  # Bring in soil moisture information
  #if using GWC
  #theta_liq  <- df$GWC/100  # GWC = Gravimetric water content
  #theta_frzn <- 0           # Not used here. TODO Needs validation.
  #if using volumetric
  theta_liq  <- df$THETA_LIQ
  theta_frzn <- df$THETA_FRZN
  
  ### Bring in mean annual temperature data
  MAT <- df$MAT
  
  # Use TSOI if 'MAT' column is not in forcing data
  if(is.null(MAT)){
    MAT <- TSOI
  }
  
  #Bring in W_SCALAR if present
  W_SCALAR = df$W_SCALAR
  
  ############################################################
  # MIMICS MODEL CODE STARTS HERE
  ############################################################
  # function calculates fMET with LIG_N if provided in input data.
  Tpars <- calc_Tpars_Conly(ANPP=ANPP, fCLAY=fCLAY, TSOI=TSOI, MAT=MAT,
                            CN=CN, LIG=LIG, LIG_N=LIG_N, AlFe=AlFe,
                            theta_liq=theta_liq, theta_frzn=theta_frzn, W_SCALAR = W_SCALAR)
  
  #debug
  #print(Tpars$VMAX)
  
  # Create arrays to hold output
  lit     <- Tpars$I
  mic     <- Tpars$I
  som     <- rep(NA, 3)
  som[1]  <- Tpars$I[1]
  som[2]  <- Tpars$I[2]
  som[3]  <- Tpars$I[1]
  CO2     <- rep(0, 2)
  
  Ty    <- c(LIT_1 = lit[1], LIT_2 = lit[2],
             MIC_1 = mic[1], MIC_2 = mic[2],
             SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3])
  
  #debug
  #print(length(Ty))
  #str(Ty)
  #print(is.numeric(Ty))
  
  # ## Set global parameters to ensure passing of variables to stode function
  # .GlobalEnv$VMAX <- Tpars$VMAX
  # .GlobalEnv$KM <- Tpars$KM
  # .GlobalEnv$fPHYS <- Tpars$fPHYS
  # .GlobalEnv$fCHEM <- Tpars$fCHEM
  # .GlobalEnv$fAVAI <- Tpars$fAVAI
  # .GlobalEnv$I <- Tpars$I
  # .GlobalEnv$tau <- Tpars$tau
  # .GlobalEnv$LITmin <- Tpars$LITmin
  # .GlobalEnv$SOMmin <- Tpars$SOMmin
  # .GlobalEnv$MICtrn <- Tpars$MICtrn
  # .GlobalEnv$desorb <- Tpars$desorb
  # .GlobalEnv$DEsorb <- Tpars$DEsorb
  # .GlobalEnv$OXIDAT <- Tpars$OXIDAT
  # .GlobalEnv$beta <- Tpars$beta
  #
  # # ------------RUN THE MODEL-------------
  
  test  <- stode(y = Ty, time = 1e7, fun = RXEQ, parms = Tpars, positive = TRUE)
  
  
  ### Calc and get MIMICS output
  MIMLIT    <- (test[[1]][[1]]+test[[1]][[2]]) * depth *1e4 / 1e6 #convert kgC/m2 from mgC/cm3 (0-30 cm)
  MIMMIC    <- (test[[1]][[3]]+test[[1]][[4]]) * depth *1e4 / 1e6
  MIM_CO    <-  test[[1]][[3]]/test[[1]][[4]]
  MIMSOC    <- sum(test[[1]]) * depth * 1e4 / 1e6
  
  table <- as.numeric(test[[1]])
  
  MIMout <- list()
  MIMout[[1]] <- as.numeric(test[[1]]) * depth *1e4 / 1e6   # MIMICS steady state C pools
  MIMout[[2]] <- Tpars                                      # MIMICS parameters used for/from simulation
  MIMout[[3]] <- df                                         # Forcing dataframe
  MIMout[[4]] <- as.numeric(test[[2]]) * depth *1e4 / 1e6   # CO2-r & CO2-k pools
  
  return(MIMout)
}


#####################
# Example use of 
#####################

# ##############################################
# #single point run
# ##############################################
# data <- data.frame(Site = 1,
#                    pGPP = 1.39,
#                    TSOI = 10.6,
#                    CLAY = 20,
#                    lig_N = 11)
# 
# MIMout_single <- MIMICS1(data[1,])
# 
# 
# ##############################################
# # Full forcing dataset run
# ##############################################
# data <- data <- read.csv("RCrk_Modelling_Data/LTER_SITE_1.csv", as.is=T)
# 
# MIMout_single <- MIMICS1(data[1,])
# 
# MIMrun <- data %>% split(1:nrow(data)) %>% map(~ MIMICS1(df=.)) %>% bind_rows()
# MIMrun <- data %>% cbind(MIMrun %>% select(-Site, -TSOI))
# 
# 
# ################################################
# # Plot SOC vs MIMSOC
# ################################################
# library(ggplot2)
# library(Metrics)
# 
# plot_data <- MIMrun
# 
# #calc SOMp turnover time
# plot_data$desorb_yr <- plot_data$desorb*24*365
# plot_data$SOMpTO <- plot_data$SOMp/plot_data$desorb_yr
# 
# r2_test <- cor.test(MIMrun$SOC, MIMrun$MIMSOC)
# r_val <- round(as.numeric(unlist(r2_test ['estimate'])),2)
# lb2 <- paste("R^2 == ", r_val)
# 
# rmse <- round(rmse(MIMrun$SOC, MIMrun$MIMSOC),2)
# 
# ggplot(plot_data, aes(x=MIMSOC, y=SOC, color=ANPP)) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
#   geom_point(size=4, alpha=0.8) +
#   geom_text(aes(label=paste0(Site)),hjust=-0.2, vjust=0.2) +
#   annotate("text", label = lb2, x = 2, y = 8.5, size = 6, colour = "black", parse=T) +
#   annotate("text", label = paste0("RMSE = ", rmse), x = 2, y = 7.4, size = 6, colour = "black") +
#   ylim(0,10) + xlim(0,10)