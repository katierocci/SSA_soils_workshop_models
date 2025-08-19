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
aV      <-rep(0.000008, 6)
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
beta    <- as.numeric(c(1.5, 1.5))# only used if tauMethod='beta', same value for both r and K
fPHYS_r <- c(0.3, 1.3)
fPHYS_K <- c(0.2, 0.8)
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(1.5e-5, -1.5)
PHYS_scalar <- c(2, -2, NA, NA, NA, NA)
FI      <- c(0.05, 0.05)
fmet_p <- c(1, 0.85, 0.013)

depth <- 20 # set soil depth
h2y        <- 24*365
MICROtoECO <- depth * 1e4 * 1e-3  # mgC/cm3 to g/m2

# fW coeeficient for Pierson-CORPSE moisture control
fW_p1 <- 1.212580 #* 0.6867031  # MSBio new
fW_p2 <- 2.748028 #* 0.6300376  # MSBio new

#Set default methods
fWmethod=1      #0-> fW=1, 1->CORPSE, 2->Calibrated, 3->water scalar from other model
historic=FALSE   #modify Vmax based on historic MAT
fixed_fMET=FALSE #calculate fMET based on litter chemistry
tauMethod='beta'  #'NPP' and 'beta' accepted
desorbMethod='clay' #'clay' and 'metal' accepted - not extensively tested!
psMethod='clay' #'clay' and 'metal' accepted - not extensively tested!

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
MIMICS_SP <- function(df){
  
  
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
  # Set lig:N value if not given
  if (is.na(LIG_N)) {
    LIG_N <- (LIG/100)/(1/(CN/2.5))
  }
  
  # Set moisture control on kinetics (fW)
  if (fWmethod==0) {
    fW=1
  } else if (fWmethod==1) {
    
    air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
    fW = (theta_liq^3 * air_filled_porosity^2.5)/0.022600567942709
    fW = max(0.05, fW)
    
  } else if (fWmethod==2) {
    
    air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
    f <- function(x, p1, p2) {x^p1 * (1-x)^p2}
    fW_p3 <- optimize(f, interval=c(0.01,1), p1=fW_p1, p2=fW_p2, maximum = T)$objective
    fW = (theta_liq^fW_p1 * air_filled_porosity^fW_p2)/fW_p3
    fW = max(0.05, fW) 
    
  } else if (fWmethod==3) {
    fW=W_SCALAR
  }
  
  # For historic MAT dependent kinetics
  if (historic==TRUE) {
    Vslope = Vslope + (MAT*0.00104)
    Vint = Vint - (MAT*0.0228)
  }
  
  
  # set fMET
  if(!fixed_fMET){
    #fMET  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * LIG_N)  
    fMET  <- ifelse(LIG_N<65, (fmet_p[1] * (fmet_p[2] - fmet_p[3] * LIG_N)), 0.01) #to account for high LIG_N that cause negtaives
  } else {
    # Fixed Fmet, ONLY for model validation against MIMICS et al. 2015 Sandbox 
    fMET <- 0.3846423
  }
  
  # Calc litter input rate from annual or daily flux then convert units
    EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4
  
  # ------------ calculate time varying parameters ---------------
  Vmax     <- exp(TSOI * (Vslope * Vslope_MULT) + (Vint * Vint_MULT)) * aV * fW   #<-- Moisture scalar applied
  Km       <- exp(TSOI * (Kslope * Kslope_MULT) + (Kint * Kint_MULT)) * aK
  
  if (tauMethod == 'NPP') {
    Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
    Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
    Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
    beta=c(1,1) # turns off beta
  } else if (tauMethod=='beta') {
    Tau_MOD1 <- 1 # turns off NPP effects on turnover
    beta=beta #[1] # use Kat's density dependent function  
    #add in another parameter that adds 1 to 
  }
  
  Tau_MOD2 <- Tau_MOD[4]                        
  
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
           tau_K[1]*exp(tau_K[2]*fMET))   
  tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT
  #print(tau)
  
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
  fAVAI    <- 1 - (fPHYS + fCHEM)
  
  if (desorbMethod == 'clay') {
    desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))
  } else if (desorbMethod=='metal') {
    fMETAL <- ifelse(AlFe<6367, AlFe/6366.6, 1)
    desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fMETAL))
  }
  
  desorb   <- desorb * desorb_MULT
  fPHYS    <- fPHYS * fPHYS_MULT
  
  if (psMethod == 'clay') {
    pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  } else if (psMethod=='metal') {
    pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fMETAL)))  #Scalar for texture effects on SOMp
  }
  
  CUE = CUE * CUE_MULT
  
  v_MOD    <- vMOD  
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD * vMOD_MULT
  KM       <- Km / (k_MOD * kMOD_MULT)
  
  I       <- array(NA, dim=2)             
  I[1]    <- (EST_LIT / depth) * fMET     
  I[2]    <- (EST_LIT / depth) * (1-fMET)
  Inputs  <- I
  
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- c(NA,NA,NA,NA,NA,NA)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  Tpars <- list(I = I, VMAX = VMAX, KM = KM, CUE = CUE,
                fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO, 
                beta=beta)
  
  # Create arrays to hold output
  lit     <- Tpars$I
  mic     <- Tpars$I
  som     <- rep(NA, 3)
  som[1]  <- Tpars$I[1]
  som[2]  <- Tpars$I[2]
  som[3]  <- Tpars$I[1]
  CO2     <- rep(0, 2)
  
  Ty    <- c( LIT_1 = lit[1], LIT_2 = lit[2], 
              MIC_1 = mic[1], MIC_2 = mic[2], 
              SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3])
  
  ## Set global parameters to allow pass to stode function
  .GlobalEnv$VMAX <- VMAX
  .GlobalEnv$KM <- KM
  .GlobalEnv$fPHYS <- fPHYS
  .GlobalEnv$fCHEM <- fCHEM
  .GlobalEnv$fAVAI <- fAVAI
  .GlobalEnv$I <- I
  .GlobalEnv$tau <- tau
  .GlobalEnv$LITmin <- LITmin
  .GlobalEnv$SOMmin <- SOMmin
  .GlobalEnv$MICtrn <- MICtrn
  .GlobalEnv$desorb <- desorb
  .GlobalEnv$DEsorb <- DEsorb
  .GlobalEnv$OXIDAT <- OXIDAT
  
  
  # Using jitter
  test  <- stode_jitter(stode_y = Ty, stode_time = 1e6, stode_fun = RXEQ, stode_parms = Tpars, stode_pos = TRUE)
  
  #print(test[[2]])
  
  # Not using jitter
  #test  <- stode(y = Ty, time = 1e6, fun = RXEQ, parms = Tpars, positive = TRUE)
  
  
  ### Calc and get MIMICS output 
  MIMLIT    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6 #convert kgC/m2 from mgC/cm3 (0-30 cm) 
  MIMMIC    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
  MIM_CO    <-  test[[1]][[3]]/test[[1]][[4]]
  MIMSOC    <- sum(test[[1]])  * depth *1e4 / 1e6   
  
  table <- as.numeric(test[[1]])
  
  #print(df)
  
  MIMout <- data.frame(Site = df$X,
                       fCLAY = fCLAY,
                       TSOI = TSOI,
                       ANPP = ANPP,
                       LIGN = df$LIG_N,
                       EST_LIT = EST_LIT,
                       MIMSOC = MIMSOC,
                       MIMMIC = MIMMIC,
                       MIMLIT = MIMLIT,
                       MIM_CO = MIM_CO,
                       desorb = as.numeric(desorb),
                       SOMpTOv = 1/(as.numeric(desorb)*24*365), #convert from per hr to per yr
                       LITm = table[1] * depth *1e4 / 1e6, #convert kgC/m2 from mgC/cm3 (0-30 cm) 
                       LITs = table[2] * depth *1e4 / 1e6,
                       MICr = table[3] * depth *1e4 / 1e6,
                       MICK = table[4] * depth *1e4 / 1e6,
                       SOMp = table[5] * depth *1e4 / 1e6,
                       SOMc = table[6] * depth *1e4 / 1e6,
                       SOMa = table[7] * depth *1e4 / 1e6,
                       JITn = test[[2]],
                       DEBUG = note
  )
  #Reset global parameters from last run
  # .GlobalEnv$VMAX <- NA
  # .GlobalEnv$KM <- NA
  # .GlobalEnv$fPHYS <- NA
  # .GlobalEnv$fCHEM <- NA
  # .GlobalEnv$fAVAI <- NA
  # .GlobalEnv$I <- NA
  # .GlobalEnv$tau <- NA
  # .GlobalEnv$LITmin <- NA
  # .GlobalEnv$SOMmin <- NA
  # .GlobalEnv$MICtrn <- NA
  # .GlobalEnv$desorb <- NA
  # .GlobalEnv$DEsorb <- NA
  # .GlobalEnv$OXIDAT <- NA
  
  #remove global variables set for stode ftn
  #rm(I, VMAX, KM, fPHYS, fCHEM, fAVAI, tau, LITmin, SOMmin, MICtrn, desorb, DEsorb, OXIDAT)
  
  return(MIMout)
}

