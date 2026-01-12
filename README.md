Code for three soil biogeochemistry models, Millennial (https://github.com/rabramoff/Millennial), MIMICS (https://github.com/katierocci/MIMICS_STODE/tree/MIMICS-tutorial), and Century (https://github.com/rabramoff/Millennial), in same directory structure to faciliate use with AfSIS dataset.


The results and findings of this analysis are presented in von Fromm SF, Rocci KS, Anuo CO, Asabere SB, Kanyiri J, Kengdo SK, Nketa KA, Mureva A, Zhang L, Abramoff RZ (2026) Transferability of process-based soil carbon models across pedological domains: lessons from sub-Saharan Africa, JGR Biogeosciences, in preparation. 


Run_Models_AfSIS.Rmd will run all three models with the AfSIS data and create some comparison plots. 

Model_Runs_Analysis.R will compare the model results of the three models for the default and fitted model versions, as well as calculating and plotting the SOC bias for each model.

StatisticalModels_RF.R will compute all random forest models and create some plots for model evaluation and interpretation.

Fit_Century_MCMC.R, Fit_MIMICS_MCMS.R, and Fit_Millennial.Rmd will run the fitted model versions.


The "ODEs" folder has the differential equations that make up Millennial (derivs_V2_MM.R), MIMICS (RXEQ.R), and Century (derivs_Century.R). 

The "parameters" folder has parameters for the default model versions: Millennial (soilpara_in_fit_with_qmax.txt), MIMICS (MIMICS_parameters_sandbox_20231129.R), and Century (soilpara_in_Century.txt), and the fitted model versions: Millennial 
(soilpara_in_fit_SSA_Millennial), MIMICS (AfSIS_MIMICS_MCMC_Parameterized.csv), and Century (Century_MCMC_out_SingleBest_120325.csv). 

The "functions" folder has other functions that are needed to run Millennial, MIMICS, and Century. For Millennial, run_functions.R defines functions for running Millennial in a transient way and to steady state. For MIMICS, calc_Tpars.R calculates final parameters to be used in RXEQ.R and MIMICS_calc_steady_state_pools.R determines steady state values. For Century, run_functions_Century.R defines functions for running Millennial in a transient way and to steady state. 

The "forcing_data" folder has the AfSIS data (afsis_ref_updated9.csv), and the code (SSA_LitterQuality.R) and data (TRY_Categorical_Traits.csv and TRY_LigN_filtered_data.csv) to derive litter quality data for each samples based on field data observations and the TRY database.
The 'model_output' folder contains the model outputs for all three models for the default and fitted model versions.

The folder 'figures' contains all the figures shown in the manuscript (including SI).

