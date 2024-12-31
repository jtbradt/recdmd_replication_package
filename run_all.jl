################################################################################
#   SCRIPT:     run_all.jl
#   AUTHOR:     Jacob Bradt (jacob.bradt@mccombs.utexas.edu)
#   NOTES:      julia script to execute all code in the replication package for 
#               "Hotelling Meets Wright: Spatial Sorting and Measurement Error 
#               in Recreation Demand Models."
################################################################################

#   Activate project
using Pkg
Pkg.activate(".")

#   Load RCall
using RCall

#   Load R packages:
R"renv::restore()"  ## Enter "y" when prompted

###################### STEP 1: Simulations ######################

# 1.1: Run simulations
include("src/julia/s_run_simulations.jl");

# 1.2: Process simulation results
R"source('src/R/s_simulation_figs.R')";

########### STEP 2: English et al (2018) Replication ###########

## STEP 2.1: Process data for first stage travel cost instruments

# 2.1.1: Construct structural oil shock time series:
R"source('src/R/e_gen_oil_vars.R')";

# 2.1.2: Process state gas tax instruments:
R"source('src/R/e_gen_gas_taxes.R')";

# 2.1.3: Construct distance matrix for use in constructing first stage cost instruments:
R"source('src/R/e_gen_dist_mat.R')";

# 2.1.4: Construct cost instruments based on CBSA census data:
R"source('src/R/e_gen_census_vars.R')";

# 2.1.5: Construct cost instruments based on commuting time data 
# (**Note: this instrument no longer used in estimates reported in paper**):
R"source('src/R/e_gen_commuting_time.R')";

# 2.1.6: Construct final first stage cost instrument objects for use in second stage 
# estimation and estimate first stage travel cost regressions for reporting in paper:
R"source('src/R/e_construct_final_inst.R')";

## STEP 2.2: Estimation of recreation demand models for English et al (2018) replication

# 2.2.1: Estimate main results in Table 3 of paper:
include("src/julia/e_est_nmnl_main.jl");

# 2.2.2: Estimate robustness checks in Table A3 of paper:
include("src/julia/e_est_nmnl_main_rc.jl");

# 2.2.3: Estimate robustness checks in Table A4 of paper:
include("src/julia/e_est_nmnl_local.jl");
