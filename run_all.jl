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

###################### STEP 1: Simulations ######################

# 1.1: Run simulations
include("src/julia/1_run_simulations.jl")

# 1.2: Process simulation results
R"source('src/R/simulations/simulation_figs.R')"