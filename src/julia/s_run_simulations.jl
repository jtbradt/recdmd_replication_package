################################################################################
#   SCRIPT:     1_run_simulations.jl
#   AUTHOR:     Jacob Bradt (jacob.bradt@mccombs.utexas.edu)
#   NOTES:      julia script to run simulations and generate data for the 
#               simulation results reported in "Hotelling Meets Wright: Spatial 
#               Sorting and Measurement Error in Recreation Demand Models."  
#               *Please note* that I set the number of simulations to 100 for 
#               the sake of time in replication. For the paper, I use 1000 
#               simulations. As a result, the results in the paper may differ
#               slightly from the results in this replication.
################################################################################

#   Import dependencies:
using DataFrames, Plots, TexTables, CSV, StatsPlots

#   Import user-written functions:
include("s_simulations.jl")

#   STEP 1: Run Simulations 1-6 -------------

#   Define constants for simulations 1-6:
nInd = 10_000;                             # Number of individuals
nSim = 1_00;                               # Number of simulations **NOTE: This is set to 100 for the sake of time in replication. For the paper, I use 1000 simulations**
nAlt = 20;                                 # Number of alternatives
nGroup = 10;                               # Number of groups of consumers
f_ξ = Normal(0.0, 1.0);                    # Distribution of idiosyncratic unobservables
rng_int = 82220;                           # Random number generator seed

#   Simulation 1: No correlation between costs and idiosyncratic unobservable:
println("Starting Simulation 1")
sim1 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, 0.0, 0.5], n_groups = nGroup, rng = MersenneTwister(rng_int));
CSV.write(
    "data/simulations/est_wtp_sim1.csv",
    vcat(
        DataFrame(sim = "sim1", est = "base", wtp = sim1[1][2,:] ./ sim1[1][1, :]),
        DataFrame(sim = "sim1", est = "groupfe", wtp = sim1[2][2,:] ./ sim1[2][1, :]),
        DataFrame(sim = "sim1", est = "cf", wtp = sim1[3][2,:] ./ sim1[3][1, :])
    ),
    writeheader = true
);
rng_int += 1;

#   Simulation 2: Sorting towards high idiosynratic preference:
println("Starting Simulation 2")
sim2 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups = nGroup, rng = MersenneTwister(rng_int));
CSV.write(
    "data/simulations/est_wtp_sim2.csv",
    vcat(
        DataFrame(sim = "sim2", est = "base", wtp = sim2[1][2,:] ./ sim2[1][1, :]),
        DataFrame(sim = "sim2", est = "groupfe", wtp = sim2[2][2,:] ./ sim2[2][1, :]),
        DataFrame(sim = "sim2", est = "cf", wtp = sim2[3][2,:] ./ sim2[3][1, :])
    ),
    writeheader = true
);
rng_int += 1;

#   Simulation 3: Sorting away from high idiosynratic preference:
println("Starting Simulation 3")
sim3 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, 1.0, 0.5], n_groups = nGroup, rng = MersenneTwister(rng_int))
CSV.write(
    "data/simulations/est_wtp_sim3.csv",
    vcat(
        DataFrame(sim = "sim3", est = "base", wtp = sim3[1][2,:] ./ sim3[1][1, :]),
        DataFrame(sim = "sim3", est = "groupfe", wtp = sim3[2][2,:] ./ sim3[2][1, :]),
        DataFrame(sim = "sim3", est = "cf", wtp = sim3[3][2,:] ./ sim3[3][1, :])
    ),
    writeheader = true
);
rng_int += 1;

#   Simulation 4: Generic additive measurement error:
println("Starting Simulation 4")
sim4 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 0.0], [5.0, 1.0, 0.0, 0.5], measurement_error = "additive", n_groups = nGroup, rng = MersenneTwister(rng_int));
CSV.write(
    "data/simulations/est_wtp_sim4.csv",
    vcat(
        DataFrame(sim = "sim4", est = "base", wtp = sim4[1][2,:] ./ sim4[1][1, :]),
        DataFrame(sim = "sim4", est = "groupfe", wtp = sim4[2][2,:] ./ sim4[2][1, :]),
        DataFrame(sim = "sim4", est = "cf", wtp = sim4[3][2,:] ./ sim4[3][1, :])
    ),
    writeheader = true
);
rng_int += 1;

#   Simulation 5: Multiplicative measurement error
println("Starting Simulation 5")
sim5 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 0.0], [5.0, 1.0, 0.0, 0.5], measurement_error = "mult", n_groups = nGroup, rng = MersenneTwister(rng_int));
CSV.write(
    "data/simulations/est_wtp_sim5.csv",
    vcat(
        DataFrame(sim = "sim5", est = "base", wtp = sim5[1][2,:] ./ sim5[1][1, :]),
        DataFrame(sim = "sim5", est = "groupfe", wtp = sim5[2][2,:] ./ sim5[2][1, :]),
        DataFrame(sim = "sim5", est = "cf", wtp = sim5[3][2,:] ./ sim5[3][1, :])
    ),
    writeheader = true
);
rng_int += 1;

#   Simulation 6: Generic additive measurement error + sorting towards high idiosynratic preference:
println("Starting Simulation 6")
sim6 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5],  measurement_error = "additive", n_groups = nGroup, rng = MersenneTwister(rng_int));
CSV.write(
    "data/simulations/est_wtp_sim6.csv",
    vcat(
        DataFrame(sim = "sim6", est = "base", wtp = sim6[1][2,:] ./ sim6[1][1, :]),
        DataFrame(sim = "sim6", est = "groupfe", wtp = sim6[2][2,:] ./ sim6[2][1, :]),
        DataFrame(sim = "sim6", est = "cf", wtp = sim6[3][2,:] ./ sim6[3][1, :])
    ),
    writeheader = true
);
rng_int += 1;

#   STEP 2: Generate sample datasets for simulations 1-6 -------------

#   Generate sample datasets for simulations 1-3:
println("Simulating example data generating processes.")
sim1_data = gen_data(nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, 0.0, 0.5], n_groups = nGroup);
sim2_data = gen_data(nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups =nGroup);
sim3_data = gen_data(nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, 1.0, 0.5], n_groups = nGroup);

#   Generate sample datasets for simulations 4-6:
sim4_data = gen_data(nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 0.0], [5.0, 1.0, 0.0, 0.5], n_groups = nGroup);
sim5_data = gen_data(nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 0.0], [5.0, 1.0, 0.0, 0.5], n_groups = nGroup);
sim6_data = gen_data(nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 0.0], [5.0, 1.0, -1.0, 0.5], n_groups =nGroup);
sim4_data_cost_obs = sim4_data[2] .+ rand(MersenneTwister(rng_int), Normal(0, 0.5), size(sim4_data[2]));
sim5_data_cost_obs = sim5_data[2] .+ sim5_data[2] .* rand(MersenneTwister(rng_int), Normal(0, 0.15), size(sim5_data[2]));
sim6_data_cost_obs = sim6_data[2] .+ rand(MersenneTwister(rng_int), Normal(0, 0.5), size(sim6_data[2]));

#   Reshape and save results:
CSV.write(
    "data/simulations/sim_data_sim1to6.csv", 
    vcat(
        DataFrame(sim = "sim1", xi = reshape(sim1_data[3], :), cost = reshape(sim1_data[2], :), cost_obs = reshape(sim1_data[2], :)),
        DataFrame(sim = "sim2", xi = reshape(sim2_data[3], :), cost = reshape(sim2_data[2], :), cost_obs = reshape(sim2_data[2], :)),
        DataFrame(sim = "sim3", xi = reshape(sim3_data[3], :), cost = reshape(sim3_data[2], :), cost_obs = reshape(sim3_data[2], :)),
        DataFrame(sim = "sim4", xi = reshape(sim4_data[3], :), cost = reshape(sim4_data[2], :), cost_obs = reshape(sim4_data_cost_obs, :)),
        DataFrame(sim = "sim5", xi = reshape(sim5_data[3], :), cost = reshape(sim5_data[2], :), cost_obs = reshape(sim5_data_cost_obs, :)),
        DataFrame(sim = "sim6", xi = reshape(sim6_data[3], :), cost = reshape(sim6_data[2], :), cost_obs = reshape(sim6_data_cost_obs, :))
    ), 
    writeheader=true
);

#   STEP 3: Robustness of estimators to different degrees of w/in group correlation of unobservables -------------

# Loop over different values governing correlation of unobservable ξ:
println("Starting Robustness Check 1")
ξ_cor = collect(0.0:0.05:1.0);
sim2_rc1_bias = zeros(size(ξ_cor, 1), 3);
sim2_rc1_tstat = zeros(size(ξ_cor, 1), 3);
for i = eachindex(ξ_cor)

    #   Run simulation 2 w/ different values for parameter governing w/in group correlation:
    sim2_temp = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, ξ_cor[i]], n_groups = nGroup, rng = MersenneTwister(rng_int));

    #   Calculate bias/tstat from each estimator in simulation 2:
    for j = 1:3
        sim2_rc1_bias[i, j] = mean((sim2_temp[j][2, :] ./ sim2_temp[j][1, :]) .+ 0.5)
        sim2_rc1_tstat[i, j] = mean((sim2_temp[j][2, :] ./ sim2_temp[j][1, :]) .+ 0.5) / (std((sim2_temp[j][2, :] ./ sim2_temp[j][1, :])) / sqrt(nSim))
    end

end;

#   Save results:
CSV.write(
    "data/simulations/est_wtp_rc1.csv",
    DataFrame(corel = repeat(ξ_cor, 3), est = repeat(["base", "groupfe", "cf"], inner = size(ξ_cor, 1)), bias = reshape(sim2_rc1_bias, :), tstat = reshape(sim2_rc1_tstat, :)),
    writeheader = true
);

#   STEP 4: Robustness of estimators to different group sizes -------------

# Define a range of values for n_groups
n_groups_values = [5, 10, 20, 25, 50, 100];

# Initialize arrays to store results
println("Starting Robustness Check 2")
sim2_rc2_bias = zeros(size(n_groups_values, 1), 3);
sim2_rc2_tstat = zeros(size(n_groups_values, 1), 3);
iter = 1
for n_groups in n_groups_values

    # Run simulation 2 w/ ξ_cor set to 0.5
    sim2_temp = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups = n_groups, rng = MersenneTwister(rng_int))

    # Calculate bias/tstat from each estimator in simulation 2
    for j in 1:3
        sim2_rc2_bias[iter, j] = mean((sim2_temp[j][2, :] ./ sim2_temp[j][1, :]) .+ 0.5)
        sim2_rc2_tstat[iter, j] = mean((sim2_temp[j][2, :] ./ sim2_temp[j][1, :]) .+ 0.5) / (std((sim2_temp[j][2, :] ./ sim2_temp[j][1, :])) / sqrt(nSim))
    end
    global iter += 1

end;

#   Save results:
CSV.write(
    "data/simulations/est_wtp_rc2.csv",
    DataFrame(n_groups = repeat(n_groups_values, 3), est = repeat(["base", "groupfe", "cf"], inner = size(n_groups_values, 1)), bias = reshape(sim2_rc2_bias, :), tstat = reshape(sim2_rc2_tstat, :)),
    writeheader = true
);

#   STEP 5: Robustness of estimators to nonlinearities in true DGP -------------

#   Simulation 2 w/ quadratic unobservables:
println("Starting Robustness Check 3")
sim2_rc_nl1 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups = nGroup, nonlinear_c = "quadratic", rng = MersenneTwister(rng_int));
rng_int += 1;
sim2_rc_nl2 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups = nGroup, nonlinear_u = "quadratic", rng = MersenneTwister(rng_int));
rng_int += 1;
sim2_rc_nl3 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups = nGroup, nonlinear_c = "quadratic", nonlinear_u = "quadratic", rng = MersenneTwister(rng_int));
rng_int += 1;

#   Simulation 2 w/ exponential unobservables:
sim2_rc_nl4 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups = nGroup, nonlinear_c = "exponential", rng = MersenneTwister(rng_int));
rng_int += 1;
sim2_rc_nl5 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups = nGroup, nonlinear_u = "exponential", rng = MersenneTwister(rng_int));
rng_int += 1;
sim2_rc_nl6 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups = nGroup, nonlinear_c = "exponential", nonlinear_u = "exponential", rng = MersenneTwister(rng_int));
rng_int += 1;

#   Combine and save results:
nl_res = [sim2_rc_nl1, sim2_rc_nl2, sim2_rc_nl3, sim2_rc_nl4, sim2_rc_nl5, sim2_rc_nl6];
nl_res_export = DataFrame()
for j = 1:6
    for i = 1:3
        append!(nl_res_export, DataFrame(sim = "sim2", est = i, nonlinear = j, wtp = nl_res[j][i][2,:] ./ nl_res[j][i][1, :]))
    end
end
CSV.write(
    "data/simulations/est_wtp_rc3.csv",
    nl_res_export,
    writeheader = true
);

#   Simulation 2 w/ quadratic Z in cost:
println("Starting Robustness Check 4")
sim2_rc_nl7 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups = nGroup, nonlinear_c_z = "quadratic", rng = MersenneTwister(rng_int));
rng_int += 1;

#   Simulation 2 w/ exponential Z in cost:
sim2_rc_nl8 = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, 1.0, -1.0, 0.5], n_groups = nGroup, nonlinear_c_z = "exponential", rng = MersenneTwister(rng_int));
rng_int += 1;

#   Combine and save results:
nl_res = [sim2_rc_nl1, sim2_rc_nl2, sim2_rc_nl3, sim2_rc_nl4, sim2_rc_nl5, sim2_rc_nl6, sim2_rc_nl7, sim2_rc_nl8];
nl_res_export = DataFrame()
for j = 7:8
    for i = 1:3
        append!(nl_res_export, DataFrame(sim = "sim2", est = i, nonlinear = j, wtp = nl_res[j][i][2,:] ./ nl_res[j][i][1, :]))
    end
end
CSV.write(
    "data/simulations/est_wtp_rc4.csv",
    nl_res_export,
    writeheader = true
);

#   STEP 6: Robustness of estimators to weak instruments -------------

# Define a range of values for relationship between Z and cost:
z_correl = [1.0 0.5 0.25 0.1 0.05 0.01];

# Initialize arrays to store results
println("Starting Robustness Check 5")
sim2_rc_weak_bias = zeros(size(z_correl, 2), 3);
sim2_rc_weak_tstat = zeros(size(z_correl, 2), 3);
sim2_rc_weak_fstats = zeros(nSim, size(z_correl, 2);)
iter = 1
for z in z_correl

    # Run simulation 2 w/ different relationships between Z and cost:
    sim2_temp = mc_simulations(nSim, nInd, nAlt, f_ξ, [1.0, 1.0, -2.0, 1.0], [5.0, z, -1.0, 0.5], n_groups = nGroup, f_stat = true, rng = MersenneTwister(rng_int))

    # Calculate bias/tstat from each estimator in simulation 2
    for j in 1:3
        sim2_rc_weak_bias[iter, j] = mean((sim2_temp[j][2, :] ./ sim2_temp[j][1, :]) .+ 0.5)
        sim2_rc_weak_tstat[iter, j] = mean((sim2_temp[j][2, :] ./ sim2_temp[j][1, :]) .+ 0.5) / (std((sim2_temp[j][2, :] ./ sim2_temp[j][1, :])) / sqrt(nSim))
    end

    #   Add F-statistics to array:
    sim2_rc_weak_fstats[:, iter] = sim2_temp[4]
    global iter += 1

end;

#   Save results:
CSV.write(
    "data/simulations/est_wtp_rc_weak.csv",
    DataFrame(z_correl = repeat(z_correl[1,:], 3), est = repeat(["base", "groupfe", "cf"], inner = size(z_correl, 2)), bias = reshape(sim2_rc_weak_bias, :), tstat = reshape(sim2_rc_weak_tstat, :)),
    writeheader = true
);
#   Save results:
CSV.write(
    "data/simulations/est_wtp_rc_weak_fstat.csv",
    DataFrame(z_correl = repeat(z_correl[1, :], inner = nSim), f_stat = reshape(sim2_rc_weak_fstats, :)),
    writeheader = true
);
