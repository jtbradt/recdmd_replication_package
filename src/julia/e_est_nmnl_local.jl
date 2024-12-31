################################################################################
#   SCRIPT:     4_est_nmnl_local.jl
#   AUTHOR:     Jacob Bradt (jacob.bradt@mccombs.utexas.edu)
#   NOTES:      julia script to run robustness estimates reported in Table A4 of 
#               "Hotelling Meets Wright: Spatial Sorting and Measurement Error 
#               in Recreation Demand Models." *Please note* that I set the 
#               number of bootstrap iterations to 20 for the sake of time in 
#               replication. For the paper, I use 200 bootstrap iterations 
#               As a result, the results in the paper may differ
#               slightly from the results in this replication.
################################################################################

#   Import dependencies:
using DataFrames, LinearAlgebra, DelimitedFiles, Plots, TexTables, CSV

#   Import user-written functions:
include("e_nmnl.jl")

#   Set number of bootstrap iterations:
n_boot = 20;        #**NOTE**: Set to 200 in paper; lower number used for speed of replication

#   STEP 1: Baseline nested logit estimates -------------

#   Import data:
tc_data = readdlm("data/dwh_replication/tc1.txt");
demo_data = readdlm("data/dwh_replication/trips_demos.txt");

#   Subset data to local trips:
local_caseids = readdlm("data/dwh_replication/respondents_local.csv", ',', header = false)[2:end, 1];
rows_to_keep = in(local_caseids).(tc_data[:, 1]);
tc_data = tc_data[rows_to_keep, :];
demo_data = demo_data[rows_to_keep, :];

#   Create choice matrix:
choi = demo_data[:, 15:97];
choi = hcat(demo_data[:, 3] .- sum(choi, dims = 2), choi);

#   Create weights vector:
wgts = demo_data[:, 3] ./ sum(demo_data[:, 3])

#   Create cost matrix, normalized by $100:
cost = tc_data[:, 4:86] ./ 100;

#   Baseline nested logit:
θ_hat1, θ_hat1_boot = est_nmnl_full(choi, cost, [-1.0, 0.8], maxiter = 500, n_boot = n_boot);

#   Save results to file:
CSV.write(
    "data/dwh_replication/est/est_nmnl_local1.csv", 
    DataFrame(θ_hat1, :auto), 
    writeheader=false
);
CSV.write(
    "data/dwh_replication/boot_est/boot_est_nmnl_local1.csv", 
    DataFrame(θ_hat1_boot, :auto), 
    writeheader=false
);

#   Import results from file, if already estimated:
#θ_hat1 = readdlm("data/dwh_replication/est/est_nmnl_local1.csv", ',', header = false)[:, 1];
#θ_hat1_boot = readdlm("data/dwh_replication/boot_est/boot_est_nmnl_local1.csv", ',', header = false);

#   Calculate standard errors:
se_hat1 = sqrt.(diag(cov(θ_hat1_boot, dims = 2)));

#   STEP 2: Baseline spill constant calibration -------------

#   Estimated (percent) reduction in visitation for Northern Gulf and 
#   Florida Panhandle regions (source: Tourganeau et al., 2017):
redu_ng = 45.224 / 100;       
redu_fp = 22.151 / 100;         

#   Construct matrix coding NG and FP sites:
spill_sites = zeros(83, 2);
spill_sites[10:35, 1] .= 1.0;        # NG sites
spill_sites[36:63, 2] .= 1.0;        # FP sites

#   Calibrate group-level adjustment constants for spill conditions
spill_adj1 = berry_contraction([θ_hat1[1]], θ_hat1[2], θ_hat1[3:end], choi, cost, spill_sites, [redu_ng, redu_fp], tol = 1e-8);

#   STEP 3: Nested logit w/ CF -------------

#   Load instruments:
inst_files = readdir("data/dwh_replication/zmat");
Z = [readdlm(string("data/dwh_replication/zmat/", inst_files[i])) for i = 1:size(inst_files, 1)];
Z_full = vcat(Z[[2, 3, 5, 6, 8, 9, 10]], [Z[8] .* Z[11], Z[9] .* Z[11], Z[10] .* Z[11]]) |> Vector{Matrix{Float64}}
Z_full = [Z_full[i][rows_to_keep, :] for i = 1:size(Z_full, 1)]

#   Nested logit w/ CF:
θ_hat2, θ_hat2_boot = est_nmnl_full(choi, cost, [-1.0, 0.0, 0.8], Z =Z_full, maxiter = 500, n_boot = n_boot);

#   Save results to file:
CSV.write(
    "data/dwh_replication/est/est_nmnl_local2.csv", 
    DataFrame(θ_hat2, :auto), 
    writeheader=false
);
CSV.write(
    "data/dwh_replication/boot_est/boot_est_nmnl_local2.csv", 
    DataFrame(θ_hat2_boot, :auto), 
    writeheader=false
);

#   Import results from file, if already estimated:
#θ_hat2 = readdlm("data/dwh_replication/est/est_nmnl_local2.csv", ',', header = false)[:, 1];
#θ_hat2_boot = readdlm("data/dwh_replication/boot_est/boot_est_nmnl_local2.csv", ',', header = false);


#   Calculate standard errors:
se_hat2 = sqrt.(diag(cov(θ_hat2_boot, dims = 2)));

#   STEP 4: Spill constant calibration for nested logit w/ CF -------------

#   Calibrate group-level adjustment constants for spill conditions
μ = first_stage(cost, Z_full, wgts[:, 1])
X_full = reduce((x,y) -> cat(x, y, dims = 3), [[cost] [μ]])
spill_adj2 = berry_contraction(θ_hat2[1:2], θ_hat2[3], θ_hat2[4:end], choi, X_full, spill_sites, [redu_ng, redu_fp], tol = 1e-8);

#   STEP 5: Calculate lost user-day value with each estimator -------------

#   Lost user-day value with baseline estimator:
wtp1 = wtp([θ_hat1[1]], θ_hat1[2], θ_hat1[3:end] .+ (spill_sites[:, 1] .* spill_adj1[1]) .+ (spill_sites[:, 2] .* spill_adj1[2]), θ_hat1[3:end], choi, cost);
ud_loss1 = ((wtp1[1] * 100) / wtp1[2]) ./ 1.7; 

#   Lost user-day value with CF estimator:
wtp2 = wtp(θ_hat2[1:2], θ_hat2[3], θ_hat2[4:end] .+ (spill_sites[:, 1] .* spill_adj2[1]) .+ (spill_sites[:, 2] .* spill_adj2[2]), θ_hat2[4:end], choi, X_full);
ud_loss2 = ((wtp2[1] * 100) / wtp2[2]) ./ 1.7;

#   STEP 6: Construct estimates table -------------

#   Create table of estimates:
@fmt Real = "{:.3f}" 
@fmt Int  = "{:,n}"
key1 = ["\$\\alpha\$", "\$\\eta\$"];
t1 = TableCol("(1)");
for (k, v, p) in zip(key1, θ_hat1[1:2], se_hat1[1:2])
    t1[k] = v, p
end;
res2_print = [θ_hat2[i] for i = [1, 3, 2]];
res2_se_print = [se_hat2[i] for i = [1, 3, 2]];
key2 = ["\$\\alpha\$", "\$\\eta\$", "\$\\lambda\$"];
t2 = TableCol("(2)");
for (k, v, p) in zip(key2, res2_print, res2_se_print)
    t2[k] = v, p
end;
t3 = TableCol("(1)");
N1 = ["Yes", FormattedNumber(ud_loss1), FormattedNumber(size(choi, 1)), FormattedNumber(size(choi, 2) - 1)];
key3 = ["Alternative Specific Constants", "Lost User Day Value (\$/day)", "\$N\$", "Sites"];
for (k, v) in zip(key3, N1)
    t3[k] = v
end;
t4 = TableCol("(2)");
N2 = ["Yes", FormattedNumber(ud_loss2), FormattedNumber(size(choi, 1)), FormattedNumber(size(choi, 2) - 1)];
for (k, v) in zip(key3, N2)
    t4[k] = v
end;
t_full = [t1 t2 
            t3 t4]
t_tex = to_tex(t_full, se_pos=:inline );
open("output/dwh_replication/dwh_second_stage_local.tex", "w") do file
    write(file, t_tex)
end;
