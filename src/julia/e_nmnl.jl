
################################################################################
#   SCRIPT:     mnl.jl
#   AUTHOR:     Jacob Bradt (jacob.bradt@mccombs.utexas.edu)
#   NOTES:      julia script containing all functions necessary to estimate 
#               recreation demand models in English et al. (2018) replication in
#               "Hotelling Meets Wright: Spatial Sorting and Measurement Error 
#               in Recreation Demand Models."
################################################################################

#   Import dependencies:
using LinearAlgebra, Optim, ForwardDiff, LineSearches, NLSolversBase, Random, Distributions, SparseArrays

#   Function: make_dummy() ----------------------------------------
#   Description: fxn to create dummy indicator matrix
#   Arguments:
#       1. data = vector of values

function make_dummy(data::AbstractVector)
    #   Get unique values:
    udata = sort(unique(data))

    #   Get indices of each value:
    inds = map(x -> findfirst(==(x), udata), data)

    #   Create dummy matrix:
    sparse(1:length(data), inds, ones(Int, length(data)))
end;


#   Function: my_ln() ----------------------------------------
#   Description: fxn to censor small values for natural log
#   Arguments:
#       1. x = vector of values

function my_ln(x)
    y = (x .< 1e-250).*1e-250 .+ (x .>= 1e-250).*x
    return log.(y)
 end;
 

#   Function: my_exp() ----------------------------------------
#   Description: fxn to censor small and large values for exponentiation
#   Arguments:
#       1. x = vector of values

 function my_exp(x)
    y = (x .< -700).*(-700) .+ (abs.(x).<=700) .* x + (x .>=700) .*(700)
    return exp.(y)
 end;
 

#   Function: nmnl_prob() ----------------------------------------
#   Description: fxn to generate nested multinomial logit probability
#   Arguments:
#       1. ch = array containing choices (wide format)
#       2. X = array containing observable data entering utility (wide format)
#       3. θ = vector containing full array of params
#       4. λ = vector containing nest parameter
#       5. δ = vector of alternative-specific constants

function nmnl_prob(
    ch::Array{Float64},
    X::Array{Float64},
    θ::Array{Float64},
    λ::Float64,
    δ::Array{Float64}
)

    #   Initialize indirect utility:
    nv = zeros(size(ch))

    #   Fill in indirect utility for inside options:
    for i = 1:size(X, 3)
        nv[:, 2:end] += X[:, :, i] .* θ[i]
    end

    #   Add ASCs for inside options:
    nv[:, 2:end] = nv[:, 2:end] .+ δ'

    #   Divide by similarity coef:
    nv = nv ./ λ

    #   Calculate unconditional choice probability:
    env = exp.(nv)
    sum_env = sum(env[:, 2:end], dims = 2)
    p_0 = 1 ./ (1 .+ sum_env .^λ)
    p_inside = (env[:, 2:end] ./ sum_env) .* ((sum_env).^ λ) ./ (1 .+ (sum_env).^ λ)
    return hcat(p_0, p_inside)

end;


#   Function: agg_shares() ----------------------------------------
#   Description: fxn to generate aggregate shares
#   Arguments:
#       1. ch = array containing chosen alternative data (wide format) 

function agg_shares(ch::Array{Float64})

    # Calculate shares based on choice matrix:
    shares = sum(ch[:, 2:end], dims = 1) ./ sum(ch)
    return shares'
    
end;


#   Function: ll_nmnl() ----------------------------------------
#   Description: fxn to generate log-likelihood of nested mnl model
#   (type agnostic to enable numerical differentiation)
#   Arguments:
#       1. ch = array containing chosen alternative data (wide format) 
#       2. X = array containing observable data entering utility (wide format)
#       3. θ = vector containing full array of params
#       4. λ = vector containing nest parameter
#       5. δ = vector of alternative-specific constants

function ll_nmnl(
    ch::Array{Float64}, 
    X::Array{Float64}, 
    θ::Array{Float64}, 
    λ::Float64, 
    δ::Array{Float64}
)

    #   Initialize indirect utility:
    nv = zeros(size(ch))

    #   Fill in indirect utility for inside options:
    for i = 1:size(X, 3)
        nv[:, 2:end] += X[:, :, i] .* θ[i]
    end

    #   Add ASCs for inside options:
    nv[:, 2:end] = nv[:, 2:end] .+ δ'

    #   Divide by similarity coef:
    nv = nv ./ λ

    #   Calculate nmnl log likelihood:
    env = exp.(nv)
    arg1 = sum(env[:, 2:end], dims = 2)
    ll = sum(ch .* nv, dims = 2) .+ sum(ch[:, 2:end], dims = 2) .* (λ - 1) .* log.(arg1) .- sum(ch, dims = 2) .* log.(1 .+ arg1 .^λ)

    return mean(ll)

end;


#   Function: ll_nmnl_grad() ----------------------------------------
#   Description: fxn to calculate gradient of log-likelihood of nested mnl model
#   Arguments:
#       1. ch = array containing chosen alternative data (wide format) 
#       2. X = array containing observable data entering utility (wide format)
#       3. θ = vector containing full array of params
#       4. λ = vector containing nest parameter
#       5. δ = vector of alternative-specific constants

function ll_nmnl_grad(
    ch::Array{Float64}, 
    X::Array{Float64}, 
    θ::Array{Float64}, 
    λ::Float64, 
    δ::Array{Float64}
)

    #   Unpack choice data:
    ch_outside = ch[:, 1]
    ch_inside = ch[:, 2:end]

    #   Initialize indirect utility:
    nv = zeros(size(ch_inside))

    #   Fill in indirect utility for inside options:
    θ_wide = reshape(θ, (1,1,length(θ)))
    nv = sum(X .* θ_wide, dims =3 )[:, :, 1]

    #   Add ASCs for inside options:
    nv = nv .+ δ'

    #   Divide by similarity coef:
    nv = nv ./ λ

    #   Calculate choice probabilities:
    env = my_exp(nv)
    arg1 = sum(env, dims = 2)
    denom1 = 1 .+ arg1.^λ
    numer1 = hcat(ones(size(nv, 1)), env .* (arg1 .^ (λ - 1)))
    prob1 = numer1 ./ denom1
    cprob1 = env ./ arg1

    # Calculate gradient of ASCs: 
    ntrips = sum(ch_inside, dims = 2)
    co = ch_outside + ntrips
    g_asc = (ch_inside ./ λ) .+ (ntrips .* (λ - 1) .* cprob1 ./ λ ) .- (co .* prob1[:, 2:end])

    # Calculate gradient of observable terms:
    g_x = sum(g_asc .* X, dims = [1,2])[1, 1, :]

    # Calculate gradient of nest param:
    ln_arg1 = my_ln(arg1)
    g_λ_arg1 = -sum((ch_inside .* nv) ./ λ) 
    g_λ_arg2 = sum(ntrips .* ln_arg1)
    g_λ_arg3 = -sum((cprob1 .* nv) .* (ntrips * (λ-1)/λ))
    g_λ_arg4 = -sum(co .* (1 .- prob1[:, 1]) .* ln_arg1)
    g_λ_arg5 = sum(co .* prob1[:, 2:end] .* nv)
    g_λ = g_λ_arg1 + g_λ_arg2 + g_λ_arg3 + g_λ_arg4 + g_λ_arg5
    
    #   Return gradient
    return vcat(g_x, g_λ, sum(g_asc, dims = 1)') ./ size(ch, 1)
    
end;

#   Function: first_stage() ----------------------------------------
#   Description: fxn to implement first stage by OLS
#   Arguments:
#       1. tc = array containing travel cost data (wide format)
#       2. Z = vector of arrays containing instruments (wide format)
#       3. wgts = vector of weights for OLS
#       4. Z_fe = vector of arrays containing choice occassion-specific fixed effects

function first_stage(
    tc::Array{Float64}, 
    Z::Vector{Matrix{Float64}}, 
    wgts::Array{Float64}; 
    Z_fe = nothing
)

    #   Get # of variables and alternatives:
    n_var = size(Z, 1)
    n_alt = size(Z[1], 2)
    n_obs = size(Z[1], 1)

    #   Reshape weights:
    wgts = repeat(wgts[:, 1], inner = n_alt)
    W = Diagonal(wgts)

    #   Construct ASC matrix:
    ascs = zeros(n_obs * n_alt, n_alt)
    for i = 1:n_obs
        ascs[((i-1)*n_alt + 1):(i*n_alt),:] .= I(n_alt)
    end

    #   Construct indicator matrices for any FE:
    if !isnothing(Z_fe)
        Z_fe = [make_dummy(repeat(Z_fe[i], inner = n_alt))[:, 2:end] for i = 1:size(Z_fe, 1)]
    end

    #   Construct matrices for OLS estimation:
    X = zeros(n_alt*n_obs, n_var)
    for i = 1:n_var
        X[:, i] .= reshape(Z[i]', n_obs*n_alt)
    end
    if !isnothing(Z_fe)
        X = hcat(X, ascs, reduce(hcat, Z_fe))
    else
        X = hcat(X, ascs)
    end
    X = Matrix(X)

    #   Reshape TC data:
    Y = reshape(tc', n_obs*n_alt)
    Y = Y

    #   Estimate OLS params:
    θ = inv(X'*W*X)*(X'*W*Y)

    #   Fit residuals
    resid = Y .- X * θ
    resid = permutedims(reshape(resid, n_alt, n_obs))

    return resid

end;



#   Function: berry_contraction_est() ----------------------------------------
#   Description: fxn to implement Berry contraction for calibration of initial 
#       alternative-specific constants
#   Arguments:
#       1. ch = array containing chosen alternative data (wide format) 
#       2. X = array containing observable data entering utility (wide format)
#       3. θ = vector containing full array of params
#       4. λ = vector containing nest parameter
#       5. δ = vector of alternative-specific constants
#       6. tol = tolerance for convergence
#       7. max_iter = maximum number of iterations

function berry_contraction_est(
    ch::Array{Float64}, 
    X::Array{Float64}, 
    θ::Array{Float64}, 
    λ::Float64, 
    δ_0::Array{Float64}; 
    tol = 1e-2,
    max_iter = 300
)

    #   Calculate observed shares:
    s_act = agg_shares(ch)

    #   Initialize values
    δ_1 = δ_0
    eps = 1.0
    iter = 1
    δ_store = zeros(size(δ_0, 1), 1)
 
    #   While loop:
    while (eps > tol) * (max_iter > iter)
 
        #   Calculate model-implied choice probabilities with grouped site correction:
        pr_model = nmnl_prob(ch, X, θ, λ,δ_0)

        #   Calculate model-implied aggregate shares:
        co = sum(ch, dims = 2)
        pr_model = pr_model .* co
        pr_model = sum(pr_model, dims = 1) ./ sum(co)

        #   Update δ:
        δ_1 = δ_0 .+ (1 - λ) * (log.(s_act) .- log.(pr_model[:, 2:end])')
        
        #   Calculate convergence criterion:
        eps = maximum(abs.(δ_1 .- δ_0))
        iter += 1

        #   Update next iteration
        if (eps > tol)
            δ_store = hcat(δ_store, δ_1)
            if iter > 50
                δ_0 = δ_1
            else
                δ_0 = 0.5 * δ_1 + 0.5 * δ_0
            end
        else
            δ_0 = δ_1
        end
    end

    return δ_1

end

#   Function: est_nmnl() ----------------------------------------
#   Description: fxn to estimate nested multinomial logit
#   Arguments:
#       1. ch = array containing chosen alternative data (wide format) 
#       2. X = array containing observable data entering utility (wide format)
#       3. θ_start = vector containing starting values for full array of params
#       4. Z = array containing instruments (wide format)
#       5. Z_fe = array containing choice occassion-specific fixed effects
#       6. gradient = boolean indicating whether to use gradient-based optimization
#       7. maxiter = maximum number of iterations
#       8. verbose = boolean indicating whether to print optimization progress
#       9. g_tol = tolerance for convergence

function est_nmnl(
    ch::Array{Float64}, 
    X::Array{Float64}, 
    θ_start::Array{Float64}; 
    Z = nothing, 
    Z_fe = nothing, 
    gradient = false, 
    maxiter = 1_000, 
    verbose = false, 
    g_tol = 1e-8
)

    #   Get parameter numbers:
    n_param = size(X, 3)
    n_asc = size(ch, 2) - 1
    δ_start = -6.0 * ones(n_asc)

    #   If instruments specified, estimate first stage:
    if !isnothing(Z)
        wgts = sum(choi, dims = 2) ./ sum(choi)
        μ = first_stage(X[:, :, 1], Z, wgts[:, 1], Z_fe = Z_fe)
        X = reduce((x,y) -> cat(x, y, dims = 3), [[X] [μ]]);
        n_param = n_param + 1
    end

    #   Update δ starting values:
    δ_start = berry_contraction_est(ch, X, θ_start[1:(end-1)], θ_start[end], δ_start, max_iter = 100)

    #   Define objective:
    function f(x)
       
        #   Objective:
        ll = ll_nmnl(ch, X, x[1:n_param], x[n_param + 1], x[(n_param + 2):end])
        return -ll

    end

    #   Define function for gradient-based optimization:
    function mnl_fg!(F,G, x)

        #   If gradient called in iteration, calculate;
        if G !== nothing
          G[1:end] = -ll_nmnl_grad(ch, X, x[1:n_param], x[n_param + 1], x[(n_param + 2):end])
        end

        #   If log-likelihood objective function called, calculate:
        if F !== nothing
          ll = ll_nmnl(ch, X, x[1:n_param], x[n_param + 1], x[(n_param + 2):end])
          return -ll
        end

    end

    #   Estimate params by MLE:
    if gradient
        res = Optim.optimize(Optim.only_fg!(mnl_fg!), vcat(θ_start, δ_start), LBFGS(linesearch = LineSearches.BackTracking()), Optim.Options(show_trace = verbose, iterations = maxiter, g_tol = g_tol))
    else
        res = Optim.optimize(f, vcat(θ_start, δ_start), Optim.Options(show_trace = verbose, iterations = maxiter, g_tol = g_tol))
    end

    #   Get final parameter estimates:
    θ_hat = Optim.minimizer(res)

    return θ_hat

end;


#   Function: vcov_boot() ----------------------------------------
#   Description: fxn to estimate nested multinomial logit
#   Arguments:
#       1. ch = array containing chosen alternative data (wide format) 
#       2. X = array containing observable data entering utility (wide format)
#       3. θ_start = vector containing starting values for full array of params
#       4. n_boot = number of bootstrap iterations
#       5. Z = array containing instruments (wide format)
#       6. Z_fe = array containing choice occassion-specific fixed effects
#       7. gradient = boolean indicating whether to use gradient-based optimization
#       8. maxiter = maximum number of iterations
#       9. g_tol = tolerance for convergence
#       10. rng = random number generator

function vcov_boot(
    ch::Array{Float64}, 
    X::Array{Float64}, 
    θ_start::Array{Float64}; 
    n_boot = 100, 
    Z = nothing, 
    Z_fe = nothing, 
    gradient = false, 
    maxiter = 1_000, 
    g_tol = 1e-8, 
    rng = MersenneTwister(82220)
)

    #   Get number of threads:
    n_threads = Threads.nthreads()

    #   Initialize matrix to store estimates:
    n_param = size(X, 3)
    n_asc = size(ch, 2) - 1
    n_obs = size(ch, 1)
    if !isnothing(Z)
        θ_boot = [zeros(n_param + n_asc + 2, 0) for i=1:n_threads]
    else
        θ_boot = [zeros(n_param + n_asc + 1, 0) for i=1:n_threads]
    end
    boot_counter = 0

    #   Generate RNGs:
    rngs = [MersenneTwister(i) for i in rand(rng, 1:10^5, n_boot)]

    #   Loop over bootstrap iterates:
    Threads.@threads for j = 1:n_boot

        #   Draw sample for  bootstrap:
        boot_sample = rand(rngs[j], 1:n_obs, n_obs)

        #   Run estimation with weights for bootstrap iterate:
        if !isnothing(Z)
            
            #   If indicators included in first stage, specify bootstrap sample of FE vector:
            if !isnothing(Z_fe)
                θ_boot_temp = est_nmnl(ch[boot_sample, :], X[boot_sample, :, :], θ_start, 
                    Z = [Z[i][boot_sample, :] for i = 1:size(Z, 1)], 
                    Z_fe = [Z_fe[i][boot_sample, :] for i = 1:size(Z_fe, 1)], 
                    verbose = false, gradient = gradient, maxiter = maxiter, g_tol = g_tol)
            else
                θ_boot_temp = est_nmnl(ch[boot_sample, :], X[boot_sample, :, :], θ_start, 
                    Z = [Z[i][boot_sample, :] for i = 1:size(Z, 1)], 
                    verbose = false, gradient = gradient, maxiter = maxiter, g_tol = g_tol)
            end
            
        else

            θ_boot_temp = est_nmnl(ch[boot_sample, :], X[boot_sample, :, :], θ_start, verbose = false, gradient = gradient, maxiter = maxiter, g_tol = g_tol)

        end

        #   Append bootstrap estimate to storage matrix:
        θ_boot[Threads.threadid()] = hcat(θ_boot[Threads.threadid()], θ_boot_temp)

        #   Print status:
        boot_counter += 1
        println("------------------------------------------------------------")
            println("******** BOOTSTRAP ", boot_counter, "/", n_boot," COMPLETE *************")
            if !isnothing(Z)
                println("       Sample est.: ", round.(θ_boot_temp[1:(n_param + 2)], digits = 3))
            else
                println("       Sample est.: ", round.(θ_boot_temp[1:(n_param + 1)], digits = 3))
            end
    end

    return(θ_boot)

end;

#   Function: est_nmnl_full() ----------------------------------------
#   Description: fxn to implement full estimation of nested multinomial logit
#   Arguments:
#       1. ch = array containing chosen alternative data (wide format) 
#       2. X = array containing observable data entering utility (wide format)
#       3. θ_start = vector containing starting values for full array of params
#       4. Z = array containing instruments (wide format)
#       5. Z_fe = array containing choice occassion-specific fixed effects
#       6. maxiter = maximum number of iterations
#       7. g_tol = tolerance for convergence
#       8. n_boot = number of bootstrap iterations

function est_nmnl_full(
    ch::Array{Float64}, 
    X::Array{Float64}, 
    θ_start::Array{Float64}; 
    Z = nothing, 
    Z_fe = nothing, 
    maxiter = 1_000, 
    g_tol = 1e-8, 
    n_boot = 100
)

    #   Main estimates:
    θ_hat = est_nmnl(ch, X, θ_start, Z = Z, Z_fe = Z_fe, maxiter = maxiter, g_tol = g_tol, verbose = false, gradient = true)

    #   Print status:
    println("------------------------------------------------------------")
    println("******** MAIN ESTIMATES COMPLETE *************")
    if !isnothing(Z)
        println("       Main est.: ", round.(θ_hat[1:(size(X, 3) + 2)], digits = 3))
    else
        println("       Main est.: ", round.(θ_hat[1:(size(X, 3) + 1)], digits = 3))
    end

    #   Bootstrap estimates:
    θ_boot = vcov_boot(ch, X, θ_start, Z = Z, Z_fe = Z_fe, n_boot = n_boot, maxiter = maxiter, g_tol = g_tol, gradient = true)

    #   Return main estimates, bootstrap estimates:
    return (θ_hat, reduce(hcat, θ_boot))

end;


#   Function: berry_contraction() ----------------------------------------
#   Description: fxn to implement Berry contraction mapping for calibration of 
#       spill condition constants
#   Arguments:
#       1. θ = vector containing full array of params
#       2. λ = vector containing nest parameter
#       3. δ_site = vector of alternative-specific constants
#       4. ch = array containing choices (wide format)
#       5. X = array containing observable data entering utility (wide format)
#       6. spill_sites = array containing indicators of group assignment (NG or FP)
#       7. spill_delta = array containing proportional reduction in visitation from 
#           spill for each group (NG and FP)

function berry_contraction(
    θ::Array{Float64}, 
    λ::Float64, 
    δ_site::Array{Float64}, 
    ch::Array{Float64}, 
    X::Array{Float64}, 
    spill_sites::Array{Float64}, 
    spill_delta::Array{Float64}; 
    tol = 1e-2, 
    max_iter = 300
)
    
    #   Calculate observed aggregate shares:
    s_g_act = vcat(sum(agg_shares(ch) .* spill_sites[:, 1]), sum(agg_shares(ch) .* spill_sites[:, 2]))
    s_g_act = s_g_act .* (1 .- spill_delta)

    #   Initialize values
    δ_g_1 = zeros(2)
    s_g_0 = zeros(size(s_g_act))
    δ_full = zeros(size(δ_site))
    eps = 1.0
    iter = 1

    #   While loop:
    while (eps > tol) * (max_iter > iter)

        #   Update w/ last iteration δ:
        δ_g_0 = δ_g_1

        #   Construct full δ:
        δ_full = δ_site .+ (spill_sites[:,1] .* δ_g_0[1]) .+ (spill_sites[:,2] .* δ_g_0[2])

        #   Calculate model-implied choice probabilities with grouped site correction:
        pr_model = nmnl_prob(ch, X, θ, λ, δ_full)

        # Calculate model-implied aggregate shares:
        co = sum(ch, dims = 2)
        pr_model = pr_model .* co
        s_0 = (sum(pr_model, dims = 1) ./ sum(co))[2:end]

        # Calculate model-implied grouped site shares:
        s_g_0 = vcat(sum(s_0 .* spill_sites[:, 1]), sum(s_0 .* spill_sites[:, 2])) 

        # Next iteration δ
        δ_g_1 = δ_g_0 .+ log.(s_g_act) .- log.(s_g_0)
        δ_g_1 = 0.5 .* δ_g_0 .+ 0.5 .* δ_g_1
        eps = maximum(abs.(δ_g_1 .- δ_g_0))
        println("Iteration ", iter, ": eps=", round(eps, digits = 3))
        iter += 1
    end

    δ_full = δ_site .+ (spill_sites[:,1] .* δ_g_1[1]) .+ (spill_sites[:,2] .* δ_g_1[2])

    return δ_g_1
    
end;

#   Function: wtp() ----------------------------------------
#   Description: fxn to calculate lost user-day value from spill
#   Arguments:
#       1. θ = vector containing full array of params
#       2. λ = vector containing nest parameter
#       3. δ = vector of alternative-specific constants
#       4. ch = array containing choices (wide format)
#       5. X = array containing observable data entering utility (wide format)

function wtp(
    θ, 
    λ, 
    δ1, 
    δ0, 
    ch, 
    X
)

    #   Initialize indirect utility:
    nv = zeros(size(ch))

    #   Fill in indirect utility for inside options:
    for i = 1:size(X, 3)
        nv[:, 2:end] += X[:, :, i] .* θ[i]
    end

    #   Add ASCs for inside options:
    nv1 = nv[:, 2:end] .+ δ1'
    nv0 = nv[:, 2:end] .+ δ0'

    #   Divide by similarity coef:
    nv1 = nv1 ./ λ
    nv0 = nv0 ./ λ

    #   Calculate expected utilities:
    env1 = exp.(nv1)
    env0 = exp.(nv0)
    eu1 = log.(1 .+ sum(env1, dims = 2).^λ)
    eu0 = log.(1 .+ sum(env0, dims = 2).^λ)

    #   Calculate WTP:
    wtp = (sum(ch, dims = 2) .* (eu1 .- eu0)) .* (1/θ[1])
    wtp = sum(wtp)

    #   Calculate change in visitation to affected sites:
    pr1 = nmnl_prob(ch, X, θ, λ, δ1)
    pr0 = nmnl_prob(ch, X, θ, λ, δ0)
    Δ_share = (pr0[:, 2:end] * spill_sites[:, 1]) .+ (pr0[:, 2:end] * spill_sites[:, 2]) .-
                (pr1[:, 2:end] * spill_sites[:, 1]) .+ (pr1[:, 2:end] * spill_sites[:, 2])
    Δ_visits = sum(Δ_share .* sum(ch, dims = 2))
    
    return (wtp, Δ_visits)
    
end;




