################################################################################
#   SCRIPT:     simulations.jl
#   AUTHOR:     Jacob Bradt (jacob.bradt@mccombs.utexas.edu)
#   NOTES:      julia functions to run simulations and generate data for the 
#               simulation results reported in "Hotelling Meets Wright: Spatial 
#               Sorting and Measurement Error in Recreation Demand Models"
################################################################################

#   Import dependencies:
using LinearAlgebra, Optim, ForwardDiff, LineSearches, NLSolversBase, Distributions, Random, Plots, StatsBase

#   Function: gen_data() ----------------------------------------
#   Description: fxn to generate choice data
#   Arguments:
#       1. n_ch = number of choice occasions
#       2. n_alt = number of alternatives
#       3. d_ξ = distribution of idiosyncratic taste
#       4. θ = vector of utility parameters
#       5. γ = vector of cost parameters
#       6. n_groups = number of groups of consumers
#       7. nonlinear_c = type of nonlinearity in cost
#       8. nonlinear_u = type of nonlinearity in utility
#       9. nonlinear_c_z = type of nonlinearity of cost instrument, Z

function gen_data(
    n_ch::Int64, 
    n_alt::Int64, 
    d_ξ::Distribution, 
    θ::Vector{Float64}, 
    γ::Vector{Float64}; 
    n_groups = 10, 
    nonlinear_c = "none", 
    nonlinear_u = "none", 
    nonlinear_c_z = "none",
    rng = MersenneTwister(82220)
)

    #   Generate exogenous variables:
    X = zeros(n_ch, n_alt)
    Z = similar(X)
    for i = eachindex(X)
        X[i] = rand(rng, Uniform(-1,1))
        Z[i] = rand(rng, Uniform(-1,1))
    end
    X_alt = reshape(repeat(rand(rng, Uniform(-1,1), n_alt), n_ch), n_alt, :)'

    #   Generate random tastes:
    ξ = similar(X)
    ε = similar(X)
    for i = eachindex(ε)
        ε[i] = rand(rng, Gumbel(0,1))
        ξ[i] = rand(rng, d_ξ)
    end

    #   Generate group-level taste:
    ξ_group = rand(rng, Normal(-0.0,1.0), n_groups)
    ξ_group_long = repeat(ξ_group, inner = (Int(n_ch/n_groups), n_alt))
    ξ_group_long[:, 1] .= 0.0

    #   Combine group-level and idiosyncratic taste:
    ξ_full = ((1 - γ[4]) .* ξ) .+ (γ[4] .* ξ_group_long)

    #   Generate travel costs:
    c = similar(X)
    for i = eachindex(c)
        if nonlinear_c == "none" && nonlinear_c_z == "none"
            c[i] = γ[1] + γ[2] * Z[i] + γ[3] * ξ_full[i] + rand(rng, Normal(0,1))
        elseif nonlinear_c == "quadratic"
            c[i] = γ[1] + γ[2] * Z[i] + γ[3] * (ξ_full[i] + ξ_full[i]^2) + rand(rng, Normal(0,1))
        elseif nonlinear_c == "exponential" 
            c[i] = γ[1] + γ[2] * Z[i] + γ[3] * exp(ξ_full[i]) + rand(rng, Normal(0,1))
        elseif nonlinear_c_z == "quadratic"
            c[i] = γ[1] + γ[2] * (Z[i] + Z[i]^2) + γ[3] * ξ_full[i] + rand(rng, Normal(0,1))
        elseif nonlinear_c_z == "exponential"
            c[i] = γ[1] + γ[2] * exp(Z[i]) + γ[3] * ξ_full[i] + rand(rng, Normal(0,1))
        end
    end

    #   Construct indirect utility:
    u = similar(X)
    for i = eachindex(u) 
        if nonlinear_u == "none" 
            u[i] = θ[1] * X_alt[i] + θ[2] * X[i] + θ[3] * c[i] + θ[4] * ξ_full[i] + ε[i]
        elseif nonlinear_u == "quadratic"
            u[i] = θ[1] * X_alt[i] + θ[2] * X[i] + θ[3] * c[i] + θ[4] * (ξ_full[i] + ξ_full[i]^2) + ε[i]
        elseif nonlinear_u == "exponential"
            u[i] = θ[1] * X_alt[i] + θ[2] * X[i] + θ[3] * c[i] + θ[4] * exp(ξ_full[i]) + ε[i]
        end
    end

    #   Identify choices:
    ch = zeros(n_ch, n_alt)
    ch[findmax(u, dims = 2)[2]] .= 1.0

    return (ch, c, ξ_full,X, Z)

end;

#   Function: mnl_prob() ----------------------------------------
#   Description: fxn to generate multinomial logit probability
#   Arguments:
#       1. ch = array containing choices (wide format)
#       2. X = array containing observable data entering utility (wide format)
#       3. θ = vector containing full array of utility params
#       4. δ = vector of alternative-specific constants
#       5. ξ_group = vector of group fixed effects

function mnl_prob(
    ch::Array{Float64}, 
    X::Array{Float64}, 
    θ::Array{Float64}; 
    δ = nothing, 
    ξ_group = nothing
)

    #   Initialize indirect utility:
    n_alt = size(ch, 2)
    nv = zeros(size(ch))

    #   Fill in indirect utility for inside options:
    for i = 1:size(X, 3)
        nv += X[:, :, i] .* θ[i]
    end

    #   Add ASCs for inside options:
    if !isnothing(ξ_group) & isnothing(δ)
        n_group = Int(size(ξ_group, 1) / (n_alt - 1))
        ξ_group = reshape(ξ_group, n_alt - 1, n_group)'
        n_ind = Int(size(ch, 1)/(n_group))
        ξ_group = repeat(ξ_group, inner = (n_ind, 1))
        nv[:, 2:end] = nv[:, 2:end] .+ ξ_group
    elseif !isnothing(δ) & isnothing(ξ_group)
        δ_full = vcat([0.0], δ)
        nv = (nv .+ δ_full')
    elseif !isnothing(δ) & !isnothing(ξ_group)
        δ_full = vcat([0.0], δ)
        n_group = size(ξ_group, 1) + 1
        n_ind = Int(size(ch, 1)/(n_group))
        ξ_group = repeat(vcat([0.0], ξ_group), inner = n_ind)
        nv = (nv .+ δ_full')
        nv[:, 2:end] = nv[:, 2:end] .+ ξ_group
    end

    #   Calculate unconditional choice probability:
    env = exp.(nv)
    sum_env = sum(env, dims = 2)
    p_inside = env ./ (sum_env)
    return p_inside

end;

#   Function: mnl_grad() ----------------------------------------
#   Description: fxn to generate multinomial logit probability
#   Arguments:
#       1. ch = array containing choices (wide format)
#       2. X = array containing observable data entering utility (wide format)
#       3. pr = array containing choice probabilities
#       4. δ = vector of alternative-specific constants
#       5. ξ_group = vector of group fixed effects

function mnl_grad(
    ch::Array{Float64}, 
    X::Array{Float64}, 
    pr::Array{Float64}; 
    δ = nothing, 
    ξ_group = nothing
)

    #   Calculate gradient wrt ASCs:
    g_asc = ch .- pr

    #   Calculate gradient wrt to observables:
    g_x = zeros(size(X, 3))
    for i = 1:size(X, 3)
        g_x[i] = sum(X[:, :, i] .* g_asc)
    end

    #   Calculate gradient for group FE, if included:
    if !isnothing(ξ_group) & isnothing(δ)

        #   Get # of groups, choice occasions per group:
        n_alt = size(ch, 2)
        n_group = Int(size(ξ_group, 1) / (n_alt - 1))
        n_ind = Int(size(ch, 1) / (n_group))

        #   Sum over groups:
        g_ξ_group = indicatormat(repeat(collect(1:1:n_group), inner = n_ind)) * g_asc[:, 2:end]
        g_ξ_group = reshape(g_ξ_group', :, 1)

        #   Construct gradient:
        g_x = vcat(g_x, g_ξ_group)

    elseif !isnothing(δ) & isnothing(ξ_group)

        #   Construct gradient:
        g_x = vcat(g_x, sum(g_asc[:, 2:end], dims = 1)')

    elseif !isnothing(δ) & !isnothing(ξ_group)

        #   Get # of groups, choice occasions per group:
        n_alt = size(ch, 2)
        n_group = Int(size(ξ_group, 1) +1 )
        n_ind = Int(size(ch, 1) / (n_group))

        #   Sum over groups:
        g_ξ_group = indicatormat(repeat(collect(1:1:n_group), inner = n_ind)) * sum(g_asc[:, 2:end], dims = 2)

        #   Construct gradient:
        g_x = vcat(g_x, sum(g_asc[:, 2:end], dims = 1)', g_ξ_group[2:end])

    end

    #   Return gradient:
    return g_x

end;


#   Function: est_mnl() ----------------------------------------
#   Description: fxn to estimate multinomial logit model via MLE
#   Arguments:
#       1. ch = array containing choices (wide format)
#       2. X = array containing observable data entering utility (wide format)
#       3. θ_start = vector of starting values for utility params
#       4. δ_start = vector of starting values for alternative-specific constants
#       5. ξ_group_start = vector of starting values for group fixed effects
#       6. max_iter = maximum number of iterations
#       7. gradient = boolean indicating whether to use gradient-based optimization
#       8. verbose = boolean indicating whether to print optimization progress

function est_mnl(
    ch::Array{Float64}, 
    X::Array{Float64}, 
    θ_start::Array{Float64}; 
    δ_start = nothing, 
    ξ_group_start = nothing, 
    max_iter = 1000, 
    gradient = true, 
    verbose = true
)

    #   Specify number of params
    n_param = size(X, 3)
    n_asc = size(ch, 2) - 1

    #   Optimization functions without group FEs:
    if isnothing(δ_start) & isnothing(ξ_group_start)
        
        #   LL function for gradient free method:
        function f1(θ) 
            pr = mnl_prob(ch, X, θ[1:n_param])
            return - sum(ch .* log.(pr))
        end

        #   LL function for gradient-based optimization:
        function mnl_fg1!(F,G, θ)
            
            #   Common probability calc:
            pr = mnl_prob(ch, X, θ[1:n_param])

            #   If gradient called in iteration, calculate;
            if G !== nothing
                G[1:end] = -mnl_grad(ch, X, pr)
            end

            #   If log-likelihood objective function called, calculate:
            if F !== nothing
                ll = sum(ch .* log.(pr))
                return -ll
            end

        end

        #   Optimization:
        if gradient == true
            res = Optim.optimize(Optim.only_fg!(mnl_fg1!), θ_start, LBFGS(), Optim.Options(iterations = max_iter, show_trace = verbose))
        else
            res = Optim.optimize(f1, θ_start, Optim.Options(iterations = max_iter))
        end

    elseif !isnothing(δ_start) &  isnothing(ξ_group_start)
        
        #   LL function for gradient free method:
        function f2(θ) 
            pr = mnl_prob(ch, X, θ[1:n_param], δ = θ[(n_param + 1):end])
            return - sum(ch .* log.(pr))
        end

        #   LL function for gradient-based optimization:
        function mnl_fg2!(F,G, θ)
            
            #   Common probability calc:
            pr = mnl_prob(ch, X, θ[1:n_param], δ = θ[(n_param + 1):end])

            #   If gradient called in iteration, calculate;
            if G !== nothing
                G[1:end] = -mnl_grad(ch, X, pr, δ = θ[(n_param + 1):end])
            end

            #   If log-likelihood objective function called, calculate:
            if F !== nothing
                ll = sum(ch .* log.(pr))
                return -ll
            end

        end

        #   Optimization:
        if gradient == true
            res = Optim.optimize(Optim.only_fg!(mnl_fg2!), vcat(θ_start, δ_start), LBFGS(), Optim.Options(iterations = max_iter, show_trace = verbose))
        else
            res = Optim.optimize(f2, vcat(θ_start, δ_start), Optim.Options(iterations = max_iter))
        end

    elseif isnothing(δ_start) & !isnothing(ξ_group_start)

        #   LL function for gradient free method:
        function f3(θ) 
            pr = mnl_prob(ch, X, θ[1:n_param], ξ_group = θ[(n_param + 1):end])
            return - sum(ch .* log.(pr))
        end

        #   LL function for gradient-based optimization:
        function mnl_fg3!(F,G, θ)
            
            #   Common probability calc:
            pr = mnl_prob(ch, X, θ[1:n_param], ξ_group = θ[(n_param + 1):end])

            #   If gradient called in iteration, calculate;
            if G !== nothing
                G[1:end] = -mnl_grad(ch, X, pr, ξ_group = θ[(n_param + 1):end])
            end

            #   If log-likelihood objective function called, calculate:
            if F !== nothing
                ll = sum(ch .* log.(pr))
                return -ll
            end

        end

        #   Optimization:
        if gradient == true
            res = Optim.optimize(Optim.only_fg!(mnl_fg3!), vcat(θ_start, ξ_group_start), LBFGS(), Optim.Options(iterations = max_iter, show_trace = verbose))
        else
            res = Optim.optimize(f3, vcat(θ_start, ξ_group_start), Optim.Options(iterations = max_iter, show_trace = verbose))
        end
    
    elseif !isnothing(δ_start) & !isnothing(ξ_group_start)

        #   LL function for gradient free method:
        function f4(θ) 
            pr = mnl_prob(ch, X, θ[1:n_param], δ = θ[(n_param + 1):(n_param + n_asc)],  ξ_group = θ[(n_param + n_asc + 1):end])
            return - sum(ch .* log.(pr))
        end

        #   LL function for gradient-based optimization:
        function mnl_fg4!(F,G, θ)
            
            #   Common probability calc:
            pr = mnl_prob(ch, X, θ[1:n_param], δ = θ[(n_param + 1):(n_param + n_asc)],  ξ_group = θ[(n_param + n_asc + 1):end])

            #   If gradient called in iteration, calculate;
            if G !== nothing
                G[1:end] = -mnl_grad(ch, X, pr, δ = θ[(n_param + 1):(n_param + n_asc)],  ξ_group = θ[(n_param + n_asc + 1):end])
            end

            #   If log-likelihood objective function called, calculate:
            if F !== nothing
                ll = sum(ch .* log.(pr))
                return -ll
            end

        end

        #   Optimization:
        if gradient == true
            res = Optim.optimize(Optim.only_fg!(mnl_fg4!), vcat(θ_start, δ_start, ξ_group_start), LBFGS(), Optim.Options(iterations = max_iter, show_trace = verbose, g_tol = 1e-24))
        else
            res = Optim.optimize(f4, vcat(θ_start, δ_start, ξ_group_start), Optim.Options(iterations = max_iter, show_trace = verbose))
        end

    end
    #   Return results:
    return Optim.minimizer(res)
    

end;


#   Function: mc_simulations() ----------------------------------------
#   Description: fxn to run Monte Carlo simulations of different estimators 
#               for different specified data generating processes
#   Arguments:
#       1. n_sim = number of simulations
#       2. n_ch = number of choice occasions
#       3. n_alt = number of alternatives
#       4. d_ξ = distribution of idiosyncratic taste
#       5. θ = vector of utility parameters
#       6. γ = vector of cost parameters
#       7. n_groups = number of groups of consumers
#       8. measurement_error = type of measurement error
#       9. nonlinear_c = type of nonlinearity in cost
#       10. nonlinear_u = type of nonlinearity in utility
#       11. nonlinear_c_z = type of nonlinearity of cost instrument, Z
#       12. f_stat = boolean indicating whether to calculate first stage F-stat

function mc_simulations(
    n_sim::Int64, 
    n_ch::Int64, 
    n_alt::Int64, 
    d_ξ::Distribution, 
    θ::Vector{Float64}, 
    γ::Vector{Float64}; 
    n_groups = 10, 
    measurement_error = nothing, 
    nonlinear_c = "none", 
    nonlinear_u = "none", 
    nonlinear_c_z = "none", 
    f_stat = false,
    rng = MersenneTwister(82220)
)

    #   Get number of threads:
    n_threads = Threads.nthreads()

    #   Initialize objects to store parameter estimates:
    θ_sim_base = [zeros(size(θ, 1) + n_alt - 3, 0) for i=1:n_threads]
    θ_sim_fe = [zeros(size(θ, 1) + n_alt + n_groups - 4, 0) for i=1:n_threads]
    θ_sim_cf = [zeros(size(θ, 1) + n_alt - 2, 0) for i=1:n_threads]

    #   If reqested, initialize object to store first stage F-stat:
    if f_stat
        f_stat_sim = [zeros(1, 0) for i=1:n_threads]
    end

    #   Initialize simulation counter:
    sim_counter = 0

    #   Initialize rngs:
    rngs = rand(rng, 1:10_000, n_sim)

    #   Loop over simulations:
    Threads.@threads for s = 1:n_sim
        
        #   Draw data:
        ch_temp, c_temp, ξ_temp, X_temp, Z_temp = gen_data(n_ch, n_alt, d_ξ, θ, γ, n_groups = n_groups, nonlinear_c = nonlinear_c, nonlinear_u = nonlinear_u, nonlinear_c_z = nonlinear_c_z, rng = MersenneTwister(rngs[s]))

        #   Measurement error?
        if !isnothing(measurement_error)
            
            #   Save real cost data:
            c_temp_real = c_temp

            #   Different measurement error models:
            if measurement_error == "additive"
                c_temp = c_temp_real .+ rand(MersenneTwister(rngs[s] + 1), Normal(0, 0.5), size(c_temp_real))
            elseif measurement_error == "mult"
                c_temp = c_temp_real .+ c_temp_real .* rand(MersenneTwister(rngs[s] + 2), Normal(0, 0.15), size(c_temp_real))
            end
            
        end

        # Estimate first stage:
        #   Convert wide items to long:
        Z_temp_long = reshape(Z_temp, :, 1)
        c_temp_long = reshape(c_temp, :, 1)

        #   Get gamma estimates:
        Z_temp_long = hcat(Z_temp_long, indicatormat(repeat(1:1:n_alt, n_ch))', indicatormat(repeat(collect(1:1:n_groups), inner = Int(n_ch/n_groups)*n_alt))')
        γ_hat = Z_temp_long \ c_temp_long

        #   Get first stage CF residuals:
        μ =  c_temp_long .- (Z_temp_long *  γ_hat)
        
        #   If requested, estimate first stage F-stat:
        if f_stat

            #   Get number of obs:
            n_fs = size(Z_temp_long, 1)

            #   Calculate variance of fs estimate:
            Z_temp_demean = Z_temp_long .- mean(Z_temp_long, dims = 1)
            σ_γ = (1/n_fs) *((1/ (n_fs - 2)) * sum((Z_temp_demean.^2 .* μ.^2))) ./ (mean(Z_temp_demean .^2).^2)

            #   Calculate first=stage F-stat:
            f_stat_sim[Threads.threadid()] = hcat(f_stat_sim[Threads.threadid()], γ_hat[1].^2 ./ σ_γ)

        end

        #   Construct second stage regressors with and w/o FS residuals:
        μ = reshape(μ, :, n_alt)
        X_full_fs_temp = reduce((x,y) -> cat(x, y, dims = 3), [[c_temp] [X_temp] [μ]])
        X_full_temp = reduce((x,y) -> cat(x, y, dims = 3), [[c_temp] [X_temp]])

        #   Estimate parameters:
        θ_sim_base[Threads.threadid()] = hcat(θ_sim_base[Threads.threadid()] , est_mnl(ch_temp, X_full_temp, ones(size(X_full_temp, 3)), δ_start = ones(n_alt - 1), verbose = false))
        θ_sim_fe[Threads.threadid()] = hcat(θ_sim_fe[Threads.threadid()], est_mnl(ch_temp, X_full_temp, ones(size(X_full_temp, 3)),δ_start = ones(n_alt - 1),  ξ_group_start = ones(n_groups - 1), verbose = false))
        θ_sim_cf[Threads.threadid()] = hcat(θ_sim_cf[Threads.threadid()], est_mnl(ch_temp, X_full_fs_temp, ones(size(X_full_fs_temp, 3)), δ_start = ones(n_alt-1), verbose = false))

        #   Report status:
        sim_counter += 1
        sim_counter_temp = Int(floor(sim_counter/(n_sim/100)))
        print("\r"*raw"-\|/"[sim_counter_temp%4+1]*"="^sim_counter_temp*" "^(100-sim_counter_temp)*"|")
        
    end

    #   Report final results:
    println("Simulations complete.")

    #   Return all requested objects:
    if f_stat
        return (reduce(hcat, θ_sim_base), reduce(hcat, θ_sim_fe), reduce(hcat, θ_sim_cf), reduce(hcat, f_stat_sim))
    else
        return (reduce(hcat, θ_sim_base), reduce(hcat, θ_sim_fe), reduce(hcat, θ_sim_cf))   
    end
    
end;
