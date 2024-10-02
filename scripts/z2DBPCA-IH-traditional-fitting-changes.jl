begin
    println("Loading packages")
    using Pkg
    #Pkg.activate(".")

    using CSV
    using DataFrames
    using Plots
    using Statistics
    using Measurements
    using LsqFit
    using LaTeXStrings
    using Interpolations
    using Random
    using DelimitedFiles
    using Alert
    using ProgressMeter

    using Distributions
    using Integrals
    using Optim
    #Pkg.status()
    println("Done loading packages")
end

begin
    class_size = 64
    SA = "traditional"
    λ₀ = 0.5
    ρ₀ = 0.3
    δλ = 0.4
    trial = 1

    path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"

    params_df = CSV.read(path * "fit_params.csv", DataFrame)

    learned = params_df[!,"learned_per_gen"]
    learned_dom = 1:length(learned)
    fit_dom = 1:0.1:length(learned)

    learned_y = learned[1:end-Int64(floor(0.75*length(learned)))]
    @. model(x, p) = p[1] * ℯ ^ (-p[2] * x) + 1
    fit = curve_fit(model, 1:length(learned_y), learned_y, [-1.0,1.0])
    coeffs = coef(fit)

    fit_vals = model(fit_dom, coeffs)

    if SA == "traditional"
            learned = learned[begin+1:end]
            learned_dom = 1:length(learned)
    end

end


begin #* Trying to use integral solvers to fit the incomplete beta function
    function solve_incomplete_beta(x, p)
        coeffs = p
        f(t, p) = real((t+0im) .^(p[1] .- 1) .* ((1 .- t)+0im) .^ (p[2] .- 1))
        _domain = (0, x)
        int_prob = IntegralProblem(f, _domain, coeffs)
        sol = solve(int_prob, QuadGKJL())
        return sol.u
    end

    function loss_function(p, data_x, data_y)
        predicted_y = [real(solve_incomplete_beta(x, p)) for x in data_x]
        return sum((predicted_y .- data_y) .^ 2)
    end

    initial_guess = [1.0, 1.0]

    result = optimize(p -> loss_function(p, 1:length(learned_y), learned_y), initial_guess)

    coeffs = result.minimizer

    p = scatter(1:length(learned), learned, label="Data", ylims=(0,1))
    p = plot!(1:length(learned), [real(solve_incomplete_beta(x, coeffs)) for x in 1:length(learned)], label="Fit")
    # p = vline!([length(learned_y)], label="75%")
end

begin
    p = scatter(learned_dom, learned, label="Data",
    yscale = :log10
    )
    p = plot!(fit_dom, fit_vals, label="Fit")
    
    fit_end_val = length(learned)-Int64(floor(0.75*length(learned)))
    p = vline!([fit_end_val], label="25%")
    
end

coeffs