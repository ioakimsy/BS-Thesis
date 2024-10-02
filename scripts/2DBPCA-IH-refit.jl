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

    if SA == "random"
        path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)-$(n_learned)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
    else
        path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
    end

    params_df = CSV.read(path * "fit_params.csv", DataFrame)

    learned = params_df[!,"learned_per_gen"]
    num_generations = size(learned)[1]
    learned_dom = 1:length(learned)
    fit_dom = 1:0.1:length(learned)

    #* if total generations < 200, just use first 10 gens. If greater, use first 10% of the data
    # diff_from_max = length(learned) < 200 ? Int64(length(learned) - 10) : Int64(floor(0.95*length(learned)))
    diff_from_max = Int64(length(learned) - 10)
    
    learned_y = learned[1:end-diff_from_max]
    generation_domain = 1:length(learned_y)

    # @. model(x,p) = p[1] * x ^ p[2]
    @. model(x,p) = 1 - ℯ ^ (-p[1] * x ^ p[2])

    fit = curve_fit(model, generation_domain, learned_y, [0.25,2.0], lower = [0.,0.], upper = [1.,5.])

    power_coeffs = coef(fit)

    fit_vals = model(fit_dom, power_coeffs)

    if SA == "traditional"
            learned = learned[begin+1:end]
            learned_dom = 2:length(learned)+1
    end
end

power_coeffs


begin
    p = plot(dpi = 300, scale = :log10, title = "$SA - $class_size - $ρ₀ - $λ₀ - $δλ - trial $trial")
    scatter!(p, learned_dom, learned)
    plot!(p, fit_dom, fit_vals)
    # if SA == "traditional"
    #     vline!(p, [1 + length(learned) < 200 ? 10 : 1 + length(learned) - Int64(floor(0.95*length(learned)))], ls = :dash, label = "5%")
    # else
    #     vline!(p, [length(learned) < 200 ? 10 : length(learned) - Int64(floor(0.95*length(learned)))], ls = :dash, label = "5%")
    # end
    vline!(p, [10], label = "10")
    vline!(p, [200], label = "200")
end

params_df.power_fit = [power_coeffs...; [missing for _ in 1:num_generations-length(power_coeffs)]]