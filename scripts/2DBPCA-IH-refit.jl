begin
    println("Loading packages")
    using Pkg
    #Pkg.activate(".")

    using CSV
    using DataFrames
    using CairoMakie
    using Statistics
    using Measurements
    using LsqFit
    using LaTeXStrings
    using Interpolations
    using Random
    using DelimitedFiles
    using Alert
    using ProgressMeter
    using ColorSchemes
    #Pkg.status()
    println("Done loading packages")
end

# begin
#     class_size = 64
#     SA = "traditional"
#     λ₀ = 0.5
#     ρ₀ = 0.9
#     δλ = 0.4
#     trial = 1

#     pts_to_fit = 10

#     if SA == "random"
#         path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)-$(n_learned)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
#     else
#         path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
#     end

#     params_df = CSV.read(path * "fit_params.csv", DataFrame)

#     learned = params_df[!,"learned_per_gen"]
#     num_generations = size(learned)[1]
#     learned_dom = 1:length(learned)
#     fit_dom = 1:0.1:length(learned)

#     #* if total generations < 200, just use first 10 gens. If greater, use first 10% of the data
#     # diff_from_max = length(learned) < 200 ? Int64(length(learned) - 10) : Int64(floor(0.95*length(learned)))
#     diff_from_max = Int64(length(learned) - pts_to_fit)
    
#     learned_y = learned[1:end-diff_from_max]
#     generation_domain = 1:length(learned_y)

#     # @. model(x,p) = p[1] * x ^ p[2]
#     # @. model(x,p) = p[3] * (x + p[1]) ^ p[2]
#     @. model(x,p) = p[2] * log(p[1] * x)

#     # fit = curve_fit(model, generation_domain, learned_y, [0.8, 0.3])

#     power_coeffs = coef(fit)

#     # power_coeffs = [1.5, 3, 5.0e-5]
#     power_coeffs = [0.5, 0.63]

#     fit_vals = model(fit_dom, power_coeffs)

#     if SA == "traditional"
#             learned = learned[begin+1:end]
#             learned_dom = 2:length(learned)+1
#     end
# end

#! Generalized logistic function
begin
    count = 0
    class_size = 64
    # SA = "inner_corner"
    λ₀ = 0.5
    ρ₀ = 0.9
    δλ = 0.4
    trial = 1

    fig = Figure(size = (1000, 1000÷sqrt(2)); dpi = 300)
    ax = Axis(fig[1,1]; 
        xlabel = "Time step", 
        ylabel = "Learned per generation", 
        title = "Learned per generation over time", 
        xscale = log10,
        # yscale = log10,
        ylabelsize = 24,
        xlabelsize = 24,
        titlesize = 32,
        subtitle = "L = $(class_size), λ₀ = $(λ₀), ρ₀ = $(ρ₀), δλ = $(δλ), trial = $(trial)",
        subtitlesize = 24
    )

    colors = ColorSchemes.seaborn_colorblind.colors

    for SA in ["traditional", "inner_corner"]

        # pts_to_fit = 10

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

        # @. model(x,p) = p[1] * x ^ p[2]
        # @. model(x,p) = p[3] * (x + p[1]) ^ p[2]

        @. model(x,p) = 1 / (1 + exp(-p[1] * (x - p[2]))) #* p[1] = steepness, p[2] = midpoint

        # @. model(x,p) = 1 / ((1 + (p[3]) * exp(-p[1] * (x - p[2]))) ^ (1/p[3])) #* p[1] = steepness, p[2] = midpoint, p[3] = shape parameter
        
        
        fit = curve_fit(model, learned_dom, learned, [0.2, 5.0, 1.0])
        fit_coeffs = coef(fit)

        print(fit_coeffs)

        # power_coeffs = [1.5, 3, 5.0e-5]
        # fit_coeffs = [0.5, 0.63]

        fit_vals = model(fit_dom, fit_coeffs)

        if SA == "traditional"
                learned = learned[begin+1:end]
                learned_dom = 2:length(learned)+1
        end
    # end

    # begin
        

        scatter!(ax, learned_dom, learned, 
            # markersize = 5, 
            label = SA, 
            color = colors[count+1]
        )        

        lines!(ax, fit_dom, fit_vals, 
            # linewidth = 2, 
            color = colors[count+1], 
            linestyle = :dash, 
            label = L"f= \frac{1}{[1 + %$(round(fit_coeffs[3],digits=4)) e^{-%$(round(fit_coeffs[1],digits=4))(t - %$(round(fit_coeffs[2],digits=4)))}]^{%$(round(fit_coeffs[3],digits=4))}}"
        )

        # vlines!(ax, [fit_coeffs[2]], linestyle = :dash, linewidth = 2, color = colors[count+3], label = "Midpoint")

        count = count + 1
    end
    Legend(fig[2,:], ax, 
        position = :lt, 
        vertical = false, 
        orientation = :horizontal, 
        halign = :center, 
        valign = :center, 
        nbanks = 2,
        labelsize = 24
    )
    fig
end

#! Power law function
begin
    count = 0
    class_size = 64
    # SA = "inner_corner"
    λ₀ = 0.5
    ρ₀ = 0.9
    δλ = 0.4
    trial = 1

    fig = Figure(size = (1000, 1000÷sqrt(2)); dpi = 300)
    ax = Axis(fig[1,1]; 
        xlabel = "Time step", 
        ylabel = "Learned per generation", 
        title = "Learned per generation over time", 
        xscale = log10,
        # yscale = log10,
        ylabelsize = 24,
        xlabelsize = 24,
        titlesize = 32,
        subtitle = "L = $(class_size), λ₀ = $(λ₀), ρ₀ = $(ρ₀), δλ = $(δλ), trial = $(trial)",
        subtitlesize = 24
    )

    colors = ColorSchemes.seaborn_colorblind.colors

    for SA in ["traditional", "inner_corner"]

        # pts_to_fit = 10

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

        @. model(x,p) = p[1] * x ^ p[2]
        # @. model(x,p) = p[3] * (x + p[1]) ^ p[2]
        # @. model(x,p) = 1 / (1 + exp(-p[1] * (x - p[2]))) #* p[1] = steepness, p[2] = midpoint
        
        fit = curve_fit(model, learned_dom, learned, [5.0, 1.0])
        fit_coeffs = coef(fit)

        # power_coeffs = [1.5, 3, 5.0e-5]
        # fit_coeffs = [0.5, 0.63]

        fit_vals = model(fit_dom, fit_coeffs)

        if SA == "traditional"
                learned = learned[begin+1:end]
                learned_dom = 2:length(learned)+1
        end
    # end

    # begin
        

        scatter!(ax, learned_dom, learned, 
            # markersize = 5, 
            label = SA, 
            color = colors[count+1]
        )        

        lines!(ax, fit_dom, fit_vals, 
            # linewidth = 2, 
            color = colors[count+1], 
            # linewidth = 3, 
            linestyle = :dash, 
            label = L"f= %$(round(fit_coeffs[1],digits=4))x^{%$(round(fit_coeffs[2],digits=4))}"
        )

        # vlines!(ax, [fit_coeffs[2]], linestyle = :dash, linewidth = 2, color = colors[count+3], label = "Midpoint")

        count = count + 1
    end
    Legend(fig[2,:], ax, 
        position = :lt, 
        vertical = false, 
        orientation = :horizontal, 
        halign = :center, 
        valign = :center, 
        nbanks = 1,
        labelsize = 24
    )    
    fig
end

begin
    class_size = 64
    SA = "inner_corner"
    λ₀ = 0.5
    ρ₀ = 0.9
    δλ = 0.4
    trial = 1

    pts_to_fit = 10

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
    diff_from_max = Int64(length(learned) - pts_to_fit)
    
    learned_y = learned[1:end-diff_from_max]
    generation_domain = 1:length(learned_y)

    # @. model(x,p) = p[1] * x ^ p[2]
    # @. model(x,p) = p[3] * (x + p[1]) ^ p[2]
    @. model(x,p) = p[2] * log(p[1] * x)

    # fit = curve_fit(model, generation_domain, learned_y, [0.8, 0.3])

    power_coeffs = coef(fit)

    # power_coeffs = [1.5, 3, 5.0e-5]
    power_coeffs = [0.5, 0.63]

    fit_vals = model(fit_dom, power_coeffs)

    if SA == "traditional"
            learned = learned[begin+1:end]
            learned_dom = 2:length(learned)+1
    end
    scatter!(p, learned_dom, learned)
end

fit_vals

params_df.power_fit = [power_coeffs...; [missing for _ in 1:num_generations-length(power_coeffs)]]
