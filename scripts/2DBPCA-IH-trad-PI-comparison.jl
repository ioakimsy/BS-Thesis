begin
    println("Loading packages")
    using Pkg
    #Pkg.activate(".")

    using CSV
    using DataFrames
    # using Plots
    using CairoMakie
    using ColorSchemes
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

    CairoMakie.activate!(px_per_unit=4)
    println("Done loading packages")
end


function read_time_series_data(class_size, SA, ρ₀, λ₀, δλ; n_learned=4)

    learned_raw = []

    # println("Reading data for $(SA) $(class_size) $(ρ₀) $(λ₀) $(δλ)")

    for trial in 1:5
        if SA == "random"
            path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)-$(n_learned)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
        else
            path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
        end

        params_df = CSV.read(path * "fit_params.csv", DataFrame)
        learned = params_df[!,"learned_per_gen"]

        push!(learned_raw, learned)
    end

    trial_1, trial_2, trial_3, trial_4, trial_5 = learned_raw

    maximum_length = maximum([length(trial_1), length(trial_2), length(trial_3), length(trial_4), length(trial_5)])

    for i in 1:5
        learned_raw[i] = append!(learned_raw[i], ones(maximum_length - length(learned_raw[i])))
    end

    learned_matrix = hcat(learned_raw...)

    learned_μ = vec(mean(learned_matrix, dims=2))
    learned_σ = vec(std(learned_matrix, dims=2))

    # learned = learned_μ .± learned_σ

    return learned_μ, learned_σ
    
    # learned_dom = 1:length(learned)
    # if SA == "traditional"
    #         learned = learned[begin+1:end]
    #         learned_dom = 2:length(learned)+1
    # end
end

function read_time_series_data_raw(class_size, SA, ρ₀, λ₀, δλ; n_learned=4)

    learned_raw = []

    # println("Reading data for $(SA) $(class_size) $(ρ₀) $(λ₀) $(δλ)")

    for trial in 1:5
        if SA == "random"
            path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)-$(n_learned)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
        else
            path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
        end

        params_df = CSV.read(path * "fit_params.csv", DataFrame)
        learned = params_df[!,"learned_per_gen"]

        push!(learned_raw, learned)
    end

    return learned_raw
    
    # learned_dom = 1:length(learned)
    # if SA == "traditional"
    #         learned = learned[begin+1:end]
    #         learned_dom = 2:length(learned)+1
    # end
end

function plot_comparison(initial_conditions)

    color = ColorSchemes.seaborn_colorblind.colors

    fig = Figure(size=(1000,1000 ÷ sqrt(2)); dpi = 300)

    ax = Axis(fig[1, 1];
        xlabel = "Time step",
        ylabel = "Fraction of learned",
        xscale = log10,
        # yscale = log10,
        # title = "Comparison between different class parameters",
        # subtitle = "Comparing: $(join(comparison, ", "))",
        # xminorgridvisible = false,
        # yminorgridvisible = false,
        # limits = (nothing, 10^3, nothing, nothing),
        xlabelsize = 28,
        ylabelsize = 28,
        titlesize = 28,
        subtitlesize = 24,
        xticklabelsize = 24,
        yticklabelsize = 24,
    )

    SA_label_dict = Dict(
        "traditional" => "T",
        "inner_corner" => "IC",
        "outer_corner" => "OC",
        "center" => "C",
        "random" => "R",
    )

    SA_symbol_dict = Dict(
        "traditional" => :circle,
        "inner_corner" => :rect,
        "outer_corner" => :hexagon,
        "center" => :utriangle,
        "random" => :diamond,
    )

    SA_color_dict = Dict(
        "traditional" => 1,
        "inner_corner" => 2,
        "outer_corner" => 3,
        "center" => 4,
        "random" => 5,
    )
    
    size_width_dict = Dict(
        32 => 1,
        48 => 2,
        64 => 3,
        96 => 4,
        128 => 5,
    )

    δλ_ls_dict = Dict(
        0.0 => :dash,
        0.2 => :dashdot,
        0.4 => :dashdotdot
    )

    rho_alpha_dict = Dict(collect(0.1:0.1:1.0)[i] => collect(range(0.5, stop=1.0, length=10))[i] for i in 1:length(collect(0.1:0.1:1.0)))

    
    # No need after changing color to correspond to SA
    # if length(initial_conditions) > 10
    #     color = ColorSchemes.tab20.colors
    # else
    #     color = ColorSchemes.seaborn_colorblind.colors
    # end

    for i in eachindex(initial_conditions)
        class_size, SA, ρ₀, λ₀, δλ = initial_conditions[i]

        learned_μ, learned_σ = read_time_series_data(class_size, SA, ρ₀, λ₀, δλ)

        if SA == "traditional"
            learned_μ = learned_μ[2:end]
            learned_μ_dom = 2:length(learned_μ)+1

            learned_σ = learned_σ[2:end]
            learned_σ_dom = 2:length(learned_σ)+1
        else
            learned_μ_dom = 1:length(learned_μ)
            learned_σ_dom = 1:length(learned_σ)
        end

        learned_min = learned_μ .- learned_σ
        learned_max = learned_μ .+ learned_σ

        band!(ax, learned_μ_dom, learned_max, learned_min;
            # whiskerwidth = 5,
            color = (color[SA_color_dict[SA]], rho_alpha_dict[ρ₀]/4),
        )

        lines!(ax, learned_μ_dom, learned_μ,
            label = "$(SA_label_dict[SA])$(class_size), ρ₀=$(ρ₀), λ=$(λ₀)±$(δλ)",
            color = (color[SA_color_dict[SA]], rho_alpha_dict[ρ₀]),
            # marker = SA_symbol_dict[SA],
            linewidth = size_width_dict[class_size],
            linestyle = δλ_ls_dict[δλ],
        )

    end

    # colsize!(fig.layout, 1, Aspect(1,1.5))

    # axislegend(ax,
    #     position = :rb,
    #     labelsize = 24,
    #     background_color = :transparent,
    #     framevisible = false
    # )
    fig

    Legend(fig[2,1], ax,
        titlesize = 24,
        labelsize = 24
    )

    return fig
end

function return_map(return_map_params) #* Comparing f(t) vs f(t-1)
    # return_map_params = [64, "traditional", 0.5, 0.5, 0.2]
    color = ColorSchemes.seaborn_colorblind.colors

    fig = Figure(size=(1000,1000÷√2); dpi = 300)

    ax = Axis(fig[1, 1];
        xlabel = L"f(t-1)",
        ylabel = L"f(t)",
        # xscale = log10,
        # yscale = log10,
        # title = "Return map of learned fraction",
        # subtitle = "L=$(return_map_params[1]), ρ₀=$(return_map_params[3]), λ=$(return_map_params[4]), δλ=$(return_map_params[5])",
        # xminorgridvisible = false,
        # yminorgridvisible = false,
        xlabelsize = 28,
        ylabelsize = 28,
        titlesize = 28,
        subtitlesize = 24,
        xticklabelsize = 24,
        yticklabelsize = 24,
    )

    ti_params = [Int(return_map_params[1]), "traditional", return_map_params[3], return_map_params[4], return_map_params[5]]
    pi_params = [Int(return_map_params[1]), "inner_corner", return_map_params[3], return_map_params[4], return_map_params[5]]

    learned_data_ti = read_time_series_data_raw(ti_params...)
    learned_data_pi = read_time_series_data_raw(pi_params...)

    for j in 1:2
        learned_data = j == 1 ? learned_data_ti : learned_data_pi
        f_t_data = []
        f_t1_data = []

        for i in eachindex(learned_data)
            f_t = learned_data[i][2:end]
            f_t_1 = learned_data[i][1:end-1]

            push!(f_t_data, f_t...)
            push!(f_t1_data, f_t_1...)
        end

        scatter!(ax, f_t1_data, f_t_data,
            label = j == 1 ? "Traditional" : "Inner corner",
            color = (color[j], 0.5),
            marker = [:circle, :rect, :hexagon, :utriangle, :diamond][j],
            markersize = 16,
        )

        if j == 1
            @. lin_model(x,p) = p[1] * x + p[2]
            lin_fit = curve_fit(lin_model, f_t1_data, f_t_data, [1.0, 0.0])
            lin_fit_params = lin_fit.param

            fit_label = L"f_t=%$(round(lin_fit_params[1], digits=3))f_{t-1}+%$(round(lin_fit_params[2], digits=3))"
            fit_y = lin_model(0:0.01:1, lin_fit_params)

            lines!(ax, 0:0.01:1, fit_y, 
            label = fit_label, 
            color = ColorSchemes.seaborn_colorblind[j],
            linestyle = :dash
            )

        elseif j == 2
            @. quad_model(x,p) = p[1] * x^2 + p[2] * x + p[3]
            quad_fit = curve_fit(quad_model, f_t1_data, f_t_data, [0.0, 1.0, 0.0])
            quad_fit_params = quad_fit.param

            fit_label = L"f_t=%$(round(quad_fit_params[1], digits=3))f_{t-1}^2+%$(round(quad_fit_params[2], digits=3))f_{t-1}+%$(round(quad_fit_params[3], digits=3))"
            fit_y = quad_model(0:0.01:1, quad_fit_params)

            lines!(ax, 0:0.01:1, fit_y, 
            label = fit_label, 
            color = ColorSchemes.seaborn_colorblind[j],
            linestyle = :dash
            )
        end        
    end
    @. logistic_model(x) = x + x * (1 - x)

    lines!(ax, 0:0.01:1, logistic_model(0:0.01:1),
        # label = L"f_t=f_{t-1}+f_{t-1}\dot(1-f_{t-1})",
        color = ColorSchemes.seaborn_colorblind[3],
        linestyle = :solid
    )
    
    lines!(ax, 0:0.01:1, 0:0.01:1,
        # label = L"y=x",
        color = ColorSchemes.seaborn_colorblind[8],
        linestyle = :solid
    )
    
    # axislegend(ax, position = :rb, labelsize = 24, background_color = :transparent, framevisible = false)
    fig

    savepath = "./output/2D-Binary-PCA-IH/analysis/plots/return-map/"
    filename = "return-map-$(return_map_params[1])-$(return_map_params[3])-$(return_map_params[4])-$(return_map_params[5])"

    save(savepath * filename * ".png", fig)
end

begin

    comparison = [
        # "size",
        # "SA",
        # "ρ₀",
        "δλ",
    ]

    sizes = in("size", comparison) ? [32, 64, 128] : [64]
    SAs = in("SA", comparison) ? ["traditional", "inner_corner", "outer_corner", "center", "random"] : [ "traditional", "inner_corner"]
    Ρs = in("ρ₀", comparison) ? [0.1, 0.5, 0.9] : [0.5]
    δλ = in("δλ", comparison) ? [0.0, 0.2, 0.4] : [0.2]

    initial_conditions = []
    for class_config in SAs, size in sizes, ρ₀ in Ρs, δλ in δλ
        push!(initial_conditions, [size, class_config, ρ₀, 0.5, δλ])
    end
    initial_conditions
end

begin #! Return map
    comparison = [
        "size",
        # "SA",
        "ρ₀",
        "δλ",
    ]

    sizes = in("size", comparison) ? [32, 64, 128] : [64]
    SAs = in("SA", comparison) ? ["traditional", "inner_corner", "outer_corner", "center", "random"] : [ "traditional", "inner_corner"]
    Ρs = in("ρ₀", comparison) ? [0.1, 0.5, 0.9] : [0.5]
    δλ = in("δλ", comparison) ? [0.0, 0.2, 0.4] : [0.2]

    initial_conditions = []
    for size in sizes, ρ₀ in Ρs, δλ in δλ
        push!(initial_conditions, [size, "_" , ρ₀, 0.5, δλ])
    end

    initial_conditions

    for i in eachindex(initial_conditions)
        return_map(initial_conditions[i])
    end
end

begin #* Comparing time series traditional with PI
    comparison_plot = plot_comparison(initial_conditions)

    savepath = "./output/2D-Binary-PCA-IH/analysis/plots/trad-PI-learned-t-comparison/"
    filename = join(comparison, "-")

    save(savepath * filename * ".png", comparison_plot)

    comparison_plot
end
