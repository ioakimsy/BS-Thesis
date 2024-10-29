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

function plot_comparison(initial_conditions)

    color = ColorSchemes.seaborn_colorblind.colors

    fig = Figure(size=(1000,1000*0.6); dpi = 300)

    ax = Axis(fig[1, 1];
        xlabel = "Time step",
        ylabel = "Fraction of learned",
        xscale = log10,
        # yscale = log10,
        title = "Comparison between different class parameters",
        subtitle = "Comparing: $(join(comparison, ", "))",
        # xminorgridvisible = false,
        # yminorgridvisible = false,
        xlabelsize = 24,
        ylabelsize = 24,
        titlesize = 24,
        subtitlesize = 16,
        xticklabelsize = 16,
        yticklabelsize = 16,
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
        48 => 1.5,
        64 => 2,
        96 => 2.5,
        128 => 3,
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

    # axislegend(ax, position = :rb)

    Legend(fig[1,2], ax,
        "Legend",
        titlesize = 16,
        labelsize = 16
    )

    return fig
end


begin

    comparison = [
        "size",
        # "SA",
        "ρ₀",
        # "δλ",
    ]

    sizes = in("size", comparison) ? [32, 64, 128] : [64]
    SAs = in("SA", comparison) ? ["traditional", "inner_corner", "outer_corner", "center", "random"] : [ "traditional", "inner_corner"]
    Ρs = in("ρ₀", comparison) ? [0.1, 0.5, 0.9] : [0.5]
    δλ = in("δλ", comparison) ? [0.0, 0.2, 0.4] : [0.2]

    initial_conditions = []
    for class_config in SAs, size in sizes, ρ₀ in Ρs, δλ in δλ
        push!(initial_conditions, [size, class_config, ρ₀, 0.5, δλ])
    end
        
    comparison_plot = plot_comparison(initial_conditions)

    savepath = "./output/2D-Binary-PCA-IH/analysis/plots/trad-PI-learned-t-comparison/"
    filename = join(comparison, "-")

    save(savepath * filename * ".png", comparison_plot)

    comparison_plot
end

