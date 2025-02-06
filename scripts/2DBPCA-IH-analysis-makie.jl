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
    using ColorSchemes

    #Pkg.status()
    println("Done loading packages")
end

function read_data(sizes, seat_configs, Ρs, δλs, n_trials; n_learned=4, λ₀=0.5, update=false)

    if update == true
        #* Initialize dataframe for the summary of the simulation data
        data = DataFrame(seat_config=String[], class_size=Int[], ρ=Float64[], λ₀=Float64[], δλ=Float64[], m=Measurement{Float64}[], ttl=Measurement{Float64}[])

        #* Reading the data from each simulation
        for seat_config in seat_configs, class_size in sizes, ρ₀ in Ρs, δλ in δλs

            m_list = []
            num_generations_list = []

            for trial in 1:n_trials
                if seat_config == "random"
                    params_df = CSV.read("./output/2D-Binary-PCA-IH/random-$(class_size)-$(n_learned)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-random-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-fit_params.csv", DataFrame)
                else
                    params_df = CSV.read("./output/2D-Binary-PCA-IH/$(seat_config)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(seat_config)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-fit_params.csv", DataFrame)
                end

                #* Adding the data per simulation to the summary
                num_generations = length(params_df[!, "learned_per_gen"]) #/ class_size^2
                m = params_df[!, "power_fit"][2] #± params_df[!, "power_fit"][4]

                push!(num_generations_list, num_generations)
                push!(m_list, m)
            end
            #* for m and t_max: x̄ ± σ/√(N)
            push!(data,
                (
                    seat_config,
                    class_size,
                    ρ₀,
                    λ₀,
                    δλ,
                    mean(m_list) ± (std(m_list) / sqrt(length(m_list))),
                    mean(num_generations_list) ± (std(num_generations_list) / sqrt(length(num_generations_list)))
                )
            )
        end
        CSV.write("./output/2D-Binary-PCA-IH/analysis/data.csv", data)

    else
        #! Need to set data type for measurments
        column_types = Dict(
            "seat_config" => String,
            "class_size" => Int,
            "ρ" => Float64,
            "λ₀" => Float64,
            "δλ" => Float64,
            "m" => Measurement{Float64},
            "ttl" => Measurement{Float64}
        )
        path = "./output/2D-Binary-PCA-IH/analysis/data.csv"
        data = DataFrame(CSV.File(path, types=column_types))
    end
    return data
end

# function plot_n_dependence(data, lengths, SAs, δλs)

begin
    sizes = [32, 48, 64, 96, 128]
    seat_configs = ["traditional", "inner_corner", "outer_corner", "center", "random"]
    Ρs = collect(0.1:0.1:1.0)
    δλs = collect(0.0:0.1:0.4)
    n_trials = 5

    data = read_data(sizes, seat_configs, Ρs, δλs, n_trials, n_learned=4, update=false)

    # lengths = [32, 48, 64, 96, 128]
    # SAs = ["traditional", "inner_corner"]
    # δλs = collect(0.0:0.2:0.4)
    # Ρs = [0.3, 0.5, 0.7]
end

begin #* Ribbon plot for different ρ₀. Ribbon for range of δλ
    lengths = [32, 48, 64, 96, 128]
    SAs = ["traditional", "inner_corner"]
    δλs = collect(0.0:0.1:0.4)
    Ρs = [0.1, 0.5, 0.9]
    # seat_config = "inner_corner"
    # δλ = 0.1

    color = ColorSchemes.seaborn_colorblind

    SA_label_dict = Dict(
        "traditional" => "T",
        "inner_corner" => "IC",
        "outer_corner" => "OC",
        "center" => "C",
        "random" => "R",
    )

    SA_ls_dict = Dict(
        "traditional" => :dash,
        "inner_corner" => :dot,
        # "outer_corner" => 3,
        # "center" => 4,
        # "random" => 5,
    )

    δλ_ls_dict = Dict(
        0.0 => :dash,
        0.2 => :dashdot,
        0.4 => :dashdotdot
    )

    rho_alpha_dict = Dict(collect(0.1:0.1:1.0)[i] => collect(range(0.5, stop=1.0, length=10))[i] for i in 1:length(collect(0.1:0.1:1.0)))

    fig = Figure(size=(1000, 1000 * 0.6); dpi=300)
    ax = Axis(fig[1, 1];
        xlabel="Class size (N)",
        # xlabel = "Class size (L)",
        ylabel="Time to learn (tₘₐₓ)",
        title="tₘₐₓ vs N for inhomogenous model",
        # title = "tₘₐₓ vs L for inhomogenous model",
        # subtitle = "ρ₀=$ρ₀",
        xlabelsize=24,
        ylabelsize=24,
        titlesize=24,
        subtitlesize=16,
        xticklabelsize=16,
        yticklabelsize=16,
        yscale=log10
    )

    # ρ₀ = 0.2
    # seat_config = "inner_corner"

    color_count = 1
    for seat_config in SAs, ρ₀ in Ρs
        ttl_maxes = []
        ttl_mins = []
        ttl_med = []
        for δλ in δλs
            subset = data[(data.seat_config.==seat_config).&(data.δλ.==δλ).&(data.ρ.∈ρ₀), :]

            ttl_max = Measurements.value.(subset.ttl) .+ Measurements.uncertainty.(subset.ttl)
            ttl_min = Measurements.value.(subset.ttl) .- Measurements.uncertainty.(subset.ttl)

            push!(ttl_maxes, ttl_max)
            push!(ttl_mins, ttl_min)

            if δλ == 0.2
                push!(ttl_med, Measurements.value.(subset.ttl))
            end

        end

        lines!(ax, sizes .^ 2, ttl_med[1],
            label="$(SA_label_dict[seat_config]), ρ₀=$ρ₀",
            linestyle=SA_ls_dict[seat_config],
            color=(color[color_count], rho_alpha_dict[ρ₀])
        )

        band!(ax, sizes .^ 2, maximum(ttl_maxes, dims=1)[1], minimum(ttl_mins, dims=1)[1],
            color=(color[color_count], rho_alpha_dict[ρ₀] / 8)
        )

        color_count = color_count + 1

        println(color_count)
    end

    Legend(fig[1, 2], ax,
        "Legend",
        titlesize=16,
        labelsize=16
    )

    save("./output/2D-Binary-PCA-IH/analysis/plots/n-t-plots/n-t-plot-ribbon.png", fig)

    fig
end

begin #* Plots for different ρ₀. Ribbon for ±σ
    lengths = [32, 48, 64, 96, 128]
    SAs = ["traditional", "inner_corner"]
    δλs = collect(0.0:0.2:0.4)
    Ρs = collect(0.1:0.1:1.0)

    color = ColorSchemes.seaborn_colorblind

    SA_label_dict = Dict(
        "traditional" => "T",
        "inner_corner" => "IC",
        "outer_corner" => "OC",
        "center" => "C",
        "random" => "R",
    )

    SA_ms_dict = Dict(
        "traditional" => :rect,
        "inner_corner" => :circle,
        "outer_corner" => 3,
        "center" => 4,
        "random" => 5,
    )

    δλ_color_dict = Dict(
        0.0 => 1,
        0.2 => 2,
        0.4 => 3
    )

    rho_alpha_dict = Dict(collect(0.1:0.1:1.0)[i] => collect(range(0.5, stop=1.0, length=10))[i] for i in 1:length(collect(0.1:0.1:1.0)))

    for ρ₀ in Ρs
        fig = Figure(size=(1000, 1000 * 0.6); dpi=300)
        ax = Axis(fig[1, 1];
            xlabel="Class size (N)",
            # xlabel = "Class size (L)",
            ylabel="Time to learn (tₘₐₓ)",
            title="tₘₐₓ vs N for inhomogenous model",
            # title = "tₘₐₓ vs L for inhomogenous model",
            subtitle="ρ₀=$ρ₀",
            xlabelsize=24,
            ylabelsize=24,
            titlesize=24,
            subtitlesize=16,
            xticklabelsize=16,
            yticklabelsize=16,
            yscale=log10,
            xscale=log10,
        )

        for seat_config in SAs
            for δλ in δλs
                subset = data[(data.seat_config.==seat_config).&(data.δλ.==δλ).&(data.ρ.==ρ₀), :]

                ttl_max = Measurements.value.(subset.ttl) .+ Measurements.uncertainty.(subset.ttl)
                ttl_min = Measurements.value.(subset.ttl) .- Measurements.uncertainty.(subset.ttl)

                scatter!(ax, subset.class_size .^ 2, Measurements.value.(subset.ttl),
                    label="$(SA_label_dict[seat_config]), δλ=0.5±$δλ ",
                    marker=SA_ms_dict[seat_config],
                    color=(color[δλ_color_dict[δλ]], 1.0),
                    markersize = 16,
                )

                # ylims!(ax, 10^1, 10^3)

                lines!(ax, subset.class_size .^ 2, Measurements.value.(subset.ttl),
                    # label="$(SA_label_dict[seat_config]), δλ=0.5±$δλ ",
                    linestyle=:dash,
                    color=(color[δλ_color_dict[δλ]], 1.0)
                )

                band!(ax, subset.class_size .^ 2, ttl_max, ttl_min,
                    color=(color[δλ_color_dict[δλ]], 1.0 / 4)
                )
            end
        end

        Legend(fig[1, 2], ax,
            "Legend",
            titlesize=16,
            labelsize=16
        )

        savepath = "./output/2D-Binary-PCA-IH/analysis/plots/n-t-plots/"
        filename = "n-t-plot-$ρ₀"

        save(savepath * filename * ".png", fig)
    end
    # fig
end

begin #* Plots for different δλ. Ribbon for ±σ
    lengths = [32, 48, 64, 96, 128]
    SAs = ["traditional", "inner_corner"]
    δλs = collect(0.0:0.1:0.4)
    Ρs = [0.1, 0.5, 0.9]

    color = ColorSchemes.seaborn_colorblind

    SA_label_dict = Dict(
        "traditional" => "T",
        "inner_corner" => "IC",
        "outer_corner" => "OC",
        "center" => "C",
        "random" => "R",
    )

    SA_ms_dict = Dict(
        "traditional" => :rect,
        "inner_corner" => :circle,
        "outer_corner" => :diamond,
        "center" => :hexagon,
        "random" => :xcross,
    )

    SA_ls_dict = Dict(
        "traditional" => :dot,
        "inner_corner" => :dash,
        # "outer_corner" => :dashdotdot,
        # "center" => :hexagon,
        # "random" => :xcross,
    )

    rho_color_dict = Dict(
        0.1 => 1,
        0.5 => 2,
        0.9 => 3
    )

    rho_alpha_dict = Dict(collect(0.1:0.1:1.0)[i] => collect(range(0.5, stop=1.0, length=10))[i] for i in 1:length(collect(0.1:0.1:1.0)))

    for δλ in δλs
        fig = Figure(size=(1000, 1000 * 0.6); dpi=300)
        ax = Axis(fig[1, 1];
            xlabel="Class size (N)",
            # xlabel = "Class size (L)",
            ylabel="Time to learn (tₘₐₓ)",
            title="tₘₐₓ vs N for inhomogenous model",
            # title = "tₘₐₓ vs L for inhomogenous model",
            subtitle="δλ=$δλ",
            xlabelsize=24,
            ylabelsize=24,
            titlesize=24,
            subtitlesize=16,
            xticklabelsize=16,
            yticklabelsize=16,
            yscale=log10,
            xscale=log2,
            limits = (nothing, nothing, 10^1, 10^3)
        )

        for seat_config in SAs
            for ρ₀ in Ρs
                subset = data[(data.seat_config.==seat_config).&(data.δλ.==δλ).&(data.ρ.==ρ₀), :]

                ttl_max = Measurements.value.(subset.ttl) .+ Measurements.uncertainty.(subset.ttl)
                ttl_min = Measurements.value.(subset.ttl) .- Measurements.uncertainty.(subset.ttl)

                scatter!(ax, subset.class_size .^ 2, Measurements.value.(subset.ttl),
                    label="$(SA_label_dict[seat_config]), ρ₀=$ρ₀ ",
                    marker=SA_ms_dict[seat_config],
                    color=(color[rho_color_dict[ρ₀]], 1.0),
                    markersize=14
                )

                lines!(ax, subset.class_size .^ 2, Measurements.value.(subset.ttl),
                    # label="$(SA_label_dict[seat_config]), δλ=0.5±$δλ ",
                    linestyle=SA_ls_dict[seat_config],
                    color=(color[rho_color_dict[ρ₀]], 1.0)
                )

                band!(ax, subset.class_size .^ 2, ttl_max, ttl_min,
                    color=(color[rho_color_dict[ρ₀]], 1.0 / 4)
                )
            end
        end

        Legend(fig[1, 2], ax,
            "Legend",
            titlesize=16,
            labelsize=16
        )

        savepath = "./output/2D-Binary-PCA-IH/analysis/plots/n-t-plots-delta/"
        filename = "n-t-plot-$δλ"

        save(savepath * filename * ".png", fig)
    end
end

