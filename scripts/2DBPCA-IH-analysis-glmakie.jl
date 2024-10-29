begin
    println("Loading packages")
    using Pkg
    #Pkg.activate(".")

    using CSV
    using DataFrames
    using GLMakie
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

# begin
#     sizes = [32, 48, 64, 96, 128]
#     seat_configs = ["traditional", "inner_corner", "outer_corner", "center", "random"]
#     Ρs = collect(0.1:0.1:1.0)
#     δλs = collect(0.0:0.1:0.4)
#     n_trials = 5
#     data = read_data(sizes, seat_configs, Ρs, δλs, n_trials, n_learned=4, update=false)

#     valid_SA = ["inner_corner", "traditional"]
#     valid_class_size = 128

#     subset = data[(data.class_size.==valid_class_size) .& in.(data.seat_config, Ref(valid_SA)), :]

#     SA_subset = groupby(subset, :seat_config)
# end

#! Generate -> save plots and animations
begin
    sizes = [32, 48, 64, 96, 128]
    seat_configs = ["traditional", "inner_corner", "outer_corner", "center", "random"]
    Ρs = collect(0.1:0.1:1.0)
    δλs = collect(0.0:0.1:0.4)
    n_trials = 5
    data = read_data(sizes, seat_configs, Ρs, δλs, n_trials, n_learned=4, update=false)

    valid_SA = ["inner_corner", "traditional"]

    for valid_class_size in sizes
        subset = data[(data.class_size.==valid_class_size) .& in.(data.seat_config, Ref(valid_SA)), :]
        SA_subset = groupby(subset, :seat_config)

        SA_label_dict = Dict(
            "traditional" => "T",
            "inner_corner" => "IC",
            "outer_corner" => "OC",
            "center" => "C",
            "random" => "R",
        )

        fig = Figure(size = (768,768);)
        ax = Axis3(fig[1, 1], 
            xlabel="ρ₀", 
            ylabel = "δλ",
            zlabel="Time to learn (tₘₐₓ)", 
            # zlabel = [],
            title="Inhomogenous Classroom Model L=$valid_class_size", 
            # subtitle="class_size = $valid_class_size",
            aspect = (1,1,1),
            # zticklabelsvisible = false,
            titlesize = 32,
            zticks = 0 : 100 : maximum(Measurements.value.(subset.ttl))
        )

        for i in 1:length(SA_subset)
            ρ_data = SA_subset[i][!, "ρ"]
            δλ_data = SA_subset[i][!, "δλ"]
            t_data = Measurements.value.(SA_subset[i][!, "ttl"])

            points = Point3f.(ρ_data,δλ_data,t_data)
            # t_data = reshape(t_data, 10, 5)'
            # t_σ = Measurements.uncertainty.(SA_subset[1][!, "ttl"])
            
            scatter!(ax, points,
                markersize = 5,
                color = ColorSchemes.seaborn_colorblind[i],
                label = "$(SA_label_dict[valid_SA[i]])",
                strokewidth = 0,
            )

            surface!(ax, ρ_data, δλ_data, t_data,
                # alpha = 0.5,
                color = fill((ColorSchemes.seaborn_colorblind[i], 1), 1:10, 1:5),
                shading = NoShading
            )

            
            # errorbars!(ax, ρ_data, δλ_data, t_data, t_σ)
        end
        ax.viewmode = :fit
        Legend(fig[1,2], ax)
        # display(fig)

        #* Animation
        println("Generating animation: $valid_class_size")
        start_angle = 1.275 * π
        n_frames = 300
        anim_savepath = "./output/2D-Binary-PCA-IH/analysis/plots/rho-dl-t-anim/"
        anim_filename = "$valid_class_size"
        record(fig, anim_savepath * anim_filename * ".mp4", 1:n_frames; framerate = 30) do frame
            ax.azimuth[] = start_angle + 2pi * frame / n_frames
        end

        #* Save the plot
        println("Generating plots: $valid_class_size")
        ax.azimuth = 1.275 * π
        plot_savepath = "./output/2D-Binary-PCA-IH/analysis/plots/rho-dl-t-plots/"
        plot_filename = "$valid_class_size"
        save(plot_savepath * plot_filename * ".png", fig)
    end
end

# #! Generate -> view and explore plot
# begin
#     sizes = [32, 48, 64, 96, 128]
#     seat_configs = ["traditional", "inner_corner", "outer_corner", "center", "random"]
#     Ρs = collect(0.1:0.1:1.0)
#     δλs = collect(0.0:0.1:0.4)
#     n_trials = 5
#     data = read_data(sizes, seat_configs, Ρs, δλs, n_trials, n_learned=4, update=false)

#     valid_SA = ["inner_corner", "traditional"]

#     valid_class_size = 128
#     _subset = data[(data.class_size.==valid_class_size) .& in.(data.seat_config, Ref(valid_SA)), :]
#     SA_subset = groupby(_subset, :seat_config)

#     SA_label_dict = Dict(
#         "traditional" => "T",
#         "inner_corner" => "IC",
#         "outer_corner" => "OC",
#         "center" => "C",
#         "random" => "R",
#     )

#     fig = Figure(size = (768,768);)
#     ax = Axis3(fig[1, 1], 
#         xlabel="ρ₀", 
#         ylabel = "δλ",
#         zlabel="Time to learn (tₘₐₓ)", 
#         # zlabel = [],
#         title="Inhomogenous Classroom Model L=$valid_class_size", 
#         # subtitle="class_size = $valid_class_size",
#         aspect = (1,1,1),
#         # zticklabelsvisible = false,
#         titlesize = 32,
#         zticks = 0 : 100 : maximum(Measurements.value.(_subset.ttl))
#     )

#     for i in 1:length(SA_subset)
#         ρ_data = SA_subset[i][!, "ρ"]
#         δλ_data = SA_subset[i][!, "δλ"]
#         t_data = Measurements.value.(SA_subset[i][!, "ttl"])

#         points = Point3f.(ρ_data,δλ_data,t_data)
#         # t_data = reshape(t_data, 10, 5)'
#         # t_σ = Measurements.uncertainty.(SA_subset[1][!, "ttl"])
        
#         scatter!(ax, points,
#             markersize = 5,
#             color = ColorSchemes.seaborn_colorblind[i],
#             label = "$(SA_label_dict[valid_SA[i]])",
#             strokewidth = 0,
#         )

#         surface!(ax, ρ_data, δλ_data, t_data,
#             # alpha = 0.5,
#             color = fill((ColorSchemes.seaborn_colorblind[i], 1), 1:10, 1:5),
#             shading = NoShading
#         )

        
#         # errorbars!(ax, ρ_data, δλ_data, t_data, t_σ)
#     end
#     ax.viewmode = :fit
#     Legend(fig[1,2], ax)

#     display(fig)
# end