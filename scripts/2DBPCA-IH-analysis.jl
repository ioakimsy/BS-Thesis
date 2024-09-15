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

    #Pkg.status()
    println("Done loading packages")
end

function read_data(sizes, seat_configs, Ρs, δλs,  n_trials; n_learned = 4, λ₀ = 0.5, update = false)

    if update == true
        #* Initialize dataframe for the summary of the simulation data
        data = DataFrame(seat_config=String[],class_size=Int[],ρ=Float64[],λ₀=Float64[], δλ = Float64[], m=Measurement{Float64}[],ttl=Measurement{Float64}[])

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
                    mean(m_list) ± (std(m_list)/sqrt(length(m_list))), 
                    mean(num_generations_list) ± (std(num_generations_list)/sqrt(length(num_generations_list)))
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
    
function plot_t_data(data, lengths, SAs, δλs)
    for seat_config in SAs, class_size in lengths
        _data = data[((data.class_size .== class_size).& (data.seat_config .== seat_config)),:]
        
        x_data = []
        y_data = []
        labels = []

        for _δλ in δλs
            _x_data = vec(_data.ρ[(_data.δλ .== _δλ), :])
            _y_data = vec((_data.ttl[(_data.δλ .== _δλ), :]))
            push!(x_data, _x_data)
            push!(y_data, _y_data)
            push!(labels, "δλ = $_δλ")
        end
    
        _plot = nothing
        for i in eachindex(x_data)
            if i == 1

                lin_interp = linear_interpolation(x_data[i], Measurements.value.(y_data[i]))

                y_fit = lin_interp(0.1:0.01:1)

                _plot = scatter(x_data[i], y_data[i], 
                    label=labels[i], 
                    xlabel="Learning coefficient (ρ)", 
                    ylabel="Time to Learn (tₘₐₓ)", 
                    title = "Inhomogenous " * seat_config * " L=$class_size",
                    yscale = :log10,
                    color = i,
                    dpi=300
                )
                _plot = plot!(0.1:0.01:1, lin_interp(0.1:0.01:1), 
                    label = false,
                    linestyle = :dash,
                    color = i
                )

            elseif i == length(x_data)
                lin_interp = linear_interpolation(x_data[i], Measurements.value.(y_data[i]))

                y_fit = lin_interp(0.1:0.01:1)

                _plot = scatter!(x_data[i], y_data[i], 
                    label=labels[i],
                    color = i
                )

                _plot = plot!(0.1:0.01:1, y_fit, 
                    label = false,
                    linestyle = :dash,
                    color = i
                )
            else
                _plot = scatter!(x_data[i], y_data[i], 
                    label=labels[i],
                    color = i
                )
            end
        end
        savefig(_plot, "./output/2D-Binary-PCA-IH/analysis/plots/t-plots/t-$(seat_config)-$(class_size).png")
    end
end

function plot_m_data(data, lengths, SAs, δλs)
    for seat_config in SAs, class_size in lengths
        _data = data[((data.class_size .== class_size).& (data.seat_config .== seat_config)),:]
        
        x_data = []
        y_data = []
        labels = []

        for _δλ in δλs
            _x_data = vec(_data.ρ[(_data.δλ .== _δλ), :])
            _y_data = vec((_data.m[(_data.δλ .== _δλ), :]))
            push!(x_data, _x_data)
            push!(y_data, _y_data)
            push!(labels, "δλ = $_δλ")
        end
    
        _plot = nothing
        for i in eachindex(x_data)
            if i == 1

                lin_interp = linear_interpolation(x_data[i], Measurements.value.(y_data[i]))

                y_fit = lin_interp(0.1:0.01:1)

                _plot = scatter(x_data[i], y_data[i], 
                    label=labels[i], 
                    xlabel="Learning coefficient (ρ)", 
                    ylabel="Learning rate (m)", 
                    title = "Inhomogenous " * seat_config * " L=$class_size",
                    #yscale = :log10,
                    color = i,
                    dpi=300
                )
                _plot = plot!(0.1:0.01:1, lin_interp(0.1:0.01:1), 
                    label = false,
                    linestyle = :dash,
                    color = i
                )

            elseif i == length(x_data)
                lin_interp = linear_interpolation(x_data[i], Measurements.value.(y_data[i]))

                y_fit = lin_interp(0.1:0.01:1)

                _plot = scatter!(x_data[i], y_data[i], 
                    label=labels[i],
                    color = i
                )

                _plot = plot!(0.1:0.01:1, y_fit, 
                    label = false,
                    linestyle = :dash,
                    color = i
                )
            else
                _plot = scatter!(x_data[i], y_data[i], 
                    label=labels[i],
                    color = i
                )
            end
        end
        savefig(_plot, "./output/2D-Binary-PCA-IH/analysis/plots/m-plots/m-$(seat_config)-$(class_size).png")
    end
end

function plot_t_data_ribbon(data, lengths, SAs, δλs)
    for class_size in lengths
        _data = data[data.class_size .== class_size, :]
        SA_count = 1
        for SA in SAs
            _data_SA = _data[(_data.seat_config .== SA), :]

            insertcols!(_data_SA, :ttl_min => Measurements.value.(_data_SA.ttl) .- Measurements.uncertainty.(_data_SA.ttl))
            
            insertcols!(_data_SA, :ttl_max => Measurements.value.(_data_SA.ttl) .+ Measurements.uncertainty.(_data_SA.ttl))

            min_ttl = combine(groupby(_data_SA, :ρ), :ttl_min => minimum => :ttl_min)[:,2]
            med_ttl = combine(groupby(_data_SA, :ρ), :ttl => median => :ttl)[:,2]
            max_ttl = combine(groupby(_data_SA, :ρ), :ttl_max => maximum => :ttl_max)[:,2]

            if SA_count == 1
                plot(Ρs, Measurements.value.(med_ttl), 
                    ribbon = (Measurements.value.(med_ttl) .- Measurements.value.(min_ttl), Measurements.value.(max_ttl) .- Measurements.value.(med_ttl)),
                    yscale = :log10,
                    fillalpha = 0.4,
                    label = SA,
                    dpi = 300,
                    title = "Inhomogenous tₘₐₓ vs ρ for L=$class_size",
                    xlabel = "Learning Coefficient (ρ)",
                    ylabel = "Time to Learn (tₘₐₓ)",
                    ylims = (10^1, 10^3),
                    legend = false,
                    ls = SA == "traditional" ? :solid : :dash
                )
            else
                plot!(Ρs, Measurements.value.(med_ttl), 
                    ribbon = (Measurements.value.(med_ttl) .- Measurements.value.(min_ttl), Measurements.value.(max_ttl) .- Measurements.value.(med_ttl)),
                    yscale = :log10,
                    fillalpha = SA == "traditional" ? 0.2 : 0.4,
                    label = SA,
                    ls = SA == "traditional" ? :solid : :dash,
                    
                )
            end

            if SA_count == length(SAs)
                savefig("./output/2D-Binary-PCA-IH/analysis/plots/rho-t-ribbon-plots/rho-t-ribbon-$(class_size).png")
            end

            SA_count += 1
        end
    end
end

function plot_dl_t_data(data, sizes; SAs = ["traditional", "inner_corner"], Ρs = [0.3, 0.5, 0.7])

    markers = [:circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]

    SA_labels = Dict(
        "outer_corner" => "OC",
        "inner_corner" => "IC",
        "center" => "C",
        "random" => "R",
        "traditional" => "T"
    )

    for class_size in sizes
        _data = filter(row -> row.class_size == class_size && row.seat_config in SAs && row.ρ in Ρs, data)
        
        # Initialize the plot
        p = plot(title="L = $class_size", 
            xlabel="Heterogeneity (δλ)", 
            ylabel="Time to learn (tₘₐₓ)", 
            dpi =300,
            yscale = :log10
        )
    
        # Loop through each seating arrangement and rho value
        for (i, SA) in enumerate(SAs)
            for ρ_i in 1:length(Ρs)
                ρ = Ρs[ρ_i]
                # Filter data for the current SA and rho
                filtered_data = filter(row -> row.seat_config == SA && row.ρ == ρ, _data)
                
                # Extract delta lambda and ttl values
                x_values = filtered_data.δλ
                y_values = filtered_data.ttl
                
                # SA Label
                SA_label = SA_labels[SA]

                # Add series to the plot
                plot!(p, x_values, y_values, 
                    label="$SA_label, ρ=$ρ", 
                    marker=markers[i],
                    color = ρ_i,
                    ls = :dash
                )
            end
        end
        savefig(p, "./output/2D-Binary-PCA-IH/analysis/plots/dl-t-plots/dl-t-$(class_size).png")
    end
end

begin
    sizes = [32,48,64,96,128]
	seat_configs = ["outer_corner", "inner_corner", "center", "random", "traditional"]
	Ρs = collect(0.1:0.1:1)
    δλs = collect(0.0:0.1:0.4)
	steady_state_tolerance = 20
	n_trials = 5
    n_learned = 4
    λ₀ = 0.5

    data = read_data(sizes, seat_configs, Ρs, δλs, n_trials; n_learned = n_learned, λ₀ = λ₀, update = false)

    # plot_t_data(data, sizes, seat_configs, δλs)
    # plot_m_data(data, sizes, seat_configs, δλs)
    # plot_t_data_ribbon(data, sizes, seat_configs, δλs)
    plot_dl_t_data(data, sizes; SAs = ["traditional", "inner_corner"], Ρs = [0.3, 0.5, 0.7])
end
