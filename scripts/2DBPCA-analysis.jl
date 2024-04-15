begin
    using Pkg
    #Pkg.activate(".")

    using CSV
    using DataFrames
    using Plots
    using Statistics
    using Measurements
    using LsqFit
    using LaTeXStrings

    #Pkg.status()
end

function read_data(sizes, seat_configs, Λs, n_trials; n_learned = 4)

    #* Initialize dataframe for the summary of the simulation data
    data = DataFrame(seat_config=String[],class_size=Int[],λ=Float64[],m=Measurement{Float64}[],ttl=Measurement{Float64}[])

    #* Reading the data from each simulation
    for seat_config in seat_configs, class_size in sizes, λ₀ in Λs
        m_list = []
        num_generations_list = []
        for trial in 1:n_trials
            if seat_config == "random"
                params_df = CSV.read("./output/2D-Binary-PCA/random-$(class_size)-$(n_learned)/$(λ₀)/trial_$(trial)/data/2DBPCA-random-$(class_size)-$(n_learned)-$(λ₀)-trial_$(trial)-fit_params.csv", DataFrame)
            else
                params_df = CSV.read("./output/2D-Binary-PCA/$(seat_config)-$(class_size)/$(λ₀)/trial_$(trial)/data/2DBPCA-$(seat_config)-$(class_size)-$(λ₀)-trial_$(trial)-fit_params.csv", DataFrame)
            end

            #* Adding the data per simulation  to the summary
            num_generations = length(params_df[!, "learned_per_gen"]) #/ class_size^2
            m = params_df[!, "power_fit"][2] #± params_df[!, "power_fit"][4]

            push!(num_generations_list, num_generations)
            push!(m_list, m)
        end
        push!(data, 
            (
                seat_config, 
                class_size, 
                λ₀, 
                mean(m_list) ± (std(m_list)/sqrt(length(m_list))), 
                mean(num_generations_list) ± (std(num_generations_list)/sqrt(length(num_generations_list)))
            )
        )
    end
    CSV.write("./output/2D-Binary-PCA/analysis/data.csv", data)
    return data
end

function plot_data(data, sizes)
    for class_size in sizes

        size_data = data[(data.class_size .== class_size),:]

        λs = []
        ms = []
        num_generations_list = []
        labels = []

        for seat_config in seat_configs
            _data = size_data[(size_data.seat_config .== seat_config),:]

            _λ = _data[:,:λ]
            _m = _data[:,:m]
            _num_generations = _data[:,:ttl]
            

            push!(λs, _λ)
            push!(ms, _m)
            push!(num_generations_list, _num_generations)
            push!(labels, "$(seat_config)")
        end

        m_plot = scatter(λs[1],ms[1], 
            label = labels[1],
            ylabel = "Characteristic variable (m)",
            xlabel = "Spread coefficient (λ)",
            title = "Classroom length: $class_size",
            #legend = :topleft,
            dpi = 300
        )

        for i in 2:length(ms)
            m_plot = scatter!(λs[i],ms[i], 
                label = labels[i],
            )
        end
        savefig(m_plot, "./output/2D-Binary-PCA/analysis/plots/m-$class_size.png")

        t_plot = scatter(λs[1],num_generations_list[1], 
            label = labels[1],
            ylabel = "Time to Learn (tₘₐₓ)",
            xlabel = "Spread coefficient (λ)",
            title = "Classroom length: $class_size",
            #legend = :topright,
            dpi = 300
        )

        for i in 2:length(ms)
            t_plot = scatter!(λs[i],num_generations_list[i], 
                label = labels[i],
            )
        end

        savefig(t_plot, "./output/2D-Binary-PCA/analysis/plots/t-$class_size.png")
    end
end

function scale_factor_analysis(data, class_configs, λs)

    markers = [:circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]

    @. power_model(x,p) = p[1] * x ^ p[2]
    data.class_size .= data.class_size .^ 2

    size_config_data = data[([data.seat_config[i] in class_configs for i in eachindex(data.seat_config)]) .& ([data.λ[i] in λs for i in 1:length(data.λ)]),:]
    
    fit_coefs = []

    color = 1

    #* Sets up plotting to have one marker for each seating arrangement.
    #* For each seating arrangement do: 
    for j in eachindex(class_configs)

        #* Sets up which seating arrangement goes first based on known list
        seat_config = class_configs[j]

        #* Only gets the seating arrangement of interest
        _data = size_config_data[[size_config_data.seat_config[i] .== seat_config for i in eachindex(size_config_data.seat_config)],:]

        if j == 1
            #* selects a series of data only from 1 lambda
            _x_data = (_data[(_data.λ .== λs[1]),:].class_size)
            _y_data = _data[(_data.λ .== λs[1]),:].ttl

            #* Calculates the best fit line: y = αx^β 
            _power_fit = curve_fit(power_model, Measurements.value.(_x_data), Measurements.value.(_y_data), [0.1,-1])
            _fit_coefs = coef(_power_fit)
            push!(fit_coefs, _fit_coefs)
            _fit_domain = minimum(_data.class_size):maximum(_data.class_size)

            #* Plots data points
            preanalysis_plot = scatter(_x_data, _y_data,
                label = seat_config * " λ = $(λs[1])",
                xlabel = "Class size (N)",
                ylabel = "Time to learn (tₘₐₓ)",
                scale = :log10,
                title = "tₘₐₓ vs N",
                markershape = markers[j],
                dpi = 300,
                legend = :outerright,
                #legend_title = "Legend",
                legend_columns = 1,
                legendfontsize = 6,
                markercolor = color,
            )

            #* Adds space in the legend
            preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)

            #* Plots the best fit line
            preanalysis_plot = plot!(_fit_domain, power_model(_fit_domain, _fit_coefs),
                label = L"y = %$(round(_fit_coefs[1], digits=3)) \cdot x^{%$(round(_fit_coefs[2], digits=3))}",
                linestyle = :dash,
                linecolor = color
            )

            #* Adds space in the legend
            preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)

            color = color + 1
        else
            #* selects a series of data only from 1 lambda
            _x_data = (_data[(_data.λ .== λs[1]),:].class_size)
            _y_data = _data[(_data.λ .== λs[1]),:].ttl

            #* Calculates the best fit line:  y = αx^β 
            _power_fit = curve_fit(power_model, Measurements.value.(_x_data), Measurements.value.(_y_data), [0.1,-1])
            _fit_coefs = coef(_power_fit)
            push!(fit_coefs, _fit_coefs)
            _fit_domain = minimum(_data.class_size):maximum(_data.class_size)

            #* Plots data points
            preanalysis_plot = scatter!((_data[(_data.λ .== λs[1]),:].class_size), _data[(_data.λ .== λs[1]),:].ttl,
            label = seat_config * " λ = $(λs[1])",
            xlabel = "Class size (N)",
            ylabel = "Time to learn (tₘₐₓ)",
            scale = :log10,
            #title = "tₘₐₓ vs N",
            markershape = markers[j],
            markercolor = color,
            )

             #* Adds space in the legend
             preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)

            #* Plots the best fit line
            preanalysis_plot = plot!(_fit_domain, power_model(_fit_domain, _fit_coefs),
                label = L"y = %$(round(_fit_coefs[1], digits=3)) \cdot x^{%$(round(_fit_coefs[2], digits=3))}",
                linestyle = :dash,
                linecolor = color,
            )

             #* Adds space in the legend
             preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)
             
            color = color + 1
        end
        
        for i in 2:length(λs)
            #* selects a series of data only from 1 lambda
            _x_data = (_data[(_data.λ .== λs[i]),:].class_size)
            _y_data = _data[(_data.λ .== λs[i]),:].ttl

            #* Calculates the best fit line: y = αx^β 
            _power_fit = curve_fit(power_model, Measurements.value.(_x_data), Measurements.value.(_y_data), [0.1,-1])
            _fit_coefs = coef(_power_fit)
            push!(fit_coefs, _fit_coefs)
            _fit_domain = minimum(_data.class_size):maximum(_data.class_size)

            #* Plots data points
            preanalysis_plot = scatter!((_data[(_data.λ .== λs[i]),:].class_size), _data[(_data.λ .== λs[i]),:].ttl,
            label =seat_config * " λ = $(λs[i])",
            markershape = markers[j],
            markercolor = color,
            )

             #* Adds space in the legend
             preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)

            #* Plots the best fit line
            preanalysis_plot = plot!(_fit_domain, power_model(_fit_domain, _fit_coefs),
                label = L"y = %$(round(_fit_coefs[1], digits=3)) \cdot x^{%$(round(_fit_coefs[2], digits=3))}",
                linestyle = :dash,
                linecolor = color
            )

             #* Adds space in the legend
             preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)

            color = color + 1
        end
    end

    savefig("./output/2D-Binary-PCA/analysis/plots/N_vs_tmax-traditional-inner_corner.png")
    return size_config_data

end

begin
    # List of parameters
    lengths = [32,48,64,96,128]
	seat_configs = ["inner_corner","outer_corner","center","random","traditional"]
	Λs = collect(0.1:0.1:1)
	steady_state_tolerance = 20
	n_trials = 5
    n_learned = 4
    
    mkpath("./output/2D-Binary-PCA/analysis")
    mkpath("./output/2D-Binary-PCA/analysis/plots")
    data = read_data(lengths, seat_configs, Λs, n_trials; n_learned = 4)
    plot_data(data,lengths)

    new_data = scale_factor_analysis(data,["inner_corner","traditional"], [0.1,0.5,0.9])
    show(new_data, allrows=true)

end