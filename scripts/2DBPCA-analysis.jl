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
        ttl_maxes = []

        for seat_config in seat_configs
            _data = size_data[(size_data.seat_config .== seat_config),:]

            _λ = _data[:,:λ]
            _m = _data[:,:m]
            _num_generations = _data[:,:ttl]

            push!(λs, _λ)
            push!(ms, _m)
            push!(num_generations_list, _num_generations)
            push!(labels, "$(seat_config)")
            push!(ttl_maxes, maximum(_num_generations))
        end

        m_plot = scatter(λs[1],ms[1], 
            label = labels[1],
            ylabel = "Class learning rate (m)",
            xlabel = "Positional Learning coefficient (ρ₀)",
            labelfontsize = 14,
            titlefontsize = 16,
            title = "Classroom length: $class_size",
            #legend = :topleft,
            dpi = 300,
            #fontfamily = "Times"
        )

        for i in 2:length(ms)

            m_plot = scatter!(λs[i],ms[i], 
                    label = labels[i],
                    title = "L= $class_size, λ=1",
                )

            if labels[i] == "traditional"
                _x_data = λs[i][1:end]
                _y_data = ms[i][1:end]
                print(typeof(_y_data))
                interp_linear = linear_interpolation(0.1:0.1:1, Measurements.value.(_y_data))

                m_plot = plot!(0.1:0.01:1, interp_linear(0.1:0.01:1),
                    label = false,
                    linestyle = :dash,
                    linecolor = i,
                )
            end 
                
        end
        savefig(m_plot, "./output/2D-Binary-PCA/analysis/plots/m-$class_size.png")

        t_plot = scatter(λs[1],num_generations_list[1], 
            label = labels[1],
            ylabel = "Time to learn (tₘₐₓ)",
            xlabel = "Positional Learning coefficient (ρ₀)",
            labelfontsize = 14,
            titlefontsize = 16,
            title = "L= $class_size, λ=1",
            yscale = :log10,
            legend = :bottomleft,
            dpi = 300,
            ylims = (1, 10^(2.5))
            #fontfamily = "Times"
        )

        for i in 2:length(ms)
            t_plot = scatter!(λs[i],num_generations_list[i], 
                label = labels[i],
            )

            if labels[i] == "traditional"
                _x_data = λs[i][1:end]
                _y_data = num_generations_list[i][1:end]
                
                interp_linear = linear_interpolation(0.1:0.1:1.0, Measurements.value.(_y_data))

                t_plot = plot!(0.1:0.01:1, interp_linear(0.1:0.01:1),
                    label = false,
                    linestyle = :dash,
                    linecolor = i,
                )

            end

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
    a = []
    b = []

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
            push!(a, _fit_coefs[1])
            push!(b, _fit_coefs[2])
            _fit_domain = minimum(_data.class_size):maximum(_data.class_size)


            #* Plots data points
            preanalysis_plot = scatter(_x_data, _y_data,
                label = seat_config * " ρ = $(λs[1])",
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
                markercolor = 1,
            )

            # #* Adds space in the legend
            # preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)

            #* Plots the best fit line
            preanalysis_plot = plot!(_fit_domain, power_model(_fit_domain, _fit_coefs),
                label = false,
                linestyle = :dash,
                linecolor = 1
            )

            # * Adds space in the legend
            # preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)

            #color = color + 1
        else
            #* selects a series of data only from 1 lambda
            _x_data = (_data[(_data.λ .== λs[1]),:].class_size)
            _y_data = _data[(_data.λ .== λs[1]),:].ttl

            #* Calculates the best fit line:  y = αx^β 
            _power_fit = curve_fit(power_model, Measurements.value.(_x_data), Measurements.value.(_y_data), [0.1,-1])
            _fit_coefs = coef(_power_fit)
            push!(fit_coefs, _fit_coefs)
            push!(a, _fit_coefs[1])
            push!(b, _fit_coefs[2])
            _fit_domain = minimum(_data.class_size):maximum(_data.class_size)

            #* Plots data points
            preanalysis_plot = scatter!((_data[(_data.λ .== λs[1]),:].class_size), _data[(_data.λ .== λs[1]),:].ttl,
            label = seat_config * " λ = $(λs[1])",
            xlabel = "Class size (N)",
            ylabel = "Time to Learn (tₘₐₓ)",
            scale = :log10,
            #title = "tₘₐₓ vs N",
            markershape = markers[j],
            markercolor = 1,
            )

            #  #* Adds space in the legend
            #  preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)

            #* Plots the best fit line
            preanalysis_plot = plot!(_fit_domain, power_model(_fit_domain, _fit_coefs),
                label = false,
                linestyle = :dash,
                linecolor = 1,
            )

             #* Adds space in the legend
            #  preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)
             
            #color = color + 1
        end
        
        for i in 2:length(λs)
            #* selects a series of data only from 1 lambda
            _x_data = (_data[(_data.λ .== λs[i]),:].class_size)
            _y_data = _data[(_data.λ .== λs[i]),:].ttl

            #* Calculates the best fit line: y = αx^β 
            _power_fit = curve_fit(power_model, Measurements.value.(_x_data), Measurements.value.(_y_data), [0.1,-1])
            _fit_coefs = coef(_power_fit)
            push!(fit_coefs, _fit_coefs)
            push!(a, _fit_coefs[1])
            push!(b, _fit_coefs[2])
            _fit_domain = minimum(_data.class_size):maximum(_data.class_size)

            #* Plots data points
            preanalysis_plot = scatter!((_data[(_data.λ .== λs[i]),:].class_size), _data[(_data.λ .== λs[i]),:].ttl,
            label =seat_config * " ρ = $(λs[i])",
            markershape = markers[j],
            markercolor = i,
            )

            #  #* Adds space in the legend
            #  preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)

            #* Plots the best fit line
            preanalysis_plot = plot!(_fit_domain, power_model(_fit_domain, _fit_coefs),
                label = false,
                linestyle = :dash,
                linecolor = i
            )

             #* Adds space in the legend
            # preanalysis_plot = scatter!([minimum(_x_data)],[minimum(_y_data)], label=" ", ms=0, mc=false, msc=false, ma = 0)

            #color = color + 1
        end
    end

    savefig("./output/2D-Binary-PCA/analysis/plots/N_vs_tmax-traditional-inner_corner.png")
    return size_config_data, a, b

end

function scale_factor_analysis1(data, class_configs, λs)
    markers = [:circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]

    
    @. power_model(x,p) = p[1] * x ^ p[2]
    data.class_size .= data.class_size .^ 2

    size_config_data = data[([data.seat_config[i] in class_configs for i in eachindex(data.seat_config)]) .& ([data.λ[i] in λs for i in 1:length(data.λ)]),:]

    fit_coefs = []
    a = []
    b = []

    #* For each seating arrangement do: 
    for (j, seat_config) in enumerate(class_configs)
        #* Only gets the seating arrangement of interest
        _data = size_config_data[size_config_data.seat_config .== seat_config,:]

        for i in 1:length(λs)
            #* selects a series of data only from 1 lambda
            _x_data = (_data[(_data.λ .== λs[i]),:].class_size)
            _y_data = _data[(_data.λ .== λs[i]),:].ttl

            #* Calculates the best fit line: y = αx^β 
            _power_fit = curve_fit(power_model, Measurements.value.(_x_data), Measurements.value.(_y_data), [0.1,-1])
            _fit_coefs = coef(_power_fit)
            push!(fit_coefs, _fit_coefs)
            push!(a, _fit_coefs[1])
            push!(b, _fit_coefs[2])
            _fit_domain = minimum(_data.class_size):maximum(_data.class_size)

            #* Plots data points
            if i == 1 && j == 1
                preanalysis_plot = scatter(_x_data, _y_data,
                    label = seat_config * " ρ = $(λs[i])",
                    xlabel = "Class size (N)",
                    ylabel = "Time to learn (tₘₐₓ)",
                    xscale = :log2,
                    yscale = :log10,
                    markershape = markers[j],
                    markercolor = i,
                    dpi = 300,
                    legend = :outerright,
                    title = "tₘₐₓ vs N",
                    #fontfamily = "Times"
                )
            else
                preanalysis_plot = scatter!(_x_data, _y_data,
                    label = seat_config * " ρ = $(λs[i])",
                    xlabel = "Class size (N)",
                    ylabel = "Time to learn (tₘₐₓ)",
                    xscale = :log2,
                    yscale = :log10,
                    markershape = markers[j],
                    markercolor = i,
                    #fontfamily = "Times"
                )
            end

            #* Plots the best fit line
            preanalysis_plot = plot!(_fit_domain, power_model(_fit_domain, _fit_coefs),
                label = false,
                linestyle = :dash,
                linecolor = i
            )
        end
    end
    
    savefig("./output/2D-Binary-PCA/analysis/plots/N_vs_tmax-traditional-inner_corner.png")
    return size_config_data, a, b
end

# Calculate R^2 value
function r_squared(y_data, y_fit)
    ȳ = mean(y_data)
    SS_tot = sum((y_data .- ȳ).^2)
    SS_res = sum((y_data .- y_fit).^2)
    return 1 - SS_res/SS_tot
end

begin
    # Define a list of parameters
    lengths = [32,48,64,96,128]  # lengths for the analysis
    seat_configs = ["inner_corner","outer_corner","center","random","traditional"]  # seating configurations
    Λs = collect(0.1:0.1:1)  # range of λ values
    steady_state_tolerance = 20  # tolerance for steady state
    n_trials = 5  # number of trials
    n_learned = 4  # number of learned
    
    # Create directories for output and plots
    mkpath("./output/2D-Binary-PCA/analysis")
    mkpath("./output/2D-Binary-PCA/analysis/plots")
    
    # Read the data using the defined parameters
    data = read_data(lengths, seat_configs, Λs, n_trials; n_learned = 4)
    
    # Plot the data
    # plot_data(data,lengths)
    
    # # Perform scale factor analysis on the data for specific configurations and λs
    # # fit coefs follow the same order as the graph legend
    new_data, as, bs = scale_factor_analysis1(data,["inner_corner","traditional"], [0.1,0.5,0.9])
    
    # # Split the 'a' coefficients for the two seating configurations
    # as_inner_corner = [as[i] for i in 1:Int(length(as)//2)]
    # as_traditional = [as[i] for i in Int(length(as)//2)+1:length(as)]
    
    # # Split the 'b' coefficients for the two seating configurations
    # bs_inner_corner = [bs[i] for i in 1:Int(length(bs)//2)]
    # bs_traditional = [bs[i] for i in Int(length(bs)//2)+1:length(bs)]
    
    # # Calculate the mean and standard error of the 'b' coefficients for the two seating configurations
    # bs_inner_corner_mean = mean(bs_inner_corner) ± std(bs_inner_corner)/sqrt(length(bs_inner_corner))
    # bs_traditional_mean = mean(bs_traditional) ± std(bs_traditional)/sqrt(length(bs_traditional))
end