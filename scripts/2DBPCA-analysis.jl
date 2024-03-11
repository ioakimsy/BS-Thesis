begin
    using Pkg
    #Pkg.activate(".")

    using CSV
    using DataFrames
    using Plots
    using Statistics
    using Measurements
    using LsqFit

    #Pkg.status()
end

function read_data(sizes, seat_configs, Λs, n_trials; n_learned = 4)

    data = DataFrame(seat_config=String[],class_size=Int[],λ=Float64[],m=Measurement{Float64}[],ttl=Measurement{Float64}[])

    for seat_config in seat_configs, class_size in sizes, λ₀ in Λs
        m_list = []
        num_generations_list = []
        for trial in 1:n_trials
            if seat_config == "random"
                params_df = CSV.read("./output/2D-Binary-PCA/random-$(class_size)-$(n_learned)/$(λ₀)/trial_$(trial)/data/2DBPCA-random-$(class_size)-$(n_learned)-$(λ₀)-trial_$(trial)-fit_params.csv", DataFrame)
            else
                params_df = CSV.read("./output/2D-Binary-PCA/$(seat_config)-$(class_size)/$(λ₀)/trial_$(trial)/data/2DBPCA-$(seat_config)-$(class_size)-$(λ₀)-trial_$(trial)-fit_params.csv", DataFrame)
            end

            num_generations = length(params_df[!, "learned_per_gen"]) / class_size^2
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

function scale_factor_analysis(data, class_config, λs)
    data.class_size .= data.class_size .^ 2
    size_data = data[(data.seat_config .== class_config) .& ([data.λ[i] in λs for i in 1:length(data.λ)]),:]
    
    preanalysis_plot = scatter((size_data[(size_data.λ .== λs[1]),:].class_size), size_data[(size_data.λ .== λs[1]),:].ttl,
    label = "λ = $(λs[1])",
    xlabel = "Class size (N)",
    ylabel = "Time to learn (tₘₐₓ)",
    scale = :log10,
    title = "$(class_config) tₘₐₓ vs N",
    )
    
    for i in 2:length(λs)
        preanalysis_plot = scatter!((size_data[(size_data.λ .== λs[i]),:].class_size), size_data[(size_data.λ .== λs[i]),:].ttl,
        label = "λ = $(λs[i])",
        )
    end

    
    @. power_model(x,p) = p[1] * x ^ p[2]
    
    for j in 1:length(λs)
        power_fit = curve_fit(power_model, 
            (size_data[(size_data.λ .== λs[j]),:].class_size), 
            Measurements.value.(size_data[(size_data.λ .== λs[j]),:].ttl),
            [0.1,-1],
        )
        α = coef(power_fit)[2]
        preanalysis_plot = plot!(minimum(size_data.class_size):maximum(size_data.class_size), power_model(minimum(size_data.class_size):maximum(size_data.class_size), coef(power_fit)), label = "α = $(round(α, digits = 5))", dpi = 300)
    end

    savefig(preanalysis_plot, "./output/2D-Binary-PCA/analysis/plots/scale_factor_preadjusted.png")
    
    #! new data
    new_data = data
    new_size_data = new_data[(new_data.seat_config .== class_config) .& ([new_data.λ[i] in λs for i in 1:length(new_data.λ)]),:]
    


    analysis_plot = scatter(new_size_data[(new_size_data.λ .== λs[1]),:].class_size, new_size_data[(new_size_data.λ .== λs[1]),:].ttl,
    label = "λ = $(λs[1])",
    xlabel = "Class size (N)",
    ylabel = "Time to learn (tₘₐₓ)",
    scale = :log10,
    title = "$(class_config) tₘₐₓ vs N"
    )
    
    for i in 2:length(λs)
        analysis_plot = scatter!(new_size_data[(new_size_data.λ .== λs[i]),:].class_size, new_size_data[(new_size_data.λ .== λs[i]),:].ttl,
        label = "λ = $(λs[i])",
        dpi = 300,
        )
    end

    savefig(analysis_plot, "./output/2D-Binary-PCA/analysis/plots/scale_factor_adjusted.png")

    CSV.write("./output/2D-Binary-PCA/analysis/scaled_data.csv", new_data)

    return new_size_data
end

begin
    # List of parameters
    lengths = [32,48,64,96,128]
	seat_configs = ["inner_corner","outer_corner","center","random"]
	Λs = collect(0.1:0.1:1)
	steady_state_tolerance = 20
	n_trials = 5
    n_learned = 4
    
    mkpath("./output/2D-Binary-PCA/analysis")
    mkpath("./output/2D-Binary-PCA/analysis/plots")
    data = read_data(lengths, seat_configs, Λs, n_trials; n_learned = 4)
    plot_data(data,lengths)

    new_data = scale_factor_analysis(data,"inner_corner", [0.1,0.5,0.9])
    # show(new_data, allrows=true)

end


new_data.class_size
new_data.class_size .= new_data.class_size .^ 2
show(new_data[(new_data.seat_config .== "inner_corner"), :], allrows=true)

