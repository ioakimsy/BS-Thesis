begin
    using Pkg
    #Pkg.activate(".")

    using CSV
    using DataFrames
    using Plots
    using Statistics

    Pkg.status()
end

function read_data(sizes, seat_configs, Λs, n_trials; n_learned = 4)
    data = DataFrame(seat_config=String[],class_size=Int[],λ=Float64[],m=Float64[],num_generations=Float64[])
    for seat_config in seat_configs, class_size in sizes, λ₀ in Λs
        m_list = []
        num_generations_list = []
        for trial in 1:n_trials
            if seat_config == "random"
                params_df = CSV.read("./output/2D-Binary-PCA/random-$(class_size)-$(n_learned)/$(λ₀)/trial_$(trial)/data/2DBPCA-random-$(class_size)-$(n_learned)-$(λ₀)-trial_$(trial)-fit_params.csv", DataFrame)
            else
                params_df = CSV.read("./output/2D-Binary-PCA/$(seat_config)-$(class_size)/$(λ₀)/trial_$(trial)/data/2DBPCA-$(seat_config)-$(class_size)-$(λ₀)-trial_$(trial)-fit_params.csv", DataFrame)
            end

            num_generations = length(params_df[!, "learned_per_gen"])
            m = params_df[!, "power_fit"][2]

            push!(num_generations_list, num_generations)
            push!(m_list, m)
        end
        push!(data, (seat_config, class_size, λ₀, mean(m_list), mean(num_generations_list)))
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

        for seat_config in seat_configs
            _data = size_data[(size_data.seat_config .== seat_config),:]

            _λ = _data[:,:λ]
            _m = _data[:,:m]
            _num_generations = _data[:,:num_generations]

            push!(λs, _λ)
            push!(ms, _m)
            push!(num_generations_list, _num_generations)
        end

        m_plot = scatter(λs,ms, 
            labels = reshape(seat_configs,1,length(seat_configs)),
            ylabel = "Characteristic variable (m)",
            xlabel = "Spread coefficient (λ)",
            title = "Classroom size: $class_size",
            legend = :topleft,
            dpi = 300
        )

        savefig(m_plot, "./output/2D-Binary-PCA/analysis/plots/m-$class_size.png")

        n_plot = scatter(λs, num_generations_list,
        labels = reshape(seat_configs,1,length(seat_configs)),
            ylabel = "Number of generations",
            xlabel = "Spread coefficient (λ)",
            title = "Classroom size: $class_size",
            legend = :topright,
            dpi = 300
        )

        savefig(n_plot, "./output/2D-Binary-PCA/analysis/plots/n-$class_size.png")
    end
end

begin
    # List of parameters
    sizes = [32,64,128]
	seat_configs = ["inner_corner","outer_corner","center","random"]
	Λs = collect(0.1:0.1:1)
	steady_state_tolerance = 20
	n_trials = 3
    n_learned = 4
    
    mkpath("./output/2D-Binary-PCA/analysis")
    mkpath("./output/2D-Binary-PCA/analysis/plots")
    data = read_data(sizes, seat_configs, Λs, n_trials; n_learned = 4)
    plot_data(data,sizes)
end