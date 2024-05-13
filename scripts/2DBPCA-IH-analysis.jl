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

function read_data(sizes, seat_configs, Ρs, δλs,  n_trials; n_learned = 4, λ₀ = 0.5)

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
    return data
end

begin
    sizes = [32,48,64,96,128]
	seat_configs = ["outer_corner", "inner_corner", "center", "random"]
	Ρs = collect(0.1:0.1:1)
    δλs = collect(0.0:0.1:0.4)
	steady_state_tolerance = 20
	n_trials = 5
    n_learned = 4
    λ₀ = 0.5

    data = read_data(sizes, seat_configs, Ρs, δλs, n_trials; n_learned = n_learned, λ₀ = λ₀)
end

