begin
    using Pkg

    using BenchmarkTools
    using ProfileView
    using Profile
    using CairoMakie
    using DataFrames
    using CSV
    using LaTeXStrings

    # List of parameters
    sizes = [128]
    seat_configs = ["random"]
    Λs = [0.25,0.5,0.75]
    steady_state_tolerance = 10
    n_trials = 3
    n_learned = 4

    set_theme!(theme_latexfonts())
end



function generate_plots(sizes::Vector{Int}, seat_configs::Vector{String},Λs::Vector{Float64}, steady_state_tolerance::Int, n_trials::Int; n_learned::Int=4)
    println("Starting to generate plots:")
    #* I tried Threads.@threads, but it doesn't seem to be very stable for this application...?
    for trial in 1:n_trials
        for seat_config in seat_configs, λ₀ in Λs, class_size in sizes#, trial in 1:n_trials

            println("$seat_config 	$λ₀ 	$class_size 	$trial")

            if seat_config == "random"
                data_df = CSV.read("./output/2D-Binary-PCA/random-$(class_size)-$(n_learned)/$(λ₀)/trial_$(trial)/data/2DBPCA-random-$(class_size)-$(n_learned)-$(λ₀)-trial_$(trial)-data.csv", DataFrame)

                params_df = CSV.read("./output/2D-Binary-PCA/random-$(class_size)-$(n_learned)/$(λ₀)/trial_$(trial)/data/2DBPCA-random-$(class_size)-$(n_learned)-$(λ₀)-trial_$(trial)-fit_params.csv", DataFrame)
            else
                data_df = CSV.read("./output/2D-Binary-PCA/$(seat_config)-$(class_size)/$(λ₀)/trial_$(trial)/data/2DBPCA-$(seat_config)-$(class_size)-$(λ₀)-trial_$(trial)-data.csv", DataFrame)

                params_df = CSV.read("./output/2D-Binary-PCA/$(seat_config)-$(class_size)/$(λ₀)/trial_$(trial)/data/2DBPCA-$(seat_config)-$(class_size)-$(λ₀)-trial_$(trial)-fit_params.csv", DataFrame)
            end

            class_size = Int(sqrt(size(data_df)[1])) #! Assumed square classroom
            num_generations = size(data_df)[2]

            generations = [data_df[!,"Generation $(i)"] for i in 1:num_generations]
            generations = [reshape(generations[i], (class_size,class_size)) for i in 1:num_generations]

            learned = params_df[!,"learned_per_gen"]
            learned_dom = 1:length(learned)
            power_coeff, power_exp = params_df[!,"power_fit"][1:2]
            fit_dom = 1:0.1:length(learned)
            fit_vals = power_coeff .* fit_dom .^ power_exp

            for i in 1:num_generations
                fig, class_ax, class = heatmap(generations[i], 
                    colormap=reverse(cgrad(:grays)),
                    colorrange=(0,1),
                    fxaa = false
                    ;
                    axis = (; title = "Classroom over time", 
                        subtitle = L"params",
                        yreversed = true, 
                        aspect = 1),
                    figure = (;dpi=300)
                )
                hidedecorations!(class_ax)

                colsize!(fig.layout, 1, Aspect(1,1))
                Colorbar(fig[1,2], class, label = "Aptitude", vertical = true)


                learned_ax = Axis(fig[2,1:2]; 
                    #title = "Learned over time", 
                    ylabel = "Fraction of learned", 
                    xlabel = "Generation number",
                    xscale = log10,
                    yscale = log10,
                    xminorgridvisible = true,
                    yminorgridvisible = true,
                    xticklabelsize = 12,
                    yticklabelsize = 12,
                    limits = (1,nothing,nothing,nothing)
                )

                learned_data = scatter!(learned_ax, learned_dom, learned,
                    markersize = 3,
                    label = "Learned students"
                )

                fit_data = lines!(learned_ax, fit_dom, fit_vals,
                    label = L"y = %$(round(power_coeff, digits = 5))x^{%$(round(power_coeff, digits = 5))}",
                    color = Cycled(2)
                )
                
                if i == 1
                    gen_line = vlines!(learned_ax, [], 
                        label = "Current generation: 1",
                        color = Cycled(3)  
                    )
                else
                    gen_line = vlines!(learned_ax, [i], 
                        label = "Current generation: $(i)",
                        color = Cycled(3)  
                    )
                end

                axislegend(learned_ax,
                    position = :lt,
                    backgroundcolor =:white,
                    framevisible = true,
                    labelsize = 12,
                    patchsize = (8,8),
                    patchlabelgap = 4,
                    #rowgap = 2
                )
                rowsize!(fig.layout, 2, Aspect(2,12))

                resize_to_layout!(fig)

                if seat_config == "random"
                    save("./output/2D-Binary-PCA/random-$(class_size)-$(n_learned)/$(λ₀)/trial_$(trial)/images/2DBPCA-random-$(class_size)-$(n_learned)-$(λ₀)-$(i).png", fig, px_per_unit = 3)
                else
                    save("./output/2D-Binary-PCA/random-$(class_size)/$(λ₀)/trial_$(trial)/images/2DBPCA-random-$(class_size)-$(λ₀)-$(i).png", fig, px_per_unit = 2)
                end
            end
        end
    end
end

@time generate_plots(sizes,seat_configs, Λs, steady_state_tolerance, n_trials; n_learned = n_learned)
