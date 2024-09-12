begin
    using Pkg
    Pkg.activate(".")

    using Plots
    using Plots.PlotMeasures
    using LaTeXStrings
    using CSV
    using DataFrames

    using BenchmarkTools
    using ProfileView
    using Profile
    using Logging

    using Alert
    using ProgressMeter
end

function generate_plots(sizes::Vector{Int}, seat_configs::Vector{String},Ρs::Vector{Float64}, δλs::Vector{Float64}, steady_state_tolerance::Int, n_trials::Int; λ₀::Float64=0.5, n_learned::Int=4, animation::Bool=false)
    default()
    l = @layout [
                a{0.75h}
                b{}
            ]
    println("Starting to generate plots:")
    max_iters = prod([length(x) for x in [sizes,seat_configs,Ρs,δλs]]) * n_trials
    prog_bar = Progress(max_iters; showspeed=true)
    #* I tried Threads.@threads, it is not thread safe
    for class_size in sizes, SA in seat_configs, ρ₀ in Ρs, δλ in δλs ,trial in 1:n_trials

        ProgressMeter.next!(prog_bar, 
                showvalues = [("Class size", class_size), ("Seat config", SA), ("ρ₀", ρ₀), ("δλ", δλ), ("Trial", trial)]
            )

        frames = []

        path = SA == "random" ? "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)-$(n_learned)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-" : "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"

        data_df = CSV.read(path * "data.csv", DataFrame)
        params_df = CSV.read(path * "fit_params.csv", DataFrame)

        λ_grid = Matrix(CSV.read(path * "lambda_grid.csv", DataFrame))
        fast = Tuple.(findall(x -> x == λ₀ + δλ, λ_grid))
        slow = Tuple.(findall(x -> x == λ₀ - δλ, λ_grid))

        λ_mult = λ_grid
        λ_mult[CartesianIndex.(fast)] .= 1
        λ_mult[CartesianIndex.(slow)] .= 0.5

        class_size = Int(sqrt(size(data_df)[1])) #! Assumed square classroom
        num_generations = size(data_df)[2]

        generations = [data_df[!,"Generation $(i)"] for i in 1:num_generations]
        generations = [reshape(generations[i], (class_size,class_size)) for i in 1:num_generations]

        learned = params_df[!,"learned_per_gen"]
        #! Check if deletable: learned_dom = 1:length(learned)
        power_coeff, power_exp = params_df[!,"power_fit"][1:2]
        fit_dom = 1:0.1:length(learned)
        fit_vals = power_coeff .* fit_dom .^ power_exp

        for i in 1:num_generations
            _generation_plot = generations[i] .* λ_mult
            class_plot = heatmap(_generation_plot, 
                #title = "Classroom over time", #\n λ=$(ρ₀) size: $(class_size)",
                titlefontvalign = :top,
                aspect_ratio = :equal,
                cbar = true,
                showaxis = false,
                c = cgrad(:Blues_9, 3, rev = false, categorical = true),
                clims = (0,1),
                colorbar_ticks = (false),
                #cbar_title = "Aptitude",
                size = (512,512),
                yflip = true,
                dpi = 300,
                axis=([], false);
                show_msg = false
            );

            learned_plot = scatter(learned, legend=:topleft,
                xlabel = "Time step",
                ylabel = "Fraction of learned",
                #title = "Learned over time",
                label = "Learned students",
                legend_font_pointsizes = 5,
                yrot = 0,
                ytickfontsize = 4,
                dpi = 300,
                leftmargin = 5mm,
                rightmargin = 5mm,
                markersize = 3,
                scale = :log10,
                guidefontsize = 9;
                show_msg = false
            );

            learned_plot = plot!(fit_dom, fit_vals, 
                label=L"y = %$(round(power_coeff, digits=5))\cdot x^{%$(round(power_exp,digits=5))}",
                lw = 1.5
            );

            learned_plot = vline!([i], label = "Current time step: $(i)", lw = 1.5);

            learned_plot = vline!([length(generations) - Int64(floor(0.5*length(learned)))], ls = :dash, label = "End of fit data")

            _plot_title = SA == "random" ? "Classroom over time\n$(SA)-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)±$(δλ)" : "Classroom over time\n$(SA)-$(class_size)-$(ρ₀)-$(λ₀)±$(δλ)"

            _plot = plot(class_plot,learned_plot, 
                layout = l,
                size = (512, 512/0.9),
                dpi = 300,
                plot_title = _plot_title,
                plot_titlefontsize = 11
            );

            if animation == true
                push!(frames,_plot)
            end

            if SA == "random"
                savefig("./output/2D-Binary-PCA-IH/random-$(class_size)-$(n_learned)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/images/2DBPCAIH-random-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)-$(δλ)-$(i).png")
            else
                savefig("./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/images/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-$(i).png")
            end # frame save end

        end # per generation plot end

        if animation == true
            anim = @animate for i in 1:length(frames)
                plot(frames[i]);
            end

            if SA == "random"
                mp4(anim, "./output/2D-Binary-PCA-IH/random-$(class_size)-$(n_learned)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/2DBPCAIH-random-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)-$(δλ)-animation.mp4", fps = 16; show_msg=false)
            else
                mp4(anim, "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-animation.mp4", fps = 16; show_msg = false)
            end # animation save end
        end
    end # params end
end # function end

begin
    # sizes = [32,48,64,96,128]
	# seat_configs = ["outer_corner", "inner_corner", "center", "random", "traditional"]
	# Ρs = collect(0.1:0.1:1)
    # δλs = collect(0.0:0.1:0.4)

    sizes = [32]
	seat_configs = ["inner_corner"]
	Ρs = [0.3]
    δλs = [0.2]

	steady_state_tolerance = 20
	n_trials = 1
    n_learned = 4
    λ₀ = 0.5

    with_logger(ConsoleLogger(Error)) do
        generate_plots(sizes, seat_configs, Ρs, δλs, steady_state_tolerance, n_trials; λ₀ = λ₀, n_learned = n_learned, animation = true)
    end
end

_size = 16
_lr = Float64.(bitrand(_size,_size))
_lr[_lr .== 1] .= 0.7
_lr[_lr .== 0] .= 0.3


_apt = Int.(bitrand(_size,_size))
_apt_plot = _apt .* _lr

_fast = (findall(x -> x == 0.7, _lr))
_slow = (findall(x -> x == 0.3, _lr))

_lr[_fast] .= 2
_lr[_slow] .= 1

heatmap(_apt_plot, 
                #title = "Classroom over time", #\n λ=$(ρ₀) size: $(class_size)",
                titlefontvalign = :top,
                aspect_ratio = :equal,
                cbar = true,
                showaxis = false,
                c = palette(:grays,rev=true),
                #colorbar_discrete_values = [0,1],
                #cbar_title = "Aptitude",
                size = (512,512),
                yflip = true,
                dpi = 300,
                axis=([], false);
                show_msg = false
            )