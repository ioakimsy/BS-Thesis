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

    using Alert
    using ProgressMeter
end

function generate_plots(lengths::Vector{Int}, Λs::Vector{Float64}, steady_state_tolerance::Int, n_trials::Int; animation::Bool=false)
    default()
    l = @layout [
                a{0.75h}
                b{}
            ]

    println("Starting to generate plots:")
    max_iters = prod([length(x) for x in [lengths,Λs]]) * n_trials
    prog_bar = Progress(max_iters; showspeed=true)

    #* I tried Threads.@threads, it is not thread safe
    for  λ₀ in Λs, class_length in lengths, trial in 1:n_trials

        ProgressMeter.next!(prog_bar, 
                showvalues = [("λ₀", λ₀), ("Class length", class_length), ("Trial", trial)]
            )

        frames = []

        data_df = CSV.read("./output/2D-Binary-PCA/traditional-$(class_length)/$(λ₀)/trial_$(trial)/data/2DBPCA-traditional-$(class_length)-$(λ₀)-trial_$(trial)-data.csv", DataFrame)
    
        params_df = CSV.read("./output/2D-Binary-PCA/traditional-$(class_length)/$(λ₀)/trial_$(trial)/data/2DBPCA-traditional-$(class_length)-$(λ₀)-trial_$(trial)-fit_params.csv", DataFrame)


        class_length = Int(sqrt(size(data_df)[1])) #! Assumed square classroom
        num_generations = size(data_df)[2]

        generations = [data_df[!,"Generation $(i)"] for i in 1:num_generations]
        generations = [reshape(generations[i], (class_length,class_length)) for i in 1:num_generations]

        learned = params_df[!,"learned_per_gen"][begin+1:end]
        #! Check if deletable: learned_dom = 1:length(learned)
        power_coeff, power_exp = params_df[!,"power_fit"][1:2]
        fit_dom = 1:0.1:length(learned)
        fit_vals = power_coeff .* fit_dom .^ power_exp

        for i in 1:num_generations

            class_plot = heatmap(generations[i], 
                #title = "Classroom over time", #\n λ=$(λ₀) size: $(class_length)",
                titlefontvalign = :top,
                aspect_ratio=:equal,
                cbar = false,
                showaxis = false,
                c = palette(:grays,rev=true),
                #colorbar_discrete_values = [0,1],
                #cbar_title = "Aptitude",
                size = (512,512),
                yflip=true,
                dpi = 300,
                axis=([], false),
            );

            learned_plot = scatter(2:length(learned), learned, legend=:topleft,
                xlabel = "Generation number",
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
                guidefontsize = 9
                #xlims = (1,num_generations),
                #ylims = (10^(-16), maximum(learned))
            );

            learned_plot = plot!(fit_dom, fit_vals, 
                label=L"y = %$(round(power_coeff, digits=5))\cdot x^{%$(round(power_exp,digits=5))}",
                lw = 1.5
            );

            learned_plot = vline!([i], label = "Current generation: $(i)", lw = 1.5);

            learned_plot = vline!([length(generations) - Int64(floor(0.75*length(learned)))], ls = :dash, label = "End of fit data")

            _plot_title = "Classroom over time\ntraditional-$(class_length)-$(λ₀)"

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

            savefig("./output/2D-Binary-PCA/traditional-$(class_length)/$(λ₀)/trial_$(trial)/images/2DBPCA-random-$(class_length)-$(λ₀)-$(i).png")
    
        end # per generation plot end

        if animation == true
            anim = @animate for i in 1:length(frames)
                plot(frames[i]);
            end

            mp4(anim, "./output/2D-Binary-PCA/traditional-$(class_length)/$(λ₀)/trial_$(trial)/2DBPCA-random-$(class_length)-$(λ₀)-animation.mp4", fps = 16; show_msg = false)

        end #animation save end
    end # params end
end # function end

begin
    #* List of parameters
    lengths = [32,48,64,96,128]
	Λs = collect(0.1:0.1:1)
	steady_state_tolerance = 20
	n_trials = 5

    @alert generate_plots(lengths, Λs, steady_state_tolerance, n_trials; animation = true)
end
