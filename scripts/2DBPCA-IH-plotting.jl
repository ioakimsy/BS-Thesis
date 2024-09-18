#=
! notes:
! 2024-07-12
! stopping because scaling issue. Will check if easier in plots. 

! 2024-07-12
! scaling issue fixed. Using Heatmap instead.
! Added Animations and Plotting
! Set up to do for all data
! Outer corner and inner corner for today

! 2024-07-16
! Generate plots and animations for traditional

! 2024-07-17
! Generate plots and animations for center
=#
begin
    using Pkg
    Pkg.activate(".")

    using CairoMakie
    using ColorSchemes
    using LaTeXStrings
    using CSV
    using DataFrames

    using BenchmarkTools
    using ProfileView
    using Profile

    using Alert
    using ProgressMeter

    CairoMakie.activate!(px_per_unit=4)
    update_theme!(colormap = :seaborn_colorblind)
    plots_default = ColorSchemes.seaborn_colorblind;
end

# begin
#     fig = Figure(;dpi=300)

#     class_ax = Axis(fig[1:2,1]; 
#         title = "Classroom over time",
#         titlesize = 24,
#         subtitle = L"params",
#         yreversed = true, 
#         limits = (0.5,size+0.5,0.5,size+0.5),
#         alignmode = Inside()
#         )
    
#     class = heatmap!(class_ax, _apt, 
#         colormap=(cgrad(:Blues_9, categorical = false)),
#         colorrange=(0,1),
#         fxaa = false
#     )
#     fast = scatter!(class_ax, _fast, label = "Fast students (λ=0.75)", marker = :star5, markersize = 16, color = :white, strokecolor = :black, strokewidth = 0.5)

#     slow = scatter!(class_ax, _slow, label = "Slow students(λ=0.25)", marker = :circle, markersize = 16, color = :white, strokecolor = :black, strokewidth = 0.5)

#     hidedecorations!(class_ax)

#     colsize!(fig.layout, 1, Aspect(2.0,1))


#     legend = Legend(fig[1,2], class_ax, "Student types";
#         tellheight = true,
#         tellwidth = true,
#         alignmode = Outside(),
#     )

#     cbar = Colorbar(fig[2,2], class, label = "Aptitude", vertical = true;
#     halign = :left,
#     alignmode = Outside(),
#     )

#     # Box(fig[2,2], color = (:red, 0.2), strokecolor = :red)
#     resize_to_layout!(fig)
#     fig
# end

# begin
#     fig = Figure(;dpi=300)

#     class_ax = Axis(fig[1:2,1]; 
#         title = "Classroom over time",
#         titlesize = 24,
#         subtitle = L"params",
#         yreversed = true, 
#         limits = (0.5,size+0.5,0.5,size+0.5),
#         alignmode = Inside()
#         )
    
#     fast = scatter!(class_ax, _fast, 
#         label = "Fast students", 
#         marker = :star5, 
#         markersize = 16, 
#         colormap = reverse(cgrad(:grays, categorical = false)),
#         color = _apt[CartesianIndex.(_fast)],
#         strokecolor = :black, 
#         strokewidth = 0.5,
#     )

#     fast = scatter!(class_ax, _slow, 
#         label = "Slow students", 
#         marker = :circle, 
#         markersize = 16, 
#         colormap = reverse(cgrad(:grays, categorical = false)),
#         color = _apt[CartesianIndex.(_slow)],
#         strokecolor = :black, 
#         strokewidth = 0.5,
#     )

#     hidedecorations!(class_ax)

#     colsize!(fig.layout, 1, Aspect(2.0,1))


#     legend = Legend(fig[1,2], class_ax, "Student types";
#         tellheight = true,
#         tellwidth = true,
#         alignmode = Outside(),
#     )

#     cbar = Colorbar(fig[2,2], class, label = "Aptitude", vertical = true;
#     halign = :left,
#     alignmode = Outside(),
#     )

#     # Box(fig[2,2], color = (:red, 0.2), strokecolor = :red)
#     resize_to_layout!(fig)
#     fig
# end

function generate_plots(sizes, seat_configs, Ρs, δλs, n_trials; n_learned = 4, λ₀ = 0.5)

    max_iters = max_iters = prod([length(x) for x in [sizes,seat_configs,Ρs,δλs]]) * n_trials

    prog_bar = Progress(max_iters; showspeed=true)

    for class_size in sizes, SA in seat_configs, ρ₀ in Ρs, δλ in δλs, trial in 1:n_trials

        ProgressMeter.next!(prog_bar, 
            showvalues = [("Seat config", SA), ("ρ₀", ρ₀), ("Class size", class_size), ("δλ", δλ), ("Trial", trial)]
        )

        path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
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
        learned_dom = 1:length(learned)
        fit_dom = 1:0.1:length(learned)
        power_coeff, power_exp = params_df[!,"power_fit"][1:2]
        fit_vals = power_coeff .* fit_dom .^ power_exp
        if SA == "traditional"
            learned = learned[begin+1:end]
            learned_dom = 1:length(learned)
        end

        # frames = []

        #! Plotting starts here

        for t in eachindex(generations)
            fig = Figure(;dpi=300)

            class_ax = Axis(fig[1,1]; 
                title = "Classroom over time",
                # titlesize = 24,
                # subtitlesize = 18,
                subtitle = "$(SA), L=$class_size, ρ₀=$(ρ₀), λ=$(λ₀)±$(δλ)",
                yreversed = true, 
                limits = (0.5,class_size+0.5,0.5,class_size+0.5),
                alignmode = Inside()
            )
            
            class = heatmap!(class_ax, generations[t] .* λ_mult, 
                colormap=(cgrad(:jblue, 3, categorical = true)),
                colorrange=(0,1),
                fxaa = false,
                limits = (0,1)
            )

            hidedecorations!(class_ax)

            colsize!(fig.layout, 1, Aspect(1,1))

            cbar = Colorbar(fig[1,2], class, 
                # label = "Aptitude", 
                vertical = true;
                halign = :left,
                alignmode = Outside(),
                ticks = ([1/6, 0.5, 5/6], ["Unlearned", "Learned (Slow)", "Learned (Fast)"]),
                ticklabelrotation = -π/2,
                ticksvisible = false,
                ticklabelsize = 10
                # limits = (0,1)
            )

            learned_ax = Axis(fig[2,1:2]; 
                xlabel = "Time step",
                ylabel = "Fraction of learned",
                xscale = log10,
                yscale = log10,
                xminorgridvisible = false,
                yminorgridvisible = false,
                xlabelsize = 12,
                ylabelsize = 10,
                xticklabelsize = 10,
                yticklabelsize = 10,
            )

            scatter!(learned_ax, learned_dom, learned,
                markersize = 3,
                label = "Learned students",
                color = Cycled(1),
                # color = :black,
            )

            lines!(learned_ax, fit_dom, fit_vals,
                label = L"y = %$(round(power_coeff, digits = 5))x^{%$(round(power_coeff, digits = 5))}",
                color = Cycled(2),
                linewidth = 1,
            )

            vlines!(learned_ax, [t], 
                label = "Current time step: $t", 
                color = Cycled(3),
                linewidth = 1,
            )

            vlines!(learned_ax, SA == "traditional" ? [length(generations) - Int64(floor(0.75*length(learned)))] : [length(generations) - Int64(floor(0.5*length(learned)))], 
                linestyle = :dash, 
                label = "End of fit data",
                lw = 1.5,
                color = Cycled(4),
                linewidth = 1,
            )

            Legend(fig[2,1], learned_ax,
                "Legend",
                titlesize = 7,
                titlegap = 0,
                alignmode = Mixed(left = 5, top = 5),
                halign = :left,
                valign = :top,
                patchsize = (15,9),
                labelsize = 7,
                rowgap = 0,
                tellheight = true,
                padding = (3,3,3,3),
                framevisible = true, 
                backgroundcolor = :white,
                patchlabelgap = 4,
                groupgap = 0,
            
            )
                
            rowsize!(fig.layout, 2, Relative(0.3))

            # Box(fig[2,1], color = (:red, 0.2), strokecolor = :red)
            resize_to_layout!(fig)
            img_path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/images/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
            # push!(frames, fig)
            save(img_path * "$t.png", fig)
        end
    end
end

function generate_animation(sizes, seat_configs, Ρs, δλs, n_trials; n_learned = 4, λ₀ = 0.5)

    max_iters = max_iters = prod([length(x) for x in [sizes,seat_configs,Ρs,δλs]]) * n_trials

    prog_bar = Progress(max_iters; showspeed=true)

    for class_size in sizes, SA in seat_configs, ρ₀ in Ρs, δλ in δλs, trial in 1:n_trials

        ProgressMeter.next!(prog_bar, 
            showvalues = [("Class size", class_size), ("Seat config", SA), ("ρ₀", ρ₀), ("δλ", δλ), ("Trial", trial)]
        )

        path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
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
        learned_dom = 1:length(learned)
        fit_dom = 1:0.1:length(learned)
        power_coeff, power_exp = params_df[!,"power_fit"][1:2]
        fit_vals = power_coeff .* fit_dom .^ power_exp
        if SA == "traditional"
            learned = learned[begin+1:end]
            learned_dom = 1:length(learned)
        end

        ts = Observable(1)

        fig = Figure(;dpi=300)

        class_ax = Axis(fig[1,1]; 
            title = "Classroom over time",
            # titlesize = 24,
            # subtitlesize = 18,
            subtitle = "$(SA), L=$class_size, ρ₀=$(ρ₀), λ=$(λ₀)±$(δλ)",
            yreversed = true, 
            limits = (0.5,class_size+0.5,0.5,class_size+0.5),
            alignmode = Inside()
        )
        
        class = heatmap!(class_ax, @lift(generations[$ts] .* λ_mult), 
            colormap=(cgrad(:jblue, 3, categorical = true)),
            colorrange=(0,1),
            fxaa = false,
            limits = (0,1)
        )

        hidedecorations!(class_ax)

        colsize!(fig.layout, 1, Aspect(1,1))

        cbar = Colorbar(fig[1,2], class, 
            # label = "Aptitude", 
            vertical = true;
            halign = :left,
            alignmode = Outside(),
            ticks = ([1/6, 0.5, 5/6], ["Unlearned", "Learned (Slow)", "Learned (Fast)"]),
            ticklabelrotation = -π/2,
            ticksvisible = false,
            ticklabelsize = 10
            # limits = (0,1)
        )

        learned_ax = Axis(fig[2,1:2]; 
            xlabel = @lift("Time step (t): $($ts)"),
            ylabel = "Fraction of learned (f)",
            xscale = log10,
            yscale = log10,
            xminorgridvisible = false,
            yminorgridvisible = false,
            xlabelsize = 12,
            ylabelsize = 10,
            xticklabelsize = 10,
            yticklabelsize = 10,
        )

        scatter!(learned_ax, learned_dom, learned,
            markersize = 3,
            label = "Learned students",
            color = Cycled(1),
            # color = :black,
        )

        lines!(learned_ax, fit_dom, fit_vals,
            label = L"f = %$(round(power_coeff, digits = 5))t^{%$(round(power_coeff, digits = 5))}",
            color = Cycled(2),
            linewidth = 1,
        )

        vlines!(learned_ax, @lift([$ts]), 
            label = "Current time step", 
            color = Cycled(3),
            linewidth = 1,
        )

        vlines!(learned_ax, SA == "traditional" ? [length(generations) - Int64(floor(0.75*length(learned)))] : [length(generations) - Int64(floor(0.5*length(learned)))], 
            linestyle = :dash, 
            label = "End of fit data",
            lw = 1.5,
            color = Cycled(4),
            linewidth = 1,
        )

        Legend(fig[2,1], learned_ax,
            "Legend",
            titlesize = 7,
            titlegap = 0,
            alignmode = Mixed(left = 5, top = 5),
            halign = :left,
            valign = :top,
            patchsize = (15,9),
            labelsize = 7,
            rowgap = 0,
            tellheight = true,
            padding = (3,3,3,3),
            framevisible = true, 
            backgroundcolor = :white,
            patchlabelgap = 4,
            groupgap = 0,
        
        )
            
        rowsize!(fig.layout, 2, Relative(0.3))

        # Box(fig[2,1], color = (:red, 0.2), strokecolor = :red)
        resize_to_layout!(fig)
        anim_path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
        
        record(fig, anim_path * "animation.mp4", eachindex(generations); framerate = 16) do t
            ts[] = t
        end
    end
end

begin
    sizes = [32,48,64,96,128]
	seat_configs = [
        # "outer_corner", 
        # "inner_corner", 
        "center", 
        # "random", 
        # "traditional",
        ]
	Ρs = collect(0.1:0.1:1)
    δλs = collect(0.0:0.1:0.4)
	steady_state_tolerance = 20
	n_trials = 5
    n_learned = 4
    λ₀ = 0.5

    generate_plots(sizes, seat_configs, Ρs, δλs, n_trials; n_learned = n_learned, λ₀ = λ₀)
    generate_animation(sizes, seat_configs, Ρs, δλs, n_trials; n_learned = n_learned, λ₀ = λ₀)
end
