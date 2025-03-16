begin
    println("Loading packages...")
    using Pkg

    using CSV
    using DataFrames
    using GLMakie
    using Statistics
    using Measurements
    using LsqFit
    using LaTeXStrings
    using Interpolations
    using ColorSchemes

    #Pkg.status()
    println("Done loading packages")
end

begin
    seconds = 0:0.1:2
    measurements = [8.2, 8.4, 6.3, 9.5, 9.1, 10.5, 8.6, 8.2, 10.5, 8.5, 7.2,
        8.8, 9.7, 10.8, 12.5, 11.6, 12.1, 12.1, 15.1, 14.7, 13.1]

    lines(seconds, measurements)
end

scatter(seconds, measurements)

lines(seconds, exp.(seconds) .+ 7)


scatter(seconds, measurements)
lines!(seconds, exp.(seconds) .+ 7)
current_figure()


begin
    f = Figure()

    ax = Axis(f[1, 1],
        title = "Experimental data and exponential fit",
        xlabel = "Time (seconds)",
        ylabel = "Value",
    )

    scatter!(ax, seconds, measurements)
    lines!(ax, seconds, exp.(seconds) .+ 7)

    f
end

begin
    f = Figure(
        fonts = (; regular = "Times New Roman", italic = "Times New Roman", bold = "Times New Roman"),
    )

    ax = Axis(f[1, 1],
        title = "Experimental data and exponential fit",
        # xlabel = L"\text{Time }(s)",
        xlabel = "Time (seconds)",
        ylabel = "Value",
        titlesize = 20,
        ylabelfont = "Arial",
        titlefont = "Comic Sans",
        xlabelfont = :regular,
    )
    scatter!(
        ax,
        seconds,
        measurements,
        color = 1:10,
        colormap = ColorSchemes.seaborn_colorblind,
        label = "Measurements",
    )
    lines!(
        ax,
        seconds,
        exp.(seconds) .+ 7,
        color = :blue,
        linestyle = :dash,
        label = L"f(x) = exp(x) + 7",
    )
    axislegend(ax, "Legend, wow!", position = :rb)
    f
end

ColorSchemes.viridis

ColorSchemes.seaborn_colorblind

begin
    @. sinwave(x, y, t) = sin(2π * (x * y) - 2π * t)

    @. gaussian2d(x, y, μx, μy, σx, σy) = exp(-((x - μx)^2 / (2σx^2) + (y - μy)^2 / (2σy^2)))

    x_in = -3:0.1:3
    y_in = -3:0.1:3

    normal = [gaussian2d(x, y, 0, 0, 1, 1) for x in x_in, y in y_in]

    sinusoid = [sinwave(x, y, 0) for x in x_in, y in y_in] .* normal

    fig = Figure()
    ax = Axis3(fig[1, 1],
        aspect = (1,1,1),
        title = "Gaussian but 2D"
    )

    ax2 = Axis3(fig[1, 2],
        aspect = (1,1,1),
        title = "Sin*normal but 2D"
    )

    sf = surface!(ax, x_in, y_in, normal, colormap = ColorSchemes.plasma,
        # alpha = 0.5
    )
    sinusoid_sf = surface!(ax2, x_in, y_in, sinusoid, colormap = ColorSchemes.viridis,
        # alpha = 0.5
    )


    Colorbar(fig[2,1], sf, label = "Intensity",
    vertical = false,
    )
    Colorbar(fig[1,3], sinusoid_sf, label = "Intensity",
    vertical = true,
    )
    

    Label(fig[0,2:3], "Super title at [0,2:3]",
        halign = :center,
        valign = :center,
        fontsize = 24,
    )

    Label(fig[-1,:], "Super² title na nakilimutan at [-1,:]",
        halign = :center,
        valign = :center,
        fontsize = 24,
    )

    Label(fig[end+1,:], "Subtitle na nakilimutan at [end+1,:]",
        halign = :center,
        valign = :center,
        fontsize = 24,
    )

    # Label(fig[1:2,-1], "i",
    #     halign = :center,
    #     valign = :center,
    #     fontsize = 24,
    #     # vertical = true
    # )

    fig
end



#! Other workshop things

a = collect(1:1000)
b = collect(1001:2000)

using BenchmarkTools

@benchmark [(a[i],b[i]) for i in eachindex(a)]
@benchmark collect(zip(a,b))

