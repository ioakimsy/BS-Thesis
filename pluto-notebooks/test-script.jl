using Pkg
Pkg.activate(".")
#Pkg.add("CairoMakie")
using CairoMakie
using BenchmarkTools
begin
    some_var = "param"
    rand_vals = rand(20,20)
    fig, ax, hm = heatmap(rand_vals;
        axis = (; title="heatmap",
            ylabel = "rows",
            xlabel = "columns $(some_var)"
        ),
        fig = (; dpi = 300),
    )
    colsize!(fig.layout, 1, Aspect(1,1))
    Colorbar(fig[1,2], hm)

    rand_sum  = sum(rand_vals, dims=1)
    ax1 = Axis(fig[2,1]; title = "help", xlabel = "x", ylabel = "y")
    lin_plot = lines!(ax1, vec(rand_sum))

    ax1.title = "a new axis axis title"
    #=
    sax = Axis(fig[2,3]; 
        title = "scatter",
        ylabel = "scatter y",
        xlabel = "scatter x",
        aspect = 1
        )
    
    colsize!(fig.layout, 3, Aspect(1,1))

    scatter_plot = scatter!(sax, rand(20), rand(20))

    ax2, bar = hist(fig[1,3], rand(100); axis =  (;title = "histogram"))
    =#
    resize_to_layout!(fig)
    fig
end

#save("./test_plot.png", fig)