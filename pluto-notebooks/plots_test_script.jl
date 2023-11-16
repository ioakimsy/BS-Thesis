using Pkg
Pkg.activate(".")
Pkg.add("BenchmarkTools")
using Plots
using BenchmarkTools

@benchmark begin
    some_var = "param"
    rand_vals = rand(20,20)

    hm = heatmap(rand_vals,
    title = "heatmap",
    xlabel = "columns $(some_var)",
    ylabel = "rows"
    )

    rand_sum = sum(rand_vals, dims = 1)
    lin_plot = plot(vec(rand_sum))

    final = plot(hm,lin_plot, layout = (2,1), dpi = 300, size=(512,1024))
end

@benchmark savefig("test_plot.png")