# Run from the repository root (paths are relative to it); see README.md.
#=
2DBPCA-IH v2: paper figures regenerated from the v2 sweep.
  * Fig. 3 — learned-fraction time series (figure3_learned_t below)
  * Fig. 5 — t_max surfaces over (ρ₀, δλ) for TI vs PI inner corner, L = 32 and 128 (figure5);
    the same with t₉₅, and t₉₅ vs t_max line plots incl. δλ = 0.45/0.49 (figure_t95_vs_tmax)

Same four panels and visual encoding as the v1 figure (2DBPCA-IH-trad-PI-comparison.jl):
  (a) lattice side L      — TI vs PI inner corner, ρ₀=0.5, δλ=0.2, L ∈ {32, 64, 128}
  (b) seating arrangement — L=64, ρ₀=0.5, δλ=0.2
  (c) coupling ρ₀         — inner corner, L=64, δλ=0.2, ρ₀ ∈ {0.1, 0.5, 0.9}
  (d) heterogeneity δλ    — inner corner, L=64, ρ₀=0.5, δλ ∈ {0, 0.2, 0.4, 0.45, 0.49}
Color = regime/SA, line width = L, opacity = ρ₀, line style = δλ. Lines are the mean over
the 20 trials; bands are ±1 SD. Runs are padded with f = 1 after absorption.

The x axis is iteration t+1 on a log₂ scale, so 2⁰ is the initial state (as in v1).
TI now starts from the same 4 seeds as PI, so its curves start at 4/N rather than 0.

=#

begin
    using Pkg
    Pkg.activate(".")
    include("2DBPCA-IH-common-v2.jl")
    using CairoMakie, ColorSchemes, Statistics
    import GLMakie #* only for the 3D surfaces (needs real depth sorting)
    const FileIO = CairoMakie.Makie.FileIO #* stacking the two Fig. 5 rows (FileIO via Makie)
    CairoMakie.activate!() #* importing GLMakie makes it the active backend; keep CairoMakie as default
end

const PLOT_DIR = "$(OUTPUT_DIR)/analysis/plots/learned-t"

const COLORS = ColorSchemes.seaborn_colorblind.colors
const SA_COLOR = Dict("TI" => 1, "inner_corner" => 2, "outer_corner" => 3, "center" => 4, "random" => 5)
const SA_LABEL = Dict("TI" => "Traditional", "inner_corner" => "Inner corner", "outer_corner" => "Outer corner",
    "center" => "Center", "random" => "Random")
const L_WIDTH = Dict(32 => 1.5, 48 => 3.0, 64 => 4.5, 96 => 6.0, 128 => 7.5)
const δλ_STYLE = Dict(0.0 => :dash, 0.2 => :dashdot, 0.4 => (:dashdotdot, :dense), 0.45 => :dashdotdot,
    0.49 => (:dashdotdot, :loose))
ρ_alpha(ρ₀) = 0.5 + 0.5 * (ρ₀ - 0.1) / 0.9

#* Series color: one color for TI, one per seating arrangement for PI
regime_color(regime, SA = "inner_corner") = COLORS[SA_COLOR[regime == "TI" ? "TI" : SA]]

#* Paper copy of a figure: legends removed (and the empty legend row/column trimmed). The legend
#* version is kept for reference; series are identified in the caption or by annotation.
strip_legends!(fig) = (foreach(delete!, filter(c -> c isa Legend, copy(fig.content))); trim!(fig.layout); fig)

#* Mean and SD of f(t) over trials, t = 0..t_end (runs padded with their final value).
function mean_curve(regime, SA, L, ρ₀, δλ; n_learned = 4, λ₀ = 0.5)
    Ts = [T for (T, _) in read_trials(param_row(regime, SA, L, ρ₀, δλ; n_learned, λ₀))]
    t_end = maximum(maximum.(Ts))
    F = reduce(hcat, [learned_fraction(T, t_end) for T in Ts])
    return vec(mean(F; dims = 2)), vec(std(F; dims = 2))
end

#* Mean t_max or t₉₅ over the 20 runs, per (regime, L, ρ₀, δλ); inner corner, 4 initial seeds
function mean_times(stat)
    runs = CSV.read(SUMMARY_PATH, DataFrame)
    runs = runs[(runs.n_learned .== 4) .& (runs.SA .== "inner_corner"), :]
    return combine(groupby(runs, [:regime, :L, :ρ₀, :δλ]), stat => mean => :t)
end

#* Matrix of mean times over (ρ₀, δλ) for one regime and lattice side, for the surface plots
surface_grid(g, regime, L, Ρs, δλs) = [only(rows_where(g; regime, L, ρ₀, δλ).t) for ρ₀ in Ρs, δλ in δλs]

#* Legend entries for the surface plots
regime_legend_entries() = ([PolyElement(color = COLORS[1]), PolyElement(color = COLORS[2])],
    [rich("Traditional"; color = COLORS[1], font = :bold), rich("Peer"; color = COLORS[2], font = :bold)])

function draw_series!(ax, μ, σ; regime, SA, L, ρ₀, δλ, label, solid = false)
    c = regime_color(regime, SA)
    x = 1:length(μ)
    band!(ax, x, clamp.(μ .- σ, 0, 1), clamp.(μ .+ σ, 0, 1); color = (c, ρ_alpha(ρ₀) / 4))
    lines!(ax, x, μ; color = (c, ρ_alpha(ρ₀)), linewidth = L_WIDTH[L],
        linestyle = solid ? :solid : δλ_STYLE[δλ], label)
end

function panel_axis(pos; xlabel = "", ylabel = "", letter)
    ax = Axis(pos; xscale = log2, xlabel, ylabel, xlabelsize = 28, ylabelsize = 28, xticklabelsize = 22,
        yticklabelsize = 22, yticks = 0:0.5:1, xticklabelsvisible = !isempty(xlabel))
    text!(ax, 0.02, 0.98; text = letter, align = (:left, :top), space = :relative, fontsize = 32)
    return ax
end

function add_legend!(ax; position = :rb)
    axislegend(ax; position, labelsize = 14, framevisible = false, backgroundcolor = :transparent,
        patchsize = (40, 12), rowgap = 0)
end

#* In-plot labels as in the published Fig. 3 (bold, colored like their series), placed in
#* data coordinates: (text, x = iteration on the log₂ axis, y, series color, ρ₀ for opacity).
label!(ax, txt, x, y, key; ρ₀ = 0.5, rotation = 0.0, align = (:left, :center), fontsize = 22) =
    text!(ax, x, y; text = txt, color = (COLORS[SA_COLOR[key]], ρ_alpha(ρ₀)), font = :bold,
        fontsize, rotation, align)

#* The four Fig. 3 panels. Series are drawn in the given order (later series on top).
#* Mean ± SD learned-fraction curve of one series, drawn with the encoding of draw_series!
function draw_mean_curve!(ax, regime, L, ρ₀, δλ; SA = "inner_corner")
    μ, σ = mean_curve(regime, SA, L, ρ₀, δλ)
    draw_series!(ax, μ, σ; regime, SA, L, ρ₀, δλ, label = "")
end

#* (a) lattice side length; ρ₀ = 0.5, δλ = 0.2
function draw_size_panel!(ax)
    for L in (32, 64, 128), regime in ("TI", "PI")
        draw_mean_curve!(ax, regime, L, 0.5, 0.2)
    end
end

#* (b) seating arrangement; L = 64, ρ₀ = 0.5, δλ = 0.2
function draw_seating_panel!(ax, SAs)
    draw_mean_curve!(ax, "TI", 64, 0.5, 0.2)
    for SA in SAs
        draw_mean_curve!(ax, "PI", 64, 0.5, 0.2; SA)
    end
end

#* (c) coupling strength; L = 64, δλ = 0.2
function draw_coupling_panel!(ax, ρs)
    for ρ₀ in ρs, regime in ("TI", "PI")
        draw_mean_curve!(ax, regime, 64, ρ₀, 0.2)
    end
end

#* (d) heterogeneity; L = 64, ρ₀ = 0.5
function draw_heterogeneity_panel!(ax)
    for δλ in (0.0, 0.2, 0.4, 0.45, 0.49), regime in ("TI", "PI")
        draw_mean_curve!(ax, regime, 64, 0.5, δλ)
    end
end

function figure3()
    fig = Figure(size = (1650, 1250))
    ax = panel_axis(fig[1, 1]; letter = "(a)")
    draw_size_panel!(ax)
    label!(ax, "L = 32, 64, 128", 2^0.5, 0.66, "TI")
    label!(ax, "L = 32", 2^2.2, 0.48, "inner_corner")
    label!(ax, "L = 64", 2^4.3, 0.55, "inner_corner")
    label!(ax, "L = 128", 2^5.55, 0.80, "inner_corner")

    ax = panel_axis(fig[2, 1]; xlabel = "Iteration number t + 1", letter = "(b)")
    draw_seating_panel!(ax, ("inner_corner", "center", "outer_corner", "random"))
    label!(ax, "Traditional", 2^1.0, 0.66, "TI")
    label!(ax, "Inner Corner", 2^2.1, 0.29, "inner_corner")
    label!(ax, "Random", 2^5.6, 0.77, "random")
    label!(ax, "Center", 2^5.3, 0.46, "center")
    label!(ax, "Outer Corner", 2^4.6, 0.14, "outer_corner")

    ax = panel_axis(fig[1, 2]; letter = "(c)")
    draw_coupling_panel!(ax, (0.9, 0.5, 0.1))
    label!(ax, "ρ = 0.9", 2^0.05, 0.68, "TI"; ρ₀ = 0.9)
    label!(ax, "ρ = 0.5", 2^2.05, 0.47, "TI"; ρ₀ = 0.5)
    label!(ax, "ρ = 0.1", 2^1.2, 0.21, "TI"; ρ₀ = 0.1)
    label!(ax, "ρ = 0.9", 2^2.6, 0.69, "inner_corner"; ρ₀ = 0.9)
    label!(ax, "ρ = 0.5", 2^4.55, 0.44, "inner_corner"; ρ₀ = 0.5)
    label!(ax, "ρ = 0.1", 2^6.25, 0.66, "inner_corner"; ρ₀ = 0.1)

    ax = panel_axis(fig[2, 2]; xlabel = "Iteration number t + 1", letter = "(d)")
    draw_heterogeneity_panel!(ax)
    label!(ax, "δλ = 0.0, 0.2", 2^0.55, 0.90, "TI"; fontsize = 20)
    label!(ax, "δλ = 0.4", 2^2.85, 0.645, "TI"; rotation = 0.6, fontsize = 18)
    label!(ax, "δλ = 0.45", 2^3.2, 0.575, "TI"; rotation = 0.45, fontsize = 18)
    label!(ax, "δλ = 0.49", 2^8.6, 0.78, "TI")
    label!(ax, "δλ = 0.0, 0.2\n0.4, 0.45", 2^1.25, 0.18, "inner_corner")
    label!(ax, "δλ = 0.49", 2^5.9, 0.52, "inner_corner")
    Label(fig[:, 0], "Fraction of learned f(t)"; fontsize = 32, rotation = π / 2, tellheight = false)
    return fig
end

function figure5(; δλs = [0.0, 0.1, 0.2, 0.3, 0.4], Ρs = round.(collect(0.1:0.1:1.0); digits = 1), stat = :t_max,
        zlabel = "Time to learn tₘₐₓ", letters = ("(a)", "(b)"), legend = true, what = "")
    GLMakie.activate!(; scalefactor = 1) #* ignore OS display scaling so the layout matches the requested size
    mean_t = mean_times(stat)

    fig = Figure(size = (1900, 950), figure_padding = (20, 20, 20, 20))
    for (j, (L, letter)) in enumerate(zip((32, 128), letters))
        ax = Axis3(fig[1, j]; xlabel = "Transmission ρ₀", ylabel = "Heterogeneity δλ", zlabel,
            aspect = (1, 1, 1), azimuth = 1.275π, viewmode = :fit, xlabelsize = 30, ylabelsize = 30, zlabelsize = 30,
            xticklabelsize = 22, yticklabelsize = 22, zticklabelsize = 22, xlabeloffset = 55, ylabeloffset = 55,
            zlabeloffset = 70)
        for regime in ("TI", "PI")
            surface!(ax, Ρs, δλs, surface_grid(mean_t, regime, L, Ρs, δλs);
                color = fill(regime_color(regime), length(Ρs), length(δλs)), shading = MultiLightShading)
        end
        Label(fig[1, j, Top()], "$letter L = $L" * (isempty(what) ? "" : ", $what"); fontsize = 40, halign = :left,
            padding = (40, 0, 0, 0))
    end
    legend && Legend(fig[1, 1], regime_legend_entries()...;
        tellwidth = false, tellheight = false, halign = :right, valign = :top, margin = (0, 40, 0, 60),
        labelsize = 34, patchsize = (55, 55), rowgap = 30, framevisible = false)
    return fig
end

#* Fig. 5 extended (revision): t_max row above a t₉₅ row. A single 2×2 GLMakie figure is taller than the
#* screen and renders squashed, so the two rows are rendered separately and stacked.
#* With legend = true, a horizontal legend strip is stacked above the panels (inside panel (a) it overlapped
#* the axes); with legend = false, the result is the clean copy for labeling in PowerPoint.
function save_figure5_with_t95(path; legend = true)
    top = figure5(; what = "last agent learned", legend = false)
    bottom = figure5(; stat = :t_95, zlabel = "Time to 95% learned t₉₅", letters = ("(c)", "(d)"), legend = false,
        what = "95% of agents learned")
    figs = Any[top, bottom]
    if legend
        strip = Figure(size = (1900, 90))
        Legend(strip[1, 1], regime_legend_entries()...;
            orientation = :horizontal, labelsize = 34, patchsize = (55, 55), colgap = 60, framevisible = false)
        pushfirst!(figs, strip)
    end
    #* render each part to a temporary PNG, then stack the images vertically
    tmp = [tempname() * ".png" for _ in figs]
    for (fig, file) in zip(figs, tmp)
        GLMakie.save(file, fig; px_per_unit = 2)
    end
    FileIO.save(path, reduce(vcat, FileIO.load.(tmp)))
    foreach(rm, tmp)
end

#* ---------------------------------------------------------------------------
#* Clean copies for re-labeling in PowerPoint: same canvas, fonts, and pixel size as the v1
#* source images (size-SA_ρ₀-δλ_quad.png, 8000×6000; {32,128}_truncated.png, 2500×2000),
#* with panel letters and axis labels but no series labels or legends.
#* ---------------------------------------------------------------------------

#* Fig. 3 as the v1 quad: 2000×1500 at 4 px/unit, paired-plot fonts, x axes linked per column
#* (so (c) shares (d)'s range, as in the published figure).
function figure3_clean()
    fig = Figure(size = (2000, 1500))
    ax_kw = (; xscale = log2, xlabelsize = 40, ylabelsize = 40, xticklabelsize = 32, yticklabelsize = 32,
        yticks = 0:0.5:1)
    letter!(ax, s) = text!(ax, 0.02, 0.98; text = s, align = (:left, :top), space = :relative, fontsize = 40)
    axa = Axis(fig[1, 1]; ax_kw..., xticksvisible = false, xticklabelsvisible = false)
    axb = Axis(fig[2, 1]; ax_kw..., xlabel = "Iteration number t + 1")
    axc = Axis(fig[1, 2]; ax_kw..., xticksvisible = false, xticklabelsvisible = false)
    axd = Axis(fig[2, 2]; ax_kw..., xlabel = "Iteration number t + 1")
    draw_size_panel!(axa)
    draw_seating_panel!(axb, ("inner_corner", "outer_corner", "center", "random"))
    draw_coupling_panel!(axc, (0.1, 0.5, 0.9))
    draw_heterogeneity_panel!(axd)
    for (ax, letter) in zip((axa, axb, axc, axd), ("(a)", "(b)", "(c)", "(d)"))
        letter!(ax, letter)
    end
    linkxaxes!(axa, axb)
    linkxaxes!(axc, axd)
    Label(fig[:, 0], "Fraction of learned f(t)"; fontsize = 40, rotation = π / 2, tellheight = false)
    return fig
end

#* One Fig. 5 panel as v1's {L}_truncated.png: 1000×800 at 2.5 px/unit, same Axis3 settings,
#* no title or legend. GLMakie (depth sorting); call with GLMakie active.
function figure5_panel(L; δλs = [0.0, 0.1, 0.2, 0.3, 0.4], Ρs = round.(collect(0.1:0.1:1.0); digits = 1), stat = :t_max,
        zlabel = "Time to learn tₘₐₓ")
    mean_t = mean_times(stat)
    grids = Dict(regime => surface_grid(mean_t, regime, L, Ρs, δλs) for regime in ("TI", "PI"))
    t_top = maximum(maximum, values(grids))
    fig = Figure(size = (1000, 800))
    ax = Axis3(fig[1, 1]; xlabel = "Transmission ρ₀", ylabel = "Heterogeneity δλ", zlabel,
        aspect = (1, 1, 1), zticks = 0:(t_top > 1000 ? 200 : 100):t_top, xlabelsize = 32, ylabelsize = 32,
        zlabelsize = 32, xlabeloffset = 60, ylabeloffset = 60, zlabeloffset = 80, xticklabelsize = 28,
        yticklabelsize = 28, zticklabelsize = 28, viewmode = :fit, azimuth = 1.275π)
    for regime in ("TI", "PI")
        surface!(ax, Ρs, δλs, grids[regime]; color = fill(regime_color(regime), length(Ρs), length(δλs)),
            shading = FastShading)
    end
    return fig
end

function save_clean_figures()
    dir = "$(OUTPUT_DIR)/analysis/plots/clean"
    mkpath(dir)
    CairoMakie.activate!()
    save("$(dir)/figure3-quad-clean.png", figure3_clean(); px_per_unit = 4)
    GLMakie.activate!(; scalefactor = 1)
    for L in (32, 128)
        GLMakie.save("$(dir)/figure5-$(L)-clean.png", figure5_panel(L); px_per_unit = 2.5)
        GLMakie.save("$(dir)/figure5-t95-$(L)-clean.png",
            figure5_panel(L; stat = :t_95, zlabel = "Time to 95% learned t₉₅"); px_per_unit = 2.5)
    end
    CairoMakie.activate!()
    #* exploration / supplementary figures: legend-free copies (same size as the legend versions)
    save("$(dir)/t95-vs-tmax-clean.png", strip_legends!(figure_t95_vs_tmax(; title = false)); px_per_unit = 2)
    save("$(dir)/seed-density-comparison-clean.png", strip_legends!(figure_seed_density()); px_per_unit = 2)
    save("$(dir)/class-split-learning-curves-clean.png", strip_legends!(figure_class_split(; title = false)); px_per_unit = 2)
    println("Clean figures saved to $(dir)")
end

#* ---------------------------------------------------------------------------
#* Seed density (reviewer R1): 4 initial seeds vs a fixed initial learned fraction 4/32²
#* ---------------------------------------------------------------------------
const Ls = [32, 48, 64, 96, 128]
n_frac(L) = round(Int, 4 * L^2 / 32^2)
const SA_ORDER = ("inner_corner", "center", "outer_corner", "random")

#* log-log slope of mean t_max vs L
loglog_slope(t) = (x = log.(Ls); y = log.(t); cov(x, y) / var(x))

function figure_seed_density()
    runs = CSV.read(SUMMARY_PATH, DataFrame)
    tmax_stats = combine(groupby(runs, [:regime, :SA, :L, :n_learned, :ρ₀, :δλ]), :t_max => mean => :t, :t_max => std => :s)
    #* (mean, SD) of t_max at one lattice side
    function stat(regime, SA, L, n_learned, ρ₀, δλ)
        match = rows_where(tmax_stats; regime, SA, L, n_learned, ρ₀, δλ)
        return match.t[1], match.s[1]
    end
    #* (mean, SD) for every L, with 4 initial seeds or with the fixed initial fraction
    series(regime, SA, ρ₀, δλ, frac) = [stat(regime, SA, L, frac ? n_frac(L) : 4, ρ₀, δλ) for L in Ls]

    fig = Figure(size = (1600, 1250))
    for (j, (ρ, δ)) in enumerate(((0.5, 0.2), (0.1, 0.4)))
        ax = Axis(fig[1, j]; xscale = log2, yscale = log10, xticks = (Ls, string.(Ls)), xlabel = "Lattice side length L",
            ylabel = j == 1 ? "Time to learn tₘₐₓ" : "", title = "($('a' + j - 1)) PI tₘₐₓ vs L, ρ₀ = $ρ, δλ = $δ",
            titlesize = 22, xlabelsize = 22, ylabelsize = 22, xticklabelsize = 18, yticklabelsize = 18)
        for frac in (false, true)
            ti_tmax = first.(series("TI", "inner_corner", ρ, δ, frac))
            lines!(ax, Ls, ti_tmax; color = (:gray40, 0.8), linewidth = 2, linestyle = frac ? :dash : :solid,
                label = frac ? "TI, fixed fraction" : "TI, 4 seeds")
        end
        for SA in SA_ORDER, frac in (false, true)
            pi_stats = series("PI", SA, ρ, δ, frac)
            tmax_mean, tmax_sd = first.(pi_stats), last.(pi_stats)
            col = COLORS[SA_COLOR[SA]]
            band!(ax, Ls, max.(tmax_mean .- tmax_sd, 1), tmax_mean .+ tmax_sd; color = (col, 0.12))
            scatterlines!(ax, Ls, tmax_mean; color = col, linewidth = 3, linestyle = frac ? :dash : :solid, markersize = 10,
                marker = frac ? :utriangle : :circle,
                label = "$(SA_LABEL[SA]), $(frac ? "fixed fraction" : "4 seeds") (γ = $(round(loglog_slope(tmax_mean); digits = 2)))")
        end
        axislegend(ax; position = :lt, labelsize = 13, framevisible = false, backgroundcolor = :transparent, rowgap = 0)
    end

    #* (c) learning curves at L = 128
    ax = Axis(fig[2, 1]; xscale = log2, xlabel = "Iteration number t + 1", ylabel = "Fraction of learned f(t)", yticks = 0:0.5:1,
        title = "(c) PI learning curves, L = 128, ρ₀ = 0.5, δλ = 0.2", titlesize = 22, xlabelsize = 22, ylabelsize = 22,
        xticklabelsize = 18, yticklabelsize = 18)
    for SA in SA_ORDER, frac in (false, true)
        μ, σ = mean_curve("PI", SA, 128, 0.5, 0.2; n_learned = frac ? n_frac(128) : 4)
        col = COLORS[SA_COLOR[SA]]
        band!(ax, 1:length(μ), clamp.(μ .- σ, 0, 1), clamp.(μ .+ σ, 0, 1); color = (col, 0.12))
        lines!(ax, 1:length(μ), μ; color = col, linewidth = 3, linestyle = frac ? :dash : :solid,
            label = "$(SA_LABEL[SA]), $(frac ? "$(n_frac(128)) seeds" : "4 seeds")")
    end
    axislegend(ax; position = :lt, labelsize = 13, framevisible = false, backgroundcolor = :transparent, rowgap = 0)

    #* (d) slope with fixed fraction vs slope with 4 seeds, every (ρ₀, δλ)
    ax = Axis(fig[2, 2]; xlabel = "Slope γ, 4 seeds", ylabel = "Slope γ, fixed seed fraction",
        title = "(d) Size scaling, all ρ₀ × δλ (PI; TI in gray)", titlesize = 22, xlabelsize = 22, ylabelsize = 22,
        xticklabelsize = 18, yticklabelsize = 18, aspect = DataAspect())
    ρs = round.(collect(0.1:0.1:1.0); digits = 1); δs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49]
    slopes(regime, SA) = [(loglog_slope(first.(series(regime, SA, ρ, δ, false))), loglog_slope(first.(series(regime, SA, ρ, δ, true))))
        for ρ in ρs, δ in δs] |> vec
    ti_slopes = slopes("TI", "inner_corner")
    scatter!(ax, first.(ti_slopes), last.(ti_slopes); color = (:gray40, 0.6), markersize = 8, label = "TI")
    for SA in SA_ORDER
        s = slopes("PI", SA)
        scatter!(ax, first.(s), last.(s); color = (COLORS[SA_COLOR[SA]], 0.75), markersize = 10, label = SA_LABEL[SA])
    end
    lines!(ax, [-0.1, 1.2], [-0.1, 1.2]; color = :black, linestyle = :dot, label = "y = x")
    limits!(ax, -0.1, 1.2, -0.1, 1.2)
    axislegend(ax; position = :lt, labelsize = 14, framevisible = false, backgroundcolor = :transparent, rowgap = 0)
    return fig
end

#* ---------------------------------------------------------------------------
#* Lattice snapshots
#* ---------------------------------------------------------------------------
#* Paper Fig. 5: lattice snapshots of one TI/PI twin pair (same seeds and λ grid) at fixed times, in the
#* v1 layout and colors: unlearned, learned slow, learned fast. labels = false gives the clean copy.
const FIG4_COLORS = [colorant"#d6e6f9", colorant"#2a87f5", colorant"#0a1a33"]

#* Snapshot category of an agent at step t (index into FIG4_COLORS)
function snapshot_state(learning_time, λ_class, t)
    0 <= learning_time <= t || return 1    #* unlearned
    return λ_class == 1 ? 3 : 2            #* learned fast / learned slow
end

function figure4_snapshots(; L = 32, ρ₀ = 0.3, δλ = 0.2, ts = (2, 8, 16, 32, 54), trial = 1, labels = true)
    fig = Figure(size = (1775, 850), figure_padding = (10, 10, 10, 10))
    for (i, regime) in enumerate(("TI", "PI"))
        T, cls = read_trial(param_row(regime, "inner_corner", L, ρ₀, δλ), trial)
        for (j, t) in enumerate(ts)
            S = snapshot_state.(T, cls, t)
            ax = Axis(fig[i, j]; aspect = DataAspect(), yreversed = true, title = (labels && i == 1) ? "t = $t" : "",
                titlesize = 40, titlefont = :regular)
            hidedecorations!(ax)
            heatmap!(ax, permutedims(S); colormap = FIG4_COLORS, colorrange = (1, 3))
        end
        labels && Label(fig[i, 0], regime == "TI" ? "Traditional" : "Peer"; rotation = π / 2, fontsize = 40, tellheight = false)
    end
    if labels
        Legend(fig[3, 1:5], [PolyElement(color = c, strokecolor = :black, strokewidth = 1) for c in FIG4_COLORS],
            ["Unlearned", "Learned, slow", "Learned, fast"]; orientation = :horizontal, labelsize = 32, patchsize = (40, 40),
            colgap = 50, framevisible = false)
    end
    colgap!(fig.layout, 8); rowgap!(fig.layout, 8)
    return fig
end

#* ---------------------------------------------------------------------------
#* Fast vs slow learners (paper Fig. 4): mean learned fraction per λ class, TI vs PI
#* ---------------------------------------------------------------------------
#* Mean over trials of the per-class learned fraction (runs padded with their final value).
function mean_class_curve(regime, SA, L, ρ₀, δλ, class; n_learned = 4, λ₀ = 0.5)
    runs = read_trials(param_row(regime, SA, L, ρ₀, δλ; n_learned, λ₀))
    t_end = maximum(maximum(T) for (T, _) in runs)
    F = reduce(hcat, [learned_fraction_by_class(T, c, class, t_end) for (T, c) in runs])
    return vec(mean(F; dims = 2))
end

function figure_class_split(; SA = "inner_corner", L = 128, Ρs = (0.5, 0.1), δλs = (0.2, 0.4, 0.49), title = true)
    fig = Figure(size = (1800, 1000))
    for (i, ρ₀) in enumerate(Ρs), (j, δλ) in enumerate(δλs)
        ax = Axis(fig[i, j]; xscale = log2, yticks = 0:0.5:1, title = "ρ₀ = $ρ₀, δλ = $δλ", titlesize = 22,
            xlabel = i == length(Ρs) ? "Iteration number t + 1" : "", ylabel = j == 1 ? "Fraction of learned f(t)" : "",
            xlabelsize = 20, ylabelsize = 20, xticklabelsize = 16, yticklabelsize = 16)
        for regime in ("TI", "PI")
            col = regime_color(regime, SA)
            μ, _ = mean_curve(regime, SA, L, ρ₀, δλ)
            f₊ = mean_class_curve(regime, SA, L, ρ₀, δλ, 1)
            f₋ = mean_class_curve(regime, SA, L, ρ₀, δλ, 2)
            lines!(ax, 1:length(μ), μ; color = (col, 0.45), linewidth = 1.5, label = "$regime whole class")
            lines!(ax, 1:length(f₊), f₊; color = col, linewidth = 3, label = "$regime fast λ = $(0.5 + δλ)")
            lines!(ax, 1:length(f₋), f₋; color = col, linewidth = 3, linestyle = :dash, label = "$regime slow λ = $(round(0.5 - δλ; digits = 2))")
        end
        axislegend(ax; position = :lt, framevisible = false, labelsize = 13, rowgap = 0)
    end
    title && Label(fig[0, :], "Fast vs slow learners: mean of $(N_TRIALS) runs, PI = $(SA_LABEL[SA]), L = $L"; fontsize = 24)
    return fig
end

#* ---------------------------------------------------------------------------
#* t₉₅ alongside t_max (mechanism finding 3). t_max is the time of the last learner, so at high δλ it
#* is a straggler statistic (the slowest of ~N/2 slow learners); t₉₅ is the time to 95% learned.
#* ---------------------------------------------------------------------------

#* Mean ± SD of t_max (solid) and t₉₅ (dashed) vs δλ, TI vs PI inner corner, 4 seeds, all δλ incl. 0.45/0.49.
function figure_t95_vs_tmax(; Ls = (32, 128), Ρs = (0.1, 0.5, 0.9), title = true)
    runs = CSV.read(SUMMARY_PATH, DataFrame)
    runs = runs[(runs.n_learned .== 4) .& (runs.SA .== "inner_corner"), :]
    times = combine(groupby(runs, [:regime, :L, :ρ₀, :δλ]), :t_max => mean => :tmax, :t_max => std => :tmax_sd,
        :t_95 => mean => :t95, :t_95 => std => :t95_sd)
    δλs = sort(unique(times.δλ))
    fig = Figure(size = (1500, 1500))
    for (i, ρ₀) in enumerate(Ρs), (j, L) in enumerate(Ls)
        ax = Axis(fig[i, j]; yscale = log10, xticks = δλs, xticklabelrotation = π / 4,
            title = "($('a' + 2(i - 1) + j - 1)) L = $L, ρ₀ = $ρ₀", titlesize = 22,
            xlabel = i == length(Ρs) ? "Heterogeneity δλ  (λ = 0.5 ± δλ)" : "", ylabel = j == 1 ? "Steps" : "",
            xlabelsize = 22, ylabelsize = 22, xticklabelsize = 16, yticklabelsize = 17)
        for regime in ("TI", "PI")
            rows = sort(rows_where(times; regime, L, ρ₀), :δλ)
            series_color = regime_color(regime)
            name = regime == "TI" ? "Traditional" : "Peer (inner corner)"
            for (col, sd, ls, mk, lab) in ((:tmax, :tmax_sd, :solid, :circle, "tₘₐₓ (last learner)"),
                    (:t95, :t95_sd, :dash, :utriangle, "t₉₅ (95% learned)"))
                band!(ax, rows.δλ, max.(rows[!, col] .- rows[!, sd], 1), rows[!, col] .+ rows[!, sd]; color = (series_color, 0.15))
                scatterlines!(ax, rows.δλ, rows[!, col]; color = series_color, linestyle = ls, marker = mk, markersize = 10,
                    linewidth = 2.5, label = "$name, $lab")
            end
        end
        i == 1 && j == 1 && axislegend(ax; position = :lt, labelsize = 15, framevisible = false, rowgap = 0)
    end
    title && Label(fig[0, :], "Time to learn: last learner (tₘₐₓ) vs 95% of the class (t₉₅)
TI vs PI, inner-corner seating, 4 initial seeds, λ₀ = 0.5"; fontsize = 26)
    return fig
end

#* Ratio t_max / t₉₅ over every (L, ρ₀, δλ) — how much of t_max is the straggler tail
function t95_ratio_table()
    runs = CSV.read(SUMMARY_PATH, DataFrame)
    runs = runs[(runs.n_learned .== 4) .& (runs.SA .== "inner_corner"), :]
    ratios = combine(groupby(runs, [:regime, :L, :δλ]), [:t_max, :t_95] => ((a, b) -> mean(a ./ b)) => :tmax_over_t95)
    return unstack(ratios, [:regime, :L], :δλ, :tmax_over_t95)
end

#! Main
begin
    mkpath(PLOT_DIR)
    save("$(OUTPUT_DIR)/analysis/plots/class-split-learning-curves.png", figure_class_split(); px_per_unit = 2)
    save("$(OUTPUT_DIR)/analysis/plots/seed-density-comparison.png", figure_seed_density(); px_per_unit = 2)
    save_clean_figures()
    fig3 = figure3()
    save("$(PLOT_DIR)/figure3-learned-t-v2.png", fig3; px_per_unit = 2)
    fig5 = figure5()
    GLMakie.save("$(OUTPUT_DIR)/analysis/plots/figure5-tmax-surface-v2.png", fig5; px_per_unit = 2)
    fig5b = figure5(; stat = :t_95, zlabel = "Time to 95% learned t₉₅")
    GLMakie.save("$(OUTPUT_DIR)/analysis/plots/figure5-t95-surface-v2.png", fig5b; px_per_unit = 2)
    #* paper Fig. 6: t_max (top) and t₉₅ (bottom) surfaces, rendered row by row and stacked
    save_figure5_with_t95("$(OUTPUT_DIR)/analysis/plots/figure5-tmax-t95-surface-v2.png")
    save_figure5_with_t95("$(OUTPUT_DIR)/analysis/plots/clean/figure5-tmax-t95-clean.png"; legend = false)
    CairoMakie.activate!()
    save("$(OUTPUT_DIR)/analysis/plots/t95-vs-tmax.png", figure_t95_vs_tmax(); px_per_unit = 2)
    #* paper Fig. 5: snapshots of one TI/PI twin pair
    save("$(OUTPUT_DIR)/analysis/plots/figure4-snapshots-v2.png", figure4_snapshots(); px_per_unit = 2)
    save("$(OUTPUT_DIR)/analysis/plots/clean/figure4-snapshots-clean.png", figure4_snapshots(; labels = false); px_per_unit = 2)
    show(t95_ratio_table(); allrows = true); println()
    CairoMakie.activate!()
    println("Saved to $(PLOT_DIR)")
end
