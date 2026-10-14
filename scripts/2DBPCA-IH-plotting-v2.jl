# ! USE BS THESIS ENVIRONMENT (file paths will break if not)
#=
2DBPCA-IH v2: paper figures regenerated from the v2 sweep.
  * Fig. 3 — learned-fraction time series (figure3_learned_t below)
  * Fig. 5 — t_max surfaces over (ρ₀, δλ) for TI vs PI inner corner, L = 32 and 128 (figure5)

Same four panels and visual encoding as the v1 figure (2DBPCA-IH-trad-PI-comparison.jl):
  (a) lattice side L      — TI vs PI inner corner, ρ₀=0.5, δλ=0.2, L ∈ {32, 64, 128}
  (b) seating arrangement — L=64, ρ₀=0.5, δλ=0.2
  (c) coupling ρ₀         — inner corner, L=64, δλ=0.2, ρ₀ ∈ {0.1, 0.5, 0.9}
  (d) heterogeneity δλ    — inner corner, L=64, ρ₀=0.5, δλ ∈ {0, 0.2, 0.4, 0.45, 0.49}
Color = regime/SA, line width = L, opacity = ρ₀, line style = δλ. Lines are the mean over
the 20 trials; bands are ±1 SD. Runs are padded with f = 1 after absorption.

The x axis is iteration t+1 on a log₂ scale, so 2⁰ is the initial state (as in v1).
TI now starts from the same 4 seeds as PI, so its curves start at 4/N rather than 0.

Also writes a v1-vs-v2 comparison for the panels the revision changes most: (d), plus
panel (a) at δλ = 0.4.
=#

begin
    using Pkg
    Pkg.activate(".")
    include("2DBPCA-IH-common-v2.jl")
    using CairoMakie, ColorSchemes, Statistics
    import GLMakie #* only for the 3D surfaces (needs real depth sorting)
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

#* Mean and SD of f(t) over trials, t = 0..t_end (runs padded with their final value).
function mean_curve(regime, SA, L, ρ₀, δλ; n_learned = 4, λ₀ = 0.5)
    row = (; regime, SA, L, n_learned, ρ₀, λ₀, δλ)
    raw = read_parquet(raw_path(row))
    Ts = [read_trial(row, trial; raw)[1] for trial in 1:N_TRIALS]
    t_end = maximum(maximum.(Ts))
    F = reduce(hcat, [learned_fraction(T, t_end) for T in Ts])
    return vec(mean(F; dims = 2)), vec(std(F; dims = 2))
end

#* v1 curves from the old per-trial fit_params.csv (learned_per_gen), padded with 1.
function mean_curve_v1(SA, L, ρ₀, δλ)
    dir = SA == "TI" ? "traditional-$(L)" : "$(SA)-$(L)"
    file = SA == "TI" ? "2DBPCAIH-traditional-$(L)" : "2DBPCAIH-$(SA)-$(L)"
    runs = [CSV.read("./output/2D-Binary-PCA-IH/$(dir)/$(ρ₀)-0.5-$(δλ)/trial_$(tr)/data/$(file)-$(ρ₀)-0.5-$(δλ)-trial_$(tr)-fit_params.csv",
        DataFrame).learned_per_gen for tr in 1:N_TRIALS]
    n = maximum(length.(runs))
    F = reduce(hcat, [[r; fill(r[end], n - length(r))] for r in runs])
    return vec(mean(F; dims = 2)), vec(std(F; dims = 2))
end

function draw_series!(ax, μ, σ; regime, SA, L, ρ₀, δλ, label, solid = false)
    c = COLORS[SA_COLOR[regime == "TI" ? "TI" : SA]]
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

function figure3()
    fig = Figure(size = (1650, 1250))
    #* (a) lattice size
    ax = panel_axis(fig[1, 1]; letter = "(a)")
    for L in (32, 64, 128), regime in ("TI", "PI")
        μ, σ = mean_curve(regime, "inner_corner", L, 0.5, 0.2)
        draw_series!(ax, μ, σ; regime, SA = "inner_corner", L, ρ₀ = 0.5, δλ = 0.2, label = "")
    end
    label!(ax, "L = 32, 64, 128", 2^0.5, 0.66, "TI")
    label!(ax, "L = 32", 2^2.2, 0.48, "inner_corner")
    label!(ax, "L = 64", 2^4.3, 0.55, "inner_corner")
    label!(ax, "L = 128", 2^5.55, 0.80, "inner_corner")
    #* (b) seating arrangement
    ax = panel_axis(fig[2, 1]; xlabel = "Iteration number", letter = "(b)")
    μ, σ = mean_curve("TI", "inner_corner", 64, 0.5, 0.2)
    draw_series!(ax, μ, σ; regime = "TI", SA = "inner_corner", L = 64, ρ₀ = 0.5, δλ = 0.2, label = "")
    for SA in ("inner_corner", "center", "outer_corner", "random")
        μ, σ = mean_curve("PI", SA, 64, 0.5, 0.2)
        draw_series!(ax, μ, σ; regime = "PI", SA, L = 64, ρ₀ = 0.5, δλ = 0.2, label = "")
    end
    label!(ax, "Traditional", 2^1.0, 0.66, "TI")
    label!(ax, "Inner Corner", 2^2.1, 0.29, "inner_corner")
    label!(ax, "Random", 2^5.6, 0.77, "random")
    label!(ax, "Center", 2^5.3, 0.46, "center")
    label!(ax, "Outer Corner", 2^4.6, 0.14, "outer_corner")
    #* (c) coupling
    ax = panel_axis(fig[1, 2]; letter = "(c)")
    for ρ₀ in (0.9, 0.5, 0.1), regime in ("TI", "PI")
        μ, σ = mean_curve(regime, "inner_corner", 64, ρ₀, 0.2)
        draw_series!(ax, μ, σ; regime, SA = "inner_corner", L = 64, ρ₀, δλ = 0.2, label = "")
    end
    label!(ax, "ρ = 0.9", 2^0.05, 0.68, "TI"; ρ₀ = 0.9)
    label!(ax, "ρ = 0.5", 2^2.05, 0.47, "TI"; ρ₀ = 0.5)
    label!(ax, "ρ = 0.1", 2^1.2, 0.21, "TI"; ρ₀ = 0.1)
    label!(ax, "ρ = 0.9", 2^2.6, 0.69, "inner_corner"; ρ₀ = 0.9)
    label!(ax, "ρ = 0.5", 2^4.55, 0.44, "inner_corner"; ρ₀ = 0.5)
    label!(ax, "ρ = 0.1", 2^6.25, 0.66, "inner_corner"; ρ₀ = 0.1)
    #* (d) heterogeneity
    ax = panel_axis(fig[2, 2]; xlabel = "Iteration number", letter = "(d)")
    for δλ in (0.0, 0.2, 0.4, 0.45, 0.49), regime in ("TI", "PI")
        μ, σ = mean_curve(regime, "inner_corner", 64, 0.5, δλ)
        draw_series!(ax, μ, σ; regime, SA = "inner_corner", L = 64, ρ₀ = 0.5, δλ, label = "")
    end
    label!(ax, "δλ = 0.0, 0.2", 2^0.55, 0.90, "TI"; fontsize = 20)
    label!(ax, "δλ = 0.4", 2^2.85, 0.645, "TI"; rotation = 0.6, fontsize = 18)
    label!(ax, "δλ = 0.45", 2^3.2, 0.575, "TI"; rotation = 0.45, fontsize = 18)
    label!(ax, "δλ = 0.49", 2^8.6, 0.78, "TI")
    label!(ax, "δλ = 0.0, 0.2\n0.4, 0.45", 2^1.25, 0.18, "inner_corner")
    label!(ax, "δλ = 0.49", 2^5.9, 0.52, "inner_corner")
    Label(fig[:, 0], "Fraction of learned"; fontsize = 32, rotation = π / 2, tellheight = false)
    return fig
end

#* v1 (thin, faded) vs v2 (thick) for PI and TI at the δλ values the revision changes most.
function figure_v1_vs_v2()
    fig = Figure(size = (1600, 650))
    for (j, regime) in enumerate(("PI", "TI"))
        ax = panel_axis(fig[1, j]; xlabel = "Iteration number", ylabel = j == 1 ? "Fraction learned" : "",
            letter = regime == "PI" ? "PI (inner corner)" : "TI")
        for (k, δλ) in enumerate((0.0, 0.4, 0.45, 0.49))
            c = Makie.wong_colors()[k]
            μ1, _ = mean_curve_v1(regime == "TI" ? "TI" : "inner_corner", 64, 0.5, δλ)
            μ2, σ2 = mean_curve(regime, "inner_corner", 64, 0.5, δλ)
            lines!(ax, 1:length(μ1), μ1; color = (c, 0.6), linewidth = 2, linestyle = :dash, label = "v1 δλ=$δλ")
            band!(ax, 1:length(μ2), clamp.(μ2 .- σ2, 0, 1), clamp.(μ2 .+ σ2, 0, 1); color = (c, 0.15))
            lines!(ax, 1:length(μ2), μ2; color = c, linewidth = 3.5, label = "v2 δλ=$δλ")
        end
        if regime == "TI" #* R3's analytical F(t) at δλ = 0.49 (dotted black)
            s0 = 4 / 64^2; p₊, p₋ = 0.5 * 0.99, 0.5 * 0.01
            t = 0:2^11
            lines!(ax, t .+ 1, s0 .+ (1 - s0) .* (1 .- 0.5 .* (1 - p₊) .^ t .- 0.5 .* (1 - p₋) .^ t);
                color = :black, linestyle = :dot, linewidth = 2, label = "F(t), δλ=0.49")
        end
        add_legend!(ax)
    end
    Label(fig[0, :], "L = 64, ρ₀ = 0.5: v1 (dashed: neighbor-λ bug, stopping rule) vs v2 (solid ± SD: learner λ, run to absorption)",
        fontsize = 20)
    return fig
end


#* Fig. 5: mean t_max over (ρ₀, δλ) for TI (blue) and PI inner corner (orange), L = 32 and 128.
#* δλ ≤ 0.4 as in the published figure (pass δλs to include 0.45/0.49). Rendered with GLMakie
#* because the two surfaces intersect and CairoMakie does not depth-sort between plots.
function figure5(; δλs = [0.0, 0.1, 0.2, 0.3, 0.4], Ρs = round.(collect(0.1:0.1:1.0); digits = 1))
    GLMakie.activate!(; scalefactor = 1) #* ignore OS display scaling so the layout matches the requested size
    sm = CSV.read(SUMMARY_PATH, DataFrame)
    sm = sm[(sm.n_learned .== 4) .& (sm.SA .== "inner_corner"), :]
    g = combine(groupby(sm, [:regime, :L, :ρ₀, :δλ]), :t_max => mean => :t)
    tgrid(regime, L) = [only(g[(g.regime .== regime) .& (g.L .== L) .& isapprox.(g.ρ₀, ρ) .& isapprox.(g.δλ, δ), :t])
        for ρ in Ρs, δ in δλs]

    fig = Figure(size = (1900, 950), figure_padding = (20, 20, 20, 20))
    for (j, (L, letter)) in enumerate(((32, "(a)"), (128, "(b)")))
        ax = Axis3(fig[1, j]; xlabel = "Transmission ρ₀", ylabel = "Heterogeneity δλ", zlabel = "Time to learn tₘₐₓ",
            aspect = (1, 1, 1), azimuth = 1.275π, viewmode = :fit, xlabelsize = 30, ylabelsize = 30, zlabelsize = 30,
            xticklabelsize = 22, yticklabelsize = 22, zticklabelsize = 22, xlabeloffset = 55, ylabeloffset = 55,
            zlabeloffset = 70)
        for regime in ("TI", "PI")
            surface!(ax, Ρs, δλs, tgrid(regime, L); color = fill(COLORS[SA_COLOR[regime == "TI" ? "TI" : "inner_corner"]], length(Ρs), length(δλs)),
                shading = MultiLightShading)
        end
        Label(fig[1, j, Top()], "$letter L = $L"; fontsize = 40, halign = :left, padding = (40, 0, 0, 0))
    end
    Legend(fig[1, 1], [PolyElement(color = COLORS[1]), PolyElement(color = COLORS[2])],
        [rich("Traditional"; color = COLORS[1], font = :bold), rich("Peer"; color = COLORS[2], font = :bold)];
        tellwidth = false, tellheight = false, halign = :right, valign = :top, margin = (0, 40, 0, 60),
        labelsize = 34, patchsize = (55, 55), rowgap = 30, framevisible = false)
    return fig
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
    axb = Axis(fig[2, 1]; ax_kw..., xlabel = "Iteration number")
    axc = Axis(fig[1, 2]; ax_kw..., xticksvisible = false, xticklabelsvisible = false)
    axd = Axis(fig[2, 2]; ax_kw..., xlabel = "Iteration number")
    for L in (32, 64, 128), regime in ("TI", "PI")
        μ, σ = mean_curve(regime, "inner_corner", L, 0.5, 0.2)
        draw_series!(axa, μ, σ; regime, SA = "inner_corner", L, ρ₀ = 0.5, δλ = 0.2, label = "")
    end
    μ, σ = mean_curve("TI", "inner_corner", 64, 0.5, 0.2)
    draw_series!(axb, μ, σ; regime = "TI", SA = "inner_corner", L = 64, ρ₀ = 0.5, δλ = 0.2, label = "")
    for SA in ("inner_corner", "outer_corner", "center", "random")
        μ, σ = mean_curve("PI", SA, 64, 0.5, 0.2)
        draw_series!(axb, μ, σ; regime = "PI", SA, L = 64, ρ₀ = 0.5, δλ = 0.2, label = "")
    end
    for ρ₀ in (0.1, 0.5, 0.9), regime in ("TI", "PI")
        μ, σ = mean_curve(regime, "inner_corner", 64, ρ₀, 0.2)
        draw_series!(axc, μ, σ; regime, SA = "inner_corner", L = 64, ρ₀, δλ = 0.2, label = "")
    end
    for δλ in (0.0, 0.2, 0.4, 0.45, 0.49), regime in ("TI", "PI")
        μ, σ = mean_curve(regime, "inner_corner", 64, 0.5, δλ)
        draw_series!(axd, μ, σ; regime, SA = "inner_corner", L = 64, ρ₀ = 0.5, δλ, label = "")
    end
    foreach(((ax, s),) -> letter!(ax, s), zip((axa, axb, axc, axd), ("(a)", "(b)", "(c)", "(d)")))
    linkxaxes!(axa, axb)
    linkxaxes!(axc, axd)
    Label(fig[:, 0], "Fraction of learned"; fontsize = 40, rotation = π / 2, tellheight = false)
    return fig
end

#* One Fig. 5 panel as v1's {L}_truncated.png: 1000×800 at 2.5 px/unit, same Axis3 settings,
#* no title or legend. GLMakie (depth sorting); call with GLMakie active.
function figure5_panel(L; δλs = [0.0, 0.1, 0.2, 0.3, 0.4], Ρs = round.(collect(0.1:0.1:1.0); digits = 1))
    sm = CSV.read(SUMMARY_PATH, DataFrame)
    sm = sm[(sm.n_learned .== 4) .& (sm.SA .== "inner_corner") .& (sm.L .== L), :]
    g = combine(groupby(sm, [:regime, :ρ₀, :δλ]), :t_max => mean => :t)
    tgrid(regime) = [only(g[(g.regime .== regime) .& isapprox.(g.ρ₀, ρ) .& isapprox.(g.δλ, δ), :t]) for ρ in Ρs, δ in δλs]
    t_top = maximum(maximum(tgrid(r)) for r in ("TI", "PI"))
    fig = Figure(size = (1000, 800))
    ax = Axis3(fig[1, 1]; xlabel = "Transmission ρ₀", ylabel = "Heterogeneity δλ", zlabel = "Time to learn tₘₐₓ",
        aspect = (1, 1, 1), zticks = 0:(t_top > 1000 ? 200 : 100):t_top, xlabelsize = 32, ylabelsize = 32,
        zlabelsize = 32, xlabeloffset = 60, ylabeloffset = 60, zlabeloffset = 80, xticklabelsize = 28,
        yticklabelsize = 28, zticklabelsize = 28, viewmode = :fit, azimuth = 1.275π)
    for regime in ("TI", "PI")
        surface!(ax, Ρs, δλs, tgrid(regime); color = fill(COLORS[SA_COLOR[regime == "TI" ? "TI" : "inner_corner"]],
            length(Ρs), length(δλs)), shading = FastShading)
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
    end
    CairoMakie.activate!()
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
    sm = CSV.read(SUMMARY_PATH, DataFrame)
    g = combine(groupby(sm, [:regime, :SA, :L, :n_learned, :ρ₀, :δλ]), :t_max => mean => :t, :t_max => std => :s)
    stat(regime, SA, L, n, ρ, δ) = (r = g[(g.regime .== regime) .& (g.SA .== SA) .& (g.L .== L) .& (g.n_learned .== n) .&
        isapprox.(g.ρ₀, ρ) .& isapprox.(g.δλ, δ), :]; (r.t[1], r.s[1]))
    series(regime, SA, ρ, δ, frac) = [stat(regime, SA, L, frac ? n_frac(L) : 4, ρ, δ) for L in Ls]

    fig = Figure(size = (1600, 1250))
    for (j, (ρ, δ)) in enumerate(((0.5, 0.2), (0.1, 0.4)))
        ax = Axis(fig[1, j]; xscale = log2, yscale = log10, xticks = (Ls, string.(Ls)), xlabel = "Lattice side L",
            ylabel = j == 1 ? "t_max (mean ± SD, 20 runs)" : "", title = "($('a' + j - 1)) PI t_max vs L, ρ₀ = $ρ, δλ = $δ",
            titlesize = 22, xlabelsize = 22, ylabelsize = 22, xticklabelsize = 18, yticklabelsize = 18)
        for frac in (false, true)
            ts = first.(series("TI", "inner_corner", ρ, δ, frac))
            lines!(ax, Ls, ts; color = (:gray40, 0.8), linewidth = 2, linestyle = frac ? :dash : :solid,
                label = frac ? "TI, fixed fraction" : "TI, 4 seeds")
        end
        for SA in SA_ORDER, frac in (false, true)
            st = series("PI", SA, ρ, δ, frac); t = first.(st); s = last.(st)
            c = COLORS[SA_COLOR[SA]]
            band!(ax, Ls, max.(t .- s, 1), t .+ s; color = (c, 0.12))
            scatterlines!(ax, Ls, t; color = c, linewidth = 3, linestyle = frac ? :dash : :solid, markersize = 10,
                marker = frac ? :utriangle : :circle,
                label = "$(SA_LABEL[SA]), $(frac ? "fixed fraction" : "4 seeds") (slope $(round(loglog_slope(t); digits = 2)))")
        end
        axislegend(ax; position = :lt, labelsize = 13, framevisible = false, backgroundcolor = :transparent, rowgap = 0)
    end

    #* (c) learning curves at L = 128
    ax = Axis(fig[2, 1]; xscale = log2, xlabel = "Iteration number", ylabel = "Fraction learned", yticks = 0:0.5:1,
        title = "(c) PI learning curves, L = 128, ρ₀ = 0.5, δλ = 0.2", titlesize = 22, xlabelsize = 22, ylabelsize = 22,
        xticklabelsize = 18, yticklabelsize = 18)
    for SA in SA_ORDER, frac in (false, true)
        μ, σ = mean_curve("PI", SA, 128, 0.5, 0.2; n_learned = frac ? n_frac(128) : 4)
        c = COLORS[SA_COLOR[SA]]
        band!(ax, 1:length(μ), clamp.(μ .- σ, 0, 1), clamp.(μ .+ σ, 0, 1); color = (c, 0.12))
        lines!(ax, 1:length(μ), μ; color = c, linewidth = 3, linestyle = frac ? :dash : :solid,
            label = "$(SA_LABEL[SA]), $(frac ? "$(n_frac(128)) seeds" : "4 seeds")")
    end
    axislegend(ax; position = :lt, labelsize = 13, framevisible = false, backgroundcolor = :transparent, rowgap = 0)

    #* (d) slope with fixed fraction vs slope with 4 seeds, every (ρ₀, δλ)
    ax = Axis(fig[2, 2]; xlabel = "log–log slope of t_max vs L, 4 seeds", ylabel = "slope, fixed seed fraction",
        title = "(d) Size scaling, all ρ₀ × δλ (PI; TI in gray)", titlesize = 22, xlabelsize = 22, ylabelsize = 22,
        xticklabelsize = 18, yticklabelsize = 18, aspect = DataAspect())
    ρs = round.(collect(0.1:0.1:1.0); digits = 1); δs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49]
    slopes(regime, SA) = [(loglog_slope(first.(series(regime, SA, ρ, δ, false))), loglog_slope(first.(series(regime, SA, ρ, δ, true))))
        for ρ in ρs, δ in δs] |> vec
    sTI = slopes("TI", "inner_corner")
    scatter!(ax, first.(sTI), last.(sTI); color = (:gray40, 0.6), markersize = 8, label = "TI")
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
#* Classroom evolution: fastest vs slowest run (by t_max) of one parameter set, side by side
#* ---------------------------------------------------------------------------
#* Cell categories: 1 unlearned fast, 2 unlearned slow, 3 learned fast, 4 learned slow, 5 seed.
const CELL_COLORS = [colorant"#dbe9f6", colorant"#b4b4b4", colorant"#08306b", colorant"#4a90d9", colorant"#e6550d"]

function cell_state(T, cls, t)
    S = similar(T, Int8)
    for i in eachindex(T)
        S[i] = T[i] == 0 ? 5 : (0 < T[i] <= t ? (cls[i] == 1 ? 3 : 4) : (cls[i] == 1 ? 1 : 2))
    end
    return S
end

function best_worst_runs(regime, SA, L, ρ₀, δλ; n_learned = 4, λ₀ = 0.5)
    row = (; regime, SA, L, n_learned, ρ₀, λ₀, δλ)
    raw = read_parquet(raw_path(row))
    runs = [read_trial(row, tr; raw) for tr in 1:N_TRIALS]
    tmax = [maximum(r[1]) for r in runs]
    ib, iw = argmin(tmax), argmax(tmax)
    return runs[ib], runs[iw], (ib, iw), tmax
end

#* Builds the figure for an observable time t; returns (fig, t_obs, t_end).
function evolution_figure(regime, SA, L, ρ₀, δλ; n_learned = 4)
    (Tb, cb), (Tw, cw), (ib, iw), tmax = best_worst_runs(regime, SA, L, ρ₀, δλ; n_learned)
    t_end = Int(maximum(Tw))
    t = Observable(0)
    fig = Figure(size = (1400, 1000))
    Label(fig[0, 1:2], "$regime, $(SA_LABEL[SA]), L = $L, n = $n_learned, ρ₀ = $ρ₀, δλ = $δλ — fastest vs slowest of $(N_TRIALS) runs " *
        "(t_max range $(minimum(tmax))–$(maximum(tmax)), median $(Int(round(median(tmax)))))"; fontsize = 20)
    for (j, (T, c, tr, name)) in enumerate(((Tb, cb, ib, "Fastest"), (Tw, cw, iw, "Slowest")))
        ax = Axis(fig[1, j]; aspect = DataAspect(), yreversed = true, title = @lift("$name run (trial $tr, t_max = $(maximum(T)))   t = $($t)"),
            titlesize = 20)
        hidedecorations!(ax)
        heatmap!(ax, @lift(cell_state(T, c, $t)); colormap = CELL_COLORS, colorrange = (1, 5))
    end
    #* Learning curves on a log₂ axis (x = t + 1, as in Fig. 3): whole class (thin), fast learners
    #* (solid), slow learners (dashed), for each run.
    ax = Axis(fig[2, 1:2]; xscale = log2, xlabel = "Iteration number (t + 1)", ylabel = "Fraction learned", yticks = 0:0.5:1,
        xlabelsize = 20, ylabelsize = 20)
    x = 1:(t_end + 1)
    for (T, c, col, name) in ((Tb, cb, COLORS[3], "fastest"), (Tw, cw, COLORS[4], "slowest"))
        f = learned_fraction(T, t_end)
        f₊ = learned_fraction_by_class(T, c, 1, t_end)
        f₋ = learned_fraction_by_class(T, c, 2, t_end)
        lines!(ax, x, f; color = (col, 0.5), linewidth = 1.5, label = "$name: whole class")
        lines!(ax, x, f₊; color = col, linewidth = 3, label = "$name: fast λ")
        lines!(ax, x, f₋; color = col, linewidth = 3, linestyle = :dash, label = "$name: slow λ")
        scatter!(ax, @lift([Point2f($t + 1, f₊[$t + 1]), Point2f($t + 1, f₋[$t + 1])]); color = col, markersize = 12)
    end
    vlines!(ax, @lift($t + 1); color = :black, linestyle = :dash)
    axislegend(ax; position = :lt, framevisible = false, nbanks = 2, labelsize = 14)
    Legend(fig[1, 3], [PolyElement(color = c) for c in CELL_COLORS],
        ["unlearned, fast λ", "unlearned, slow λ", "learned, fast λ", "learned, slow λ", "initial seed"]; framevisible = false)
    rowsize!(fig.layout, 2, Relative(0.25))
    return fig, t, t_end
end

#* Long runs are subsampled to at most max_frames evenly spaced steps (title shows the true t).
function animate_evolution(regime, SA, L, ρ₀, δλ; n_learned = 4, framerate = 8, max_frames = 300, strip_times = nothing)
    dir = "$(OUTPUT_DIR)/analysis/plots/evolution"
    mkpath(dir)
    name = "$(regime)-$(SA)-$(L)-$(n_learned)-$(ρ₀)-$(δλ)"
    fig, t, t_end = evolution_figure(regime, SA, L, ρ₀, δλ; n_learned)
    frames = t_end + 1 <= max_frames ? (0:t_end) : unique(round.(Int, range(0, t_end; length = max_frames)))
    record(fig, "$(dir)/$(name).mp4", frames; framerate) do τ
        t[] = τ
    end
    #* Still strip of selected frames (for quick viewing without a video player)
    times = strip_times === nothing ? round.(Int, range(0, t_end; length = 6)) : strip_times
    strip = Figure(size = (300 * length(times), 640))
    (Tb, cb), (Tw, cw), _, _ = best_worst_runs(regime, SA, L, ρ₀, δλ; n_learned)
    for (k, τ) in enumerate(times), (j, (T, c, nm)) in enumerate(((Tb, cb, "fastest"), (Tw, cw, "slowest")))
        ax = Axis(strip[j, k]; aspect = DataAspect(), yreversed = true, title = "$nm, t = $τ", titlesize = 16)
        hidedecorations!(ax)
        heatmap!(ax, cell_state(T, c, τ); colormap = CELL_COLORS, colorrange = (1, 5))
    end
    save("$(dir)/$(name)-strip.png", strip)
    println("Saved $(dir)/$(name).mp4 ($(length(frames)) frames for t = 0..$(t_end)) and strip")
end

#* ---------------------------------------------------------------------------
#* Fast vs slow learners: mean learned fraction per λ class, TI vs PI (exploration)
#* ---------------------------------------------------------------------------
#* Mean over trials of the per-class learned fraction (runs padded with their final value).
function mean_class_curve(regime, SA, L, ρ₀, δλ, class; n_learned = 4, λ₀ = 0.5)
    row = (; regime, SA, L, n_learned, ρ₀, λ₀, δλ)
    raw = read_parquet(raw_path(row))
    runs = [read_trial(row, tr; raw) for tr in 1:N_TRIALS]
    t_end = maximum(maximum(r[1]) for r in runs)
    return vec(mean(reduce(hcat, [learned_fraction_by_class(T, c, class, t_end) for (T, c) in runs]); dims = 2))
end

function figure_class_split(; SA = "inner_corner", L = 64, Ρs = (0.5, 0.1), δλs = (0.2, 0.4, 0.49))
    fig = Figure(size = (1800, 1000))
    for (i, ρ₀) in enumerate(Ρs), (j, δλ) in enumerate(δλs)
        ax = Axis(fig[i, j]; xscale = log2, yticks = 0:0.5:1, title = "ρ₀ = $ρ₀, δλ = $δλ", titlesize = 22,
            xlabel = i == length(Ρs) ? "Iteration number" : "", ylabel = j == 1 ? "Fraction learned" : "",
            xlabelsize = 20, ylabelsize = 20, xticklabelsize = 16, yticklabelsize = 16)
        for regime in ("TI", "PI")
            col = COLORS[SA_COLOR[regime == "TI" ? "TI" : SA]]
            μ, _ = mean_curve(regime, SA, L, ρ₀, δλ)
            f₊ = mean_class_curve(regime, SA, L, ρ₀, δλ, 1)
            f₋ = mean_class_curve(regime, SA, L, ρ₀, δλ, 2)
            lines!(ax, 1:length(μ), μ; color = (col, 0.45), linewidth = 1.5, label = "$regime whole class")
            lines!(ax, 1:length(f₊), f₊; color = col, linewidth = 3, label = "$regime fast λ = $(0.5 + δλ)")
            lines!(ax, 1:length(f₋), f₋; color = col, linewidth = 3, linestyle = :dash, label = "$regime slow λ = $(round(0.5 - δλ; digits = 2))")
        end
        axislegend(ax; position = :lt, framevisible = false, labelsize = 13, rowgap = 0)
    end
    Label(fig[0, :], "Fast vs slow learners: mean of $(N_TRIALS) runs, PI = $(SA_LABEL[SA]), L = $L"; fontsize = 24)
    return fig
end

#! Main
begin
    mkpath(PLOT_DIR)
    save("$(OUTPUT_DIR)/analysis/plots/class-split-learning-curves.png", figure_class_split(); px_per_unit = 2)
    save("$(OUTPUT_DIR)/analysis/plots/seed-density-comparison.png", figure_seed_density(); px_per_unit = 2)
    save_clean_figures()
    fig3 = figure3()
    save("$(PLOT_DIR)/figure3-learned-t-v2.png", fig3; px_per_unit = 2)
    figc = figure_v1_vs_v2()
    save("$(PLOT_DIR)/learned-t-v1-vs-v2.png", figc; px_per_unit = 2)
    fig5 = figure5()
    GLMakie.save("$(OUTPUT_DIR)/analysis/plots/figure5-tmax-surface-v2.png", fig5; px_per_unit = 2)
    CairoMakie.activate!()
    println("Saved to $(PLOT_DIR)")
end
