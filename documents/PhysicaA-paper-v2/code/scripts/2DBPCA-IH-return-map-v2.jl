# Run from the repository root (paths are relative to it); see README.md.
#=
2DBPCA-IH v2 return-map fits (paper §3.4, Fig. 6, Table 3), from the v2 sweep.

Same functional forms as v1 (2DBPCA-IH-return-map-stats.jl), with h(1) = 1 imposed for the
simulations because f = 1 is the absorbing state:
  TI: h(f) = f + a₀(1 − f)              → constant gain g = a₀
  PI: h(f) = f + (1 − f)(a₀ − a₂ f)     → linear gain g(f) = a₀ − a₂ f

Uncertainty (reviewer R1): consecutive return-map points from one trajectory are correlated, so
a single fit pooled over all runs badly understates the coefficient uncertainty (its standard
error assumes every point is independent). Here every run is fitted separately; the 20 runs
are independent (own seeds), so
  mean ± SE (= SD/√n) is the uncertainty of the typical coefficient (for Table 3), and
  SD is the run-to-run spread (for comparisons against individual empirical classes).
The pooled fit and its naive standard error are also written, for comparison.

PI uses the inner-corner seating arrangement, as in the paper. Fits use each run's own points
f_t for t = 0..t_max (no padding after absorption).

Outputs (output/2D-Binary-PCA-IH-v2/analysis/):
  return-map-per-trial.csv  one row per (regime, L, ρ₀, δλ, trial): coefficients, R², n_points
  return-map-summary.csv    one row per (regime, L, ρ₀, δλ): mean/SD/SE/median of each
                            coefficient over trials, plus the pooled fit and its naive SE
  analysis/plots/return-map/  Fig. 6 option 1 (student-level empirical cloud, as published) and
                            option 2 (class-level empirical comparison); clean Fig. 6 panels at the
                            v1 source dimensions go to analysis/plots/clean/

Empirical data (Roxas; available from the authors on request, not included here): place
data/empirical/{in,out}_8_<k>.txt in the repository root to include them. The files hold
the first and second assessment scores of the 64 students of 8×8 class k (k = 1..14). The model's f
is a class-level learned fraction, so the class-level analogue is (mean pre, mean post) per class.
The model cannot produce score decreases (learning is irreversible by design).
=#

begin
    using Pkg
    Pkg.activate(".")
    include("2DBPCA-IH-common-v2.jl")
    using LsqFit, Statistics
    using CairoMakie, ColorSchemes, Distributions
end

ti_model(f, p) = f .+ p[1] .* (1 .- f)
pi_model(f, p) = f .+ (1 .- f) .* (p[1] .- p[2] .* f)

r_squared(y, ŷ) = 1 - sum((y .- ŷ) .^ 2) / sum((y .- mean(y)) .^ 2)

#* Paper copy of a figure: legends removed (and the empty legend row/column trimmed). The legend
#* version is kept for reference; series are identified in the caption or by annotation.
strip_legends!(fig) = (foreach(delete!, filter(c -> c isa Legend, copy(fig.content))); trim!(fig.layout); fig)

#* (f_{t-1}, f_t) pairs from one learning-time grid
function return_map_pairs(T)
    f = learned_fraction(T)
    return f[1:end-1], f[2:end]
end

#* Return-map pairs of every trial of a parameter row
trial_pairs(row) = [return_map_pairs(T) for (T, _) in read_trials(row)]

function fit_return_map(regime, x, y)
    model, p0 = regime == "TI" ? (ti_model, [0.5]) : (pi_model, [0.0, 1.0])
    fit = curve_fit(model, x, y, p0)
    se = try stderror(fit) catch; fill(NaN, length(p0)) end
    return fit.param, se, r_squared(y, model(x, fit.param))
end

coef_names(regime) = regime == "TI" ? ["a0"] : ["a0", "a2"]

function return_map_fits(; Ls = [32, 64, 128], Ρs = [0.1, 0.5, 0.9], δλs = [0.0, 0.2, 0.4, 0.45, 0.49],
        SA = "inner_corner", λ₀ = 0.5)
    per_trial = DataFrame()
    summary = DataFrame()
    for regime in ("TI", "PI"), L in Ls, ρ₀ in Ρs, δλ in δλs
        pairs = trial_pairs(param_row(regime, SA, L, ρ₀, δλ; λ₀))
        n_coefs = length(coef_names(regime))

        #* per-trial fits
        coefs = Matrix{Float64}(undef, N_TRIALS, n_coefs)
        r2s = zeros(N_TRIALS)
        for (tr, (x, y)) in enumerate(pairs)
            p, _, r2 = fit_return_map(regime, x, y)
            coefs[tr, :] = p
            r2s[tr] = r2
            a2 = regime == "TI" ? NaN : p[2]
            push!(per_trial, (; regime, SA, L, ρ₀, δλ, trial = tr, n_points = length(x), r_squared = r2, a0 = p[1], a2);
                cols = :union)
        end

        #* pooled fit (v1 method) for comparison
        pooled_x = reduce(vcat, first.(pairs)); pooled_y = reduce(vcat, last.(pairs))
        pooled_coefs, pooled_se, pooled_r2 = fit_return_map(regime, pooled_x, pooled_y)

        summary_row = Dict{Symbol, Any}(:regime => regime, :SA => SA, :L => L, :ρ₀ => ρ₀, :δλ => δλ, :n_trials => N_TRIALS,
            :mean_r_squared => mean(r2s), :pooled_n_points => length(pooled_x), :pooled_r_squared => pooled_r2)
        for (k, name) in enumerate(["a0", "a2"])
            #* TI has no a₂: its a2 columns are NaN
            stats = if k <= n_coefs
                trial_coefs = coefs[:, k]
                (mean(trial_coefs), std(trial_coefs), std(trial_coefs) / sqrt(N_TRIALS), median(trial_coefs), pooled_coefs[k], pooled_se[k])
            else
                ntuple(_ -> NaN, 6)
            end
            for (suffix, value) in zip(("mean", "sd", "se", "median", "pooled", "pooled_naive_se"), stats)
                summary_row[Symbol("$(name)_$(suffix)")] = value
            end
        end
        push!(summary, summary_row; cols = :union)
    end
    return per_trial, summary
end

#* ---------------------------------------------------------------------------
#* Plots
#* ---------------------------------------------------------------------------
const COLORS = ColorSchemes.seaborn_colorblind.colors
const EMP_DIR = "./data/empirical" #* optional; empirical overlays and the class-level figure are skipped without it
#* Fig. 6 grid: rows = unfavorable / favorable conditions (slowest / fastest corner of the δλ ≤ 0.4 sweep), columns = L.
const FIG6_LS = (32, 128)
const FIG6_ROWS = (("Unfavorable conditions", 0.1, 0.4), ("Favorable conditions", 0.9, 0.0))
const FIG6_CASES = Tuple((L, ρ₀, δλ) for (_, ρ₀, δλ) in FIG6_ROWS for L in FIG6_LS)

#* One score per line (blank lines skipped)
read_scores(path) = parse.(Float64, filter(!isempty, strip.(readlines(path))))

#* (first assessment, second assessment) scores of each of the 14 classes; empty if the data are absent
function read_empirical()
    isfile("$(EMP_DIR)/in_8_1.txt") || return Tuple{Vector{Float64}, Vector{Float64}}[]
    return [(read_scores("$(EMP_DIR)/in_8_$(k).txt"), read_scores("$(EMP_DIR)/out_8_$(k).txt")) for k in 1:14]
end

empirical_model(f, p) = p[1] .* f .^ 2 .+ p[2] .* f .+ p[3]
nitta(f) = f .+ f .* (1 .- f) #* Nitta (2019): g(f) = f

function mean_coefs(summary, regime, L, ρ₀, δλ)
    rows = rows_where(summary; regime, L, ρ₀, δλ)
    return regime == "TI" ? [rows.a0_mean[1]] : [rows.a0_mean[1], rows.a2_mean[1]]
end

#* One Fig. 6 panel (v1 style: 1000×707 canvas). Simulation points pooled over the 20 runs for
#* display; dashed lines are the constrained fits with per-run mean coefficients.
function return_map_panel!(ax, L, ρ₀, δλ, summary, emp, emp_p)
    fs = 0:0.01:1
    for (j, regime) in enumerate(("TI", "PI"))
        pairs = trial_pairs(param_row(regime, "inner_corner", L, ρ₀, δλ))
        scatter!(ax, reduce(vcat, first.(pairs)), reduce(vcat, last.(pairs)); color = (COLORS[j], 0.5),
            marker = j == 1 ? :circle : :rect, markersize = 16)
        coefs = mean_coefs(summary, regime, L, ρ₀, δλ)
        lines!(ax, fs, (regime == "TI" ? ti_model : pi_model)(fs, coefs); color = COLORS[j], linestyle = :dash, linewidth = 3)
    end
    lines!(ax, fs, nitta(fs); color = COLORS[3], linewidth = 3)
    lines!(ax, fs, fs; color = COLORS[8], linewidth = 3)
    isempty(emp) && return
    scatter!(ax, reduce(vcat, first.(emp)), reduce(vcat, last.(emp)); color = (COLORS[5], 0.3), marker = :diamond, markersize = 16)
    lines!(ax, fs, empirical_model(fs, emp_p); color = COLORS[5], linestyle = :dash, linewidth = 3)
end

v1_axis(pos; kw...) = Axis(pos; xlabel = "Fraction of learned f(t − 1)", ylabel = "Fraction of learned f(t)", xlabelsize = 32, ylabelsize = 32,
    xticklabelsize = 28, yticklabelsize = 28, kw...)

#* Option 1: Fig. 6 with the student-level empirical cloud (as published), now 2×2: rows are the
#* unfavorable and favorable conditions, columns are L = 32 and 128, so each row compares sizes at fixed parameters.
function figure6_option1(summary, emp, emp_p)
    fig = Figure(size = (2000, 1560))
    letters = "abcd"
    for (i, (case, ρ₀, δλ)) in enumerate(FIG6_ROWS), (j, L) in enumerate(FIG6_LS)
        k = 2(i - 1) + j
        ax = v1_axis(fig[i, j]; title = "($(letters[k])) $case: L = $L, ρ₀ = $ρ₀, λ = 0.5 ± $δλ", titlesize = 28,
            xlabel = i == 2 ? "Fraction of learned f(t − 1)" : "", xticklabelsvisible = i == 2)
        return_map_panel!(ax, L, ρ₀, δλ, summary, emp, emp_p)
    end
    Legend(fig[3, :], [MarkerElement(color = COLORS[1], marker = :circle, markersize = 16),
            MarkerElement(color = COLORS[2], marker = :rect, markersize = 16),
            MarkerElement(color = COLORS[5], marker = :diamond, markersize = 16),
            LineElement(color = COLORS[3], linewidth = 3), LineElement(color = COLORS[8], linewidth = 3),
            LineElement(color = :black, linestyle = :dash, linewidth = 3)],
        ["TI (20 runs)", "PI, inner corner (20 runs)", "Roxas: individual students (14 classes, 64 students each)",
            "Nitta: g(f) = f", "y = x", "fits (simulation: constrained, per-run mean; Roxas: free quadratic)"];
        orientation = :horizontal, nbanks = 2, labelsize = 22, framevisible = false)
    return fig
end

#* Option 2: class-level empirical comparison.
#*  (a) class return map: (mean pre, mean post) for each of the 14 Roxas classes (64 students each).
#*  (b) Hake gain vs pre-score per class, against the gain *shapes* the model predicts: TI constant,
#*      PI g ∝ a₀ − a₂ f. Amplitudes are fitted (one free scale each), because one simulation step
#*      is not one assessment interval (for small gains, k steps ≈ k × the one-step gain).
#* Bootstrap intervals (resampling each class's 64 students) are computed for the caption only; they
#* made the panel too busy to plot.
#* Class 9 is tested as an outlier (Grubbs, robust z) and excluded from the fits; the statistics
#* are returned as caption text.
function class_bootstrap(a, b; n_boot = 5000, rng = StableRNG(2026))
    n = length(a)
    pre_b = zeros(n_boot); g_b = zeros(n_boot)
    for k in 1:n_boot
        idx = rand(rng, 1:n, n)
        p, q = mean(a[idx]), mean(b[idx])
        pre_b[k] = p; g_b[k] = (q - p) / (1 - p)
    end
    return quantile(pre_b, (0.025, 0.975)), quantile(g_b, (0.025, 0.975))
end

#* Two-sided Grubbs test for the single most extreme value
function grubbs(x; α = 0.05)
    n = length(x)
    G = maximum(abs.(x .- mean(x))) / std(x)
    t = quantile(TDist(n - 2), 1 - α / (2n))
    G_crit = (n - 1) / sqrt(n) * sqrt(t^2 / (n - 2 + t^2))
    return G, G_crit
end
robust_z(x, i) = (x[i] - median(x)) / (1.4826 * median(abs.(x .- median(x))))

function figure6_option2(summary, emp)
    pre = [mean(a) for (a, _) in emp]; post = [mean(b) for (_, b) in emp]
    g = (post .- pre) ./ (1 .- pre)
    n_dec = [count(b .< a) for (a, b) in emp]
    out = argmin(g) #* class 9
    keep = setdiff(1:length(g), out)
    cis = [class_bootstrap(a, b) for (a, b) in emp]
    fs = 0:0.01:1
    marker_kw = (; marker = :diamond, markersize = 24)
    open_kw = (; color = :white, strokecolor = COLORS[5], strokewidth = 3)

    fig = Figure(size = (2000, 820))
    ax = v1_axis(fig[1, 1]; xlabel = "Class mean, first assessment", ylabel = "Class mean, second assessment",
        title = "(a) Class-level return map: one point per Roxas class", titlesize = 28, limits = (0, 1, 0, 1))
    lines!(ax, fs, fs; color = COLORS[8], linewidth = 3, label = "y = x")
    lines!(ax, fs, nitta(fs); color = COLORS[3], linewidth = 3, label = "Nitta: g(f) = f")
    scatter!(ax, pre[keep], post[keep]; color = COLORS[5], marker_kw..., label = "Class mean (Roxas): 14 classes, 64 students each")
    scatter!(ax, pre[[out]], post[[out]]; open_kw..., marker_kw..., label = "Class $(out): $(n_dec[out]) of its 64 students scored lower")
    #* free quadratic over the class means (class 9 excluded), drawn over the data range only;
    #* its f² coefficient a₂ is the class-level analogue of the published student-level fit
    qfit = curve_fit(empirical_model, pre[keep], post[keep], [0.0, 1.0, 0.0])
    a2c, a2c_se = qfit.param[1], stderror(qfit)[1]
    fq = range(minimum(pre[keep]), maximum(pre[keep]); length = 100)
    lines!(ax, fq, empirical_model(fq, qfit.param); color = COLORS[5], linestyle = :dash, linewidth = 3,
        label = "Quadratic fit, 13 classes: a₂ = $(round(a2c; digits = 2)) ± $(round(a2c_se; digits = 2))")
    axislegend(ax; position = :lt, labelsize = 20, framevisible = false)

    ax = v1_axis(fig[1, 2]; xlabel = "Class mean, first assessment", ylabel = "Class Hake gain g",
        title = "(b) Gain vs pre-score: model shapes, amplitude fitted", titlesize = 28)
    hlines!(ax, [0]; color = (:black, 0.3))
    scatter!(ax, pre[keep], g[keep]; color = COLORS[5], marker_kw..., label = "Class mean (Roxas): Hake g of each class")
    scatter!(ax, pre[[out]], g[[out]]; open_kw..., marker_kw..., label = "Class $(out): outlier, excluded from fits")
    shapes = Any[("TI: constant gain", f -> one(f), COLORS[1], :solid)]
    #* all four Fig. 6 cases (unfavorable/favorable × L = 32/128); the shapes nearly coincide because a₀ ≪ |a₂|
    for ((L, ρ₀, δλ), ls) in zip(FIG6_CASES, (:dash, :dot, :dashdot, (:dot, :loose)))
        a0, a2 = mean_coefs(summary, "PI", L, ρ₀, δλ)
        push!(shapes, ("PI: g ∝ a₀ − a₂f (L = $L, ρ₀ = $ρ₀, λ = 0.5 ± $δλ)", f -> a0 - a2 * f, COLORS[2], ls))
    end
    #* least-squares amplitude c of the model c · shape(f) through the origin
    amplitude(y, shape) = sum(y .* shape) / sum(shape .^ 2)
    r2_lines = String[]
    for (name, s, col, ls) in shapes
        shape_kept = s.(pre[keep]); amplitude_kept = amplitude(g[keep], shape_kept)
        lines!(ax, fs, amplitude_kept .* s.(fs); color = col, linestyle = ls, linewidth = 3, label = name)
        shape_all = s.(pre); amplitude_all = amplitude(g, shape_all)
        push!(r2_lines, "  $name: R² = $(round(r_squared(g[keep], amplitude_kept .* shape_kept); digits = 2)) without class $(out) " *
            "($(round(r_squared(g, amplitude_all .* shape_all); digits = 2)) with it)")
    end
    axislegend(ax; position = :lt, labelsize = 18, framevisible = false)

    G, G_crit = grubbs(g)
    caption = """
    Figure 6 option 2 — numbers for the caption / text

    Data: 14 classes (Roxas), 64 students each; class value = mean score; Hake gain g = (post − pre)/(1 − pre).
    R² = 1 − SS_res/SS_tot over the class gains (13 classes, class 9 excluded): the share of the class-to-class
    variation in g explained by each model shape, with one fitted amplitude per shape (a constant line can only
    sit at the mean gain, so TI's R² is 0 by construction).

    Fitted shapes (amplitude fitted; only the shape is compared):
    $(join(r2_lines, "\n"))
    Correlation between pre-score and gain: r = $(round(cor(pre[keep], g[keep]); digits = 2)) without class $(out), $(round(cor(pre, g); digits = 2)) with it.
    Class-level free quadratic (panel a, class $(out) excluded): a₂ = $(round(a2c; digits = 3)) ± $(round(a2c_se; digits = 3)).

    Class $(out) as an outlier:
      Hake gain $(round(g[out]; digits = 3)); the other 13 classes range $(round(minimum(g[keep]); digits = 3)) to $(round(maximum(g[keep]); digits = 3)).
      Grubbs test (two-sided, α = 0.05, n = 14): G = $(round(G; digits = 2)) vs critical $(round(G_crit; digits = 2)) → $(G > G_crit ? "outlier" : "not an outlier").
      Robust z (median/MAD) of its gain: $(round(robust_z(g, out); digits = 1)).
      Students scoring lower on the second assessment: $(n_dec[out])/64, vs $(minimum(n_dec[keep]))–$(maximum(n_dec[keep])) in the other classes (robust z = $(round(robust_z(Float64.(n_dec), out); digits = 1))).
      Its 95% bootstrap interval for g (5000 resamples of its 64 students): $(round.(cis[out][2]; digits = 3)), entirely below zero;
      it is the only class whose interval lies entirely below zero.
    """
    return fig, caption
end

function save_return_map_figures(summary)
    emp = read_empirical()
    emp_p = nothing #* student-level free quadratic
    if isempty(emp)
        @info "Empirical data not found in $(EMP_DIR); plotting simulations only"
    else
        pre, post = reduce(vcat, first.(emp)), reduce(vcat, last.(emp))
        emp_p = curve_fit(empirical_model, pre, post, [0.0, 1.0, 0.0]).param
    end
    dir = "$(OUTPUT_DIR)/analysis/plots/return-map"; mkpath(dir)
    save("$(dir)/figure6-option1-student-level.png", figure6_option1(summary, emp, emp_p); px_per_unit = 2)
    if !isempty(emp) #* class-level empirical figure (paper Fig. 9)
        fig2, caption = figure6_option2(summary, emp)
        save("$(dir)/figure6-option2-class-level.png", fig2; px_per_unit = 2)
        mkpath("$(OUTPUT_DIR)/analysis/plots/clean")
        save("$(OUTPUT_DIR)/analysis/plots/clean/figure6-option2-class-level-clean.png",
            strip_legends!(first(figure6_option2(summary, emp))); px_per_unit = 2)
        write("$(dir)/figure6-option2-caption.txt", caption)
        print(caption)
    end
    #* clean panels at the v1 source dimensions (1000×707 at 4 px/unit = 4000×2828)
    cdir = "$(OUTPUT_DIR)/analysis/plots/clean"; mkpath(cdir)
    for (L, ρ₀, δλ) in FIG6_CASES #* all four panels of the 2×2 option-1 grid
        fig = Figure(size = (1000, 1000 ÷ √2))
        return_map_panel!(v1_axis(fig[1, 1]), L, ρ₀, δλ, summary, emp, emp_p)
        save("$(cdir)/figure6-$(L)-$(ρ₀)-0.5-$(δλ)-clean.png", fig; px_per_unit = 4)
    end
    println("Saved return-map figures to $(dir) and clean panels to $(cdir)")
end

#* ---------------------------------------------------------------------------
#* Supplementary (reviewers R1, R3, §4): coefficients vs L at fixed (ρ₀, δλ), inner corner, for the
#* main set (4 seeds) and the seed-density set (fixed initial fraction 4/32², n = round(4L²/32²)).
#* δλ ≤ 0.4 only: the PI quadratic breaks down at 0.45/0.49.
#* ---------------------------------------------------------------------------
const VS_L_LS = [32, 48, 64, 96, 128]
n_frac(L) = round(Int, 4 * L^2 / 32^2)
loglog_slope(x, y) = (lx = log.(x); ly = log.(y); cov(lx, ly) / var(lx))

function return_map_vs_L(; Ρs = [0.1, 0.5, 0.9], δλs = [0.0, 0.2, 0.4], λ₀ = 0.5)
    out = DataFrame()
    for regime in ("TI", "PI"), seeds in ("4 seeds", "fixed fraction"), L in VS_L_LS, ρ₀ in Ρs, δλ in δλs
        n_learned = seeds == "4 seeds" ? 4 : n_frac(L)
        row = param_row(regime, "inner_corner", L, ρ₀, δλ; n_learned, λ₀)
        #* one row of fitted coefficients per trial
        coefs = reduce(vcat, [permutedims(first(fit_return_map(regime, x, y))) for (x, y) in trial_pairs(row)])
        a2 = regime == "TI" ? fill(NaN, N_TRIALS) : coefs[:, 2]
        push!(out, (; regime, seeds, L, n_learned = row.n_learned, ρ₀, δλ,
            a0_mean = mean(coefs[:, 1]), a0_se = std(coefs[:, 1]) / sqrt(N_TRIALS),
            a2_mean = mean(a2), a2_se = std(a2) / sqrt(N_TRIALS)))
    end
    return out
end

#* Rows = δλ; columns = PI a₀, PI −a₂ (log–log), TI a₀/ρ₀λ₀. Color = ρ₀; solid = 4 seeds, dashed = fixed fraction.
#* Error bars are ±SE over the 20 per-run fits (mostly smaller than the markers). Returns (fig, slopes).
function figure_coefs_vs_L(vsL; λ₀ = 0.5, title = true)
    Ρs = sort(unique(vsL.ρ₀)); δλs = sort(unique(vsL.δλ))
    sel(regime, seeds, ρ₀, δλ) = sort(rows_where(vsL; regime, seeds, ρ₀, δλ), :L)
    slopes = DataFrame()
    fig = Figure(size = (1700, 1500))
    cols = (("PI: a₀", rows -> rows.a0_mean, rows -> rows.a0_se, true), ("PI: −a₂", rows -> -rows.a2_mean, rows -> rows.a2_se, true),
        ("TI: a₀ / ρ₀λ₀", rows -> rows.a0_mean ./ (rows.ρ₀ .* λ₀), rows -> rows.a0_se ./ (rows.ρ₀ .* λ₀), false))
    for (i, δλ) in enumerate(δλs), (j, (name, val, err, logy)) in enumerate(cols)
        ax = Axis(fig[i, j]; xscale = log2, yscale = logy ? log10 : identity, xticks = (VS_L_LS, string.(VS_L_LS)),
            title = "$name,  δλ = $δλ", titlesize = 22, xlabel = i == length(δλs) ? "Lattice side length L" : "",
            xlabelsize = 22, xticklabelsize = 17, yticklabelsize = 17)
        logy || ylims!(ax, 0.4, 1.1)
        for (k, ρ₀) in enumerate(Ρs), seeds in ("4 seeds", "fixed fraction")
            rows = sel(j == 3 ? "TI" : "PI", seeds, ρ₀, δλ)
            y = val(rows); errors = err(rows)
            if logy && any(y .<= 0)
                @warn "non-positive $(name) at ρ₀ = $ρ₀, δλ = $δλ ($seeds); skipped on log axis"
                continue
            end
            s = logy ? loglog_slope(rows.L, y) : NaN
            logy && push!(slopes, (; coefficient = name, seeds, ρ₀, δλ, loglog_slope = s))
            errorbars!(ax, rows.L, y, errors; color = COLORS[k], whiskerwidth = 6)
            scatterlines!(ax, rows.L, y; color = COLORS[k], linewidth = 2.5, linestyle = seeds == "4 seeds" ? :solid : :dash,
                marker = seeds == "4 seeds" ? :circle : :utriangle, markersize = 11,
                label = "ρ₀ = $ρ₀, $seeds" * (logy ? " (slope $(round(s; digits = 2)))" : ""))
        end
        axislegend(ax; position = j == 3 ? :rb : :lb, labelsize = 13, framevisible = false, rowgap = 0)
    end
    title && Label(fig[0, :], "Return-map coefficients vs lattice size at fixed ρ₀ and δλ (λ = 0.5 ± δλ)
Inner-corner seating; 4 seeds (solid) or fixed seed fraction 4/32² (dashed)";
        fontsize = 26)
    return fig, slopes
end

function save_coefs_vs_L()
    dir = "$(OUTPUT_DIR)/analysis"
    vsL = return_map_vs_L()
    CSV.write("$(dir)/return-map-vs-L.csv", vsL)
    fig, slopes = figure_coefs_vs_L(vsL)
    CSV.write("$(dir)/return-map-vs-L-slopes.csv", slopes)
    save("$(dir)/plots/return-map/coefficients-vs-L.png", fig; px_per_unit = 2)
    save("$(dir)/plots/clean/coefficients-vs-L-clean.png", strip_legends!(figure_coefs_vs_L(vsL; title = false)[1]); px_per_unit = 2)
    show(slopes; allrows = true); println()
end


#! Main
begin
    per_trial, summary = return_map_fits()
    dir = "$(OUTPUT_DIR)/analysis"
    CSV.write("$(dir)/return-map-per-trial.csv", per_trial)
    cols = [:regime, :L, :ρ₀, :δλ, :a0_mean, :a0_se, :a0_sd, :a0_pooled, :a0_pooled_naive_se,
        :a2_mean, :a2_se, :a2_sd, :a2_pooled, :a2_pooled_naive_se, :mean_r_squared, :pooled_n_points]
    summary = select(summary, cols, Not(cols))
    CSV.write("$(dir)/return-map-summary.csv", summary)
    println("Wrote $(nrow(per_trial)) per-trial fits and $(nrow(summary)) summary rows to $(dir)")
    save_return_map_figures(summary)
    save_coefs_vs_L()
end
