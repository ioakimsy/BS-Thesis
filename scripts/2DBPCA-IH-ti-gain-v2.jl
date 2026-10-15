# ! USE BS THESIS ENVIRONMENT (file paths will break if not)
#=
TI per-step gain vs reviewer R3's exact broadcast benchmark (paper §3.4, Eq. 6; review synthesis §11).

Under TI every agent learns independently with p± = ρ₀(λ₀ ± δλ) per step, so the unlearned
fraction is exactly
  u(t) = (1 − s₀) [w (1 − p₊)ᵗ + (1 − w)(1 − p₋)ᵗ]       (s₀ = seed fraction, w = fast share ≈ ½)
and the per-step normalized gain (Hake gain of one step) is
  g(t) = [f(t+1) − f(t)] / [1 − f(t)] = 1 − u(t+1)/u(t),
the unlearned-pool average of p₊ and p₋. It starts at ρ₀λ₀ (pool half fast, half slow) and decays
to p₋ = ρ₀(λ₀ − δλ) as the fast learners are used up. Eq. 6 (constant gain a₀) is exact only at δλ = 0.

Simulation gain: all (run, step) pairs are binned by f(t); a bin's gain is
Σ newly learned / Σ unlearned over its pairs (ratio estimator, so near-absorbed steps with few
unlearned agents don't dominate).

Outputs (output/2D-Binary-PCA-IH-v2/analysis/):
  ti-gain-drift.csv                       per (L, ρ₀, δλ): fitted a₀ vs the exact gain's range,
                                          the f-interval where the exact gain is within ±10% of a₀,
                                          and the largest simulation − exact gain difference
  plots/return-map/ti-gain-vs-exact.png   normalized gain vs f: simulation, exact curve, fitted a₀
=#

begin
    using Pkg
    Pkg.activate(".")
    include("2DBPCA-IH-common-v2.jl")
    using Statistics
    using CairoMakie, ColorSchemes
end

const COLORS = ColorSchemes.seaborn_colorblind.colors
const λ₀ = 0.5
const W = 0.5 #* expected fast share (λ classes are Bernoulli(½))

#* Exact unlearned fraction and per-step gain at (real) time t
exact_u(t, ρ₀, δλ, s₀) = (1 - s₀) * (W * (1 - ρ₀ * (λ₀ + δλ))^t + (1 - W) * (1 - ρ₀ * (λ₀ - δλ))^t)
exact_gain(t, ρ₀, δλ, s₀) = 1 - exact_u(t + 1, ρ₀, δλ, s₀) / exact_u(t, ρ₀, δλ, s₀)

#* Exact (f, g) curve, sampled finely in t until 1 − f < u_min
function exact_curve(ρ₀, δλ, s₀; u_min = 1e-4)
    fs, gs = Float64[], Float64[]
    t = 0.0
    dt = 0.02 / ρ₀
    while exact_u(t, ρ₀, δλ, s₀) > u_min
        push!(fs, 1 - exact_u(t, ρ₀, δλ, s₀)); push!(gs, exact_gain(t, ρ₀, δλ, s₀))
        t += dt
    end
    return fs, gs
end

#* Simulation gain binned by f (ratio estimator over all runs)
function simulated_gain(L, ρ₀, δλ; nbins = 40)
    row = (; regime = "TI", SA = "inner_corner", L, n_learned = 4, ρ₀, λ₀, δλ)
    raw = read_parquet(raw_path(row))
    new = zeros(nbins); unl = zeros(nbins); fsum = zeros(nbins); npair = zeros(nbins)
    for tr in 1:N_TRIALS
        f = learned_fraction(read_trial(row, tr; raw)[1])
        for t in 1:length(f)-1
            b = clamp(floor(Int, f[t] * nbins) + 1, 1, nbins)
            new[b] += f[t+1] - f[t]; unl[b] += 1 - f[t]
            fsum[b] += f[t]; npair[b] += 1
        end
    end
    keep = npair .> 0
    return fsum[keep] ./ npair[keep], new[keep] ./ unl[keep], npair[keep]
end

function gain_drift(summary; Ls = [32, 64, 128], Ρs = [0.1, 0.5, 0.9], δλs = [0.0, 0.2, 0.4, 0.45, 0.49], tol = 0.10)
    out = DataFrame()
    for L in Ls, ρ₀ in Ρs, δλ in δλs
        s₀ = 4 / L^2
        r = summary[(summary.regime .== "TI") .& (summary.L .== L) .& isapprox.(summary.ρ₀, ρ₀) .& isapprox.(summary.δλ, δλ), :]
        a0 = r.a0_mean[1]
        fs, gs = exact_curve(ρ₀, δλ, s₀)
        within = abs.(gs .- a0) ./ a0 .<= tol
        f_lo = any(within) ? minimum(fs[within]) : NaN
        f_hi = any(within) ? maximum(fs[within]) : NaN
        sf, sg, n = simulated_gain(L, ρ₀, δλ)
        ok = n .>= 20 #* bins with enough pairs to compare
        sim_err = maximum(abs.(sg[ok] .- exact_gain_at.(sf[ok], ρ₀, δλ, s₀)) ./ (ρ₀ * λ₀))
        push!(out, (; L, ρ₀, δλ, rho_lambda0 = ρ₀ * λ₀, a0_fit = a0, a0_over_rho_lambda0 = a0 / (ρ₀ * λ₀),
            gain_start = gs[1], gain_end = ρ₀ * (λ₀ - δλ), start_over_end = gs[1] / (ρ₀ * (λ₀ - δλ)),
            max_rel_dev_from_a0 = maximum(abs.(gs .- a0) ./ a0),
            f_within_10pct_lo = f_lo, f_within_10pct_hi = f_hi,
            max_sim_minus_exact_over_rho_lambda0 = sim_err))
    end
    return out
end

#* Exact gain at a given f (invert f(t) by bisection; f is increasing in t)
function exact_gain_at(f, ρ₀, δλ, s₀)
    u = 1 - f
    u >= 1 - s₀ && return exact_gain(0.0, ρ₀, δλ, s₀)
    lo, hi = 0.0, 1.0
    while exact_u(hi, ρ₀, δλ, s₀) > u; hi *= 2; end
    for _ in 1:60
        mid = (lo + hi) / 2
        exact_u(mid, ρ₀, δλ, s₀) > u ? (lo = mid) : (hi = mid)
    end
    return exact_gain((lo + hi) / 2, ρ₀, δλ, s₀)
end

function ti_gain_figure(summary; L = 128, Ρs = [0.1, 0.5, 0.9], δλs = [0.0, 0.2, 0.4, 0.49])
    s₀ = 4 / L^2
    fig = Figure(size = (2200, 900))
    markers = (:circle, :rect, :utriangle)
    for (j, δλ) in enumerate(δλs)
        ax = Axis(fig[1, j]; title = "δλ = $(δλ)", titlesize = 26, ylabel = j == 1 ? "Per-step gain / ρ₀λ₀" : "",
            xticklabelsize = 20, yticklabelsize = 20, ylabelsize = 24, limits = (0, 1, 0, 1.12))
        axr = Axis(fig[2, j]; xlabel = "Learned fraction f(t)", ylabel = j == 1 ? "(simulation − exact) / ρ₀λ₀" : "",
            xticklabelsize = 20, yticklabelsize = 20, xlabelsize = 24, ylabelsize = 22, limits = (0, 1, -0.06, 0.06))
        hidexdecorations!(ax; grid = false, ticks = false)
        hlines!(axr, [0]; color = (:black, 0.4))
        for (i, ρ₀) in enumerate(Ρs)
            norm = ρ₀ * λ₀
            ef, eg = exact_curve(ρ₀, δλ, s₀)
            lines!(ax, ef, eg ./ norm; color = COLORS[i], linewidth = 3, label = "Exact (R3), ρ₀ = $(ρ₀)")
            a0 = summary[(summary.regime .== "TI") .& (summary.L .== L) .& isapprox.(summary.ρ₀, ρ₀) .& isapprox.(summary.δλ, δλ), :a0_mean][1]
            hlines!(ax, [a0 / norm]; color = COLORS[i], linestyle = :dash, linewidth = 2.5, label = "Fitted constant a₀ (Eq. 6), ρ₀ = $(ρ₀)")
            sf, sg, n = simulated_gain(L, ρ₀, δλ)
            ok = n .>= 20
            scatter!(ax, sf[ok], sg[ok] ./ norm; color = (COLORS[i], 0.7), marker = markers[i], markersize = 12,
                label = "Simulation (20 runs), ρ₀ = $(ρ₀)")
            scatter!(axr, sf[ok], (sg[ok] .- exact_gain_at.(sf[ok], ρ₀, δλ, s₀)) ./ norm; color = (COLORS[i], 0.7),
                marker = markers[i], markersize = 12)
        end
        j == 1 && (global leg_ax = ax)
    end
    Legend(fig[1:2, length(δλs)+1], leg_ax; merge = true, labelsize = 18, framevisible = false)
    rowsize!(fig.layout, 2, Relative(0.3))
    Label(fig[0, 1:length(δλs)], "TI per-step gain vs R3's exact broadcast result (L = $L, λ₀ = 0.5)"; fontsize = 28)
    return fig
end

#! Main
begin
    dir = "$(OUTPUT_DIR)/analysis"
    summary = CSV.read("$(dir)/return-map-summary.csv", DataFrame)
    drift = gain_drift(summary)
    CSV.write("$(dir)/ti-gain-drift.csv", drift)
    show(drift[drift.L .== 128, [:ρ₀, :δλ, :a0_over_rho_lambda0, :start_over_end, :max_rel_dev_from_a0,
        :f_within_10pct_lo, :f_within_10pct_hi, :max_sim_minus_exact_over_rho_lambda0]]; allrows = true)
    println()
    println("Max |a0(L) − a0(128)| / a0(128) across L: ",
        maximum(maximum(abs.(g.a0_fit .- g.a0_fit[g.L .== 128][1]) ./ g.a0_fit[g.L .== 128][1])
            for g in groupby(drift, [:ρ₀, :δλ])))
    mkpath("$(dir)/plots/return-map")
    save("$(dir)/plots/return-map/ti-gain-vs-exact.png", ti_gain_figure(summary); px_per_unit = 2)
end
